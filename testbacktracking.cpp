#include <utility>
#include <chrono>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
// #include <seqan3/alignment/all.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>

#include <seqan3/search/all.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>   // for using the bi_fm_index
#include <seqan3/search/fm_index/fm_index.hpp>      // for using the fm_index
#include <seqan3/search/algorithm/search.hpp>
#include <type_traits>

#include <range/v3/view/drop_exactly.hpp>

#include <seqan3/range/concept.hpp>
#include <seqan3/search/algorithm/detail/search_common.hpp>

#include <seqan3/core/metafunction/pre.hpp>
#include <seqan3/search/algorithm/detail/search_trivial.hpp>
#include <seqan3/search/fm_index/concept.hpp>

#include <mutex>
#include <ctime>


using namespace seqan3;
using namespace std::string_literals; // for using the ""s string literal
using namespace seqan3::search_cfg;

namespace seqan3::test
{

struct search_param
{
    uint8_t total, substitution, insertion, deletion;
};

} // namespace seqan3::detail


namespace seqan3::test
{

template <bool abort_on_hit, typename query_t, typename cursor_t, typename delegate_t>
inline bool search_trivial(cursor_t cur, query_t & query, typename cursor_t::size_type const query_pos,
                           search_param const error_left, delegate_t && delegate) noexcept(noexcept(delegate))
{
    // Exact case (end of query sequence or no errors left)
    if (query_pos == std::ranges::size(query) || error_left.total == 0)
    {
        // If not at end of query sequence, try searching the remaining suffix without any errors.
        if (query_pos == std::ranges::size(query) || cur.extend_right(std::view::drop(query, query_pos)))
        {
            delegate(cur);
            return true;
        }
    }
    // Approximate case
    else
    {
        // Insertion
        if (error_left.insertion > 0)
        {
            search_param error_left2{error_left};
            error_left2.insertion--;
            error_left2.total--;

            // always perform a recursive call. Abort recursion if and only if recursive call found a hit and
            // abort_on_hit is set to true.
            if (search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left2, delegate) && abort_on_hit)
                return true;
        }

        // Do not allow deletions at the beginning of the query sequence
        if (((query_pos > 0 && error_left.deletion > 0) || error_left.substitution > 0) && cur.extend_right())
        {
            do
            {
                // Match (when error_left.substitution > 0) and Mismatch
                if (error_left.substitution > 0)
                {
                    bool delta = cur.last_char() != query[query_pos];
                    search_param error_left2{error_left};
                    error_left2.total -= delta;
                    error_left2.substitution -= delta;

                    if (search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left2, delegate) && abort_on_hit)
                        return true;
                }

                // Deletion (Do not allow deletions at the beginning of the query sequence.)
                if (query_pos > 0)
                {
                    // Match (when error_left.substitution == 0)
                    if (error_left.substitution == 0 && cur.last_char() == query[query_pos])
                    {
                        if (search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left, delegate) &&
                            abort_on_hit)
                        {
                            return true;
                        }
                    }

                    // Deletions at the end of the sequence are not allowed. This cannot happen: when the algorithm
                    // arrives here, it cannot be at the end of the query and since deletions do not touch the query
                    // (i.e. increase query_pos) it won't be at the end of the query after the deletion.
                    if (error_left.deletion > 0)
                    {
                        search_param error_left2{error_left};
                        error_left2.total--;
                        error_left2.deletion--;

                        if (search_trivial<abort_on_hit>(cur, query, query_pos, error_left2, delegate) && abort_on_hit)
                            return true;
                    }
                }
            } while (cur.cycle_back());
        }
        else
        {
            // Match (when error_left.substitution == 0)
            if (cur.extend_right(query[query_pos]))
            {
                if (search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left, delegate) && abort_on_hit)
                    return true;
            }
        }
    }

    return false;
}

template <bool abort_on_hit, typename index_t, typename query_t, typename delegate_t>
inline void search_trivial(index_t const & index, query_t & query, search_param const error_left,
                           delegate_t && delegate) noexcept(noexcept(delegate))
{
    search_trivial<abort_on_hit>(index.begin(), query, 0, error_left, delegate);
}


enum class ErrorCode : std::uint8_t {
    LAST_MATCH = 0, LAST_INSERTION = 1, LAST_NOTHING = 2, LAST_DELETION = 3, LAST_MISMATCH = 4
};

template <bool abort_on_hit, typename query_t, typename cursor_t, typename delegate_t>
inline bool my_search_trivial(cursor_t cur, query_t & query, typename cursor_t::size_type const query_pos,
                           search_param const error_left, ErrorCode memory,
                           delegate_t && delegate) noexcept(noexcept(delegate))
{
    // Exact case (end of query sequence or no errors left)
    if (query_pos == std::ranges::size(query) || error_left.total == 0)
    {
        // If not at end of query sequence, try searching the remaining suffix without any errors.
        if (query_pos == std::ranges::size(query) || cur.extend_right(std::view::drop(query, query_pos)))
        {
            delegate(cur);
            return true;
        }
    }
    // Approximate case
    else
    {
        // Do not use insertion if current char is the same as previous   AAAA not AAAA
        //                                                                IAAA not AIAA
        // Cannot check for match since the insertion alignment was not enumerated as is the case for deletion
        // Check if cursor did a step in the index.
        bool const do_insertion = cur.query_length() > 0 ? cur.last_char() != query[query_pos] : true;

        // Insertion
        // Do not allow insertions after an deletion or mismatch.
        if (do_insertion && (memory <= ErrorCode::LAST_NOTHING
            || error_left.substitution == 0) && error_left.insertion > 0)
        {
            search_param error_left2{error_left};
            error_left2.insertion--;
            error_left2.total--;

            // always perform a recursive call. Abort recursion if and only if recursive call found a hit and
            // abort_on_hit is set to true.
            if (my_search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left2, ErrorCode::LAST_INSERTION, delegate) && abort_on_hit)
                return true;
        }

        // Do not allow deletions at the beginning of the query sequence
        if (query_pos > 0 && ((error_left.deletion > 0) || error_left.substitution > 0) && cur.extend_right())
        {
            do
            {
                //(Do not allow deletions or mismatches at the beginning of the query sequence.)
//                 if (query_pos > 0)
//                 {
                    // Match (when error_left.substitution > 0) and Mismatch
                    if (error_left.substitution > 0)
                    {
                        bool delta = cur.last_char() != query[query_pos];
                        search_param error_left2{error_left};
                        error_left2.total -= delta;
                        error_left2.substitution -= delta;

                        ErrorCode errorCode = (delta) ? ErrorCode::LAST_MISMATCH : ErrorCode::LAST_MATCH;

                        if (my_search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left2, errorCode, delegate) && abort_on_hit)
                            return true;
                    }

                    // Deletion
                    // Match (when error_left.substitution == 0)
                    if (error_left.substitution == 0 && cur.last_char() == query[query_pos])
                    {
                        if (my_search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left, ErrorCode::LAST_MATCH, delegate) &&
                            abort_on_hit)
                        {
                            return true;
                        }
                    }

                    // Deletions at the end of the sequence are not allowed. This cannot happen: when the algorithm
                    // arrives here, it cannot be at the end of the query and since deletions do not touch the query
                    // (i.e. increase query_pos) it won't be at the end of the query after the deletion.
                    // Do not allow deletions after an insertion.
                    if ((memory == ErrorCode::LAST_MATCH || memory == ErrorCode::LAST_DELETION || error_left.substitution == 0) &&
                        error_left.deletion > 0 && cur.last_char() != query[query_pos])
                    {
                        search_param error_left2{error_left};
                        error_left2.total--;
                        error_left2.deletion--;
                        if (my_search_trivial<abort_on_hit>(cur, query, query_pos, error_left2, ErrorCode::LAST_DELETION, delegate) && abort_on_hit)
                            return true;
//                     }
                }
            } while (cur.cycle_back());
        }
        else
        {
            // Match (when error_left.substitution == 0)
            if (cur.extend_right(query[query_pos]) || query_pos == 0)
            {
                if (my_search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left, ErrorCode::LAST_MATCH, delegate) && abort_on_hit)
                    return true;
            }
        }
    }

    return false;
}

template <bool abort_on_hit, typename index_t, typename query_t, typename delegate_t>
inline void my_search_trivial(index_t const & index, query_t & query, search_param const error_left,
                           delegate_t && delegate) noexcept(noexcept(delegate))
{
    my_search_trivial<abort_on_hit>(index.begin(), query, 0, error_left, ErrorCode::LAST_NOTHING, delegate);
}


}

// template<typename Tunit>
struct Timer{
    double defaultTime = 0.0;
    double myTime = 0.0;
//     std::chrono::time_point<std::chrono::high_resolution_clock> start;
//     std::chrono::time_point<std::chrono::high_resolution_clock> end;

    void addDefault(double t)
    {
//         std::chrono::duration<double> elapsedT = end - start;
        defaultTime += t;
    }

    void addMy(double t)
    {
//         std::chrono::duration<double> elapsedT = end - start;
        myTime += t;
    }
};


template <typename index_t>
void my_search(index_t const & index, /*queries_t */auto & queries, uint8_t errors, auto & randomtext, Timer & timer, bool verbose = false
){
    test::search_param max_errorNormal{errors, errors, errors, errors};
    test::search_param max_error{errors, errors, errors, errors};

    std::vector<typename index_t::cursor_type> internal_hits;
    std::vector<typename index_t::cursor_type> my_internal_hits;
//     std::vector<int> internal_hits(5);
    auto internal_delegate = [&internal_hits, &max_error](auto const & cur)
    {
//         std::cout << "found something\n";
        internal_hits.push_back(cur);
//         for (auto const & pos : cur.locate())
//              debug_stream << pos << ' ';
    };

    auto my_internal_delegate = [&my_internal_hits, &max_error](auto const & cur)
    {
        my_internal_hits.push_back(cur);
//         for (auto const & pos : cur.locate())
//             debug_stream << pos << ' ';
    };

    using hit_t = std::conditional_t<index_t::is_collection,
                                         std::pair<typename index_t::size_type, typename index_t::size_type>,
                                         typename index_t::size_type>;

    std::vector<std::vector<hit_t> > hits(queries.size());
    std::vector<std::vector<hit_t> > myhits(queries.size());


    std::chrono::time_point<std::chrono::high_resolution_clock> start, end, startmy, endmy;
//     std::chrono::time_point<std::chrono::high_resolution_clock> end;

    //normal
    start = std::chrono::high_resolution_clock::now();
    int k = 0;
    for (auto const query : queries)
    {
//         std::cout << "Searching\n";
        test::search_trivial<false>(index, query, max_errorNormal, internal_delegate);

        if(verbose)
            std::cout << "Number of backtracked positions: " << internal_hits.size() << "\n";

        /*
        if(internal_hits.size() == 0){
            std::cout << k << "\n";
            debug_stream << "Something went wrong?\n";
            debug_stream << randomtext << "\n";
            for(int i = 0; i < errors; ++i)
                debug_stream << " ";
            debug_stream << query << "\n";
            exit(1);
        }*/
        for (auto const & cur : internal_hits)
        {
            for (auto const & text_pos : cur.locate())
                hits[k].push_back(text_pos);
            std::sort(hits[k].begin(), hits[k].end());
            hits[k].erase(std::unique(hits[k].begin(), hits[k].end()), hits[k].end());
        }
        internal_hits.clear();
        ++k;
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    timer.addDefault(elapsed.count());
    std::cout << "Searched normal\n";
    startmy = std::chrono::high_resolution_clock::now();
    int k2 = 0;
    auto startDefault = std::chrono::high_resolution_clock::now();
    for (auto const query : queries)
    {
        test::my_search_trivial<false>(index, query, max_error, my_internal_delegate);
        if(verbose)
            std::cout << "Number of backtracked positions: " << my_internal_hits.size()  << "\t My\n";
/*
        if(my_internal_hits.size() == 0){
            std::cout << k2 << "\n";
            debug_stream << "Something went wrong? My Version\n";
            debug_stream << randomtext << "\n";
            for(int i = 0; i < errors; ++i)
                debug_stream << " ";
            debug_stream << query << "\n";
            exit(1);
        }*/
        for (auto const & cur : my_internal_hits)
        {
            for (auto const & text_pos : cur.locate()){
                myhits[k2].push_back(text_pos);
            }
            std::sort(myhits[k2].begin(), myhits[k2].end());
            myhits[k2].erase(std::unique(myhits[k2].begin(), myhits[k2].end()), myhits[k2].end());

        }
        my_internal_hits.clear();
        ++k2;
    }
    endmy = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedmy = endmy - startmy;
    timer.addMy(elapsedmy.count());
/*
    std::cout << "Checking Occurrences\n";
    for(int l = 0; l < hits.size(); ++l)
    {
        if(hits[l].size() != myhits[l].size()){
            for(auto & pos : hits[l])
            debug_stream << pos << "\n";
            std::cout << "myHits\n\n";
            for(auto & pos : myhits[l])
                debug_stream << pos << "\n";
            exit(0);
        }
        for(int m = 0; m < hits[l].size(); ++m)
        {
            if(hits[l][m] != myhits[l][m])
            {
                std::cout << "Something went wrong at read: " << l << " occ: " << m << "\n";
                debug_stream  << hits[l] << "\t myhits: " << myhits[l] << "\n";
                exit(0);
            }

        }
    }*/

    hits.clear();
    myhits.clear();
}

template <Alphabet alphabet_t>
auto generate_sequence_seqan3(size_t const len = 500,
                              size_t const variance = 0,
                              size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, alphabet_size_v<alphabet_t> - 1);
    std::uniform_int_distribution<size_t> dis_length(len - variance, len + variance);

    std::vector<alphabet_t> sequence;

    size_t length = dis_length(gen);
    for (size_t l = 0; l < length; ++l)
        sequence.push_back(alphabet_t{}.assign_rank(dis_alpha(gen)));

    return sequence;
}

template<Alphabet alphabet_t>
void mutate_insertion(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed = 0){
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, alphabet_size_v<alphabet_t> - 1);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(seq) - overlap);
    alphabet_t cbase;
    seq.insert(seq.begin() + random_pos(gen), alphabet_t{}.assign_rank(dis_alpha(gen)));
}

template<Alphabet alphabet_t>
void mutate_deletion(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed = 0){
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(seq) - overlap);
    seq.erase(seq.begin() + random_pos(gen));
}

template<Alphabet alphabet_t>
void mutate_substitution(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed = 0){
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha_short(0, alphabet_size_v<alphabet_t> - 2);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(seq) - overlap);
    alphabet_t & cbase = seq[random_pos(gen)];
    uint8_t crank = to_rank(cbase);
    uint8_t rrank = dis_alpha_short(gen);
    if (rrank >=  crank)
        ++rrank;
    cbase.assign_rank(rrank);
}

template<Alphabet alphabet_t>
auto generate_reads(std::vector<alphabet_t> & ref,
                    size_t const number_of_reads,
                    size_t const read_length,
                    size_t const simulated_errors,
                    float const prob_insertion,
                    float const prob_deletion,
                    size_t const seed = 0)
{
    std::vector<std::vector<alphabet_t> > reads;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> seeds (0, SIZE_MAX);
    std::uniform_int_distribution<size_t> random_pos(simulated_errors, std::ranges::size(ref) - read_length - simulated_errors);
    std::uniform_real_distribution<double> probability(0.0, 1.0);
    for (size_t i = 0; i < number_of_reads; ++i){
        size_t rpos = random_pos(gen);
        std::vector<alphabet_t> read_tmp{ref.begin() + rpos,
            ref.begin() + rpos + read_length + simulated_errors};
        for (size_t j = 0; j < simulated_errors; ++j)
        {
            double prob = probability(gen);
            //Substitution
            if (prob_insertion + prob_deletion < prob)
            {
                mutate_substitution(read_tmp, simulated_errors, seeds(gen));
            }
            //Insertion
            else if (prob_insertion < prob)
            {
                mutate_insertion(read_tmp, simulated_errors, seeds(gen));
            }
            //Deletion
            else
            {
                mutate_deletion(read_tmp, simulated_errors, seeds(gen));
            }
        }
        read_tmp.erase(read_tmp.begin() + read_length, read_tmp.end());
        reads.push_back(read_tmp);
    }
    return reads;
}

int main(int argc, char * argv[])
{
    
//     int simulatedErrors = 2;
    int searchErrors = 3;
//     uint32_t olength = 600;
    uint32_t olength = 100;
//     uint32_t olength = 40;
    float probI = 0.18;
    float probD = 0.18;
    int length = olength;
    int iterationText = 1;
//     int numberofReads = 20;
    int numberofReads = 60;
//     int numberofReads = 250;
    int indexSize = 100'000;
    length += searchErrors * 2 + indexSize;

    std::chrono::duration<double> elapseduni;
    std::chrono::duration<double> elapsedbi;
//     srand (time(NULL));
    for(uint8_t sE = 0; sE <= searchErrors; ++sE){
        Timer timer;
        std::cout << "simulated Errors " << (int)sE << "\n";
        for(int t = 0; t < iterationText; ++t){
//             dna4_vector randomtext{};
//             generateText(randomtext, length);
            dna4_vector randomtext = generate_sequence_seqan3<seqan3::dna4>(length);

            fm_index<dna4_vector> index{randomtext};
            std::cout << "built Index\n";

            std::vector<dna4_vector> reads = generate_reads(randomtext, numberofReads, olength, sE, probI, probD);
        //     std::span text_view{std::data(query) + searchErrors, olength + searchErrors};
        //     dna4_vector read{text_view.begin(), text_view.end()};
            std::cout << "simulated Reads\n";
            my_search(index, reads, searchErrors, randomtext, timer);

            
            std::cout << "uni Search\n";
            uint8_t et = searchErrors;
            configuration cfg = max_error{total{et}, deletion{et}, insertion{et}, substitution{et}};
            auto startUni = std::chrono::high_resolution_clock::now();
            auto resultsuni = search(index, reads, cfg);
            auto endUni = std::chrono::high_resolution_clock::now();
            elapseduni = endUni - startUni;
            
            std::cout << "bi Search\n";
            bi_fm_index<dna4_vector> biIndex{randomtext};
//             configuration cfg = max_error{total{sE}, deletion{sE}, insertion{sE}, substitution{sE}};
            auto startBi = std::chrono::high_resolution_clock::now();
            auto results = search(biIndex, reads, cfg);
            auto endBi = std::chrono::high_resolution_clock::now();
            elapsedbi = endBi - startBi;

        }

        std::cout << "Default Time: \t\t" << timer.defaultTime << "\n";
        std::cout << "My Time: \t\t" << timer.myTime << "\n";
        std::cout << "Uni Search Time: \t" << elapseduni.count() << "\n";
        std::cout << "Bi Search Time: \t" << elapsedbi.count() << "\n";//elapsedmy.count()
    }

    indexSize = 0;
    probI = 0.33;
    probD = 0.33;
    olength = 50;
    length = olength + searchErrors * 2;
    numberofReads = 100;

    std::cout << "Checking\n";
    bool did_not_find_every_read{false};
    for(size_t t = 0; t < 50; ++t){
        for(uint8_t e = searchErrors; e <= searchErrors; ++e){

            dna4_vector randomtext = generate_sequence_seqan3<seqan3::dna4>(length);
            fm_index<dna4_vector> index{randomtext};
            std::vector<dna4_vector> reads = generate_reads(randomtext, numberofReads, olength, e, probI, probD);

            uint8_t tmpe = searchErrors;
            test::search_param max_error{tmpe, tmpe, tmpe, tmpe};
            std::vector<typename fm_index<dna4_vector>::cursor_type> my_internal_hits;
            auto my_internal_delegate = [&my_internal_hits, &max_error](auto const & cur)
            {
                my_internal_hits.push_back(cur);
        //         for (auto const & pos : cur.locate())
        //             debug_stream << pos << ' ';
            };

//             std::cout << "simulated Errors " << (int)e << "\n";
//             std::cout << "search Errors " << searchErrors << "\n";
            for(auto & read : reads)
            {
                test::my_search_trivial<false>(index, read, max_error, my_internal_delegate);
//                 auto results = search(index, read, cfg);
//                 std::cout << "Size: " << my_internal_hits.size() << "\n";
                if(my_internal_hits.size() == 0){
                    debug_stream << "Print sequences\n";
                    debug_stream << randomtext << "\n";
                    for(int i = 0; i < searchErrors; ++i)
                        debug_stream << " ";
                    debug_stream << read << "\n";
                    did_not_find_every_read = true;
                    return 1;
                }
                my_internal_hits.clear();
            }
        }
    }

    std::cout << "fin\n";
    return 0;
}
