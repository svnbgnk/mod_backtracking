#include <utility>
#include <chrono>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
// #include <seqan3/alignment/all.hpp>
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

#include <seqan3/search/algorithm/detail/search_trivial.hpp>
#include <seqan3/search/fm_index/concept.hpp>

#include <mutex>
#include <ctime>


#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/algorithm/all.hpp>

#include <seqan3/std/filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <fstream>                                  // for writing/reading files
#include <cereal/archives/binary.hpp>               // for storing/loading indices via cereal
#include <seqan3/range/view/convert.hpp>

using namespace seqan3;

template <typename alphabet_t>
auto generate_sequence(size_t const len = 500,
                              size_t const variance = 0,
                              size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, alphabet_size<alphabet_t> - 1);
    std::uniform_int_distribution<size_t> dis_length(len - variance, len + variance);

    std::vector<alphabet_t> sequence;

    size_t length = dis_length(gen);
    for (size_t l = 0; l < length; ++l)
        sequence.push_back(alphabet_t{}.assign_rank(dis_alpha(gen)));

    return sequence;
}

struct options
{
    size_t const sequence_length;
    bool const has_repeats;
    size_t const number_of_reads;
    size_t const read_length;
    double const prob_insertion;
    double const prob_deletion;
    uint8_t const simulated_errors;
    uint8_t const searched_errors;
    uint8_t const strata;
    double const stddev{0};
    uint32_t repeats{20};
};

template <typename alphabet_t>
void mutate_insertion(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed)
{
    std::mt19937_64 gen{seed};
    std::uniform_int_distribution<uint8_t> dis_alpha{0, alphabet_size<alphabet_t> - 1};
    std::uniform_int_distribution<size_t> random_pos{0, std::ranges::size(seq) - overlap};
    seq.insert(std::ranges::begin(seq) + random_pos(gen), alphabet_t{}.assign_rank(dis_alpha(gen)));
}

template <typename alphabet_t>
void mutate_deletion(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed)
{
    std::mt19937_64 gen{seed};
    std::uniform_int_distribution<size_t> random_pos{0, std::ranges::size(seq) - overlap};
    seq.erase(std::ranges::begin(seq) + random_pos(gen));
}

template <typename alphabet_t>
void mutate_substitution(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed)
{
    std::mt19937_64 gen{seed};
    std::uniform_int_distribution<uint8_t> dis_alpha_short{0, alphabet_size<alphabet_t> - 2};
    std::uniform_int_distribution<size_t> random_pos{0, std::ranges::size(seq) - overlap};
    alphabet_t & cbase = seq[random_pos(gen)];
    uint8_t crank = to_rank(cbase);
    uint8_t rrank = dis_alpha_short(gen);
    if (rrank >= crank)
        ++rrank;
    cbase.assign_rank(rrank);
}

template <typename alphabet_t>
std::vector<std::vector<alphabet_t>> generate_reads(std::vector<alphabet_t> const & ref,
                                                    size_t const number_of_reads,
                                                    size_t const read_length,
                                                    uint8_t const simulated_errors_,
                                                    double const prob_insertion,
                                                    double const prob_deletion,
                                                    double const stddev = 0,
                                                    size_t const seed = 0)
{
    std::vector<std::vector<alphabet_t>> reads;
    std::mt19937_64 gen{seed};
    std::uniform_int_distribution<size_t> seeds{0};
    std::normal_distribution<> dis_error_count{static_cast<double>(simulated_errors_), stddev};
    std::uniform_real_distribution<double> probability{0.0, 1.0};

    for (size_t i = 0; i < number_of_reads; ++i)
    {
        // TODO avoid mutating the same position multiple times
        // simulate concrete error number or use normal distribution
        uint8_t simulated_errors = (stddev == 0) ? simulated_errors_ :
                                                   std::abs(std::round(dis_error_count(gen)));
        std::uniform_int_distribution<size_t> random_pos{0, std::ranges::size(ref) - read_length - simulated_errors};
        size_t rpos = random_pos(gen);
        std::vector<alphabet_t> read_tmp{std::ranges::begin(ref) + rpos,
                                         std::ranges::begin(ref) + rpos + read_length + simulated_errors};

        for (uint8_t j = 0; j < simulated_errors; ++j)
        {
            double prob = probability(gen);
            // Substitution
            if (prob_insertion + prob_deletion < prob)
            {
                mutate_substitution(read_tmp, simulated_errors, seeds(gen));
            }
            // Insertion
            else if (prob_insertion < prob)
            {
                mutate_insertion(read_tmp, simulated_errors, seeds(gen));
            }
            // Deletion
            else
            {
                mutate_deletion(read_tmp, simulated_errors, seeds(gen));
            }
        }
        read_tmp.erase(std::ranges::begin(read_tmp) + read_length, std::ranges::end(read_tmp));
        reads.push_back(read_tmp);
    }
    return reads;
}

template <typename alphabet_t>
std::vector<alphabet_t> generate_repeating_sequence(size_t const template_length = 5000,
                                                    size_t const repeats = 20,
                                                    double const template_fraction = 1,
                                                    size_t const seed = 0)
{
    std::vector<alphabet_t> seq_template = generate_sequence<alphabet_t>(template_length, 0, seed);

    // copy substrings of length len from seq_template mutate and concatenate them
    size_t len = std::round(template_length * template_fraction);
    uint8_t simulated_errors = 5;
    len = (len + simulated_errors  > template_length) ? template_length - simulated_errors : len;

    return generate_reads(seq_template, repeats, len, simulated_errors, 0.15, 0.15) | view::persist | std::view::join;
}

//============================================================================
//  undirectional; trivial_search, collection, dna4, all-mapping
//============================================================================
/*
void unidirectional_search_all_collection(benchmark::State & state, options && o)
{
    size_t set_size = 10;
    std::vector<std::vector<seqan3::dna4>> collection;
    std::vector<std::vector<seqan3::dna4>> reads;
    for (size_t i = 0; i < set_size; ++i)
    {
        collection.push_back(generate_sequence<seqan3::dna4>(o.sequence_length, 0, 0));
        std::vector<std::vector<seqan3::dna4>> seq_reads = generate_reads(collection.back(), o.number_of_reads,
                                                                          o.read_length, o.simulated_errors,
                                                                          o.prob_insertion, o.prob_deletion,
                                                                          o.stddev);
        std::move(std::ranges::begin(seq_reads), std::ranges::end(seq_reads), std::back_inserter(reads));
    }

    fm_index index{collection};
    configuration cfg = search_cfg::max_error{search_cfg::total{o.searched_errors}};

    for (auto _ : state)
        auto results = search(reads, index, cfg);
}

//============================================================================
//  undirectional; trivial_search, single, dna4, all-mapping
//============================================================================*/




int main(int argc, char * argv[])
{


    options o{2'000'0, false, 200, 100, 0.18, 0.18, 0, 2, 2, 1.75};
/*
    std::string tmp_dir = "tmp_dir/";

    std::string text{"Garfield the fat cat without a hat."};

    {
        fm_index index{text};
        std::ofstream os{"index.file", std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }

    {
        fm_index<text_layout::single> index; // we need to tell the index that we work on a single text before loading
        std::ifstream is{"index.file", std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
    }*/

    /*
    if(!std::filesystem::exists("bench.fasta"))
    {
        std::vector<seqan3::dna4> ref_c = //(o.has_repeats) ?generate_repeating_sequence<seqan3::dna4>(2 * o.sequence_length / o.repeats, o.repeats, 0.5, 0) :
                                            generate_sequence<seqan3::dna4>(o.sequence_length, 0, 0);


//         std::vector<seqan3::dna4> ref_c{"ACGTTTTTAACCGGT"_dna4};
        std::cout << "Built Index and save\n";

        sequence_file_output fout{(tmp_dir + "bench.fasta")};

        fout.emplace_back((ref_c | view::convert<dna5>), "randomSeq");

        fm_index index_c{ref_c};

        std::cout << "Index built\n";
        std::ofstream os{(tmp_dir + "bench_index"), std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index_c);
    }

    fm_index<text_layout::single> index;
    std::cout << "Load Index";
    std::ifstream is{(tmp_dir + "index.file"), std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(index);
    fm_index index{"ACGTTTTTAACCGGT"_dna4};

    std::cout << "Load Text";
    std::vector<seqan3::dna4> ref;
    sequence_file_input fin{(tmp_dir + "bench.fasta")};

    for (auto & rec : fin)
    {
//         debug_stream << "ID:  " << get<field::ID>(rec) << '\n';
        ref = get<field::SEQ>(rec) | view::convert<dna4>;
        // a quality field also exists, but is not printed, because we know it's empty for FastA files.
    }
*/

    std::vector<seqan3::dna4> ref = //(o.has_repeats) ?generate_repeating_sequence<seqan3::dna4>(2 * o.sequence_length / o.repeats, o.repeats, 0.5, 0) :
                                            generate_sequence<seqan3::dna4>(o.sequence_length, 0, 0);

    fm_index index{ref};



    std::cout << "Index built\n";
    std::vector<std::vector<seqan3::dna4>> reads = generate_reads(ref, o.number_of_reads * 10, o.read_length,
                                                                  o.simulated_errors, o.prob_insertion,
                                                                  o.prob_deletion, /*o.stddev*/0);

    std::vector<std::vector<seqan3::dna4>> many_reads = generate_reads(ref, o.number_of_reads*30, o.read_length,
                                                                  o.simulated_errors, o.prob_insertion,
                                                                  o.prob_deletion, o.stddev);

    configuration cfg = search_cfg::max_error{search_cfg::total{o.searched_errors}};

    configuration cfg2 = search_cfg::itv_threshold{std::make_pair(2, 13)} | search_cfg::max_error{search_cfg::total{o.searched_errors}};

    configuration cfg10 = search_cfg::itv_threshold{std::make_pair(10, 13)} | search_cfg::max_error{search_cfg::total{o.searched_errors}};

    configuration cfg10_15 = search_cfg::itv_threshold{std::make_pair(10, 13)} | search_cfg::max_error{search_cfg::total{o.searched_errors}};

    configuration cfg10_11 = search_cfg::itv_threshold{std::make_pair(10, 7)} | search_cfg::max_error{search_cfg::total{o.searched_errors}};

    configuration cfg2_13 = search_cfg::itv_threshold{std::make_pair(2, 10)} | search_cfg::max_error{search_cfg::total{o.searched_errors}};

    configuration cfg0 = search_cfg::itv_threshold{std::make_pair(0, 13)} | search_cfg::max_error{search_cfg::total{o.searched_errors}};
//     debug_stream << "Text: " << ref << "\n";
//     debug_stream << reads[0] << "\n\n";




    auto start = std::chrono::high_resolution_clock::now();
    auto results4 = search(reads, index, ref, cfg0);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    debug_stream << "DTime " << elapsed.count() << "\n";

    start = std::chrono::high_resolution_clock::now();
    auto results = search(reads, index, ref, cfg);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "DPTime " << elapsed.count() << "\n";
//         debug_stream  << "results " << results << "\n";
/*
    start = std::chrono::high_resolution_clock::now();
    auto results2 = search(reads, index, ref, cfg2);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "Time " << elapsed.count() << "\n";
*/

    start = std::chrono::high_resolution_clock::now();
    auto results3 = search(reads, index, ref, cfg10);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "Time " << elapsed.count() << "\n";


    start = std::chrono::high_resolution_clock::now();
    auto results32 = search(reads, index, ref, cfg10_15);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "Time " << elapsed.count() << "\n";

    start = std::chrono::high_resolution_clock::now();
    auto results33 = search(reads, index, ref, cfg10_11);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "Time " << elapsed.count() << "\n";

    start = std::chrono::high_resolution_clock::now();
    auto results34 = search(reads, index, ref, cfg2_13);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "2Time " << elapsed.count() << "\n";

/*
    configuration cfgHE0 = search_cfg::itv_threshold{0} | search_cfg::max_error{search_cfg::total{o.searched_errors}, search_cfg::substitution{o.searched_errors}, search_cfg::insertion{1}, search_cfg::deletion{0}};

    start = std::chrono::high_resolution_clock::now();
    auto resultsHE0 = search(many_reads, index, ref, cfgHE0);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "D1 Insertion 3 S Time " << elapsed.count() << "\n";

    configuration cfgHE1 = search_cfg::itv_threshold{1} | search_cfg::max_error{search_cfg::total{o.searched_errors}, search_cfg::substitution{o.searched_errors}, search_cfg::insertion{1}, search_cfg::deletion{0}};

    start = std::chrono::high_resolution_clock::now();
    auto resultsHE1 = search(many_reads, index, ref, cfgHE1);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "1 Insertion 3 S Time " << elapsed.count() << "\n";

    configuration cfgH2E1 = search_cfg::itv_threshold{1} | search_cfg::max_error{search_cfg::total{o.searched_errors}, search_cfg::substitution{o.searched_errors}, search_cfg::insertion{2}, search_cfg::deletion{1}};

    start = std::chrono::high_resolution_clock::now();
    auto resultsH2E1 = search(reads, index, ref, cfgH2E1);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "2I 1D 3S Time " << elapsed.count() << "\n";

    configuration cfgH2E0 = search_cfg::itv_threshold{0} | search_cfg::max_error{search_cfg::total{o.searched_errors}, search_cfg::substitution{o.searched_errors}, search_cfg::insertion{2}, search_cfg::deletion{1}};

    start = std::chrono::high_resolution_clock::now();
    auto resultsH2E0 = search(reads, index, ref, cfgH2E0);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "D 2I1D 3S Time " << elapsed.count() << "\n";

    configuration cfgH0 = search_cfg::itv_threshold{0} | search_cfg::max_error{search_cfg::total{o.searched_errors}, search_cfg::substitution{o.searched_errors}, search_cfg::insertion{0}, search_cfg::deletion{0}};

    start = std::chrono::high_resolution_clock::now();
    auto resultsH0 = search(many_reads, index, ref, cfgH0);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "DTime " << elapsed.count() << "\n";

    configuration cfgH1 = search_cfg::itv_threshold{1} | search_cfg::max_error{search_cfg::total{o.searched_errors}, search_cfg::substitution{o.searched_errors}, search_cfg::insertion{0}, search_cfg::deletion{0}};

    start = std::chrono::high_resolution_clock::now();
    auto resultsH1 = search(many_reads, index, ref, cfgH1);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "TH 1 Time " << elapsed.count() << "\n";

    configuration cfgH5 = search_cfg::itv_threshold{5} | search_cfg::max_error{search_cfg::total{o.searched_errors}, search_cfg::substitution{o.searched_errors}, search_cfg::insertion{0}, search_cfg::deletion{0}};

    start = std::chrono::high_resolution_clock::now();
    auto resultsH5 = search(many_reads, index, ref, cfgH5);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "TH 5 Time " << elapsed.count() << "\n";


    std::cout << "aminoacid\n";

    std::vector<seqan3::aa27> refAA = generate_sequence<seqan3::aa27>(o.sequence_length, 0, 0);

    fm_index indexAA{refAA};
    std::cout << "Index built\n";
    std::vector<std::vector<seqan3::aa27>> readsAA = generate_reads(refAA, o.number_of_reads, o.read_length,
                                                                  o.simulated_errors, o.prob_insertion,
                                                                  o.prob_deletion, o.stddev);

    start = std::chrono::high_resolution_clock::now();
    auto results5 = search(readsAA, indexAA, refAA, cfg0);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "Default Time " << elapsed.count() << "\n";

    start = std::chrono::high_resolution_clock::now();
    auto results6 = search(readsAA, indexAA, refAA, cfg);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    debug_stream << "Time " << elapsed.count() << "\n";*/

    std::cout << "fin\n";


    return 0;
}
