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
#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/range/view/pairwise_combine.hpp>


using namespace seqan3;
using namespace std::string_literals; // for using the ""s string literal
// using namespace seqan3::search_cfg;



int main(int argc, char * argv[])
{

//     std::vector vec{"ACGTGACTGACT"_dna4,
//                     "ACGAAGACCGAT"_dna4};
    std::vector<dna4> /*const */seq1 = "CGTGAACTGACTTTTTTTTTTTTTT"_dna4;
    std::vector<dna4> /*const */seq2 = "CGTGACTGACT"_dna4;

    auto slice1 = seq1 | view::slice(0, 15);
    auto slice2 = seq2 | view::slice(0, 15);

    debug_stream << "slice1 " << slice1 << "\t" << "size: " << std::ranges::size(slice1) << "\n";
    debug_stream << "slice2 " << slice2 << "\t" << "size: " << std::ranges::size(slice2) << "\n";

// seqan3::free_ends_first
    front_end_first fef{std::false_type()};
    back_end_first bef{std::false_type()};
        // Configure the alignment kernel.
    auto config = align_cfg::edit |
                    align_cfg::aligned_ends{end_gaps{fef, bef}}|
                    align_cfg::result{seqan3::with_alignment};
//     auto filter_v = std::view::filter([](auto && res) { return res.score() >= -6;});
    for (auto const & res : align_pairwise(std::tie(slice1, slice2), config))
    {
        debug_stream << "Score: " << res.score() << '\n';
         debug_stream << "End: " << res.back_coordinate().first << '\n';
         debug_stream << "Alignment: \n" << res.alignment() << '\n';

    }




    std::cout << "Finished\n";
    return 0;
}
