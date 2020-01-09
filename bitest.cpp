#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/std/filesystem>

#include <utility>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/pairwise_combine.hpp>

using namespace seqan3;
int main()
{


    dna4_vector seq1 = "TTACGTACGGACTAGCTACAACATTACGGACTAC"_dna4;
    dna4_vector seq2 = "GGACGACATGACGTACGACTTTACGTACGACTAGC"_dna4;
    dna4_vector seq3 = "GGACGACATGATTTTCGACTTTACGTACGACTAGCCC"_dna4;
    dna4_vector seq4 = "TTACGTACGGA"_dna4;
    dna4_vector seq5 = "GGACGACATGA"_dna4;

    std::vector vec{"ACGTGAACTGACT"_dna4,
                    "ACGAAGACCGAT"_dna4,
                    "ACGTGACTGACT"_dna4,
                    "AGGTACGAGCGACACT"_dna4};

    std::vector<std::tuple<dna4_vector, dna4_vector> > jobs;
    jobs.push_back(std::tie(seq1, seq2));
    jobs.push_back(std::tie(seq1, seq3));
    jobs.push_back(std::tie(seq4, seq5));
    jobs.push_back(std::tie(seq4, seq5));


    auto config  = align_cfg::edit /*| align_cfg::aligned_ends{free_ends_first}*/ | align_cfg::result{seqan3::with_alignment};
//     auto results = align_pairwise(std::tie(text_view, seq), config);

    for (auto const & res : align_pairwise(/*views::pairwise_combine(vec)*/jobs, config))
    {
        debug_stream << "Score: " << res.score() << '\n';
        debug_stream << "Begin: " << res.front_coordinate() << '\n';
        debug_stream << "End: " << res.back_coordinate() << '\n';
        debug_stream << "Alignment: \n" << res.alignment() << '\n';
    }

}
