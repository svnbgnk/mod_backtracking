/*
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



#include <seqan3/range/concept.hpp>
#include <seqan3/search/algorithm/detail/search_common.hpp>

#include <seqan3/core/metafunction/pre.hpp>
#include <seqan3/search/algorithm/detail/search_trivial.hpp>
#include <seqan3/search/fm_index/concept.hpp>

#include <mutex>
#include <ctime>
*/

#include <seqan3/io/alignment_file/all.hpp>
#include <range/v3/view/drop_exactly.hpp>
#include <fstream>
#include <numeric>
#include <range/v3/view/split.hpp>
#include <seqan3/argument_parser/all.hpp> // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>   // our custom output stream
#include <seqan3/std/charconv>            // includes std::from_chars
#include <seqan3/std/filesystem>          // use std::filesystem::path
#include <bitset>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/std/ranges> // std::ranges::copy

#include <utility>

#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <map>
#include <seqan3/range/views/all.hpp>
#include <range/v3/range/conversion.hpp>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
#include <seqan3/alphabet/cigar/cigar_op.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>

using namespace seqan3;
using namespace std::string_literals; // for using the ""s string literal


int main(int argc, char * argv[])
{
    using seqan3::operator""_tag;
    using seqan3::get;
    using seqan3::operator""_cigar_op;
    std::filesystem::path bamPath;
    std::filesystem::path referencePath;
    int stop = 0;
    bool verbose{};



    argument_parser parser{"CheckAlignments", argc, argv};        // initialise myparser

    parser.info.short_description = "Check Alignments.";
    parser.info.version = "1.0.0";
    parser.add_option(bamPath, 'b', "bam", "Alignmenfile.");
    parser.add_option(referencePath, 'r', "ref", "Reference file.");
    parser.add_option(stop, 's', "stop", "Reference file.");
    parser.add_flag(verbose, 'v', "verbose", "Print Alignments read info.");

    // ... add information, options, flags and positional options
    try
    {
         parser.parse();                                          // trigger command line parsing
    }
    catch (parser_invalid_argument const & ext)                     // catch user errors
    {
        debug_stream << "Parser error  " << ext.what() << "\n"; // customise your error message
        return -1;
    }


    debug_stream << "BAM:" << bamPath << "\n";
    debug_stream << "Ref:" << referencePath << "\n";

/*
    using seqan3::get;
    seqan3::sequence_file_input fin{referencePath};
    using record_type = decltype(fin)::record_type;
    std::vector<record_type> reference{};
    std::ranges::copy(fin, std::ranges::back_inserter(reference));
*/

    std::vector<std::string> refIds;
    std::vector<std::vector<dna5>> redSeq;

    sequence_file_input reference_in{referencePath};
//     typedef std::map<std::string, unsigned>  IDMAP;
//     IDMAP idMap;

    int kk = 0;
    for (auto & [seq, id, qual] : reference_in)
    {
        refIds.push_back(std::move(id));
        redSeq.push_back(std::move(seq));
//         idMap[refIds.back()] = kk;
//         debug_stream << "!!!!!!!!!!!!!!!!!!!" << refIds.back() << "\n";
//         ++kk;
    }


    alignment_file_input finalign{bamPath, fields<field::ID, field::SEQ, field::FLAG, field::REF_ID, field::REF_OFFSET, field::CIGAR, field::TAGS/*, field::ALIGNMENT*/>{}};


/*
    seqan3::sam_tag_dictionary dict{};         // initialise empty dictionary
    dict.get<"NM"_tag>() = 3;          // set SAM tag 'NM' to 3 (integer type)
    dict.get<"CO"_tag>() = "comment";  // set SAM tag 'CO' to "comment" (string type)
    auto nm = dict.get<"NM"_tag>();    // get SAM tag 'NM' (note: type is int32_t)
    auto co = dict.get<"CO"_tag>();    // get SAM tag 'CO' (note: type is std::string)
    seqan3::debug_stream << nm << '\n';      // will print '3'
    seqan3::debug_stream << co << '\n';      // will print "comment"*/

    int nalignmentcounter = 0;
    int wrongalignmentcounter = 0;
    int k = 0;
    for (auto & record : finalign)
    {
        ++k;
        if(stop <= k && stop > 0)
            break;
        auto & id = get<seqan3::field::ID>(record);
        int flags = get<seqan3::field::FLAG>(record);
        auto flag =  std::bitset<16>(flags);
        auto & seq = get<seqan3::field::SEQ>(record);
        auto & cigar = get<seqan3::field::CIGAR>(record);
//         auto & alignment = get<seqan3::field::ALIGNMENT>(record);
        auto & refid = get<seqan3::field::REF_ID>(record);
//         auto & offset = get<seqan3::field::OFFSET>(record);
        auto & offsetref = get<seqan3::field::REF_OFFSET>(record);

        seqan3::sam_tag_dictionary dict = get<seqan3::field::TAGS>(record);
        auto nm = dict.get<"NM"_tag>();

        int endPos = 0;
        int soft = 0;
        for(auto & cha : cigar){
            if(get<1>(cha) == 'M'_cigar_op)
                endPos += get<0>(cha);
            if(get<1>(cha) == 'I'_cigar_op)
                ;
            if(get<1>(cha) == 'D'_cigar_op)
                endPos += get<0>(cha);
            if(get<1>(cha) == 'S'_cigar_op)
                soft += get<0>(cha);
        }


        if(verbose){
    //         auto & alignment = get<seqan3::field::ALIGNMENT>(record);
            debug_stream << id << "\n";
            debug_stream << flag.to_string() << "\n";
    //         debug_stream << seq << "\n";
            debug_stream << refid << "\n";
            debug_stream << offsetref << "\t" << endPos << "\n\n";
    //         debug_stream << cigar << "\n\n";
        }


        int refidd = *refid;
        int offsett = *offsetref;


/*
        std::cout << "refidd  " << refidd << "\n";
        std::cout << "off  " << offsett << "\n";*/
        if(redSeq[refidd].size() <= offsett + endPos)
        {
            std::cerr << "WARNING: Seq to short. Skipping\n";
            continue;
        }
//         std::span text_view{std::data(redSeq[refidd]) + offsetref, seq.size() + 1};
        std::span text_view{std::data(redSeq[refidd]) + offsett, static_cast<unsigned>(endPos)};
//         std::span text_view{std::data(redSeq[0]) + offsetref, seq.size() + 1};

        if(verbose){
            debug_stream << text_view << "\n";
            debug_stream << seq << "\n";
        }

        bool noN = false;
        for(int i = 0; i < text_view.size(); ++i){
            if(text_view[i] == 'N'_dna5){
                noN = true;
                ++nalignmentcounter;
                break;
                exit(0);
            }
        }



/*
        //check for rc
        if(flag[4]){
            std::cout << "RC\n";
            seq = seq | seqan3::views::complement | std::views::reverse; //std::ranges::views::reverse;// std::views::reverse;
        }
*/




        // Configure the alignment kernel.
        auto config  = align_cfg::edit /*| align_cfg::aligned_ends{free_ends_first}*/ | align_cfg::result{seqan3::with_alignment};
//         auto config = align_cfg::mode{global_alignment}| align_cfg::edit | align_cfg::result{with_alignment} |align_cfg::aligned_ends{free_ends_first};


        // Invoke the pairwise alignment which returns a lazy range over alignment results.
        auto results = align_pairwise(std::tie(text_view, seq), config);
        auto & res = *results.begin();
        if(verbose){
            std::cout << "Reported errors: "<< nm << " + " << soft << "\n";
            debug_stream << "Score: " << res.score() << "\n";
            debug_stream << "Cigar: " << cigar << "\n";
            debug_stream << "Alignment: " << res.alignment() << "\n\n";
        }



        if(nm + soft != -1*(res.score()))
        {
            ++wrongalignmentcounter;
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
//             exit(0);
        }

//         debug_stream << "Alignment from cigar \n" << get<0>(alignment) << "\n" << get<1>(alignment) << '\n'; // Now you can print the whole alignment!
//         debug_stream << cigar << "\n";


    }



    std::cout << "Alignments containing Ns: " << nalignmentcounter << "\n";
    std::cout << "Alignments with wrong error count: " << wrongalignmentcounter << "\n";


    std::cout << "fin\n";
    return 0;
}
