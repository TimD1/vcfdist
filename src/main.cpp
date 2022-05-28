#include "vcf.h"
#include "print.h"
#include "globals.h"
#include "fasta.h"
#include "bed.h"
#include "dist.h"

Globals g;
std::vector< std::string > type_strs = { "REF", "SUB", "INS", "DEL", "GRP"};
std::vector< std::string > gt_strs = {
".|.", "0|0", "0|1", "0|2", "1|0", "1|1", "1|2", "2|0", "2|1", "2|2", "?|?" };
std::vector< std::string > region_strs = {"OUTSIDE  ", "INSIDE   ", "BORDER   ", "OTHER CTG"};
std::vector< std::string > aln_strs = {"CAL1-HAP1", "CAL1-HAP2", "CAL2-HAP1", "CAL2-HAP2"};
 
int main(int argc, char **argv) {

    // parse and store command-line args
    g.parse_args(argc, argv);
    fastaData ref(g.ref_fasta_fp);

    vcfData calls(g.calls_vcf_fp);
    /* edit_dist_realign(&calls, &ref); */

    vcfData truth(g.truth_vcf_fp);
    /* edit_dist_realign(&truth, &ref); */

    edit_dist(&calls, &truth, &ref);

    return EXIT_SUCCESS;
}
