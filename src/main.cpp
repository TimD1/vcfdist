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
std::vector< std::string > reg_strs = {"OUTSIDE", "INSIDE", "BORDER"};
 
int main(int argc, char **argv) {

    // parse and store command-line args
    g.parse_args(argc, argv);
    fastaData ref(g.ref_fasta_fp);

    bedData bed("./test.bed");
    printf("%s\n\n", reg_strs[bed.contains("chr10", 110, 150)].data());  // OUTSIDE
    printf("%s\n\n", reg_strs[bed.contains("chr20", 0, 50)].data());     // OUTSIDE
    printf("%s\n\n", reg_strs[bed.contains("chr20", 50, 150)].data());   // BORDER
    printf("%s\n\n", reg_strs[bed.contains("chr20", 50, 250)].data());   // BORDER
    printf("%s\n\n", reg_strs[bed.contains("chr20", 140, 150)].data());  // INSIDE
    printf("%s\n\n", reg_strs[bed.contains("chr20", 160, 250)].data());  // BORDER
    printf("%s\n\n", reg_strs[bed.contains("chr20", 160, 600)].data());  // BORDER
    printf("%s\n\n", reg_strs[bed.contains("chr20", 10, 600)].data());   // BORDER
    printf("%s\n\n", reg_strs[bed.contains("chr20", 200, 203)].data());  // OUTSIDE
    printf("%s\n\n", reg_strs[bed.contains("chr20", 204, 400)].data());  // BORDER
    printf("%s\n\n", reg_strs[bed.contains("chr20", 350, 550)].data());  // BORDER
    printf("%s\n\n", reg_strs[bed.contains("chr20", 450, 550)].data());  // BORDER
    printf("%s\n\n", reg_strs[bed.contains("chr20", 550, 550)].data());  // OUTSIDE

    printf("%s\n\n", reg_strs[bed.contains("chr20", 200, 205)].data());  // OUTSIDE
    printf("%s\n\n", reg_strs[bed.contains("chr20", 200, 206)].data());  // BORDER
    printf("%s\n\n", reg_strs[bed.contains("chr20", 199, 200)].data());  // INSIDE
    printf("%s\n\n", reg_strs[bed.contains("chr20", 200, 201)].data());  // OUTSIDE
    printf("%s\n\n", reg_strs[bed.contains("chr20", 199, 201)].data());  // BORDER
    printf("%s\n\n", reg_strs[bed.contains("chr20", 204, 205)].data());  // OUTSIDE
    printf("%s\n\n", reg_strs[bed.contains("chr20", 205, 206)].data());  // INSIDE
    printf("%s\n\n", reg_strs[bed.contains("chr20", 204, 206)].data());  // BORDER

    /* vcfData calls(g.calls_vcf_fp); */
    /* edit_dist_realign(&calls, &ref); */

    /* vcfData truth(g.truth_vcf_fp); */
    /* edit_dist_realign(&truth, &ref); */

    /* /1* edit_dist(&calls, &truth, &ref); *1/ */

    return EXIT_SUCCESS;
}
