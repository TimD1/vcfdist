#include "vcf.h"
#include "print.h"
#include "globals.h"
#include "fasta.h"
#include "bed.h"
#include "dist.h"

Globals g;
 
int main(int argc, char **argv) {

    // parse and store command-line args
    g.parse_args(argc, argv);
    fastaData ref(g.ref_fasta_fp);

    vcfData calls(g.calls_vcf_fp);
    ed_align(&calls, &ref);

    /* vcfData truth(g.truth_vcf_fp); */
    /* ed_align(&truth, &ref); */

    return EXIT_SUCCESS;
}
