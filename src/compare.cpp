#include <stdio.h>
#include "vcf.h"
#include "vcfutils.h"
 
int main(int argc, char **argv) {
    if (argc != 2) {
        printf("supply VCF");
        return 1;
    }

    
    htsFile * inf = bcf_open(argv[1], "r");
    if (inf == NULL) {
        return EXIT_FAILURE;
    }
        
    // read header
    bcf_hdr_t *hdr = bcf_hdr_read(inf);
 
    // struc for storing each record
    bcf1_t *rec = bcf_init();
    if (rec == NULL) {
        return EXIT_FAILURE;
    }
 
    while (bcf_read(inf, hdr, rec) == 0) {
        if (bcf_is_snp(rec)) {
            printf("Number of alleleles: %lu\n", (unsigned long)rec->n_allele);
            printf("Chr name is %s\n", bcf_hdr_id2name(hdr, rec->rid));
            printf("Position for this SNP is %li\n", rec->pos);
            printf("REF:%s ALT:%s\n", rec->d.allele[0], rec->d.allele[1]);
            // iterating over alleles
            for (int i=0; i<rec->n_allele; ++i)
                printf("%s\n", rec->d.allele[i]);
        } else {
            continue;
        }
    }
 
    bcf_hdr_destroy(hdr);
    bcf_close(inf);
    bcf_destroy(rec);
    return EXIT_SUCCESS;
}
