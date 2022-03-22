#include <assert.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "edit_dist.h"
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *ks1, *ks2;
	ketopt_t o = KETOPT_INIT;
	int c, s, report_cigar = 0;
	int32_t n_cigar;
	uint32_t *cigar = 0;

	while ((c = ketopt(&o, argc, argv, 1, "c", 0)) >= 0) {
		if (c == 'c') report_cigar = 1;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: vcfdist [options] <in1.fa> <in2.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -c    report CIGAR\n");
		return 1;
	}

	fp1 = gzopen(argv[o.ind+0], "r");
	fp2 = gzopen(argv[o.ind+1], "r");
	assert(fp1 && fp2);
	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);
	kseq_read(ks1);
	kseq_read(ks2);

    if (report_cigar) {
        cigar = edit_dist_cigar(ks1->seq.l, ks1->seq.s, ks2->seq.l, ks2->seq.s, &s, &n_cigar);
    } else {
        uint8_t *mem = (uint8_t*)malloc(edit_dist_bufsize(ks1->seq.l, ks2->seq.l));
        s = edit_dist(ks1->seq.l, ks1->seq.s, ks2->seq.l, ks2->seq.s, mem);
        free(mem);
    }

	if (!report_cigar) {
		printf("%s\t%s\t%d\n", ks1->name.s, ks2->name.s, s);
	} else {
		int32_t i;
		printf("%s\t%s\t%d\t", ks1->name.s, ks2->name.s, s);
		for (i = 0; i < n_cigar; ++i)
			printf("%d%c", cigar[i]>>4, "MIDNSHP=X"[cigar[i]&0xf]);
		putchar('\n');
	}
	free(cigar);

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}
