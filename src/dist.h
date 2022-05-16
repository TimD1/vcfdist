#ifndef _DIST_H_
#define _DIST_H_

#include "fasta.h"
#include "vcf.h"

#define PTR_NONE 0
#define PTR_UP   1
#define PTR_LEFT 2
#define PTR_BOTH 3
#define PTR_DIAG 4

int ed_align(const vcfData* vcf, const fastaData* const ref);

#endif