#ifndef _DIST_H_
#define _DIST_H_

#include "fasta.h"
#include "vcf.h"
#include "cluster.h"

#define PTR_NONE 0
#define PTR_UP   1
#define PTR_LEFT 2
#define PTR_BOTH 3
#define PTR_DIAG 4
#define PTR_SWAP 8
#define LEFT_PATH 16
#define RIGHT_PATH 32
#define PTR_SUB 64
#define PTR_DONE 128
#define PTR_NEXT 256

#define CAL1_HAP1 0
#define CAL1_HAP2 1
#define CAL2_HAP1 2
#define CAL2_HAP2 3

#define ERRTYPES 6
#define ERRTYPE_TP 0 // true positive
#define ERRTYPE_FP 1 // false positive
#define ERRTYPE_FN 2 // false negative
#define ERRTYPE_PP 3 // partial positive (reduces ED, but not TP)
#define ERRTYPE_PE 4 // phase error (0|1 -> 1|0)
#define ERRTYPE_GE 5 // genotype error (0|1 -> 1|1)

vcfData edit_dist_realign(const vcfData* vcf, const fastaData* const ref);
clusterData edit_dist(vcfData* calls, vcfData* truth, fastaData* ref);

#endif
