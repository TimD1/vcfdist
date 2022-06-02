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

#define CAL1_HAP1 0
#define CAL1_HAP2 1
#define CAL2_HAP1 2
#define CAL2_HAP2 3

vcfData edit_dist_realign(const vcfData* vcf, const fastaData* const ref);
clusterData edit_dist(vcfData* calls, vcfData* truth, fastaData* ref);

#endif
