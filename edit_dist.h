#ifndef VCF_DIST_H
#define VCF_DIST_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// CIGAR definitions
#define M 0
#define I 1
#define D 2
#define N 3
#define S 4
#define H 5
#define P 6
#define E 7
#define X 8

/**
 * Find the edit distance between two sequences
 *
 * @param tl         target sequence length
 * @param ts         target sequence
 * @param ql         query sequence length
 * @param qs         query sequence
 * @param mem        temporary memory of lv_ed_bufsize(tl, ql) bytes
 *
 * @return edit distance
 */
int32_t edit_dist(int32_t tl, const char *ts, int32_t ql, const char *qs, uint8_t *mem);

uint32_t *edit_dist_cigar(int32_t tl, const char *ts, 
        int32_t ql, const char *qs, int32_t *score, int32_t *n_cigar);

int32_t edit_dist_bufsize(int32_t tl, int32_t ql);

#ifdef __cplusplus
}
#endif

#endif
