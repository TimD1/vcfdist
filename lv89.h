#ifndef LV89_ALN_H
#define LV89_ALN_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Find the edit distance between two sequences
 *
 * @param tl         target sequence length
 * @param ts         target sequence
 * @param ql         query sequence length
 * @param qs         query sequence
 * @param is_global  false to ignore gaps at the end of the target
 * @param mem        temporary memory of lv_ed_bufsize(tl, ql) bytes
 *
 * @return edit distance
 */
int32_t lv_ed(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t is_global, uint8_t *mem);

int32_t lv_ed_semi(int32_t tl, const char *ts, int32_t ql, const char *qs, uint8_t *mem);

uint32_t *lv_ed_semi_cigar(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t *score, int32_t *n_cigar);

int32_t lv_ed_bufsize(int32_t tl, int32_t ql);

#ifdef __cplusplus
}
#endif

#endif
