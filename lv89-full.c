#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "lv89.h"

typedef struct {
	int32_t d, k; // diagonal and furthest-reaching index
} wf_diag_t;



/// extend current wavefront, following matches along diagonal
static int32_t wf_extend(int32_t tl, const char *ts, int32_t ql, const char *qs, 
        int32_t n, wf_diag_t *a, int32_t is_global)
{
	int32_t j;
	for (j = 0; j < n; ++j) {
		wf_diag_t *p = &a[j];
		int32_t k = p->k;
		int32_t max_k = (ql - p->d < tl? ql - p->d : tl) - 1;
		const char *ts_ = ts + 1, *qs_ = qs + p->d + 1;
		uint64_t cmp = 0;

        // compare 8 values at once for efficiency
		while (k + 7 < max_k) {
			uint64_t x = *(uint64_t*)(ts_ + k); // warning: unaligned memory access
			uint64_t y = *(uint64_t*)(qs_ + k);
			cmp = x ^ y;
			if (cmp == 0) k += 8;
			else break;
		}
		if (cmp)
			k += __builtin_ctzl(cmp) >> 3; // on x86, this is done via the BSR instruction: https://www.felixcloutier.com/x86/bsr
		else if (k + 7 >= max_k)
			while (k < max_k && *(ts_ + k) == *(qs_ + k)) // use this for generic CPUs. It is slightly faster than the unoptimized version
				++k;

        // if we're done, stop early
		if (k + p->d == ql - 1 && (!is_global || k == tl - 1)) return -1;

        // update furthest-reaching index
		p->k = k;
	}
	return 0;
}



static int32_t wf_next(int32_t n, const wf_diag_t *a, wf_diag_t *b)
{
    /// current wavefront (a) used to calculate next (b), later copied over
	int32_t j;

    // compute start of b (insert, mismatch)
    // ---------------
    // | a[0] | a[1] |
    // ---------------
    // | b[0] | b[1] |
    // ---------------
	b[0].d = a[0].d - 1;
	b[0].k = a[0].k + 1;
	b[1].d = a[0].d;
	b[1].k = (n == 1 || a[0].k > a[1].k? a[0].k : a[1].k) + 1;

    // compute middle of b (insert, delete, mismatch)
    // --------------------------
    // | a[j-1] | a[j] | a[j+1] |
    // --------------------------
    //                 | b[j+1] |
    // --------------------------
	for (j = 1; j < n - 1; ++j) {
		int32_t k = a[j-1].k;
		k = k > a[j+1].k + 1? k : a[j+1].k + 1;
		k = k > a[j].k + 1? k : a[j].k + 1;
		b[j+1].d = a[j].d, b[j+1].k = k;
	}
    
    // compute end of b (delete, mismatch)
	if (n >= 2) {
		b[n].d = a[n-1].d;
		b[n].k = a[n-2].k > a[n-1].k + 1? a[n-2].k : a[n-1].k + 1;
	}
	b[n+1].d = a[n-1].d + 1;
	b[n+1].k = a[n-1].k;

    // two more possible rows
	return n + 2;
}



static int32_t wf_step(int32_t tl, const char *ts, int32_t ql, const char *qs, 
        int32_t n, wf_diag_t *a, uint64_t *mask, int32_t is_global)
{
	int32_t j, l, m, x;
	wf_diag_t *b = a + n + 2; // temporary array

	// wfa_extend
	x = wf_extend(tl, ts, ql, qs, n, a, is_global);
	if (x < 0) return -1; // we're done

	// wfa_next
    // x is index into a of WF cells used so far
    // m is index into b of WF cells in next wave generated so far
    // this loop is necessary since we drop out-of-bound cells later
	for (x = 0, m = 0, l = 1; l <= n; ++l)
		if (l == n || a[l].d != a[l-1].d + 1)
			m += wf_next(l - x, &a[x], &b[m]), x = l;
            // m += l-x+2, x=l

	// drop out-of-bound cells
	for (j = 0, n = 0; j < m; ++j) {

        // diagonals can be negative, force positive
		int32_t y = b[j].d + tl;

        // boolean mask packed into uint64_t
        // if this diag has been masked, skip (previously out-of-bounds)
		if (mask[y>>6] & 1ULL<<(y&0x3f)) continue;

        // only copy in-bounds cells
		if (b[j].d + b[j].k < ql && b[j].k < tl)
			a[n++] = b[j];

        // if not in-bounds, mask this diagonal
		else mask[y>>6] |= 1ULL << (y&0x3f);
	}

    // return wavefront size (some dropped)
	return n;
}



// create buffer for landau-vishkin edit distance calculation
int32_t lv_ed_bufsize(int32_t tl, int32_t ql)
{
    // packed uint64_t boolean mask for valid diagonals
	int32_t mask = (tl + ql + 1 + 63) / 64 * 8;

    // two arrays of int32_t (d,k) for current/next wavefront
    int32_t a = (tl + ql + 1) * 2 * 8;
    return mask + a;
    /* return (tl + ql + 1 + 63) / 64 * 8 + (tl + ql + 1) * 2 * 8; */
}

// mem should be at least (tl+ql)*16 long
int32_t lv_ed(int32_t tl, const char *ts, int32_t ql, const char *qs, 
        int32_t is_global, uint8_t *mem)
{
    // n is the number of cells in current wavefront
	int32_t s = 0, n = 1;

    // mem is concatenated array (mask, a)
	uint64_t *mask = (uint64_t*)mem;
	wf_diag_t *a;

    // zero mem of size tl+ql+1, rounded to next multiple of 64 -> mask
	memset(mask, 0, (tl + ql + 1 + 63) / 64 * sizeof(*mask));

    // a is array of L_{d,e} pairs renamed (d,k)
	a = (wf_diag_t*)(mask + (tl + ql + 1 + 63) / 64);

    // d = diagonal
    // k = offset of furthest-reaching point
	a[0].d = 0, a[0].k = -1;
	while (1) {
		n = wf_step(tl, ts, ql, qs, n, a, mask, is_global);
		if (n < 0) break;
		++s;
	}

    // s is number of steps/edits
	return s;
}
