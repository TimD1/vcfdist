#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "edit_dist.h"

// wavefront diagonal info
typedef struct {
	int32_t d;  // diagonal
	int32_t k:30, p:2; // wavefront distance index, backtrace pointer
} wf_diag_t;

// per-wavefront traceback info
typedef struct {
	int32_t n, d0; // cells in wavefront, first diagonal of wavefront
	uint64_t *a; // TODO: wavefront distance indices??? 
} wf_tb1_t;

// per-alignment traceback info
typedef struct {
	int32_t m, n; // max size of a, current size of a
	wf_tb1_t *a; // array of wavefronts
} wf_tb_t;



// extend traceback data structure, allocating more memory and updating size
static void wf_tb_add(wf_tb_t *tb, int32_t m)
{
    // expand array of alignment traceback info if necessary
	if (tb->n == tb->m) {
		tb->m = tb->m + (tb->m>>1) + 4; // m = m + m/2 + 4
		tb->a = (wf_tb1_t*)realloc(tb->a, tb->m * sizeof(*tb->a));
	}

    // generate zero-ed tb 
	wf_tb1_t *p;
	p = &tb->a[tb->n++];
	p->n = m;
	p->a = (uint64_t*)calloc((m + 31) / 32, sizeof(uint64_t)); // 2-bit pointer
}



static int32_t wf_step(wf_tb_t *tb, int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t n, wf_diag_t *a)
{
	int32_t j;
	wf_tb1_t *q;
	wf_diag_t *b = a + n + 2; // temporary array

	// wfa_extend
	for (j = 0; j < n; ++j) {
		wf_diag_t *p = &a[j];
		int32_t k = p->k, max_k;
		const char *ts_, *qs_;
		uint64_t cmp = 0;
		if (k >= tl || k + p->d >= ql) continue;
		max_k = (ql - p->d < tl? ql - p->d : tl) - 1;
		ts_ = ts + 1;
		qs_ = qs + p->d + 1;

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

        // if we're done, stop
		if (k + p->d == ql - 1 && k == tl - 1) return -1;

        // update furthest-reaching index
		p->k = k;
	}

	// wfa_next
    
    // compute start of b (insert, mismatch)
	b[0].d = a[0].d - 1;
	b[0].p = 1;
	b[0].k = a[0].k + 1;
	b[1].d = a[0].d;
	b[1].p =  n == 1 || a[0].k > a[1].k? 0 : 1;
	b[1].k = (n == 1 || a[0].k > a[1].k? a[0].k : a[1].k) + 1;

    // compute middle of b (insert, delete, mismatch)
	for (j = 1; j < n - 1; ++j) {
		int32_t k = a[j-1].k, p = -1;
		p = k > a[j+1].k + 1? p : 1;
		k = k > a[j+1].k + 1? k : a[j+1].k + 1;
		p = k > a[j].k + 1? p : 0;
		k = k > a[j].k + 1? k : a[j].k + 1;
		b[j+1].d = a[j].d, b[j+1].k = k, b[j+1].p = p;
	}

    // compute end of b (delete, mismatch)
	if (n >= 2) {
		b[n].d = a[n-1].d;
		b[n].p = a[n-2].k > a[n-1].k + 1? -1 : 0;
		b[n].k = a[n-2].k > a[n-1].k + 1? a[n-2].k : a[n-1].k + 1;
	}
	b[n+1].d = a[n-1].d + 1;
	b[n+1].p = -1;
	b[n+1].k = a[n-1].k;

	// keep traceback information
	wf_tb_add(tb, n + 2); // add new wavefront to alignment info
	q = &tb->a[tb->n - 1];
	q->d0 = b[0].d; //store first diagonal

    // store backtrace pointers
    // TODO: why +1? any checks for non-zero pointers?
	for (j = 0; j < n + 2; ++j)
		q->a[j>>5] |= (uint64_t)(b[j].p + 1) << (j&0x1f)*2;

	// save results to a, two more rows now
	memcpy(a, b, (n + 2) * sizeof(*a));
	return n + 2;
}

typedef struct {
	int32_t m, n; // max cigar len, current cigar len
	uint32_t *cigar;
} wf_cigar_t;

static void wf_cigar_push(wf_cigar_t *c, int32_t op, int32_t len)
{
    // 4-bit op, 28-bit len = 268M (pretty reasonable)
	if (c->n && op == (c->cigar[c->n-1]&0xf)) { // extend len of most recent op
		c->cigar[c->n-1] += len<<4;
	} else {
		if (c->n == c->m) { // allocate more space in cigar
			c->m = c->m + (c->m>>1) + 4;
			c->cigar = (uint32_t*)realloc(c->cigar, c->m * sizeof(*c->cigar));
		}
        // add op and len to cigar
		c->cigar[c->n++] = len<<4 | op;
	}
}



static uint32_t *wf_traceback(int32_t tl, const char *ts, 
        int32_t ql, const char *qs, wf_tb_t *tb, int32_t *n_cigar)
{
    // initialize cigar and pointers
	wf_cigar_t cigar = {0,0,0};
	int32_t i = ql-1, k = tl-1, s = tb->n - 1; // query, template ptrs, score

	for (;;) {

        // slide back along diag (opposite of wfa_extend)
		int32_t k0 = k, j, pre;
		while (i >= 0 && k >= 0 && qs[i] == ts[k])
			--i, --k;

        // push matches to cigar
		if (k0 - k > 0)	
			wf_cigar_push(&cigar, E, k0 - k);

        // end if we're done
		if (i < 0 || k < 0) break;

        // score should always be positive
		if (s < 0) fprintf(stderr, "i=%d, k=%d, s=%d\n", i, k, tb->n);
		assert(s >= 0);

        // get index of current cell in wavefront
		j = i - k - tb->a[s].d0;
		assert(j < tb->a[s].n); // check in range

        // get predecessor pointer
		pre = (int32_t)(tb->a[s].a[j>>5] >> (j&0x1f)*2 & 0x3) - 1;

        // update CIGAR
		if (pre == 0) {
			wf_cigar_push(&cigar, X, 1);
			--i, --k;
		} else if (pre < 0) {
			wf_cigar_push(&cigar, I, 1);
			--i;
		} else {
			wf_cigar_push(&cigar, D, 1);
			--k;
		}

        // go to previous wavefront
		--s;
	}
    // push remaining INSs
	if (i > 0) wf_cigar_push(&cigar, I, i);
    // push remaining DELs
	else if (k > 0) wf_cigar_push(&cigar, D, k);

    // reverse cigar
	for (i = 0; i < cigar.n>>1; ++i) {
		uint32_t t = cigar.cigar[i];
		cigar.cigar[i] = cigar.cigar[cigar.n - i - 1];
		cigar.cigar[cigar.n - i - 1] = t;
	}

    // return cigar and length
	*n_cigar = cigar.n;
	return cigar.cigar;
}



// mem should be at least (tl+ql)*16 long
uint32_t *edit_dist_cigar(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t *score, int32_t *n_cigar)
{
	int32_t s = 0, n = 1, i;
	wf_diag_t *a;
	uint32_t *cigar;
	wf_tb_t tb = {0,0,0};
	a = (wf_diag_t*)malloc(2 * (tl + ql + 1) * sizeof(*a));
	a[0].d = 0, a[0].k = -1;
	while (1) {
		n = wf_step(&tb, tl, ts, ql, qs, n, a);
		if (n < 0) break;
		++s;
	}
	free(a);
	cigar = wf_traceback(tl, ts, ql, qs, &tb, n_cigar);
	for (i = 0; i < tb.n; ++i) free(tb.a[i].a);
	free(tb.a);
	*score = s;
	return cigar;
}
