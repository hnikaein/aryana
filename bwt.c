/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "utils.h"
#include "bwt.h"
#include "kvec.h"
#include "bwtinline.h"

void bwt_gen_cnt_table(bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		bwt->cnt_table[i] = x;
	}
}

// bwt->bwt and bwt->occ must be precalculated
void bwt_cal_sa(bwt_t *bwt, int intv)
{
	//--------------------------------Milad
	//fprintf(stderr, "\n in bwt_cal_sa\n");
	//--------------------------------END

	bwtint_t isa, sa, i; // S(isa) = sa

	xassert(bwt->bwt, "bwt_t::bwt is not initialized.");

	if (bwt->sa) free(bwt->sa);
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) / intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	if (bwt->sa == 0) {
		fprintf(stderr, "[%s] Fail to allocate %.3fMB memory. Abort!\n", __func__, bwt->n_sa * sizeof(bwtint_t) / 1024.0/1024.0);
		abort();
	}
	// calculate SA value
	isa = 0; sa = bwt->seq_len;
	for (i = 0; i < bwt->seq_len; ++i) {
		if (isa % intv == 0) bwt->sa[isa/intv] = sa;
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if (isa % intv == 0) bwt->sa[isa/intv] = sa;
	bwt->sa[0] = (bwtint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k)
{
	bwtint_t sa = 0;
	while (k % bwt->sa_intv != 0) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	return sa + bwt->sa[k/bwt->sa_intv] % (bwt->seq_len + 1);
}



//--------------------------------Milad

int bwt_match_exact2(const bwt_t *bwt, int len, const ubyte_t *str)
{
	bwtint_t k, l, ok, ol;
	int i;
	k = 0; l = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) break; // no match
	}
	if (k > l) return 0; // no match
	return l - k + 1;
}


//--------------------------------END


//Aryan
/*
int bwt_match_limit_maxocc(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end, bwtint_t *limit, bwtint_t maxocc,int minmatch,int offset,
	int *nums, bwtint_t ups[], bwtint_t downs[], bwtint_t limits[])
{
	bwtint_t k, l, ok, ol,bk=0,bl=0;
	int i;

	*limit=len;
	k = 0; l = bwt->seq_len;
	bwtint_t lastocc=0;
	int last_index;
	for (i = len - 1; i >= 0; --i) {
		if ((len-1)-i>=minmatch && ((lastocc==0 && l-k+1<=maxocc) || (lastocc!=0 && l-k+1!=lastocc && last_index-i>=offset)))
		{
			last_index=i;
			lastocc=l-k+1;
			downs[*nums]=k;
			ups[*nums]=l;
			limits[*nums]=(bwtint_t)(len-i-1);
			(*nums)++;
		}
		ubyte_t c = str[i];
		if (c > 3)
		{
			(*limit)=(bwtint_t)(len-i-1);
			if (sa_begin) *sa_begin = k;
			if (sa_end)   *sa_end = l;

			return 0; // no match
		}
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l)
		{
			*limit=len-i-1;
			if (sa_begin) *sa_begin = bk;
			if (sa_end)   *sa_end = bl;
			return 0; // no match
		}
		bk=k;
		bl=l;
	}
	if (sa_begin) *sa_begin = k;
	if (sa_end)   *sa_end = l;
	if (k > l) return 0; // no match
	return l - k + 1;
}*/


// Tries to find the exact matches of a read (str) having length (len) 
// Returns the interval of the identified exact matches in bwt table, from (sa_begin) to (sa_end)
// If no exact matches are found, it returns the maximum length of a suffix of (str) with exact match (limit)
int bwt_match_limit(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end, bwtint_t *limit)
{
	bwtint_t k, l, ok, ol,bk=0,bl=0;
	int i;

	*limit=len;
	k = 0; l = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3)
		{
			(*limit)=(bwtint_t)(len-i-1);
			if (sa_begin) *sa_begin = k;
			if (sa_end)   *sa_end = l;

			return 0; // no match
		}
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l)
		{
			*limit=len-i-1;
			if (sa_begin) *sa_begin = bk;
			if (sa_end)   *sa_end = bl;
			return 0; // no match
		}
		bk=k;
		bl=l;
	}
	if (sa_begin) *sa_begin = k;
	if (sa_end)   *sa_end = l;
	if (k > l) return 0; // no match
	return l - k + 1;
}

// Same as bwt_match_limit, but finds the exact matches for the reverse complement of (str)
int bwt_match_limit_rev(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end, bwtint_t *limit)
{
	bwtint_t k, l, ok, ol;
	int i;

	*limit=len;
	k = 0; l = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[len-i-1];
		if (c > 3)
		{
			(*limit)=(bwtint_t)(len-i-1);
			return 0; // no match
		}
		c = 3 - c;
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l)
		{
			*limit=len-i-1;
			break; // no match
		}
	}
	if (k > l) return 0; // no match
	if (sa_begin) *sa_begin = k;
	if (sa_end)   *sa_end = l;
	return l - k + 1;
}

//end

int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end)
{
	bwtint_t k, l, ok, ol;
	int i;
	k = 0; l = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) break; // no match
	}
	if (k > l) return 0; // no match
	if (sa_begin) *sa_begin = k;
	if (sa_end)   *sa_end = l;
	return l - k + 1;
}

int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0)
{
	int i;
	bwtint_t k, l, ok, ol;
	k = *k0; l = *l0;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // there is an N here. no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) return 0; // no match
	}
	*k0 = k; *l0 = l;
	return l - k + 1;
}

/*********************
 * Bidirectional BWT *
 *********************/

void bwt_extend(const bwt_t *bwt, const bwtintv_t *ik, bwtintv_t ok[4], int is_back)
{
	bwtint_t tk[4], tl[4];
	int i;
	bwt_2occ4(bwt, ik->x[!is_back] - 1, ik->x[!is_back] - 1 + ik->x[2], tk, tl);
	for (i = 0; i != 4; ++i) {
		ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i];
		ok[i].x[2] = tl[i] - tk[i];
	}
	ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->x[2] - 1 >= bwt->primary);
	ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2];
	ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2];
	ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2];
}

static void bwt_reverse_intvs(bwtintv_v *p)
{
	if (p->n > 1) {
		int j;
		for (j = 0; j < p->n>>1; ++j) {
			bwtintv_t tmp = p->a[p->n - 1 - j];
			p->a[p->n - 1 - j] = p->a[j];
			p->a[j] = tmp;
		}
	}
}

int bwt_smem1(const bwt_t *bwt, int len, const uint8_t *q, int x, bwtintv_v *mem, bwtintv_v *tmpvec[2])
{
	int i, j, c, ret;
	bwtintv_t ik, ok[4];
	bwtintv_v a[2], *prev, *curr, *swap;

	mem->n = 0;
	if (q[x] > 3) return x + 1;
	kv_init(a[0]); kv_init(a[1]);
	prev = tmpvec[0]? tmpvec[0] : &a[0];
	curr = tmpvec[1]? tmpvec[1] : &a[1];
	bwt_set_intv(bwt, q[x], ik);
	ik.info = x + 1;

	for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
		if (q[i] > 3) break;
		c = 3 - q[i];
		bwt_extend(bwt, &ik, ok, 0);
		if (ok[c].x[2] != ik.x[2]) // change of the interval size
			kv_push(bwtintv_t, *curr, ik);
		if (ok[c].x[2] == 0) break; // cannot be extended
		ik = ok[c]; ik.info = i + 1;
	}
	if (i == len) kv_push(bwtintv_t, *curr, ik); // push the last interval if we reach the end
	bwt_reverse_intvs(curr); // s.t. smaller intervals visited first
	ret = curr->a[0].info; // this will be the returned value
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= -1; --i) { // backward search for MEMs
		if (q[i] > 3) break;
		c = i < 0? 0 : q[i];
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			bwtintv_t *p = &prev->a[j];
			bwt_extend(bwt, p, ok, 1);
			if (ok[c].x[2] == 0 || i == -1) { // keep the hit if reaching the beginning or not extended further
				if (curr->n == 0) { // curr->n to make sure there is no longer matches
					if (mem->n == 0 || i + 1 < mem->a[mem->n-1].info>>32) { // skip contained matches
						ik = *p; ik.info |= (uint64_t)(i + 1)<<32;
						kv_push(bwtintv_t, *mem, ik);
					}
				} // otherwise the match is contained in another longer match
			}
			if (ok[c].x[2] && (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2])) {
				ok[c].info = p->info;
				kv_push(bwtintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}
	bwt_reverse_intvs(mem); // s.t. sorted by the start coordinate

	if (tmpvec[0] == 0) free(a[0].a);
	if (tmpvec[1] == 0) free(a[1].a);
	return ret;
}

int bwt_smem(const bwt_t *bwt, int len, const uint8_t *q, bwtintv_v *mem, bwtintv_v *tmpvec[3])
{
	int x = 0, i;
	bwtintv_v a[3], *tvec[2], *mem1;
	kv_init(a[0]); kv_init(a[1]); kv_init(a[2]); // no memory allocation here
	tvec[0] = tmpvec[0]? tmpvec[0] : &a[0];
	tvec[1] = tmpvec[1]? tmpvec[1] : &a[1];
	mem1    = tmpvec[2]? tmpvec[2] : &a[2];
	mem->n = 0;
	do {
		x = bwt_smem1(bwt, len, q, x, mem1, tvec);
		for (i = 0; i < mem1->n; ++i)
			kv_push(bwtintv_t, *mem, mem1->a[i]);
	} while (x < len);
	if (tmpvec[0] == 0) free(a[0].a);
	if (tmpvec[1] == 0) free(a[1].a);
	if (tmpvec[2] == 0) free(a[2].a);
	return mem->n;
}
