#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "utils.h"
#include "bwt.h"
#include "kvec.h"
//#include "bwtinline.h"


int bwt_match_exact2(const bwt_t *bwt, int len, const ubyte_t *str) {
    bwtint_t k, l, ok, ol;
    int i;
    k = 0;
    l = bwt->seq_len;
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


// Tries to find the exact matches of a read (str) having length (len)
// Returns the interval of the identified exact matches in bwt table, from (sa_begin) to (sa_end)
// If no exact matches are found, it returns the maximum length of a suffix of (str) with exact match (limit)
int
bwt_match_limit(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end, bwtint_t *limit) {
    bwtint_t k, l, ok, ol, bk = 0, bl = 0;
    int i;

    *limit = len;
    k = 0;
    l = bwt->seq_len;
    for (i = len - 1; i >= 0; --i) {
        ubyte_t c = str[i];
        if (c > 3) {
            (*limit) = (bwtint_t) (len - i - 1);
            if (sa_begin) *sa_begin = k;
            if (sa_end) *sa_end = l;

            return 0; // no match
        }
        bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
        k = bwt->L2[c] + ok + 1;
        l = bwt->L2[c] + ol;
        if (k > l) {
            *limit = len - i - 1;
            if (sa_begin) *sa_begin = bk;
            if (sa_end) *sa_end = bl;
            return 0; // no match
        }
        bk = k;
        bl = l;
    }
    if (sa_begin) *sa_begin = k;
    if (sa_end) *sa_end = l;
    if (k > l) return 0; // no match
    return l - k + 1;
}

// Same as bwt_match_limit, but finds the exact matches for the reverse complement of (str)
int bwt_match_limit_rev(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end,
                        bwtint_t *limit) {
    bwtint_t k, l, ok, ol;
    int i;

    *limit = len;
    k = 0;
    l = bwt->seq_len;
    for (i = len - 1; i >= 0; --i) {
        ubyte_t c = str[len - i - 1];
        if (c > 3) {
            (*limit) = (bwtint_t) (len - i - 1);
            return 0; // no match
        }
        c = 3 - c;
        bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
        k = bwt->L2[c] + ok + 1;
        l = bwt->L2[c] + ol;
        if (k > l) {
            *limit = len - i - 1;
            break; // no match
        }
    }
    if (k > l) return 0; // no match
    if (sa_begin) *sa_begin = k;
    if (sa_end) *sa_end = l;
    return l - k + 1;
}
