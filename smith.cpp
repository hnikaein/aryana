#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include "aryana_args.h"
extern "C" char getNuc(uint64_t place, uint64_t *reference, uint64_t seq_len);
#include "bwa2.h"
#include "smith.h"

const int mgi = 3; // Maximum length of gap interval to check. Be careful with increasing this.
const char delC = 'a';
const char insC = 'n';
const char matC = 'z';
const char misC = 'Z';

extern "C" int
smith_waterman(aryana_args *options, uint64_t match_start, uint64_t match_end, uint64_t index_start, uint64_t index_end,
               char *cigar, int head, const ubyte_t *read, int len,
               int *mismatch_num, uint64_t seq_len, int **d, char **arr, char *tmp_cigar, uint64_t *reference,
               ignore_mismatch_t ignore, ubyte_t *qual);



int max(int q, int p) {
    return q < p ? p : q;
}

static inline int gap(int len, int go, int ge) {
    return (len <= 0) ? 0 : (len - 1) * ge + go;
}

static inline int insertion(char c) {
    return (c >= insC && c <= insC + mgi - 1);
}

static inline int deletion(char c) {
    return (c >= delC && c <= delC + mgi - 1);
}

static inline int mismatch(ubyte_t qual, int mp) {
    double factor = (74 - (int) qual) / 40;
    int qual_cost = (mp - 1) * factor;
    return mp - qual_cost;
}

int
smith_waterman(aryana_args *options, uint64_t match_start, uint64_t match_end, uint64_t index_start, uint64_t index_end,
               char *cigar, int head, const ubyte_t *read, int len, int *mismatch_num, uint64_t seq_len, int **d,
               char **arr, char *tmp_cigar, uint64_t *reference, ignore_mismatch_t ignore, ubyte_t *qual) {
    // d: the array that stores dynamic programming penalty table
    // arr: the array that stores the best strategy for each cell of dynamic programming table
    // match_start, match_end: positions of the alignment region in the read
    // index_start, index_end: positions of the alignment region in the reference
    // head: the number of characters already filled in cigar, e.g. the position to add new string to cigar
    // return value: the new number of characters filled in cigar (e.g. the new position to add further strings to cigar)
    int mp = options->mismatch_penalty, go = options->gap_open_penalty, ge = options->gap_ext_penalty, ms = options->match_score;
    if ((match_end - match_start == 0) && (index_end - index_start == 0)) return head;    // When both strings are empty
    if (match_end - match_start == 0)
        return head + snprintf(cigar + head, 1000, "%dD", (int) (index_end - index_start)); // When one of them is empty
    if (index_end - index_start == 0)
        return head + snprintf(cigar + head, 1000, "%dI", (int) (match_end - match_start)); // When one of them is empty
    int off = max(2 * (abs((signed) (index_end - index_start) - (signed) (match_end - match_start))), 10);
    //assert(off <= 98);
    if (off > 98) {
        //fprintf(stderr,"off too long\n");
        head += snprintf(cigar + head, 10, "%d",
                         abs((signed) (index_end - index_start) - (signed) (match_end - match_start)));
        if (match_end - match_start > index_end - index_start) {
            head += snprintf(cigar + head, 10, "I");
            head += snprintf(cigar + head, 10, "%"PRIu64"M", index_end - index_start);
        } else {
            head += snprintf(cigar + head, 10, "D");
            head += snprintf(cigar + head, 10, "%"PRIu64"M", match_end - match_start);
        }
        return head;
    }
    off *= 2;
    if (off > 98)
        off = 100;


    int i = 0, j = 0;
    for (i = off / 2 + 1; i < off; i++) {
        if (match_start == 0 || match_start == len) // zero deletion penalty for the beginning of the read
            d[0][i] = 0;
        else
            d[0][i] = gap(i - off / 2, go, ge);
        arr[0][i] = delC;
    }
    d[0][off / 2] = 0;
    int cur_off = 0, best_pen = INT_MAX;
    for (i = 1; i <= match_end - match_start; i++) // The position in read
        for (j = 0; j < off; j++)                // The real position - read position + off/2
        {
            int real_off = j - off / 2;
            if (i < -real_off) continue;
            int ref_i = i + real_off;
            if (ref_i == 0) {
                d[i][j] = gap(i, go, ge);
                arr[i][j] = insC;
                continue;
            }

            d[i][j] = d[i - 1][j];
            arr[i][j] = matC; // match

            if (index_start + ref_i - 1 >= seq_len) {
                d[i][j] += (insertion(arr[i - 1][j]) ? ge : go);
                arr[i][j] = insC;
            } else {
                char gc = getNuc(index_start + ref_i - 1, reference, seq_len), rc = read[match_start + i - 1];
                if (gc != rc) {
                    if (ignore == ignore_none || (ignore == ignore_CT && (gc != 1 || rc != 3)) ||
                        (ignore == ignore_GA && (gc != 2 || rc != 0))) {
                        d[i][j] += mismatch(qual[i], mp);
                        arr[i][j] = misC; // mismatch
                    }
                } else
                    d[i][j] -= ms;
            }
            int k;
            for (k = 1; k <=
                        mgi; k++) { // This is not an optimized dynamic programming here, it is a simple heuristic just to trigger the gaps open
                // For an optimal implementation we need to have two "d" and "d'" tables, one for "no-gap" and the other for "gap already open" cases.
                // The alternative way would be to increase the upperbound 3 on k, which results in quadratic time.
                if (k <= j && i >= off / 2 - j + k) {
                    int pen = d[i][j - k] + (deletion(arr[i][j - k]) ? (k * ge) : (go + (k - 1) * ge));
                    if (d[i][j] > pen) {
                        d[i][j] = pen;
                        arr[i][j] = delC + k - 1;
                    }
                }

                if (k <= i && j + k < off) {
                    int pen = d[i - k][j + k] + (insertion(arr[i - k][j + k]) ? (k * ge) : ((k - 1) * ge + go));
                    if (d[i][j] > pen) {
                        d[i][j] = pen;
                        arr[i][j] = insC + k - 1;
                    }
                }
            }
            if (i == match_end - match_start &&
                d[i][j] < best_pen) { // To find the best "local" alignment, by omitting possible deletions from the end
                best_pen = d[i][j];
                cur_off = j;
            }
        }
    int row, col;
    if (options->debug >= 4) {
        fprintf(stderr, "    ");
        for (col = 0; col < off; col++) fprintf(stderr, "%4d", col - off / 2);
        for (row = 0; row <= (match_end - match_start); row++) {
            fprintf(stderr, "\n%4d", row);
            for (col = 0; col < off; col++) fprintf(stderr, "%4d", d[row][col]);
        }
        fprintf(stderr, "\n\n");
        fprintf(stderr, "    ");
        for (col = 0; col < off; col++) fprintf(stderr, "%4d", col - off / 2);
        for (row = 0; row <= (match_end - match_start); row++) {
            fprintf(stderr, "\n%4d", row);
            for (col = 0; col < off; col++) {
                if (!arr[row][col]) arr[row][col] = ' ';
                fprintf(stderr, "%4c", arr[row][col]);
            }
        }
        fprintf(stderr, "\n\n");
    }
    if (match_end != len) cur_off = (index_end - index_start) - (match_end - match_start) + off / 2;
    int trimmed = (index_end - index_start) - (match_end - match_start) + off / 2 - cur_off;
    int cur_i = match_end - match_start;
    best_pen = d[cur_i][cur_off];
    if (options->debug >= 3)
        fprintf(stderr, "smith_waterman best penalty: %d, cur_off: %d, real_off: %d, off: %d, trimmed: %d\n", best_pen,
                cur_off, cur_off - off / 2, off, trimmed);
    char *tail = tmp_cigar;
    char last = matC;
    while (1) {
        int ref_i = cur_i + cur_off - off / 2;
        if (cur_i == 0) {
            memset(tail, 'D', ref_i);
            tail += ref_i;
            break;
        }
        if (ref_i == 0) {
            memset(tail, 'I', cur_i);
            tail += cur_i;
            break;
        }
        char c = arr[cur_i][cur_off];
        if (c == misC) (*mismatch_num)++;

        if (insertion(c)) {
            int num = c - insC + 1;
            cur_off += num;
            cur_i -= num;
            memset(tail, 'I', num);
            tail += num;
        } else if (deletion(c)) {
            int num = c - delC + 1;
            cur_off -= num;
            memset(tail, 'D', num);
            tail += num;
        } else if (c == matC || c == misC) {
            cur_i--;
            *(tail++) = 'M';
        } else {
            fprintf(stderr, "Error: Invalid value in smith_waterman table\n");
            break;
        }
    }
    *tail = 0;
    int ct = 0;
    last = *(--tail);
    for (; tail >= tmp_cigar; tail--) {
        if (*tail == last)
            ct++;
        else {
            head += snprintf(cigar + head, 10, "%d%c", ct, last);
            last = *tail;
            ct = 1;
        }
    }
    head += snprintf(cigar + head, 10, "%d%c", ct, last);
    if (trimmed != 0) head += snprintf(cigar + head, 10, "%dD", trimmed);
    return head;
}

