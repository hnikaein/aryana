#include "aryana_args.h"
#include "bwt.h"
//#include "kbwt.h"
//#include "ttest.h"

void aligner(bwt_t *const bwt, int len, ubyte_t *seq, bwtint_t level, hash_element * table, int *best, int best_size, int *best_found, aryana_args *args);

void create_cigar(hash_element *best, char *cigar, int len, const ubyte_t *seq, uint64_t seq_len,int **d, char **arr, char * tmp_cigar, penalty_t * penalty, uint64_t * reference);

void showerr(int len, const ubyte_t *seq);

void show(int len, const ubyte_t *seq);
