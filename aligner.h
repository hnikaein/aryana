#include "aryana_args.h"
#include "bwt.h"
//#include "kbwt.h"
//#include "ttest.h"

void aligner(bwt_t *const bwt, int len, ubyte_t *seq, bwtint_t level, hash_element *table, int *best, int best_size,
             int *best_found, aryana_args *args, uint64_t *reference);

void
create_cigar(aryana_args *options, hash_element *best, char *cigar, int len, const ubyte_t *seq, const ubyte_t *qual,
             uint64_t seq_len, int **d, char **arr, char *tmp_cigar, penalty_t *penalty, uint64_t *reference,
             ignore_mismatch_t ignore);

void showerr(int len, const ubyte_t *seq);

void show(int len, const ubyte_t *seq);
