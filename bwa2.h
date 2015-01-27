#include "bwtaln.h"
#include "aligner.h"
//#include "aryana_args.h"

char getNuc(uint64_t place, uint64_t seq_len);
void bwa_aln_core2(aryana_args *args);
void bwa_aln_single(const char *prefix, const char *fn_fa);
