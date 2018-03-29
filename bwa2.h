

typedef struct {
    int mismatch_num, gap_open_num, gap_ext_num, penalty;
    double mapq;
} penalty_t;



char getNuc(uint64_t place, uint64_t * reference, uint64_t seq_len);
void bwa_aln_core2(aryana_args *args);
void bwa_aln_single(const char *prefix, const char *fn_fa);
