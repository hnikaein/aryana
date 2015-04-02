int bwt_match_exact2(const bwt_t *bwt, int len, const ubyte_t *str);
int bwt_match_limit(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end, bwtint_t *limit);
int bwt_match_limit_rev(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end, bwtint_t *limit);
