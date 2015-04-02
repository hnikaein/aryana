static inline uint32_t bwtl_occ(const bwtl_t *bwt, uint32_t k, uint8_t c)
{
    uint32_t n, b;
    if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
    if (k == (uint32_t)(-1)) return 0;
    if (k >= bwt->primary) --k; // because $ is not in bwt
    n = bwt->occ[k/16<<2|c];
    b = bwt->bwt[k/16] & ~((1U<<((15-(k&15))<<1)) - 1);
    n += (bwt->cnt_table[b&0xff] + bwt->cnt_table[b>>8&0xff]
          + bwt->cnt_table[b>>16&0xff] + bwt->cnt_table[b>>24]) >> (c<<3) & 0xff;
    if (c == 0) n -= 15 - (k&15); // corrected for the masked bits
    return n;
}

static inline void bwtl_occ4(const bwtl_t *bwt, uint32_t k, uint32_t cnt[4])
{
    uint32_t x, b;
    if (k == (uint32_t)(-1)) {
        memset(cnt, 0, 16);
        return;
    }
    if (k >= bwt->primary) --k; // because $ is not in bwt
    memcpy(cnt, bwt->occ + (k>>4<<2), 16);
    b = bwt->bwt[k>>4] & ~((1U<<((~k&15)<<1)) - 1);
    x = bwt->cnt_table[b&0xff] + bwt->cnt_table[b>>8&0xff]
        + bwt->cnt_table[b>>16&0xff] + bwt->cnt_table[b>>24];
    x -= 15 - (k&15);
    cnt[0] += x&0xff;
    cnt[1] += x>>8&0xff;
    cnt[2] += x>>16&0xff;
    cnt[3] += x>>24;
}

static inline void bwtl_2occ4(const bwtl_t *bwt, uint32_t k, uint32_t l, uint32_t cntk[4], uint32_t cntl[4])
{
    bwtl_occ4(bwt, k, cntk);
    bwtl_occ4(bwt, l, cntl);
}

