static inline int __occ_aux(uint64_t y, int c) {
    // reduce nucleotide counting to bits counting
    y = ((c & 2) ? y : ~y) >> 1 & ((c & 1) ? y : ~y) & 0x5555555555555555ull;
    // count the number of 1s in y
    y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
    return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

static inline bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c) {
    bwtint_t n, l, j;
    uint32_t *p;

    if (k == bwt->seq_len) return bwt->L2[c + 1] - bwt->L2[c];
    if (k == (bwtint_t)(-1)) return 0;
    if (k >= bwt->primary) --k; // because $ is not in bwt

    // retrieve Occ at k/OCC_INTERVAL
    n = ((bwtint_t * )(p = bwt_occ_intv(bwt, k)))[c];
    p += sizeof(bwtint_t); // jump to the start of the first BWT cell

    // calculate Occ up to the last k/32
    j = k >> 5 << 5;
    for (l = k / OCC_INTERVAL * OCC_INTERVAL; l < j; l += 32, p += 2)
        n += __occ_aux((uint64_t) p[0] << 32 | p[1], c);

    // calculate Occ
    n += __occ_aux(((uint64_t) p[0] << 32 | p[1]) & ~((1ull << ((~k & 31) << 1)) - 1), c);
    if (c == 0) n -= ~k & 31; // corrected for the masked bits

    return n;
}

// an analogy to bwt_occ() but more efficient, requiring k <= l
static inline void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol) {
    bwtint_t _k, _l;
    _k = (k >= bwt->primary) ? k - 1 : k;
    _l = (l >= bwt->primary) ? l - 1 : l;
    if (_l / OCC_INTERVAL != _k / OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
        *ok = bwt_occ(bwt, k, c);
        *ol = bwt_occ(bwt, l, c);
    } else {
        bwtint_t m, n, i, j;
        uint32_t *p;
        if (k >= bwt->primary) --k;
        if (l >= bwt->primary) --l;
        n = ((bwtint_t * )(p = bwt_occ_intv(bwt, k)))[c];
        p += sizeof(bwtint_t);
        // calculate *ok
        j = k >> 5 << 5;
        for (i = k / OCC_INTERVAL * OCC_INTERVAL; i < j; i += 32, p += 2)
            n += __occ_aux((uint64_t) p[0] << 32 | p[1], c);
        m = n;
        n += __occ_aux(((uint64_t) p[0] << 32 | p[1]) & ~((1ull << ((~k & 31) << 1)) - 1), c);
        if (c == 0) n -= ~k & 31; // corrected for the masked bits
        *ok = n;
        // calculate *ol
        j = l >> 5 << 5;
        for (; i < j; i += 32, p += 2)
            m += __occ_aux((uint64_t) p[0] << 32 | p[1], c);
        m += __occ_aux(((uint64_t) p[0] << 32 | p[1]) & ~((1ull << ((~l & 31) << 1)) - 1), c);
        if (c == 0) m -= ~l & 31; // corrected for the masked bits
        *ol = m;
    }
}

#define __occ_aux4(bwt, b)                                                                                      \
        ((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]             \
         + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

static inline void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4]) {
    bwtint_t l, j, x;
    uint32_t *p;
    if (k == (bwtint_t)(-1)) {
        memset(cnt, 0, 4 * sizeof(bwtint_t));
        return;
    }
    if (k >= bwt->primary) --k; // because $ is not in bwt
    p = bwt_occ_intv(bwt, k);
    memcpy(cnt, p, 4 * sizeof(bwtint_t));
    p += sizeof(bwtint_t);
    j = k >> 4 << 4;
    for (l = k / OCC_INTERVAL * OCC_INTERVAL, x = 0; l < j; l += 16, ++p)
        x += __occ_aux4(bwt, *p);
    x += __occ_aux4(bwt, *p & ~((1U << ((~k & 15) << 1)) - 1)) - (~k & 15);
    cnt[0] += x & 0xff;
    cnt[1] += x >> 8 & 0xff;
    cnt[2] += x >> 16 & 0xff;
    cnt[3] += x >> 24;
}

// an analogy to bwt_occ4() but more efficient, requiring k <= l
static inline void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4]) {
    bwtint_t _k, _l;
    _k = (k >= bwt->primary) ? k - 1 : k;
    _l = (l >= bwt->primary) ? l - 1 : l;
    if (_l / OCC_INTERVAL != _k / OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
        bwt_occ4(bwt, k, cntk);
        bwt_occ4(bwt, l, cntl);
    } else {
        bwtint_t i, j, x, y;
        uint32_t *p;
        if (k >= bwt->primary) --k; // because $ is not in bwt
        if (l >= bwt->primary) --l;
        p = bwt_occ_intv(bwt, k);
        memcpy(cntk, p, 4 * sizeof(bwtint_t));
        p += sizeof(bwtint_t);
        // prepare cntk[]
        j = k >> 4 << 4;
        for (i = k / OCC_INTERVAL * OCC_INTERVAL, x = 0; i < j; i += 16, ++p)
            x += __occ_aux4(bwt, *p);
        y = x;
        x += __occ_aux4(bwt, *p & ~((1U << ((~k & 15) << 1)) - 1)) - (~k & 15);
        // calculate cntl[] and finalize cntk[]
        j = l >> 4 << 4;
        for (; i < j; i += 16, ++p) y += __occ_aux4(bwt, *p);
        y += __occ_aux4(bwt, *p & ~((1U << ((~l & 15) << 1)) - 1)) - (~l & 15);
        memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
        cntk[0] += x & 0xff;
        cntk[1] += x >> 8 & 0xff;
        cntk[2] += x >> 16 & 0xff;
        cntk[3] += x >> 24;
        cntl[0] += y & 0xff;
        cntl[1] += y >> 8 & 0xff;
        cntl[2] += y >> 16 & 0xff;
        cntl[3] += y >> 24;
    }
}
