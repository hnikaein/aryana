//YA HAGh
#include "bwt.h"

#define max_refs 10000

typedef enum {
    sf_paired = 0x1,
    sf_pair_mapped = 0x2,
    sf_unmapped = 0x4,
    sf_mate_unmapped = 0x8,
    sf_reverse = 0x10,
    sf_mate_reverse = 0x20,
    sf_first_mate = 0x40,
    sf_second_mate = 0x80,
    sf_not_primary = 0x100,
    sf_fail_quality = 0x200,
    sf_pcr = 0x400,
    sf_supplement = 0x800
} sam_flags;

int load_reference(const char *filename);

typedef struct {
    char *rname;
    bwtint_t pos;
} genomic_location;

int sam_generator(char *buffer, char *qname, int flag, uint32_t mapq, bwtint_t index, char *cigar, bwtint_t inext2,
                  char *cigar2, ubyte_t *seq, ubyte_t *quality, int len, int len2, bwt_t *bwt, bwtint_t *offset,
                  bwtint_t offInd);

char name[max_refs][1000];

//int sam_headers(char * buffer,char ** name, long long*  offset, int size);
int sam_headers(char *buffer, bwtint_t *offset, int size, int buf_size);

int get_chr(bwtint_t index, bwtint_t seq_len, bwtint_t *offset, bwtint_t offInd);
