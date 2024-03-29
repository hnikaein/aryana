#include <stdint.h>
#include "bwt.h"
#include "const.h"

#ifndef hash_header
#define hash_header
typedef struct {
    uint64_t place; //index of 2l
    int value;
    uint64_t level;
//	uint64_t width;
    uint64_t index; //real matched index
    //char cigar[100];
    //uint64_t left; //
    //uint64_t old_left;
    uint64_t last; //last seed that updtates left
    uint64_t groupid;
    uint64_t match_start[SEED_NUMS_PER_READ];// start positions of the seeds in read
    uint64_t match_index[SEED_NUMS_PER_READ];// start positions of the seeds in reference
    uint64_t matched[SEED_NUMS_PER_READ]; // lengths of the seeds
    int parts; // Number of seeds matched to a tag
    int mate;
} hash_element;

typedef struct {
    int paired;
    hash_element *pair_table;
    uint64_t min_dis, max_dis;
    uint64_t len1, len2;
    int pairID;
    char ori[3];
    int strict;
} pair_opt;


void add(bwt_t *const bwt, uint64_t place, uint64_t value, uint64_t level, uint64_t index, int *best, int best_size,
         int *best_found, hash_element *table,
         uint64_t read_start, uint64_t read_size, uint64_t read_index, int whole_read_size, uint64_t groupid);

void reset_hash(hash_element *table);

uint64_t uminus(uint64_t x, uint64_t y);

#endif
