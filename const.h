#define seq_num_per_file_read 7
#define MAX_READ_LEN (1000 + 50)
#define HASH_TABLE_SIZE 200
#define SEED_NUMS_PER_READ 100
#define MAX_SAM_LINE (4 * MAX_READ_LEN + 2000)
#define MAX_CIGAR_SIZE (4 * MAX_READ_LEN)
#define MAX_POTENTS 200
#define MAX_THREADS_COUNT 100
#define sleep_time 1
#define true 1
#define false 0
#define maxPenalty 100000000
#define max_chrom_num 10000
#define max_sam_header_size 1000000

#define max(a, b) (((a)>(b))?(a):(b))

#define RAND_SEED (14)
#define MAX_READ_LENGTH 10000
#define MAX_READ_NAME_LENGTH 1000
#define MAX_CHR_NAME_LENGTH 1000
#define MAX_SAM_LINE_LENGTH (4*(MAX_READ_LENGTH)+(MAX_READ_NAME_LENGTH)+(2*MAX_CHR_NAME_LENGTH)+1000)

#define ME_BUFFER_SIZE (2 * MAX_READ_LENGTH + 1) // Should be higher than maximum read length
#define BS_GENOMES_COUNT 5
#define UNALIGNED_PENALTY  1000000
#define DISCORD_PENALTY  100000