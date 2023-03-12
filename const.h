#define read_size 600000
#define seq_num_per_read 10 // xxx 5000
#define max_sam_line 3*MAX_READ_SIZE+1000
#define MAX_CIGAR_SIZE 6000
#define output_buffer (1 << 9)
#define sleep_time (1)
#define true 1
#define false 0
#define sleep_time (1)
#define maxPenalty (100000000)
#define max_chrom_num (10000)
#define max_sam_header_size (1000000)

#define max(a, b) (((a)>(b))?(a):(b))

#define RAND_SEED (14)
#define MAX_READ_LENGTH 10000
#define MAX_READ_NAME_LENGTH 1000
#define MAX_CHR_NAME_LENGTH 1000
#define MAX_SAM_LINE_LENGTH (4*(MAX_READ_LENGTH)+(MAX_READ_NAME_LENGTH)+(2*MAX_CHR_NAME_LENGTH)+1000)

#define BS_GENOMES_COUNT 5
#define UNALIGNED_PENALTY  1000000
#define DISCORD_PENALTY  100000