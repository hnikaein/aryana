#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "utils.h"
#include "bwt.h"

char* samName, *referenceName;
FILE *samFile;

int main(int argc, char *argv[]) {
//	if (argc < 3) {
//		fprintf(stderr, "Need more inputs\n");
//		return -1;
//	}
	referenceName = "hgTest.fa";
	fprintf(stderr, "salam");
	ref_read(referenceName);
	samName = "a.sam";
	samFile = fopen(samName, "r");
	char line[1000];
	int header = 1;
	char *qname, *rname, *cigar, *rnext, *pnext, *seq_string, *quality_string;
	int flag;
	uint64_t pos;
	uint32_t mapq;
	long long int tlen;
	qname = malloc(100 * sizeof(char));
	rname = malloc(100 * sizeof(char));
	cigar = malloc(200 * sizeof(char));
	rnext = malloc(100 * sizeof(char));
	pnext = malloc(100 * sizeof(char));
	seq_string = malloc(1000 * sizeof(char));
	quality_string = malloc(500 * sizeof(char));
	while (fgets(line, 1000, samFile) != NULL) {
		if (!header) {
			sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, rname, &pos,&mapq, cigar,rnext,pnext, &tlen,seq_string,quality_string);
			fprintf(stderr, "%s\t%s\n", cigar, seq_string);
		} else {
			if (line[0] == '@')
				fprintf(stderr, "header\t");
			else {
				header = 0;
				fprintf(stderr, line);
				sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, rname, &pos,&mapq, cigar,rnext,pnext, &tlen,seq_string,quality_string);
				fprintf(stderr, "salam");
				fprintf(stderr, "\n%s\t%lld\t%s\n", cigar, tlen, seq_string);
			}
		}
	}
}
uint64_t * reference;
uint32_t  reference_size;
uint32_t  reference_reminder;

int ref_read(char * file_name){
	fprintf(stderr, "salam %s", file_name);
	fprintf(stderr, "inside ref_read with %s\n", file_name);
	struct stat file_info;
	if(stat(file_name , &file_info) == -1){
		fprintf(stderr, "Could not get the information of file %s\nplease make sure the file exists\n", file_name);
		return -1;
	}
	int fd = open(file_name , O_RDONLY);
	//FILE *fd;
	//fd = xopen(file_name, "rb");
	if(fd == -1){
		fprintf(stderr, "Could not open the file %s\nplease make sure the file exists\n", file_name);
		return -1;
	}
	off_t file_size_bytes = file_info.st_size;
	reference_size = ceil ( ((double)file_size_bytes) / (double)(sizeof(uint64_t)) );
	fprintf(stderr, "reference_size = %u\n", reference_size);
	reference_reminder = file_size_bytes % sizeof(uint64_t) ;
	//reference = new base64 [ reference_size ];
	reference = (uint64_t *)malloc(reference_size * sizeof(uint64_t));
	memset ( reference , 0 , reference_size * sizeof(uint64_t) );
	size_t read_size2 = 0;//there is a read_size defined above
	size_t signal;
	size_t total_size = (file_size_bytes);
	unsigned char *uc_buffer = (unsigned char *)(reference);
	int counter=0;

	do{
		signal = read ( fd , (void *)uc_buffer , total_size - read_size2 );
		//signal = fread((void *)uc_buffer, )
		if ( signal == -1 )
		{
			fprintf(stderr, "Error: while writing to file\n");
			if ( close(fd) == -1 )
				fprintf(stderr, "Error: while closing file\n");
			return -1;
		}
		counter++;
		read_size2 += signal;
		uc_buffer += signal;
	}
	while ( read_size2 < total_size );
	if ( close(fd) == -1 )
	{
		fprintf(stderr, "Unable to close the file\n");
		return -1;
	}

	return 0;
}

char getNuc(uint64_t place, uint64_t seq_len){
	int rev=0;
	if(place > (seq_len / 2))
	{
		place = (seq_len / 2) - (place - (seq_len / 2))-1;
		rev=1;
	}
	uint64_t block=place/(sizeof(bwtint_t)*4);
	int offset=place%(sizeof(bwtint_t)*4);
	uint64_t mask=3;
	mask=mask & (reference[block] >> (2*offset));
	if (rev==1)
		mask=3-mask;
	return mask;//atom[mask];
}
