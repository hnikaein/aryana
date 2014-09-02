#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <limits.h>

#include "utils.h"
#include "bwt.h"
//#include "align_bisulfite.h"

#include <stdio.h>
#include <stdlib.h>

#define maxGenomeSize 1e5
#define maxChromosomeNum 1000
//const long maxGenomeSize = 4e9;
//const int maxChromosomeNum = 1000;

const int interval_side = 3000;
long chromPos[maxChromosomeNum];
unsigned long gs;
char chromName[maxChromosomeNum][100];
char * genome;
int chromNum, xChromosome = -1, yChromosome = -1;
char * genomeFile, *annotationFile, *outputFile;
uint64_t islandStarts[(long) maxGenomeSize];
uint64_t islandEnds[(long) maxGenomeSize];
long islandsNum = 0;
int highPenalty = 0, medPenalty = 0, lowPenalty = 0;
long penalties;
long readNum;
long readPenalties[4];

char *referenceName, *annotationFile;
char *samNames[4];
FILE *samFiles[4];


struct ChrIslands
{
	int chrNum;
	int islandsNum;
	uint64_t islandStarts[(long) maxGenomeSize];
	uint64_t islandEnds[(long) maxGenomeSize];
    
};

struct ChrIslands chrIslands[maxChromosomeNum] ;

int count =0; ////////////////////////////////
int main(int argc, char *argv[]) {
        //printf("heu1");
	if (argc < 6) {
		fprintf(stderr, "Need more inputs\n");
		return -1;
	}

	static struct option long_options[] = {
		{ "samfiles", required_argument, 0, 's' },
		{ "reference", required_argument, 0, 'x' },
		{ "CpG-islands", required_argument, 0, 'c' },
		{ "penalties", required_argument, 0, 'p' }
	//	{0, 0, 0, 0}
			};
	int option_index = 0;
	int c;
	while ((c = getopt_long(argc, argv, "x:s:c:p:", long_options, &option_index))
			>= 0) {
            //printf("heu0");
		switch (c) {
		case 0:
			break;
		case 1:
			break;
		case 'x':
			referenceName = (char *) malloc(strlen(optarg));
			strcpy(referenceName, optarg);
			break;
		case 's':
			samNames[0] = (char *) malloc(strlen(optarg));
			strcpy(samNames[0], optarg);
			samNames[1] = (char *) malloc(strlen(argv[optind]));
			strcpy(samNames[1], argv[optind]);
			samNames[2] = (char *) malloc(strlen(argv[optind + 1]));
			strcpy(samNames[2], argv[optind + 1]);
			samNames[3] = (char *) malloc(strlen(argv[optind + 2]));
			strcpy(samNames[3], argv[optind + 2]);
			optind = optind + 3;
			break;
		case 'c':
			annotationFile = (char *) malloc(strlen(optarg));
			strcpy(annotationFile, optarg);
			break;
		case 'p':
			lowPenalty = atoi(optarg);
			medPenalty = atoi(argv[optind]);
			highPenalty = atoi(argv[optind + 1]);
			optind = optind + 2;
			break;
		}
	}
   // printf("heu1");
	ReadCpGIslands(annotationFile);
    int p=0;
    for(p;p<100;p++)
        fprintf(stderr,"start:  %lld    end : %lld number:  %d \n",chrIslands[4].islandStarts[p],chrIslands[4].islandEnds[p],chrIslands[4].chrNum);
        //printf("heu2");
	//fprintf(stderr, "salam");
	ref_read(referenceName);
       // printf("heu3");
	char * command = "sort -k1 ";
	char cmd_pointer[strlen(samNames[0]) + strlen(command) + 30];
	strcpy(cmd_pointer, command);
	strcat(cmd_pointer, samNames[0]);
	strcat(cmd_pointer, " > samFile1.sam");
	system(cmd_pointer);

	char cmd_pointer2[strlen(samNames[1]) + strlen(command) + 1];
	strcpy(cmd_pointer2, command);
	strcat(cmd_pointer2, samNames[1]);
	strcat(cmd_pointer2, " > samFile2.sam");
	system(cmd_pointer2);

	char cmd_pointer3[strlen(samNames[2]) + strlen(command) + 1];
	strcpy(cmd_pointer3, command);
	strcat(cmd_pointer3, samNames[2]);
	strcat(cmd_pointer3, " > samFile3.sam");
	system(cmd_pointer3);

	char cmd_pointer4[strlen(samNames[3]) + strlen(command) + 1];
	strcpy(cmd_pointer4, command);
	strcat(cmd_pointer4, samNames[3]);
	strcat(cmd_pointer4, " > samFile4.sam");
	system(cmd_pointer4);

	samFiles[0] = fopen("samFile1.sam", "r");
	samFiles[1] = fopen("samFile2.sam", "r");
	samFiles[2] = fopen("samFile3.sam", "r");
	samFiles[3] = fopen("samFile4.sam", "r");
	char line[1000];
	int header = 1;
	char *qname, *rnext, *pnext, *seq_string, *quality_string;
	char *rname[4], *cigar[4];
	int flag, i;
	uint64_t pos[4];
	uint32_t mapq[4];
	long long int tlen;
	for (i = 0; i < 4; i++) {
		rname[i] = malloc(100 * sizeof(char));
		cigar[i] = malloc(200 * sizeof(char));
	}
	qname = malloc(100 * sizeof(char));
	rnext = malloc(100 * sizeof(char));
	pnext = malloc(100 * sizeof(char));
	seq_string = malloc(1000 * sizeof(char));
	quality_string = malloc(500 * sizeof(char));
	int chosen[4];
    int j=0;
    for (j; j < 4; j++)
		chosen[j] = 0;


	int stop = 0;
	while (1 && !stop) {
		for (i = 0; i < 4 && !stop; i++) {
			if (fgets(line, 1000, samFiles[i]) == NULL) {
				stop = 1;
				break;
			}
			while (line[0] == '@'){
				if(fgets(line, 1000, samFiles[i]) == NULL){
					stop = 1;
					break;
				}
			}
			if(stop)
				break;
			sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, rname[i], &pos[i],&mapq[i], cigar[i],rnext,pnext, &tlen,seq_string,quality_string);
			//fprintf(stderr, "AAA %s\n", qname);
            //printf("cigar : %s \n",cigar[i]);
			readCigar(cigar[i], pos[i], seq_string, i ,ChromIndex(rname[i]));
		}
		if(stop)
			break;
		int min = min_penalty();
		chosen[min]++;
		fprintf(stdout, "%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, flag, rname[min], pos[min],mapq[min], cigar[min],rnext,pnext, tlen,seq_string,quality_string);
        int j=0;
        for (j; j < 4; j++)
            readPenalties[j]=0;
	}
	
    for (j=0; j < 4; j++)
		fprintf(stderr, "CHOSEN %d: %d\n", j, chosen[j]);


//	while (fgets(line, 1000, samFile1) != NULL) {
//		if (!header) {
//			sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, rname, &pos,&mapq, cigar,rnext,pnext, &tlen,seq_string,quality_string);
//			//fprintf(stderr, "%s\t%s\n", cigar, seq_string);
//			//printf("%s\n",cigar);
//			readCigar(cigar, pos, seq_string, readNumber);
//			readNumber++;
//		} else {
//			if (line[0] == '@')
//				fprintf(stderr, "header\t");
//			else {
//				header = 0;
//				//fprintf(stderr, line);
//				sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, rname, &pos,&mapq, cigar,rnext,pnext, &tlen,seq_string,quality_string);
//				//fprintf(stderr, "salam");
//				//fprintf(stderr, "\n%s\t%lld\t%s\n", cigar, tlen, seq_string);
//				readCigar(cigar, pos, seq_string, readNumber);
//				readNumber++;
//			}
//		}
//	}
//
//	printf("hhh : %d", isInIsland(711539));
//    for(i=0;i<100;i++){
//        printf("sssstart: %" PRIu64 "     eeeend : %" PRIu64 " ", islandStarts[i],islandEnds[i]);
//    }
}
uint64_t * reference;
uint32_t reference_size;
uint32_t reference_reminder;

int ref_read(char * file_name){
	//fprintf(stderr, "salam %s", file_name);
	fprintf(stderr, "inside ref_read with %s\n", file_name);
	struct stat file_info;
	if (stat(file_name, &file_info) == -1) {
		fprintf(stderr,
				"Could not get the information of file %s\nplease make sure the file exists\n",
				file_name);
		return -1;
	}
	int fd = open(file_name, O_RDONLY);
	//FILE *fd;
	//fd = xopen(file_name, "rb");
	if (fd == -1) {
		fprintf(stderr,
				"Could not open the file %s\nplease make sure the file exists\n",
				file_name);
		return -1;
	}
	off_t file_size_bytes = file_info.st_size;
	//fprintf(stderr, "size: %d\n", file_size_bytes);
	reference_size = ceil(
			((double) file_size_bytes) / (double) (sizeof(uint64_t)));
	fprintf(stderr, "reference_size = %u\n", reference_size);
	reference_reminder = file_size_bytes % sizeof(uint64_t);
	//reference = new base64 [ reference_size ];
	reference = (uint64_t *) malloc(reference_size * sizeof(uint64_t));
	memset(reference, 0, reference_size * sizeof(uint64_t));
	size_t read_size2 = 0; //there is a read_size defined above
	size_t signal;
	size_t total_size = (file_size_bytes);
	unsigned char *uc_buffer = (unsigned char *) (reference);
	int counter = 0;

	do {
		signal = read(fd, (void *) uc_buffer, total_size - read_size2);
		//signal = fread((void *)uc_buffer, )
		if (signal == -1) {
			fprintf(stderr, "Error: while writing to file\n");
			if (close(fd) == -1)
				fprintf(stderr, "Error: while closing file\n");
			return -1;
		}
		counter++;
		read_size2 += signal;
		uc_buffer += signal;
	} while (read_size2 < total_size);
	if (close(fd) == -1) {
		fprintf(stderr, "Unable to close the file\n");
		return -1;
	}
	return 0;
}


int min_penalty(){
	int i = 0,j;
    if(count++<50)
        for(j=0;j<4;j++)
            fprintf(stderr, "Penalties %d : %lld \n",j,readPenalties[j]);
	long min = readPenalties[0];
	if(min > readPenalties[1]){
		min = readPenalties[1];
		i = 1;
	}
	if(min > readPenalties[2]){
		min = readPenalties[2];
		i = 2;
	}
	if(min > readPenalties[3]){
		min = readPenalties[3];
		i = 3;
	}
	return i;
}
char getNuc(uint64_t place, uint64_t seq_len) {
	char atomic[4] = { 'A', 'C', 'G', 'T' };
	int rev = 0;
	//fprintf(stderr, "AAAA %" PRIu64 " %" PRIu64 "\n", place, seq_len);
	// if(place > (seq_len / 2))
// 	{
// 		place = (seq_len / 2) - (place - (seq_len / 2))-1;
// 		rev=1;
// 	}
	uint64_t block = place / (sizeof(bwtint_t) * 4);
	//fprintf(stderr, "BBBB %" PRIu64 "\n", block);
	int offset = place % (sizeof(bwtint_t) * 4);
	//fprintf(stderr, "CCCC %d\n", offset);
	uint64_t mask = 3;
	//mask = mask & (reference[block] >> (2 * offset));
	//fprintf(stderr, "DDDD %" PRIu64 "\n", mask);
	if (rev == 1)
		mask = 3 - mask;
	//return atomic[mask];
	return 'A';
}

inline void ToLower(char * s) {
	int i = 0;
	for (i = 0; i < strlen(s); i++)
		if (s[i] >= 'A' && s[i] <= 'Z')
			s[i] += 'a' - 'A';
}

int ChromIndex(char * chrom) {
	ToLower(chrom);
	if (strstr(chrom, "chr") == chrom)
		chrom += 3;
	if (strstr(chrom, "ch") == chrom)
		chrom += 2;
	if (chrom[0] >= '0' && chrom[0] <= '9' && strlen(chrom) < 3)
		return atoi(chrom) - 1;
	if (strlen(chrom) > 1)
		return -1;
	if (chrom[0] == 'x')
		return xChromosome;
	else if (chrom[0] == 'y')
		return yChromosome;
	else
		return -1;
}

int compare_function(const void *a, const void *b) {
	int *x = (uint64_t *) a;
	int *y = (uint64_t *) b;
	return *x - *y;
}

void ReadCpGIslands(char * annotationFile) {
	// cerr << "Processing CpG island locations from file: " <<  annotationFile << endl;
	// ifstream f(annotationFile);
    //    if (! f.is_open()) {
    //       // cerr << "Error: CpG island locations file not found or could not be opened" << endl;
    //        exit(1);
    //    }
	fprintf(stderr, "salam");
	FILE* file = fopen(annotationFile, "r"); /* should check the result */
	fprintf(stderr, "salam2");
	char fLine[10000], chrom[10], strand[10];
	uint64_t wStart, wEnd, chr;
	int bin;
    
	long index = 0;
	fgets(fLine, sizeof(fLine), file);
	fgets(fLine, sizeof(fLine), file);
	wStart = 0;
	sscanf(fLine, "%s %" PRIu64 " %" PRIu64 " ", chrom, &wStart, &wEnd);
	fprintf(stderr, "%s \n",fLine);
	long chrIndex =0;
	int stop = 0;
	while (!stop) {
		char *copy = (char *)malloc(strlen(chrom));
		memcpy(copy, chrom,strlen(chrom));
		index = 0;
		while(!strcmp(copy, chrom)){
            
			if (!wStart)
				continue;
            
			// cerr << chrom << '\t' << wStart << '\t' << wEnd << endl;
			//islandStarts[index] = wStart;
			chrIslands[chrIndex].islandStarts[index] = wStart;
			chrIslands[chrIndex].islandEnds[index] = wEnd;
			printf("starts: %" PRIu64 "     ends : %" PRIu64 "", chrIslands[chrIndex].islandStarts[index],chrIslands[chrIndex].islandEnds[index]);
			chrIslands[chrIndex].chrNum = ChromIndex(chrom);
			chrIslands[chrIndex].islandsNum++;
			index++;
			if(NULL==fgets(fLine, sizeof(fLine), file)){
				stop =1;
				break;
			}
			sscanf(fLine, "%s %" PRIu64 " %" PRIu64 " ", chrom, &wStart, &wEnd);
			
		}
		chrIndex++;
		index=0;
        
	}
}

void setPenalties(int p1, int p2, int p3) {
	highPenalty = p1;
	medPenalty = p2;
	lowPenalty = p3;
	return;
}

int isInIsland(uint64_t ref_i , int chr) {

	//printf("island refindex : %" PRIu64 "\n",ref_i);
    int i;
    int chr2;
    for(i=0;i<maxChromosomeNum ; i++)
        if(chrIslands[i].chrNum == chr){
            chr2 = i;
            break;
        }
	uint64_t first = 0, last = chrIslands[chr2].islandsNum - 1;
	uint64_t middle = (first + last) / 2;
	int isInIsland = 0;
	while (first <= last) {
		//printf("first :%" PRIu64 "  last :%" PRIu64 "",islandStarts[middle],islandEnds[middle]);
		if (chrIslands[chr2].islandStarts[middle] >= ref_i && ref_i <= chrIslands[chr2].islandEnds[middle]) {

			isInIsland = 1;
			return 1;
		} else if (chrIslands[chr2].islandEnds[middle] < ref_i)
			first = middle + 1;
		else
			last = middle - 1;

		middle = (first + last) / 2;
	}
	return 0;
}

void CalcPenalties(uint64_t ref_i, char read, uint64_t seq_len, long readNum,int chr) {
	//printf("   7salam");
	//printf("read : %c   refrence: %c \n ", read,getNuc(ref_i, seq_len));
	//printf("read : %c   ", read);
	char atomic[4] = { 'A', 'C', 'G', 'T' };
	//printf("read : %c    ref : %c  refindex : %" PRIu64 "\n",read,getNuc(ref_i,seq_len) ,ref_i);
	if (read == 'A' || read == 'G') {
		if (getNuc(ref_i, seq_len) != read)
			readPenalties[readNum] += highPenalty;
	} else { //read = C or T

		if (read == 'T' && getNuc(ref_i, seq_len) == 'C') {
			if (getNuc(ref_i + 1, seq_len) == 'G') { // in the CpG context
				if (isInIsland(ref_i,chr)) { // in CpG and also island
					readPenalties[readNum] += medPenalty;
				} else {
					readPenalties[readNum] += highPenalty;
				}
			} else { // out of CpG context

				readPenalties[readNum] += lowPenalty;
			}
		} else if (read == 'C' && getNuc(ref_i, seq_len) == 'C') {
			if (getNuc(ref_i + 1, seq_len) == 'G') { // in the CpG context
				int temp = isInIsland(ref_i,chr);
				if (temp == 1) { // in CpG and also island
					readPenalties[readNum] += medPenalty;
				} else {
					readPenalties[readNum] += lowPenalty;
				}
			} else { // out of CpG context
				readPenalties[readNum] += highPenalty;
			}

		} else if (read == 'C' && getNuc(ref_i, seq_len) == 'T')
			readPenalties[readNum] += highPenalty;

	}

}
void readCigar(char * cigar, uint64_t ref_i, char *seq_string, long readNum,int chr) {
	//fprintf(stderr, "salam\n");
	int pos = 0;
	int value = 0;
	uint64_t ref_index = ref_i;
	long read_index = 0;
	char alignType;
	//printf("   %s\n", seq_string);
    //printf("cigar:   %s\n", cigar);
	while (1) {
		if (!isdigit(cigar[pos])) {
			//printf("   1salam\n");
			if (value > 0) {

				//printf("value:   %d",value);
				if (cigar[pos] == 'm') {
					int j;
					for (j = 0; j < value; j++) {
						//printf("   71salam\n");
						CalcPenalties(++ref_index, seq_string[read_index++],
								reference_size, readNum,chr);
                        //          if(strstr(cigar,"79m1d9m1d5m1i6m"))
                        //                  fprintf(stderr,"   penalties for cigar : %lld   ",readPenalties[readNum]);
					}
				} else if (cigar[pos] == 'd') {
					ref_index += value;
                    readPenalties[readNum] += highPenalty * value;

				} else if (cigar[pos] == 'i'){
					read_index += value;
                    readPenalties[readNum] += highPenalty * value;

                }

			}
            else if(cigar[pos]=='*'){
                readPenalties[readNum] += LONG_MAX;
                break;
            }
//			printf("*");
			if (cigar[pos] == 0)
				break;
//            else{
//                printf("   3salam");
//                alignType = cigar[pos];
//                printf("  %c",alignType);
//            }
			value = 0;
            } else {
                value = value * 10 + cigar[pos] - '0';
            }
            pos++;
	}
    if(count <20)
        fprintf(stderr, "read : %s \n, cigar : %s \n , penalties : %lld \n",seq_string,cigar,readNum);
}
    
