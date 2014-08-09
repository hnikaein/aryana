#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <inttypes.h>

#include "utils.h"
#include "bwt.h"
//#include "align_bisulfite.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define maxGenomeSize 4e9
#define maxChromosomeNum 1000
//const long maxGenomeSize = 4e9;
//const int maxChromosomeNum = 1000;

const int interval_side = 3000;
long chromPos[maxChromosomeNum];
unsigned long gs;
char chromName[maxChromosomeNum][100];
char * genome;
int chromNum, xChromosome = -1, yChromosome = -1;
char * genomeFile, * annotationFile, * outputFile;



uint64_t islandStarts[(long)maxGenomeSize];
uint64_t islandEnds[(long)maxGenomeSize];
long islandsNum = 0;
int highPenalty=0,medPenalty=0,lowPenalty=0;
long penalties;
long readNum;
long *readPenalties;


char *samName, *referenceName, *annotationFile;
FILE *samFile;

int main(int argc, char *argv[]) {
	if (argc < 6) {
		fprintf(stderr, "Need more inputs\n");
		return -1;
	}

	
	static struct option long_options[] =
	{
		{"samfiles", required_argument, 0, 's'},
		{"reference", required_argument, 0, 'x'},
		{"CpG-islands", required_argument, 0, 'c'}
	//	{0, 0, 0, 0}
	};
	int option_index = 0;
	int c;
	while((c = getopt_long(argc, argv, "x:s:c:", long_options, &option_index)) >= 0){
		switch(c){
			case 0:
				break;
			case 1:
				break;
			case 'x':
				referenceName = (char *)malloc(strlen(optarg));
				strcpy(referenceName, optarg);
				break;
			case 's':
				samName = (char *)malloc(strlen(optarg));
				strcpy(samName, optarg);
				break;
			case 'c':
				annotationFile= (char *)malloc(strlen(optarg));
				strcpy(annotationFile, optarg);
				break;
		}
	}
	
    int i;
    ReadCpGIslands(annotationFile);
    

	fprintf(stderr, "salam");
	ref_read(referenceName);
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
	
	struct stat file_info;
	if(stat(referenceName , &file_info) == -1){
		fprintf(stderr, "Could not get the information of file %s\nplease make sure the file exists\n", file_name);
		return -1;
	}
	off_t file_size_bytes = file_info.st_size;
	readNum = (ceil((double)file_size_bytes/800)+200);
    readPenalties = malloc(readNum * sizeof(long));
    memset (readPenalties, 0, sizeof (long) * readNum);
	long readNumber= 0;
    
	while (fgets(line, 1000, samFile) != NULL) {
		if (!header) {
			sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, rname, &pos,&mapq, cigar,rnext,pnext, &tlen,seq_string,quality_string);
			//fprintf(stderr, "%s\t%s\n", cigar, seq_string);
            //printf("%s\n",cigar);
            readCigar(cigar,pos,seq_string,readNumber);
	        readNumber++;
		} else {
			if (line[0] == '@')
				fprintf(stderr, "header\t");
			else {
				header = 0;
				//fprintf(stderr, line);
				sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, rname, &pos,&mapq, cigar,rnext,pnext, &tlen,seq_string,quality_string);
				//fprintf(stderr, "salam");
				//fprintf(stderr, "\n%s\t%lld\t%s\n", cigar, tlen, seq_string);
				readCigar(cigar,pos,seq_string,readNumber);
				readNumber++;
			}
		}
	}
    
    
    printf("hhh : %d", isInIsland(711539));
//    for(i=0;i<100;i++){
//        printf("sssstart: %" PRIu64 "     eeeend : %" PRIu64 " ", islandStarts[i],islandEnds[i]);
//    }
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

	int i =0;
	for (i=0; i< 10; i++)
		fprintf(stderr, "%d\t", reference[i]);
	return 0;
}

char getNuc(uint64_t place, uint64_t seq_len){
   char atomic[4] = { 'A' , 'C' , 'G' , 'T'};
	int rev=0;
	fprintf(stderr, "AAAA %d %d\n", place, seq_len);
	// if(place > (seq_len / 2))
// 	{
// 		place = (seq_len / 2) - (place - (seq_len / 2))-1;
// 		rev=1;
// 	}
	uint64_t block=place/(sizeof(bwtint_t)*4);
	fprintf(stderr, "BBBB %d\n", block);
	int offset=place%(sizeof(bwtint_t)*4);
	fprintf(stderr, "CCCC %d\n", offset);
	uint64_t mask=3;
	mask=mask & (reference[block] >> (2*offset));
	fprintf(stderr, "DDDD %d\n", mask);
	if (rev==1)
		mask=3-mask;
	return atomic[mask];
}

inline void ToLower(char * s) {
    int i=0;
    for (i = 0; i < strlen(s); i++)
        if (s[i] >= 'A' && s[i] <= 'Z') s[i] += 'a' - 'A';
}

int ChromIndex(char * chrom) {
    ToLower(chrom);
    if (strstr(chrom, "chr") == chrom) chrom += 3;
    if (strstr(chrom, "ch") == chrom) chrom += 2;
	if (chrom[0] >= '0' && chrom[0] <= '9' && strlen(chrom) < 3)
		return atoi(chrom) - 1;
    if (strlen(chrom) > 1) return -1;
	if (chrom[0] == 'x') return xChromosome;
	else
	    if (chrom[0] == 'y') return yChromosome;
		else return -1;
}



int compare_function(const void *a,const void *b) {
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
   //	f.getline(fLine, sizeof(fLine)); // First row
   	long index = 0;
    while (fgets(fLine, sizeof(fLine), file)) {
		//fLine[0] = 0;
//        f.getline(fLine, sizeof(fLine));
        
		//if (! fLine[0]) continue;
		// cerr << fLine << endl;
//        if(islandsNum==0){
//            islandsNum++;
//            continue;
//        }
		wStart = 0;
        sscanf(fLine, "%s %" PRIu64 " %" PRIu64 " ", chrom, &wStart, &wEnd);

		if (! wStart) continue;
		// cerr << chrom << '\t' << wStart << '\t' << wEnd << endl;
		islandStarts[index] = wStart;
		islandEnds[index] = wEnd;
        //printf("starts: %" PRIu64 "     ends : %" PRIu64 "", islandStarts[index],islandEnds[index]);
        chr = ChromIndex(chrom);
        islandsNum++;
        index++;
    }
//    qsort(islandStarts,sizeof(islandStarts)/sizeof(*islandStarts),sizeof(*islandStarts),compare_function);
//    qsort(islandEnds,maxGenomeSize,sizeof(uint64_t),compare_function);
    //islandsNum--;
}

void setPenalties(int p1,int p2,int p3){
	highPenalty = p1;
	medPenalty = p2;
	lowPenalty = p3;
	return;
}

int isInIsland(uint64_t ref_i){
    
    printf("island refindex : %" PRIu64 "\n",ref_i);
    
	uint64_t first =0,last = islandsNum-1;
	uint64_t middle = (first+last)/2;
	int isInIsland = 0;
	while( first <= last )
	{
            printf("first :%" PRIu64 "  last :%" PRIu64 "",islandStarts[middle],islandEnds[middle]);
		if ( islandStarts[middle] >= ref_i && ref_i <= islandEnds[middle]){
			
			isInIsland = 1;
			return 1;
		}
        else if ( islandEnds[middle] < ref_i )
            first = middle + 1;
        else
            last = middle - 1;
        
        middle = (first + last)/2;
	}
	return 0;
}

void CalcPenalties(uint64_t ref_i , char read ,uint64_t seq_len , long readNum){
    //printf("   7salam");
     printf("read : %c   ",read);
    char atomic[4] = { 'A' , 'C' , 'G' , 'T'};
    //printf("read : %c    ref : %c  refindex : %" PRIu64 "\n",read,getNuc(ref_i,seq_len) ,ref_i);
	if(read == 'A' || read =='G'){
		if(getNuc(ref_i,seq_len) != read )
			readPenalties[readNum] += highPenalty;
	}
	else{ //read = C or T
        
		if(read == 'T' && getNuc(ref_i,seq_len) == 'C'){
			if(getNuc(ref_i + 1,seq_len) =='G'){ // in the CpG context
				if(isInIsland(ref_i) ){ // in CpG and also island
					readPenalties[readNum] += medPenalty;
				}
				else{
					readPenalties[readNum] += highPenalty;
				}
			}
			else{ // out of CpG context
                
				readPenalties[readNum] += lowPenalty;
			}
		}
		else if(read =='C' && getNuc(ref_i,seq_len) == 'C'){
			if(getNuc(ref_i + 1,seq_len) =='G'){ // in the CpG context
				int temp = isInIsland(ref_i);
				if(temp == 1 ){ // in CpG and also island
					readPenalties[readNum] += medPenalty;
				}
				else{
					readPenalties[readNum] += lowPenalty;
				}
			}
			else{ // out of CpG context
				readPenalties[readNum] += highPenalty;
			}
            
            
		}
		else if(read == 'C' && getNuc(ref_i,seq_len) == 'T')
            readPenalties[readNum] += highPenalty;
        
        
	}
       
}
void readCigar(char * cigar,uint64_t ref_i,char *seq_string,long readNum){
    fprintf(stderr, "salam");
    int pos = 0;
    int value = 0;
    uint64_t ref_index = ref_i;
    long read_index = 0;
    char alignType ;
                  printf("   %s",seq_string);
    while ( 1 ) {
        if ( !isdigit(cigar[pos]) ) {
            	//printf("   1salam");
            if ( value > 0 ){
  
                //printf("value:   %d",value);
                if(cigar[pos] == 'm'){
                    int j;
                    for (j=0; j < value; j++) {
                        //printf("   71salam");
                        CalcPenalties(++ref_index ,seq_string[read_index++], reference_size , readNum );
                    }
                }
                else if (cigar[pos] == 'd') {
                    ref_index += value;        
                }
                else if (cigar[pos] == 'i')
                    read_index += value;
                else{
                    printf("*");
                    break;
                }
 
            }
            printf("*");
            if ( cigar[pos]==0 )
                break;
//            else{
//                printf("   3salam");
//                alignType = cigar[pos];
//                printf("  %c",alignType);
//            }
            value = 0;
        } else {
            value = value * 10 + cigar[pos]-'0';
        }
        pos++;
    }
}
