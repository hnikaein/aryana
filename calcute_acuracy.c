//
//  calcute_acuracy.c
//  
//
//  Created by Maryam Rabiee on 8/21/14.
//
//

#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>


char * run(const char* cmd){
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[262144];
    char data[100];
    char *result;
    int dist=0;
    int size;
    //TIME_START
    while(!feof(pipe)) {
        size=(int)fread(buffer,1,262144, pipe); //cout<<buffer<<" size="<<size<<endl;
        memcpy(data,buffer,100);
    }
    //TIME_PRINT_
    fprintf(stdout,"%s \n",data);
    pclose(pipe);
    return data;
}

int main(int argc, char *argv[]) {
    
    
    if ( argc != 2 ){
        /* We print argv[0] assuming it is the program name */
        printf( "usage: %s filename", argv[0] );
    }
    else
    {
        // We assume argv[1] is a filename to open
        FILE *samFile = fopen( argv[1], "r" );
        /* fopen returns 0, the NULL pointer, on failure */
        if ( samFile == 0 )
            printf( "Could not open file\n" );
    
//    FILE *samFile;
//    samFile = fopen(samName, "r");
    long long badAlignedReads = 0;
    long long notAlignedReads = 0;
    long long readNum = 0;
    char line[1000];
    int header = 1;
    char *qname, *rnext, *pnext, *seq_string, *quality_string,*rname, *cigar;
    int flag, i;
    uint64_t pos;
    uint32_t mapq;
    long long int tlen;
    int min;
        int repeated = 0;
    rname = malloc(100 * sizeof(char));
    cigar = malloc(200 * sizeof(char));
    qname = malloc(100 * sizeof(char));
    rnext = malloc(100 * sizeof(char));
    pnext = malloc(100 * sizeof(char));
    seq_string = malloc(1000 * sizeof(char));
    quality_string = malloc(500 * sizeof(char));

    int stop = 0;
    while (1) {
            if (fgets(line, 1000, samFile) == NULL) 
                break;
            while (line[0] == '@'){
                if(fgets(line, 1000, samFile) == NULL){
                    stop = 1;
                    break;
                }
            }
            if(stop)
                break;
            readNum++;
            long long penalty;
            sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\t%d\t%lld\n",qname, &flag, rname, &pos,&mapq, cigar,rnext,pnext, &tlen,seq_string,quality_string,&min,&penalty);
        char *copy = (char *)malloc(strlen(qname) + 1);
        strcpy(copy, qname);
        char* tokens=strtok(qname, ":");
            tokens = strtok(NULL, ":");
            
            char *first,*second;
            first = strtok(tokens, "-");
            second = strtok(NULL, "-");
    
            //printf( "f: %s   sec:%s \n",first,second  );
            uint64_t start, end;
            start = strtoll(first,NULL,10);
            end = strtoll(second,NULL,10);
            //printf( "first: %" PRIu64 "   sec:%" PRIu64 " \n",start,end  );
            if(start-20 <= pos && pos <=end+20);
                //exact aligning
            else if(!strstr(cigar,"*")){
                badAlignedReads += 1;
                char str[80];
                strcpy (str,rname);
                strcat (str,":");
                char buffer [33];
                sprintf(buffer,"%d",pos);
                strcat (str,buffer);
                strcat (str,"-");
                sprintf(buffer,"%d",pos+99);
                strcat (str,buffer);
                fprintf(stdout,"%s \n",str);
                
                char * command = "samtools faidx /home/a.sharifi/Aryana_BS/Resources/Resources/Human/Index/hg19.fa ";
                char cmd_pointer[strlen(command) + 60];
                strcpy(cmd_pointer, command);
                strcat(cmd_pointer, str);
                //system(cmd_pointer);
                char * data=run(cmd_pointer);
                //1043_>chr8:6522050-6522149
                tokens=strtok(copy, ">");
                tokens = strtok(NULL, ">");
                char cmd_pointer2[strlen(command) + 60];
                strcpy(cmd_pointer2, command);
                strcat(cmd_pointer2, tokens);
                char * data2=run(cmd_pointer2);
                if(strcmp(data, data2)){
                    fprintf(stdout,"hooooooraaaaa \n");
                    repeated++;
                }

                
                fprintf(stdout,"aligned with: %d \t line :%s\t penalty:%ld\n",min+1,line,penalty);

		}
           // else
             //   fprintf(stdout,"%s\n",line);
    //readCigar(cigar, pos, seq_string, i,start ,end);

            //fprintf(stdout, "%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n \n",qname, flag, rname, pos,mapq, cigar,rnext,pnext, tlen,seq_string,quality_string);

           // fprintf(stdout,"\n cigar %s \n ",cigar);
            if(strstr(cigar,"*"))
                notAlignedReads++;
    }
    float accuracy = ((float)badAlignedReads/readNum)*100.0 ;
    float accuracy2 = ((float)notAlignedReads/readNum)*100.0 ;
    fprintf(stdout, "number of not aligned reads: %lld \n number of reads with wrong alignment : %lld \n total reads : %lld \n percentage of bad aligned reads :%10f \n percentage of not aligned reads :%10f \n",notAlignedReads, badAlignedReads, readNum , accuracy,accuracy2);
        
        fprintf(stdout, "reapeted: %d \n",repeated);
    fclose(samFile);
    }
}
