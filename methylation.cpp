//
//  methylation.cpp
//  aryana
//
//  Created by Maryam Rabiee on 9/20/14.
//  Copyright (c) 2014 Maryam Rabiee. All rights reserved.
//


#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <vector>
#include <math.h>

#include <iostream>
#include <fstream>
#include "methylation.h"

using namespace std;



#define maxGenomeSize 1e5
#define maxChromosomeNum 1000
#define READ_LENGHT 100
#define READ_SAM_BY 5
unsigned long gs;
char chromName[maxChromosomeNum][100];
int chromNum;

char * reference;
uint32_t reference_size;
uint32_t reference_reminder;
FILE *samFile;
int firstPointer = 0 , secondPointer=0;

struct Line{
    char line[1000];
    char chr[30];
    uint64_t start , pos;
   // char *qname, *seq_string,*rname, *cigar;
    char cigar[300];
    char cigar2[500];
    char rname[200];
    string seq_string;
    char strand;
    
};
//struct line lines[1000];
std::vector<Line> lines;
struct Chrom
{
    char * chrName;
    long chrStart;
    long chrNum;
};

struct Cytosine
{
    long pos;
    char chr[30];
    int methylated;
    int unmethylated ;
    
};
std::vector<Cytosine> cytosines;
struct Chrom chrom[maxChromosomeNum] ;


int ChromIndex(char * chr) {
	int i = 0;
	for(i; i < maxChromosomeNum; i++){
		if(chrom[i].chrName == NULL)
			break;
		if(strcmp(chrom[i].chrName, chr) == 0)
			return i;
	}
	return -1;
    
}

int checkGAorCT(Line line){
    int GA =0;
    int CT = 0 , i = 0;
    int refPos , readPos = 0;
    int ref = chrom[ChromIndex(line.chr)].chrStart +line.pos;
    cout<<"GA1 "<<GA<<endl;
    cout << "cigar2"<<line.cigar2<<endl;
    while (line.cigar2[i] != '\0') {
        //cout<<"GA2    "<<GA;


        if(line.cigar2[i] == 'i'){
            readPos++;
            i++;
            continue;
        }
        if(line.cigar2[i] == 'd'){
            i++;
            ref++;
            continue;
        }
        if(line.cigar2[i] == 'm'){
            
            cout<<"  "<<line.seq_string[readPos]<<i<<reference[ref]<<" ";
            if(line.seq_string[readPos] == 'A'){
                if(reference[ref] == 'G')
                    GA++;
            }
            if(line.seq_string[readPos] == 'T'){
                if(reference[ref] == 'C')
                    CT++;
            }
            i++;
            ref++;
            readPos++;
        }
    }
    cout<<"cigar2[i]  :"<<line.cigar2[i]<<" "<<i<<endl;
    cout<<"GA  :"<<GA<<endl;
    cout<<"CT   :"<<CT<<endl;


    if(GA > CT + 5)
        return 16;
    if(GA < CT - 5)
        return 4;
    else
        return -1;
    
}

int checkGAorCT2(Line line){
    int GA =0;
    int CT = 0 , i = 0;
    int refPos , readPos = 0;
    int ref = chrom[ChromIndex(line.chr)].chrStart +line.pos;
    
    cout<<"GA1 "<<ref<<endl;
    cout << "cigar2"<<line.cigar2<<endl;
    while (i < line.seq_string.size()) {
            cout<<"  "<<line.seq_string[i]<<i<<reference[ref+i]<<" ";
            if(line.seq_string[i] == 'A'){
                if(reference[ref+i] == 'G')
                    GA++;
            }
            if(line.seq_string[i] == 'T'){
                if(reference[ref+i] == 'C')
                    CT++;
            }
            i++;

    }
    cout<<"cigar2[i]  :"<<line.cigar2[i]<<" "<<i<<endl;
    cout<<"GA  :"<<GA<<endl;
    cout<<"CT   :"<<CT<<endl;
    
    
    if(GA > CT)
        return 16;
    else if(GA < CT)//////////////
        return 4;
    else
        return -1;
    
}


void computeMethylation(){
    

    cout << "heyyyy1"<<endl;
    if(lines.size() == 0){
        return;
    }
    for(int k=0;k<cytosines.size();k++)
        cout <<"cytosin  "<< cytosines[k].pos<<"   "<<cytosines[k].methylated<<"   "<<cytosines[k].chr<<endl;
    
    int lastchecked = lines[0].pos;

    if(cytosines.size()!=0){
        cout << "heyyyy888   "<<strcmp(cytosines[cytosines.size()-1].chr,lines[0].chr)<<endl;
        cout << lines[0].pos <<"    "<<lines[0].seq_string <<"   "<<lines[0].chr<<"  "<<cytosines[cytosines.size()-1].chr<<endl;
        lastchecked= cytosines[cytosines.size()-1].pos;
        if (lines[0].pos > lastchecked || strcmp(cytosines[cytosines.size()-1].chr,lines[0].chr)) {
                cout << "heyyyy77"<<lines[0].pos<<" "<<lastchecked<<cytosines[cytosines.size()-1].chr<<endl;
            lastchecked = lines[0].pos;
        }
    }

    cout << "lastcheck  :"<<lastchecked<<endl;
    int refPos = chrom[ChromIndex(lines[0].chr)].chrStart;
    cout << "size :"<<lines[0].seq_string<<"   "<< lines[0].chr<<endl;
    
    for(int i = lastchecked + 1; i< (lines[0].seq_string.size()+lines[0].pos) ;i++){
        if(reference[refPos + i] == 'C'){
            //cout<<refPos + i<<"     "<<reference[refPos + i]<<"i:   "<<i<<endl;
           // cout << "rrrr"<<lines[0].pos <<"    "<<lines[0].seq_string <<"   "<<lines[0].chr<<endl;
            setPointer(i , lines[0].chr, secondPointer);
            cytosines.push_back(Cytosine());
            strcpy(cytosines[cytosines.size()-1].chr, lines[0].chr);
            cout<<lines[0].chr<<endl;
            cytosines[cytosines.size()-1].pos = i + lines[0].pos ;
            cytosines[cytosines.size()-1].methylated = 0;
            cytosines[cytosines.size()-1].unmethylated = 0;
            cout<<"scn "<<secondPointer<<endl;
            
            
            for(int j=0 ; j<  secondPointer ; j++){
                //cout<<"seq_string[ i ].pos]  "<<lines[j].seq_string[i - lines[j].pos+1]<<endl;
                if ((lines[0].pos + i - lines[j].pos) >= 0 && lines[j].seq_string[i - lines[j].pos+1] == 'C') {                  
                    cytosines[cytosines.size()-1].methylated++;                   
                }
                else if(lines[j].seq_string[lines[0].pos + i - lines[j].pos] == 'T' && (lines[0].pos + i - lines[j].pos) >= 0){
                    cytosines[cytosines.size()-1].unmethylated++;
                }
            }
            for(int k=0;k<cytosines.size();k++)
                cout <<"cytosin2  "<< cytosines[k].pos<<"   "<<cytosines[k].methylated<<"   "<<cytosines[k].chr<<endl;
            
        }
    }
    cout<<"sizeeeeee:       "<<lines.size()<<endl;
    lines.erase(lines.begin());
    cout<<"sizeeeeee:       "<<lines.size()<<endl;
    computeMethylation();
}

void convertRead(Line &line){
    //cout <<"before :  "<<line.seq_string<<endl;
    int c=0;
    int j = 0;
    for(int i=0;i< line.seq_string.size();i++){
        if (line.cigar2[i] == 'i') {
            line.seq_string.erase(j,1);
            //cout <<"after:  "<<line.seq_string<<endl;
        }
        else if (line.cigar2[i] == 'd') {
            line.seq_string.insert(j,"M");
            j++;
        }
        else
            j++;
        
    }
    cout <<"after:  "<<line.seq_string<<endl;
}

void reverseRead(Line &line){////
    int i,j=0 ;
    char * copy = new char[lines[0].seq_string.size()+1];
//    strcpy(copy , line.seq_string);
    line.seq_string.copy(copy, line.seq_string.size(),0);
    for (i = lines[0].seq_string.size(); i > 0 ; i--){
        if(copy[i] == 'A'){
            line.seq_string[j]='T';
        }
        if(copy[i] == 'T'){
            line.seq_string[j]='A';
        }
        if(copy[i] == 'C'){
            line.seq_string[j]='G';
        }
        if(copy[i] == 'G'){
            line.seq_string[j]='C';
        }
        j++;
    }    
}

char * convertCigar(char * cigar){
    char * cigar2 = new char[500];
    int pos = 0;
    int value = 0, j , index =0;;
    long read_index = 0;
    char alignType;
    cout << "cigar in convert"<<cigar<<endl;
    while (1) {
		if (!isdigit(cigar[pos])) {
			if (value > 0) {
                
				if (cigar[pos] == 'm') 
					for (j = 0; j < value; j++) {
                        cigar2[index] = 'm';
                        index++;
                    }
                else if (cigar[pos] == 'd')
                    for (j = 0; j < value; j++) {
                        cigar2[index] = 'd';
                        index++;
                    }
                else if (cigar[pos] == 'i')
                    for (j = 0; j < value; j++) {
                        cigar2[index] = 'i';
                        index++;
                    }
                
            
            }
            if (cigar[pos] == 0)
                break;
            value = 0;
        } else {
            value = value * 10 + cigar[pos] - '0';
        }
        pos++;
    }
    cout << "cigar2"<<cigar2<<endl;
    return cigar2;
}

int readSamFile(FILE * samFile){
    int index = 0;
    char line[1000];
    int header = 1;
    char *qname, *rnext, *pnext, *seq_string, *quality_string,*rname, *cigar;
    int flag, i;
    uint64_t pos;
    uint32_t mapq;
    long long int tlen;
    
    rname = new char[100];
    cigar = new char[200];
    qname = new char[100];
    rnext = new char[100];
    pnext = new char[100];
    seq_string = new char[1000];
    quality_string = new char[500];
    
    int count_to=0;
    int stop = 0;
    while (count_to< READ_SAM_BY) {
        if (fgets(line, 1000, samFile) == NULL){
            stop = -1;
            break;
        }
            
        while (line[0] == '@'){
            if(fgets(line, 1000, samFile) == NULL){
                stop = 1;
                break;
            }
        }
        cout<<line<<endl;
        if(stop)
            break;
        struct Line temp;
        
        sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, temp.chr, &(temp.pos) ,&mapq, temp.cigar,rnext,pnext, &tlen,seq_string,quality_string);
        string str(seq_string);
        temp.seq_string = str;
        
        
        char * cigar2 = new char[500];
        cout<<line<<endl;
        
        if(cigar[0] != '*'){
            cigar2 = convertCigar(temp.cigar);
            
            strcpy(temp.cigar2, cigar2);
            cout << "hey  "<<temp.cigar2<<endl;
        }
        cout<<line<<"      2"<<endl;

        
        convertRead(temp);
        cout<<"heyyyyyyyyyyyyyyyyy              "<<temp.seq_string<<"    "<<temp.chr<<endl;
    
        int result = checkGAorCT2(temp);
        cout << result<<endl;
        if (result == -1) {
            continue;//////////////////////////////////////////
        }
        if(result == 4)
            temp.strand = '+';
        else{
            temp.strand = '-';
            //reverseRead(temp);/////////////////////////////////////////////////////////////////
        }
        cout<<line<<"      3"<<endl;
        count_to++;
        //1000000_chr21:30778456-30778555
        //968_>chr2:24442257-24442356

        index++;
        if (count_to==READ_SAM_BY) {
            //computeMethylation(lines[count_to_tousand].start,lines[count_to_tousand].end , lines[count_to_tousand].seq_string);
            count_to=0;
        }
        lines.push_back(temp);
        cout<<lines[0].pos<<"      3"<<endl;
        
        
    }//while
    if(stop == -1)
        return -1;
    return 1;

}

void readSam_Compute(){
    
    int result = readSamFile(samFile);
    cout<<result <<"llll        "<<lines.size()<<endl;
    while (result == 1) {
        computeMethylation();
        result = readSamFile(samFile);
    }
    computeMethylation();
}


void setPointer(int pos, char * chr , int lineStart){
    cout<<"in set pointer ,pos:"<<pos<<"start:"<<lineStart<<endl;
    int i = lineStart;
    while(i < lines.size() && lines[i].pos <= pos && i < lineStart + READ_SAM_BY && !strcmp(lines[i].chr,chr)){
        if (lines[i].seq_string.empty()) {
            break;
        }
        secondPointer++;
        i++;
    }
    if (i - lineStart == READ_SAM_BY && lines[i].pos <= pos && !strcmp(lines[i].chr,chr)) {
        readSamFile(samFile);
        setPointer(pos, chr, lineStart+READ_SAM_BY);
        
    }
    if (secondPointer>lines.size()) {
        secondPointer = lines.size();
    }
}


int ReadGenome(char * genomeFile) {
    fprintf(stderr, "Allocating memory...\n");
    struct stat file_info;
    if (stat(genomeFile, &file_info) == -1) {
        fprintf(stderr,
                "Could not get the information of file %s\nplease make sure the file exists\n",
                genomeFile);
        return -1;
    }
    FILE *fp;
    fp = fopen(genomeFile, "r");
    off_t file_size_bytes = file_info.st_size;
    reference_size = ceil(((double) file_size_bytes) / (double) (sizeof(char)));
    reference = (char *) malloc(reference_size * sizeof(char));
    memset(reference, 0, reference_size * sizeof(char));
	gs = 0;
	chromNum = 0;
    //fprintf(stderr, "Reading genome...\n");
    char fLine[10000];
    while (! feof(fp)) {
        // fprintf(stderr, "Reading genome1...\n");
		int n = fscanf(fp, "%s\n", fLine);
		if (n == EOF) break;
		n = strlen(fLine);
        // fprintf(stderr, "Reading genome2 :n %d ...\n",n);
		if (fLine[0] == '>') {
			//chromPos[chromNum++] = gs;
            chrom[chromNum++].chrStart = gs;
            char * temp = fLine;
            temp++;
			int lenght = strlen(temp);
			chrom[chromNum-1].chrName = (char *) malloc(lenght * sizeof(char));
            strcpy(chrom[chromNum-1].chrName, temp);
            //fprintf(stderr, "Reading genome3...\n");
			if (chromNum >= 1)
			    fprintf(stderr, "chrName: %s, chrStart: %ld\n",chrom[chromNum-1].chrName, chrom[chromNum-1].chrStart);
			fprintf(stderr, "%s\n",fLine);
			//strcpy(chromName[chromNum - 1], fLine);
		} else {
            //            ToUpper(fLine);
			memcpy(reference+gs, fLine, n);
			gs += n;
		}
	}
    //chromPos[chromNum] = gs;
    // fprintf(stderr, "Lenght: %ld\n",chromPos[chromNum] - chromPos[chromNum - 1]);
    fclose(fp);
}




long methylated = 0;
//void computeMethylation(){
//    int j=0;
//    int lastchecked ;
//    
//    while(j < 1000){
//        
//        
//        int pos = lines[j].pos;
//        char *chrom = lines[j].chr;
//        pch=strchr(lines[j].seq_string,'C');
//        int refPos2 = lines[j].pos + chrom[ChromIndex(lines[j].rname)].chrStart;
//        int refPos = lines[j].pos + chrom[ChromIndex(lines[j].rname)].chrStart +
//        pch - lines[j].seq_string + 1;
//        lastchecked = j;
//        if(reference[refPos] == 'C' && reference[refPos+1] == 'G'){
//            methylated ++;
//        }
//        else{
//            outOfcontext++;
//            pch=strchr(lines[j].seq_string,'C');
//            if(pch == NULL){
//                j = lastchecked + 1;
//                continue;
//            }
//            
//        }
//        while (!strcmp(line[j].chr, chrom) && (pos + READ_LENGHT) >= lines[j].pos && j < 1000) {
//            char * pch;
//            pch2=strchr(lines[j].seq_string,'C');               ////////////handle lowercase too
//            if(pch2 =pch)
//                methylfolan++;
//            else
//                j++;
//            
//            //            int refPos2 = lines[j].pos + chrom[ChromIndex(lines[j].rname)].chrStart;
//            //            int loc = pch - lines[j].seq_string+1;
//            
//            
//        }
//        
//    }
//    
//}

//void computeMethylation2(char *chr){
//    int i;
//    
//    int index = ChromIndex(chr);
//    int start = chrom[index].chrStart;
//    int end = chrom[index].
//}

int main(int argc, char *argv[]) {
    char *referenceName, *annotationFile;
    
    if ( argc != 3 ){
        /* We print argv[0] assuming it is the program name */
        printf( "usage: %s filename", argv[0] );
    }
    else
    {
        ReadGenome(argv[1]);
        for(int k = 0 ; k<200;k++)
            cout<<reference[k]<<" ";
        samFile = fopen( argv[2], "r" );
        if ( samFile == 0 )
            printf( "Could not open file\n" );
        
        //readSamFile(samFile);
        readSam_Compute();
//        fclose(samFile);
//        int i=0;
//        while (i<100) {
//            cerr << lines[i].pos << "\t"<<lines[i].cigar<<endl;
//            i++;
//        }
        
        for(int k=0;k<cytosines.size();k++)
            cout <<"cytosin  "<< cytosines[k].pos<<"   "<<cytosines[k].methylated<<"   "<<cytosines[k].chr<<endl;
    }
    
}//main

