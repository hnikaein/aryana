//
//  methylation.h
//  aryana
//
//  Created by Maryam Rabiee on 9/20/14.
//  Copyright (c) 2014 Maryam Rabiee. All rights reserved.
//

#ifndef __aryana__methylation__
#define __aryana__methylation__

#include <iostream>
//int checkGAorCT(Line line);
//void computeMethylation();
//void convertGenome(Line &line);
char * convertCigar(char * cigar);
int readSamFile(FILE * samFile);
void setPointer(int pos, char * chr);
void removeLine();
#endif /* defined(__aryana__methylation__) */


