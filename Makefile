CC=			gcc
CXX=		g++
CFLAGS=		-Wall -Wno-unused-function -Ofast
CXXFLAGS=	-Wall -Wno-unused-function -Ofast
OBJS=		QSufSort.o bwt_gen.o utils.o bwt.o bwtaln.o bwa2.o bwtgap.o sam.o hash.o smith.o aligner.o fa2bin.o \
			is.o bntseq.o bwtindex.o ksw.o stdaln.o simple_dp.o \
			bwaseqio.o bwase.o bwape.o kstring.o cs2nt.o \
			bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwt_lite.o \
			bwtsw2_chain.o bamlite.o bwtsw2_pair.o bwt2.o bwa.o
PROG=		aryana
INCLUDES=	
LIBS=		-lm -lz -lpthread
SUBDIRS=	. 
debug:		CFLAGS=-Wall -Wno-unused-function -DDEBUG -g3 -O0
debug:		CXXFLAGS=-Wall -Wno-unused-function -DDEBUG -g3 -O0

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

all:	$(PROG) convert_genomes align_bs methyl_extract read_simul SamAnalyzer fastaseq bwtcheck SamToNormWig
debug:	all

aryana:$(OBJS) aryana_main.o
		$(CC) $(CFLAGS) $(OBJS) aryana_main.o -o aryana $(LIBS)

convert_genomes:
		$(CXX) $(CXXFLAGS) convert_genomes.cpp -o convert_genomes

align_bs:
		$(CXX) $(CXXFLAGS)  align_bs.cpp -o align_bs 

methyl_extract:
		$(CXX) $(CXXFLAGS) methyl_extract.cpp -o methyl_extract

read_simul:	
		$(CXX) $(CXXFLAGS) read_simul.cpp -o read_simul

SamAnalyzer:
		$(CXX) $(CXXFLAGS) -std=c++11 SamAnalyzer.cpp -o SamAnalyzer

SamToNormWig:
		$(CXX) $(CXXFLAGS) SamToNormWig.cpp -o SamToNormWig

fastaseq:
		$(CXX) $(CXXFLAGS) fastaseq.cpp -o fastaseq

bwtcheck:
		$(CXX) $(CXXFLAGS) bwtcheck.cpp -o bwtcheck

QSufSort.o:QSufSort.h

bwt.o:bwt.h
bwt2.o:bwt.h
bwtaln.o:bwt.h bwtaln.h kseq.h bwa2.h
bwt1away.o:bwt.h bwtaln.h
bwt2fmv.o:bwt.h
bntseq.o:bntseq.h
bwtgap.o:bwtgap.h bwtaln.h bwt.h
aligner.o: aligner.h bwt.h hash.h smith.h
fa2bin.o: fa2bin.c
hash.o: hash.h
bwa2.o: bwa2.h sam.h aligner.h bwt.h
sam.o: sam.h bwt.h
smith.o: smith.h bwt.h
bwtsw2_core.o:bwtsw2.h bwt.h bwt_lite.h stdaln.h
bwtsw2_aux.o:bwtsw2.h bwt.h bwt_lite.h stdaln.h
bwtsw2_main.o:bwtsw2.h
bwa.o: bntseq.h bwa.h bwt.h ksw.h utils.h kstring.h malloc_wrap.h kvec.h
aryana_main.o: aryana_main.h aryana_args.h bwt.h bwtaln.h kseq.h bwa2.h

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a aryana align_bs methyl_extract convert_genomes read_simul SamAnalyzer bwtcheck fastaseq SamToNormWig
