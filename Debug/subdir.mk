################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../QSufSort.o \
../aligner.o \
../aryana_main.o \
../bamlite.o \
../bntseq.o \
../bwa2.o \
../bwape.o \
../bwase.o \
../bwaseqio.o \
../bwt.o \
../bwt_gen.o \
../bwt_lite.o \
../bwtaln.o \
../bwtgap.o \
../bwtindex.o \
../bwtio.o \
../bwtmisc.o \
../bwtsw2_aux.o \
../bwtsw2_chain.o \
../bwtsw2_core.o \
../bwtsw2_main.o \
../bwtsw2_pair.o \
../cs2nt.o \
../fa2bin.o \
../fastmap.o \
../hash.o \
../is.o \
../kstring.o \
../ksw.o \
../main.o \
../sam.o \
../simple_dp.o \
../smith.o \
../stdaln.o \
../utils.o 

C_SRCS += \
../QSufSort.c \
../aligner.c \
../aryana_args.c \
../aryana_main.c \
../bamlite.c \
../bntseq.c \
../bwa2.c \
../bwape.c \
../bwase.c \
../bwaseqio.c \
../bwt.c \
../bwt_gen.c \
../bwt_lite.c \
../bwtaln.c \
../bwtgap.c \
../bwtindex.c \
../bwtio.c \
../bwtmisc.c \
../bwtsw2_aux.c \
../bwtsw2_chain.c \
../bwtsw2_core.c \
../bwtsw2_main.c \
../bwtsw2_pair.c \
../cigar.c \
../cigars.c \
../cs2nt.c \
../fa2bin.c \
../fastmap.c \
../hash.c \
../is.c \
../kbwt.c \
../kstring.c \
../ksw.c \
../main.c \
../sam.c \
../simple_dp.c \
../smith.c \
../stdaln.c \
../utils.c 

OBJS += \
./QSufSort.o \
./aligner.o \
./aryana_args.o \
./aryana_main.o \
./bamlite.o \
./bntseq.o \
./bwa2.o \
./bwape.o \
./bwase.o \
./bwaseqio.o \
./bwt.o \
./bwt_gen.o \
./bwt_lite.o \
./bwtaln.o \
./bwtgap.o \
./bwtindex.o \
./bwtio.o \
./bwtmisc.o \
./bwtsw2_aux.o \
./bwtsw2_chain.o \
./bwtsw2_core.o \
./bwtsw2_main.o \
./bwtsw2_pair.o \
./cigar.o \
./cigars.o \
./cs2nt.o \
./fa2bin.o \
./fastmap.o \
./hash.o \
./is.o \
./kbwt.o \
./kstring.o \
./ksw.o \
./main.o \
./sam.o \
./simple_dp.o \
./smith.o \
./stdaln.o \
./utils.o 

C_DEPS += \
./QSufSort.d \
./aligner.d \
./aryana_args.d \
./aryana_main.d \
./bamlite.d \
./bntseq.d \
./bwa2.d \
./bwape.d \
./bwase.d \
./bwaseqio.d \
./bwt.d \
./bwt_gen.d \
./bwt_lite.d \
./bwtaln.d \
./bwtgap.d \
./bwtindex.d \
./bwtio.d \
./bwtmisc.d \
./bwtsw2_aux.d \
./bwtsw2_chain.d \
./bwtsw2_core.d \
./bwtsw2_main.d \
./bwtsw2_pair.d \
./cigar.d \
./cigars.d \
./cs2nt.d \
./fa2bin.d \
./fastmap.d \
./hash.d \
./is.d \
./kbwt.d \
./kstring.d \
./ksw.d \
./main.d \
./sam.d \
./simple_dp.d \
./smith.d \
./stdaln.d \
./utils.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


