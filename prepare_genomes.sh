#! /bin/bash

if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters: <reference genome> <position of CpG islands file> <output folder>"
	exit
fi

echo "Start of convert_genomes"
./convert_genomes -g $1 -a $2 -o $3
echo "End of convert_genomes"
echo "Start of bwa"
./bwa index $3originalGenome.fa &
./bwa index $3BisulfiteGenomeIslandConsideredCT.fa &
./bwa index $3BisulfiteGenomeContextConsideredCT.fa &
./bwa index $3BisulfiteGenomeCompleteCT.fa &
./bwa index $3BisulfiteGenomeIslandConsideredGA.fa &
./bwa index $3BisulfiteGenomeContextConsideredGA.fa &
./bwa index $3BisulfiteGenomeCompleteGA.fa &
./bwa fa2bin $3originalGenome.fa &
./bwa fa2bin $3BisulfiteGenomeIslandConsideredCT.fa &
./bwa fa2bin $3BisulfiteGenomeContextConsideredCT.fa &
./bwa fa2bin $3BisulfiteGenomeCompleteCT.fa &
./bwa fa2bin $3BisulfiteGenomeIslandConsideredGA.fa &
./bwa fa2bin $3BisulfiteGenomeContextConsideredGA.fa &
./bwa fa2bin $3BisulfiteGenomeCompleteGA.fa &
echo "bwa is running in the background"
exit
