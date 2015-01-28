#! /bin/bash

if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters: <reference genome> <position of CpG islands file> <output folder>"
	exit
fi
DIR=$(dirname $(readlink -f $0));
echo "Start of convert_genomes"
"$DIR/convert_genomes" -g "$1" -a "$2" -o "$3"
echo "End of convert_genomes"
echo "Start of bwa"
"$DIR/bwa" index "$3originalGenome.fa" &
"$DIR/bwa" index "$3BisulfiteGenomeIslandConsideredCT.fa" &
#"$DIR/bwa" index "$3BisulfiteGenomeContextConsideredCT.fa" &
"$DIR/bwa" index "$3BisulfiteGenomeCompleteCT.fa" &
"$DIR/bwa" index "$3BisulfiteGenomeIslandConsideredGA.fa" &
#"$DIR/bwa" index "$3BisulfiteGenomeContextConsideredGA.fa" &
"$DIR/bwa" index "$3BisulfiteGenomeCompleteGA.fa" &
"$DIR/bwa" fa2bin "$3originalGenome.fa" &
"$DIR/bwa" fa2bin "$3BisulfiteGenomeIslandConsideredCT.fa" &
#"$DIR/bwa" fa2bin "$3BisulfiteGenomeContextConsideredCT.fa" &
"$DIR/bwa" fa2bin "$3BisulfiteGenomeCompleteCT.fa" &
"$DIR/bwa" fa2bin "$3BisulfiteGenomeIslandConsideredGA.fa" &
#"$DIR/bwa" fa2bin "$3BisulfiteGenomeContextConsideredGA.fa" &
"$DIR/bwa" fa2bin "$3BisulfiteGenomeCompleteGA.fa" &
echo "bwa is running in the background"
exit
