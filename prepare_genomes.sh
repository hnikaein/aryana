#! /bin/bash

if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters: <reference genome> <position of CpG islands file> <output folder>"
	exit
fi
DIR=$(dirname $(readlink -f $0));
mkdir "$3"
echo "Start of convert_genomes"
"$DIR/convert_genomes" -g "$1" -a "$2" -o "$3"
echo "End of convert_genomes"
echo "Start of aryana"
"$DIR/aryana" index "$3/originalGenome.fa" &
"$DIR/aryana" index "$3/BisulfiteGenomeIslandConsideredCT.fa" &
#"$DIR/aryana" index "$3/BisulfiteGenomeContextConsideredCT.fa" &
"$DIR/aryana" index "$3/BisulfiteGenomeCompleteCT.fa" &
"$DIR/aryana" index "$3/BisulfiteGenomeIslandConsideredGA.fa" &
#"$DIR/aryana" index "$3/BisulfiteGenomeContextConsideredGA.fa" &
"$DIR/aryana" index "$3/BisulfiteGenomeCompleteGA.fa" &
"$DIR/aryana" fa2bin "$3/originalGenome.fa" &
"$DIR/aryana" fa2bin "$3/BisulfiteGenomeIslandConsideredCT.fa" &
#"$DIR/aryana" fa2bin "$3/BisulfiteGenomeContextConsideredCT.fa" &
"$DIR/aryana" fa2bin "$3/BisulfiteGenomeCompleteCT.fa" &
"$DIR/aryana" fa2bin "$3/BisulfiteGenomeIslandConsideredGA.fa" &
#"$DIR/aryana" fa2bin "$3/BisulfiteGenomeContextConsideredGA.fa" &
"$DIR/aryana" fa2bin "$3/BisulfiteGenomeCompleteGA.fa" &
echo "aryana is running in the background"
exit
