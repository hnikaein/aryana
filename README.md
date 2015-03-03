src
===

The source code of the Aryana next generation sequencing (NGS) aligner.

Installation
============

Download Aryana from github and run "make"

Running Aryana
==============

Create genome index and bin file from the genome

<code>
./aryana index [reference genome (fasta)]

./aryana fa2bin [reference genome (fasta)]
</code>

Alignment of single end reads

<code>
./aryana [-x Index] [-U input (fastq)]

Optional arguments:

	-p INT	Number of threads

	-P INT	Number of potential candidates to select for more precise alignment

	--seed INT	Fixed length for the seeds
</code>

Alignment of paired end reads

<code>
./aryana [-x Index] [-1 input1 (fastq)] [-2 input2 (fastq)] 

Optional arguments:

	--fr/--ff/--rf	Refers to orientation of the pairs. /forward-reverse/forward-forward/reverse-forward

	-I INT	Minimum distance between reads pair

	-X INT	Maximum distance between reads pair
</code>

Aryana Bisulfite Sequencing
===========================

Alignment of Bisulfite treated reads.

Installation
============
<code>
	make methyl
</code>

Running Aryana_BS
=================

Creating the converted reference genomes and the index files

<code>
	./prepare_genomes.sh <reference genome> <position of CpG islands file> <output folder>
</code>

Running an script for alignment of the reads to different converted reference genomes, selecting the best alignment, and computing the methylation ratios for each cytosine

<code>
	./aryana_bs <reference genome> <reference index folder> <CpG islands file> <input fastq file> <output file, without extensions> [ar="additional arguments to aryana"] [me="additional arguments to methyl_extract"]
</code>

Or the following script for the paired-end reads

<code>
    ./aryana_bsp <reference genome> <reference index folder> <CpG islands file> <input fastq file 1>  <input fastq file 2> <output file, without extensions> [ar="additional arguments to aryana"] [me="additional arguments to methyl_extract"]
</code>

** It is not needed to run "prepare_genomes.sh" for every alignment. You just need to run it once on your reference to be able to run "aryana_bs".

Example:

<code>

	./prepare_genomes hg19.fa cpg_island_hg19.txt BS_Genomes/

	./aryana_bs hg19.fa BS_Genomes/ cpg_islands_hg19.txt reads.fastq result ar="-p 10"

</code>

Analyzing output
================

The BisAnalyzer receives the SAM file of aligining Bis-Seq reads to the genome, produces a table containing the minimum edit distance between each read
and the aligned position, considering difference scenarios the read can be produced: Allowing either C->T or G->A conversion, checking PCR-read possibility, 
and considering the positive or negative strand.

If the reads are produces through BisSimul program, the BisAnalyzer can also analyze each read against the real genomic position the read has been simulated from.

The output is a tab-delimited table including several columns: the read name (ReadName), chromosome name and chromosomal position of the alignment according to SAM file (AlnChr and  AlnPos), the minimum edit distance - penalty - between the read and the aligned sequence based on the CIGAR sequence (AlnPen), and the same information for the real read location (RealChr, RealPos and RealPen).

<code>
    ./BisAnalyzer  -g <reference genome, mandatory> -i <alignment SAM file> -s (use this argument if the reads are generated via BisSimul)
</code>

