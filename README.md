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

** It is not needed to run "prepare_genomes.sh" for every alignment. You just need to run it once on your reference to be able to run "aryana_bs".

Example:

<code>

	./prepare_genomes hg19.fa cpg_island_hg19.txt BS_Genomes/

	./aryana_bs hg19.fa BS_Genomes/ cpg_islands_hg19.txt reads.fastq result ar="-p 10"

</code>
