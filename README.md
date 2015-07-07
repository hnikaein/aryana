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
./aryana {arguments}
	-x, --index 		<reference genome index>
	-i, --input 		<reads file (fastq format)>

Optional arguments:
	-o, --output		<alignment file (SAM format)>
	-t, --threads 		<threads number> 
	-c, --candidates	<num of alignment candidates, default=10>
	-s, --seed		<fixed length of seed sequence, default=dynamic>
	-D 			<debug info level, default=0>
	-O			(keep the order of reads in output the same as input)
    -e          <number of exact matches to select>

</code>

Alignment of paired end reads

<code>
./aryana {arguments}
	-x, --index 		<reference genome index>
	-1, --first 		<reads file 1 (fastq format)>
	-2, --second 		<reads file 2 (fastq format)>

Optional arguments:
	--fr, --ff, --rf	(orientation of paired ends, --fr:forward-reverse, --ff=forward-forward, --rf=reverse-forward)
	-m, --min 		<min distance between pair reads>
	-M, --max  		<max distance between pair reads>
	-d, --no-discordant 	(do not print discordants reads)
</code>

Making Index for Bisulfite Sequencing
=====================================


Creating the converted reference genomes and the index files

<code>
./prepare_genomes.sh <reference genome> <position of CpG islands file> <output folder>
</code>

The first 3 columns of the CpG-islands file should contain chromosome name, start position and end position of each CpG island, respectively. 
The chromosome names should match the reference genome. No header line should be provided. 

Aligning Bisulfite Sequencing Reads Using Aryana
================================================

For this purpose you should replace -x or --index argument of aryana with -b or --bisulfite. All the other arguments of aryana apply here.
<code>
./aryana {arguments}
	-b, --bisulfite 	<bisulfite reference genome index, the output folder given to prepare_genomes.sh>
</code>

For simplicity, there is a script for alignment of the reads to different converted reference genomes, selecting the best alignment, and computing the methylation ratios for each cytosine

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

The SamAnalyzer receives the SAM file of aligining either normal DNA or Bis-Seq reads to the genome, produces a table containing the minimum edit distance between each read
and the aligned position. For Bis-Seq reads, it considers difference scenarios the read can be produced: Allowing either C->T or G->A conversion, checking PCR-read possibility, 
and considering the positive or negative strand.

If the reads are produces through BisSimul program, or the names of the reads include the true genomic location, in the same way that BisSimul does, the SamAnalyzer can also analyze 
each read against the real genomic position the read has been simulated from.

The output is a tab-delimited table including several columns: the read name (ReadName), chromosome name and chromosomal position of the alignment according to SAM file (AlnChr and  AlnPos), the minimum edit distance - penalty - between the read and the aligned sequence based on the CIGAR sequence (AlnPen), and the same information for the real read location (RealChr, RealPos and RealPen) if the -s argument is used.

If -i and -o arguments are not used, the standard input and standard output will be used.

<code>
./SamAnalyzer {arguments}
 	-g 		<reference genome, mandatory> 
	-i 		<alignment SAM file> 
	-o 		<output tabular file>  
    	-S 		(print SAM, CIGAR and alignment sequences in output)
    	-b 		(input is Bis-Seq, default: normal DNA-seq) 
    	-s 		(use this argument if the reads are generated via BisSimul, or the read names include the true position in a similar way)
</code>

Here is a description about the SamAnalyzer results:


ReadName: name of read, as appears in FASTQ and/or SAM file
Aligned: shows wether the read is aligned to the reference genome (1) or not (0)
AlnChr: The chromosome to which the read is aligned, or NA if the read is not aligned
AlnPos: The position in the chromosome to which the read is aligned
AlnMismatch: The number of mismatches between the read sequence and alignment position
AlnGapOpen: The number of gap intervals between the read and the alignment position. 
                       Any number of consecutive insertions in read or reference are considered as one gap interval.
AlnGapExt: The total number of nucleotides extending the gaps. For an interval of N insertions, the gap extension number is N-1. 


The following columns are reported if the reads are simulated using the read_simul of aryana:
RealChr: The true chromosome from which the read is simulated
RealPos: The true position in the chromosome from which the read is simulated
RealMismatch: The number of mismatches between the true genomic sequence and the simulated read, possibly due to simulated SNPs or sequencing error

The following columns are added if -s argument is used:
SamSeq: The nucleotide sequence reported in SAM file
CIGAR: The CIGAR sequence reported in SAM file, that tells how 
AlnSeq: The genomic sequence to which the read is aligned (i.e. located in AlnChr,AlnPos)


The following columns are added if -s argument is used and the reads are simulated using read_simul:
RealSeq: The genomic sequence from which the read is simulated

To find out how good each read is aligned we might use Aligned, AlnMismatch, AlnGapOpen and AlnGapExt columns. We can define a penalty, for instance Penalty = 5 * (AlnMismatch + AlnGapOpen) + 3 * AlnGapExt. This will make higher penalty for mismatches and start of the gaps, but lower penalty for gap extensions, since it's more likely for an opened gap to be extended.

A simpler approach would be just to consider AlnMismatch.
