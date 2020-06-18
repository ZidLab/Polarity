#!/bin/sh
dest=$HOME'/Zid-Lab/Analysis/'
input=$dest'Zid_old_data/'
cur=$dest'Inputs/'
# Set file names and paths for use in bowtie and data extraction
mapDir='Maps/'
bowDir='Bowtie/'
# Set index basenames
rdnaDir='Indexes/rDNA'
geneDir='Indexes/genome'
ecoliDir='Indexes/e_coli'
# Set bowtie output filenames
rDNA=$dest$bowDir'rDNA-'
coding=$dest$bowDir'coding-'
map='aligned-'
dchrom='Chrom/'
trimf='trim/'
sortedGeneDir='SortedFiles/'
trimpy='trimPolyA.py'
fct='FeatureCounts.py'
# Create bowtie argument filepaths
arg1=$dest$rdnaDir
arg4=$dest$geneDir


for f in $input/*.fastq; do 
	file=`basename $f`;
	echo $file;
	python $cur$trimpy $dest$trimf$file $f;
	#echo $dest$mapDir$file;	
	export PATH=$PATH:/opt/biotools/bowtie/bin/
        bowtie -a --best --strata $arg1 $dest$trimf$file $rDNA$file --un $coding$file
        bowtie -m 1 --best --strata $arg4 $coding$file $dest$mapDir$file
	python $cur$fct $dest$sortedGeneDir $dest$mapDir$file $dest$dchrom 
done
