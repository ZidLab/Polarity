#!/bin/sh
#PBS -q hotel
#PBS -N run_brian_bk6
#PBS -l nodes=1:ppn=6,walltime=10:00:00
#PBS -o run_brian_bk6.o
#PBS -e run_brian_bk6.e
#PBS -V
#PBS -M shz334@ucsd.edu
#PBS -m abe
dest=$HOME'/Zid-Lab/Analysis/'
input=$dest'Kyle/Kyle_3_7_17/'
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
sortedGeneDir='SortedFiles/Kyle_data_sorted/Kyle_3_7_17_sorted/'
trimpy='trimPolyA.py'
fct='FeatureCounts.py'
fctsam='FeatureCounts_sam.py'
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
for f in $input/*.sam; do
	file=`basename $f`;
        echo $file;
	python $cur$fctsam $dest$sortedGeneDir $f $dest$dchrom
done
