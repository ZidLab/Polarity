#!/bin/sh
filename=$1
echo $filename
dest=$HOME'/Zid-Lab/Analysis/'
input=$dest'Anna/longmers/'
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
fct='FeatureCounts_sam.py'
# Create bowtie argument filepaths
arg1=$dest$rdnaDir
arg2=$rDNA$filename
arg3=$coding$filename
arg4=$dest$geneDir
#arg5=$dest$mapDir$filename
arg5=$input$filename

#python $cur$trimpy $dest$trimf$filename $input$filename
#export PATH=$PATH:/opt/biotools/bowtie/bin/
#module load bowtie
#cd /opt/biotools/bowtie/bin/
#bowtie -a --best --strata $arg1 $dest$trimf$filename $arg2 --un $arg3
#bowtie -m 1 --best --strata $arg4 $arg3 $arg5
#cd $cur
python $cur$fct $dest$sortedGeneDir $arg5 $dest$dchrom
