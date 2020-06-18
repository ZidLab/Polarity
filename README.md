#This is README file for polarity analysis

##For every input/output folder, remember to append '/' after the folder name

**In Inputs folder:**

Run script *Seq_analysis.sh* for sequencing data (The script *Seq_analysis_shihui.sh* is a test script which only takes 1 file in)

InputFolder: Folder contains .fasq files(or .sam files)

OuputFolder: SortedFiles/

**If the input files are in .fastq format (execute the first for loop in Seq-analysis.sh)**

1. *trimPolyA.py*

*parameters:*
input: .fastq files
output: trim/

*features:*
trim poly A tails

2. run bowtie (below are example file path)

bowtie -a --best --strata $HOME/Indexes/rDNA/ $HOME/trim/$file $HOME/Bowtie/rDNA-$file --un $HOME/Bowtie/coding-$file

bowtie -m 1 --best --strata $HOME/Indexes/genome $HOME/Bowtie/coding-/$file $dest$/Maps/$file

3. *FeatureCounts.py:* 

*parameters:*
input: Maps/, Chrom/
output: SortedFiles/XXX/

*features:*
Take in .sam files, converted in .csv format, stored in SortedFiles/

**If the input files are in .sam format(execute the second for loop in Seq-analysis.sh)**

1. *FeatureCounts_sam.py*

*parameters:*
input: .sam files, Chrom/
output: sortedFiles/

*features:*
Take in .sam files, converted in .csv format, stored in SortedFiles/

**Other Scripts**

*AlignSortCount.py*: align and sort the a.a. sequence

*MochiView.py*: count a.a. sequence

*countRibo.py*: count a.a. sequence

*Features.py*: similar with FeatureCounts, without format conversion

*polarity.py*: 

*parameters:*
input: sorted files in .csv format
output: Polarity/, in .csv format

*features:*
count polarity per gene
