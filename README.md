### ProTree_pipe.sh

ProTree is a simple (and not very clever!) automated pipeline for generating phylogenetic trees based on protein data sets.

It should suit quick and simple phylogenetic inferences without much input, however if there are any desired modifications,
you may do so within this script at their respective command.

The script should be run within a directory containing desired protein fastas (labelled with '.faa' suffix).

### quick start - with all protein sets in current directory, ending in .faa
>  ./Protree_pipe.sh num_threads

### dependencies

conda install orthofinder

conda install -c bioconda gblocks

conda install -c bioconda mafft

faSomeRecords  -  https://github.com/santiagosnchez/faSomeRecords

fasta_to_phylip.py  -  https://github.com/audy/bioinformatics-hacks/blob/master/bin/fasta-to-phylip

prottest3 - https://github.com/ddarriba/prottest3
  - prottest3 also requires java! (tested on "openjdk 11.0.1 2018-10-16 LTS")

conda install -c genomedk raxml-ng 

### how to run

0. activate the conda environment and install the necessary software

e.g   source activate "env"

1. copy all protein sets to compare into a single directory 
  - it is helpful to rename each with an organism/specific identifier to rename the proteins within the file. 
  - The name must also be short, where Gblocks will complain of protein fasta headers being too long and disregarding the
  - alignnment.

###	- ALL FASTAS MUST END IN " .faa ".

2. execute via : ' ProTree_pipe.sh "threads to use" ' 

### overview of pipeline
 
1. single copy orthologue inference from Orthofinder
2. alignment of single copy orthologues with mafft
3. trimming of poorly aligned regions using Gblocks and alignments concatenated
4. conversion to phyllip format and model testing using prottest3 (Chooses BIC model)
5. raxml-ng with 10 randomized parsimony starting trees and 1000 bootstraps

### to do!

1. fix the tar'ing of final working directory.
2. allow continuing from specific steps, rather than re-do pipeline all over again.
