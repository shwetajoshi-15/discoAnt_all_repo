# discoAnt
The repository is under-construction.

1. git clone repository
2. cd discoAnt
3. bash discoAnt_setup.sh
4. conda activate discoAnt.env
5. cd /path/to/folder/discoAnt/programs/cDNA_Cupcake
6. python setup.py build
7. python setup.py install
8. Submit discoAnt.sh in the appropriate submission script

# What does discoAnt do?

## This script performs the following steps
## 1. Quanlity control for the fasta files 
## 1.a. Generates a file with no. of reads
## 1.b. average read length
## 1.c. read length in each barcode/sample 
## 2. Align the sample fasta files to reference genome
## 3. Merge the primary alignments
## 4. Generate transcript files (GTF) for the merged alignments
## 4. Classify the transcripts based on known annotations 
## 6. Create a metagene based on the transcripts
## 7. Align the samples fasta files to the metagene
## 8. Create transcript counts based on the alignments to the metagene
## 9. Normalise the counts
## 10. Generate a heatmap and PCA plot


