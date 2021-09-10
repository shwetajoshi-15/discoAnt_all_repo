# discoAnt
The repository is under-construction.

## Prepare FASTA files in a folder
Date_of_sequencing.barcode01.fa (e.g. 2019_12_16.barcode01.fa)


1. git clone repository
2. cd discoAnt
3. Update discoAnt_params.txt with - Gene name, start and end coordinates and strands
4. Update discoAnt_params.txt with - path to FASTA folder and discoAnt folder
5. When running the pipeline for the first time - bash discoAnt_setup.sh
6. conda activate discoAnt.env
7. cd /path/to/folder/discoAnt/programs/cDNA_Cupcake
8. python setup.py build
9. python setup.py install
10. Submit discoAnt.sh in the appropriate submission script

## What does discoAnt do?

1. Quanlity control for the fasta files 
1.a. Generates a file with no. of reads
1.b. average read length
1.c. read length in each barcode/sample 
2. Align (and correcting) the sample fasta files to reference genome
3. Merge the primary alignments
4. Generate transcript files (GTF) for the merged alignments
4. Classify the transcripts based on known annotations 
6. Create a metagene based on the transcripts
7. Align the samples fasta files to the metagene
9. Create transcript counts based on the alignments to the metagene (under construction)
10. Normalise the counts (under construction)
11. Generate a heatmap and PCA plot (under construction)


