#!/bin/bash

source discoAnt_params.txt

## Installing programs

https://github.com/mortazavilab/TranscriptClean.git $PROGRAMS/.
git clone https://github.com/ConesaLab/SQANTI3.git $PROGRAMS/.
git clone https://github.com/Magdoll/cDNA_Cupcake.git $PROGRAMS/.


python $PROGRAMS/cDNA_Cupcake/setup.py build
python $PROGRAMS/cDNA_Cupcake/setup.py install


## GENCODE genome reference and annotation

	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.p13.genome.fa.gz $REF_HG38/.
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz $REF_HG38/.

## Editing the reference and annotation files

	gunzip $REF_HG38/GRCh38.p13.genome.fa.gz
	gunzip $REF_HG38/gencode.v34.transcripts.fa.gz
	gunzip $REF_HG38/gencode.v34.annotation.gtf.gz

	cut -d " " -f 1 $REF_HG38/GRCh38.p13.genome.fa > $REF_HG38/GRCh38.p13.genome_edit.fa
	cut -d "|" -f 1 $REF_HG38/gencode.v34.transcripts.fa > $REF_HG38/gencode.v34.transcripts_edit.fa

## Downlaoding reference filed for SQANTI3 annotation

	wget http://reftss.clst.riken.jp/datafiles/3.3/human/refTSS_v3.3_human_coordinate.hg38.bed.gz $REF_HG38/.
	wget https://raw.githubusercontent.com/Magdoll/images_public/master/SQANTI2_support_data/human.polyA.list.txt $REF_HG38/.
	wget https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz $REF_HG38/.

## Downlaoding reference filed for SQANTI3 annotation

	gunzip $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed.gz 
	gunzip $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed.gz
	
	sed -i 's/^/chr/' $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed
