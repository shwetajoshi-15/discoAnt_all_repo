#/bin/bash!

### This version of the pipeline (as of 13th October 2021) is based on Gencode v35

## This script performs the following steps
## 1. Quality control for the fasta files
## 1.a. Generates a file with no. of reads
## 1.b. average read length
## 1.c. read length in each barcode/sample
## 2. Align the sample fasta files to reference genome
## 3. Filter out off-target alignments
## 4. Merge alignments from all samples
## 5. Correct the alignments based on reference genome
## 6. Generate transcript files (GTF) for the merged alignments
## 7. Classify the transcripts based on known annotations
## 8.a. Create a metagene based on known and novel exons
## 8.b. Create a GTF file for the metagene
## 9. Align the samples fasta files to the metagene
## 10. Generate counts based on the alignments to the metagene
## 11. Merge the transcript annotations and counts in "results"
## 12. Generate a heatmap and PCA plot



source discoAnt_params.txt
export PYTHONPATH="$PROGRAMS/cDNA_Cupcake/sequence/:$PYTHONPATH"

echo "Making folders"

mkdir -p $DISCOANT/$GENE
mkdir -p $DISCOANT/$GENE/fasta_stats
mkdir -p $DISCOANT/$GENE/minimap2
mkdir -p $DISCOANT/$GENE/minimap2_target
mkdir -p $DISCOANT/$GENE/transcriptclean
mkdir -p $DISCOANT/$GENE/stringtie
mkdir -p $DISCOANT/$GENE/stringtie/metagene
mkdir -p $DISCOANT/$GENE/stringtie/clean_metagene
mkdir -p $DISCOANT/$GENE/sqanti3
mkdir -p $DISCOANT/$GENE/minimap2_metagene
mkdir -p $DISCOANT/$GENE/minimap2_metagene/metagene
mkdir -p $DISCOANT/$GENE/minimap2_metagene/clean_metagene
mkdir -p $DISCOANT/$GENE/minimap2_metagene_transcriptome
mkdir -p $DISCOANT/$GENE/minimap2_metagene_transcriptome/metagene
mkdir -p $DISCOANT/$GENE/minimap2_metagene_transcriptome/clean_metagene
mkdir -p $DISCOANT/$GENE/featurecounts
mkdir -p $DISCOANT/$GENE/salmon
mkdir -p $DISCOANT/$GENE/results

##########                                          ##########
##########  1. Quality control for the fasta files  ##########
##########                                          ##########

echo "Extracting the number of reads and read lengths"

	basename -s .fa $FASTA/*.fa | sed 's/^.*bar/bar/g' > $DISCOANT/$GENE/fasta_stats/barcodes.txt

	file=$DISCOANT/$GENE/fasta_stats/$GENE_stats.COMPLETED

	if [[ ! -f $DISCOANT/$GENE/fasta_stats/"$GENE"_stats.COMPLETED ]]
	then

	for filename in $FASTA/*.fa
	do
   	base=$(basename $filename .fa)
   	echo "On sample : $base"

   	cat $FASTA/${base}.fa | grep -c "^>" >> $DISCOANT/$GENE/fasta_stats/"$GENE"_tmp1.txt
   	paste $DISCOANT/$GENE/fasta_stats/barcodes.txt $DISCOANT/$GENE/fasta_stats/"$GENE"_tmp1.txt | column -s $'\t' -t > $DISCOANT/$GENE/fasta_stats/"$GENE"_total_number_of_reads.txt

   	mkdir -p $DISCOANT/$GENE/fasta_stats/"$GENE"_read_lengths
   	bioawk -c fastx '{print length($seq)}' $FASTA/${base}.fa > $DISCOANT/$GENE/fasta_stats/"$GENE"_read_lengths/${base}.txt

   	awk '{/>/&&++a||b+=length()}END{print b/a}' $FASTA/${base}.fa >> $DISCOANT/$GENE/fasta_stats/"$GENE"_tmp2.txt
   	paste $DISCOANT/$GENE/fasta_stats/barcodes.txt $DISCOANT/$GENE/fasta_stats/"$GENE"_tmp2.txt | column -s $'\t' -t > $DISCOANT/$GENE/fasta_stats/"$GENE"_mean_read_lengths.txt

 	rm $DISCOANT/$GENE/fasta_stats/*tmp* && rm $DISCOANT/$GENE/fasta_stats/barcodes.txt && touch $DISCOANT/$GENE/fasta_stats/"$GENE"_stats.COMPLETED
	done

	else
	echo "FASTA stats are present in $DISCOANT/$GENE/fasta_stats"
	fi


##########                                                      ##########
########## 2. Align the sample fasta files to reference genome  ##########
##########                                                      ##########


echo "minimap2 - Mapping fasta files to genome"

	if [[ ! -f $DISCOANT/$GENE/minimap2/"$GENE"_sorted_pri_align.COMPLETED ]]
	then

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	minimap2 -ax splice $REF_HG38/GRCh38.p13.genome_edit.fa $FASTA/${base}.fa > $DISCOANT/$GENE/minimap2/${base}.sam
	samtools view -S -h -b $DISCOANT/$GENE/minimap2/${base}.sam | samtools sort - > $DISCOANT/$GENE/minimap2/${base}_sorted.bam
	samtools view -h -F 2308 $DISCOANT/$GENE/minimap2/${base}_sorted.bam | samtools sort - > $DISCOANT/$GENE/minimap2/${base}_pri_sorted.bam
	
	done
	touch $DISCOANT/$GENE/minimap2/"$GENE"_sorted_pri_align.COMPLETED

	else
	echo "Alignments to the genome are present in $DISCOANT/$GENE/minimap2"
	fi

##########                                      ##########
########## 3. Filter out off-target alignments  ##########
##########                                      ##########
### This filering is based on the the start and end site cooordinates provided in the GENE_params.txt file. 

echo "Extracting alignments to the target region/gene"

	if [[ ! -f $DISCOANT/$GENE/minimap2_target/"$GENE"_sorted_pri_tar_align.COMPLETED  ]]
	then

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"
   	samtools index $DISCOANT/$GENE/minimap2/${base}_pri_sorted.bam
	samtools view -h $DISCOANT/$GENE/minimap2/${base}_pri_sorted.bam "$CHR:$GENE_START-$GENE_END" > $DISCOANT/$GENE/minimap2_target/${base}_pri_tar_sorted.bam
	done
	touch $DISCOANT/$GENE/minimap2_target/"$GENE"_sorted_pri_tar_align.COMPLETED

	else
	echo "Filtered primary alignments are present in $DISCOANT/$GENE/minimap2_target"
	fi

##########                                      ##########
########## 4. Merge alignments from all samples ##########
##########                                      ##########

echo "Merging minimap2 primary alignments"

	if [[ ! -f $DISCOANT/$GENE/minimap2_target/"$GENE"_sorted_pri_tar_merged_align.COMPLETED ]]
	then

	samtools merge -f $DISCOANT/$GENE/minimap2_target/"$GENE"_merged.bam $DISCOANT/$GENE/minimap2_target/*_pri_tar_sorted.bam
	samtools sort $DISCOANT/$GENE/minimap2_target/"$GENE"_merged.bam > $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam
	samtools index $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam
	samtools view -h -o $DISCOANT/$GENE/minimap2_target/"$GENE"_merged.sam $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam

	touch $DISCOANT/$GENE/minimap2_target/"$GENE"_sorted_pri_tar_merged_align.COMPLETED

	else
	echo "Merged alignments are present in $DISCOANT/$GENE/minimap2_target"
	fi

##########                                                     ##########
########## 5. Correct the alignments based on reference genome ##########
##########                                                     ##########
### This correction step uses the reference genome. 

echo "Correcting the merged primary alignments with TranscriptClean"

	if [[ ! -f $DISCOANT/$GENE/transcriptclean/"$GENE"_merged_clean_sorted.COMPLETED ]]
	then

	python $PROGRAMS/TranscriptClean/TranscriptClean.py \
	--sam $DISCOANT/$GENE/minimap2_target/"$GENE"_merged.sam \
	--primaryOnly --genome $REF_HG38/GRCh38.p13.genome_edit.fa \
	--outprefix $DISCOANT/$GENE/transcriptclean/$GENE

	samtools view -S -h -b $DISCOANT/$GENE/transcriptclean/"$GENE"_clean.sam | samtools sort - > $DISCOANT/$GENE/transcriptclean/"$GENE"_clean_sorted.bam

	touch $DISCOANT/$GENE/transcriptclean/"$GENE"_merged_clean_sorted.COMPLETED

	else
	echo "TranscriptClean-ed alignments are present in $DISCOANT/$GENE/transcriptclean"
	fi

##########                                                              ##########
########## 6. Generate transcript files (GTF) for the merged alignments ##########
##########                                                              ##########

echo "Constructing transcripts based on the corrected alignments"

	if [[ ! -f $DISCOANT/$GENE/stringtie/"$GENE"_transcripts_pre_clean.COMPLETED ]]
	then

	stringtie -L -t $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam \
	-G $REF_HG38/gencode.v35.annotation.gtf \
	-o $DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE.gtf

	touch $DISCOANT/$GENE/stringtie/"$GENE"_transcripts_pre_clean.COMPLETED

	else
	echo "Transcripts constructed (pre correction) are present in $DISCOANT/$GENE/stringtie"
	fi

	if [[ ! -f $DISCOANT/$GENE/stringtie/"$GENE"_transcripts_post_clean.COMPLETED ]]
	then

	stringtie -L -t $DISCOANT/$GENE/transcriptclean/"$GENE"_clean_sorted.bam \
	-G $REF_HG38/gencode.v35.annotation.gtf \
	-o $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE.gtf

	touch $DISCOANT/$GENE/stringtie/"$GENE"_transcripts_post_clean.COMPLETED

	else
	echo "Transcripts constructed (post correction) are present in $DISCOANT/$GENE/stringtie"
	fi

##########                                                         ##########
########## 7.a Classify the transcripts based on known annotations ##########
##########                                                         ##########

echo "Annotating the transcripts with SQANTI3"

	if [[ ! -f $DISCOANT/$GENE/sqanti3/"$GENE"_sqanti3.COMPLETED ]]
	then
	python $PROGRAMS/SQANTI3-1.3/sqanti3_qc.py \
	--gtf $DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE.gtf \
	$REF_HG38/gencode.v35.annotation.gtf $REF_HG38/GRCh38.p13.genome_edit.fa \
	--cage_peak $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed \
	--polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/human.polyA.list.txt \
	-d $DISCOANT/$GENE/sqanti3 -o "$GENE"

	cat $DISCOANT/$GENE/sqanti3/"$GENE"_classification.txt | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$15,$17,$30,$31,$37,$40,$41,$42,$43}' > $DISCOANT/$GENE/sqanti3/"$GENE"_sqanti_matrix.txt

	touch $DISCOANT/$GENE/sqanti3/"$GENE"_sqanti3.COMPLETED

	else
	echo "Annotated transcripts (pre correction) are present in $DISCOANT/$GENE/sqanti3"
	fi

echo "Annotating the transcripts with SQANTI3"

	if [[ ! -f $DISCOANT/$GENE/sqanti3/"$GENE"_sqanti3_clean.COMPLETED ]]
	then

	python $PROGRAMS/SQANTI3-1.3/sqanti3_qc.py \
	--gtf $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE.gtf \
	$REF_HG38/gencode.v35.annotation.gtf $REF_HG38/GRCh38.p13.genome_edit.fa \
	--cage_peak $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed \
	--polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/human.polyA.list.txt \
	-d $DISCOANT/$GENE/sqanti3 -o "$GENE"_clean

	cat $DISCOANT/$GENE/sqanti3/"$GENE"_clean_classification.txt | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$15,$17,$30,$31,$37,$40,$41,$42,$43}' > $DISCOANT/$GENE/sqanti3/"$GENE"_sqanti_matrix.txt

	touch $DISCOANT/$GENE/sqanti3/"$GENE"_sqanti3_clean.COMPLETED

	else
	echo "Annotated transcripts (post correction) are present in $DISCOANT/$GENE/sqanti3"
	fi


##########                                                       ##########
########## 8.a. Create a metagene based on known and novel exons ##########
##########                                                       ##########
### This step uses the stringtie transcripts to create an improved transcriptome (metagene). 
### This section of the script borrows from TAQLoRe.

echo "Creating a metagene FASTA"

	if [[ ! -f $DISCOANT/$GENE/stringtie/"$GENE"_metagene.COMPLETED ]]
	then

	cat $DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE.gtf | awk '{ if ($3 == "exon") { print } }' | awk 'BEGIN{OFS="\t"}{print $1, $4, $5, ".", ".", -v "$STRAND"}' > $DISCOANT/$GENE/stringtie/"$GENE"_all_exons.bed

	bedtools getfasta -s -fi $REF_HG38/GRCh38.p13.genome_edit.fa \
	-bed $DISCOANT/$GENE/stringtie/"$GENE"_all_exons.bed -fo $DISCOANT/$GENE/stringtie/"$GENE"_meta_gene_exons.fa

	echo ">meta_gene_$GENE" > $DISCOANT/$GENE/stringtie/"$GENE"_meta_gene.fa && \
	cat $DISCOANT/$GENE/stringtie/"$GENE"_meta_gene_exons.fa | grep -v "^>" | tr -d '\n' >> $DISCOANT/$GENE/stringtie/"$GENE"_meta_gene.fa && \
	echo >> $DISCOANT/$GENE/stringtie/"$GENE"_meta_gene.fa

	touch $DISCOANT/$GENE/stringtie/"$GENE"_metagene.COMPLETED

	else
	echo "Metagene (pre correction transcripts) is present in $DISCOANT/$GENE/stringtie"
	fi

## check if the strand variable is assigned correctly
	if [[ ! -f $DISCOANT/$GENE/stringtie/"$GENE"_metagene_clean.COMPLETED ]]
	then
	cat $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE.gtf | awk '{ if ($3 == "exon") { print } }' | awk 'BEGIN{OFS="\t"}{print $1, $4, $5, ".", ".", -v "$STRAND"}' > $DISCOANT/$GENE/stringtie/"$GENE"_clean_all_exons.bed

	bedtools getfasta -s -fi $REF_HG38/GRCh38.p13.genome_edit.fa \
	-bed $DISCOANT/$GENE/stringtie/"$GENE"_clean_all_exons.bed -fo $DISCOANT/$GENE/stringtie/"$GENE"_clean_meta_gene_exons.fa

	echo ">meta_gene_$GENE" > $DISCOANT/$GENE/stringtie/"$GENE"_clean_meta_gene.fa && \
	cat $DISCOANT/$GENE/stringtie/"$GENE"_clean_meta_gene_exons.fa | grep -v "^>" | tr -d '\n' >> $DISCOANT/$GENE/stringtie/"$GENE"_clean_meta_gene.fa && \
	echo >> $DISCOANT/$GENE/stringtie/"$GENE"_clean_meta_gene.fa

	touch $DISCOANT/$GENE/stringtie/"$GENE"_metagene_clean.COMPLETED

	else
	echo "Metagene (post correction transcripts) is present in $DISCOANT/$GENE/stringtie"
	fi

##########                                         ##########
########## 8.b. Create a GTF file for the metagene ##########
##########                                         ##########
### The generation of counts requires a reference GTF. 
### The stringtie is modified by adding metagene coordinates for quantification. 

echo "Modifying stringtie GTF to create a metagene GTF"

	if [[ ! -f $DISCOANT/$GENE/stringtie/"$GENE"_metagene_gtf.COMPLETED ]]
	then
	cat $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE.gtf | awk '{ if ($3 == "exon") { print $4 } }' > $DISCOANT/$GENE/stringtie/"$GENE"_clean_exon_start.txt
	cat $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE.gtf | awk '{ if ($3 == "exon") { print $5 } }' > $DISCOANT/$GENE/stringtie/"$GENE"_clean_exon_end.txt

	Rscript $PROGRAMS/exon_coord_conversion.R --col1 $DISCOANT/$GENE/stringtie/"$GENE"_clean_exon_start.txt \
	--col2 $DISCOANT/$GENE/stringtie/"$GENE"_clean_exon_end.txt \
	--out $DISCOANT/$GENE/stringtie/"$GENE"_clean_metagene_exon_coord.txt

	cat $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE.gtf | \
	awk '{ if ($3 == "exon") { print $9,$10,$11,$12 } }' > $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE_modified_tmp1.gtf

	cat $DISCOANT/$GENE/stringtie/"$GENE"_clean_metagene_exon_coord.txt | \
	awk -v GENE="$GENE" -v STRAND="$STRAND" 'BEGIN{{OFS="\t"}}{ print "meta_gene_"GENE, "Stringtie", "exon", $1, $2, ".", STRAND, "."}' > $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE_modified_tmp2.gtf

 	paste $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE_modified_tmp2.gtf \
	$DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE_modified_tmp1.gtf > $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE_modified.gtf

	touch $DISCOANT/$GENE/stringtie/"$GENE"_clean_metagene_gtf.COMPLETED

	else
	echo "Metagene is present in $DISCOANT/$GENE/stringtie"
	fi

	if [[ ! -f $DISCOANT/$GENE/stringtie/"$GENE"_metagene_gtf.COMPLETED ]]
	then
	cat $DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE.gtf | awk '{ if ($3 == "exon") { print $4 } }' > $DISCOANT/$GENE/stringtie/"$GENE"_exon_start.txt
	cat $DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE.gtf | awk '{ if ($3 == "exon") { print $5 } }' > $DISCOANT/$GENE/stringtie/"$GENE"_exon_end.txt

	Rscript $PROGRAMS/exon_coord_conversion.R --col1 $DISCOANT/$GENE/stringtie"$GENE"_exon_start.txt \
	--col2 $DISCOANT/$GENE/stringtie/"$GENE"_exon_end.txt \
	--out $DISCOANT/$GENE/stringtie/"$GENE"_metagene_exon_coord.txt

	cat $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE.gtf | \
	awk '{ if ($3 == "exon") { print $9,$10,$11,$12 } }' > $DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE_modified_tmp1.gtf

	cat $DISCOANT/$GENE/stringtie/"$GENE"_metagene_exon_coord.txt | \
	awk -v GENE="$GENE" -v STRAND="$STRAND" 'BEGIN{{OFS="\t"}}{ print "meta_gene_"GENE, "Stringtie", "exon", $1, $2, ".", STRAND, "."}' > $DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE_modified_tmp2.gtf

 	paste $DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE_modified_tmp2.gtf \
	$DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE_modified_tmp1.gtf > $DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE_modified.gtf

	rm $DISCOANT/$GENE/stringtie/*tmp* && touch $DISCOANT/$GENE/stringtie/"$GENE"_metagene_gtf.COMPLETED
	else
	echo "Metagene is present in $DISCOANT/$GENE/stringtie"
	fi

##########                                                   ##########
########## 9.a Align the samples fasta files to the metagene ##########
##########                                                   ##########

echo "Mapping sample FASTA to the metagene"

	if [[ ! -f $DISCOANT/$GENE/minimap2_metagene/"$GENE"_sorted_pri_align_clean.COMPLETED ]]
	then
	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	minimap2 -ax splice $DISCOANT/$GENE/stringtie/"$GENE"_clean_meta_gene.fa $FASTA/${base}.fa > $DISCOANT/$GENE/minimap2_metagene/${base}_clean.sam
	samtools view -S -h -b $DISCOANT/$GENE/minimap2_metagene/${base}_clean.sam | samtools sort - > $DISCOANT/$GENE/minimap2_metagene/${base}_clean_sorted.bam
	samtools view -h -F 2308 $DISCOANT/$GENE/minimap2_metagene/${base}_clean_sorted.bam | samtools sort - > $DISCOANT/$GENE/minimap2_metagene/${base}_clean_pri_sorted.bam
	samtools index $DISCOANT/$GENE/minimap2_metagene/${base}_clean_pri_sorted.bam
	done
	touch $DISCOANT/$GENE/minimap2_metagene/"$GENE"_sorted_pri_align_clean.COMPLETED
	else
	echo "Alignments to the metagene are present in $DISCOANT/$GENE/minimap2_metagene"
	fi

	if [[ ! -f $DISCOANT/$GENE/minimap2_metagene/metagene/"$GENE"_sorted_pri_align.COMPLETED ]]
	then
	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	minimap2 -ax splice $DISCOANT/$GENE/stringtie/"$GENE"_meta_gene.fa $FASTA/${base}.fa > $DISCOANT/$GENE/minimap2_metagene/metagene/${base}.sam
	samtools view -S -h -b $DISCOANT/$GENE/minimap2_metagene/metagene/${base}.sam | samtools sort - > $DISCOANT/$GENE/minimap2_metagene/metagene/${base}_sorted.bam
	samtools view -h -F 2308 $DISCOANT/$GENE/minimap2_metagene/metagene/${base}_sorted.bam | samtools sort - > $DISCOANT/$GENE/minimap2_metagene/metagene/${base}_pri_sorted.bam
	samtools index $DISCOANT/$GENE/minimap2_metagene/metagene/${base}_pri_sorted.bam
	done
	touch $DISCOANT/$GENE/minimap2_metagene/metagene/"$GENE"_sorted_pri_align.COMPLETED
	else
	echo "Alignments to the metagene are present in $DISCOANT/$GENE/minimap2_metagene"
	fi

##########                                                                   ##########
########## 9.b Align the samples fasta files to the metagene (transcriptome) ##########
##########                                                                   ##########

echo "Mapping sample FASTA to the metagene"

	if [[ ! -f $DISCOANT/$GENE/minimap2_metagene_transcriptome/"$GENE"_sorted_pri_align_clean.COMPLETED ]]
	then
  for filename in $FASTA/*.fa
  do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	minimap2 -ax splice $DISCOANT/$GENE/sqanti3/"$GENE"_clean_corrected.fasta $FASTA/${base}.fa > $DISCOANT/$GENE/minimap2_metagene_transcriptome/${base}_clean.sam
	samtools view -S -h -b $DISCOANT/$GENE/minimap2_metagene_transcriptome/${base}_clean.sam | samtools sort - > $DISCOANT/$GENE/minimap2_metagene_transcriptome/${base}_clean_sorted.bam
	samtools view -h -F 2308 $DISCOANT/$GENE/minimap2_metagene_transcriptome/${base}_clean_sorted.bam | samtools sort - > $DISCOANT/$GENE/minimap2_metagene_transcriptome/${base}_clean_pri_sorted.bam
	samtools index $DISCOANT/$GENE/minimap2_metagene_transcriptome/${base}_clean_pri_sorted.bam
	done
	touch $DISCOANT/$GENE/minimap2_metagene_transcriptome/"$GENE"_sorted_pri_align_clean.COMPLETED
	else
	echo "Alignments to the metagene are present in $DISCOANT/$GENE/minimap2_metagene_transcriptome/clean_metagene"
	fi

	if [[ ! -f $DISCOANT/$GENE/minimap2_metagene_transcriptome/metagene/"$GENE"_sorted_pri_align.COMPLETED ]]
	then
	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	minimap2 -ax splice $DISCOANT/$GENE/sqanti3/"$GENE"_corrected.fasta  $FASTA/${base}.fa > $DISCOANT/$GENE/minimap2_metagene_transcriptome/metagene/${base}.sam
	samtools view -S -h -b $DISCOANT/$GENE/minimap2_metagene_transcriptome/metagene/${base}.sam | samtools sort - > $DISCOANT/$GENE/minimap2_metagene_transcriptome/metagene/${base}_sorted.bam
	samtools view -h -F 2308 $DISCOANT/$GENE/minimap2_metagene_transcriptome/metagene/${base}_sorted.bam | samtools sort - > $DISCOANT/$GENE/minimap2_metagene_transcriptome/metagene/${base}_pri_sorted.bam
	samtools index $DISCOANT/$GENE/minimap2_metagene_transcriptome/metagene/${base}_pri_sorted.bam
	done
	touch $DISCOANT/$GENE/minimap2_metagene_transcriptome/metagene/"$GENE"_sorted_pri_align.COMPLETED
	else
	echo "Alignments to the metagene are present in $DISCOANT/$GENE/minimap2_metagene_transcriptome"
	fi


##########                                                              ##########
########## 10.a Generate counts based on the alignments to the metagene ##########
##########                                                              ##########

echo "Generating read counts for the metagene alignments"

	if [[ ! -f $DISCOANT/$GENE/featurecounts/"$GENE"_metagene_counts.COMPLETED ]]
	then
	featureCounts -L -g transcript_id -a $DISCOANT/$GENE/stringtie/metagene/"$GENE"_STRINGTIE_modified.gtf \
	-o $DISCOANT/$GENE/featurecounts/"$GENE"_counts.txt $DISCOANT/$GENE/minimap2_metagene/metagene/*_pri_sorted.bam

	cat $DISCOANT/$GENE/featurecounts/"$GENE"_counts.txt | cut -f2,3,4,5,6 --complement | awk 'FNR > 2' > $DISCOANT/$GENE/featurecounts/"$GENE"_counts_matrix.txt
	cat $DISCOANT/$GENE/"$GENE"_samplenames.txt | tr '\n' '\t' | awk 'FNR > 0' - $DISCOANT/$GENE/featurecounts/"$GENE"_counts_matrix.txt > $DISCOANT/$GENE/featurecounts/"$GENE"_counts_matrix_samplenames.txt

	touch $DISCOANT/$GENE/featurecounts/"$GENE"_metagene_counts.COMPLETED
	else
	echo "Counts matrix is present in $DISCOANT/$GENE/featurecounts"
	fi

	if [[ ! -f $DISCOANT/$GENE/featurecounts/"$GENE"_metagene_counts_clean.COMPLETED ]]
	then
	featureCounts -L -g transcript_id -a $DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE_modified.gtf \
	-o $DISCOANT/$GENE/featurecounts/"$GENE"_clean_counts.txt $DISCOANT/$GENE/minimap2_metagene/*_clean_pri_sorted.bam

	cat $DISCOANT/$GENE/featurecounts/"$GENE"_clean_counts.txt | cut -f2,3,4,5,6 --complement | awk 'FNR > 2' > $DISCOANT/$GENE/featurecounts/"$GENE"_clean_counts_matrix.txt
	cat $DISCOANT/$GENE/"$GENE"_samplenames.txt | tr '\n' '\t' | awk 'FNR > 0' - $DISCOANT/$GENE/featurecounts/"$GENE"_clean_counts_matrix.txt > $DISCOANT/$GENE/featurecounts/"$GENE"_clean_counts_matrix_samplenames.txt

	touch $DISCOANT/$GENE/featurecounts/"$GENE"_metagene_counts_clean.COMPLETED
	else
	echo "Counts matrix is present in $DISCOANT/$GENE/featurecounts"
	fi

##########                                                                              ##########
########## 10.b Generate counts based on the alignments to the metagene (transcriptome) ##########
##########                                                                              ##########

echo "Quantifying transcripts with salmon"

	if [[ ! -f $DISCOANT/$GENE/salmon/"$GENE"_salmon_clean.COMPLETED ]]
	then
	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	salmon quant -t $DISCOANT/$GENE/sqanti3/"$GENE"_clean_corrected.fasta -l A -a $DISCOANT/$GENE/minimap2_metagene_transcriptome/${base}_clean_pri_sorted.bam -o $DISCOANT/$GENE/salmon/${base}_clean
        
	done
	touch $DISCOANT/$GENE/salmon/"$GENE"_salmon_clean.COMPLETED
	else
	echo "Alignments to the metagene are present in $DISCOANT/$GENE/salmon/clean_metagene"
	fi

	if [[ ! -f $DISCOANT/$GENE/salmon/"$GENE"_salmon.COMPLETED ]]
	then
	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	salmon quant -t $DISCOANT/$GENE/sqanti3/"$GENE"_corrected.fasta -l A -a $DISCOANT/$GENE/minimap2_metagene_transcriptome/metagene/${base}_pri_sorted.bam -o $DISCOANT/$GENE/salmon/${base}

	done
	touch $DISCOANT/$GENE/salmon/"$GENE"_salmon.COMPLETED
	else
	echo "Alignments to the metagene are present in $DISCOANT/$GENE/salmon/metagene"
	fi


	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	cut -d\t -f4 $DISCOANT/$GENE/salmon/${base}/quant.sf | awk 'FNR > 1' | paste - | column -s $'\t' -t >> $DISCOANT/$GENE/salmon/"$GENE"_salmon_TPM_tmp1.txt
	cut -d\t -f5 $DISCOANT/$GENE/salmon/${base}/quant.sf | awk 'FNR > 1' | paste - | column -s $'\t' -t >> $DISCOANT/$GENE/salmon/"$GENE"_salmon_counts_tmp1.txt

	done

##########                                                              ##########
########## 11. Merge the transcript annotations and counts in "results" ##########
##########                                                              ##########

echo "Merging SQANTI3 annotations with counts matrix"

	paste $DISCOANT/$GENE/sqanti3/"$GENE"_sqanti_matrix.txt $DISCOANT/$GENE/featurecounts/"$GENE"_counts_matrix_samplenames.txt \
	| column -s $'\t' -t > $DISCOANT/$GENE/results/$GENE_annotated_transcript_counts.txt
