#/bin/bash!

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


source discoAnt_params.txt
export PYTHONPATH="$PROGRAMS/cDNA_Cupcake/sequence/:$PYTHONPATH"

echo "Making folders"

mkdir -p $DISCOANT/$GENE
mkdir -p $DISCOANT/$GENE/fasta_stats
mkdir -p $DISCOANT/$GENE/minimap2
mkdir -p $DISCOANT/$GENE/minimap2_target
mkdir -p $DISCOANT/$GENE/transcriptclean
mkdir -p $DISCOANT/$GENE/stringtie
mkdir -p $DISCOANT/$GENE/sqanti3
mkdir -p $DISCOANT/$GENE/minimap2_metagene

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


echo "Merging minimap2 primary alignments"

if [[ ! -f $DISCOANT/$GENE/minimap2_target/"$GENE"_sorted_pri_tar_merged_align.COMPLETED ]]

then

   samtools merge $DISCOANT/$GENE/minimap2_target/"$GENE"_merged.bam $DISCOANT/$GENE/minimap2_target/*_pri_tar_sorted.bam 
   samtools sort $DISCOANT/$GENE/minimap2_target/"$GENE"_merged.bam > $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam
   samtools index $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam 
   samtools view -h -o $DISCOANT/$GENE/minimap2_target/"$GENE"_merged.sam $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam

   touch $DISCOANT/$GENE/minimap2_target/"$GENE"_sorted_pri_tar_merged_align.COMPLETED
   
else 

echo "Merged alignments are present in $DISCOANT/$GENE/minimap2_target"
	
fi


echo "Correcting the merged primary alignments with TranscriptClean"

if [[ ! -f $DISCOANT/$GENE/transcriptclean/"$GENE"_merged_clean_sorted.COMPLETED ]]

then

   python $PROGRAMS/TranscriptClean/TranscriptClean.py \
   --sam $DISCOANT/$GENE/minimap2/"$GENE"_merged.sam \
   --primaryOnly --genome $REF_HG38/GRCh38.p13.genome_edit.fa \
   --outprefix $DISCOANT/$GENE/transcriptclean/$GENE

   samtools view -S -h -b $DISCOANT/$GENE/transcriptclean/"$GENE"_clean.sam | samtools sort - > $DISCOANT/$GENE/transcriptclean/"$GENE"_clean_sorted.bam

   touch $DISCOANT/$GENE/transcriptclean/"$GENE"_merged_clean_sorted.COMPLETED

else

echo "TranscriptClean-ed alignments are present in $DISCOANT/$GENE/transcriptclean"
	
fi


echo "Constructing transcripts based on the corrected alignments"

if [[ ! -f $DISCOANT/$GENE/stringtie/"$GENE"_transcripts_pre_clean.COMPLETED ]]

then

   stringtie -t $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam \
   -G $REF_HG38/gencode.v35.annotation.gtf \
   -A $DISCOANT/$GENE/stringtie/"$GENE"_gene_abund.tab \
   -C $DISCOANT/$GENE/stringtie/"$GENE"_cov_ref.gtf \
   -o $DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE.gtf
   
   touch $DISCOANT/$GENE/stringtie/"$GENE"_transcripts_pre_clean.COMPLETED

else

echo "Transcripts constructed pre correction are present in $DISCOANT/$GENE/stringtie"
	
fi

if [[ ! -f $DISCOANT/$GENE/stringtie/"$GENE"_transcripts_post_clean.COMPLETED ]]

then

   stringtie -t $DISCOANT/$GENE/transcriptclean/"$GENE"_clean_sorted.bam \
   -G $REF_HG38/gencode.v35.annotation.gtf \
   -A $DISCOANT/$GENE/stringtie/"$GENE"_clean_gene_abund.tab \
   -C $DISCOANT/$GENE/stringtie/"$GENE"_clean_cov_ref.gtf \
   -o $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE.gtf
   
   touch $DISCOANT/$GENE/stringtie/"$GENE"_transcripts_post_clean.COMPLETED

else

echo "Transcripts constructed post correction are present in $DISCOANT/$GENE/stringtie"
	
fi

echo "Annotating the transcripts with SQANTI3"

   python $PROGRAMS/SQANTI3/sqanti3_qc.py \
   $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE.gtf \
   $REF_HG38/gencode.v35.annotation.gtf $REF_HG38/GRCh38.p13.genome_edit.fa \
   --cage_peak $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed \
   --polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/human.polyA.list.txt \
   -d $DISCOANT/$GENE/sqanti3 -o "$GENE"_clean
   
      python $PROGRAMS/SQANTI3/sqanti3_qc.py \
   $DISCOANT/$GENE/stringtie/"$GENE"_STRINGTIE.gtf \
   $REF_HG38/gencode.v35.annotation.gtf $REF_HG38/GRCh38.p13.genome_edit.fa \
   --cage_peak $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed \
   --polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/human.polyA.list.txt \
   -d $DISCOANT/$GENE/sqanti3 -o "$GENE"

echo "Creating a metagene FASTA"
   ## check if the strand variable is assigned correctly

   cat $DISCOANT/$GENE/stringtie/"$GENE"_clean_STRINGTIE.gtf | awk '{ if ($3 == "exon") { print } }' | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5, ".", ".", "$STRAND"}' > $DISCOANT/$GENE/stringtie/"$GENE"_clean_all_exons.bed

   bedtools getfasta -s -fi $REF_HG38/GRCh38.p13.genome_edit.fa -bed $DISCOANT/$GENE/stringtie/"$GENE"_clean_all_exons.bed -fo $DISCOANT/$GENE/stringtie/"$GENE"_clean_meta_gene_exons.fa

   echo ">meta_gene_$GENE" > $DISCOANT/$GENE/stringtie/"$GENE"_clean_meta_gene.fa && \
   cat $DISCOANT/$GENE/stringtie/"$GENE"_clean_meta_gene_exons.fa | grep -v "^>" | tr -d '\n' >> $DISCOANT/$GENE/stringtie/"$GENE"_clean_meta_gene.fa && \
   echo >> $DISCOANT/$GENE/stringtie/"$GENE"_clean_meta_gene.fa 
 
echo "Mapping sample FASTA to the metagene"

   for filename in $FASTA/*.fa

   do
   base=$(basename $filename .fa)
   echo "On sample : $base" 

   minimap2 -ax splice $DISCOANT/$GENE/stringtie/"$GENE"_meta_gene.fa $FASTA/${base}.fa > $DISCOANT/$GENE/minimap2_metagene/${base}.sam 

   done
