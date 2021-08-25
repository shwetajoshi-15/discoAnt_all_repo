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


source $DISCOANT/discoAnt_params.txt

echo "Making folders"

mkdir -p $DISCOANT/$GENE/fasta_stats
mkdir -p $DISCOANT/$GENE/minimap2
mkdir -p $DISCOANT/$GENE/transcriptclean
mkdir -p $DISCOANT/$GENE/stringtie
mkdir -p $DISCOANT/$GENE/sqanti3
mkdir -p $DISCOANT/$GENE/minimap2_metagene

echo "minimap2 - Mapping fasta files to genome"

   for filename in $DISCOANT/$GENE/fasta/downsample_500/*.fa

   do

   base=$(basename $filename .fa)
   echo "On sample : $base"

   minimap2 -ax splice $REF_HG38/GRCh38.p13.genome_edit.fa \
   $DISCOANT/$GENE/fasta/downsample_500/${base}.fa > $DISCOANT/$GENE/minimap2/${base}.sam
   samtools view -S -h -b $DISCOANT/$GENE/minimap2/${base}.sam | samtools sort - > $DISCOANT/$GENE/minimap2/${base}_sorted.bam
   samtools view -h -F 2308 $DISCOANT/$GENE/minimap2/${base}_sorted.bam | samtools sort - > $DISCOANT/$GENE/minimap2/${base}_pri_sorted.bam 
   
   done

echo "Merging minimap2 primary alignments"

   samtools merge $DISCOANT/$GENE/minimap2/$GENE_merged.bam $DISCOANT/$GENE/minimap2/*_pri_sorted.bam 
   samtools sort $DISCOANT/$GENE/minimap2/$GENE_merged.bam > $DISCOANT/$GENE/minimap2/$GENE_merged_sorted.bam
   samtools index $DISCOANT/$GENE/minimap2/$GENE_merged_sorted.bam 

echo "Correecting the merged primary alignments with TranscriptClean"

   python $PROGRAMS/TranscriptClean/TranscriptClean.py \
   --sam $DISCOANT/$GENE/minimap2/$GENE_merged.sam \
   --primaryOnly --genome $REF_HG38/GRCh38.p13.genome_edit.fa \
   --outprefix $DISCOANT/$GENE/transcriptclean/$GENE

   samtools view -S -h -b $DISCOANT/$GENE/transcriptclean/$GENE_clean.sam | samtools sort - > $DISCOANT/$GENE/transcriptclean/$GENE_clean_sorted.bam

echo "Constructing transcripts based on the corrected alignments"

   stringtie -t $DISCOANT/$GENE/minimap2/$GENE_merged_sorted.bam \
   -G $REF_HG38/gencode.v35.annotation.gtf \
   -A $DISCOANT/$GENE/stringtie/$GENE_gene_abund.tab \
   -C $DISCOANT/$GENE/stringtie/$GENE_cov_ref.gtf \
   -o $DISCOANT/$GENE/stringtie/$GENE_STRINGTIE.gtf

echo "Annotating the transcripts with SQANTI3"

   python $PROGRAMS/SQANTI3/sqanti3_qc.py \
   $DISCOANT/$GENE/stringtie/$GENE_clean_STRINGTIE.gtf \
   $REF_HG38/gencode.v35.annotation.gtf $REF_HG38/GRCh38.p13.genome_edit.fa \
   --cage_peak $REF_HG38/refTSS_v3.1_human_coordinate.hg38.bed \
   --polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/polyA.list \
   -o $DISCOANT/$GENE/sqanti/$GENE_clean

echo "Creating a metagene FASTA"
   ## check if the strand variable is assigned correctly

   cat $DISCOANT/$GENE/stringtie/$GENE_clean_STRINGTIE_edited.gtf | awk '{ if ($3 == "exon") { print } }' | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5, ".", ".", "$STRAND"}' > $DISCOANT/$GENE/stringtie/$GENE_clean_all_exons.bed

   bedtools getfasta -s -fi $REF_HG38/GRCh38.p13.genome_edit.fa -bed $DISCOANT/$GENE/stringtie/$GENE_clean_all_exons.bed -fo $DISCOANT/$GENE/stringtie/$GENE_clean_meta_gene_exons.fa


   echo ">meta_gene_$$GENE" > $DISCOANT/$GENE/stringtie/$GENE_clean_meta_gene.fa && \
   cat $DISCOANT/$GENE/stringtie/$GENE_clean_meta_gene_exons.fa | grep -v "^>" | tr -d '\n' >> $DISCOANT/$GENE/stringtie/$GENE_clean_meta_gene.fa && \
   echo >> $DISCOANT/$GENE/stringtie/$GENE_clean_meta_gene.fa 

echo "Mapping sample FASTA to the metagene"

   for filename in $DISCOANT/$GENE/fasta/downsample_500/*.fa

   do
   base=$(basename $filename .fa)
   echo "On sample : $base" 

   minimap2 -ax splice $DISCOANT/$GENE/stringtie/$GENE_clean_meta_gene.fa $DISCOANT/$GENE/fasta/downsample_500/${base}.fa > $DISCOANT/$GENE/minimap2_metagene/${base}.sam 

   done

