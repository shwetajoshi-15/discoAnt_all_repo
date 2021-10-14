#/bin/bash!

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


source discoAnt_params_PLCG1_R1.txt
export PYTHONPATH="$PROGRAMS/cDNA_Cupcake/sequence/:$PYTHONPATH"

echo "Making folders"

mkdir -p $DISCOANT/$GENE
mkdir -p $DISCOANT/$GENE/fasta_stats
mkdir -p $DISCOANT/$GENE/minimap2
mkdir -p $DISCOANT/$GENE/minimap2_target
mkdir -p $DISCOANT/$GENE/flair
mkdir -p $DISCOANT/$GENE/sqanti3
mkdir -p $DISCOANT/$GENE/flair_metagene_minimap2
mkdir -p $DISCOANT/$GENE/flair_metagene_salmon


##########                                          ##########
##########  1. Quality control for the fasta files  ##########
##########                                          ##########

echo "Extracting the number of reads and read lengths"

#basename -s .fa $FASTA/*.fa | sed 's/^.*bar/bar/g' > $DISCOANT/$GENE/fasta_stats/barcodes.txt

#file=$DISCOANT/$GENE/fasta_stats/$GENE_stats.COMPLETED

#if [[ ! -f $DISCOANT/$GENE/fasta_stats/"$GENE"_stats.COMPLETED ]]
#then

#        for filename in $FASTA/*.fa
#        do
#        base=$(basename $filename .fa)
#        echo "On sample : $base"

#        cat $FASTA/${base}.fa | grep -c "^>" >> $DISCOANT/$GENE/fasta_stats/"$GENE"_tmp1.txt
#        paste $DISCOANT/$GENE/fasta_stats/barcodes.txt $DISCOANT/$GENE/fasta_stats/"$GENE"_tmp1.txt | column -s $'\t' -t > $DISCOANT/$GENE/fasta_stats/"$GENE"_total_number_of_reads.txt

#        mkdir -p $DISCOANT/$GENE/fasta_stats/"$GENE"_read_lengths
#        bioawk -c fastx '{print length($seq)}' $FASTA/${base}.fa > $DISCOANT/$GENE/fasta_stats/"$GENE"_read_lengths/${base}.txt

#        awk '{/>/&&++a||b+=length()}END{print b/a}' $FASTA/${base}.fa >> $DISCOANT/$GENE/fasta_stats/"$GENE"_tmp2.txt
#        paste $DISCOANT/$GENE/fasta_stats/barcodes.txt $DISCOANT/$GENE/fasta_stats/"$GENE"_tmp2.txt | column -s $'\t' -t > $DISCOANT/$GENE/fasta_stats/"$GENE"_mean_read_lengths.txt

#rm $DISCOANT/$GENE/fasta_stats/*tmp* && rm $DISCOANT/$GENE/fasta_stats/barcodes.txt && touch $DISCOANT/$GENE/fasta_stats/"$GENE"_stats.COMPLETED
#        done

#else
#echo "FASTA stats are present in $DISCOANT/$GENE/fasta_stats"
#fi


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

echo "FLAIR - correct and collapse"

cat $FASTA/*.fa > /g/data/tv28/$GENE/all.fa

python $PROGRAMS/flair/bin/bam2Bed12.py -i $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam > $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bed

python $PROGRAMS/flair/flair.py correct -q $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bed -c $REF_HG38/chrom_sizes.txt -f $REF_HG38/gencode.v35.annotation.gtf -g $REF_HG38/GRCh38.p13.genome_edit.fa -o $DISCOANT/$GENE/flair/$GENE

python $PROGRAMS/flair/flair.py collapse -r /g/data/tv28/$GENE/all.fa -q $DISCOANT/$GENE/flair/"$GENE"_all_corrected.psl -f $REF_HG38/gencode.v35.annotation.gtf -g $REF_HG38/GRCh38.p13.genome_edit.fa -o $DISCOANT/$GENE/flair/"$GENE"_collapse_500 --temp_dir $DISCOANT/$GENE/flair -s 500

python $PROGRAMS/flair/flair.py quantify -r $DISCOANT/$GENE/reads_manifest_"$GENE".tsv -i $DISCOANT/$GENE/flair/"$GENE"_collapse_500.isoforms.fa -o $DISCOANT/$GENE/flair/"$GENE"_quant --tpm

echo "Annotated transcripts are present in $DISCOANT/$GENE/sqanti3"

if [[ ! -f $DISCOANT/$GENE/sqanti3/"$GENE"_sqanti3.COMPLETED ]]
then
        python $PROGRAMS/SQANTI3-1.3/sqanti3_qc.py \
        --gtf $DISCOANT/$GENE/flair/"$GENE"_collapse_500.isoforms.gtf \
        $REF_HG38/gencode.v35.annotation.gtf $REF_HG38/GRCh38.p13.genome_edit.fa \
        --cage_peak $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed \
        --polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/human.polyA.list.txt \
        -d $DISCOANT/$GENE/sqanti3 -o "$GENE"

        touch $DISCOANT/$GENE/sqanti3/"$GENE"_sqanti3.COMPLETED

else
echo "Annotated transcripts are present in $DISCOANT/$GENE/sqanti3"
fi


for filename in $DISCOANT/$GENE/minimap2/*_pri_sorted.bam
        do
        base=$(basename $filename _pri_sorted.bam)
        echo "On sample : $base"

        minimap2 -ax map-ont $DISCOANT/$GENE/flair/"$GENE"_collapse_500.isoforms.fa $FASTA/${base}.fa > $DISCOANT/$GENE/flair_metagene_minimap2/${base}.sam
        samtools view -S -h -b $DISCOANT/$GENE/flair_metagene_minimap2/${base}.sam | samtools sort - > $DISCOANT/$GENE/flair_metagene_minimap2/${base}_sorted.bam
        samtools view -h -F 2308 $DISCOANT/$GENE/flair_metagene_minimap2/${base}_sorted.bam | samtools sort - > $DISCOANT/$GENE/flair_metagene_minimap2/${base}_pri_sorted.bam
        samtools index $DISCOANT/$GENE/flair_metagene_minimap2/${base}_pri_sorted.bam
        done

for filename in $DISCOANT/$GENE/minimap2/*_pri_sorted.bam
        do
        base=$(basename $filename _pri_sorted.bam)
        echo "On sample : $base"

        salmon quant -t $DISCOANT/$GENE/flair/"$GENE"_collapse_500.isoforms.fa -l A -a $DISCOANT/$GENE/flair_metagene_minimap2/${base}_pri_sorted.bam -o $DISCOANT/$GENE/flair_metagene_salmon/${base}

        done

