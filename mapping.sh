#!/bin/bash

# Check if three arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 thread input output"
    exit 1
fi

thread=$1
input_bam=$2
output_dir=$3

ref_dir="/home/czlab/Downloads/reference_genome/hg38"
ref_genome="${ref_dir}/hg38.fa"
ref_genes_bed="${ref_dir}/gencode.v44.annotation.bed"
ref_genes_gtf="${ref_dir}/gencode.v44.annotation.gtf"
ref_chrom_sizes="${ref_dir}/chr_sizes"

# prepare fastq file for minimap2
tmp_fastq="/tmp/tmp.fastq"
mapped_bam="${output_dir}/mapped.bam"
primary_align="${output_dir}/primary_alignments.bam"
spt_align="${output_dir}/supplementary_alignments.bam"
primary_bed="${output_dir}/primary_alignments.bed"

samtools fastq -@ $thread $input_bam > $tmp_fastq

# Run minimap2 alignment
minimap2 -ax splice -uf --secondary=no --MD -t $thread --junc-bed $ref_genes_bed $ref_genome $tmp_fastq \
    | samtools view -b --no-PG -t $ref_chrom_sizes - \
    | samtools sort -@ 2 --no-PG - > $mapped_bam

rm $tmp_fastq

# Index the sorted BAM file
samtools -b -F 0x800 $mapped_bam > $primary_align
samtools -b -f 0x800 $mapped_bam > $spt_align
samtools index $primary_align
bedtools bamtobed -i $primary_align > $primary_bed

# assign genes to reads
featureCounts -a $ref_genes_gtf -L -o "${output_dir}/gene_assigned" -R CORE -g gene_name $mapped_bam


#featurecounts reports error for extreme long reads output in bam format
#gene_assign_tmp="${mapped_bam}.featureCounts"
#gene_assign="${output_dir}/gene_assign.bam"
#gene_tsv="${output_dir}/gene_assign.tsv"
#featureCounts -a $ref_genes_gtf -L -o "${output_dir}/gene_assigned" -R BAM -g gene_name $mapped_bam
#samtools sort -bS -@ 2 --no-PG $gene_assign_tmp > $gene_assign
#samtools index $gene_assign
#rm $gene_assign_tmp
#samtools view $gene_assign | \
#awk 'BEGIN{OFS="\t"; print "read_id", "gene"} \
#{geneName="-"; if($NF ~ /^XT:Z:/){split($NF, arr, ":"); \
#geneName=arr[3]}; print $1, geneName}' > $gene_tsv


