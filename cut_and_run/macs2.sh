#!/bin/bash

#trimming
trim_galore --paired hTERT_IgG_rep1_R1.fastq.gz hTERT_IgG_rep1_R2.fastq.gz
trim_galore --paired hTERT_GR_rep1_R1.fastq.gz hTERT_GR_rep1_R2.fastq.gz


#hTERT_IgG_rep1
#mapping
bowtie2 --very-sensitive  -p 8 -x /home/users/ymuto/GRCh38_noalt_as/GRCh38_noalt_as --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -1 hTERT_IgG_rep1_R1_val_1.fq.gz -2 hTERT_IgG_rep1_R2_val_2.fq.gz -S hTERT_IgG_rep1.sam

#sort&index
samtools view -bS hTERT_IgG_rep1.sam -@ 8 | samtools sort > hTERT_IgG_rep1.bam
samtools index -@ 8 hTERT_IgG_rep1.bam

#Visualization
bamCoverage -b hTERT_IgG_rep1.bam -o hTERT_IgG_rep1.bw

#hTERT_GR_rep1
#mapping
bowtie2 --very-sensitive  -p 8 -x /home/users/ymuto/GRCh38_noalt_as/GRCh38_noalt_as --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -1 hTERT_GR_rep1_R1_val_1.fq.gz -2 hTERT_GR_rep1_R2_val_2.fq.gz -S hTERT_GR_rep1.sam

#sort&index
samtools view -bS hTERT_GR_rep1.sam -@ 8 | samtools sort > hTERT_GR_rep1.bam
samtools index -@ 8 hTERT_GR_rep1.bam

#Visualization
bamCoverage -b hTERT_GR_rep1.bam -o hTERT_GR_rep1.bw

#peak calling
macs2 callpeak -t hTERT_GR_rep1.bam -c hTERT_IgG_rep1.bam --name=hTERT_GR_rep1 --gsize=hs

#Intersection of peaks among replicates to generate consensus peak list
bedtools intersect -a hTERT_GR_rep1_peaks.narrowPeak -b hTERT_GR_rep2_peaks.narrowPeak > hTERT_GR_consensus.bed
