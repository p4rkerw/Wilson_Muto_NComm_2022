#!/usr/bin/env Rscript
# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v /mnt/g/diabneph:$HOME/project \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.0

# to run interactively on the RIS compute1 cluster
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

# note: annotatr and org.Hs.eg.db bioconductor packagages need to be installed

library(openxlsx)
library(dplyr)
library(here)
library(tibble)
library(Signac)
library(stringr)
library(annotatr)
library(org.Hs.eg.db)
library(data.table)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


# read in cell-specific dar chipseeker annotations
dar <- read.csv(here("project","analysis","dkd","markers","dar_chipseeker_anno.csv"), row.names=1) %>%
  dplyr::filter(abs_log2FC > 0.1)
dar.gr <- StringToGRanges(dar$peak)
dar.gr$peakloc <- dar$peakloc

# load databases
db <- c("hg38_genes_promoters","hg38_enhancers_fantom")
anno <- build_annotations(genome = 'hg38', annotations=db)

# annotate peaks
dar_anno <- annotate_regions(dar.gr, annotations=anno, ignore.strand=TRUE)

# convert to df 
dar_anno.df <- as.data.frame(dar_anno)

# create a peak column for the df
# some peaks have multiple annotations
dar_anno.df <- dplyr::mutate(dar_anno.df, peak = paste0(seqnames, "-", start, "-", end))

# count peaks mapping to enhancers
enh <- dplyr::filter(dar_anno.df, annot.type == "hg38_enhancers_fantom")
n_distinct(enh$peak) # n=245

# count peaks mapping to promoters
prom <- dplyr::filter(dar_anno.df, annot.type == "hg38_genes_promoters")
n_distinct(prom$peak) # n=496

# calculate mean distance of enh to nearest TSS
# annotate the list of GRanges DAR for each cell type
enh_dist <- left_join(data.frame(peak=enh$peak), dar, by="peak") %>%
  distinct(peak, dist, peakloc)
