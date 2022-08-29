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
#
# to run interactively on the RIS compute1 cluster:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

# to run detached:
# git clone https://github.com/p4rkerw/dkd $SCRATCH1/dkd
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log_step6_ccan.out Rscript $SCRATCH1/Wilson_Muto_NComm_2022/snATAC_prep/step6_find_ccan.R

library(Signac) # 1.3.0
library(Seurat) # 4.0.3
library(SeuratWrappers) # 0.3.0
library(ggplot2) # 3.3.5
library(patchwork) # 1.1.1
library(cicero) # 1.3.4.11
library(here) # 1.0.1

atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step5_chromVAR.rds"))
DefaultAssay(atacAggr) <- "peaks"

# convert to cds monocle/cicero format
atac.cds <- as.cell_data_set(atacAggr)
atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates=reducedDims(atac.cds)$UMAP)

# remove alt contigs
genomes <- seqlengths(atacAggr)[1:25]

# convert chr sizes to df
genome.df <- data.frame("chr" = names(genomes), "length" = genomes)

# run cicero
conns <- run_cicero(atac.cicero, genomic_coords = genome.df, sample_num=100)

# generate ccan
ccans <- generate_ccans(conns)

# add to seurat obj
links <- ConnectionsToLinks(conns=conns, ccans=ccans)
Links(atacAggr) <- links

saveRDS(atacAggr, file=here(atac_aggr_prep,"step6_ccan.rds"), compress=FALSE)
