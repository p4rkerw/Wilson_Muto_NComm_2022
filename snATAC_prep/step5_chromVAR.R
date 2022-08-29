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
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log_step5.out Rscript $SCRATCH1/Wilson_Muto_NComm_2022/snATAC_prep/step5_chromVAR13.R

library(Signac) # 1.3.0
library(Seurat) # 4.0.3
library(JASPAR2020) # # 0.99.10
library(TFBSTools) # 1.39.0
library(BSgenome.Hsapiens.UCSC.hg38) # 1.4.3
library(patchwork) # 1.1.1
library(motifmatchr) # 1.14.0
library(here) # 1.0.1
library(chromVAR) # 1.14.0
library(future) # 1.21.0
library(openxlsx) # 4.2.4
set.seed(1234)

atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step4_m2anno.rds"))
DefaultAssay(atacAggr) <- "peaks"

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
atacAggr <- AddMotifs(
  object = atacAggr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

atacAggr <- RegionStats(
  object = atacAggr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c("-", "-")
)

# compute motif activities using chromvar
atacAggr <- RunChromVAR(
  object = atacAggr,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(atacAggr, here(atac_aggr_prep,"step5_chromVAR.rds"), compress=FALSE)
