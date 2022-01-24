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
# to run interactively on the RIS compute1 cluster
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=256GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

# to run detached:
# git clone https://github.com/p4rkerw/dkd $SCRATCH1/dkd
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -R 'rusage[mem=256GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log.find_gene_enh.out Rscript $SCRATCH1/dkd/analysis/find_gene_enh.R

library(openxlsx) # 4.1.4
library(Signac) # 1.3.0
library(Seurat) # 4.0.3
library(EnsDb.Hsapiens.v86)
library(dplyr) # 1.0.7
library(here) # 1.0.1
library(future) # 1.21.0

# read in aggregated snATAC object
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan.rds"))
Idents(atacAggr) <- "celltype"
DefaultAssay(atacAggr) <- "peaks"

# create output directory
markers <- here("project","analysis","dkd","markers")
dir.create(here(markers), showWarnings = FALSE)
idents <- levels(atacAggr)

atacAggr <- LinkPeaks(atacAggr,
	peak.assay="peaks",
	expression.assay="IMPRNA",
	expression.slot="data",
	gene.coords=NULL,
	distance=200000,
	min.distance=NULL,
	min.cells=10,
	method="pearson",
	genes.use=NULL,
	n_sample=200,
	pvalue_cutoff = 0.05,
	score_cutoff = 0.05,
	verbose=TRUE)

saveRDS(atacAggr, here(atac_aggr_prep, "step7_links.rds"), compress=FALSE)
