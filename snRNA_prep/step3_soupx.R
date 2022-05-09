# this script will estimate the ambient RNA contamination fraction in aggregated snRNA dataset

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

install.packages("SoupX")
library(SoupX)
library(Seurat) # 4.0.3
library(ggplot2) # 3.3.5
library(harmony) # 0.1.0
library(DoubletFinder) # 2.0.3
library(dplyr) # 1.0.7
library(stringr) # 1.4.0
library(tibble) # 3.1.3
library(here) # 1.0.1
library(future) # 1.21.0
set.seed(1234)

# create output and plots directory
dir.create(here("project","analysis","dkd","rna_aggr_prep","plots"), recursive=TRUE, showWarnings=FALSE)
rna_aggr_prep <- here("project","analysis","dkd","rna_aggr_prep")

# read in annotated object
rnaAggr <- readRDS(here(rna_aggr_prep, "step2_anno.rds"))

# define cellranger-atac aggregation input dir
aggr_input_dir <- here("project","analysis","dkd","cellranger_rna_aggr","outs")

# read in raw and filtered counts
raw.matrix <- Read10X_h5(here(aggr_input_dir, "raw_feature_bc_matrix.h5"))
filt.matrix <- Read10X_h5(here(aggr_input_dir, "filtered_feature_bc_matrix.h5"))

# create a soup channel
soup.channel  <- SoupChannel(raw.matrix, filt.matrix)

# create a seurat object and cluster
srat  <- CreateSeuratObject(counts = filt.matrix)
srat    <- SCTransform(srat, verbose = TRUE)
srat    <- RunPCA(srat, verbose = TRUE)
srat    <- RunUMAP(srat, dims = 1:30, verbose = TRUE)
srat    <- FindNeighbors(srat, dims = 1:30, verbose = TRUE)
srat    <- FindClusters(srat, verbose = T)

# update soup channel with annotated clusters
meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
head(meta)

# estimate contamination
soup.channel  <- autoEstCont(soup.channel) # 0.04

# do the analysis a second time only including the annotated barcodes in the filtered matrix
filt.matrix <- filt.matrix[,colnames(filt.matrix) %in% rownames(rnaAggr@meta.data)]

# create a soup channel
soup.channel  <- SoupChannel(raw.matrix, filt.matrix)

# create a seurat object and cluster
srat  <- CreateSeuratObject(counts = filt.matrix)
srat    <- SCTransform(srat, verbose = TRUE)
srat    <- RunPCA(srat, verbose = TRUE)
srat    <- RunUMAP(srat, dims = 1:30, verbose = TRUE)
srat    <- FindNeighbors(srat, dims = 1:30, verbose = TRUE)
srat    <- FindClusters(srat, verbose = T)

# update soup channel with annotated clusters
meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
head(meta)

# estimate contamination
soup.channel  <- autoEstCont(soup.channel)




