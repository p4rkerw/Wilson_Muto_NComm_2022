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
# bsub -Is -G compute-parkerw -R 'rusage[mem=256GB]' -n 10 -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

# to run detached:
# git clone https://github.com/p4rkerw/dkd $SCRATCH1/dkd
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -R 'rusage[mem=256GB]' -n 10 -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log.find_deg.out Rscript $SCRATCH1/Wilson_Muto_NComm_2022/analysis/find_deg.R

library(Seurat) 
library(EnsDb.Hsapiens.v86)
library(openxlsx) 
library(here) #
library(future) # 

# wrapper functions for FindMarkers 
CelltypeMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding DEG for: ",cluster))
  deg <- FindMarkers(seurat_aggregate, 
                     ident.1 = cluster,    
                     min.pct = 0.2) # find all cluster-specific degs
  return(deg)
}

CompareMarkers <- function(cluster, seurat_aggregate, test, ref, meta_group) {
  print(paste0("Comparing DEG for: ",cluster))
  group <- dplyr::select(seurat_aggregate@meta.data, all_of(meta_group))[,1]
  seurat_aggregate@meta.data$celltype.stim <- paste0(seurat_aggregate@meta.data$celltype,"_", group)
  Idents(seurat_aggregate) <- "celltype.stim"
  deg <- tryCatch(FindMarkers(seurat_aggregate, 
                     ident.1 = paste0(cluster,"_", test),  
                     ident.2 = paste0(cluster, "_", ref),
                     logfc.threshold = 0),
                  error=function(e) NULL) # find all celltype specific degs that differ between disease and control
  return(deg)
}

# enable parallel processing
plan("multiprocess", workers = 10) 
options(future.globals.maxSize = 256000 * 1024^2) # for 256 Gb RAM
plan()

# load annotated snRNA object and normalize the RNA assay before finding deg
rnaAggr <- readRDS(here("project","analysis","dkd","rna_aggr_prep","step2_anno.rds"))
DefaultAssay(rnaAggr) <- "RNA"
rnaAggr <- NormalizeData(rnaAggr)
Idents(rnaAggr) <- "celltype"

# FindMarkers and write to an xlsx file with default parameters
markers <- here("project","analysis","dkd","markers")
dir.create(here(markers), showWarnings = FALSE)
idents <- levels(rnaAggr@meta.data$celltype)

list.disease.deg <- lapply(idents, function(x) {CompareMarkers(x, seurat_aggregate = rnaAggr, test=1, ref=0, meta_group="diabetes")})                          
write.xlsx(list.disease.deg, file = here(markers,"deg.celltype.diab_vs_ctrl.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

list.cluster.deg <- lapply(idents, function(x) {CelltypeMarkers(x, seurat_aggregate = rnaAggr)})
write.xlsx(list.cluster.deg, file = here(markers,"deg.celltype.markers.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)
                        
# identify deg between PT and PT_VCAM1
deg <- FindMarkers(rnaAggr, ident.1 = "PTVCAM1", ident.2 = "PT", logfc.threshold = 0)
write.xlsx(deg, file = here(markers,"deg.PT_vs_PTVCAM1.markers.xlsx"), rowNames = T, overwrite=TRUE)

