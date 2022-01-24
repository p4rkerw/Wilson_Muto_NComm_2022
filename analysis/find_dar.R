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
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

# to run detached:
# git clone https://github.com/p4rkerw/dkd $SCRATCH1/dkd
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log.find_dar.out Rscript $SCRATCH1/dkd/analysis/find_dar.R

library(openxlsx) # 4.1.4
library(Signac) # 1.3.0
library(Seurat) # 4.0.3
library(EnsDb.Hsapiens.v86)
library(dplyr) # 1.0.7
library(here) # 1.0.1
library(future) # 1.21.0
library(stringr)

# wrapper functions for FindMarkers 
CelltypeMarkers <- function(cluster, seurat_aggregate, anno) {
  print(paste0("Finding DAR for: ", cluster))
  dar <- FindMarkers(seurat_aggregate, 
                     ident.1 = cluster,    
                     test.use = 'LR', 
                     latent.vars = "macs2_FRiP",
                     logfc.threshold = 0,
                     min.pct = 0.2) # find all cluster-specific peaks that are present in at least 20% of specified cell type
  cf <- ClosestFeature(seurat_aggregate, regions=rownames(dar), annotation=anno)
  return(cbind(dar, gene=cf$gene_name, distance=cf$distance))
}

CompareMarkers <- function(cluster, seurat_aggregate, test, ref, meta_group, anno) {
  print(paste0("Comparing DAR for: ", cluster))
  group <- dplyr::select(seurat_aggregate@meta.data, all_of(meta_group))[,1]
  seurat_aggregate@meta.data$celltype.stim <- paste0(seurat_aggregate@meta.data$celltype,"_", group)
  Idents(seurat_aggregate) <- "celltype.stim"
  dar <- tryCatch(FindMarkers(seurat_aggregate, 
                     ident.1 = paste0(cluster,"_", test),  
                     ident.2 = paste0(cluster, "_", ref),
                     test.use = 'LR', 
                     latent.vars = "macs2_FRiP",
                     logfc.threshold = 0),
                     error = function(e) NULL) # find all celltype specific dars that differ between groups
  cf <- ClosestFeature(seurat_aggregate, regions=rownames(dar), annotation=anno)
  return(cbind(dar, gene=cf$gene_name, distance=cf$distance))
}

# enable parallel processing
plan("multiprocess", workers = 4) 
options(future.globals.maxSize = 128000 * 1024^2) # for 128 Gb RAM
plan()

# read in aggregated snATAC object
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan.rds"))
Idents(atacAggr) <- "celltype"
DefaultAssay(atacAggr) <- "peaks"

# set up annotation for ClosestFeature 
gene.ranges <- genes(EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(gene.ranges),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(gene.ranges) <- ucsc.levels
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

# create output directory
markers <- here("project","analysis","dkd","markers")
dir.create(here(markers), showWarnings = FALSE)
idents <- levels(atacAggr)

# compare celltype DAR for diabetes vs. control and diabetes groups
# diabetes: 1=true, 0=false
list.disease.dar <- lapply(idents, function(x) {CompareMarkers(x, seurat_aggregate = atacAggr, anno=gene.ranges, test=1, ref=0, meta_group="diabetes")})
write.xlsx(list.disease.dar, file = here(markers,"dar.macs2.celltype.diab_vs_ctrl.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# identify celltype DAR
plan("sequential")
list.celltype.dar <- lapply(idents, function(x) {CelltypeMarkers(x, seurat_aggregate = atacAggr, anno=gene.ranges)})
write.xlsx(list.celltype.dar, file = here(markers,"dar.macs2.celltype.markers.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# identify DAR between PCT and PT_VCAM1
dar <- FindMarkers(atacAggr, ident.1 = "PT_VCAM1", ident.2 = "PCT", logfc.threshold = 0)
cf <- ClosestFeature(atacAggr, regions=rownames(dar), annotation=gene.ranges)
df <- cbind(dar, gene=cf$gene_name, distance=cf$distance)
write.xlsx(df, file = here(markers,"dar.macs2.PCT_vs_PTVCAM1.markers.xlsx"), rowNames = T, overwrite=TRUE)
