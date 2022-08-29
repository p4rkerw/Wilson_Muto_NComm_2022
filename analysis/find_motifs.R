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
# bsub -G compute-parkerw -R 'rusage[mem=256GB]' -n 10 -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log.find_chromvar.out Rscript $SCRATCH1/Wilson_Muto_NComm_2022/analysis/find_motifs.R

library(Signac) # 1.3.0
library(Seurat) # 4.0.3
library(JASPAR2020) # 0.99.10
library(TFBSTools) # 1.30.0
library(BSgenome.Hsapiens.UCSC.hg38) # 1.4.3
library(patchwork) # 1.1.1
library(here) # 1.0.1
library(future) # 1.21.0
library(openxlsx) # 4.2.4
library(dplyr) # 1.0.7
set.seed(1234)

# identify enriched motifs using a cell type ident and list of dar peaks from an xlsx file
GetEnrichedMotifs <- function(ident, seurat_aggregate, dar_file) { 
  print(paste0("Finding enriched motifs in dar for: ", ident))
  # read in previously computed dar and filter by pval
  dar <- read.xlsx(here(markers, dar_file), sheet = ident)
  colnames(dar)[1] <- "peak"
  
  # order sig dar by pval
    sig_peaks <- dar[dar$p_val < 0.05, ] %>%
    arrange(-desc(p_val))
    sig_peaks <- sig_peaks$peak

    # find background peaks
    open.peaks <- AccessiblePeaks(seurat_aggregate, idents = ident)

    # eliminate any sig peaks that do not meet 'open peak' criteria as defined by AccessiblePeaks()
    # these are peaks present in a small fraction of cells
    sig_peaks <- sig_peaks[sig_peaks %in% open.peaks]

    meta.feature <- GetAssayData(seurat_aggregate, assay = "peaks", slot = "meta.features")
    peaks.matched <- tryCatch(MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[sig_peaks, ],
    n = 50000
    ), error = function(e) return(NULL)) 

    enriched.motifs <- tryCatch(FindMotifs(seurat_aggregate, features=sig_peaks, background = peaks.matched), error = function(e) return(NULL)) 

    # sort by pvalue
    enriched.motifs <- dplyr::arrange(enriched.motifs, desc(-pvalue))
  return(enriched.motifs)
}

GetChromvarActivities <- function(ident, seurat_aggregate, motif) {
  print(paste0("Finding chromVAR activities for: ", ident))
  DefaultAssay(seurat_aggregate) <- 'chromvar'
  chromvar <- FindMarkers(seurat_aggregate, 
                     ident.1 = ident,
                     test.use = 'LR',
                     latent.vars = "macs2_FRiP",
                     logfc.threshold = 0,
                     only.pos=TRUE) # find all cluster-specific chromVAR motif activities
  motif_lookup <- rownames(chromvar)
  motif_names <- sapply(motif_lookup, function(x) motif[[x]])
  return(cbind(chromvar, gene = motif_names))
}

CompareChromvarActivities <- function(ident, seurat_aggregate, test, ref, meta_group, motif) {
  print(paste0("Finding chromVAR activities for: ",ident))
  seurat_aggregate <- seurat_aggregate
  DefaultAssay(seurat_aggregate) <- "chromvar"
  
  group <- dplyr::select(seurat_aggregate@meta.data, all_of(meta_group))[,1]
  seurat_aggregate@meta.data$celltype.stim <- paste0(seurat_aggregate@meta.data$celltype,"_", group)
  Idents(seurat_aggregate) <- "celltype.stim"
  chromvar <- tryCatch(FindMarkers(seurat_aggregate, 
                     ident.1 = paste0(ident,"_", test),  
                     ident.2 = paste0(ident, "_", ref),
                     test.use = 'LR', 
                     latent.vars = "macs2_FRiP",
                     logfc.threshold = 0.10),
                     error = function(e) NULL) # find all celltype specific motif activities that differ between groups

  motif_lookup <- rownames(chromvar)
  motif_names <- sapply(motif_lookup, function(x) motif[[x]])
  return(cbind(chromvar, gene = motif_names))
}


###########################################################################

# enable parallel processing
plan("multiprocess", workers = 10) 
options(future.globals.maxSize = 256000 * 1024^2) # for 256 Gb RAM
plan()

# read in aggregated snATAC object
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan.rds"))
Idents(atacAggr) <- "celltype"
DefaultAssay(atacAggr) <- "peaks"

# define markers directory previously computed files with find_dar.R
markers <- here("project","analysis","dkd","markers")

# assign idents
idents <- levels(atacAggr)

motif.names <- GetMotifData(object = atacAggr, assay = "peaks", slot = "motif.names")

# find cell-specific dar motifs
list.celltype.motifs <- lapply(idents, function(x) GetEnrichedMotifs(x, seurat_aggregate = atacAggr, dar_file = "dar.macs2.celltype.markers.xlsx"))
write.xlsx(list.celltype.motifs, file = here(markers,"motifs_in_dar.macs2.celltype.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# find diabetes cell-specific dar motif enrichment
list.cluster.enriched_motif <- lapply(idents, function(x) GetEnrichedMotifs(x, seurat_aggregate = atacAggr, dar_file = "dar.macs2.celltype.diab_vs_ctrl.xlsx"))
write.xlsx(list.cluster.enriched_motif, file = here(markers,"motifs_in_dar.macs2.celltype.diab_vs_ctrl.xlsx"), sheetName = idents, rowNames = F, overwrite=TRUE)
                               
# cell-specific chromVAR activities                               
list.cluster.chromvar <- lapply(idents, function(x) GetChromvarActivities(x, seurat_aggregate = atacAggr, motif=motif.names))
write.xlsx(list.cluster.chromvar, file = here(markers,"chromvar.macs2.celltype.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

# cell-specific chromVAR activities for diabetes vs. control
list.cluster.dam <- lapply(idents, function(x) CompareChromvarActivities(x, seurat_aggregate = atacAggr, test=1, ref=0, meta_group="diabetes", motif=motif.names))
write.xlsx(list.cluster.dam, file = here(markers,"chromvar.macs2.celltype.diab_vs_ctrl.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

