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
# bsub -Is -G compute-parkerw -R 'rusage[mem=256GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

# to run detached:
# git clone https://github.com/p4rkerw/dkd $SCRATCH1/dkd
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -R 'rusage[mem=256GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log.find_footprint.out Rscript $SCRATCH1/dkd/analysis/find_footprint.R

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
library(stringr)
library(ggplot2)
set.seed(1234)

# read in aggregated snATAC object
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan.rds"))
Idents(atacAggr) <- "celltype"
DefaultAssay(atacAggr) <- "peaks"

# define markers directory previously computed files with find_dar.R
markers <- here("project","analysis","dkd","markers")

# create a footprints directory
footprints <- here("project","analysis","dkd","footprints")
dir.create(here(footprints), showWarnings=FALSE)

# find footprints in peaks
atacAggr <- Footprint(atacAggr, genome = BSgenome.Hsapiens.UCSC.hg38, motif.name = c("NR3C1","NR3C2","RELA","HNF4A"), in.peaks = TRUE)
p1 <- PlotFootprint(atacAggr, features = "NR3C1") + ggtitle("NR3C1 footprint in all cell types")
p2 <- PlotFootprint(atacAggr, features = "NR3C1", idents = "PCT") + ggtitle("NR3C1 footprint in PCT")
p3 <- PlotFootprint(atacAggr, features = "NR3C1", group.by = "diabetes") + ggtitle("NR3C1 footprint in all cell types with and without T2DM")

saveRDS(atacAggr, here(atac_aggr_prep,"step7_footprint.rds"), compress=FALSE)

# idents argument doesn't work in PlotFootprint when used with grouping var
atacAggr@meta.data$celltype.stim <- paste0(atacAggr@meta.data$celltype,"_", atacAggr@meta.data$diabetes)
Idents(atacAggr) <- "celltype.stim"
p4 <- PlotFootprint(atacAggr, features = "NR3C1", idents = c("PCT_0","PCT_1")) + ggtitle("NR3C1 PCT footprint with and without T2DM")
p5 <- PlotFootprint(atacAggr, features = "NR3C2", idents = c("PCT_0","PCT_1")) + ggtitle("NR3C2 PCT footprint with and without T2DM")
p6 <- PlotFootprint(atacAggr, features = "RELA", idents = c("PCT_0","PCT_1")) + ggtitle("RELA PCT footprint with and without T2DM")
p7 <- PlotFootprint(atacAggr, features = "HNF4A", idents = c("PCT_0","PCT_1")) + ggtitle("HNF4A PCT footprint with and without T2DM")

p8 <- PlotFootprint(atacAggr, features = "NR3C1", idents = c("DCT1_0","DCT1_1")) + ggtitle("NR3C1 DCT1 footprint with and without T2DM")
p9 <- PlotFootprint(atacAggr, features = "NR3C2", idents = c("DCT1_0","DCT1_1")) + ggtitle("NR3C2 DCT1 footprint with and without T2DM")

p10 <- PlotFootprint(atacAggr, features = "NR3C1", idents = c("PC_0","PC_1")) + ggtitle("NR3C1 PC footprint with and without T2DM")
p11 <- PlotFootprint(atacAggr, features = "NR3C2", idents = c("PC_0","PC_1")) + ggtitle("NR3C2 PC footprint with and without T2DM")

pdf(here(footprints, "footprint.inpeaks.pdf"))
print(list(p1,p2,p3,p4,p5,p6,p7, p8,p9,p10,p11))
dev.off()
