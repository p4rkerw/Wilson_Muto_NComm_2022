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
# bsub -Is -G compute-parkerw -R 'rusage[mem=64GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

library(here)
library(Signac)
library(stringr)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(tibble)
library(ComplexHeatmap)
library(Seurat)

figures <- here("project","analysis","dkd","figures")
dir.create(here(figures), showWarnings=FALSE)

#################################################
# supplemental figure number peaks by cell type
files <- list.files(here("project","analysis","dkd","atac_aggr_prep","macs2"), pattern = "narrowPeak", full.names = TRUE)
idents <- str_split(basename(files), pattern="_", simplify = TRUE)[,1]

peaks <- lapply(seq(files), function(index){
  file <- files[index]
  ident <- idents[index]
  # qval FDR < 0.05
  df <- fread(file) %>%
    dplyr::filter(V10 > 1.3) %>% 
    as.data.frame()
  df$ident <- ident
  return(df)
}) %>% bind_rows()

atacAggr <- readRDS(here("project","analysis","dkd","atac_aggr_prep","step5_chromVAR.rds"))
num_cells <- table(atacAggr@meta.data$celltype) %>% as.data.frame()
colnames(num_cells)[1] <- "celltype"
colnames(num_cells)[2] <- "num_cells"

# change TCELL and BCELL to LYMPH
# and change FIB_VSMC_MC to MESFIB
# remove PT_CD36
metadata <- atacAggr@meta.data %>% as.data.frame()
metadata <- dplyr::mutate(metadata, celltype_update = ifelse(celltype == "TCELL" | celltype == "BCELL","LYMPH",as.character(celltype)))
metadata <- dplyr::filter(metadata, celltype != "PT_CD36")
num_cells <- table(metadata$celltype_update) %>% 
  as.data.frame() %>%
  rename(celltype = Var1)

num_peaks <- table(peaks$ident) %>%
  as.data.frame()
colnames(num_peaks)[1] <- "celltype"
colnames(num_peaks)[2] <- "num_peaks"
num_peaks <- dplyr::mutate(num_peaks, celltype = ifelse(celltype == "MESFIB","FIB_VSMC_MC",as.character(celltype)))
num_peaks <- dplyr::mutate(num_peaks, celltype = ifelse(celltype == "PTVCAM1","PT_VCAM1",as.character(celltype)))

toplot <- merge(num_peaks, num_cells, by = "celltype")

levels <- c("PCT","PST","PT_VCAM1","PEC","ATL",
            "TAL1","TAL2","DCT1","DCT2","PC",
            "ICA","ICB","ENDO","PODO","FIB_VSMC_MC",
            "LYMPH","MONO")
toplot$celltype <- factor(toplot$celltype, levels=levels)
pdf(here(figures,"sfigure1.pdf"), width=10, height=6)
ggplot(toplot, aes(x=Freq, y=num_peaks, color=celltype, label=celltype)) +
  geom_point() +
  geom_text_repel(aes(label=celltype)) +
  theme_bw() +
  xlab("Number of Cells") +
  ylab("Number of Peaks")
dev.off()
###########################################################
# supplemental figure dotplot of atacAggr gene activity and imputed RNA markers
levels(atacAggr) <- rev(levels(atacAggr))
features <- c("SLC34A1","LRP2","SLC5A1","SLC5A2","HAVCR1",
              "VCAM1","CD36","CFH","SLC12A1","CLDN16","SLC12A3",
              "AQP2","SLC26A7","SLC26A4","NPHS2","FLT1",
              "PIEZO2","COL1A2","PTPRC","CD3E","MS4A1","CSF1R")
p1 <- DotPlot(atacAggr, assay="RNA", features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #820x600
p2 <- DotPlot(atacAggr, assay="IMPRNA", features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #820x600
pdf(here(figures,"sfigure2.pdf"), width=10, height=6)
list(p1,p2)
dev.off()
###########################################################3
# supplemental figure dotplot of rnaAggr gene expression
rnaAggr <- readRDS(here("project","analysis","dkd","rna_aggr_prep","step2_anno.rds"))
levels(rnaAggr) <- rev(levels(rnaAggr))
features <- c("SLC34A1","LRP2","SLC5A1","SLC5A2","HAVCR1",
              "VCAM1","CFH","SLC12A1","CLDN16","SLC12A3",
              "AQP2","SLC26A7","SLC26A4","NPHS2","FLT1",
              "PIEZO2","COL1A2","PTPRC","CD3E","MS4A1","CSF1R")
p1 <- DotPlot(rnaAggr, assay="SCT", features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #820x600
pdf(here(figures,"sfigure3.pdf"), width=10, height=6)
print(p1)
dev.off()
#################################################################
# draw ATP1B1 gene model peak coverage
# edit afterwards in inskcape / adobe illustrator etc.
library(Signac)
library(Seurat)
library(here)
figures <- here("project","analysis","dkd","figures")
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan.rds"))
DefaultAssay(atacAggr) <- "peaks"

# this will print the coverage plot by itself
Idents(atacAggr) <- paste0(atacAggr@meta.data$celltype,"_",atacAggr@meta.data$diabetes)
p1 <- CoveragePlot(atacAggr, ident=c("PCT_0","PCT_1"), region = "ATP1B1", peaks=FALSE, links=FALSE)
pdf(here(figures, "sfigure_ATP1B1_coverage.pdf"))
print(p1)
dev.off()

# # add multiple tracks together (have to use genomic coordinates rather than gene name here)
p2 <- LinkPlot(atacAggr, region="chr1-169106683-169135009", min.cutoff=0.4)
p3 <- CombineTracks(list(p2, p3))

# this will print coverage plot with the links
pdf(here(figures, "sfigure_ATP1B1.pdf"))
print(p7)
dev.off()



