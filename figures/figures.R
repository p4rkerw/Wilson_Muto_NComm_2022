#!/usr/bin/env Rscript
# to run locally:
SCRATCH1=/mnt/g/scratch
docker run -it \
--workdir $HOME \
-v /mnt/g/diabneph:$HOME/project \
-v $HOME:$HOME \
-v $SCRATCH1:$SCRATCH1 \
-e SCRATCH1="/mnt/g/scratch" \
p4rkerw/sctools:R4.1.0

# to run interactively on the RIS compute1 cluster
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=256GB]' -n 10 -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

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

#######################
# figure snATAC-seq UMAP
# read in aggregated snATAC object
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan13.rds"))
Idents(atacAggr) <- "celltype"
DefaultAssay(atacAggr) <- "peaks"
p1 <- DimPlot(atacAggr, label=FALSE) + NoLegend()
png(here(figures,"figure1A.png"), height=600, width=600)
print(p1)
dev.off()

# figure snATAC-seq Gene Activity
marker.genes <- c("SLC34A1","SLC5A1","SLC5A2","VCAM1", # PT and PT-VCAM1+ markers
                  "CD36", # PT_CD36
                  "CFH", # PEC
                  "S100A2", #ATL
                  "SLC12A1", # TAL NKCC2
                  "CLDN10", #MTAL
                  "CLDN16", #CTAL
                  "SLC12A3", # DCT1 and DCT2 NCC
                  "AQP2", # PC
                  "SLC4A1", # ICA
                  "SLC26A4", # ICB
                  "NPHS1",# PODO
                  "FLT1", # ENDO
                  "PDGFRB", # VSMC, MES, JGA,
                  "PTPRC","CD3E","MS4A1", # Lymphocytes
                  "FCGR3A","CD14","CSF1R") # Monocyte / Macrophage
# update levels
levels(atacAggr) <- rev(c("PCT","PST","PT_VCAM1","PT_CD36","PEC",
  "ATL","TAL1","TAL2","DCT1","DCT2",
  "PC","ICA","ICB","PODO","ENDO",
  "FIB_VSMC_MC","TCELL","BCELL","MONO"))
                  
DefaultAssay(atacAggr) <- "RNA"
p2 <- DotPlot(atacAggr, features=marker.genes, cols=c("lightgrey","darkred")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
png(here(figures,"figure1B.png"), height=600, width=800)
print(p2)
dev.off()


####################################################
# figure upset plot cell-specific DAR intersection for dkd vs. control
file <- here("analysis","combined_adv","markers","dar.macs2.celltype.markers.xlsx")
idents <- getSheetNames(file)

dar.ls <- lapply(levels(atacAggr), function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    rownames_to_column(var = "peak")
  list <- df$peak
  return(list)
})
names(dar.ls) <- levels(atacAggr)

# reorder
dar.ls <- dar.ls[rev(levels(atacAggr))]

# remove non-kidney cell types
id_remove <- names(dar.ls) %in% c("BCELL","TCELL","MONO")
dar.ls <- dar.ls[!id_remove]

mat <- list_to_matrix(dar.ls)

# cell-specific dar shared across cell types
m1 <- make_comb_mat(mat)
# filter by intersection size and set size
# only include intersections with a size of greater than 20, but include individual cell types
# identify combinations with only a single ident by summing all the occurrences of 1 in the string
df <- strsplit(names(comb_size(m1)), split="") %>%
  as.data.frame() 
df <- sapply(df, as.numeric)
colsums <- colSums(df)

m1 <- m1[comb_size(m1) > 5 | colsums == 1]
# m1 <- m1[colsums == 1]
# m1 <- m1[comb_size(m1) > 10]

ht <- draw(UpSet(m1, set_order=names(dar.ls)))
od <- column_order(ht)
cs <- comb_size(m1)
p3 <- draw(UpSet(m1, set_order=names(dar.ls)))
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})

# find peaks in all 4 cell types
# all4_peaks <- extract_comb(m1, "1111")
# dar <- lapply(idents, function(ident){
#   df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
#     dplyr::filter(p_val_adj < 0.05) %>%
#     rownames_to_column(var = "peak")
#   df$celltype <- ident
#   return(df)
# }) %>% bind_rows()
# all4 <- dar[dar$peak %in% all4_peaks, ]

######################################
# plot panther db GO enrichment for DAR in PCT diabetes vs. control
panther <- fread(here("analysis","dkd","panther","panther_pct_dar.macs2.diab_vs_ctrl.txt"), skip=11, header=TRUE) #skip white space and database annotation
barplot(panther$`upload_1 (fold Enrichment)`, xlab="GO Biological Process", ylab = "Pathway Fold Enrichment", main="Pathway enrichment for PCT DAR in diabetes")
######################################
# figure Upset snATAC-seq DAR diabetes group comparison
files <- c("dar.celltype.early_vs_ctrl.xlsx",
           "dar.celltype.adv_vs_ctrl.xlsx",
           "dar.celltype.adv_vs_early.xlsx")

return_dar <- function(file) {
  file <- here("analysis","combined_adv","markers", file)
  idents <- getSheetNames(file)
  dar <- lapply(idents, 
                function(ident){
                                df <- read.xlsx(file, sheet = ident, rowNames = TRUE)
                                df <- rownames_to_column(df, var="peak")
                                df$ident <- ident
                                return(df)
          }) %>% bind_rows
  return(dar)
}

dar.adv_vs_ctrl <- return_dar("dar.celltype.adv_vs_ctrl.xlsx") %>% filter(p_val_adj < 0.05)
dar.early_vs_ctrl <- return_dar("dar.celltype.early_vs_ctrl.xlsx") %>% filter(p_val_adj < 0.05)
dar.adv_vs_early <- return_dar("dar.celltype.adv_vs_early.xlsx") %>% filter(p_val_adj < 0.05)

ls <- list(early_vs_ctrl=dar.early_vs_ctrl$peak,
           adv_vs_ctrl=dar.adv_vs_ctrl$peak,
           adv_vs_early=dar.adv_vs_early$peak)

mat <- list_to_matrix(ls)

m1 <- make_comb_mat(mat)
ht <- draw(UpSet(m1, set_order =  c("early_vs_ctrl","adv_vs_ctrl","adv_vs_early")))
od <- column_order(ht)
cs <- comb_size(m1)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})

# find DAR unique to early DKD
unique.early_vs_ctrl <- extract_comb(m1, "100")
genes <- dar.early_vs_ctrl[dar.early_vs_ctrl$peak %in% unique.early_vs_ctrl, ]$gene
write.csv(genes, "G:/downloads/genes.csv")

# find DAR unique to advanced DKD
unique.adv_vs_early <- extract_comb(m1, "001")
genes <- dar.adv_vs_early[dar.adv_vs_early$peak %in% unique.adv_vs_early, ]$gene
write.csv(genes, "G:/downloads/genes.csv")

# find DAR that decrease in early DKD and further decrease in advanced DKD
decreased.early_vs_ctrl <- dar.early_vs_ctrl[dar.early_vs_ctrl$avg_logFC < 0,]$peak
decreased.adv_vs_early <- dar.adv_vs_early[dar.adv_vs_early$avg_logFC < 0,]$peak

overlap <- intersect(decreased.early_vs_ctrl, decreased.adv_vs_early)
genes <- dar.early_vs_ctrl[dar.early_vs_ctrl$peak %in% overlap, ]$gene
write.csv(genes, "G:/downloads/genes.csv")

##########################################
return_deg <- function(file) {
  file <- here("analysis","combined_adv","markers", file)
  idents <- getSheetNames(file)
  dar <- lapply(idents, 
                function(ident){
                  df <- read.xlsx(file, sheet = ident, rowNames = TRUE)
                  df <- rownames_to_column(df, var="peak")
                  df$ident <- ident
                  return(df)
                }) %>% bind_rows
  return(dar)
}

files <- c("deg.celltype.early_vs_ctrl.xlsx",
           "deg.celltype.adv_vs_ctrl.xlsx",
           "deg.celltype.adv_vs_early.xlsx")

# snRNAseq analysis DEGs
# find degs upregulated in early DKD that are not present in the broader comparison
deg.early_vs_ctrl <- return_deg("deg.celltype.early_vs_ctrl.xlsx") %>% filter(p_val_adj < 0.05)
deg.adv_vs_ctrl <- return_deg("deg.celltype.adv_vs_ctrl.xlsx") %>% filter(p_val_adj < 0.05)
deg.adv_vs_early <- return_deg("deg.celltype.adv_vs_early.xlsx") %>% filter(p_val_adj < 0.05)









