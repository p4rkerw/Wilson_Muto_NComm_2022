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
rna_aggr_prep <- here("project","analysis","dkd","rna_aggr_prep")
rnaAggr <- readRDS(here(rna_aggr_prep,"step2_anno.rds"))
Idents(rnaAggr) <- "celltype"
p1 <- DimPlot(rnaAggr, label=FALSE) + NoLegend()
png(here(figures,"figure2A.png"), height=600, width=600)
print(p1)
dev.off()

levels(rnaAggr) <- rev(levels(rnaAggr))


# figure snATAC-seq Gene Activity
marker.genes <- c("SLC34A1","VCAM1", # PT and PT-VCAM1+ markers
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
                  "ACTA2",
                  "PTPRC") # Leukocytes

####################################################
# volcano plot of PT degs
# this file does not have a lfc threshold so all points can be displayed
pt.deg <- read.xlsx(here("project","analysis","dkd","markers","deg.celltype.diab_vs_ctrl.xlsx"), sheet = "PT", rowNames=TRUE) %>%
  rownames_to_column(var = "gene")
pt.deg <- dplyr::filter(pt.deg, p_val_adj < 0.05) %>%
  dplyr::mutate(color = ifelse(abs(avg_log2FC) > 0.25,1,0))

# modify outliers with absolute LFC > 1.5 and set to 1.5
# this will make volcano plot more narrow
# approx three points meet this criteria
pt.deg <- dplyr::mutate(pt.deg, update_log2FC = ifelse(avg_log2FC > 1.5, 1.5, avg_log2FC)) %>%
  dplyr::mutate(update_log2FC = ifelse(avg_log2FC < -1.5, -1.5, avg_log2FC))

pdf(here(figures, "figure2b.pdf"))
ggplot(pt.deg, aes(x=update_log2FC, y=-log10(p_val_adj), color=color)) +
 geom_point() +
 theme_bw() +
 guides(color="none") +
 geom_vline(xintercept=0.25) +
 geom_vline(xintercept=-0.25)
dev.off()

####################################################
# upset plot
# read in degs
file <- here("project","analysis","dkd","markers","deg.celltype.diab_vs_ctrl.xlsx")
idents <- getSheetNames(file)
deg.ls <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::filter(abs(avg_log2FC) > 0.25) %>%
    rownames_to_column(var = "gene")
  list <- df$gene
  return(list)
})
names(deg.ls) <- idents

# remove non-kidney cell types
id_remove <- names(deg.ls) %in% c("LEUK")
deg.ls <- deg.ls[!id_remove]

mat <- list_to_matrix(deg.ls)

# cell-specific dar shared across cell types
m1 <- make_comb_mat(mat)
# filter by intersection size and set size
# only include intersections with a size of greater than 20, but include individual cell types
# identify combinations with only a single ident by summing all the occurrences of 1 in the string
df <- strsplit(names(comb_size(m1)), split="") %>%
  as.data.frame() 
df <- sapply(df, as.numeric)
colsums <- colSums(df)

# keep intersection set size of at least 10 and retain all groups with at least 1 deg
m1 <- m1[comb_size(m1) > 10 | colsums == 1]

ht <- draw(UpSet(m1, set_order=names(deg.ls)))
od <- column_order(ht)
cs <- comb_size(m1)
p3 <- draw(UpSet(m1, set_order=names(deg.ls)))
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
}) # export as pdf 6x7" figure2c.pdf

######################################
# panther pathways for proximal tubule deg
file <- here("project","analysis","dkd","markers","deg.celltype.diab_vs_ctrl.xlsx")
idents <- c("PT")

deg.df <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
    dplyr::filter(abs_log2FC > 0.25) %>%
    rownames_to_column(var = "peak")
  df$celltype <- ident
  return(df)
}) %>% bind_rows()

genes <- unique(deg.df$gene) # n=607
# enter gene list into pantherdb at http://geneontology.org/ and export results using filename below

# plot panther db GO enrichment for DAR in PCT diabetes vs. control
panther <- fread(here("project","analysis","dkd","panther","panther_pt_deg.diab_vs_ctrl.txt"), skip=11, header=TRUE) #skip white space and database annotation

# remove any pathways with negative enrichment or NA values
panther$'upload_1 (fold Enrichment)' <- as.numeric(as.character(panther$'upload_1 (fold Enrichment)'))
panther <- na.omit(panther)

# barplot(panther$`upload_1 (fold Enrichment)`, xlab="GO Biological Process", ylab = "Pathway Fold Enrichment", main="Proximal tubule DEG in diabetes")
panther$order <- seq(nrow(panther))
panther$color <- ifelse(panther$order <=25, "Top25", "Not_Top25")

ggplot(panther, aes(x=order, y=`upload_1 (fold Enrichment)`, fill=color, color=color)) + 
  geom_bar(show.legend=FALSE,stat="identity") +
  ylab("Pathway fold enrichment") +
  xlab("GO Biological Process") +
  theme_bw() # export 6x6" pdf figure2d.pdf

######################################
# read in aggregated snATAC object
library(Seurat)
library(Signac)
library(here)
library(openxlsx)
library(dplyr)
library(tibble)

figures <- here("project","analysis","dkd","figures")
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan.rds"))
Idents(atacAggr) <- "celltype"
DefaultAssay(atacAggr) <- "peaks"
Idents(atacAggr) <- paste0(atacAggr@meta.data$celltype,"_",atacAggr@meta.data$diabetes)

# cell-specific DAR intersection for dkd vs.control
file <- here("project","analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx")
dar.pct <- read.xlsx(file, sheet="PCT", rowNames = TRUE) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  rownames_to_column(var = "peak")

dar.gr <- StringToGRanges(dar.pct$peak)
plot.gr <- StringToGRanges("chr20-57560110-57578121")
dar.gr <- join_overlap_intersect(dar.gr, plot.gr)

# no dmr intersect with this region

# pck1 with 1k upstream and 10k downstream flank
cp <- CoveragePlot(atacAggr, ident=c("PCT_0","PCT_1"), region = "chr20-57560110-57578121", peaks=TRUE, links=FALSE)
ccan.plot <- LinkPlot(atacAggr, region="chr20-57560110-57578121", min.cutoff=0.4)
dar.plot <- PeakPlot(atacAggr, region = "chr20-57560110-57578121", peaks=dar.gr)
plot <- CombineTracks(list(cp,dar.plot, ccan.plot))

pdf(here(figures, "figure2e.pdf"))
print(plot)
dev.off()

################
# draw gluconeogenesis genes vln plot
rna_aggr_prep <- here("project","analysis","dkd","rna_aggr_prep")
rnaAggr <- readRDS(here(rna_aggr_prep, "step2_anno.rds"))
pdf(here(figures,"figure2f.pdf"))
VlnPlot(rnaAggr, ident = "PT", features = c("PCK1","ALDOB","FBP1","G6PC"), ncol = 2, group.by = "diabetes", pt.size = 0)
dev.off()








