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
features <- c("SLC34A1","LRP2","SLC5A1","SLC5A2","HAVCR1","PROM1",
              "VCAM1","CD36","CFH","SLC12A1","CLDN16","SLC12A3",
              "AQP2","SLC26A7","SLC26A4","NPHS2","FLT1",
              "PIEZO2","COL1A2","PTPRC","CD3E","MS4A1","CSF1R")
p1 <- DotPlot(atacAggr, assay="RNA", features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("") +
  ggtitle("A) snATAC-seq Gene Activity")
p2 <- DotPlot(atacAggr, assay="IMPRNA", features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("") +
  ggtitle("B) snATAC-seq Imputed RNA Expression")
pdf(here(figures,"sfigure3.pdf"), width=10, height=10)
grid.arrange(p1,p2)
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
pdf(here(figures,"sfigure8.pdf"), width=10, height=6)
print(p1)
dev.off()
#################################################################
# supplemental figure PT_VCAM1 proportion by donor in snRNA-seq and snATAC-seq
rnaAggr <- readRDS(here("project","analysis","dkd","rna_aggr_prep","step2_anno.rds"))
df <- rnaAggr@meta.data %>%
  dplyr::select(orig.ident, celltype) %>%
  group_by(orig.ident) %>%
  dplyr::mutate(total_cells = n()) %>%
  group_by(orig.ident, celltype) %>%
  dplyr::mutate(total_celltype = n()) %>%
  dplyr::mutate(prop_celltype = total_celltype / total_cells) %>%
  distinct(orig.ident, celltype, prop_celltype) %>%
  dplyr::filter(celltype == "PTVCAM1") %>%
  arrange(orig.ident)

p1 <- df %>%
  ggplot(aes(x=orig.ident, y=prop_celltype, fill=orig.ident)) +
  geom_bar(stat = "identity") +
  xlab("") +
  ylab("Proportion PT_VCAM1") +
  ggtitle("A) snRNA-seq") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")

atacAggr <- readRDS(here("project","analysis","dkd","atac_aggr_prep","step5_chromVAR.rds"))
df <- atacAggr@meta.data %>%
  dplyr::select(orig.ident, celltype) %>%
  group_by(orig.ident) %>%
  dplyr::mutate(total_cells = n()) %>%
  group_by(orig.ident, celltype) %>%
  dplyr::mutate(total_celltype = n()) %>%
  dplyr::mutate(prop_celltype = total_celltype / total_cells) %>%
  distinct(orig.ident, celltype, prop_celltype) %>%
  dplyr::filter(celltype == "PT_VCAM1") %>%
  arrange(orig.ident)

p2 <- df %>%
  ggplot(aes(x=orig.ident, y=prop_celltype, fill=orig.ident)) +
  geom_bar(stat = "identity") +
  xlab("") +
  ylab("Proportion PT_VCAM1") +
  ggtitle("B) snATAC-seq") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")

library(gridExtra)
pdf(here(figures,"sfigure4.pdf"), width=6, height=10)
grid.arrange(p1,p2)
dev.off()
###########################################
# draw ATP1B1 gene model peak coverage
# edit afterwards in inskcape / adobe illustrator etc.
library(Signac)
library(Seurat)
library(here)
library(dplyr)
library(tibble)
library(openxlsx)
library(plyranges)
library(data.table)

figures <- here("project","analysis","dkd","figures")
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan.rds"))
DefaultAssay(atacAggr) <- "peaks"
Idents(atacAggr) <- paste0(atacAggr@meta.data$celltype,"_",atacAggr@meta.data$diabetes)

# cell-specific DAR for dkd vs.control
file <- here("project","analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx")
idents <- getSheetNames(file)
dar.df <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
    dplyr::filter(abs_log2FC > 0.1) %>%
    rownames_to_column(var = "peak")
  if(nrow(df) > 0) { df$celltype <- ident }
    return(df)
}) %>% bind_rows()

dar.gr <- StringToGRanges(dar.df$peak)

# flank the gene coordinates with window
plot.gr <- StringToGRanges("chr1-169106683-169135009")
window <- 20000
start(plot.gr) <- start(plot.gr) - window
end(plot.gr) <- end(plot.gr) + window

file <- here("project","cut_and_run","hTERT_RPTEC","GR","GR_NT","hTERT_GR_consensus.bed")
GRpeak.gr <- fread(file) %>%
  dplyr::rename(chrom = V1, start = V2, end = V3) %>%
  makeGRangesFromDataFrame()
GRpeak.gr <- join_overlap_intersect(GRpeak.gr, plot.gr)

# read in DMR
file <- here("project","analysis","dkd","methylation","ST19_intersection_DMR_with_DAR_and_GR_cut_and_run.xlsx")
dmr.gr <- read.xlsx(file, sheet = "ALL_DMR") %>%
  distinct(seqnames, start, end) %>%
  makeGRangesFromDataFrame()

dmr.gr <- join_overlap_intersect(dmr.gr, plot.gr)

# make dmr a little wider to visualize in coverageplot
dmr_flank.gr <- Extend(dmr.gr, upstream=100, downstream=100)
dmr.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=dmr_flank.gr, color="blue")

# intersect plotting region with dar
select.gr <- join_overlap_intersect(dar.gr, plot.gr)

# make plots
cp <- CoveragePlot(atacAggr, ident=c("PCT_0","PCT_1"), region = plot.gr, peaks=TRUE, links=FALSE)
peaks.plot <- PeakPlot(atacAggr, region = plot.gr)
dar.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=select.gr)
gr.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=GRpeak.gr)
ccan.plot <- LinkPlot(atacAggr, region= plot.gr, min.cutoff=0.4)

# combine tracks
plot <- CombineTracks(list(cp, peaks.plot, dar.plot, gr.plot, dmr.plot, ccan.plot))

# this will print coverage plot with the links
pdf(here(figures, "sfigure_PT_ATP1B1.pdf"))
print(plot)
dev.off()

cp <- CoveragePlot(atacAggr, ident=c("TAL1_0","TAL1_1"), region = plot.gr, peaks=TRUE, links=FALSE)
ccan.plot <- LinkPlot(atacAggr, region = plot.gr, min.cutoff=0.4)
dar.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=select.gr)
gr.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=GRpeak.gr)
dmr.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=dmr_flank.gr, color="blue")
plot <- CombineTracks(list(cp, dar.plot, gr.plot, dmr.plot, ccan.plot))

# this will print coverage plot with the links
pdf(here(figures, "sfigure_TAL1_ATP1B1.pdf"))
print(plot)
dev.off()

cp <- CoveragePlot(atacAggr, ident=c("TAL2_0","TAL2_1"), region = plot.gr, peaks=TRUE, links=FALSE)
ccan.plot <- LinkPlot(atacAggr, region = plot.gr, min.cutoff=0.4)
dar.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=select.gr)
plot <- CombineTracks(list(cp,dar.plot, gr.plot, ccan.plot))

# this will print coverage plot with the links
pdf(here(figures, "sfigure_TAL2_ATP1B1.pdf"))
print(plot)
dev.off()
#############################################
# pairwise comparison of INSR DAR between control and DKD snATAC-seq donors
figures <- here("project","analysis","dkd","figures")
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan.rds"))
DefaultAssay(atacAggr) <- "peaks"

# create lists of control and dkd donors
library_ids <- unique(atacAggr@meta.data$orig.ident)
controls <- library_ids[grepl("Control", library_ids)] %>% sort()
dkd <- library_ids[grepl("DN", library_ids)] %>% sort()

# create a df of unique pairwise combos
pairwise <- expand.grid(controls, dkd)

# extract PCT barcodes from metadata
pct <- atacAggr@meta.data %>%
  dplyr::filter(celltype == "PCT") %>%
  dplyr::select(barcode, orig.ident)

# INSR DAR region
region <- "chr19-7196798-7198626"

# compare each combo for the INSR DAR
res <- lapply(1:nrow(pairwise), function(row) {
  control <- pairwise[row, 1]
  control_cells <- pct[pct$orig.ident == control,]$barcode
  
  dkd <- pairwise[row, 2]
  dkd_cells <- pct[pct$orig.ident == dkd,]$barcode
  
  mark <- FindMarkers(atacAggr, ident.2 = control_cells, ident.1 = dkd_cells, logfc.threshold = 0, features = region)
  return(cbind(mark, control = control, dkd = dkd))
  }) %>% bind_rows()

# create a geom tile plot indicating which pairwise comparisons were significant (red text) and which direction the change was in (blue fill decreased)
p1 <- res %>%
  dplyr::mutate(fill = ifelse(avg_log2FC < 0, 1, 0)) %>%
  dplyr::mutate(textcolor = ifelse(p_val_adj < 0.05, "red", "black")) %>%
  dplyr::mutate(label = round(avg_log2FC, 2)) %>%
  ggplot(aes(x=control, y=dkd, fill=fill, label=label)) + 
    geom_tile() +
    geom_text(aes(label=label, color=textcolor)) +
    scale_color_manual(values = c('black','red')) +
    scale_fill_gradient(low = "white",high = "lightblue") +
    xlab("") +
    ylab("") +
    ggtitle(paste0("A) Pairwise comparison for INSR region: ", region),
            subtitle = "Blue Fill = Decreased accessibility, Red Text = padj < 0.05") +
    NoLegend()

pt <- subset(rnaAggr, celltype %in% c("PT","PTVCAM1"))
exp <- AverageExpression(pt, features = "INSR", group.by = "orig.ident")$SCT %>%
  as.data.frame() 
exp <- exp %>%
  pivot_longer(cols = colnames(exp))
p2 <- exp %>%
  ggplot(aes(x = name, y = value, color = name, fill = name)) + 
  geom_bar(stat = "identity") +
  NoLegend() +
  ggtitle("B) Average INSR expression in proximal tubule") +
  xlab("") +
  ylab("Average INSR expression") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf(here(figures, "sfigure_INSR.pdf"), width = 6, height = 10)
grid.arrange(p1, p2)
dev.off()

#############################################
# ALDOB
file <- here("project","analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx")
idents <- getSheetNames(file)
dar.df <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
    rownames_to_column(var = "peak")
  if(nrow(df) > 0) { df$celltype <- ident }
    return(df)
}) %>% bind_rows()
dar.gr <- StringToGRanges(dar.df$peak)

# flank the gene coordinates with window
plot.gr <- StringToGRanges("chr9-101421439-101449664")
start(plot.gr) <- start(plot.gr) - 20000
end(plot.gr) <- end(plot.gr) + 100000

# intersect plotting region with dar
select.gr <- join_overlap_intersect(dar.gr, plot.gr)

gr.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=GRpeak.gr)

# read in DMR
file <- here("project","analysis","dkd","methylation","ST19_intersection_DMR_with_DAR_and_GR_cut_and_run.xlsx")
dmr.gr <- read.xlsx(file, sheet = "ALL_DMR") %>%
  distinct(seqnames, start, end) %>%
  makeGRangesFromDataFrame()

dmr.gr <- join_overlap_intersect(dmr.gr, plot.gr)

# make dmr a little wider to visualize in coverageplot
dmr_flank.gr <- Extend(dmr.gr, upstream=100, downstream=100)
dmr.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=dmr_flank.gr, color="blue")

# gene plot
cp <- CoveragePlot(atacAggr, ident=c("PCT_0","PCT_1"), region = plot.gr, peaks=TRUE, links=FALSE)
ccan.plot <- LinkPlot(atacAggr, region = plot.gr, min.cutoff=0.7)
dar.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=select.gr)
plot <- CombineTracks(list(cp,dar.plot, gr.plot, dmr.plot, ccan.plot))

# this will print coverage plot with the links
pdf(here(figures, "sfigure_ALDOB.pdf"))
print(plot)
dev.off()

#############################################
# G6PC
file <- here("project","analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx")
idents <- getSheetNames(file)
dar.df <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
    rownames_to_column(var = "peak")
  if(nrow(df) > 0) { df$celltype <- ident }
    return(df)
}) %>% bind_rows()
dar.gr <- StringToGRanges(dar.df$peak)

# flank the gene coordinates with window
plot.gr <- StringToGRanges("chr17-42900799-42914438")
start(plot.gr) <- start(plot.gr) - 100000
end(plot.gr) <- end(plot.gr) + 100000

# intersect plotting region with dar
select.gr <- join_overlap_intersect(dar.gr, plot.gr)

gr.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=GRpeak.gr)

# read in DMR
file <- here("project","analysis","dkd","methylation","ST19_intersection_DMR_with_DAR_and_GR_cut_and_run.xlsx")
dmr.gr <- read.xlsx(file, sheet = "ALL_DMR") %>%
  distinct(seqnames, start, end) %>%
  makeGRangesFromDataFrame()

dmr.gr <- join_overlap_intersect(dmr.gr, plot.gr)

# make dmr a little wider to visualize in coverageplot
dmr_flank.gr <- Extend(dmr.gr, upstream=100, downstream=100)
dmr.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=dmr_flank.gr, color="blue")

# gene plot
cp <- CoveragePlot(atacAggr, ident=c("PCT_0","PCT_1"), region = plot.gr, peaks=TRUE, links=FALSE)
ccan.plot <- LinkPlot(atacAggr, region = plot.gr, min.cutoff=0.4)
dar.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=select.gr)
plot <- CombineTracks(list(cp,dar.plot, gr.plot, dmr.plot, ccan.plot))

# this will print coverage plot with the links
pdf(here(figures, "sfigure_G6PC.pdf"))
print(plot)
dev.off()

#############################################
# FBP1
file <- here("project","analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx")
idents <- getSheetNames(file)
dar.df <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
    rownames_to_column(var = "peak")
  if(nrow(df) > 0) { df$celltype <- ident }
    return(df)
}) %>% bind_rows()
dar.gr <- StringToGRanges(dar.df$peak)

# flank the gene coordinates with window
plot.gr <- StringToGRanges("chr9-94603141-94640249")
start(plot.gr) <- start(plot.gr) - 50000
end(plot.gr) <- end(plot.gr) + 50000

# intersect plotting region with dar
select.gr <- join_overlap_intersect(dar.gr, plot.gr)

# intersect plotting region with GR
gr.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=GRpeak.gr)

# read in DMR
file <- here("project","analysis","dkd","methylation","ST19_intersection_DMR_with_DAR_and_GR_cut_and_run.xlsx")
dmr.gr <- read.xlsx(file, sheet = "ALL_DMR") %>%
  distinct(seqnames, start, end) %>%
  makeGRangesFromDataFrame()

dmr.gr <- join_overlap_intersect(dmr.gr, plot.gr)

# gene plot
cp <- CoveragePlot(atacAggr, ident=c("PCT_0","PCT_1"), region = plot.gr, peaks=TRUE, links=FALSE)
ccan.plot <- LinkPlot(atacAggr, region = plot.gr, min.cutoff=0.4)
dar.plot <- PeakPlot(atacAggr, region = plot.gr, peaks=select.gr)
plot <- CombineTracks(list(cp,dar.plot,gr.plot, dmr.plot, ccan.plot))

# this will print coverage plot with the links
pdf(here(figures, "sfigure_FBP1.pdf"))
print(plot)
dev.off()

###################
# draw atacqc plots
p1 <- VlnPlot(atacAggr, features = "reads_count", group.by="orig.ident", pt.size = 0) +
  ggtitle("A) snATAC-seq total reads per cell") +
  xlab("")
p2 <- VlnPlot(atacAggr, features = "peak_region_fragments", group.by = "orig.ident", pt.size = 0) + 
  ggtitle("B) snATAC-seq peak region fragments per cell") +
  xlab("")
# p3 see snATAC_prep step1_amulet.pdf
p4 <- DimPlot(atacAggr, split.by = "orig.ident", ncol = 6) + 
  ggtitle("D) snATAC-seq Sample Integration")
pdf(here(figures, "sfig12.pdf"), width=15, height=10)
grid.arrange(p1,p2,p4, ncol = 2)
dev.off()

############
# draw rnaqc plots
p1 <- VlnPlot(rnaAggr, features = "nCount_RNA", group.by="orig.ident", pt.size = 0) +
  ggtitle("A) snRNA-seq number RNA counts per cell") +
  xlab("")
p2 <- VlnPlot(rnaAggr, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0) + 
  ggtitle("B) snRNA-seq number features per cell") +
  xlab("")
# p3 see snRNA_prep step1_doublets.pdf
p4 <- DimPlot(rnaAggr, split.by = "orig.ident", ncol = 6) + 
  ggtitle("D) snRNA-seq Sample Integration")
# p5, p6 see snRNA_prep step3 soupx_clusters.pdf
pdf(here(figures, "sfig13.pdf"), width = 20, height = 10)
grid.arrange(p1,p2,p4, ncol = 2)
dev.off()

