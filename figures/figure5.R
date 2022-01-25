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
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

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

# # read in TF motif enrichment for diabetes vs. control
file <- here("project","analysis","dkd","markers","motifs_in_dar.macs2.celltype.diab_vs_ctrl.xlsx")
atac_idents <- getSheetNames(file)
motifs.df <- lapply(atac_idents, function(ident) {
  df <- read.xlsx(file, sheet=ident)
  df$celltype <- ident
  return(df)
  }) %>% bind_rows

# reorder the cell type levels and filter by padj
motifs.df$celltype <- factor(motifs.df$celltype, levels=atac_idents)
motifs.sub <- dplyr::filter(motifs.df, pvalue < 0.05)

###########################
# read in cell-specific degs
file <- here("project","analysis","dkd","markers","deg.celltype.diab_vs_ctrl.xlsx")
rna_idents <- getSheetNames(file)
deg.df <- lapply(rna_idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    rownames_to_column(var = "gene")
  df$celltype <- ident
  return(df)
}) %>% bind_rows()

# update deg.df so rna celltypes match the atac data
new_atac_idents <- c("PT","DCT1","PT","TAL1","DCT2","TAL2","PTVCAM1","ATL","ENDO","ICA","PC","ICB","PEC","FIB","LEUK","LEUK","PODO","LEUK","PT")
new_atac_celltype <- plyr::mapvalues(motifs.sub$celltype, from = atac_idents, to = new_atac_idents)

# plot deg on y-axis vs. motif enrichment on x-axis for differentially expressed TF that also 
# so a change in motif enrichment
motifs.sub$celltype <- new_atac_celltype
colnames(motifs.sub)[8] <- "gene"
deg_motifs <- deg.df[deg.df$gene %in% motifs.sub$gene,]
deg_motifs <- left_join(deg_motifs, motifs.sub, by=c("gene","celltype"))

# filter out repeat motifs
deg_motifs <- dplyr::arrange(deg_motifs, motif, p_val_adj) 
deg_motifs <- dplyr::distinct(deg_motifs, motif, celltype, .keep_all=TRUE)
p3 <- ggplot(deg_motifs, aes(x=fold.enrichment, y=avg_log2FC, color=celltype)) + geom_point()

# label selected motifs
motifs.label <- c("NR3C1","NR3C2","HIF1A")
deg_motifs <- dplyr::mutate(deg_motifs, label = ifelse(gene %in% motifs.label, 1, 0)) %>%
  dplyr::mutate(alpha = ifelse(label == 1, 1, 0.5)) # make background more transparent
deg_motifs <- dplyr::arrange(deg_motifs, label)
p3 <- ggplot(deg_motifs, aes(x=fold.enrichment, y=avg_log2FC, color=gene, alpha=alpha, label=ifelse(label==1,celltype,""))) +
 geom_text() +
 geom_point() +
 scale_color_manual(values = c("NR3C1" = "red", "NR3C2" = "blue", "HIF1A" = "black")) +
 theme_bw() +
 guides(alpha = "none")
pdf(here(figures, "figure5a.pdf"))
print(p3)
dev.off()

###########################################
# figure 3b,c
# visualize NR3C1/2 chromvar motif activity
p16a <- VlnPlot(atacAggr, assay="chromvar", features = "MA0113.3", group.by="celltype", split.by="diabetes", pt.size=0) +
 ggtitle("NR3C1 chromVAR motif activity")
png(here(figures, "figure5b1.png"))
print(p16a)
dev.off()

p18 <- VlnPlot(atacAggr, assay="chromvar", features = "MA0101.1", group.by="celltype", split.by="diabetes", pt.size=0) +
 ggtitle("REL chromVAR motif activity")
png(here(figures, "figure5b4.png"))
print(p18) 
dev.off()

pdf(here(figures, "figure5b4.pdf"))
print(p18) 
dev.off()

###########################################
# fkbp5 coverageplot
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan.rds"))
Idents(atacAggr) <- "celltype"
DefaultAssay(atacAggr) <- "peaks"

# calculate threshold for peaks open in at least x% of cells
threshold <- table(atacAggr@meta.data$celltype == "PCT")[[2]]/20

# identify accessible peaks for celltype
peaks.gr <- AccessiblePeaks(atacAggr, assay="peaks", ident="PCT", min.cells=threshold) %>%
  StringToGRanges()


Idents(atacAggr) <- paste0(atacAggr@meta.data$celltype,"_",atacAggr@meta.data$diabetes)

# read in the GR cut and run peaks and intersect with plot region
library(plyranges)
file <- here("project","cut_and_run","hTERT_RPTEC","GR","GR_NT","hTERT_GR_consensus.bed")
GRpeak.gr <- fread(file) %>%
  dplyr::rename(chrom = V1, start = V2, end = V3) %>%
  makeGRangesFromDataFrame()
GRpeak.gr <- join_overlap_intersect(GRpeak.gr, plot.gr)

# read cell-specific DAR dkd vs. control
file <- here("project","analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx")
dar.pct <- read.xlsx(file, sheet="PCT", rowNames = TRUE) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  rownames_to_column(var = "peak")

dar.gr <- StringToGRanges(dar.pct$peak)
region <- "chr6-35573585-35750000"
plot.gr <- StringToGRanges(region)
dar.gr <- join_overlap_intersect(dar.gr, plot.gr)
peaks_to_plot.gr <- join_overlap_intersect(peaks.gr, plot.gr)

# fkbp5 with 1k upstream and 10k downstream flank
cp <- CoveragePlot(atacAggr, ident=c("PCT_0","PCT_1"), region = region, peaks=TRUE, links=FALSE)

peaks.plot <- PeakPlot(atacAggr, region=region, peaks=peaks_to_plot.gr)
dar.plot <- PeakPlot(atacAggr, region = region, peaks=dar.gr)
# tcre.plot <- PeakPlot(atacAggr, region = region, peaks=tcre.gr)
gr.plot <- PeakPlot(atacAggr, region = region, peaks=GRpeak.gr)
ccan.plot <- LinkPlot(atacAggr, region= region, min.cutoff=0.5)

plot <- CombineTracks(list(cp,peaks.plot, dar.plot, gr.plot, ccan.plot))

pdf(here(figures, "figure5_fkbp5.pdf"))
print(plot)
dev.off()

##################################
# draw fkbp5 pt violin plot
rna_aggr_prep <- here("project","analysis","dkd","rna_aggr_prep")
rnaAggr <- readRDS(here(rna_aggr_prep, "step2_anno.rds"))
pdf(here(figures,"figure5f.pdf"))
VlnPlot(rnaAggr, idents = "PT", features = c("FKBP5"), group.by = "diabetes", pt.size = 0)
VlnPlot(rnaAggr, idents = "PT", features = c("NR3C1"), group.by = "diabetes", pt.size = 0)
VlnPlot(rnaAggr, idents = c("PT","PTVCAM1","TAL1","DCT1","PC","ENDO","PODO","LEUK"), features = "FKBP5", split.by = "diabetes", pt.size = 0)
VlnPlot(rnaAggr, features = "FKBP5", split.by = "diabetes", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  coord_flip()
dev.off()
