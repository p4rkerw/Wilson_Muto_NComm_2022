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
library(plyranges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

figures <- here("project","analysis","dkd","figures")
dir.create(here(figures), showWarnings=FALSE)

##############################################
# tcre histogram distance to TSS and barplot
# read in the tcre
tcre <- fread(here("project","analysis","dkd","scafe","scafe_pool","annotate","pool","bed","pool.CRE.annot.bed.gz")) %>%
  dplyr::mutate(peak = paste0(V1,"-",V2,"-",V3))
tcre.gr <- StringToGRanges(tcre$peak)
               
# annotate
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene             
peakAnno <- annotatePeak(tcre.gr, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)
dist_tcre_tss <- peakAnno@anno@elementMetadata@listData$distanceToTSS

# extract peak locations and take first word in string 
peakloc <- peakAnno@anno@elementMetadata@listData$annotation
peakloc <- word(peakloc, 1)

df <- as.data.frame(peakAnno@anno@ranges)
df$seqnames <- as.data.frame(peakAnno@anno@seqnames)
df$dist <- dist_tcre_tss
df$peakloc <- peakloc
df$peakloc <- recode_factor(df$peakloc, 
  Downstream = "Promoter")

# draw a histogram for dist from TSS
limit <- 2 * sd(df$dist)
p1 <- ggplot(df, aes(x=dist)) + 
  geom_histogram() + 
  xlim(-limit, limit) +
  theme_bw() +
  xlab("Distance to Nearest Protein-Coding TSS") +
  ylab("Number of tCRE") +
  theme(text = element_text(size = 16))

pdf(here(figures, "figure3a.pdf"))
print(p1)
dev.off()

# change factor levels of tcre locations
df$peakloc <- factor(df$peakloc, levels = c("Promoter","Intron","Exon","Distal","3'","5'"))

# draw barplot of tcre locations
toplot <- group_by(df, peakloc) %>% 
  mutate(peak = paste0(seqnames$value,"-",start,"-",end)) %>%
  mutate(count = n_distinct(peak)) %>% 
  as.data.frame() %>% 
  distinct(peakloc, .keep_all = TRUE)
p1 <- ggplot(toplot, aes(x=peakloc, y=count, fill=peakloc)) + 
  geom_bar(stat="identity", show.legend=FALSE) +
  xlab("") +
  ylab("Number of tCRE") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
       text = element_text(size = 16))

pdf(here(figures, "figure3b.pdf"))
print(p1)
dev.off()

############################################
library(ChIPpeakAnno)
# read in the tcre
tcre <- fread(here("project","analysis","dkd","scafe","scafe_pool","annotate","pool","bed","pool.CRE.annot.bed.gz")) %>%
  dplyr::mutate(peak = paste0(V1,"-",V2,"-",V3))
tcre.gr <- StringToGRanges(unique(tcre$peak))

atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atac_peaks.gr <- fread(here(atac_aggr_prep,"step2_peaks.gr")) %>%
  makeGRangesFromDataFrame()

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
dar.peaks <- unique(dar.df$peak)
dar.gr <- StringToGRanges(dar.peaks, sep = c("-","-"))
dar.gr <- sort(dar.gr)

# draw venn diagram
makeVennDiagram(list(atac_peaks.gr, tcre.gr, dar.gr), 
                      NameOfPeaks = c("atac","tcre","dar"), 
                      fill = c("cornflowerblue", "lightgreen", "yellow"),
                      alpha = 0.50)
# number of overlapping peaks does not appear to be correct so use the following nums
tcre_int_dar <- join_overlap_intersect(tcre.gr, dar.gr) # n=494
tcre_int_atac <- join_overlap_intersect(tcre.gr, atac_peaks.gr) # n=23,048

pdf(here(figures, "figure3c.pdf"))
makeVennDiagram(list(atac_peaks.gr, tcre.gr, dar.gr), 
                      NameOfPeaks = c("atac","tcre","dar"), 
                      fill = c("cornflowerblue", "lightgreen", "yellow"),
                      alpha = 0.50)
dev.off()
         
############################################
# vcam1 gene model for PCT
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan.rds"))
Idents(atacAggr) <- "celltype"
DefaultAssay(atacAggr) <- "peaks"

# calculate threshold for peaks open in at least x% of cells
threshold <- table(atacAggr@meta.data$celltype == "PCT")[[2]]/50

# identify accessible peaks for celltype
peaks.gr <- AccessiblePeaks(atacAggr, assay="peaks", ident="PCT", min.cells=threshold) %>%
  StringToGRanges()

# designate plot region around vcam1
region <- "chr1-100642257-100740029"
plot.gr <- StringToGRanges(region)

# read in the GR cut and run peaks and intersect with plot region
file <- here("project","cut_and_run","hTERT_RPTEC","GR","GR_NT","GR_NT_peaks.narrowPeak")
GRpeak.gr <- fread(file) %>%
  dplyr::rename(chrom = V1, start = V2, end = V3) %>%
  makeGRangesFromDataFrame()
GRpeak.gr <- join_overlap_intersect(GRpeak.gr, plot.gr)

# read DAR for PCT vs PT_VCAM1
file <- here("project","analysis","dkd","markers","dar.macs2.PCT_vs_PTVCAM1.markers.xlsx")
dar <- read.xlsx(file, rowNames = TRUE) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  rownames_to_column(var = "peak")
dar.gr <- StringToGRanges(dar$peak)
dar.gr <- join_overlap_intersect(dar.gr, plot.gr)

# read in the tCRE peaks and intersect with plot region
file <- here("project","analysis","dkd","scafe","scafe_pool","annotate","pool","bed","pool.CRE.annot.bed.gz")
tcre.gr <- fread(file) %>%
  dplyr::rename(chrom = V1, start = V2, end = V3) %>%
  makeGRangesFromDataFrame()
tcre.gr <- join_overlap_intersect(tcre.gr, plot.gr)

file <- here("project","analysis","dkd","methylation","ST19_intersection_DMR_with_DAR_and_GR_cut_and_run.xlsx")
dmr.gr <- read.xlsx(file, sheet = "ALL_DMR") %>%
  distinct(seqnames, start, end, dmr) %>%
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)

dmr.gr <- join_overlap_intersect(dmr.gr, plot.gr)

# make dmr a little wider to visualize in coverageplot
dmr_flank.gr <- Extend(dmr.gr, upstream=100, downstream=100)

# intersect peaks with plotting region
peaks_to_plot.gr <- join_overlap_intersect(peaks.gr, plot.gr)

# vcam1 region plot
cp <- CoveragePlot(atacAggr, 
                   ident=c("PCT","PT_VCAM1"), 
                   region = region, 
                   peaks=TRUE, 
                   links=FALSE)

peaks.plot <- PeakPlot(atacAggr, region=region, peaks=peaks_to_plot.gr)
tcre.plot <- PeakPlot(atacAggr, region = region, peaks=tcre.gr)
dar.plot <- PeakPlot(atacAggr, region = region, peaks=dar.gr)
dmr.plot <- PeakPlot(atacAggr, region=region, peaks=dmr_flank.gr, color="blue")
gr.plot <- PeakPlot(atacAggr, region = region, peaks=GRpeak.gr)
ccan.plot <- LinkPlot(atacAggr, region= region, min.cutoff=0.2)
plot <- CombineTracks(list(cp, peaks.plot, tcre.plot, dar.plot, gr.plot, dmr.plot, ccan.plot))

pdf(here(figures, "figure3_vcam1.pdf"))
print(plot)
dev.off()

######################################
# draw vcam1 violin plots
rna_aggr_prep <- here("project","analysis","dkd","rna_aggr_prep")
rnaAggr <- readRDS(here(rna_aggr_prep, "step2_anno.rds"))
pdf(here(figures,"figure3f.pdf"))
VlnPlot(rnaAggr, idents = c("PT","PTVCAM1"), features = c("VCAM1"), pt.size = 0)
dev.off()
