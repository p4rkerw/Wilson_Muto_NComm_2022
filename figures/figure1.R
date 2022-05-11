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
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RColorBrewer)

figures <- here("project","analysis","dkd","figures")
dir.create(here(figures), showWarnings=FALSE)

#######################
# figure snATAC-seq UMAP
# read in aggregated snATAC object
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan.rds"))
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
p2 <- DotPlot(atacAggr, features=marker.genes, cols=c("lightgrey","darkgreen")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf(here(figures,"figure1B.pdf"))
print(p2)
dev.off()
####################################################
# figure relation of DAR to nearest TSS
# annotate location of cell-specific DAR for diabetes vs. control
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(dplyr)
library(tibble)
library(RColorBrewer)
library(openxlsx)
library(Signac)
library(stringr)
library(ggplot2)
library(here)

file <- here("project","analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx")
idents <- getSheetNames(file)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
dar.df <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
    rownames_to_column(var = "peak")
  if(nrow(df) > 0) { df$celltype <- ident }
    return(df)
}) %>% bind_rows()

dar.gr <- StringToGRanges(dar.df$peak, sep = c("-","-"))

# annotate the list of GRanges DAR for each cell type
peakAnno <- annotatePeak(dar.gr, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)

p6 <- plotAnnoPie(peakAnno)#total DAR in the dataset
p7 <- plotAnnoBar(peakAnno) #celltype-specific analysis
p8 <- plotDistToTSS(peakAnno)

plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno, title="",ylab="",xlab="")

dist_dar_tss <- peakAnno@anno@elementMetadata@listData$distanceToTSS

# extract peak locations and take first word in string 
peakloc <- peakAnno@anno@elementMetadata@listData$annotation
peakloc <- word(peakloc, 1)

df <- data.frame(peak=dar.df$peak, dist=dist_dar_tss, avg_log2FC=dar.df$avg_log2FC, peakloc=peakloc, abs_log2FC = dar.df$abs_log2FC)
df$peakloc <- recode_factor(df$peakloc, 
  Downstream = "Promoter", 
  Distal = "Distal_Intergenic")

# save a csv of the dar annotations
write.csv(df, file=here("project","analysis","dkd","markers","dar_chipseeker_anno.csv"))

# set factor levels for peak location
df$peakloc <- factor(df$peakloc, levels=c("Promoter","Intron","Exon","Distal_Intergenic","3'","5'"))

# modify any points that are greater than 100kb from TSS to fit in plot
df <- dplyr::mutate(df, distmod = ifelse(abs(dist_dar_tss) > 100000, sign(dist_dar_tss) * 100000, dist_dar_tss))

# create layers to plot DAR meeting threshold with different fill color
toplot1 <- dplyr::filter(df, abs_log2FC > 0.1)
toplot2 <- dplyr::filter(df, abs_log2FC < 0.1)

top_layer <- ggplot() + 
  geom_point(data = toplot1, aes(x=distmod, y=avg_log2FC, color=peakloc)) +
  geom_hline(yintercept=0.1) +
  geom_hline(yintercept=-0.1) +
  theme_bw() +
  xlab("Distance from TSS") +
  ylab("Average log fold change") +
  scale_colour_brewer('DAR location', palette = "Set1")
p1 <- top_layer + geom_point(data = toplot2, size=1, aes(x=distmod, y=avg_log2FC, fill="black"))

# plot the DAR that dont meet threshold on the bottom
p1$layers <- rev(p1$layers)
p1

# in inkscape select the plot area and then Path > Union to flatten geom_point objects
png(here(figures, "figure1b.png"))
p1
dev.off()

pdf(here(figures, "figure1b.pdf"))
p1
dev.off()


####################################################
# figure upset plot cell-specific DAR intersection for dkd vs. control
library(here)
library(openxlsx)
library(ComplexHeatmap)
library(dplyr)
library(tibble)

file <- here("project","analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx")
idents <- getSheetNames(file)

dar.ls <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
    dplyr::filter(abs_log2FC > 0.1) %>%
    rownames_to_column(var = "peak")
  list <- df$peak
  return(list)
})
names(dar.ls) <- idents

# remove non-kidney cell types
id_remove <- names(dar.ls) %in% c("BCELL","TCELL","MONO")
dar.ls <- dar.ls[!id_remove]

levels <- c("PCT","PST","PT_VCAM1","PT_PROM1","PT_CD36","PEC",
  "ATL","TAL1","TAL2","DCT1","DCT2",
  "PC","ICA","ICB","PODO","ENDO",
  "FIB_VSMC_MC")

# reorder the levels
dar.ls <- dar.ls[levels]

mat <- list_to_matrix(dar.ls)

# cell-specific dar shared across cell types
m1 <- make_comb_mat(mat)
# filter by intersection size and set size
# identify combinations with only a single ident by summing all the occurrences of 1 in the string
df <- strsplit(names(comb_size(m1)), split="") %>%
  as.data.frame() 
df <- sapply(df, as.numeric)
colsums <- colSums(df)

# only plot intersections with size greater than 5 unless there's only 1 peak in the group
m1 <- m1[comb_size(m1) > 5 | colsums == 1]

ht <- draw(UpSet(m1, set_order=names(dar.ls)))
od <- column_order(ht)
cs <- comb_size(m1)
p3 <- draw(UpSet(m1, set_order=names(dar.ls)))
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
}) # export to pdf 6x6in figure1c.pdf

# return DAR present in more than one cell type
dar.df <- group_by(dar.df, peak) %>% 
  dplyr::mutate(multi_dar = n_distinct(celltype)) %>% 
  arrange(-multi_dar) %>%
  filter(p_val_adj < 0.05, abs_log2FC > 0.1) %>%
  as.data.frame()
write.csv(here(figures), "multi_dar.csv")

######################################
# plot panther db GO enrichment for DAR in PCT and PST diabetes vs. control
library(data.table)
library(openxlsx)
library(dplyr)

file <- here("project","analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx")
idents <- c("PCT","PST")

dar.df <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
    dplyr::filter(abs_log2FC > 0.1) %>%
    rownames_to_column(var = "peak")
  df$celltype <- ident
  return(df)
}) %>% bind_rows()

genes <- unique(dar.df$gene) # n=514
write.csv(genes, file=here(figures, "ptdargenes.csv"))
# enter gene list into pantherdb at http://geneontology.org/ and export results using filename below

# read in panther results
panther <- fread(here("project","analysis","dkd","panther","panther_allpt_dar.macs2.diab_vs_ctrl.txt"), skip=11, header=TRUE) #skip white space and database annotation
# barplot(panther$`upload_1 (fold Enrichment)`, xlab="GO Biological Process", ylab = "Pathway Fold Enrichment", main="Proximal tubule DAR in diabetes")
panther$order <- seq(nrow(panther))
panther$color <- ifelse(panther$order <=25, "Top25", "Not_Top25")

ggplot(panther, aes(x=order, y=`upload_1 (fold Enrichment)`, fill=color, color=color)) + 
  geom_bar(show.legend=FALSE,stat="identity") +
  ylab("Pathway fold enrichment") +
  xlab("GO Biological Process") +
  theme_bw() # export 6x6 pdf figure1d.pdf

######################################
# draw INSR gene model peak coverage
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
p5 <- CoveragePlot(atacAggr, ident=c("PCT_0","PCT_1"), region = "INSR", peaks=FALSE, links=FALSE)
pdf(here(figures, "figure1e_coverage.pdf"))
print(p5)
dev.off()

# # add multiple tracks together (have to use genomic coordinates rather than gene name here)
p6 <- LinkPlot(atacAggr, region="chr19-7112255-7293931", min.cutoff=0.4)
p7 <- CombineTracks(list(p5, p6))

# this will print coverage plot with the links
pdf(here(figures, "figure1e.pdf"))
print(p7)
dev.off()
##################################
# draw insr pt violin plot
rna_aggr_prep <- here("project","analysis","dkd","rna_aggr_prep")
rnaAggr <- readRDS(here(rna_aggr_prep, "step2_anno.rds"))
pdf(here(figures,"figure1f.pdf"))
VlnPlot(rnaAggr, ident = "PT", features = "INSR", group.by = "diabetes", pt.size = 0)
dev.off()

