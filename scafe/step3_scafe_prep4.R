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
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log.scafe_prep4.out Rscript $SCRATCH1/dkd/scafe/scafe_prep4.R


library(Seurat) # 4.0.3
library(SeuratObject)
library(ggplot2) # 3.3.5
library(harmony) # 0.1.0
library(dplyr) # 1.0.7
library(stringr) # 1.4.0
library(tibble) # 3.1.3
library(here) # 1.0.1
library(openxlsx)
set.seed(1234)

# create output and plots directory
dir.create(here("project","analysis","dkd","scafe","scafe_pool","plots"), recursive=TRUE, showWarnings=FALSE)
rna_aggr_prep <- here("project","analysis","dkd","scafe","scafe_pool","count")

library_ids=c("Control_4","Control_5","DN_4","DN_5")
obj.ls <- lapply(seq(library_ids), function(index) {
	library_id <- library_ids[index]
	counts <- Read10X(here(rna_aggr_prep, library_id, "matrix"))
	obj <- CreateSeuratObject(counts)
	obj@meta.data$orig.ident <- library_id
	return(obj)
	})

rnaAggr <- merge(obj.ls[[1]], c(obj.ls[[2]], obj.ls[[3]], obj.ls[[4]]))

# fix barcodes
# extract GEM groups from individual barcodes using string split and the suffix integer
# use the GEM groups to assign sample origin (control vs. diabetes) from the aggregation.csv metadata
gemgroup <- sapply(strsplit(rownames(rnaAggr@meta.data), split="_"), "[[", 2) %>% as.numeric()
barcodes <- sapply(strsplit(rownames(rnaAggr@meta.data), split="-"), "[[", 1)
bc_update <- paste0(barcodes,"-", gemgroup)

# control vs diabetes
group <- c(0,0,1,1) # diabetes=1 control=0
orig.ident <- rnaAggr@meta.data$orig.ident
rnaAggr@meta.data$diabetes <- plyr::mapvalues(orig.ident, from = library_ids, to=group)

# draw qc plots
pdf(here("project","analysis","dkd","scafe","scafe_pool","plots","step1_qc4.pdf"))
  VlnPlot(object = rnaAggr, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0)
  VlnPlot(object = rnaAggr, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "orig.ident",pt.size=0)
dev.off()

barcodes <- read.csv(here("project","analysis","dkd","barcodes","rna_barcodes.csv"))
bc_suffix <- lapply(seq(library_ids), function(index) {
  ident <- library_ids[index]
  df <- dplyr::filter(barcodes, orig.ident == ident)
  df <- dplyr::mutate(df, bc_suffix = paste0(barcode,"-1_",index))
  return(df$bc_suffix)
  }) %>% unlist()

# filter the aggregated dataset for annotated barcodes
rnaAggr <- rnaAggr[, rownames(rnaAggr@meta.data) %in% bc_suffix]

# add cell type anno from snRNA
bc_filter <- barcodes[barcodes$orig.ident %in% library_ids, ]
rownames(bc_filter) <- bc_suffix
rnaAggr <- AddMetaData(rnaAggr, bc_filter)

# visualize doublets before proceeding with snRNA preprocessing
rnaAggr <- SCTransform(rnaAggr, vars.to.regress = c("nCount_RNA","nFeature_RNA"), verbose = TRUE)
rnaAggr <- RunPCA(rnaAggr, verbose = TRUE)
# ElbowPlot(rnaAggr, ndims = 50) # to determine number of dimensions for clustering
rnaAggr <- RunHarmony(rnaAggr, "orig.ident", plot_convergence = TRUE, assay.use="SCT")
rnaAggr <- FindNeighbors(rnaAggr, dims = 1:30, verbose = TRUE, reduction = "harmony")
rnaAggr <- FindClusters(rnaAggr, verbose = TRUE, resolution = 0.6)
rnaAggr <- RunUMAP(rnaAggr, dims = 1:30, verbose = TRUE, reduction = "harmony")

Idents(rnaAggr) <- "seurat_clusters"
p1 <- DimPlot(rnaAggr, reduction = "umap", label = TRUE) +
  ggtitle("scafe Seurat Clustering with Harmony")

Idents(rnaAggr) <- "celltype"  
p2 <- DimPlot(rnaAggr, reduction = "umap", label = TRUE) +
  ggtitle("scafe Seurat Clustering with Harmony")
p3 <- DimPlot(rnaAggr, reduction = "umap", label = TRUE, group.by="orig.ident") +
  ggtitle("scafe Seurat Clustering with Harmony")
p4 <- DimPlot(rnaAggr, reduction = "umap", label = TRUE, group.by="diabetes") +
  ggtitle("scafe Seurat Clustering with Harmony")      

print("Drawing UMAP clusters")
pdf(here("project","analysis","dkd","scafe","scafe_pool","plots","step1_clusters.pdf"))
  tryCatch(print(p1), error=function(e) NULL)
  tryCatch(print(p2), error=function(e) NULL)
  tryCatch(print(p3), error=function(e) NULL)
  tryCatch(print(p4), error=function(e) NULL)
dev.off()


# save the preprocessed rna file
saveRDS(rnaAggr, file = here(rna_aggr_prep,"scafe_prep.rds"))

# wrapper functions for FindMarkers 
CelltypeMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding DEG for: ",cluster))
  deg <- FindMarkers(seurat_aggregate, 
                     ident.1 = cluster,    
                     min.pct = 0.2) # find all cluster-specific degs
  return(deg)
}

CompareMarkers <- function(cluster, seurat_aggregate, test, ref, meta_group) {
  print(paste0("Comparing DEG for: ",cluster))
  group <- dplyr::select(seurat_aggregate@meta.data, all_of(meta_group))[,1]
  seurat_aggregate@meta.data$celltype.stim <- paste0(seurat_aggregate@meta.data$celltype,"_", group)
  Idents(seurat_aggregate) <- "celltype.stim"
  deg <- tryCatch(FindMarkers(seurat_aggregate, 
                     ident.1 = paste0(cluster,"_", test),  
                     ident.2 = paste0(cluster, "_", ref),
                     logfc.threshold = 0.25),
                  error=function(e) NULL) # find all celltype specific degs that differ between disease and control
  return(deg)
}

levels(rnaAggr) <- c("PT","PTVCAM1","PEC","ATL","TAL1",
           "TAL2","DCT1","DCT2","PC","ICA",
           "ICB","PODO","ENDO","MES","FIB",
           "LEUK")
rnaAggr@meta.data$celltype <- rnaAggr@active.ident

# FindMarkers and write to an xlsx file with default parameters
markers <- here("project","analysis","dkd","scafe","markers")
dir.create(here(markers), showWarnings = FALSE)
idents <- levels(rnaAggr@meta.data$celltype)

list.disease.deg <- lapply(idents, function(x) {CompareMarkers(x, seurat_aggregate = rnaAggr, test=1, ref=0, meta_group="diabetes")})                          
write.xlsx(list.disease.deg, file = here(markers,"tcre4.celltype.diab_vs_ctrl.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)

list.cluster.deg <- lapply(idents, function(x) {CelltypeMarkers(x, seurat_aggregate = rnaAggr)})
write.xlsx(list.cluster.deg, file = here(markers,"tcre4.celltype.markers.xlsx"), sheetName = idents, rowNames = T, overwrite=TRUE)       

# identify DAR between PCT and PT_VCAM1
tcre <- FindMarkers(rnaAggr, ident.1 = "PTVCAM1", ident.2 = "PT", logfc.threshold = 0) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    rownames_to_column(var = "peak") %>%
    dplyr::mutate(peak = str_sub(peak, 1, nchar(peak) - 2))
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tcre.gr <- StringToGRanges(tcre$peak, sep = c("-","-"))
tcre.gr$avg_log2FC <- tcre$avg_log2FC
tcre.gr$p_val_adj <- tcre$p_val_adj
# set up annotation 
gene.ranges <- genes(EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(gene.ranges),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(gene.ranges) <- ucsc.levels
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')
tcre.gr <- join_overlap_intersect(tcre.gr, gene.ranges)
tcre.df <- as.data.frame(tcre.gr)
write.csv(tcre.df, file = here(markers,"tcre.macs2.PT_vs_PTVCAM1.markers.csv"), row.names = FALSE, quote = FALSE)		  
		  
# annotate location of cell-specific DAR for diabetes vs. control
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
file <- here(markers,"tcre4.celltype.markers.xlsx")
tcre.df <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    rownames_to_column(var = "peak")
  # remove strand orientation
  df <- dplyr::mutate(df, peak_update = str_sub(peak, 1, nchar(peak) - 2))
  df$peak <- df$peak_update
  return(df)
}) %>% bind_rows()

library(Signac)
tcre.gr <- StringToGRanges(tcre.df$peak, sep = c("-","-"))

# annotate the list of GRanges DAR for each cell type
library(ChIPseeker)
peakAnno <- annotatePeak(tcre.gr, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)

dist_tcre_tss <- peakAnno@anno@elementMetadata@listData$distanceToTSS

# extract peak locations and take first word in string 
peakloc <- peakAnno@anno@elementMetadata@listData$annotation
peakloc <- word(peakloc, 1)

df <- data.frame(peak=tcre.df$peak, dist=dist_tcre_tss, avg_log2FC=tcre.df$avg_log2FC, peakloc=peakloc)
df$peakloc <- recode_factor(df$peakloc, 
  Downstream = "Promoter", 
  Distal = "Distal_Intergenic")

# set factor levels for peak location
df$peakloc <- factor(df$peakloc, levels=c("Promoter","Intron","Exon","Distal_Intergenic","3'","5'"))

# modify any points that are greater than 100kb from TSS to fit in plot
df <- dplyr::mutate(df, distmod = ifelse(abs(dist_tcre_tss) > 100000, sign(dist_tcre_tss) * 100000, dist_tcre_tss))

# in inkscape select the plot area and then Path > Union to flatten geom_point objects
library(RColorBrewer)
png(here("project","analysis","dkd","scafe","scafe_pool","plots","tcre.png"))
ggplot(df, aes(x=distmod, y=avg_log2FC, color=peakloc)) + 
  geom_point() +
  geom_hline(yintercept=0.1) +
  geom_hline(yintercept=-0.1) +
  theme_bw() +
  xlab("Distance from TSS") +
  ylab("Average log fold change") +
  scale_colour_brewer('tCRE location', palette = "Set1")
dev.off()
