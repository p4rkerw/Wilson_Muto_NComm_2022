#!/usr/bin/env Rscript
# this script will eliminate doublets from an aggregated snRNA object prior to preprocessing
#
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
# to run interactively on the RIS compute1 cluster:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

# to run detached:
# git clone https://github.com/p4rkerw/dkd $SCRATCH1/dkd
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log_step1.out Rscript $SCRATCH1/dkd/snRNA_prep/step1_prep.R

library(Seurat) # 4.0.3
library(ggplot2) # 3.3.5
library(harmony) # 0.1.0
library(DoubletFinder) # 2.0.3
library(dplyr) # 1.0.7
library(stringr) # 1.4.0
library(tibble) # 3.1.3
library(here) # 1.0.1
library(future) # 1.21.0
set.seed(1234)

# create output and plots directory
dir.create(here("project","analysis","dkd","rna_aggr_prep","plots"), recursive=TRUE, showWarnings=FALSE)
rna_aggr_prep <- here("project","analysis","dkd","rna_aggr_prep")

# define cellranger-atac aggregation input dir
aggr_input_dir <- here("project","analysis","dkd","cellranger_rna_aggr","outs")

# load aggregated snRNAseq data from cellranger aggregate matrix and create seurat objects
counts <- Read10X_h5(here(aggr_input_dir,"filtered_feature_bc_matrix.h5"))
aggcsv <- read.csv(here(aggr_input_dir,"aggregation.csv"))

# load aggregated snRNAseq data from cellranger aggregate matrix and create seurat objects
rnaAggr <- CreateSeuratObject(counts = counts, min.cells = 10, min.features = 500, 
                              project = "RNA")
# extract GEM groups from individual barcodes using string split and the suffix integer
# use the GEM groups to assign sample origin (control vs. diabetes) from the aggregation.csv metadata
gemgroup <- sapply(strsplit(rownames(rnaAggr@meta.data), split="-"), "[[", 2) %>% as.numeric()
current.gemgroups <- seq(1, length(unique(aggcsv$library_id)))
orig.ident <- factor(aggcsv$library_id, levels=aggcsv$library_id) 
sampleID <- plyr::mapvalues(gemgroup, from = current.gemgroups, to = as.character(aggcsv$library_id))
rnaAggr <- AddMetaData(object=rnaAggr, metadata=data.frame(orig.ident=sampleID, row.names=rownames(rnaAggr@meta.data)))
rnaAggr <- PercentageFeatureSet(rnaAggr, pattern = "^MT-", col.name = "percent.mt")
rnaAggr <- PercentageFeatureSet(rnaAggr, pattern = "^RPL", col.name = "percent.rpl")
rnaAggr <- PercentageFeatureSet(rnaAggr, pattern = "^RPS", col.name = "percent.rps")

# control vs diabetes
rnaAggr@meta.data$diabetes <- plyr::mapvalues(gemgroup, from = current.gemgroups, to=as.character(aggcsv$diabetes))

# draw qc plots
pdf(here(rna_aggr_prep,"plots","step1_qc.pdf"))
  VlnPlot(object = rnaAggr, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0)
  VlnPlot(object = rnaAggr, features = c("percent.mt","percent.rps", "percent.rpl"), ncol = 2, pt.size=0, y.max=1)
  VlnPlot(object = rnaAggr, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "orig.ident",pt.size=0)
  VlnPlot(object = rnaAggr, features = c("percent.mt","percent.rps", "percent.rpl"), ncol = 2, group.by = "orig.ident", pt.size=0, y.max=1)
dev.off()

# filter the aggregated dataset for low quality cells
rnaAggr <- subset(rnaAggr, subset = nFeature_RNA > 500 # use same filter params as aggr
                  & nFeature_RNA < 5000 
                  & nCount_RNA < 16000 
                  & percent.mt < 0.5 
                  & percent.rps < 0.3 
                  & percent.rpl < 0.3)

# Doublet removal with the assumption that doublets represent 6% of cells.
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------

FindDoublets <- function(library_id, seurat_aggregate) {
  rnaAggr <- seurat_aggregate
  seurat_obj <- subset(rnaAggr, idents = library_id)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- RunPCA(seurat_obj)
  # ElbowPlot(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
  # DimPlot(seurat_obj)

  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list_kidney <- paramSweep_v3(seurat_obj, PCs = 1:20, sct = F, num.cores=10)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  pK <- bcmvn_kidney %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  seurat_doublets <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = pK,
                               nExp = round(0.05*length(seurat_obj@active.ident)), 
                               reuse.pANN = FALSE, sct = F)
 
  # create doublet groupings and visualize results
  DF.class <- names(seurat_doublets@meta.data) %>% str_subset("DF.classifications")
  pANN <- names(seurat_doublets@meta.data) %>% str_subset("pANN")
  
  p1 <- ggplot(bcmvn_kidney, aes(x=pK, y=BCmetric)) +
    geom_bar(stat = "identity") + 
    ggtitle(paste0("pKmax=",pK)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2 <- DimPlot(seurat_doublets, group.by = DF.class)
  p3 <- FeaturePlot(seurat_doublets, features = pANN)
  
  outFile <- here(rna_aggr_prep,"plots", paste0("step1_doublets.",library_id,".pdf"))
  pdf(outFile)
    print(p1) # need to use print() when drawing pdf in a function call
    print(p2)
    print(p3)
  dev.off()
  
  # create a df of barcodes and doublet designations
  df_doublet_barcodes <- as.data.frame(cbind(rownames(seurat_doublets@meta.data), seurat_doublets@meta.data[[DF.class]]))
  return(df_doublet_barcodes)
}

# take an aggregated snRNA seurat object and a list of library_id to find doublets. return a df of doublet barcodes
# send DimPlot and FeaturePlot of doublets for each library to plots dir
Idents(rnaAggr) <- "orig.ident"
list.doublet.bc <- lapply(orig.ident, function(x) {FindDoublets(x, seurat_aggregate = rnaAggr)})
doublet_id <- list.doublet.bc %>%
  bind_rows() %>%
  dplyr::rename("doublet_id" = "V2") %>%
  tibble::column_to_rownames(var = "V1") # this is the barcode column
table(doublet_id) # quantify total doublet vs. singlet calls (expect ~6% doublets)
  
# add doublet calls to aggregated snRNA object as doublet_id in meta.data slot
rnaAggr <- AddMetaData(rnaAggr, doublet_id)

# enable parallel processing via future package
plan("multiprocess", workers = 10) 
options(future.globals.maxSize = 100000 * 1024^2) # for 100 Gb RAM
plan()

# visualize doublets before proceeding with snRNA preprocessing
rnaAggr <- SCTransform(rnaAggr, vars.to.regress = c("percent.mt","percent.rpl", "percent.rps","nCount_RNA"), verbose = TRUE)
rnaAggr <- RunPCA(rnaAggr, verbose = TRUE)
# ElbowPlot(rnaAggr, ndims = 50) # to determine number of dimensions for clustering
rnaAggr <- RunHarmony(rnaAggr, "orig.ident", plot_convergence = TRUE, assay.use="SCT")
rnaAggr <- FindNeighbors(rnaAggr, dims = 1:30, verbose = TRUE, reduction = "harmony")
rnaAggr <- FindClusters(rnaAggr, verbose = TRUE, resolution = 0.6)
rnaAggr <- RunUMAP(rnaAggr, dims = 1:30, verbose = TRUE, reduction = "harmony")


rnaAggr@meta.data$doublet_viz <- ifelse(rnaAggr@meta.data$doublet_id == "Singlet",0,1)

Idents(rnaAggr) <- "seurat_clusters"
p1 <- DimPlot(rnaAggr, reduction = "umap", label = TRUE) +
  ggtitle("snRNA Seurat Clustering with Harmony Including Doublets")
p2 <- FeaturePlot(rnaAggr, features = "doublet_viz", order=TRUE) + ggtitle("snRNA-seq DoubletFinder")
pdf(here(rna_aggr_prep,"plots","step1_doublets.pdf"))
print(list(p1,p2))
dev.off()

# save an aggregated snRNA object with doublets
# saveRDS(rnaAggr, here(rna_aggr_prep,"step1_doublets.rds"))

# subset the object for singlets to do additional processing
Idents(rnaAggr) <- "doublet_id"
rnaAggr <- subset(rnaAggr, idents = "Singlet")

# Regress out the mitochondrial reads and nCount_RNA
rnaAggr <- SCTransform(rnaAggr, vars.to.regress = c("percent.mt","percent.rpl", "percent.rps","nCount_RNA"), verbose = TRUE)
rnaAggr <- RunPCA(rnaAggr, verbose = TRUE)
# ElbowPlot(rnaAggr, ndims = 50) # to determine number of dimensions for clustering
rnaAggr <- RunHarmony(rnaAggr, "orig.ident", plot_convergence = TRUE, assay.use="SCT")
rnaAggr <- FindNeighbors(rnaAggr, dims = 1:24, verbose = TRUE, reduction = "harmony")
rnaAggr <- FindClusters(rnaAggr, verbose = TRUE, resolution = 0.8, future.seed=TRUE)
rnaAggr <- RunUMAP(rnaAggr, dims = 1:24, verbose = TRUE, reduction = "harmony")

# visualize the clustering
Idents(rnaAggr) <- "seurat_clusters"
p2 <- DimPlot(rnaAggr, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("snRNA Seurat Clustering with Harmony No Doublets") +
  NoLegend()

marker.genes <- c("CUBN","HAVCR1","SLC5A1","SLC5A2","VCAM1", # PT and PTVCAM1+ markers
                  "CFH", # PEC
                  "SLC12A1", # TAL NKCC2
                  "CLDN10", #MTAL
                  "CLDN16", #CTAL
                  "S100A2", #ATL
                  "SLC12A3","TRPM6", # DCT1 and DCT2 NCC
                  "SCNN1G","TRPV5", # DCT2/CNT ENaC
                  "CALB1", # CNT
                  "AQP2", # PC
                  "ATP6V0D2", # ICA and ICB
                  "SLC4A1","SLC26A7", # ICA
                  "SLC26A4", # ICB
                  "NPHS1","NPHS2", # PODO
                  "PECAM1","FLT1", # ENDO 
                  "IGFBP5","IGFBP7", # PTC and AVR
                  "PLVAP", # PTC and AVR https://www.nature.com/articles/s41467-019-12872-5
                  "EHD3", # GEC
                  "SLC6A6","SLC14A1","AQP1", # EA and DVR
                  "NOS1", # MD
                  "ITGA8","PDGFRB","MEIS2","PIEZO2","REN", # MES and JGA
                  "ACTA2","CALD1", # FIB
                  "PROX1","FLT4","PDPN", # Lymphatics
                  "PTPRC","CD3E","MS4A1", # Lymphocytes
                  "FCGR3A","CD14","CSF1R") # Monocyte / Macrophage

print("Drawing UMAP clusters")
pdf(here(rna_aggr_prep,"plots","step1_clusters.pdf"))
  tryCatch(print(p2), error=function(e) NULL)
  tryCatch(print(p1), error=function(e) NULL)
dev.off()

print("Drawing UMAP markers")
pdf(here(rna_aggr_prep,"plots","step1_markers.pdf")) 
    lapply(marker.genes, function(gene) {
      tryCatch({plot <- FeaturePlot(rnaAggr, features=gene, reduction="umap")
                        AugmentPlot(plot) # downsample
                }, warning=function(w) return(NULL)) 
    })
dev.off()

# save the preprocessed rna file
saveRDS(rnaAggr, file = here(rna_aggr_prep,"step1_prep.rds"))


        
        
        
        
        
        
        
        
        
        
        
        
