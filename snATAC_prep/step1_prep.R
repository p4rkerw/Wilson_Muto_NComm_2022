#!/usr/bin/env Rscript
# this script will preprocess aggregated snATACseq data from kidney cortex samples
# counted and aggregated by cellranger-atac v2.0 without library normalization
#
# to run locally
# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v /mnt/g/diabneph:$HOME/project \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.0

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
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log.atac.step1.out Rscript $SCRATCH1/Wilson_Muto_NComm_2022/snATAC_prep/step1_prep.R

library(Signac) # 1.3.0
library(Seurat) # 4.0.3
library(GenomeInfoDb) # 1.28.1
library(harmony) # 0.1.0
library(EnsDb.Hsapiens.v86)
library(ggplot2) # 3.3.5
library(patchwork) # 1.1.1
library(tibble) # 3.1.3
library(dplyr) # 1.0.7
library(here) # 1.0.1
library(data.table) # 1.14.0

start <- Sys.time()

# create output and plots directory
dir.create(here("project","analysis","dkd","atac_aggr_prep","plots"), recursive=TRUE, showWarnings=FALSE)
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")

# define cellranger-atac aggregation input dir
aggr_input_dir <- here("project","analysis","dkd","cellranger_atac_aggr_13","outs")

# define cellranger-atac count dir
# individual counts are subfolders labeled with corresponding library_id
count_input_dir <- here("project","cellranger_atac_counts","version_2.0")

# load aggregated snATACseq data and create a seurat object
counts <- Read10X_h5(here(aggr_input_dir,"filtered_peak_bc_matrix.h5"))
metadata <- read.csv(here(aggr_input_dir ,"singlecell.csv"), header = TRUE, row.names = 1)
aggcsv <- read.csv(here(aggr_input_dir,"aggregation_csv.csv"), header = TRUE, row.names = 1)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = here(aggr_input_dir,"fragments.tsv.gz"),
  min.cells = 10,
  min.features = 200
)

atacAggr <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# Add sample information to the metadata of the Seurat object
gemgroup <- sapply(strsplit(rownames(atacAggr@meta.data), split="-"), "[[", 2) 
current.gemgroups <- seq(length(rownames(aggcsv))) # no. gemgroups is no. samples
orig.ident <- rownames(aggcsv)
sampleID <- plyr::mapvalues(gemgroup, from = current.gemgroups, to = orig.ident)
atacAggr <- AddMetaData(object=atacAggr, metadata=data.frame(orig.ident=sampleID, 
                                                             row.names=rownames(atacAggr@meta.data)))

# update metadata columns for disease groups
gemgroup <- sapply(strsplit(rownames(atacAggr@meta.data), split="-"), "[[", 2) 
current.gemgroups <- seq(n_distinct(gemgroup)) # no. gemgroups is no. samples
diabetes <- c(0,0,0,0,0,1,1,1,1,1,0,1,1) # 0=control, 1=diabetes

meta_group <- plyr::mapvalues(gemgroup, from = current.gemgroups, to = diabetes)
atacAggr   <- AddMetaData(object=atacAggr , metadata=data.frame(diabetes=meta_group, 
                                                             row.names=rownames(atacAggr@meta.data)))
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(atacAggr) <- annotations

# calculate nucleosome signal
atacAggr <- NucleosomeSignal(object = atacAggr)

# Add the metadata for mitochondrial fragments from individual snATACseq counts into the aggregated dataset 
# Collect total fragment number data from each of the original CellRangerATAC datasets.
metaqc <- 
  mclapply(current.gemgroups, 
         function(gemgroup) {
           sampleName <- rownames(aggcsv)[gemgroup]
           file <- here(count_input_dir,sampleName,"outs","singlecell.csv")
           df <- read.csv(file, header=TRUE, row.names=1)
           rownames(df) <- paste(substr(rownames(df), 1, 16), gemgroup, sep = "-") # change gemgroup to reflect sample order
           df <- tibble::rownames_to_column(df, var = "barcode")
           return(df)
         }, mc.cores=10) %>%
  bind_rows() %>%
  tibble::column_to_rownames(var = "barcode") %>%
  dplyr::select(c("total","mitochondrial"))
atacAggr <- AddMetaData(atacAggr,metaqc)
remove(metaqc)

# compute TSS enrichment score per cell
# store TSS matrix (fast=FALSE)
atacAggr <- TSSEnrichment(object = atacAggr, fast = FALSE)

# QC metrics and filtering
atacAggr$pct_reads_in_peaks <- atacAggr$peak_region_fragments / atacAggr$passed_filters * 100 # %fragments in peaks
atacAggr$mito_ratio <- atacAggr$mitochondrial / atacAggr$total # %fragments mapping to mitochondrial genome
atacAggr$nucleosome_group <- ifelse(atacAggr$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
atacAggr$high.tss <- ifelse(atacAggr$TSS.enrichment > 2, 'High', 'Low')

# draw QC plots
# Note: pct_read_in_peaks histogram does not look normally distributed and is left-skewed
# this may represent low quality nuclei OR fewer peaks called in less common cell types (eg. podocytes, leukocytes etc)
# examine this metric again after cell-specific macs2 peaks
# similarly, mito_ratio is higher in podocytes than other cell types. it may be helpful to filter by mito_ratio
# within cell types rather than using a single threshold for all cell types
pdf(here(atac_aggr_prep,"plots","step1_qc.pdf"))
VlnPlot(
  object = atacAggr,
  features = c('pct_reads_in_peaks','peak_region_fragments','mito_ratio'), 
  pt.size = 0,
  ncol = 2) + NoLegend()
p1 <- qplot(atacAggr@meta.data$pct_reads_in_peaks, geom="histogram")
p2 <- qplot(atacAggr@meta.data$peak_region_fragments, geom="histogram")
p3 <- qplot(atacAggr@meta.data$mito_ratio, geom="histogram")
print(list(p1,p2,p3))
FragmentHistogram(object = atacAggr, group.by = 'nucleosome_group')
TSSPlot(atacAggr, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
dev.off()

# write the object to file
# saveRDS(atacAggr, file = here(atac_aggr_prep,"step1_prep.rds"))
# atacAggr <- readRDS(here(atac_aggr_prep,"step1_qc.rds"))
atacAggr_sub <- subset(
  x = atacAggr,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 20000 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# perform normalization and dimensional reduction of filtered snATAC object
atacAggr_sub <- RunTFIDF(atacAggr_sub)
atacAggr_sub <- FindTopFeatures(atacAggr_sub, min.cutoff = 'q0')
atacAggr_sub <- RunSVD(atacAggr_sub)

# show that first dimension is related to depth of sequencing
pdf(here(atac_aggr_prep, "plots","step1_depthCor.pdf"))
DepthCor(atacAggr_sub)
dev.off()

# ElbowPlot(atacAggr_sub, ndim = 40) # select number of dimensions for UMAP embedding
atacAggr_sub <- RunHarmony(atacAggr_sub, group.by.vars = "orig.ident", reduction = "lsi", assay.use = "peaks", project.dim = FALSE)
atacAggr_sub <- FindNeighbors(object = atacAggr_sub, reduction = "harmony", dims = 2:30)
atacAggr_sub <- FindClusters(object = atacAggr_sub, verbose = TRUE, algorithm = 1) # Louvain algorithm
atacAggr_sub <- RunUMAP(object = atacAggr_sub, reduction = "harmony", dims = 2:30)

# plot the original clustering prior to label transfer / doublet removal
pdf(here(atac_aggr_prep,"plots","step1_cluster.pdf"))
p1 <- DimPlot(atacAggr_sub) + ggtitle("Harmony snATAC Prior to AMULET Doublet Removal")
AugmentPlot(p1)
dev.off()

# read in AMULET designated atac doublet metrics and filter barcodes with FDR < 0.05
amulet.df <- lapply(seq(orig.ident), function(index) {
  df <- fread(here(atac_aggr_prep, "doublets", orig.ident[index], "MultipletProbabilities.txt")) %>%
    as.data.frame()
  df$barcode_update <- paste(substr(df$barcode, 1, 16), index, sep = "-") # change barcode suffix to reflect gemgroup
  df <- dplyr::select(df, barcode_update, "p-value", "q-value") %>%
    dplyr::rename(amulet_pval = "p-value") %>%
    dplyr::rename(amulet_qval = "q-value") %>%
    dplyr::mutate(doublet = ifelse(amulet_qval < 0.05, 1, 0))
  return(df)
}) %>% 
bind_rows() %>%
column_to_rownames(var = "barcode_update")

# add amulet annotation to object. 1=doublet, 0=not a doublet
atacAggr_sub <- AddMetaData(atacAggr_sub, amulet.df)

# count doublets that were not already filtered and visualize
table(atacAggr_sub@meta.data$doublet)
pdf(here(atac_aggr_prep,"plots","step1_amulet.pdf"))
p1 <- FeaturePlot(atacAggr_sub, features="doublet", reduction="umap", order=TRUE) + ggtitle("Amulet Doublets")
AugmentPlot(p1)
dev.off()

# remove doublets
atacAggr_sub <- subset(
  x = atacAggr_sub,
  subset = doublet == 0
)

# recluster without doublets
atacAggr_sub <- FindNeighbors(object = atacAggr_sub, reduction = "harmony", dims = 2:30)
atacAggr_sub <- FindClusters(object = atacAggr_sub, verbose = TRUE, algorithm = 1) # Louvain algorithm
atacAggr_sub <- RunUMAP(object = atacAggr_sub, reduction = "harmony", dims = 2:30)

# visualize without doublets
pdf(here(atac_aggr_prep,"plots","step1_cluster_nodb.pdf"))
p1 <- DimPlot(atacAggr_sub) + ggtitle("Harmony snATAC After AMULET Doublet Removal")
AugmentPlot(p1)
dev.off()

# create gene activity matrix
gene.activities <- GeneActivity(atacAggr_sub)
atacAggr_sub[['RNA']] <- CreateAssayObject(counts = gene.activities)
atacAggr_sub <- NormalizeData(
  object = atacAggr_sub,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atacAggr_sub$nCount_RNA)
)

# perform label transfer from previously annotated snRNAseq object
rna_aggr_prep <- here("project","analysis","dkd","rna_aggr_prep")
rnaAggr = readRDS(here(rna_aggr_prep,"step2_anno.rds"))
Idents(rnaAggr) <- "celltype"  

transfer.anchors <- FindTransferAnchors(
  reference = rnaAggr,
  query = atacAggr_sub,
  reduction = 'cca',
  normalization.method = 'SCT',
  reference.assay = 'SCT',
  query.assay = 'RNA'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rnaAggr$celltype,
  weight.reduction = atacAggr_sub[['lsi']],
  dims = 2:30
)

atacAggr_sub <- AddMetaData(object = atacAggr_sub, metadata = predicted.labels)

# visualize label transfer results
pdf(here(atac_aggr_prep,"plots","step1_label.pdf"))
p1 <- DimPlot(atacAggr_sub, label=TRUE, group.by="predicted.id") + ggtitle("snATAC predicted celltypes")
AugmentPlot(p1)
p2 <- DimPlot(atacAggr_sub, label=TRUE, group.by="seurat_clusters") + ggtitle("snATAC seurat clusters")
AugmentPlot(p2)
dev.off()

# visualize marker genes using gene activity
marker.genes <- c("CUBN","HAVCR1","SLC5A1","SLC5A2","VCAM1","PROM1", # PT and PT-VCAM1+ markers
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

DefaultAssay(atacAggr_sub) <- "RNA"
print("Drawing UMAP markers")
pdf(here(atac_aggr_prep,"plots","step1_markers.pdf")) 
    lapply(marker.genes, function(gene) {
      tryCatch({plot <- FeaturePlot(atacAggr_sub, features=gene, reduction="umap")
                        AugmentPlot(plot) # downsample
                }, warning=function(w) return(NULL)) 
    })
dev.off()

saveRDS(atacAggr_sub, file=here(atac_aggr_prep, "step1_prep.rds"), compress=FALSE)
print("Step 1 finished in: ")
Sys.time() - start
