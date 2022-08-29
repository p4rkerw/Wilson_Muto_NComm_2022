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
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log.rna.step3.out Rscript $SCRATCH1/Wilson_Muto_NComm_2022/snATAC_prep/step3_count_m2peaks.R

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
library(stringr) # 1.4.0
library(future) # 1.21.0
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)

# extract metadata to add to new object 
atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep, "step1_prep.rds"))
metadata <- atacAggr@meta.data
rm(atacAggr)

# set fragment file
aggr_input_dir <- here("project","analysis","dkd","cellranger_atac_aggr_13","outs")
fragment.path = here(aggr_input_dir,"fragments.tsv.gz")
barcodes <- rownames(metadata)

# read in macs2 peaks
peaks.gr <- fread(here(atac_aggr_prep,"step2_peaks.gr")) %>%
  makeGRangesFromDataFrame()

# enable parallel processing via future package
plan("multiprocess", workers = 6) 
options(future.globals.maxSize = 128000 * 1024^2) # for 128 Gb RAM
plan()

# create a new barcode-peak matrix using macs2 peaks
fragments <- CreateFragmentObject(
  fragment.path,
  cells = barcodes,
  validate.fragments = TRUE,
  verbose = TRUE
)

counts <- FeatureMatrix(
  fragments=fragments,
  features=peaks.gr,
  cells=barcodes
)

# cleanup
rm(fragments)

# Create new object and add metadata
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = here(aggr_input_dir,"fragments.tsv.gz"),
  min.cells = 10,
  min.features = 200
)

# cleanup
rm(counts)

atacAggr_macs2 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
rm(chrom_assay)

# compute the fraction of fragments within macs2 peaks (FRiP) to use as a latent
# variable in differential accessibility testing
fragment.counts <- CountFragments(fragment.path, cells=barcodes)
fragment.counts <- column_to_rownames(fragment.counts, var = "CB")
atacAggr_macs2 <- AddMetaData(atacAggr_macs2, fragment.counts)
rm(fragment.counts)

# the frequency_count metadata column is total fragments per cell
atacAggr_macs2 <- FRiP(
  object = atacAggr_macs2,
  assay = 'peaks',
  total.fragments = "frequency_count",
  col.name = 'macs2_FRiP'
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(atacAggr_macs2) <- annotations
rm(annotations)

# switch plan to avoid random number generation in parallel
plan("sequential") 
options(future.globals.maxSize = 128000 * 1024^2) # for 128 Gb RAM
plan()

# perform normalization and dimensional reduction of aggregated snATACseq object
atacAggr_macs2 <- RunTFIDF(atacAggr_macs2)
atacAggr_macs2 <- FindTopFeatures(atacAggr_macs2, min.cutoff = 'q0')
atacAggr_macs2 <- RunSVD(atacAggr_macs2)
# ElbowPlot(atacAggr_macs2, ndim = 40) # select number of dimensions for UMAP embedding
atacAggr_macs2 <- RunHarmony(atacAggr_macs2, group.by.vars = "orig.ident", reduction = "lsi", assay.use = "peaks", project.dim = FALSE)
atacAggr_macs2 <- FindNeighbors(object = atacAggr_macs2, reduction = "harmony", dims = 2:30)
atacAggr_macs2 <- FindClusters(object = atacAggr_macs2, verbose = TRUE, algorithm = 1) # Louvain algorithm
atacAggr_macs2 <- RunUMAP(object = atacAggr_macs2, reduction = "harmony", dims = 2:30)

# create gene activity matrix
gene.activities <- GeneActivity(atacAggr_macs2)
atacAggr_macs2[['RNA']] <- CreateAssayObject(counts = gene.activities)
atacAggr_macs2 <- NormalizeData(
  object = atacAggr_macs2,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atacAggr_macs2$nCount_RNA)
)
rm(gene.activities)

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

# plot clustering with macs2 peaks
pdf(here(atac_aggr_prep,"plots","step3_m2cluster.pdf"))
p1 <- DimPlot(atacAggr_macs2, label=TRUE) + ggtitle("snATAC Clustering with macs2 peaks")
AugmentPlot(p1)
dev.off()

DefaultAssay(atacAggr_macs2) <- "RNA"
print("Drawing UMAP markers")
pdf(here(atac_aggr_prep,"plots","step3_m2markers.pdf")) 
    lapply(marker.genes, function(gene) {
      tryCatch({plot <- FeaturePlot(atacAggr_macs2, features=gene, reduction="umap")
                        AugmentPlot(plot) # downsample
                }, warning=function(w) return(NULL)) 
    })
dev.off()

saveRDS(atacAggr_macs2, file=here(atac_aggr_prep,"step3_m2peaks.rds"), compress=FALSE)
