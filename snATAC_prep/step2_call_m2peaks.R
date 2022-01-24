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
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log_step2.out Rscript $SCRATCH1/dkd/snATAC_prep/step2_call_m2peaks.R

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
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table) # 1.14.0
library(future) # 1.21.0
library(tidyr)


atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step1_prep.rds"))

# compute correlation matrix between label transfer predicted.id using variable peaks
DefaultAssay(atacAggr) <- "peaks"
atacAggr <- FindTopFeatures(atacAggr, assay="peaks")
variable_peaks <- VariableFeatures(atacAggr)
av.access <- AverageExpression(atacAggr, features=variable_peaks, group.by="predicted.id", assays="peaks")
cor.access <- as.data.frame(cor(av.access$peaks))
cor.df <- rownames_to_column(cor.access) %>%
gather(celltype, value, -rowname)

# plot the correlation results
pdf(here(atac_aggr_prep,"plots","step2_correlation.pdf"))
p1 <- ggplot(cor.df, aes(x = rowname, y = celltype, fill = value)) + geom_tile()
print(p1)
dev.off()


# loop through each cluster and return barcodes for cells with predicted.id that do not match the most
# abundant celltype in the cluster && are negatively correlated with that cluster
bc_anno <- lapply(unique(atacAggr@meta.data$seurat_clusters), function(cluster.sel){
    
    meta <- rownames_to_column(atacAggr@meta.data, var="barcode")

    # cells that are predicted to be a different cell type within cluster
    meta <- dplyr::filter(meta, seurat_clusters == cluster.sel) %>%
      add_count(predicted.id, name="num_predicted.id_in_cluster", sort=TRUE) %>%
      mutate(most_abundant_in_cluster = predicted.id[1]) %>%
      mutate(not_most_abundant_cell = ifelse(predicted.id != most_abundant_in_cluster, 1, 0)) 

    correlated_celltypes <- dplyr::filter(cor.df, celltype == unique(meta$most_abundant_in_cluster)) %>%
      dplyr::filter(value > 0) 
    meta$related_celltype <- ifelse(meta$predicted.id %in% correlated_celltypes$rowname, 1, 0)
    meta$filter_barcode <- ifelse(meta$related_celltype == 0 & meta$not_most_abundant_cell == 1, 1, 0)

    # print cell types filtered from cluster
    print(paste0("Filtering for cluster: ", cluster.sel))
    print(table(meta$predicted.id, meta$filter_barcode))

    rownames(meta) <- meta$barcode
    return(meta)
}) %>% bind_rows() 

atacAggr <- AddMetaData(atacAggr, bc_anno)

# visualize predicted doublets
pdf(here(atac_aggr_prep, "plots","step2_doublets.pdf"))
p1 <- FeaturePlot(atacAggr, features = "filter_barcode", order=TRUE)
AugmentPlot(p1)
dev.off()

# remove predicted doublets / low quality cells
atacAggr <- subset(atacAggr, subset = filter_barcode == 0)

# recluster
plan("sequential")
atacAggr <- FindNeighbors(object = atacAggr, reduction = "harmony", dims = 2:30)
atacAggr <- FindClusters(object = atacAggr, verbose = TRUE, algorithm = 1) # Louvain algorithm
atacAggr <- RunUMAP(object = atacAggr, reduction = "harmony", dims = 2:30)

# plots
pdf(here(atac_aggr_prep,"plots","step2_clusternodb.pdf"))
p1 <- DimPlot(atacAggr, label=TRUE, group.by="predicted.id") + ggtitle("snATAC predicted celltypes")
AugmentPlot(p1)
p2 <- DimPlot(atacAggr, label=TRUE, group.by="seurat_clusters") + ggtitle("snATAC seurat clusters")
AugmentPlot(p2)
dev.off()

# visualize marker genes using gene activity
marker.genes <- c("CUBN","HAVCR1","SLC5A1","SLC5A2","VCAM1", # PT and PT-VCAM1+ markers
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

DefaultAssay(atacAggr) <- "RNA"
print("Drawing UMAP markers")
pdf(here(atac_aggr_prep,"plots","step2_markers.pdf")) 
    lapply(marker.genes, function(gene) {
      tryCatch({plot <- FeaturePlot(atacAggr, features=gene, reduction="umap")
                        AugmentPlot(plot) # downsample
                }, warning=function(w) return(NULL)) 
    })
dev.off()

# annotate the object
atacAggr <- RenameIdents(atacAggr,
  '0' = 'PCT',
  '1' = 'PCT',
  '2' = 'PCT',
  '3' = 'DCT1',
  '4' = 'TAL1',
  '5' = 'PST',
  '6' = 'TAL2',
  '7' = 'PC',
  '8' = 'PTVCAM1',
  '9' = 'TAL1',
  '10' = 'ENDO',
  '11' = 'DCT2',
  '12' = 'ICA',
  '13' = 'ATL',
  '14' = 'ICB',
  '15' = 'PEC',
  '16' = 'PST',
  '17' = 'LYMPH',
  '18' = 'PCT',
  '19' = 'MESFIB',
  '20' = 'MONO',
  '21' = 'PODO')

atacAggr@meta.data$celltype <- Idents(atacAggr)

pdf(here(atac_aggr_prep,"plots","step2_anno.pdf"))
p1 <- DimPlot(atacAggr, label=TRUE, group.by="celltype") + ggtitle("snATAC annotated celltypes")
AugmentPlot(p1)
p2 <- DimPlot(atacAggr, label=TRUE, group.by="seurat_clusters") + ggtitle("snATAC seurat clusters")
AugmentPlot(p2)
dev.off()

# enable parallel processing via future package
plan("multiprocess", workers = 8) 
options(future.globals.maxSize = 128000 * 1024^2) # for 128 Gb RAM
plan()

# call macs2 peaks with celltype annotations
print("Calling macs2 peaks")
dir.create(here(atac_aggr_prep,"macs2"), showWarnings=FALSE)
DefaultAssay(atacAggr) <- "peaks"
peaks.gr <- CallPeaks(
  object = atacAggr,
  group.by = "celltype",
  macs2.path = "/usr/local/bin/macs2",
  outdir = here(atac_aggr_prep,"macs2"),
  cleanup = FALSE
)

# do not keep peaks in alt contigs
# this will make the final object compatible with the chromvar workflow
chroms.keep <- seqnames(BSgenome.Hsapiens.UCSC.hg38)[1:25]
peaks.gr <- peaks.gr[seqnames(peaks.gr) %in% chroms.keep]

# write peaks to file
fwrite(as.data.frame(peaks.gr), file=here(atac_aggr_prep,"step2_peaks.gr"), sep="\t")

# save the cellranger rds
saveRDS(atacAggr, file=here(atac_aggr_prep,"step2_cranno.rds"), compress=FALSE)
