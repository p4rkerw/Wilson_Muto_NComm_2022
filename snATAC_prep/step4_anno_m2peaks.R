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
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log.rna.step4.out Rscript $SCRATCH1/dkd/snATAC_prep/step4_anno_m2peaks.R

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
library(tidyr)

atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr_macs2 <- readRDS(here(atac_aggr_prep,"step3_m2peaks.rds"))
DefaultAssay(atacAggr_macs2) <- "peaks"

# perform label transfer from previously annotated snRNAseq object
# to generate an imputed RNA assay in the macs2 object
rna_aggr_prep <- here("project","analysis","dkd","rna_aggr_prep")
rnaAggr = readRDS(here(rna_aggr_prep,"step2_anno.rds"))
Idents(rnaAggr) <- "celltype" # use the low resolution clustering for label transfer

# calculate cell-specific mean and sd for macs2_FRiP and filter bc that are 2 stdev below mean
cell_thresholds <- group_by(atacAggr_macs2@meta.data, predicted.id) %>%
                  dplyr::mutate(limit_macs2_FRiP = mean(macs2_FRiP) - 2*sd(macs2_FRiP)) %>%
                  dplyr::mutate(filter_macs2_FRiP = ifelse(macs2_FRiP < limit_macs2_FRiP, 1, 0)) %>%
                  as.data.frame()
atacAggr_macs2$filter_macs2_FRiP <- cell_thresholds$filter_macs2_FRiP

# create qc plots for macs2_FRiP
pdf(here(atac_aggr_prep,"plots","step4_qc.pdf"))
p1 <- DimPlot(atacAggr_macs2, label=TRUE, group.by="predicted.id") + ggtitle("snATAC predicted celltypes")
AugmentPlot(p1)
p2 <- DimPlot(atacAggr_macs2, label=TRUE, group.by="seurat_clusters") + ggtitle("snATAC seurat clusters")
AugmentPlot(p2)
FeaturePlot(atacAggr_macs2, features="macs2_FRiP")
FeaturePlot(atacAggr_macs2, features="mito_ratio")
FeaturePlot(atacAggr_macs2, features="filter_macs2_FRiP", order=TRUE)
dev.off()

# remove barcodes that are in the bottom 2.5% for macs2_FRiP for individual cell types
atacAggr_macs2 <- subset(atacAggr_macs2, subset = filter_macs2_FRiP == 0) # 1=filter, 0=keep   

# recluster
atacAggr_macs2 <- FindNeighbors(object = atacAggr_macs2, reduction = "harmony", dims = 2:30)
atacAggr_macs2 <- FindClusters(object = atacAggr_macs2, verbose = TRUE, algorithm = 1) # Louvain algorithm
atacAggr_macs2 <- RunUMAP(object = atacAggr_macs2, reduction = "harmony", dims = 2:30)

# redo label transfer celltype prediction with macs2 peaks
transfer.anchors <- FindTransferAnchors(
  reference = rnaAggr,
  query = atacAggr_macs2,
  reduction = 'cca',
  normalization.method = 'SCT',
  reference.assay = 'SCT',
  query.assay = 'RNA'
)

# update predicted cell type metadata
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rnaAggr$celltype,
  weight.reduction = atacAggr_macs2[['lsi']],
  dims = 2:30
)

atacAggr_macs2 <- AddMetaData(atacAggr_macs2, predicted.labels)

# create plots 
pdf(here(atac_aggr_prep,"plots","step4_label.pdf"))
p1 <- DimPlot(atacAggr_macs2, label=TRUE, group.by="predicted.id") + ggtitle("snATAC predicted celltypes")
AugmentPlot(p1)
p2 <- DimPlot(atacAggr_macs2, label=TRUE, group.by="seurat_clusters") + ggtitle("snATAC seurat clusters")
AugmentPlot(p2)
dev.off()

# compute correlation matrix between predicted.id groups to identify unrelated cell types
# eliminate cells with a predicted.id that does not match the most abundant predicted.id within a cluster (ie its unrelated)
DefaultAssay(atacAggr_macs2) <- "peaks"
atacAggr_macs2 <- FindTopFeatures(atacAggr_macs2, assay="peaks")
variable_peaks <- VariableFeatures(atacAggr_macs2)
av.access <- AverageExpression(atacAggr_macs2, features=variable_peaks, group.by="predicted.id", assays="peaks")
cor.access <- as.data.frame(cor(av.access$peaks))
cor.df <- rownames_to_column(cor.access) %>%
gather(celltype, value, -rowname)

# plot the correlation results
pdf(here(atac_aggr_prep,"plots","step4_correlation.pdf"))
p1 <- ggplot(cor.df, aes(x = rowname, y = celltype, fill = value)) + geom_tile()
print(p1)
dev.off()

# loop through each cluster and return barcodes for cells with a predicted.id that does not match the most
# abundant celltype in the cluster && are negatively correlated with that cluster (ie they are unrelated cell types. eg. ENDO-PT)
bc_anno <- lapply(unique(atacAggr_macs2@meta.data$seurat_clusters), function(cluster.sel){
    
    meta <- rownames_to_column(atacAggr_macs2@meta.data, var="barcode")

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

atacAggr_macs2 <- AddMetaData(atacAggr_macs2, bc_anno)

# visualize predicted doublets
pdf(here(atac_aggr_prep, "plots","step4_doublets.pdf"))
p1 <- FeaturePlot(atacAggr_macs2, features = "filter_barcode", order=TRUE)
AugmentPlot(p1)
dev.off()

# remove predicted doublets / low quality cells
atacAggr_macs2 <- subset(atacAggr_macs2, subset = filter_barcode == 0) # 0=keep, 1=filter

# recluster
atacAggr_macs2 <- FindNeighbors(object = atacAggr_macs2, reduction = "harmony", dims = 2:30)
atacAggr_macs2 <- FindClusters(object = atacAggr_macs2, verbose = TRUE, algorithm = 1) # Louvain algorithm
atacAggr_macs2 <- RunUMAP(object = atacAggr_macs2, reduction = "harmony", dims = 2:30)

# plots
pdf(here(atac_aggr_prep,"plots","step4_clusternodb.pdf"))
p1 <- DimPlot(atacAggr_macs2, label=TRUE, group.by="predicted.id") + ggtitle("snATAC predicted celltypes")
AugmentPlot(p1)
p2 <- DimPlot(atacAggr_macs2, label=TRUE, group.by="seurat_clusters") + ggtitle("snATAC seurat clusters")
AugmentPlot(p2)
dev.off()
           
# redo label transfer celltype prediction after bc filtering
transfer.anchors <- FindTransferAnchors(
  reference = rnaAggr,
  query = atacAggr_macs2,
  reduction = 'cca',
  normalization.method = 'SCT',
  reference.assay = 'SCT',
  query.assay = 'RNA'
)

# update predicted cell type metadata
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rnaAggr$celltype,
  weight.reduction = atacAggr_macs2[['lsi']],
  dims = 2:30
)

# impute "RNA" values using snRNA dataset as a reference to create pseudomultiomics cells
refdata <- GetAssayData(
  object = rnaAggr, 
  assay = "RNA", 
  slot = "data"
)

imputation <- TransferData(
  anchorset = transfer.anchors, 
  refdata = refdata, 
  weight.reduction = atacAggr_macs2[["lsi"]],
  dims = 2:30 
)
atacAggr_macs2[["IMPRNA"]] <- imputation

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
                  "ITGA8","PDGFRB", # VSMC, MES, JGA,
                  "MEIS2","PIEZO2","REN", # JGA, (mesangium is REN -ve) 
                  "CNN1", # VSMC
                  "ACTA2","CALD1", # Mesangium / Fibroblast markers
                  "PDGFRA","GATA3", # Fibroblasts / MES
                  "PROX1","FLT4","PDPN", # Lymphatics
                  "PTPRC","CD3E","MS4A1", # Lymphocytes
                  "SDC1", # CD138 plasma cells
                  "FCGR3A","CD14","CSF1R", # Monocyte / Macrophage
                  "KRT7","KRT20","UPK2") # Non-specific epithelial markers and urothelium

DefaultAssay(atacAggr_macs2) <- "IMPRNA"
print("Drawing UMAP markers")
pdf(here(atac_aggr_prep,"plots","step4_impmarkers.pdf")) 
    lapply(marker.genes, function(gene) {
      tryCatch({plot <- FeaturePlot(atacAggr_macs2, features=gene, reduction="umap")
                        AugmentPlot(plot) # downsample
                }, warning=function(w) return(NULL)) 
    })
dev.off()

# # plots
pdf(here(atac_aggr_prep,"plots","step4_cluster2.pdf"))
p1 <- DimPlot(atacAggr_macs2, label=TRUE, group.by="predicted.id") + ggtitle("snATAC predicted celltypes")
AugmentPlot(p1)
p2 <- DimPlot(atacAggr_macs2, label=TRUE, group.by="seurat_clusters") + ggtitle("snATAC seurat clusters")
AugmentPlot(p2)
dev.off()


atacAggr_macs2 <- RenameIdents(atacAggr_macs2,
  '0' = 'PCT',
  '1' = 'PCT',
  '2' = 'DCT1',
  '3' = 'PST',
  '4' = 'TAL1',
  '5' = 'PCT',
  '6' = 'DCT2',
  '7' = 'TAL2',
  '8' = 'PT_VCAM1',
  '9' = 'TAL1',
  '10' = 'ATL',
  '11' = 'ENDO',
  '12' = 'ICA',
  '13' = 'PC',
  '14' = 'ICB',
  '15' = 'PEC',
  '16' = 'FIB_VSMC_MC',
  '17' = 'MONO',
  '18' = 'PT_PROM1',
  '19' = 'TCELL',
  '20' = 'PODO',
  '21' = 'BCELL',
  '22' = 'PT_CD36')

# update levels
levels(atacAggr_macs2) <- c("PCT","PST","PT_VCAM1","PT_CD36","PEC","ATL",
          "TAL1","TAL2","DCT1","DCT2","PC",
          "ICA","ICB","PODO","ENDO","FIB_VSMC_MC",
          "TCELL","BCELL","MONO")


atacAggr_macs2@meta.data$celltype <- Idents(atacAggr_macs2)


pdf(here(atac_aggr_prep,"plots","step4_anno.pdf"))
p1 <- DimPlot(atacAggr_macs2, label=TRUE, group.by="celltype") + ggtitle("snATAC annotated celltypes")
AugmentPlot(p1)
p2 <- DimPlot(atacAggr_macs2, label=TRUE, group.by="predicted.id") + ggtitle("snATAC predicted.id")
AugmentPlot(p2)
dev.off()

# update metadata columns for disease groups
gemgroup <- sapply(strsplit(rownames(atacAggr_macs2@meta.data), split="-"), "[[", 2) 
current.gemgroups <- seq(n_distinct(gemgroup)) # no. gemgroups is no. samples
diabetes <- c(0,0,0,0,0,1,1,1,1,1,0,1,1) # 0=control, 1=diabetes

meta_group <- plyr::mapvalues(gemgroup, from = current.gemgroups, to = diabetes)
atacAggr_macs2   <- AddMetaData(object=atacAggr_macs2 , metadata=data.frame(diabetes=meta_group, 
                                                             row.names=rownames(atacAggr_macs2@meta.data)))

# create a list of cell type barcodes and write to file
dir.create(here("project","analysis","dkd","barcodes"), showWarnings=FALSE)
bc <- atacAggr_macs2@meta.data %>%
  dplyr::select(barcode, celltype, orig.ident) %>%
  dplyr::mutate(barcode = str_split(barcode, pattern="-", simplify=TRUE)[,1]) 
write.table(bc, file=here("project","analysis","dkd","barcodes","atac_barcodes.csv"), sep=",", row.names=FALSE, quote=FALSE)

DefaultAssay(atacAggr_macs2) <- "peaks"
saveRDS(atacAggr_macs2, file=here(atac_aggr_prep,"step4_m2anno.rds"), compress=FALSE)
