#!/usr/bin/env Rscript
# this script will annotate aggregated snRNA object 
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
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' -o $SCRATCH1/log.rna.step2.out Rscript $SCRATCH1/Wilson_Muto_NComm_2022/snRNA_prep/step2_anno.R

library(Seurat) # 4.0.3
library(ggplot2) # 3.3.5
library(dplyr) # 1.0.7
library(here) # 1.0.1
library(tibble) # 3.1.3
library(stringr)
set.seed(1234)

rna_aggr_prep <- here("project","analysis","dkd","rna_aggr_prep")
rnaAggr <- readRDS(here(rna_aggr_prep,"step1_prep.rds"))

########## Annotation of the clusters for data integration with snATAC dataset ##########
Idents(rnaAggr) <- "seurat_clusters"
rnaAggr[["orig.clusters"]] <- Idents(object = rnaAggr) # stash cluster idents prior to annotation
new.cluster.ids <- c("PT","PT","DCT1","TAL1","PC",
                     "DCT2","TAL2","ICA","PTVCAM1","PODO",
                     "PEC","DCT2","TAL1","ENDO","ATL",
                     "ICB","ENDO","MES","FIB","ENDO","LEUK")
names(new.cluster.ids) <- levels(rnaAggr)
rnaAggr <- RenameIdents(rnaAggr, new.cluster.ids)

# reorder the idents and save celltype annotations in celltype slot metadata
levels(rnaAggr) <- c("PT","PTVCAM1","PEC","ATL","TAL1",
					 "TAL2","DCT1","DCT2","PC","ICA",
					 "ICB","PODO","ENDO","MES","FIB",
					 "LEUK")
rnaAggr@meta.data$celltype <- rnaAggr@active.ident

# reset celltype as primary ident before saving and plotting
Idents(rnaAggr) <- "celltype"

# redraw umap and dotplot with reordered idents
p1 <- DimPlot(rnaAggr, reduction = "umap", label = TRUE) +
  ggtitle("snRNA Annotated Celltypes")

# draw pdf plots for before and after annotation
pdf(here(rna_aggr_prep,"plots","step2_anno.pdf"))
list(p1)
dev.off()

# create a list of cell type barcodes and write to file
dir.create(here("project","analysis","dkd","barcodes"), showWarnings=FALSE)
bc <- rownames_to_column(rnaAggr@meta.data, var = "barcode") %>%
  dplyr::select(barcode, celltype, orig.ident) %>%
  dplyr::mutate(barcode = str_split(barcode, pattern="-", simplify=TRUE)[,1]) 
write.table(bc, file=here("project","analysis","dkd","barcodes","rna_barcodes.csv"), sep=",", row.names=FALSE, quote=FALSE)

# save preprocessed rnaAggr file
print("Saving annotated snRNAseq object as step2_anno.rds in:")
saveRDS(rnaAggr, file = here(rna_aggr_prep,"step2_anno.rds"), compress=FALSE)
