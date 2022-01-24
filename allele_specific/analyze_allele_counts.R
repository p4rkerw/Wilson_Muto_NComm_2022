# run interactive
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

# run detached
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -n 10 -J "allele_mods" -o "$SCRATCH1/log.allele_mods" -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/sctools:R4.1.0)' \
# Rscript $SCRATCH1/dkd/allele_specific/analyze_allele_counts.R

# RUN LOCAL INTERACTIVE R terminal
# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v /mnt/g/diabneph:$HOME/project \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.0

library(data.table)
library(dplyr)
library(GenomicRanges)
library(plyranges)
library(tidyr)
library(stringr)
library(Matrix)
library(ggplot2)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Seurat)
library(future.apply)
library(tibble)
library(here)

# load aggregated atac dataset with peak-gene links precomputed
atacAggr <- readRDS(here("project","analysis","dkd","atac_aggr_prep","step7_links.rds"))

# load imputed RNA values
imprna.mat <- GetAssayData(object = atacAggr, assay = "IMPRNA", slot = "data") # this assay contains the imputed RNA values
imprna.mat <- as.data.frame(imprna.mat)

# load links
links.gr <- Links(atacAggr)

# allele specific counts dir
counts_dir <- here("project","analysis","dkd","wasp_atac","counts")

# read in single cell allele counts
samp_array <- c("Control_1","Control_2","Control_3","Control_4","Control_5",
                "DN_1","DN_2","DN_3","DN_4","DN_5","Control_6")
barcode_suffix_array <- c(1,2,3,4,5,6,7,8,9,10,11)
counts <- lapply(seq(samp_array), function(index){
  samp <- samp_array[index]
  file <- here(counts_dir,paste0(samp,".counts.single_cell.phased.table"))
  counts <- fread(file)
  counts$sample <- samp
  counts$barcode <- gsub("-1","", counts$barcode)
  counts$barcode <- paste0(counts$barcode,"-",barcode_suffix_array[index])
  return(counts)
}) %>% bind_rows()

# add cell type barcodes to df
bc <- data.frame(barcode=rownames(atacAggr@meta.data), celltype=atacAggr@meta.data$celltype)
counts <- left_join(counts, bc, by="barcode")

# filter proximal tubule cell types
counts <- dplyr::filter(counts, celltype %in% c("PCT","PST"))

# convert count matrix to granges
counts.bed <- data.frame(chrom=counts$contig, start=counts$position, end=counts$position)
counts.gr <- makeGRangesFromDataFrame(counts.bed)
counts.gr$variant_id <- counts$variant_id
counts.gr$GT <- counts$GT
counts.gr$refAllele <- counts$refAllele
counts.gr$altAllele <- counts$altAllele
counts.gr$refCount <- counts$refCount
counts.gr$altCount <- counts$altCount
counts.gr$sample <- counts$sample
counts.gr$barcode <- counts$barcode
rm(counts)
rm(counts.bed)

# annotate each peak with the linked gene
overlap.gr <- join_overlap_intersect(counts.gr, links.gr)
df <- as.data.frame(overlap.gr) %>%
  dplyr::rename(peak_id = peak)

# filter the count matrix
sub <- dplyr::group_by(df, peak_id) %>% mutate(total_ref_count = sum(refCount)) 
sub <- dplyr::group_by(sub, peak_id) %>% mutate(total_alt_count = sum(altCount)) 
sub <- dplyr::mutate(sub, ref=ifelse(refCount > 0, 1, 0), alt=ifelse(altCount > 0, 1, 0))
sub <- dplyr::mutate(sub, total_peak_counts = total_alt_count + total_ref_count)
sub <- mutate(sub, VAF=total_alt_count/(total_ref_count + total_alt_count))
sub <- dplyr::filter(sub, total_peak_counts > 20, total_ref_count > 5, total_alt_count > 5)

# add sample specific alt and total peak counts (not aggregated by peak)
sub <- dplyr::group_by(sub, sample, peak_id) %>% 
  mutate(sample_ref_count = sum(refCount)) %>% 
  mutate(sample_alt_count = sum(altCount))  
sub <- dplyr::mutate(sub, sample_peak_counts = sample_alt_count + sample_ref_count)

# enable parallel processing via future package
plan("multicore", workers=6) 
options(future.globals.maxSize = 640000 * 1024^2) # for 64 Gb RAM
plan()

# make histograms
toplot <- distinct(sub, peak_id, VAF)
hist(toplot$VAF, main="ATAC Peak Allele Fraction", xlab = "ATAC Peak Allele Fraction")

toplot <- distinct(sub, peak_id, total_ref_count, total_peak_counts, total_alt_count, VAF)
hist(toplot$VAF)

dir.create(here("project","analysis","dkd","allele_specific"), showWarnings=FALSE)
pdf(here("project","analysis","dkd","allele_specific","asca_vaf.pdf"))
ggplot(toplot, aes(x=total_ref_count, y=total_alt_count)) + geom_point()
dev.off()

# find binomial p values for aggregated peaks across samples
agg_pvals <- future_lapply(seq(nrow(toplot)), function(x){
    total_ref <- toplot[x, ]$total_ref_count
    total_count <- toplot[x, ]$total_peak_counts
    binom.test(total_ref, total_count, 0.5)$p.value
    })

toplot$pval <- agg_pvals

# visualize which peaks pass binomial test in scatter plot
toplot <- mutate(toplot, binomial_test=ifelse(pval < 0.05, "p < 0.05", "NS"))
pdf(here("project","analysis","dkd","allele_specific","asca_vaf2.pdf"))
ggplot(toplot, aes(x=total_ref_count, y=total_alt_count)) + geom_point(aes(color=binomial_test))
dev.off()

# save df to file
df <- as.data.frame(toplot)

# adjust binomial pval for multiple comparisons
df$padj <- p.adjust(df$pval, method="hochberg", n = length(df$pval))
fwrite(df, file=here("project","analysis","dkd","allele_specific","allele_binomial.csv"), row.names=FALSE, quote=FALSE, sep=",", col.names=TRUE)

# binarize the ref and alt allele counts
allele_peak_counts <- dplyr::select(sub, peak_id, sample, barcode, gene, ref, alt)
allele_peak_counts <- ungroup(allele_peak_counts)

# add variables for diabetes and sex to matrix
# update metadata columns for disease groups
diabetes <- c(0,0,0,0,0,1,1,1,1,1,0) # 0=control, 1=diabetes
sex_male <- c(1,1,0,1,0,1,1,0,1,0,0)
diabetes_group <- plyr::mapvalues(allele_peak_counts$sample, from = samp_array, to = diabetes)
allele_peak_counts$diabetes <- diabetes_group
sex_group <- plyr::mapvalues(allele_peak_counts$sample, from = samp_array, to = sex_male)
allele_peak_counts$sex <- sex_group

# save the allele counts to file
fwrite(allele_peak_counts, file=here("project","analysis","dkd","allele_specific","allele_counts.csv"), row.names=FALSE, quote=FALSE, sep=",", col.names=TRUE)

# filter the imputed RNA matrix for barcodes containing allele counts
imprna_allele.mat <- imprna.mat[colnames(imprna.mat) %in% unique(allele_peak_counts$barcode)]
rm(imprna.mat)
imprna_allele.mat <- rownames_to_column(imprna_allele.mat, var = "gene")
fwrite(imprna_allele.mat, file=here("project","analysis","dkd","allele_specific","imprna_allele.csv"), row.names=FALSE, quote=FALSE, sep=",", col.names=TRUE)
