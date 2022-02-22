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
library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


figures <- here("project","analysis","dkd","figures")
dir.create(here(figures), showWarnings=FALSE)

####################################
# cell-specific partitioned heritability
files <- list.files(here("project","analysis","dkd","ldsc","partition","celltype_markers"), pattern = "*.cell_type_results.txt", full.names = TRUE)
res <- lapply(files, function(file) {
  df <- fread(file)
  df$trait <- str_split(basename(file), pattern = "_", simplify = TRUE)[,1]
  df$padj <- p.adjust(df$Coefficient_P_value, method = "BH")
  return(df)
  }) %>% bind_rows()

# plot each gwas trait result
for(gwas_trait in unique(res$trait)) {
  toplot <- dplyr::filter(res, trait == gwas_trait)

  # reorder the idents
  idents <- c("PCT","PST","PT_VCAM1","PT_CD36","PEC","ATL","TAL1","TAL2","DCT1","DCT2","PC","ICA","ICB","PODO","ENDO","FIB_VSMC_MC","TCELL","BCELL","MONO")
  toplot$celltype <- factor(toplot$Name, levels = rev(idents))

  pval_threshold = -log10(0.05)

  # adjust pval and calculate pval threshold
  p1 <- ggplot(toplot, aes(x=celltype, y=-log10(padj), fill=ifelse(padj < 0.05,1,0))) + 
    geom_bar(stat="identity") + 
    geom_hline(yintercept = pval_threshold) +
    RotatedAxis() +
    ggtitle(gwas_trait) +
    theme_bw() +
    guides(fill = "none") +
    xlab("") +
    ylim(0,4) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1),
          text = element_text(size = 20),
          legend.title = element_blank()) +
    coord_flip()

  pdf(here(figures, paste0(gwas_trait,".celltype.pdf")))
  print(p1)
  dev.off()
}
##################################
# cell-specific DAR partitioned heritability
files <- list.files(here("project","analysis","dkd","ldsc","partition","celltype_diabetes_dar"), pattern = "*.cell_type_results.txt", full.names = TRUE)
res <- lapply(files, function(file) {
  df <- fread(file)
  df$trait <- str_split(basename(file), pattern = "_", simplify = TRUE)[,1]
  df$padj <- p.adjust(df$Coefficient_P_value, method = "BH")
  return(df)
  }) %>% bind_rows()

# remove the merged dar category
res <- dplyr::filter(res, Name != "merge")

# plot each gwas trait result
for(gwas_trait in unique(res$trait)) {
  toplot <- dplyr::filter(res, trait == gwas_trait)

  # reorder the idents
  idents <- c("PCT","PST","PT_VCAM1","PT_CD36","PEC","ATL","TAL1","TAL2","DCT1","DCT2","PC","ICA","ICB","PODO","ENDO","FIB_VSMC_MC","TCELL","BCELL","MONO")
  toplot$celltype <- factor(toplot$Name, levels = rev(idents))

  pval_threshold = -log10(0.05)

  # adjust pval and calculate pval threshold
  p1 <- ggplot(toplot, aes(x=celltype, y=-log10(padj), fill=ifelse(padj < 0.05,1,0))) + 
    geom_bar(stat="identity") + 
    geom_hline(yintercept = pval_threshold) +
    RotatedAxis() +
    ggtitle(gwas_trait) +
    theme_bw() +
    guides(fill = "none") +
    xlab("") +
    ylim(0,4) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1),
          text = element_text(size = 20),
          legend.title = element_blank()) +
    coord_flip()

  pdf(here(figures, paste0(gwas_trait,".celltype_dar.pdf")))
  print(p1)
  dev.off()
}
#############################
# also see allele_specific/analyze_allele_counts.R asca_vaf2.pdf
df <- fread(here("project","analysis","dkd","allele_specific","allele_binomial.csv"))
png(here(figures, "figure4a.png"))
ggplot(df, aes(x=total_ref_count, y=total_alt_count)) + geom_point(aes(color=binomial_test)) + theme_bw()
dev.off()
#############################
# annotate peaks with an allele-specific effect 
df <- fread(here("project","analysis","dkd","allele_specific","glme_results.csv"))

# load annotation database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# filter by padj
df$padj <- p.adjust(df$df$p.value_exp, method = "BH")

# convert to granges
peaks.gr <- StringToGRanges(df$peak, sep = c("-","-"))

# annotate the peaks
peakAnno <- annotatePeak(peaks.gr, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)

# grab distance to nearest TSS
dist_dar_tss <- peakAnno@anno@elementMetadata@listData$distanceToTSS

# extract peak locations and take first word in string 
peakloc <- peakAnno@anno@elementMetadata@listData$annotation
peakloc <- word(peakloc, 1)

df <- data.frame(peak=df$peak, dist=dist_dar_tss, Estimate=df$estimate_exp, peakloc=peakloc, gene=df$gene, pval=df$p.value_exp, padj=df$padj)

# this will convert very large or very small estimates to +Inf or -Inf
# these are likely outliers and are removed from dataset
df <- dplyr::mutate(df, log_odds = log(exp(Estimate)))
df[is.infinite(df$log_odds),]$log_odds <- NA
df <- df[!is.na(df$log_odds),]

# add probability
df <- dplyr::mutate(df, odds = exp(Estimate)) %>%
  mutate(prob = odds / (1 + odds))
                   
# recode peak locations
df$peakloc <- recode_factor(df$peakloc, 
  Downstream = "Promoter", 
  Distal = "Distal_Intergenic")

# set factor levels for peak location
df$peakloc <- factor(df$peakloc, levels=c("Promoter","Intron","Exon","Distal_Intergenic","3'","5'"))

# modify any points that are greater than 100kb from TSS to fit in plot
df <- dplyr::mutate(df, distmod = ifelse(abs(dist) > 100000, sign(dist) * 100000, dist))

# draw histogram excluding outliers with large logit
limit <- 2 * sd(df$log_odds)
p0 <- ggplot(df, aes(x=log_odds)) + geom_histogram() + xlim(-limit, limit)

# draw a density plot and inspect for outliers
pdensity <- ggplot(df, aes(x=log_odds)) + geom_density() + xlim(-limit, limit)

# assign limit
limit <- 0.2

# create layers to plot DAR meeting threshold with different fill color
toplot <- dplyr::filter(df, padj < 0.05)

p1 <- ggplot() + 
  geom_point(data = toplot, aes(x=distmod, y=log_odds, color=peakloc)) +
  theme_bw() +
  geom_hline(yintercept=0) +
  xlab("Distance from TSS") +
  ylab("Log Odds (ALT Allele)") +
  scale_colour_brewer('DAR location', palette = "Set1") +
  ylim(-limit, limit)

# create layers to plot peaks meeting adjusted pval 
toplot1 <- dplyr::filter(df, pval < 0.05)
toplot2 <- dplyr::filter(df, pval >= 0.05)

top_layer <- ggplot() + 
  geom_point(data = toplot1, aes(x=distmod, y=log_odds, color=peakloc)) +
  geom_hline(yintercept=0) +
  theme_bw() +
  xlab("Distance from TSS") +
  ylab("Log Odds (ALT Allele)") +
  scale_colour_brewer('Peak Location', palette = "Set1") +
  ylim(-limit, limit)
p3 <- top_layer + geom_point(data = toplot2, size=1, aes(x=distmod, y=log_odds, fill="black"))

# plot the peaks that dont meet padj on bottom
p3$layers <- rev(p3$layers)
p3

# in inkscape select the plot area and then Path > Union to flatten geom_point objects
pdf(here(figures, "figure6d.pdf"))
p0
p1
p2
p3
dev.off()

##################################
# allele-specific PT partitioned heritability
files <- list.files(here("project","analysis","dkd","ldsc","partition","allele_specific"), pattern = "*.cell_type_results.txt", full.names = TRUE)
res <- lapply(files, function(file) {
  df <- fread(file)
  df$trait <- str_split(basename(file), pattern = "_", simplify = TRUE)[,1]
  df$padj <- p.adjust(df$Coefficient_P_value, method = "BH")
  return(df)
  }) %>% bind_rows()

# remove the merged dar category
res <- dplyr::filter(res, Name != "merge")

# plot each gwas trait result
for(gwas_trait in unique(res$trait)) {
  toplot <- dplyr::filter(res, trait == gwas_trait)

  # reorder the idents
  models <- c("asca","glme_results","glme_dkd_results","glme_dkdint_results")
  toplot$model <- factor(toplot$Name, levels = models)

  pval_threshold = -log10(0.05)

  # adjust pval and calculate pval threshold
  p1 <- ggplot(toplot, aes(x=model, y=-log10(padj), fill=ifelse(padj < 0.05,1,0))) + 
    geom_bar(stat="identity") + 
    geom_hline(yintercept = pval_threshold) +
    RotatedAxis() +
    ggtitle(gwas_trait) +
    theme_bw() +
    guides(fill = "none") +
    xlab("") +
    ylim(0,4) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1),
          text = element_text(size = 20),
          legend.title = element_blank()) +
    coord_flip()

  pdf(here(figures, paste0(gwas_trait,".pt_allele_specific.pdf")))
  print(p1)
  dev.off()
}

############################################
# read in the glme dkd interaction results
allele.df <- fread(here("project","analysis","dkd","allele_specific","glme_dkdint_results.csv"))
allele.df$padj <- p.adjust(allele.df$p.value_exp, method = "BH")
allele.df <- dplyr::filter(allele.df, p.value_exp < 0.05)

# convert to granges
allele.gr <- StringToGRanges(unique(allele.df$peak), sep = c("-","-"))

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

# filter for PT-specific DAR
dar.df <- dplyr::filter(dar.df, celltype %in% c("PCT","PST"))
dar.peaks <- unique(dar.df$peak)
dar.gr <- StringToGRanges(dar.peaks, sep = c("-","-"))
dar.gr <- sort(dar.gr)

# double check the intersection numbers
allele_int_dar <- join_overlap_intersect(allele.gr, dar.gr) # 104

# draw venn diagram
makeVennDiagram(list(allele.gr, dar.gr), 
                      NameOfPeaks = c("allele","dar"), 
                      fill = c("cornflowerblue", "lightgreen"),
                      alpha = 0.50)
pdf(here(figures, "figure7e.pdf"))
makeVennDiagram(list(allele.gr, dar.gr), 
                      NameOfPeaks = c("allele","dar"), 
                      fill = c("cornflowerblue", "lightgreen"),
                      alpha = 0.50)
dev.off()






