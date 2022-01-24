#!/usr/bin/env Rscript

# # export docker volumes
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"

# # to run interactive
# bsub -Is -G compute-parkerw -R 'rusage[mem=16GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

# SCRATCH1=/mnt/g/scratch
# docker run -it \
# -v /mnt/g/diabneph:$HOME/project \
# -v /mnt/g/reference:$HOME/reference \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# --workdir $HOME \
# p4rkerw/sctools:R4.1.0

library(here)
library(data.table)
library(GenomicRanges)
library(Signac)
library(openxlsx)
library(dplyr)

# set output directory
bed_dir <- here("project","analysis","dkd","ldsc","bed")
dir.create(bed_dir, showWarnings=FALSE)

# prepare cell-specific bed files. If there are less than 2000 sig regions include the top 2000 regions sorted first by pval
# and then by avg_log2FC > 0 (ie regions that are preferentially open in a cell type)
file <- here("project","analysis","dkd","markers","dar.macs2.celltype.markers.xlsx")
idents <- getSheetNames(file)
dar.grls <- lapply(idents, function(ident) {
	print(ident)
	df <- read.xlsx(file, sheet = ident)
	colnames(df)[1] <- "peak"
	df <- dplyr::filter(df, avg_log2FC > 0) %>%
		  dplyr::arrange(p_val, -avg_log2FC)
	# take top 2000 peaks if there are not enough peaks that meet padj threshold
	# if there are less than 2000 peaks with avg_log2FC > 0 then take all the peaks
	if(nrow(df[df$p_val_adj < 0.05,]) < 2000 && nrow(df) > 2000) {
		print("Taking top 2000 peaks")
		df <- df[1:2000, ] # else take all the rows
	} else if (nrow(df[df$p_val_adj < 0.05,]) > 2000) {
		# if there are more than 2000 peaks that meet padj only take those peaks
		print("Taking all peaks that meet padj threshold")
		df <- dplyr::filter(df, p_val_adj < 0.05)
	} else if (nrow(df) < 2000) {
		print("Taking all available peaks")
	}
	gr <- StringToGRanges(df$peak)
	# convert to 0-based coord
	start(gr) <- start(gr) - 1
	gr$p_val <- df$p_val
	gr$p_val_adj <-	df$p_val_adj
	gr$avg_log2FC <- df$avg_log2FC
	gr <- sort(gr)
	return(gr)
	})
# export gr as bed files
lapply(seq(idents), function(index) {
	df <- as.data.frame(dar.grls[[index]])
	fwrite(df, file=here(bed_dir, paste0(idents[index],".dar.macs2.celltype.markers.bed")),
	 quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	})

###############################################
# prepare diabetes cell-specific DAR from xlsx file
# select top 2000 cell-specific DAR sorted by padj and avg_log2FC
# here we will include DAR that show decreased accessibility
file <- here("project","analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx")
idents <- getSheetNames(file)
dar.grls <- lapply(idents, function(ident) {
	print(ident)
	df <- read.xlsx(file, sheet = ident)
	colnames(df)[1] <- "peak"
	df <- dplyr::arrange(df, p_val, avg_log2FC)
	# take top 2000 peaks if there are not enough peaks that meet padj threshold
	# if there are less than 2000 peaks that meet avg_log2FC > 0 then take all the peaks
	if(nrow(df[df$p_val_adj < 0.05,]) < 2000 && nrow(df) > 2000) {
		print("Taking top 2000 peaks")
		df <- df[1:2000, ] # else take all the rows
	} else if (nrow(df[df$p_val_adj < 0.05,]) > 2000) {
		# if there are more than 2000 peaks that meet padj only take those peaks
		print("Taking all peaks that meet padj threshold")
		df <- dplyr::filter(df, p_val_adj < 0.05)
	} else if (nrow(df) < 2000) {
		print("Taking all available peaks")
	}
	gr <- StringToGRanges(df$peak)
	# convert to 0-based coord
	start(gr) <- start(gr) - 1
	gr$p_val <-	df$p_val
	gr$p_val_adj <- df$p_val_adj
	gr$avg_log2FC <- df$avg_log2FC
	gr <- sort(gr)
	return(gr)
	})
# write cell-specific dar to bed files
lapply(seq(idents), function(index) {
	df <- as.data.frame(dar.grls[[index]]) 
	fwrite(df, file=here(bed_dir, paste0(idents[index],".dar.macs2.celltype.diab_vs_ctrl.bed")),
	 quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	})
# here we will take all unique cell-specific DAR that meet padj < 0.05
merge.dar.gr <- lapply(seq(dar.grls), function(index) {
	df <- as.data.frame(dar.grls[[index]]) %>%
		dplyr::filter(p_val_adj < 0.05) %>%
		dplyr::mutate(peak = paste0(seqnames,"-",start,"-",end))
	return(df$peak)	
	}) %>% 
	unlist() %>%
	unique() %>%
	StringToGRanges() %>%
	sort()
merge.dar.df <- as.data.frame(merge.dar.gr) 
fwrite(merge.dar.df, file=here(bed_dir, "merge.dar.macs2.celltype.diab_vs_ctrl.bed"),
	quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# prepare allele-specific effect bed from glmer model
# this will be a combination of genomic ranges that show allele-specific chromatin accessibility
# and ranges that show an allele-specic effect on a nearby gene
asca_file <- here("project","analysis","dkd","allele_specific","allele_binomial.csv")
asca.df <- fread(asca_file) %>%
		dplyr::filter(padj < 0.05)
asca.gr <- StringToGRanges(asca.df$peak_id, sep=c(":","-")) %>% 
	sort()
asca.df <- as.data.frame(asca.gr) %>%
	dplyr::mutate(start = start - 1) # convert to 0-based coordinates
fwrite(asca.df, file=here(bed_dir,"asca.bed"),
	 quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

glmer_file <- here("project","analysis","dkd","allele_specific","glme_dkdint_results.csv")
glmer.df <- fread(glmer_file) %>%
		dplyr::filter(p.value_exp < 0.05)
glmer.gr <- StringToGRanges(glmer.df$peak, sep=c(":","-")) %>% 
	sort()
glmer.df <- as.data.frame(glmer.gr) %>%
	dplyr::mutate(start = start - 1) # convert to 0-based coordinates
fwrite(glmer.df, file=here(bed_dir,"glme_dkdint_results.bed"),
	 quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

glmer_file <- here("project","analysis","dkd","allele_specific","glme_dkd_results.csv")
glmer.df <- fread(glmer_file) %>%
		dplyr::filter(p.value_exp < 0.05)
glmer.gr <- StringToGRanges(glmer.df$peak, sep=c(":","-")) %>% 
	sort()
glmer.df <- as.data.frame(glmer.gr) %>%
	dplyr::mutate(start = start - 1) # convert to 0-based coordinates
fwrite(glmer.df, file=here(bed_dir,"glme_dkd_results.bed"),
	 quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

glmer_file <- here("project","analysis","dkd","allele_specific","glme_results.csv")
glmer.df <- fread(glmer_file) %>%
		dplyr::filter(p.value_exp < 0.05)
glmer.gr <- StringToGRanges(glmer.df$peak, sep=c(":","-")) %>% 
	sort()
glmer.df <- as.data.frame(glmer.gr) %>%
	dplyr::mutate(start = start - 1) # convert to 0-based coordinates
fwrite(glmer.df, file=here(bed_dir,"glme_results.bed"),
	 quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)











