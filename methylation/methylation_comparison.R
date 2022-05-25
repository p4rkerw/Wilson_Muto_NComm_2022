# this script compares publicly-available differentially methylated loci to the cell-specific DAR for DKD vs. control samples
library(tidyr)
library(dplyr)
library(openxlsx)
library(GenomicRanges)
library(plyranges)
library(tibble)
library(rtracklayer)
library(Signac)
library(data.table)


# import liftover chains
ch19 = import.chain("G:/diabneph/analysis/dkd/methylation/hg19ToHg38.over.chain")
ch18 = import.chain("G:/diabneph/analysis/dkd/methylation/hg18ToHg38.over.chain")

# read in DAR
file <- here("project","analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx")
idents <- getSheetNames(file)
dar.df <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    rownames_to_column(var = "peak")  
  if(nrow(df) > 0) { df$celltype <- ident }
    return(df)
}) %>% bind_rows()

dar.df$overlap_atac_peak <- dar.df$peak
dar.df <- tidyr::separate(dar.df, col = peak, into = c("chrom","start","end"), sep = "-")
dar.gr <- makeGRangesFromDataFrame(dar.df, keep.extra.columns=TRUE)

# read in GR CUT&RUN peaks for bulk kidney
cut.gr <- fread(here("analysis","dkd","cut_and_run","kidney_GR_peaks.narrowPeak")) %>%
  dplyr::rename(chrom = V1, start = V2, end = V3) %>%
  dplyr::mutate(gr_peak = paste0(chrom,"-",start,"-",end)) %>%
  dplyr::select(chrom, start, end, gr_peak) %>%
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)

# Kidney cytosine methylation changes improve renal function decline estimation in patients with diabetic kidney disease
# hg19
# PMID: 31165727
# Supplemental Table 1
# Description: Top 65 probes significantly associated with interstitial fibrosis and passed independent replication. (a) Association determined by linear regression models adjusted for age, gender, race, diabetes, hypertension, batch, bisulfite conversion, and degree of lymphocytic infiltrate on histology.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6549146/bin/41467_2019_10378_MOESM4_ESM.xlsx

# read in table and convert col to hg19 coordinates
study_id <- "pmid31165727"
phenotype <- "interstitial_fibrosis_dkd"
dmr.df <- read.xlsx(here("analysis","dkd","methylation","41467_2019_10378_MOESM4_ESM.xlsx")) %>%
  dplyr::mutate(dmr = paste0(Chr,"-",Position,"-",Position)) %>%
  dplyr::select(dmr)
dmr.df$study_id <- study_id
dmr.df$phenotype <- phenotype
dmr.df <- tidyr::separate(dmr.df, col = dmr, into = c("chrom","start","end"), sep = "-", remove=FALSE)
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_flank.gr <- Extend(dmr38.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over1.gr <- join_overlap_intersect(dar.gr, dmr38.gr)
over1_flank.gr <- join_overlap_intersect(dar.gr, dmr38_flank.gr)

# overlap with GR cut and run peaks
cut1.gr <- join_overlap_intersect(cut.gr, dmr38.gr)
cut1_flank.gr <- join_overlap_intersect(cut.gr, dmr38_flank.gr)

# Kidney cytosine methylation changes improve renal function decline estimation in patients with diabetic kidney disease
# hg19
# PMID: 31165727
# File Name: Supplementary Data 3
# Description: Top 471 probes that improve model of kidney function declineusing weighted regression. (a) Model is a weighted linear regression model of adjusted eGFR slope (weight = inverse variance of adjusted eGFR slope). Base model includes variables: baseline eGFR, Diabetes, and Age (base model AIC = 206). When methylation level of probe is added to base model, the following variables are also added: methylation batch, and bisulfite conversion efficiency. (b) Association determined by linear regression models adjusted for age, gender, race, diabetes, hypertension, batch, bisulfite conversion, and degree of lymphocytic infiltrate on histology.
# read in table and convert col to hg19 coordinates
study_id <- "pmid31165727"
phenotype <- "renal_function_decline_dkd"
dmr.df <- read.xlsx(here("analysis","dkd","methylation","41467_2019_10378_MOESM6_ESM.xlsx")) %>%
  dplyr::mutate(dmr = paste0(chr,"-",Position,"-",Position)) %>%
  dplyr::select(dmr)
dmr.df$study_id <- study_id
dmr.df$phenotype <- phenotype
dmr.df <- tidyr::separate(dmr.df, col = dmr, into = c("chrom","start","end"), sep = "-", remove=FALSE)
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_flank.gr <- Extend(dmr38.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over2.gr <- join_overlap_intersect(dar.gr, dmr38.gr)
over2_flank.gr <- join_overlap_intersect(dar.gr, dmr38_flank.gr)

# overlap with GR cut and run peaks
cut2.gr <- join_overlap_intersect(cut.gr, dmr38.gr)
cut2_flank.gr <- join_overlap_intersect(cut.gr, dmr38_flank.gr)

# Assessment of differentially methylated loci in individuals with end-stage kidney disease attributed to diabetic kidney disease: an exploratory study
# PMID: 33933144
# hg19
# Supplementary Table ST3: Top-ranked differentially methylated CpG sites for matched individuals with T1DM-ESKD (n=107) vs. T1DM (n=107): FDR p≤x10-8 (Analysis 1)
study_id <- "pmid33933144"
phenotype <- "eskd_t1dm"
dmr.df <- read.xlsx(here("analysis","dkd","methylation","13148_2021_1081_MOESM1_ESM.xlsx"), sheet = "ST3", startRow = 2) %>%
    dplyr::mutate(dmr = paste0(CHR,"-",MAPINFO,"-",MAPINFO)) %>%
    dplyr::select(dmr)
dmr.df$study_id <- study_id
dmr.df$phenotype <- phenotype
dmr.df <- tidyr::separate(dmr.df, col = dmr, into = c("chrom","start","end"), sep = "-", remove=FALSE) %>%
  dplyr::filter(chrom != "NA") %>%
  mutate(chrom = paste0("chr",chrom)) %>%
  mutate(dmr = paste0("chr",dmr)) 
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_flank.gr <- Extend(dmr38.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over3.gr <- join_overlap_intersect(dar.gr, dmr38.gr)
over3_flank.gr <- join_overlap_intersect(dar.gr, dmr38_flank.gr)

# overlap with GR cut and run peaks
cut3.gr <- join_overlap_intersect(cut.gr, dmr38.gr)
cut3_flank.gr <- join_overlap_intersect(cut.gr, dmr38_flank.gr)

# Assessment of differentially methylated loci in individuals with end-stage kidney disease attributed to diabetic kidney disease: an exploratory study
# PMID: 33933144
# hg19
# Supplementary Table ST4:Top-ranked differentially methylated CpG sites for matched individuals with T1DM-ESKD (n=107) vs. T1DM (n=107): FDR p≤x10-8 and FC ±2 (Analysis 1)
study_id <- "pmid33933144"
phenotype <- "eskd_t1dm_fc2"
dmr.df <- read.xlsx(here("analysis","dkd","methylation","13148_2021_1081_MOESM1_ESM.xlsx"), sheet = "ST4", startRow = 2) %>%
    dplyr::mutate(dmr = paste0(CHR,"-",MAPINFO,"-",MAPINFO)) %>%
    dplyr::select(dmr)
dmr.df$study_id <- study_id
dmr.df$phenotype <- phenotype
dmr.df <- tidyr::separate(dmr.df, col = dmr, into = c("chrom","start","end"), sep = "-", remove=FALSE) %>%
  dplyr::filter(chrom != "NA") %>%
  mutate(chrom = paste0("chr",chrom)) %>%
  mutate(dmr = paste0("chr",dmr)) 
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_flank.gr <- Extend(dmr38.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over4.gr <- join_overlap_intersect(dar.gr, dmr38.gr)
over4_flank.gr <- join_overlap_intersect(dar.gr, dmr38_flank.gr)

# overlap with GR cut and run peaks
cut4.gr <- join_overlap_intersect(cut.gr, dmr38.gr)
cut4_flank.gr <- join_overlap_intersect(cut.gr, dmr38_flank.gr)

# Assessment of differentially methylated loci in individuals with end-stage kidney disease attributed to diabetic kidney disease: an exploratory study
# PMID: 33933144
# hg19
# Supplementary Table ST11: Top-ranked differentially methylated CpG sites for individuals with T1DM-ESKD - either dialysis or transplant (n=107) vs. T1DM (n=253): FDR p≤x10-8 (Analysis 2)
study_id <- "pmid33933144"
phenotype <- "eskd_t1dm_dialysis_transplant"
dmr.df <- read.xlsx(here("analysis","dkd","methylation","13148_2021_1081_MOESM1_ESM.xlsx"), sheet = "ST11", startRow = 2) %>%
    dplyr::mutate(dmr = paste0(CHR,"-",MAPINFO,"-",MAPINFO)) %>%
    dplyr::select(dmr)
dmr.df$study_id <- study_id
dmr.df$phenotype <- phenotype
dmr.df <- tidyr::separate(dmr.df, col = dmr, into = c("chrom","start","end"), sep = "-", remove=FALSE) %>%
  dplyr::filter(chrom != "NA") %>%
  mutate(chrom = paste0("chr",chrom)) %>%
  mutate(dmr = paste0("chr",dmr)) 
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_flank.gr <- Extend(dmr38.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over5.gr <- join_overlap_intersect(dar.gr, dmr38.gr)
over5_flank.gr <- join_overlap_intersect(dar.gr, dmr38_flank.gr)

# overlap with GR cut and run peaks
cut5.gr <- join_overlap_intersect(cut.gr, dmr38.gr)
cut5_flank.gr <- join_overlap_intersect(cut.gr, dmr38_flank.gr)

# Cytosine methylation changes in enhancer regions of core pro-fibrotic genes characterize kidney fibrosis development
# PMID: 24098934
# hg18
# Supplemental Table 2 List of differentially methylated regions with p-value methylation level, genomic locus and nearest annotated transcript
study_id <- "pmid24098934"
phenotype <- "interstitial_fibrosis_ckd"
dmr.df <- read.xlsx(here("analysis","dkd","methylation","gb-2013-14-10-r108-S3.xlsx"), startRow = 3) %>%
  dplyr::select(chr, start, end) %>%
  dplyr::mutate(dmr = paste0("chr",chr,"-",start,"-",end))
dmr.df$study_id <- study_id
dmr.df$phenotype <- phenotype
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38.gr <- liftOver(dmr.gr, ch18) %>% unlist()
dmr38_flank.gr <- Extend(dmr38.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over6.gr <- join_overlap_intersect(dar.gr, dmr38.gr)
over6_flank.gr <- join_overlap_intersect(dar.gr, dmr38_flank.gr)

# overlap with GR cut and run peaks
cut6.gr <- join_overlap_intersect(cut.gr, dmr38.gr)
cut6_flank.gr <- join_overlap_intersect(cut.gr, dmr38_flank.gr)

# Systematic integrated analysis of genetic and epigenetic variation in diabetic kidney disease
# PMID: 33144501
# hg19
# Supplemental Table 3: DNA methylation sites associated with albuminuria
study_id <- "pmid33144501"
phenotype <- "dkd_albuminuria"
dmr.df <- read.xlsx(here("analysis","dkd","methylation","pnas.2005905117.sd03.xlsx"), startRow = 3) %>%
  dplyr::select(chromosome, `Position.(b37)`) %>%
  dplyr::rename(chrom = chromosome) %>%
  dplyr::rename(start = `Position.(b37)`) %>%
  dplyr::mutate(end = start) %>%
  dplyr::mutate(dmr = paste0("chr",chrom,"-",start,"-",end))
dmr.df$study_id <- study_id
dmr.df$phenotype <- phenotype
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_flank.gr <- Extend(dmr38.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over7.gr <- join_overlap_intersect(dar.gr, dmr38.gr)
over7_flank.gr <- join_overlap_intersect(dar.gr, dmr38_flank.gr)

# overlap with GR cut and run peaks
cut7.gr <- join_overlap_intersect(cut.gr, dmr38.gr)
cut7_flank.gr <- join_overlap_intersect(cut.gr, dmr38_flank.gr)

# Systematic integrated analysis of genetic and epigenetic variation in diabetic kidney disease
# PMID: 33144501
# hg19
# Supplemental Table 3: DNA methylation sites associated with albuminuria
study_id <- "pmid33144501"
phenotype <- "dkd_eGFR"
dmr.df <- read.xlsx(here("analysis","dkd","methylation","pnas.2005905117.sd04.xlsx"), startRow = 3) %>%
  dplyr::select(chromosome, `Position.(b37)`) %>%
  dplyr::rename(chrom = chromosome) %>%
  dplyr::rename(start = `Position.(b37)`) %>%
  dplyr::mutate(end = start) %>%
  dplyr::mutate(dmr = paste0("chr",chrom,"-",start,"-",end))
dmr.df$study_id <- study_id
dmr.df$phenotype <- phenotype
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_flank.gr <- Extend(dmr38.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over8.gr <- join_overlap_intersect(dar.gr, dmr38.gr)
over8_flank.gr <- join_overlap_intersect(dar.gr, dmr38_flank.gr)

# overlap with GR cut and run peaks
cut8.gr <- join_overlap_intersect(cut.gr, dmr38.gr)
cut8_flank.gr <- join_overlap_intersect(cut.gr, dmr38_flank.gr)

# Systematic integrated analysis of genetic and epigenetic variation in diabetic kidney disease
# PMID: 33144501
# hg19
# Supplemental Table 3: DNA methylation sites associated with albuminuria
study_id <- "pmid33144501"
phenotype <- "dkd_eGFR"
dmr.df <- read.xlsx(here("analysis","dkd","methylation","pnas.2005905117.sd05.xlsx"), startRow = 3) %>%
  dplyr::select(chromosome, `Position.(b37)`) %>%
  dplyr::rename(chrom = chromosome) %>%
  dplyr::rename(start = `Position.(b37)`) %>%
  dplyr::mutate(end = start) %>%
  dplyr::mutate(dmr = paste0("chr",chrom,"-",start,"-",end))
dmr.df$study_id <- study_id
dmr.df$phenotype <- phenotype
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_flank.gr <- Extend(dmr38.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over9.gr <- join_overlap_intersect(dar.gr, dmr38.gr)
over9_flank.gr <- join_overlap_intersect(dar.gr, dmr38_flank.gr)

# overlap with GR cut and run peaks
cut9.gr <- join_overlap_intersect(cut.gr, dmr38.gr)
cut9_flank.gr <- join_overlap_intersect(cut.gr, dmr38_flank.gr)

# DNA hypermethylation and DNA hypomethylation is present at different loci in chronic kidney disease
# PMID: 24253112
# hg19
# Supplemental Table 2: Top ranked markers in genes of interest where multiple CpG’s affected per gene
# table reformatted from pdf
study_id <- "pmid24253112"
phenotype <- "ckd"
dmr.df <- read.csv(here("analysis","dkd","methylation","pmid24253112.st2.csv"), header=FALSE) %>%
  dplyr::rename(position = V3) %>%
  dplyr::rename(gene = V2) %>%
  tidyr::separate(col = position, into = c("chrom","start"), sep = c(":")) %>%
  dplyr::mutate(end = start) %>%
  dplyr::select(chrom, start, end) %>%
  dplyr::mutate(dmr = paste0("chr",chrom,"-",start,"-",end))
dmr.df$study_id <- study_id
dmr.df$phenotype <- phenotype
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_flank.gr <- Extend(dmr38.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over10.gr <- join_overlap_intersect(dar.gr, dmr38.gr)
over10_flank.gr <- join_overlap_intersect(dar.gr, dmr38_flank.gr)

# overlap with GR cut and run peaks
cut10.gr <- join_overlap_intersect(cut.gr, dmr38.gr)
cut10_flank.gr <- join_overlap_intersect(cut.gr, dmr38_flank.gr)

# compile the results
dmr.compile.df <- lapply(list(over1.gr, over2.gr, over3.gr, over4.gr, over5.gr, over6.gr, over7.gr, over8.gr, over9.gr, over10.gr), function(gr) {
  tmp <- as.data.frame(gr)
  }) %>% bind_rows() %>% arrange(seqnames, start)

cut.compile.df <- lapply(list(cut1.gr, cut2.gr, cut3.gr, cut4.gr, cut5.gr, cut6.gr, cut7.gr, cut8.gr, cut9.gr, cut10.gr), function(gr) {
  tmp <- as.data.frame(gr)
  }) %>% bind_rows() %>% arrange(seqnames, start)

# compile the flank results
dmr.compile.flank.df <- lapply(list(over1_flank.gr, over2_flank.gr, over3_flank.gr, over4_flank.gr, over5_flank.gr, over6_flank.gr, over7_flank.gr,
                                    over8_flank.gr, over9_flank.gr, over10_flank.gr), function(gr) {
  tmp <- as.data.frame(gr)
  }) %>% bind_rows() %>% arrange(seqnames, start)

cut.compile.flank.df <- lapply(list(cut1_flank.gr, cut2_flank.gr, cut3_flank.gr, cut4_flank.gr, cut5_flank.gr, cut6_flank.gr, cut7_flank.gr,
                                    cut8_flank.gr, cut9_flank.gr, cut10_flank.gr), function(gr) {
  tmp <- as.data.frame(gr)
  }) %>% bind_rows() %>% arrange(seqnames, start)
