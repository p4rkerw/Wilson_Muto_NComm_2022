# this script compares publicly-available differentially methylated loci to the cell-specific DAR for DKD vs. control samples
library(tidyr)
library(dplyr)
library(openxlsx)
library(GenomicRanges)
library(plyranges)
library(rtracklayer)


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

dar.df <- tidyr::separate(dar.df, col = peak, into = c("chrom","start","end"), sep = "-")
dar.gr <- makeGRangesFromDataFrame(dar.df, keep.extra.columns=TRUE)

# Kidney cytosine methylation changes improve renal function decline estimation in patients with diabetic kidney disease
# hg19
# PMID: 31165727
# Supplemental Table 1
# Description: Top 65 probes significantly associated with interstitial fibrosis and passed independent replication. (a) Association determined by linear regression models adjusted for age, gender, race, diabetes, hypertension, batch, bisulfite conversion, and degree of lymphocytic infiltrate on histology.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6549146/bin/41467_2019_10378_MOESM4_ESM.xlsx

# read in table and convert col to hg19 coordinates
dmr.df <- read.xlsx(here("analysis","dkd","methylation","41467_2019_10378_MOESM4_ESM.xlsx")) %>%
  dplyr::mutate(peak = paste0(Chr,"-",Position,"-",Position)) %>%
  dplyr::select(peak, Methylation.probe)
dmr.df <- tidyr::separate(dmr.df, col = peak, into = c("chrom","start","end"), sep = "-")
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()

# overlap with cell-specific DAR
over1.gr <- join_overlap_intersect(dar.gr, dmr38.gr)

# Kidney cytosine methylation changes improve renal function decline estimation in patients with diabetic kidney disease
# hg19
# PMID: 31165727
# File Name: Supplementary Data 3
# Description: Top 471 probes that improve model of kidney function declineusing weighted regression. (a) Model is a weighted linear regression model of adjusted eGFR slope (weight = inverse variance of adjusted eGFR slope). Base model includes variables: baseline eGFR, Diabetes, and Age (base model AIC = 206). When methylation level of probe is added to base model, the following variables are also added: methylation batch, and bisulfite conversion efficiency. (b) Association determined by linear regression models adjusted for age, gender, race, diabetes, hypertension, batch, bisulfite conversion, and degree of lymphocytic infiltrate on histology.
# read in table and convert col to hg19 coordinates
dmr.df <- read.xlsx(here("analysis","dkd","methylation","41467_2019_10378_MOESM6_ESM.xlsx")) %>%
  dplyr::mutate(peak = paste0(chr,"-",Position,"-",Position)) %>%
  dplyr::select(peak)
dmr.df <- tidyr::separate(dmr.df, col = peak, into = c("chrom","start","end"), sep = "-")
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()

# overlap with cell-specific DAR
over2.gr <- join_overlap_intersect(dar.gr, dmr38.gr)

# Assessment of differentially methylated loci in individuals with end-stage kidney disease attributed to diabetic kidney disease: an exploratory study
# PMID: 33933144
# hg19
# Supplementary Table ST3: Top-ranked differentially methylated CpG sites for matched individuals with T1DM-ESKD (n=107) vs. T1DM (n=107): FDR p≤x10-8 (Analysis 1)
dmr.df <- read.xlsx(here("analysis","dkd","methylation","13148_2021_1081_MOESM1_ESM.xlsx"), sheet = "ST3", startRow = 2) %>%
    dplyr::mutate(peak = paste0(CHR,"-",MAPINFO,"-",MAPINFO)) %>%
    dplyr::select(peak, CpG.ID, UCSC_RefGene_Name)
dmr.df <- tidyr::separate(dmr.df, col = peak, into = c("chrom","start","end"), sep = "-") %>%
  dplyr::filter(chrom != "NA")
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()

# overlap with cell-specific DAR
over3.gr <- join_overlap_intersect(dar.gr, dmr38.gr)

# Assessment of differentially methylated loci in individuals with end-stage kidney disease attributed to diabetic kidney disease: an exploratory study
# PMID: 33933144
# hg19
# Supplementary Table ST4:Top-ranked differentially methylated CpG sites for matched individuals with T1DM-ESKD (n=107) vs. T1DM (n=107): FDR p≤x10-8 and FC ±2 (Analysis 1)
dmr.df <- read.xlsx(here("analysis","dkd","methylation","13148_2021_1081_MOESM1_ESM.xlsx"), sheet = "ST4", startRow = 2) %>%
    dplyr::mutate(peak = paste0(CHR,"-",MAPINFO,"-",MAPINFO)) %>%
    dplyr::select(peak, CpG.ID)
dmr.df <- tidyr::separate(dmr.df, col = peak, into = c("chrom","start","end"), sep = "-") %>%
  dplyr::filter(chrom != "NA")
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()

# overlap with cell-specific DAR
over4.gr <- join_overlap_intersect(dar.gr, dmr38.gr)

# Assessment of differentially methylated loci in individuals with end-stage kidney disease attributed to diabetic kidney disease: an exploratory study
# PMID: 33933144
# hg19
# Supplementary Table ST11: Top-ranked differentially methylated CpG sites for individuals with T1DM-ESKD - either dialysis or transplant (n=107) vs. T1DM (n=253): FDR p≤x10-8 (Analysis 2)
dmr.df <- read.xlsx(here("analysis","dkd","methylation","13148_2021_1081_MOESM1_ESM.xlsx"), sheet = "ST11", startRow = 2) %>%
    dplyr::mutate(peak = paste0(CHR,"-",MAPINFO,"-",MAPINFO)) %>%
    dplyr::select(peak, CpG.ID)
dmr.df <- tidyr::separate(dmr.df, col = peak, into = c("chrom","start","end"), sep = "-") %>%
  dplyr::filter(chrom != "NA")
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38.gr <- liftOver(dmr.gr, ch19) %>% unlist()

# overlap with cell-specific DAR
over5.gr <- join_overlap_intersect(dar.gr, dmr38.gr)

# Cytosine methylation changes in enhancer regions of core pro-fibrotic genes characterize kidney fibrosis development
# PMID: 24098934
# hg18
# Supplemental Table 2 List of differentially methylated regions with p-value methylation level, genomic locus and nearest annotated transcript
dmr.df <- read.xlsx(here("analysis","dkd","methylation","gb-2013-14-10-r108-S3.xlsx"), startRow = 3) %>%
  dplyr::select(chr, start, end)
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38.gr <- liftOver(dmr.gr, ch18) %>% unlist()

# overlap with cell-specific DAR
over6.gr <- join_overlap_intersect(dar.gr, dmr38.gr)




