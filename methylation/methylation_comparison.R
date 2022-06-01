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
library(EnsDb.Hsapiens.v86)
library(stringr)


# import liftover chains
ch19 = import.chain("G:/diabneph/analysis/dkd/methylation/hg19ToHg38.over.chain")
ch18 = import.chain("G:/diabneph/analysis/dkd/methylation/hg18ToHg38.over.chain")

# read in DAR
file <- here("analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx")
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
bulk.gr <- read.xlsx(here("analysis","dkd","cut_and_run","ST11_kidney_GR_peaks.xlsx"), colNames=FALSE) %>%
  dplyr::rename(chrom = X1, start = X2, end = X3) %>%
  dplyr::mutate(gr_peak = paste0(chrom,"-",start,"-",end)) %>%
  dplyr::select(chrom, start, end, gr_peak) %>%
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)

bulk.gr <- keepStandardChromosomes(bulk.gr, pruning.mode = 'coarse')

# annotate with nearest gene
gene.ranges <- genes(EnsDb.Hsapiens.v86)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(gene.ranges),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(gene.ranges) <- ucsc.levels
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

nearest_feature <- distanceToNearest(bulk.gr, subject = gene.ranges)
feature_hits <- gene.ranges[subjectHits(x = nearest_feature)]
df <- as.data.frame(x = mcols(x = feature_hits))
bulk.gr$gene <- df$symbol

# read in GR CUT&RUN peaks for DEX stimulated rptec
rptec.gr <- read.xlsx(here("analysis","dkd","cut_and_run","ST15_hTERT_GR_consensus.xlsx"), colNames=FALSE) %>%
  dplyr::rename(chrom = X1, start = X2, end = X3) %>%
  dplyr::mutate(gr_peak = paste0(chrom,"-",start,"-",end)) %>%
  dplyr::select(chrom, start, end, gr_peak) %>%
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)

rptec.gr <- keepStandardChromosomes(rptec.gr, pruning.mode = 'coarse')

# annotate with nearest gene
nearest_feature <- distanceToNearest(rptec.gr, subject = gene.ranges)
feature_hits <- gene.ranges[subjectHits(x = nearest_feature)]
df <- as.data.frame(x = mcols(x = feature_hits))
rptec.gr$gene <- df$symbol

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
dmr38_1.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_1_flank.gr <- Extend(dmr38_1.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over1.gr <- join_overlap_intersect(dar.gr, dmr38_1.gr)
over1_flank.gr <- join_overlap_intersect(dar.gr, dmr38_1_flank.gr)

# overlap with GR cut and run peaks
bulk1.gr <- join_overlap_intersect(bulk.gr, dmr38_1.gr)
bulk1_flank.gr <- join_overlap_intersect(bulk.gr, dmr38_1_flank.gr)
rptec1.gr <- join_overlap_intersect(rptec.gr, dmr38_1.gr)
rptec1_flank.gr <- join_overlap_intersect(rptec.gr, dmr38_1_flank.gr)

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
dmr38_2.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_2_flank.gr <- Extend(dmr38_2.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over2.gr <- join_overlap_intersect(dar.gr, dmr38_2.gr)
over2_flank.gr <- join_overlap_intersect(dar.gr, dmr38_2_flank.gr)

# overlap with GR cut and run peaks
bulk2.gr <- join_overlap_intersect(bulk.gr, dmr38_2.gr)
bulk2_flank.gr <- join_overlap_intersect(bulk.gr, dmr38_2_flank.gr)
rptec2.gr <- join_overlap_intersect(rptec.gr, dmr38_2.gr)
rptec2_flank.gr <- join_overlap_intersect(rptec.gr, dmr38_2_flank.gr)

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
dmr38_3.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_3_flank.gr <- Extend(dmr38_3.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over3.gr <- join_overlap_intersect(dar.gr, dmr38_3.gr)
over3_flank.gr <- join_overlap_intersect(dar.gr, dmr38_3_flank.gr)

# overlap with GR cut and run peaks
bulk3.gr <- join_overlap_intersect(bulk.gr, dmr38_3.gr)
bulk3_flank.gr <- join_overlap_intersect(bulk.gr, dmr38_3_flank.gr)
rptec3.gr <- join_overlap_intersect(rptec.gr, dmr38_3.gr)
rptec3_flank.gr <- join_overlap_intersect(rptec.gr, dmr38_3_flank.gr)

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
dmr38_4.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_4_flank.gr <- Extend(dmr38_4.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over4.gr <- join_overlap_intersect(dar.gr, dmr38_4.gr)
over4_flank.gr <- join_overlap_intersect(dar.gr, dmr38_4_flank.gr)

# overlap with GR cut and run peaks
bulk4.gr <- join_overlap_intersect(bulk.gr, dmr38_4.gr)
bulk4_flank.gr <- join_overlap_intersect(bulk.gr, dmr38_4_flank.gr)
rptec4.gr <- join_overlap_intersect(rptec.gr, dmr38_4.gr)
rptec4_flank.gr <- join_overlap_intersect(rptec.gr, dmr38_4_flank.gr)

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
dmr38_5.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_5_flank.gr <- Extend(dmr38_5.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over5.gr <- join_overlap_intersect(dar.gr, dmr38_5.gr)
over5_flank.gr <- join_overlap_intersect(dar.gr, dmr38_5_flank.gr)

# overlap with GR cut and run peaks
bulk5.gr <- join_overlap_intersect(bulk.gr, dmr38_5.gr)
bulk5_flank.gr <- join_overlap_intersect(bulk.gr, dmr38_5_flank.gr)
rptec5.gr <- join_overlap_intersect(rptec.gr, dmr38_5.gr)
rptec5_flank.gr <- join_overlap_intersect(rptec.gr, dmr38_5_flank.gr)

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
dmr38_6.gr <- liftOver(dmr.gr, ch18) %>% unlist()
dmr38_6_flank.gr <- Extend(dmr38_6.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over6.gr <- join_overlap_intersect(dar.gr, dmr38_6.gr)
over6_flank.gr <- join_overlap_intersect(dar.gr, dmr38_6_flank.gr)

# overlap with GR cut and run peaks
bulk6.gr <- join_overlap_intersect(bulk.gr, dmr38_6.gr)
bulk6_flank.gr <- join_overlap_intersect(bulk.gr, dmr38_6_flank.gr)
rptec6.gr <- join_overlap_intersect(rptec.gr, dmr38_6.gr)
rptec6_flank.gr <- join_overlap_intersect(rptec.gr, dmr38_6_flank.gr)

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
dmr38_7.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_7_flank.gr <- Extend(dmr38_7.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over7.gr <- join_overlap_intersect(dar.gr, dmr38_7.gr)
over7_flank.gr <- join_overlap_intersect(dar.gr, dmr38_7_flank.gr)

# overlap with GR cut and run peaks
bulk7.gr <- join_overlap_intersect(bulk.gr, dmr38_7.gr)
bulk7_flank.gr <- join_overlap_intersect(bulk.gr, dmr38_7_flank.gr)
rptec7.gr <- join_overlap_intersect(rptec.gr, dmr38_7.gr)
rptec7_flank.gr <- join_overlap_intersect(rptec.gr, dmr38_7_flank.gr)

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
dmr38_8.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_8_flank.gr <- Extend(dmr38_8.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over8.gr <- join_overlap_intersect(dar.gr, dmr38_8.gr)
over8_flank.gr <- join_overlap_intersect(dar.gr, dmr38_8_flank.gr)

# overlap with GR cut and run peaks
bulk8.gr <- join_overlap_intersect(bulk.gr, dmr38_8.gr)
bulk8_flank.gr <- join_overlap_intersect(bulk.gr, dmr38_8_flank.gr)
rptec8.gr <- join_overlap_intersect(rptec.gr, dmr38_8.gr)
rptec8_flank.gr <- join_overlap_intersect(rptec.gr, dmr38_8_flank.gr)

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
dmr38_9.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_9_flank.gr <- Extend(dmr38_9.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over9.gr <- join_overlap_intersect(dar.gr, dmr38_9.gr)
over9_flank.gr <- join_overlap_intersect(dar.gr, dmr38_9_flank.gr)

# overlap with GR cut and run peaks
bulk9.gr <- join_overlap_intersect(bulk.gr, dmr38_9.gr)
bulk9_flank.gr <- join_overlap_intersect(bulk.gr, dmr38_9_flank.gr)
rptec9.gr <- join_overlap_intersect(rptec.gr, dmr38_9.gr)
rptec9_flank.gr <- join_overlap_intersect(rptec.gr, dmr38_9_flank.gr)

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
  dplyr::mutate(dmr = paste0(chrom,"-",start,"-",end))
dmr.df$study_id <- study_id
dmr.df$phenotype <- phenotype
dmr.gr <- makeGRangesFromDataFrame(dmr.df, keep.extra.columns=TRUE)
seqlevelsStyle(dmr.gr) <- "UCSC"
dmr38_10.gr <- liftOver(dmr.gr, ch19) %>% unlist()
dmr38_10_flank.gr <- Extend(dmr38_10.gr, upstream=1000, downstream=1000)

# overlap with cell-specific DAR
over10.gr <- join_overlap_intersect(dar.gr, dmr38_10.gr)
over10_flank.gr <- join_overlap_intersect(dar.gr, dmr38_10_flank.gr)

# overlap with GR cut and run peaks
bulk10.gr <- join_overlap_intersect(bulk.gr, dmr38_10.gr)
bulk10_flank.gr <- join_overlap_intersect(bulk.gr, dmr38_10_flank.gr)
rptec10.gr <- join_overlap_intersect(rptec.gr, dmr38_10.gr)
rptec10_flank.gr <- join_overlap_intersect(rptec.gr, dmr38_10_flank.gr)

# compile the results
dmr.compile.df <- lapply(list(over1.gr, over2.gr, over3.gr, over4.gr, over5.gr, over6.gr, over7.gr, over8.gr, over9.gr, over10.gr), function(gr) {
  tmp <- as.data.frame(gr)
  }) %>% bind_rows() %>% arrange(seqnames, start)

bulk.compile.df <- lapply(list(bulk1.gr, bulk2.gr, bulk3.gr, bulk4.gr, bulk5.gr, bulk6.gr, bulk7.gr, bulk8.gr, bulk9.gr, bulk10.gr), function(gr) {
  tmp <- as.data.frame(gr)
  }) %>% bind_rows() %>% arrange(seqnames, start)

rptec.compile.df <- lapply(list(rptec1.gr, rptec2.gr, rptec3.gr, rptec4.gr, rptec5.gr, rptec6.gr, rptec7.gr, rptec8.gr, rptec9.gr, rptec10.gr), function(gr) {
  tmp <- as.data.frame(gr)
  }) %>% bind_rows() %>% arrange(seqnames, start)

# compile all dmr
dmr.compile.df <- lapply(list(dmr38_1.gr, dmr38_2.gr, dmr38_3.gr, dmr38_4.gr, dmr38_5.gr, dmr38_6.gr, dmr38_7.gr,
                                    dmr38_8.gr, dmr38_9.gr, dmr38_10.gr), function(gr) {
  tmp <- as.data.frame(gr)
  }) %>% bind_rows() %>% arrange(seqnames, start)

# compile the flank results
dar.compile.flank.df <- lapply(list(over1_flank.gr, over2_flank.gr, over3_flank.gr, over4_flank.gr, over5_flank.gr, over6_flank.gr, over7_flank.gr,
                                    over8_flank.gr, over9_flank.gr, over10_flank.gr), function(gr) {
  tmp <- as.data.frame(gr)
  }) %>% bind_rows() %>% arrange(seqnames, start)

bulk.compile.flank.df <- lapply(list(bulk1_flank.gr, bulk2_flank.gr, bulk3_flank.gr, bulk4_flank.gr, bulk5_flank.gr, bulk6_flank.gr, bulk7_flank.gr,
                                    bulk8_flank.gr, bulk9_flank.gr, bulk10_flank.gr), function(gr) {
  tmp <- as.data.frame(gr)
  }) %>% bind_rows() %>% arrange(seqnames, start)

rptec.compile.flank.df <- lapply(list(rptec1_flank.gr, rptec2_flank.gr, rptec3_flank.gr, rptec4_flank.gr, rptec5_flank.gr, rptec6_flank.gr, rptec7_flank.gr,
                                    rptec8_flank.gr, rptec9_flank.gr, rptec10_flank.gr), function(gr) {
  tmp <- as.data.frame(gr)
  }) %>% bind_rows() %>% arrange(seqnames, start)
