library(GenomicRanges)
library(plyranges)
library(data.table)
library(openxlsx)
library(stringr)
library(dplyr)

htert_gr_sites.gr <- fread("G:/diabneph/cut_and_run/GR_sigpeaks.bed")  %>%
  rename(chrom=V1, start=V2, end=V3) %>%
  makeGRangesFromDataFrame()

dar_ptc_diabetes <- read.xlsx("G:/diabneph/analysis/combined_adv/markers/dar.celltype.diab_vs_ctrl.xlsx",
                     sheet="PCT") 
dar_ptc_diabetes <- dar_ptc_diabetes[dar_ptc_diabetes$p_val_adj < 0.05, ]

dar_ptc_diabetes.gr <- dar_ptc_diabetes[,1] %>%
  str_split(pattern="-", simplify=TRUE) %>%
  as.data.frame() %>%
  rename(chrom=V1, start=V2, end=V3) %>%
  makeGRangesFromDataFrame()


dar_ptc_diabetes_with_gr.gr <- join_overlap_intersect(htert_gr_sites.gr, dar_ptc_diabetes.gr)

htert_bulk_atac_peaks.gr <- fread("G:/diabneph/bulkATAC/hTERT_ATAC_rep1/hTERT_rep1_peaks.narrowPeak") %>%
  filter(V9 > -log10(0.05)) %>%
  select(V1, V2, V3) %>%
  rename(chrom=V1, start=V2, end=V3) %>%
  makeGRangesFromDataFrame()


pct_snrna_atac_peaks.gr <- fread("G:/diabneph/analysis/combined_adv/narrowPeak/PCT_peaks.narrowPeak") %>%
  filter(V9 > -log10(0.05)) %>%
  select(V1, V2, V3) %>%
  rename(chrom=V1, start=V2, end=V3) %>%
  makeGRangesFromDataFrame()

primary_rptec_bulk_atac_peaks.gr <- fread("G:/diabneph/bulkATAC/Prim_ATAC_rep1/Primary_rep1_peaks.narrowPeak") %>%
  filter(V9 > -log10(0.05)) %>%
  select(V1, V2, V3) %>%
  rename(chrom=V1, start=V2, end=V3) %>%
  makeGRangesFromDataFrame()

overlap_pct_htert_peaks.gr <- join_overlap_intersect(htert_bulk_atac_peaks.gr, pct_snrna_atac_peaks.gr)

overlap_pct_primary_peaks.gr <- join_overlap_intersect(primary_rptec_bulk_atac_peaks.gr, pct_snrna_atac_peaks.gr)

##################################PT_VCAM1 comparison
ptvcam1_snrna_atac_peaks.gr <- fread("G:/diabneph/analysis/combined_adv/narrowPeak/PTVCAM1_peaks.narrowPeak") %>%
  filter(V9 > -log10(0.05)) %>%
  select(V1, V2, V3) %>%
  rename(chrom=V1, start=V2, end=V3) %>%
  makeGRangesFromDataFrame()


overlap_ptvcam1_htert_peaks.gr <- join_overlap_intersect(htert_bulk_atac_peaks.gr, ptvcam1_snrna_atac_peaks.gr)

overlap_ptvcam1_primary_peaks.gr <- join_overlap_intersect(primary_rptec_bulk_atac_peaks.gr, ptvcam1_snrna_atac_peaks.gr)

dar_ptvcam_diabetes <- read.xlsx("G:/diabneph/analysis/combined_adv/markers/dar.celltype.diab_vs_ctrl.xlsx",
                              sheet="PTVCAM1")

dar_ptvcam_diabetes <- dar_ptvcam_diabetes[dar_ptvcam_diabetes$p_val_adj < 0.05, ]


dar_ptvcam_diabetes.gr <- dar_ptvcam_diabetes[,1] %>% 
  str_split(pattern="-", simplify=TRUE) %>%
  as.data.frame() %>%
  rename(chrom=V1, start=V2, end=V3) %>%
  makeGRangesFromDataFrame()


dar_ptvcam_diabetes_with_gr.gr <- join_overlap_intersect(htert_gr_sites.gr, dar_ptvcam_diabetes.gr)




