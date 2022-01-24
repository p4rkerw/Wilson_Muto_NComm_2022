# figure upset plot cell-specific DAR intersection with GR cut and run
library(here)
library(openxlsx)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(data.table)
library(Signac)
library(GenomicRanges)
library(Signac)
library(stringr)

figures <- here("project","analysis","dkd","figures")

file <- here("project","analysis","dkd","markers","dar.macs2.celltype.markers.xlsx")
idents <- getSheetNames(file)

dar.df <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
    dplyr::filter(abs_log2FC > 0.1) %>%
    rownames_to_column(var = "peak")
  if(nrow(df) > 0) {
  df$celltype <- ident
  }
  return(df)
}) %>% bind_rows()

dar.bed <- str_split(unlist(dar.df$peak), pattern="-", simplify=TRUE) %>% as.data.frame()
colnames(dar.bed) <- c("chrom","start","end")
dar.bed$value <- dar.df$avg_log2FC
dar.bed$start <- as.numeric(dar.bed$start)
dar.bed$end <- as.numeric(dar.bed$end)

dar.bed <- dplyr::arrange(dar.bed, chrom, start, end)
dar.up.bed <- dplyr::filter(dar.bed, value > 0)
dar.down.bed <- dplyr::filter(dar.bed, value < 0)

gre.bed <- fread(here("project","cut_and_run","bulk_kidney","kidney_GR_peaks.narrowPeak")) %>%
  select(V1, V2, V3, V5)
colnames(gre.bed) <- c("chrom","start","end","value")

# #####################
pdf(here(figures, "figure4a.pdf"))
circos.initializeWithIdeogram(plotType = "labels")
circos.genomicRainfall(dar.up.bed, cex=0.1, col=c("darkgreen"))
circos.genomicDensity(dar.up.bed, baseline=0, col = "black", bg.border="white")
circos.genomicDensity(gre.bed, baseline=0, col = "purple", bg.border="white")
dev.off()
#####################
# upset plot for cell-specific DAR and GR cut and run
library(here)
library(openxlsx)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(data.table)
library(Signac)
library(GenomicRanges)

file <- here("project","analysis","dkd","markers","dar.macs2.celltype.markers.xlsx")
idents <- getSheetNames(file)

dar.ls <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(avg_log2FC > 0.1) %>%
    rownames_to_column(var = "peak")
  return(df$peak)
})
names(dar.ls) <- idents

levels <- c("PCT","PST","PT_VCAM1","PT_CD36","PEC",
            "ATL","TAL1","TAL2","DCT1","DCT2",
            "PC","ICA","ICB","PODO","ENDO",
            "FIB_VSMC_MC","BCELL","TCELL","MONO")

# reorder the levels
dar.ls <- dar.ls[levels]

# create a list of all dar and convert to gr
dar.gr <- unlist(dar.ls) %>% StringToGRanges() %>% sort()
dar_peaks <- unique(unlist(dar.ls)) %>% sort()

# read in kidney bulk GR cut and run sites
gre.df <- fread(here("project","cut_and_run","bulk_kidney","kidney_GR_peaks.narrowPeak")) %>%
  mutate(peak = paste0(V1,"-",V2,"-",V3))
gre.gr <-  StringToGRanges(gre.df$peak, sep=c("-","-"))

# intersect gre.gr with dar.gr to identify gr binding sites in atac peaks
# reassign overlapping gr peaks as the dar peaks they intersect with
# this will allow quantification of number of dar peaks with a gr binding site
# append with gr binding sites that do not overlap an atac peak
hits <- findOverlaps(query=gre.gr, subject=dar.gr)
gre_hits.gr <- dar.gr[subjectHits(hits)] # atac peaks with a GR binding site
gre_nohits.gr <- gre.gr[-queryHits(hits)] # GR binding sites that do not overlap an ATAC peak
gre_hitpeaks.df <- data.frame(seqnames=seqnames(gre_hits.gr), ranges=ranges(gre_hits.gr)) %>%
  mutate(peak = paste0(seqnames,"-",ranges.start,"-",ranges.end))
gre_nohitpeaks.df <- data.frame(seqnames=seqnames(gre_nohits.gr), ranges=ranges(gre_nohits.gr)) %>%
  mutate(peak = paste0(seqnames,"-",ranges.start,"-",ranges.end))

# combine the list of gr peaks that do and do not intersect an atac peak
gre_allhits.df <- bind_rows(gre_hitpeaks.df, gre_nohitpeaks.df)

# create a GR peaks
gre_peaks.ls <- list(gre_allhits.df$peak)

# filter the dar to only include peaks with a GR binding site
dar_with_gre.ls <- lapply(seq_along(dar.ls), function(index) {
  peaks <- dar.ls[[index]]
  peaks <- peaks[peaks %in% gre_allhits.df$peak]
  return(peaks)
})
names(dar_with_gre.ls) <- names(dar.ls)

ls <- c(GR=gre_peaks.ls, dar_with_gre.ls)

mat <- list_to_matrix(ls)

# cell-specific dar shared across cell types
m1 <- make_comb_mat(mat)
# filter by intersection size and set size
# identify combinations with only a single ident by summing all the occurrences of 1 in the string
df <- strsplit(names(comb_size(m1)), split="") %>%
  as.data.frame() 
df <- sapply(df, as.numeric)
colsums <- colSums(df)

# only plot intersections with size greater than 5 unless there's only 1 peak in the group
m1 <- m1[comb_size(m1) > 10 | colsums == 1]

ht <- draw(UpSet(m1, set_order=names(ls)))
od <- column_order(ht)
cs <- comb_size(m1)
pdf(here(figures, "figure4b.pdf"))
p3 <- draw(UpSet(m1, set_order=names(ls)))
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
}) # export to pdf 6x6in 
dev.off()
