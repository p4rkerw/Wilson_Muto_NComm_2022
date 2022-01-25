# to run interactively on the RIS compute1 cluster
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=64GB]' -q general-interactive -a 'docker(p4rkerw/sctools:R4.1.0)' /bin/bash

# this script will annotate cicero ccan with the ChIPSeeker database to determine what genomic
# regions are linked by the connections and create circos plots

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(openxlsx)
library(circlize)
library(dplyr)
library(tibble)
library(data.table)
library(Signac)
library(stringr)
library(RColorBrewer)
library(here)
library(plyranges)

atac_aggr_prep <- here("project","analysis","dkd","atac_aggr_prep")
atacAggr <- readRDS(here(atac_aggr_prep,"step6_ccan.rds"))

# designate figures dir
figures <- here("project","analysis","dkd","figures")

# read in cell-specific dar for PCT
dar.df <- read.xlsx(here("project","analysis","dkd","markers","dar.macs2.celltype.diab_vs_ctrl.xlsx"), sheet = "PCT", rowNames = TRUE) %>%
  rownames_to_column(var = "peak") %>%
  filter(p_val < 0.05)

dar.gr <- StringToGRanges(dar.df$peak)

# extract links from aggregated object
links.gr <- Links(atacAggr)
# saveRDS(links.gr, here(figures, "links.rds"))

# split the links.gr into start and end points
links.df <- as.data.frame(links.gr)
links.start.gr <- data.frame(chrom=links.df$seqnames, start=links.df$start, end=links.df$start) %>% makeGRangesFromDataFrame()
links.end.gr <- data.frame(chrom=links.df$seqnames, start=links.df$end, end=links.df$end) %>% makeGRangesFromDataFrame()

# add back score and ccan group metadata
links.start.gr$score <- links.gr$score
links.start.gr$group <- links.gr$group
links.end.gr$score <- links.gr$score
links.end.gr$group <- links.gr$group

# extract peaks from aggregated object and overlap with links
# the links are actually the peak midpoints (and not the whole peak)
peaks.gr <- atacAggr[["peaks"]]@ranges
# saveRDS(peaks.gr, here(figures, "peaks.rds"))

# overlap peaks.gr with the ccan links
hits.start <- findOverlaps(peaks.gr, links.start.gr)
peaks.start.gr <- peaks.gr[queryHits(hits.start)]
hits.end <- findOverlaps(peaks.gr, links.end.gr)
peaks.end.gr <- peaks.gr[queryHits(hits.end)]

# add back metadata
peaks.start.gr$score <- links.start.gr$score
peaks.start.gr$group <- links.start.gr$group
peaks.end.gr$score <- links.end.gr$score
peaks.end.gr$group <- links.end.gr$group

# add a metadata column indicating whether the peak is a dar
dar.start.hits <- findOverlaps(dar.gr, peaks.start.gr)
peaks.start.gr$dar <- ifelse(seq(1:length(peaks.start.gr)) %in% subjectHits(dar.start.hits),1,0) 
dar.end.hits <- findOverlaps(dar.gr, peaks.end.gr)
peaks.end.gr$dar <- ifelse(seq(1:length(peaks.end.gr)) %in% subjectHits(dar.end.hits),1,0) 

# identify row id for start and end that contain a dar and overlap
peaks.start.id <- peaks.start.gr$dar == 1
peaks.end.id <- peaks.end.gr$dar == 1

peaks.start.dar.gr <- peaks.start.gr[peaks.start.id | peaks.end.id]
peaks.end.dar.gr <- peaks.end.gr[peaks.start.id | peaks.end.id]

# load glucocorticoid receptor binding sites
gre.gr <- fread(here("project","cut_and_run","hTERT_RPTEC","GR","GR_NT","hTERT_GR_consensus.bed")) %>%
  dplyr::rename(chrom = V1, start = V2, end = V3) %>%
  makeGRangesFromDataFrame()

# intersect start and end ccan with glucocorticoid receptor binding sites
start.gre.hits <- findOverlaps(gre.gr, peaks.start.dar.gr)
peaks.start.dar.gr$gre <- ifelse(seq(1:length(peaks.start.dar.gr)) %in% subjectHits(start.gre.hits),1,0) 

end.gre.hits <- findOverlaps(gre.gr, peaks.end.dar.gr)
peaks.end.dar.gr$gre <- ifelse(seq(1:length(peaks.end.dar.gr)) %in% subjectHits(end.gre.hits),1,0) 

# annotate the start and end coords with the peak location
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
start_anno <- annotatePeak(peaks.start.dar.gr, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)
end_anno <- annotatePeak(peaks.end.dar.gr, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)

# group utr, intron, exon, downstream into gene body
start_loc <- start_anno@anno$annotation 
start_loc <- str_split(start_loc, pattern = " ", simplify = TRUE)[,1]
start_loc <- str_replace(start_loc, pattern="3'", replacement="GeneBody")
start_loc <- str_replace(start_loc, pattern="5'", replacement="GeneBody")
start_loc <- str_replace(start_loc, pattern="Exon", replacement="GeneBody")
start_loc <- str_replace(start_loc, pattern="Downstream", replacement="GeneBody")
start_loc <- str_replace(start_loc, pattern="Intron", replacement="GeneBody")
start_loc <- str_replace(start_loc, pattern="Distal", replacement="Intergenic")

end_loc <- end_anno@anno$annotation 
end_loc <- str_split(end_loc, pattern = " ", simplify = TRUE)[,1]
end_loc <- str_replace(end_loc, pattern="3'", replacement="GeneBody")
end_loc <- str_replace(end_loc, pattern="5'", replacement="GeneBody")
end_loc <- str_replace(end_loc, pattern="Exon", replacement="GeneBody")
end_loc <- str_replace(end_loc, pattern="Downstream", replacement="GeneBody")
end_loc <- str_replace(end_loc, pattern="Intron", replacement="GeneBody")
end_loc <- str_replace(end_loc, pattern="Distal", replacement="Intergenic")

# create df with peak1 and peak2 locations
df <- cbind(start_loc, end_loc, peaks.start.dar.gr$gre, peaks.end.dar.gr$gre) %>%
  as.data.frame()
colnames(df) <- c("Peak1","Peak2","Peak1_GRE","Peak2_GRE")

# count the number of times each combo occurs
counts <- count(df, Peak1, Peak2, Peak1_GRE, Peak2_GRE) %>% 
    as.data.frame()
toplot <- counts

# generate circos plots of predicted chromatin-chromatin interactions and save to plots directory
# convert to an adjacency list with a value indicating the number of connections between
# each of the unique genomic location pairs
unique_combos <- !duplicated(t(apply(toplot, 1, sort)))
toplot <- toplot[unique_combos, ]
toplot <- dplyr::mutate(toplot, from=paste0(Peak1,"_",Peak1_GRE), to=paste0(Peak2,"_", Peak2_GRE))

# sort by GRE
#toplot <- dplyr::arrange(toplot, desc(gre))
toplot <- data.frame(from=toplot$from, to=toplot$to, n=toplot$n)

grid.col = brewer.pal(n_distinct(toplot$from), "Paired")
names(grid.col) <- unique(toplot$from)

# group all GRE together
# order the group sectors by whether or not they contain a GRE
toplot <- dplyr::mutate(toplot, order = str_split(from, pattern="_", simplify=TRUE)[,2]) %>%
  dplyr::arrange(order)

# assign grouping variable
group <- unique(unlist(toplot$from))
nm = unique(unlist(toplot$from))
names(group) <- nm

toplot$from <- factor(toplot$from, levels=group)  

pdf(here(figures,"circos_gre_dar.pdf"))
circos.clear()
par(mar=c(0.5,0.5,0.5,0.5))
circos.par(gap.after = 5)
chordDiagram(toplot,
           grid.col = grid.col,
           annotationTrack = "grid",
           preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(toplot))))))

# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter,
              CELL_META$ylim[1],
              CELL_META$sector.index,
              facing = "clockwise",
              niceFacing = TRUE,
              adj = c(0, 0.5))
 }, bg.border = NA) # here set bg.border to N
dev.off()







