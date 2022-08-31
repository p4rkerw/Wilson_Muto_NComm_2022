# **Multimodal single cell sequencing implicates chromatin accessibility and genetic background in diabetic kidney disease progression**
__*Parker C. Wilson, *Yoshiharu Muto, Haojia Wu, Anil Karihaloo, Sushrut S. Waikar, Benjamin D. Humphreys__  
*These authors contributed equally  

This manuscript has been provisionally accepted at Nature Communications. In the meantime,
please cite the [preprint](https://www.biorxiv.org/content/10.1101/2022.01.28.478204v1)
```
Wilson PC, Muto Y, Wu H, Karihaloo A, Waikar SS, Humphreys BD.
Multimodal single cell sequencing of human diabetic kidney disease implicates chromatin accessibility and genetic background in disease progression
bioRxiv. 2022 Jan 1
doi: https://doi.org/10.1101/2022.01.28.478204
```
The code associated with this publication has been deposited in Zenodo
[![DOI](https://zenodo.org/badge/451634843.svg)](https://zenodo.org/badge/latestdoi/451634843)



<a href="https://zenodo.org/account/settings/github/repository/p4rkerw/Wilson_Muto_NComm_2022"><img src="https://zenodo.org/badge/451634843.svg" alt="DOI"></a>

The bioRxiv preprint can be found [here](https://www.biorxiv.org/content/10.1101/2022.01.28.478204v1)

Single cell sequencing data generated for this manuscript (snRNA-seq: 1 control, 2 DKD; snATAC-seq: 1 control, 7 DKD) and cellranger v4.0 / cellranger-atac v2.0 gene and peak count matrices for all analyzed samples (snRNA-seq: 6 control, 5 DKD; snATAC-seq: 6 control, 7 DKD) can be found in GEO. </br>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195460

CUT&RUN and bulk ATAC data generated for this manuscript can be found in a second repository: </br>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195443

Sequencing data generated for previously-published kidney snATAC (5 control) and snRNA libraries (5 control, 3 DKD) can be found at the links below:
Note that there are two different sequencing configurations for the snATAC libraries. Samples 1-3 have an I1 index file, whereas samples 4-5 do not. The I1 index files are not needed to run cellranger-atac count on samples 4-5. Index sequences used for demultplexing are printed in the header lines of R1 and R2. R1 and R3 contain the paired-end reads and R2 contains the cell barcode information.  <br/>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151302 <br/>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131882 <br/>
https://data.humancellatlas.org/explore/projects/2af52a13-65cb-4973-b513-39be38f2df3f

[snATAC barcodes](https://github.com/p4rkerw/Wilson_Muto_NComm_2022/blob/main/barcodes/atac_barcodes.csv) and [snRNA barcodes](https://github.com/p4rkerw/Wilson_Muto_NComm_2022/blob/main/barcodes/rna_barcodes.csv) used for the analysis can be found in this github repository

Finalized R objects can be visualized interactively or downloaded from the cellxgene [website](https://cellxgene.cziscience.com/collections/b3e2c6e3-9b05-4da9-8f42-da38a664b45b). The snATAC object only includes the gene activity "RNA" assay and does not have a raw or normalized peak count matrix.


Welcome to our GitHub repository!  
Here you will find analysis scripts for our manuscript where we compare control kidney to diabetic kidney disease (DKD) samples using paired snRNAseq and snATACseq. Please contact the corresponding author, Dr. Benjamin Humphreys, with questions or comments.  
<br/>
Thanks,  
Parker and Yoshi
<br/><br/>
![alt text](http://humphreyslab.com/wp-content/uploads/2015/12/favicon-H.jpg)  
Visit the Humphrey's lab website:   
www.humphreyslab.com  
<br/>
Check out our interactive datasets with Kidney Interactive mulTiomics (KIT):  
http://humphreyslab.com/SingleCell/
<br/><br/>
Find us on Twitter: 
<br/>
  <a href="https://twitter.com/parkercwilson?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false"> @parkercwilson</a>
  <a href="https://twitter.com/YoshiharuMuto?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false"> @YoshiharuMuto</a>
  <a href="https://twitter.com/HumphreysLab?ref_src=twsrc%5Etfw" class="twitter-follow-button" data-show-count="false"> @HumphreysLab</a>
<br/><br/>
Find us on Docker Hub:  
[p4rkerw@dockerhub](https://hub.docker.com/search?q=p4rkerw&type=image)
<br/>

**snRNA preprocessing and quality control workflow**
1. Align and count each snRNA library (cellranger/cellranger_rna_count.sh)
Libraries were generated from a nuclear dissociation and aligned to refdata-gex-GRCh38-2020-A, which can be downloaded from the 10X genomics website: https://support.10xgenomics.com/.

2. Aggregate the six snRNA libraries using the cellranger/rna_aggr.csv file (cellranger/cellranger_rna_aggr.sh)

3. Run standard Seurat QC on the aggregated snRNA data and eliminate doublet barcodes identified by DoubletFinder (snRNA_prep/step1_prep.R)

4. Make final cell annotations and generate the final snRNA object (snRNA_prep/step2_anno.R) - The snRNA_prep/step2_anno.rds object is considered the final snRNA object for downstream analysis.

5. (optional) Correct for ambient RNA contamination (snRNA_prep/step3_soupx.R) - This step was not included in the final workflow because we had a very low ambient RNA contamination.

**snATAC preprocessing and quality control workflow**  
1. Align and count each ATAC library (cellranger/cellranger_atac_count.sh)  
Libraries were generated from a nuclear dissociation and aligned to refdata-cellranger-arc-A-2.0.0 which can be downloaded from the 10X genomics website: https://support.10xgenomics.com/. 

2. Aggregate the seven snATAC libraries using the cellranger/atac_aggr.csv file (cellranger/cellranger_atac_aggr.sh)

3. Run AMULET on individual snATAC libraries to identify doublet barcodes (snATAC_prep/step0_amulet.sh)

4. Run standard Signac QC on the aggregated snATAC data, eliminate doublet barcodes identified by AMULET, and transfer labels from the final snRNA object to predict cell types. (snATAC_prep/step1_prep.R)

5. Filter barcodes using the predicted.id thresholds from Seurat label transfer and call cell-specific peaks using the Signac MACS2 wrapper (snATAC_prep/step2_call_m2peaks.R)

6. Create a new peak-cell matrix with the MACS2 peaks using the Signac FeatureMatrix function (snATAC_prep/step3_count_m2peaks.R)

7. Perform additional filtering based on the fragments of reads in MACS2 peaks and clustering of snATAC cell types (snATAC_prep/step4_anno_m2peaks.R). The annotations in this step are considered the final barcode annotations for the downstream analysis.

8. Run chromVAR on the final snATAC object (snATAC_prep/step5_chromVAR.R)

9. Identify cis-coaccessibility networks in the final snATAC object (snATAC_prep/step6_ccan.R)



**snRNA and snATAC standard analysis workflow**

1. Find cell-specific genes and differentially expressed genes (DEG) in DKD in the snRNA dataset (analysis/find_deg.R)  

2. Find cell-specific peaks and differentially accessible chromatin (DAR) in DKD in the snATAC dataset (analysis/find_dar.R)  
find_atac_dar.R

3. Find cell-specific motifs and differential chromVAR activities in DKD in the snATAC dataset (analysis/find_motifs.R)

4. Identify peak-gene links using the Signac LinkPeaks function (analysis/find_gene_enh.R)

5. Perform footprinting analysis for NR3C1, NR3C2, REL, and HNF4A (analysis/find_footprints.R)

6. Annotate DAR using ChIPseeker and the Fantom promoter / enhancer database (analysis/dar_enh_anno.R)

7. Visualize GR interactions with cell-specific DAR using circlize (analysis/circlize_ccans_gre.R)


**Allele-specific workflow**

1. Use SALSA to genotype, phase, correct for mapping bias, and generate single cell allele-specific counts for all libraries with both snRNA and snATAC (6 control, 5 DKD) - (allele_specific/run_salsa.sh)

2. Generate an allele-specific count matrix for each predicted peak-gene combination (allele_specific/analyze_allele_counts.R)

3. Model the allele-specific count matrix using a glmer with lme4 (allele_specific/model_allele_counts.R)


**Bulk ATAC workflow**

0-7. Align and deduplicate bulk ATAC before calling MACS2 bulk ATAC peaks (bulk_atac/*)


**CUT&RUN workflow**

1. Align fastq with bowtie2 and call peaks with MACS2 (cut_and_run/*)


**Bulk RNA-seq workflow**

1. Download publicly-available bulk RNA-seq for human DKD (bulk_rnaseq/sra_download_fastq.sh)

2. Count the fastq with salmon (bulk_rnaseq/salmon_count.sh)

3. Analyze the counts with DESeq2 (bulk_rnaseq/find_bulk_degs.R)


**Partitioned heritability workflow**

1. Generate bed files using the cell-specific peaks, cell-specific DAR that change in DKD, and allele-specific effect peaks (ldsc/step0_prep_bed.R)

2. Munge publicly-available GWAS statistics into ldsc format (ldsc/step1_munge.sh)

3. Generate custom annotations and ld scores for each of the bed files (ldsc/step2a_annoscore.sh, ldsc/step2b_annoscore.sh)

4. Partition heritability using the cell-type-specific workflow and estimate proportion of heritability for each annotation (ldsc/step3_partition.sh)


**Transcribed cis-regulatory element workflow with SCAFE**

1. Generate tCRE counts for each snRNA library (scafe/step1_scafe.sh)

2. Pool the tCRE into consensus tCRE (scafe/step2_scafe_pool.sh)

3. Process the aggregated tCRE and identify differential tCRE between conditions (scafe/step3_scafe_prep4.R)


**Methylation workflow**

1. Download publicly-available differentially methylated regions and intersect with DAR (methylation/methylation_comparison.R)


**Figures**

Scripts for generating figures in the manuscript.



