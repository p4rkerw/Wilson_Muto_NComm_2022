#!/bin/bash

# useful resource https://github.com/Nealelab/UKBB_ldsc_scripts
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"

# # # to run interactive
# bsub -Is -G compute-parkerw -R 'rusage[mem=64GB]' -q general-interactive -a 'docker(p4rkerw/ldsc:1.0)' /bin/bash

# run interactive local
# SCRATCH1=/mnt/g/scratch
# docker run -it \
# -v /mnt/g/diabneph:$HOME/project \
# -v /mnt/g/reference:$HOME/reference \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# --workdir $HOME \
# p4rkerw/ldsc:latest

# to run detached
# bsub -G compute-parkerw -J "partition" -o "$SCRATCH1/log.partition" -R 'rusage[mem=32GB]' -q general -a 'docker(p4rkerw/ldsc:latest)' \
# bash $SCRATCH1/dkd/ldsc/step3_partition.sh

export PATH=$PATH:/opt/conda/envs/ldsc/bin:/opt/conda/condabin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
source activate ldsc

# download reference files
wget -P $SCRATCH1 https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz
tar -C $SCRATCH1 -xvzf $SCRATCH1/1000G_Phase3_baseline_ldscores.tgz

wget -P $SCRATCH1 https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baseline_v1.2_ldscores.tgz
tar -C $SCRATCH1 -xvzf $SCRATCH1/1000G_Phase3_baseline_v1.2_ldscores.tgz

wget -P $SCRATCH1/1000G_Phase3_baselineLD_v2.2_ldscores https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
tar -C $SCRATCH1/1000G_Phase3_baselineLD_v2.2_ldscores -xvzf $SCRATCH1/1000G_Phase3_baselineLD_v2.2_ldscores.tgz

wget -P $SCRATCH1 https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz
tar -C $SCRATCH1 -xvzf $SCRATCH1/weights_hm3_no_hla.tgz

wget -P $SCRATCH1 https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_frq.tgz
tar -C $SCRATCH1 -xvzf $SCRATCH1/1000G_Phase3_frq.tgz

# function for calculating partitioned heritability of a single annotation without controlling
# for background peaks using the v1.2 baselineLD scores
function run_ldsc_partition1.2() {
	annotation_dir=$1
	results_dir=$2
	celltype=$3
	gwas_trait=$4
	python /opt/ldsc/ldsc.py \
	--h2 project/analysis/dkd/ldsc/sumstats/${gwas_trait}_munge.sumstats.gz \
	--ref-ld-chr $annotation_dir/${celltype}.,$SCRATCH1/baseline_v1.2/baseline. \
	--w-ld-chr $SCRATCH1/weights_hm3_no_hla/weights. \
	--overlap-annot \
	--frqfile-chr $SCRATCH1/1000G_Phase3_frq/1000G.EUR.QC. \
	--print-coefficients \
	--out $results_dir/${gwas_trait}_base12_${celltype}
}
export -f run_ldsc_partition1.2

# function for calculating partitioned heritability of a single annotation without controlling
# for background peaks using the v2.2 baselineLD scores
# this model is recommended for estimating proportion of heritability
function run_ldsc_partition() {
	annotation_dir=$1
	results_dir=$2
	celltype=$3
	gwas_trait=$4
	python /opt/ldsc/ldsc.py \
	--h2 project/analysis/dkd/ldsc/sumstats/${gwas_trait}_munge.sumstats.gz \
	--ref-ld-chr $annotation_dir/${celltype}.,$SCRATCH1/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD. \
	--w-ld-chr $SCRATCH1/weights_hm3_no_hla/weights. \
	--overlap-annot \
	--frqfile-chr $SCRATCH1/1000G_Phase3_frq/1000G.EUR.QC. \
	--print-coefficients \
	--out $results_dir/${gwas_trait}_base2.2_${celltype}
}
export -f run_ldsc_partition

# analyze cell-specific peaks obtained from Signac FindMarkers function
seurat_celltypes=(ATL DCT1 DCT2 ENDO ICA ICB BCELL TCELL FIB_VSMC_MC MONO PC PCT PEC PODO PST PT_VCAM1 PT_PROM1 PT_CD36 TAL1 TAL2)
gwas_traits=(eGFR CKD MICRO URINENA)
results_dir=project/analysis/dkd/ldsc/partition/celltype_markers
annotation_dir=project/analysis/dkd/ldsc/annoscore/celltype_markers
mkdir -p $results_dir	
parallel -j 6 run_ldsc_partition $annotation_dir $results_dir ::: ${seurat_celltypes[@]} ::: ${gwas_traits[@]} 
# run_ldsc_partition $annotation_dir $results_dir $seurat_celltype $gwas_trait

# analyze differentially accessible cell-specific peaks in diabetes obtained by signac findmarkers
# seurat_celltypes=(ATL DCT1 DCT2 ENDO ICA ICB BCELL TCELL FIB_VSMC_MC MONO PC PCT PEC PODO PST PT_VCAM1 PT_PROM1 PT_CD36 TAL1 TAL2 merge)
gwas_traits=(eGFR CKD MICRO URINENA)
results_dir=project/analysis/dkd/ldsc/partition/celltype_diabetes_dar
annotation_dir=project/analysis/dkd/ldsc/annoscore/celltype_diabetes_dar
mkdir -p $results_dir	
parallel -j 6 run_ldsc_partition $annotation_dir $results_dir ::: ${seurat_celltypes[@]} ::: ${gwas_traits[@]} 

# analyze allele-specific peaks
allele_annos=(glme_dkdint_results glme_dkd_results glme_results asca)
gwas_traits=(eGFR CKD MICRO URINENA)
results_dir=project/analysis/dkd/ldsc/partition/allele_specific
annotation_dir=project/analysis/dkd/ldsc/annoscore/allele_specific
mkdir -p $results_dir	
parallel -j 6 run_ldsc_partition $annotation_dir $results_dir ::: ${allele_annos[@]} ::: ${gwas_traits[@]} 
# run_ldsc_partition $annotation_dir $results_dir ${allele_annos[@]} ${gwas_traits[@]} 

##################################################################
# cell-type-specific partitioned heritability controlling for background peaks
##################################################################
# create an ldcts file with locations of cell-specific annotations
rm $SCRATCH1/partition.ldcts
annotation_dir=project/analysis/dkd/ldsc/annoscore/celltype_markers
control_anno_dir=project/analysis/dkd/ldsc/annoscore/celltype_markers
results_dir=project/analysis/dkd/ldsc/partition/celltype_markers
mkdir -p $results_dir	

seurat_celltypes=(ATL DCT1 DCT2 ENDO ICA ICB BCELL TCELL FIB_VSMC_MC MONO PC PCT PEC PODO PST PT_VCAM1 PT_PROM1 PT_CD36 TAL1 TAL2)
for seurat_celltype in ${seurat_celltypes[*]}; do
echo -e "${seurat_celltype}\t${annotation_dir}/${seurat_celltype}.,${control_anno_dir}/allpeaks." >> $SCRATCH1/partition.ldcts
done

# function for computing annotation specific partitioned heritability while controlling
# for background peaks using the baseline 1.2 LD weights
# this model is recommended for prioritizing annotations wrt a summary statistic
# ie. which cell type or annotation is "most relevant" for a specific function or trait
function run_ldsc_cts() {
	results_dir=$1
	ldcts_file=$2
	gwas_trait=$3
	python /opt/ldsc/ldsc.py \
	--h2-cts project/analysis/dkd/ldsc/sumstats/${gwas_trait}_munge.sumstats.gz \
	--ref-ld-chr $SCRATCH1/baseline_v1.2/baseline. \
	--w-ld-chr $SCRATCH1/weights_hm3_no_hla/weights. \
	--ref-ld-chr-cts $ldcts_file \
	--overlap-annot \
	--print-coefficients \
	--frqfile-chr $SCRATCH1/1000G_Phase3_frq/1000G.EUR.QC. \
	--out $results_dir/${gwas_trait}_base1.2
}
export -f run_ldsc_cts
gwas_traits=(eGFR CKD MICRO URINENA)
parallel -j 12 run_ldsc_cts $results_dir $SCRATCH1/partition.ldcts ::: ${gwas_traits[@]} 

# cell-specific analyses for dar in diabetes
rm $SCRATCH1/dar.partition.ldcts
annotation_dir=project/analysis/dkd/ldsc/annoscore/celltype_diabetes_dar
control_anno_dir=project/analysis/dkd/ldsc/annoscore/celltype_markers
results_dir=project/analysis/dkd/ldsc/partition/celltype_diabetes_dar
mkdir -p $results_dir	

seurat_celltypes=(ATL DCT1 DCT2 ENDO ICA ICB BCELL TCELL FIB_VSMC_MC MONO PC PCT PEC PODO PST PT_VCAM1 PT_PROM1 PT_CD36 TAL1 TAL2 merge)
for seurat_celltype in ${seurat_celltypes[*]}; do
echo -e "${seurat_celltype}\t${annotation_dir}/${seurat_celltype}.,${control_anno_dir}/allpeaks." >> $SCRATCH1/dar.partition.ldcts
done
gwas_traits=(eGFR CKD MICRO URINENA)
parallel -j 12 run_ldsc_cts $results_dir $SCRATCH1/dar.partition.ldcts ::: ${gwas_traits[@]} 

# allele-specific analyses
rm $SCRATCH1/allele.partition.ldcts
annotation_dir=project/analysis/dkd/ldsc/annoscore/allele_specific
control_anno_dir=project/analysis/dkd/ldsc/annoscore/celltype_markers
results_dir=project/analysis/dkd/ldsc/partition/allele_specific
mkdir -p $results_dir	

allele_annos=(asca glme_dkdint_results glme_dkd_results glme_results)
for allele_anno in ${allele_annos[*]}; do
echo -e "${allele_anno}\t${annotation_dir}/${allele_anno}.,${control_anno_dir}/allpeaks." >> $SCRATCH1/allele.partition.ldcts
done
gwas_traits=(eGFR CKD MICRO URINENA)
parallel -j 12 run_ldsc_cts $results_dir $SCRATCH1/allele.partition.ldcts ::: ${gwas_traits[@]} 

# workflow for summarizing results and accompanying paper
# https://kevinlkx.github.io/analysis_pipelines/sldsc_example_neuron_ATACseq.html
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7773145/

