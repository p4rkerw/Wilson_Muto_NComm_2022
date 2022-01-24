#!/bin/bash
# this script will compute annotations and estimate ld scores for cell-specific ATAC peaks
# and cell-specific DAR in diabetes

# export docker volumes
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"

# # to run interactive
# bsub -Is -G compute-parkerw -R 'rusage[mem=64GB]' -q general-interactive -a 'docker(p4rkerw/ldsc:1.0)' /bin/bash

# # run interactive local
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
# export docker volumes
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"
# seurat_celltypes=(merge ATL DCT1 DCT2 ENDO ICA ICB BCELL TCELL FIB_VSMC_MC MONO PC PCT PEC PODO PST PT_VCAM1 PT_CD36 TAL1 TAL2)
# for seurat_celltype in ${seurat_celltypes[@]}; do
# bsub -G compute-parkerw -J "ldsc_${seurat_celltype}" -o "$SCRATCH1/log.ldsc_${seurat_celltype}" -R 'rusage[mem=64GB]' -q general -a 'docker(p4rkerw/ldsc:latest)' \
# bash $SCRATCH1/dkd/ldsc/step2a_anno_score.sh $seurat_celltype
# done

############################
# liftover from hg38 to hg19
function lift_over {
	output_bed=$(basename $1)
	cut -f1-3 $1 > /tmp/input.bed
	liftOver /tmp/input.bed /opt/hg38ToHg19.over.chain.gz /tmp/bed19/$output_bed /tmp/unlifted/$output_bed
}
export -f lift_over

# make cell-specific annotation files for autosomal chromosomes
function run_ldsc_annoscore() {
	output_dir=$1
	output_prefix=$2
	bed=$3
	chrom=$4
	awk -v a="chr${chrom}" '$1 == a' $bed > /tmp/$output_prefix.$chrom.bed	

	# make annotations using 100kb window
	python /opt/ldsc/make_annot.py \
	--bed-file /tmp/$output_prefix.$chrom.bed \
	--bimfile /opt/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chrom.bim \
	--windowsize 100000 \
	--annot-file $output_dir/$output_prefix.$chrom.annot.gz

	# generate ld scores
	python /opt/ldsc/ldsc.py \
	--l2 \
	--bfile /opt/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chrom \
	--ld-wind-cm 1 \
	--annot $output_dir/$output_prefix.$chrom.annot.gz \
	--thin-annot \
	--out $output_dir/$output_prefix.$chrom \
	--print-snps $SCRATCH1/list.txt
}
export -f run_ldsc_annoscore
###########################
# assign positional arg to cell type
seurat_celltype=$1

# update path and activate ldsc conda environ
export PATH=$PATH:/opt/conda/envs/ldsc/bin:/opt/conda/condabin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
source activate ldsc

# note that the baseline model used the below snp list to generate annotations (and not the hapmap3 list)
wget -N -P $SCRATCH1 https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/list.txt

# create temporary directories
mkdir /tmp/bed19 2> /dev/null
mkdir /tmp/unlifted 2> /dev/null

# chromosome array for iteration
chroms=($(seq 1 22))

# process celltype markers
bed_file=project/analysis/dkd/ldsc/bed/$seurat_celltype.dar.macs2.celltype.markers.bed
lift_over $bed_file
output_dir=project/analysis/dkd/ldsc/annoscore/celltype_markers
mkdir -p $output_dir
bed19_file=/tmp/bed19/$seurat_celltype.dar.macs2.celltype.markers.bed
output_prefix=$seurat_celltype
parallel -j 20 run_ldsc_annoscore $output_dir $output_prefix $bed19_file ::: ${chroms[@]}	

# process DAR files
bed_file=project/analysis/dkd/ldsc/bed/$seurat_celltype.dar.macs2.celltype.diab_vs_ctrl.bed
lift_over $bed_file
output_dir=project/analysis/dkd/ldsc/annoscore/celltype_diabetes_dar
mkdir -p $output_dir
bed19_file=/tmp/bed19/$seurat_celltype.dar.macs2.celltype.diab_vs_ctrl.bed
output_prefix=$seurat_celltype
parallel -j 20 run_ldsc_annoscore $output_dir $output_prefix $bed19_file ::: ${chroms[@]}	
