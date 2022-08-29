#!/bin/bash
# this script will compute annotations and estimate ld scores for regions identified using the allele-specific models in proximal tubule

# # export docker volumes
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
# p4rkerw/ldsc:1.0

# to run detached
# bsub -G compute-parkerw -J "ldsc_${celltype}" -o "$SCRATCH1/log.ldsc_${celltype}" -R 'rusage[mem=32GB]' -q general -a 'docker(p4rkerw/ldsc:1.0)' \
# bash $SCRATCH1/Wilson_Muto_NComm_2022/ldsc/step2b_anno_score.sh 

export PATH=$PATH:/opt/conda/envs/ldsc/bin:/opt/conda/condabin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
source activate ldsc

# liftover from hg38 to hg19
function lift_over {
	output_bed=$(basename $1)
	liftOver $1 /opt/hg38ToHg19.over.chain.gz /tmp/bed19/$output_bed /tmp/unlifted/$output_bed
}
mkdir /tmp/bed19 2> /dev/null
mkdir /tmp/unlifted 2> /dev/null
bed_dir=project/analysis/dkd/ldsc/bed/
bed_files=($bed_dir/asca.bed $bed_dir/glme_dkdint_results.bed $bed_dir/glme_dkd_results.bed $bed_dir/glme_results.bed)
export -f lift_over
parallel lift_over ::: ${bed_files[@]} 

# note that the baseline model used the below snp list to generate annotations (and not the hapmap3 list)
wget -N -P $SCRATCH1 https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/list.txt

# make cell-specific annotation files for autosomal chromosomes
function run_ldsc_annoscore() {
	output_dir=$1
	output_prefix=$2
	bed=$3
	chrom=$4
	awk -v a="chr${chrom}" '$1 == a' $bed > /tmp/$output_prefix.$chrom.bed	

	# make annotations
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


# create chromosome array for iteration
chroms=($(seq 1 22))

# create output directory
output_dir=project/analysis/dkd/ldsc/annoscore/allele_specific
mkdir -p $output_dir

# asca
bed19_file=/tmp/bed19/asca.bed
output_prefix=asca
parallel -j 20 run_ldsc_annoscore $output_dir $output_prefix $bed19_file ::: ${chroms[@]}

# glme_results
bed19_file=/tmp/bed19/glme_results.bed
output_prefix=glme_results
parallel -j 20 run_ldsc_annoscore $output_dir $output_prefix $bed19_file ::: ${chroms[@]}

# glme_dkd_results
bed19_file=/tmp/bed19/glme_dkd_results.bed
output_prefix=glme_dkd_results
parallel -j 20 run_ldsc_annoscore $output_dir $output_prefix $bed19_file ::: ${chroms[@]}

# glme_dkdint_results
bed19_file=/tmp/bed19/glme_dkdint_results.bed
output_prefix=glme_dkdint_results
parallel -j 20 run_ldsc_annoscore $output_dir $output_prefix $bed19_file ::: ${chroms[@]}

# all peaks control for all features included in final seurat object
tail -n+2 project/analysis/dkd/atac_aggr_prep/step2_peaks.gr | awk '{print $1,$2-1,$3,$4}' > /tmp/allpeaks.bed
lift_over /tmp/allpeaks.bed
bed19_file=/tmp/bed19/allpeaks.bed
output_prefix=allpeaks
output_dir=project/analysis/dkd/ldsc/annoscore/celltype_markers
parallel -j 20 run_ldsc_annoscore $output_dir $output_prefix $bed19_file ::: ${chroms[@]}	
