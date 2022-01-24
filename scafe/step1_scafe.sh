#!/bin/bash
# https://github.com/chung-lab/scafe

# docker pull p4rkerw/scafe:latest

# SCRATCH1=/mnt/g/scratch
# docker run -it \
# -v /mnt/g/diabneph:$HOME/project \
# -v /mnt/g/reference:$HOME/reference \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/scafe:latest

# # export docker volumes
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"

# # to run interactive
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/scafe:latest)' /bin/bash

# library_ids=(Control_4 Control_5 DN_4 DN_5)
# for library_id in ${library_ids[*]}
# do
# bsub -G compute-parkerw -J "scafe_${library_id}" -o "$SCRATCH1/log.scafe_${library_id}" -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/scafe:latest)' \
# bash $SCRATCH1/dkd/scafe/scafe.sh $library_id
# done

export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/SCAFE/scripts:$PATH

library_id=$1
threads=10

# download the hg38 ref provided by scafe
# scafe.download.resources.genome --genome=hg38.gencode_v32

# grab cell barcodes
barcodes=$HOME/project/analysis/dkd/barcodes/rna_barcodes.csv

# filter csv by library_id
awk -F ',' 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}}{ print $(f["barcode"]),$(f["celltype"]),$(f["orig.ident"]) }' $barcodes|\
awk -v a="${library_id}" 'BEGIN { OFS="," } {if($3 == a) {print $1"-1"}}' > /tmp/barcodes.$library_id.csv

# count individual libs
rna_count_dir=/home/parkerw/project/cellranger_rna_counts/version_4.0
scafe.workflow.sc.solo \
--run_bam_path $rna_count_dir/$library_id/outs/possorted_genome_bam.bam \
--run_cellbarcode_path /tmp/barcodes.$library_id.csv \
--genome hg38.gencode_v32 \
--run_tag $library_id \
--run_outDir $SCRATCH1/scafe \
--max_thread $threads

