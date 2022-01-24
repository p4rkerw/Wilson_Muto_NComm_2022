#!/bin/bash
# this script will aggregate snAtacseq output using cellranger-atac 

# to run locally:
# SCRATCH1=/mnt/g/scratch
# reference=/mnt/g/reference
# docker run -it \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v /mnt/g/diabneph/cellranger_atac_counts/version_2.0:$HOME/atac_counts \
# -v /mnt/g/reference:$HOME/reference \
# -v /mnt/g/diabneph/github_repository/dkd:$HOME/github_repository \
# -v $reference:$HOME/reference \
# -v $SCRATCH1:$SCRATCH1 \
# -e HOME=$HOME \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/cellranger-atac:2.0 /bin/bash 

# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph/cellranger_atac_counts/version_2.0:$HOME/atac_counts \
# $SCRATCH1/dkd:$HOME/github_repository \
# $SCRATCH1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/cellranger-atac:2.0)' /bin/bash

# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph/cellranger_atac_counts/version_2.0:$HOME/atac_counts \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/cellranger-atac:2.0)' bash $SCRATCH1/dkd/cellranger/cellranger_atac_aggr.sh

ln -s $SCRATCH1 scratch1
cd scratch1

# clone the github repo
# git clone https://github.com/p4rkerw/dkd

# download reference if not present in $SCRATCH1
ref=refdata-cellranger-arc-GRCh38-2020-A-2.0.0
if [ ! -d $HOME/reference/$ref ]; then
	wget https://cf.10xgenomics.com/supp/cell-atac/$ref.tar.gz &&
	tar -C reference -xvzf $ref.tar.gz &&
	rm $ref.tar.gz
fi

# update file paths in aggregation csv on the fly with sed to match docker mount path
# csv file takes absolute paths 
cellranger_atac_counts_dir=$HOME/atac_counts
sed "s|cellranger_atac_counts_dir|${cellranger_atac_counts_dir}|g" $HOME/github_repository/cellranger/atac_aggr.csv > /tmp/atac_aggr.csv

# aggregate the output
cellranger-atac aggr \
--id=cellranger_atac_aggr \
--csv=/tmp/atac_aggr.csv \
--reference=$HOME/reference/$ref \
--normalize=none \
--nosecondary \
--localcores=8 \
--localmem=128


