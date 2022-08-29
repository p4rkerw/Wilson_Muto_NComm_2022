#! /bin/bash

# SCRATCH1=/mnt/g/scratch
# reference=/mnt/g/reference
# docker run -it \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v /mnt/g/diabneph/cellranger_atac_counts/version_2.0:$HOME/atac_counts \
# -v /mnt/g/reference:$HOME/reference \
# -v /mnt/g/diabneph/github_repository/dkd:$HOME/github_repository \
# -v /mnt/g/diabneph/atacFastq:$HOME/atac_fastq \
# -v $reference:$HOME/reference \
# -v $SCRATCH1:$SCRATCH1 \
# -e HOME=$HOME \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/cellranger-atac:2.0 /bin/bash 

# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph/cellranger_atac_counts/version_2.0:$HOME/outs \
# $STORAGE1/atacFastq:$HOME/fastq \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/cellranger-atac:2.0)' /bin/bash 

# assign library_id to positional arg
library_id=$1 # eg. Control_1

ln -s $SCRATCH1 scratch1
cd scratch1

# clone the github repo
# git clone https://github.com/p4rkerw/Wilson_Muto_NComm_2022

# download reference if not present in $SCRATCH1
ref=refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
if [ ! -d $HOME/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 ]; then
	wget https://cf.10xgenomics.com/supp/cell-atac/$ref &&
	tar -xvzf $ref
fi

# count the fastq files
ln -s $SCRATCH1/cellranger_atac_counts/version_2.0 counts
mkdir -p $SCRATCH1/cellranger_atac_counts/version_2.0
cd counts
cellranger-atac count \
--id=$library_id \
--reference=$HOME/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--fastqs=$HOME/atac_fastq/$library_id \
--sample=$library_id \
--localcores=8 \
--localmem=128


