#!/bin/bash
# this script will detect doublets in individual snATACseq 10X libraries

# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run \
# --workdir $HOME \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -v /mnt/g/diabneph:$HOME/project \
# -v /mnt/g/diabneph/cellranger_atac_counts/version_2.0:$HOME/atac_counts \
# -v /mnt/g/scratch:$HOME/output \
# -e SCRATCH1="/mnt/g/scratch" \
# -it p4rkerw/amulet:1.0 /bin/bash

# to run LSF interactive
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph/cellranger_atac_counts/version_2.0:$HOME/atac_counts \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=64GB]' -q general-interactive -a 'docker(p4rkerw/amulet:1.0)' /bin/bash

# to run LSF detached
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph/cellranger_atac_counts/version_2.0:$HOME/atac_counts \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# library_ids=(Control_1 Control_2 Control_3 Control_4 Control_5 Control_6 DN_1 DN_2 DN_3 DN_4 DN_5 DN_6 DN_7)
# for library_id in ${library_ids[*]}
# do
# bsub -Is -G compute-parkerw -R 'rusage[mem=64GB]' -q general-interactive -a 'docker(p4rkerw/amulet:1.0)' $SCRATCH1/Wilson_Muto_NComm_2022/snATAC_prep/step0_amulet.sh $library_id
# done

# update path o/w LSF environ wont see correct python version
export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH

library_id=$1
outputdir=$HOME/project/analysis/dkd/atac_aggr_prep/doublets/$library_id
mkdir -p $outputdir
bash /opt/AMULET/AMULET.sh \
atac_counts/$library_id/outs/possorted_bam.bam \
atac_counts/$library_id/outs/singlecell.csv \
/opt/AMULET/human_autosomes.txt \
/opt/AMULET/ENCFF356LFX.bed \
$outputdir \
/opt/AMULET \
--forcesorted


