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

# bsub -G compute-parkerw -J "scafe_pool" -o "$SCRATCH1/log.scafe_pool" -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/scafe:latest)' \
# bash $SCRATCH1/Wilson_Muto_NComm_2022/scafe/scafe_pool.sh

export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/SCAFE/scripts:$PATH

# pool individual libs
scafe.workflow.sc.pool \
--lib_list_path $SCRATCH1/dkd/scafe/lib_list_path.txt \
--genome hg38.gencode_v32 \
--run_tag pool \
--run_outDir $SCRATCH1/scafe_pool \
--max_thread $threads

