#! /bin/bas
# this script will aggregate cellranger count output

# run detached:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $HOME/Wilson_Muto_NComm_2022/cellranger:$HOME/github_repository \
# $STORAGE1/diabneph/cellranger_rna_counts:$HOME/rna_counts \
# $STORAGE1/diabneph/analysis/combined_adv:$HOME/outs \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -R 'rusage[mem=128GB]' -q general -o log.out -a 'docker(p4rkerw/cellranger:4.0)' bash Wilson_Muto_NComm_2022/cellranger/cellranger_rna_aggr.sh 

# to run interactive:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $HOME/Wilson_Muto_NComm_2022/cellranger:$HOME/github_repository \
# $STORAGE1/diabneph/cellranger_rna_counts/version_4.0:$HOME/rna_counts \
# $STORAGE1/diabneph/analysis/combined_adv:$HOME/outs \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/cellranger:4.0)' /bin/bash 

# the cellranger count output
ln -s $HOME/outs ~/cellranger_output
cd ~/cellranger_output

# aggregate the output
cellranger aggr \
--id=cellranger_rna_aggr \
--csv=$HOME/github_repository/crRnaAggr_combined_adv.csv \
--normalize="none" \
--nosecondary
