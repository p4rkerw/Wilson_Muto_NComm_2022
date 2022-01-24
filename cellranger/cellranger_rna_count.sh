#! /bin/bash

# to run interactive:
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph/cellranger_rna_counts/version_4.0:$HOME/outs \
# $STORAGE1/rnaFastq:$HOME/fastq \
# $STORAGE1:$STORAGE1 \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/cellranger:4.0)' /bin/bash

# assign library_id to positional arg
library_id=$1 # eg. Control_1

# change to scratch directory
ln -s $SCRATCH1 scratch1
cd scratch1

# download and build reference if not present in $SCRATCH1
if [ ! -d GRCh38-2020-A.premrna ]; then
	wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
	tar -xvzf refdata-gex-GRCh38-2020-A.tar.gz

	# create custom reference that will count intronic reads
	awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ print; $3="exon"; $9 = \
		gensub("(transcript_id\\s\"{0,1})([^;\"]+)(\"{0,1});", "\\1\\2_premrna\\3;", "g", $9); print; next}{print}' \
	      refdata-gex-GRCh38-2020-A/genes/genes.gtf > GRCh38-2020-A.premrna.gtf

	cellranger mkref \
	--genome="GRCh38-2020-A.premrna" \
	--fasta="refdata-gex-GRCh38-2020-A/fasta/genome.fa" \
	--genes="GRCh38-2020-A.premrna.gtf" \
	--nthreads=16 \
	--memgb=128
fi

# count the fastq files
ln -s $SCRATCH1/cellranger_rna_counts/version_4.0 counts
mkdir -p $SCRATCH1/cellranger_rna_counts/version_4.0
cd counts
cellranger count \
--id=$library_id \
--transcriptome=GRCh38-2020-A.premrna/ \
--fastqs=$HOME/fastq/$library_id \
--nosecondary \
--localcores=16 \
--localmem=128

# move to storage
mv -r $library_id $HOME/outs