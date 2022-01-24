#!/bin/bash
# this script will run the salsa workflow in an LSF HPC environ

# clone salsa repo to scratch dir
git -C $SCRATCH1 clone https://github.com/p4rkerw/SALSA

# download and unpack cellranger-atac ref to scratch (if not already in $STORAGE1/reference)
wget -P $SCRATCH1 https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
tar -C $SCRATCH1 -xvzf $SCRATCH1/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
rm $SCRATCH1/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz

# libraries to analyze
library_ids=(Control_1 Control_2 Control_3 Control_4 Control_5 DN_1 DN_2 DN_3 DN_4 DN_5 Control_6)

# export docker volumes
export LSF_DOCKER_VOLUMES="$HOME:$HOME \
$STORAGE1/diabneph:$HOME/project \
$STORAGE1/reference:$HOME/reference \
$SCRATCH1:$SCRATCH1"

# to run interactive
# bsub -Is -G compute-parkerw -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/salsa:latest)' /bin/bash

# to run detached on the RIS compute1 cluster:
#######GENOTYPE ATAC SAMPLES##############
modality=atac
for library_id in ${library_ids[*]}
do
bsub -G compute-parkerw -J "gatk_${modality}_${library_id}" -o "$SCRATCH1/log.gatk_${modality}_${library_id}" -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/salsa:latest)' \
bash $SCRATCH1/SALSA/step1_gatk_genotype.sh \
--inputbam project/cellranger_atac_counts/version_2.0/$library_id/outs/possorted_bam.bam \
--library_id $library_id \
--reference $SCRATCH1/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--gatk_bundle reference/gatk \
--outputdir project/analysis/dkd/vcf/atac_genotype \
--outputvcf $library_id.$modality.vcf.gz \
--modality $modality \
--threads 8
done

#######GENOTYPE RNA SAMPLES##############
# genotype rna samples
modality=rna
for library_id in ${library_ids[*]}
do
bsub -G compute-parkerw -J "gatk_${modality}_${library_id}" -o "$SCRATCH1/log.gatk_${modality}_${library_id}" -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/salsa:latest)' \
bash $SCRATCH1/SALSA/step1_gatk_genotype.sh \
--inputbam project/cellranger_rna_counts/version_4.0/$library_id/outs/possorted_genome_bam.bam \
--library_id $library_id \
--reference reference/refdata-gex-GRCh38-2020-A \
--gatk_bundle reference/gatk \
--outputdir  project/analysis/dkd/vcf/rna_genotype \
--outputvcf $library_id.$modality.vcf.gz \
--modality $modality \
--threads 8
done

######MERGE GENOTYPES######
joint_outputdir=project/analysis/dkd/vcf/joint_genotype
for library_id in ${library_ids[*]}
do
bsub -G compute-parkerw -J "merge_${library_id}" -o "$SCRATCH1/log.merge_${library_id}" -R 'rusage[mem=32GB]' -q general -a 'docker(p4rkerw/salsa:latest)' \
bash $SCRATCH1/SALSA/step2_merge_geno.sh \
--library_id $library_id \
--vcfone  project/analysis/dkd/vcf/atac_genotype/$library_id.atac.vcf.gz \
--vcftwo  project/analysis/dkd/vcf/rna_genotype/$library_id.rna.vcf.gz \
--outputvcf $library_id.pass.joint.vcf.gz \
--outputdir $joint_outputdir \
--threads 10
done

# ensure that the downloaded vcf have the same contig style by adding chr to the chromosome names
# for SNV and INDEL
# export PATH=/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH
# mkdir -p reference/phasing/biallelic_SNV_and_INDEL/ucsc
# for i in {1..22} X;do echo "${i} chr${i}";done > /tmp/rename_chrm.txt
# for VCF in $(ls reference/phasing/biallelic_SNV_and_INDEL/*.vcf.gz); do
# echo $VCF
# bcftools annotate $VCF --threads $threads --rename-chrs /tmp/rename_chrm.txt -Oz -o reference/phasing/biallelic_SNV_and_INDEL/ucsc/$(basename $VCF)
# bcftools index --threads $threads --tbi reference/phasing/biallelic_SNV_and_INDEL/ucsc/$(basename $VCF)
# done

#######PHASE VCF#########
for library_id in ${library_ids[*]}
do
bsub -G compute-parkerw -J "phase_${library_id}" -o "$SCRATCH1/log.phase_${library_id}" -R 'rusage[mem=32GB]' -q general -a 'docker(p4rkerw/salsa:latest)' \
bash $SCRATCH1/SALSA/step3_phase_vcf.sh \
--library_id $library_id \
--inputvcf project/analysis/dkd/vcf/joint_genotype/$library_id.pass.joint.vcf.gz \
--outputdir project/analysis/dkd/vcf/phasing \
--outputvcf $library_id.pass.joint.hcphase.vcf.gz \
--phasingref reference/phasing/biallelic_SNV_and_INDEL/ucsc \
--hcphase \
--snvindel \
--verbose \
--reproduce 
done

#######ANNOTATE VCF#########
for library_id in ${library_ids[*]}
do
bsub -G compute-parkerw -J "anno_${library_id}" -o "$SCRATCH1/log.anno_${library_id}" -R 'rusage[mem=32GB]' -q general -a 'docker(p4rkerw/salsa:latest)' \
bash $SCRATCH1/SALSA/step4_gatk_anno_vcf.sh \
--library_id $library_id \
--inputvcf project/analysis/dkd/vcf/phasing/$library_id.pass.joint.hcphase.vcf.gz \
--outputdir project/analysis/dkd/vcf/funcotation \
--outputvcf $library_id.pass.joint.hcphase.funco.vcf.gz \
--reference $SCRATCH1/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--funcotation reference/funcotator_dataSources.v1.6.20190124g \
--output_table $library_id.pass.joint.hcphase.formatted.csv \
--modality atac \
--threads 20
done

#######FILTER RNA BAM WITH BARCODES######
for library_id in ${library_ids[*]}
do
bsub -G compute-parkerw -J "bcfr_${library_id}" -o "$SCRATCH1/log.bcfr_${library_id}" -R 'rusage[mem=32GB]' -q general -a 'docker(p4rkerw/salsa:latest)' \
bash $SCRATCH1/SALSA/step5_filterbam.sh \
--library_id $library_id \
--validate \
--inputbam project/cellranger_rna_counts/version_4.0/$library_id/outs/possorted_genome_bam.bam \
--modality rna \
--barcodes project/analysis/dkd/barcodes/rna_barcodes.csv \
--outputdir project/analysis/dkd/wasp_rna \
--outputbam $library_id.bcfilter.bam \
--threads 20
done

#######FILTER ATAC BAM WITH BARCODES######
wasp_atac_dir=project/analysis/dkd/wasp_atac
mkdir $wasp_atac_dir
for library_id in ${library_ids[*]}
do
bsub -G compute-parkerw -J "bcfa_${library_id}" -o "$SCRATCH1/log.bcfa_${library_id}" -R 'rusage[mem=32GB]' -q general -a 'docker(p4rkerw/salsa:latest)' \
bash $SCRATCH1/SALSA/step5_filterbam.sh \
--library_id $library_id \
--validate \
--inputbam project/cellranger_atac_counts/version_2.0/$library_id/outs/possorted_bam.bam \
--modality atac \
--barcodes project/analysis/dkd/barcodes/atac_barcodes.csv \
--outputdir project/analysis/dkd/wasp_atac \
--outputbam $library_id.bcfilter.bam \
--threads 20
done

######WASP WITH RNA BAM########
for library_id in ${library_ids[*]}
do
bsub -G compute-parkerw -J "waspr_${library_id}" -o "$SCRATCH1/log.waspr_${library_id}" -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/salsa:latest)' \
bash $SCRATCH1/SALSA/step6_wasp.sh \
--inputvcf project/analysis/dkd/vcf/funcotation/$library_id.pass.joint.hcphase.funco.vcf.gz  \
--inputbam project/analysis/dkd/wasp_rna/$library_id.bcfilter.bam \
--outputdir project/analysis/dkd/wasp_rna \
--outputbam $library_id.hcphase.wasp.bam \
--genotype rna \
--stargenome reference/refdata-gex-GRCh38-2020-A/star \
--library_id $library_id \
--modality rna \
--isphased \
--threads 20
done

######WASP WITH ATAC BAM########
for library_id in ${library_ids[*]}
do
rm $SCRATCH1/log.waspa_${library_id}
bsub -G compute-parkerw -J "waspa_${library_id}" -o "$SCRATCH1/log.waspa_${library_id}" -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/salsa:latest)' \
bash $SCRATCH1/SALSA/step6_wasp.sh \
--inputvcf project/analysis/dkd/vcf/funcotation/$library_id.pass.joint.hcphase.funco.vcf.gz \
--inputbam project/analysis/dkd/wasp_atac/$library_id.bcfilter.bam \
--outputdir project/analysis/dkd/wasp_atac \
--outputbam $library_id.hcphase.wasp.bam \
--atacref $SCRATCH1/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--genotype atac \
--library_id $library_id \
--modality atac \
--isphased \
--threads 20
done

######COUNT WASP RNA BAM#########
# mkdir project/analysis/dkd/wasp_rna/counts
for library_id in ${library_ids[*]}
do
bsub -G compute-parkerw -J "rcnt_${library_id}" -o "$SCRATCH1/log.rcnt_${library_id}" -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/salsa:latest)' \
bash $SCRATCH1/SALSA/step7_gatk_alleleCount.sh \
--inputvcf project/analysis/dkd/vcf/funcotation/$library_id.pass.joint.hcphase.funco.vcf.gz \
--inputbam project/analysis/dkd/wasp_rna/$library_id.hcphase.wasp.bam \
--outputdir project/analysis/dkd/wasp_rna/counts \
--barcodes project/analysis/dkd/barcodes/rna_barcodes.csv \
--genotype joint_genotype \
--library_id $library_id \
--modality rna \
--reference reference/refdata-gex-GRCh38-2020-A \
--pseudobulk_counts \
--single_cell_counts \
--celltype_counts \
--isphased \
--threads 20
done

######COUNT WASP ATAC BAM########
mkdir project/analysis/dkd/wasp_atac/counts
for library_id in ${library_ids[*]}
do
bsub -G compute-parkerw -J "acnt_${library_id}" -o "$SCRATCH1/log.acnt_${library_id}" -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/salsa:latest)' \
bash $SCRATCH1/SALSA/step7_gatk_alleleCount.sh \
--inputvcf project/analysis/dkd/vcf/funcotation/$library_id.pass.joint.hcphase.funco.vcf.gz \
--inputbam project/analysis/dkd/wasp_atac/$library_id.hcphase.wasp.bam \
--outputdir project/analysis/dkd/wasp_atac/counts \
--barcodes project/analysis/dkd/barcodes/atac_barcodes.csv \
--genotype joint_genotype \
--library_id $library_id \
--modality atac \
--reference $SCRATCH1/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--pseudobulk_counts \
--single_cell_counts \
--celltype_counts \
--isphased \
--threads 20
done
