# # export docker volumes
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $STORAGE1/reference:$HOME/reference \
# $SCRATCH1:$SCRATCH1"

# # to run interactive
# bsub -Is -G compute-parkerw -R 'rusage[mem=32GB]' -q general-interactive -a 'docker(p4rkerw/ldsc:1.0)' /bin/bash

# SCRATCH1=/mnt/g/scratch
# docker run -it \
# -v /mnt/g/diabneph:$HOME/project \
# -v /mnt/g/reference:$HOME/reference \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# --workdir $HOME \
# p4rkerw/ldsc:latest

# lots of traits located here: https://alkesgroup.broadinstitute.org/UKBB/ 
# and here https://alkesgroup.broadinstitute.org/sumstats_formatted/

wget -P $SCRATCH1 https://ckdgen.imbi.uni-freiburg.de/files/Stanzick2021/metal_eGFR_meta1.TBL.map.annot.gc.gz
gunzip $SCRATCH1/metal_eGFR_meta1.TBL.map.annot.gc.gz

# use the following header
# MarkerName chr position Allele1 Allele2 Effect P-value n
# Note: do not use the original "MarkerName" column without reformatting it. Solution is to substitute it for the RSID in column 15
cat $SCRATCH1/metal_eGFR_meta1.TBL.map.annot.gc | awk 'BEGIN { OFS="\t" } {print $15,$11,$12,$2,$3,$6,$13,$4}' > /tmp/sumstats.txt

# fix header for p-value
{ echo "MarkerName chr position Allele1 Allele2 Effect P-value n"; tail -n+2 /tmp/sumstats.txt; } > /tmp/header_sumstats.txt

# activate ldsc conda environ
export PATH=$PATH:/opt/conda/envs/ldsc/bin:/opt/conda/condabin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
source activate ldsc

wget -P $SCRATCH1 https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bzip2 -d $SCRATCH1/w_hm3.snplist.bz2

# munge egfr sumstats using ldsc
python /opt/ldsc/munge_sumstats.py \
--sumstats /tmp/header_sumstats.txt \
--merge-alleles $SCRATCH1/w_hm3.snplist \
--out $SCRATCH1/eGFR_munge \
--chunksize 500000

# ckd gwas statistics https://www.nature.com/articles/s41588-019-0407-x
wget -P $SCRATCH1 https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/CKD_overall_ALL_JW_20180223_nstud30.dbgap.txt.gz
gzip -dc $SCRATCH1/ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/CKD_overall_ALL_JW_20180223_nstud30.dbgap.txt.gz | sed 's/n_total_sum/n/g' > /tmp/ckd_gwas.txt
# munge t2d sumstats using ldsc
python /opt/ldsc/munge_sumstats.py \
--sumstats /tmp/ckd_gwas.txt \
--merge-alleles $SCRATCH1/w_hm3.snplist \
--out $SCRATCH1/CKD_munge \
--chunksize 500000

# microalbuminuria overall trans-ethnic
wget -P $SCRATCH1 https://ckdgen.imbi.uni-freiburg.de/files/Teumer2019/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz
gzip -dc $SCRATCH1/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz | \
	sed 's/n_total_sum/n/g' | \
	awk 'BEGIN { OFS="\t" } {print $3,$4,$5,$7,$9,$10}' > /tmp/micro_gwas.txt
python /opt/ldsc/munge_sumstats.py \
--sumstats /tmp/micro_gwas.txt \
--merge-alleles $SCRATCH1/w_hm3.snplist \
--out $SCRATCH1/MICRO_munge \
--chunksize 500000

# urine sodium
# ALLELE1 is effect allele "A2" and ALLELE2 is ref allele and N was obtained from the study sample size
# munge_sumstats will auto-detect A1 as the ref allele
wget -P $SCRATCH1 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008647/PazokiR_prePMID_Sodium.GWAS.csv
cat $SCRATCH1/PazokiR_prePMID_Sodium.GWAS.csv | sed 's/,/\t/g' | sed 's/P_BOLT_LMM_INF/P/' | sed 's/ALLELE0/A1/' | sed 's/ALLELE1/A2/' | \
	awk 'NR==1{$(NF+1)="n"} NR>1{$(NF+1)="446237"}1' > /tmp/urinena_gwas.txt

python /opt/ldsc/munge_sumstats.py \
--sumstats /tmp/urinena_gwas.txt \
--merge-alleles $SCRATCH1/w_hm3.snplist \
--out project/analysis/dkd/ldsc/sumstats/URINENA_munge \
--chunksize 500000
