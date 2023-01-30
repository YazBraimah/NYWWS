#!/usr/bin/env bash


### Copy new BAM files to "ALL_BAMS directory"

echo -e "\n\t____________________________________\n"
echo -e "\n\t--- Updating BAM files directory ---\n"
echo -e "\n\t____________________________________\n"

cp -u /home/yahmed/GCS/su_nywws_test_bucket/nymc/202*ptrim.bam \
/home/yahmed/GCS/su_nywws_test_bucket/rochester/nywws/202*ptrim.bam \
/home/yahmed/GCS/su_nywws_test_bucket/suny_upstate/202*ptrim.bam \
/home/yahmed/GCS/su_nywws_test_bucket/wadsworth/202*ptrim.bam \
/home/yahmed/NYWWS/BAM_DATA/ALL_BAMS/


### Update the sample BAM JSON file

echo -e "\n\t________________________________\n"
echo -e "\n\t--- Updating BAM sample file ---\n"
echo -e "\n\t________________________________\n"
/home/yahmed/NYWWS/scripts/make_bam_json_samples.py /home/yahmed/NYWWS/BAM_DATA/ALL_BAMS/202*ptrim.bam


### Run Snakemake pipeline

snakemake --use-conda -j 50


### copy results folders to GCS su_nywws_test_bucket
