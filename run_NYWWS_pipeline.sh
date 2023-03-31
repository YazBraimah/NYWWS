#!/usr/bin/env bash


# Parameters
# ----------

# Name of the Google Cloud Storage bucket for rclone
RCLONE_BUCKET_NAME="gcs:su_nywws_test_bucket"

# This file lists the directories in the GCS to download BAM files from
GCS_DIRS_FILE="config/remote_dirs.txt"

# Folder to download BAM files to
BAM_FOLDER="BAM_DATA"

# Name of main conda environment of the pipeline
CONDA_ENV="nywws"

# File with the pipeline parameters
PIPELINE_CONFIG="config/pipeline_parameters.yml"

# Upload files after running the pipeline?
# UPLOAD_FILES=yes or no


# Download BAM files
# ------------------

# Change "rclone copy" to "rclone sync" to also delete local BAMs
# that are deleted in the bucket.
cat ${GCS_DIRS_FILE} | while read -r DIR
do
    SOURCE="${RCLONE_BUCKET_NAME}/${DIR}"
    DEST="${BAM_FOLDER}/${DIR}"
    mkdir -p ${DEST}
    rclone copy \
	   --progress \
	   ${SOURCE} ${DEST} \
	   --include "/202*.ptrim.bam"
done


# Run the pipeline
# ----------------

# Run the pipeline with the conda environment
conda run -n ${CONDA_ENV} snakemake \
      -j 50 \
      --use-conda \
      --configfile ${PIPELINE_CONFIG}


# Upload results
# --------------
# 
# mkdir /home/yahmed/GCS/su_nywws_test_bucket/SU_Results/Cumulative_Results/"$(date +"%d-%m-%Y")"
# cp -r /home/yahmed/NYWWS/Cumulative_Results/BAM/ \
# /home/yahmed/NYWWS/Cumulative_Results/Coverage/ \
# /home/yahmed/NYWWS/Cumulative_Results/fastQC/ \
# /home/yahmed/NYWWS/Cumulative_Results/Freyja/ \
# /home/yahmed/NYWWS/Cumulative_Results/iVar/ \
# /home/yahmed/NYWWS/Cumulative_Results/MultiQC/ \
# /home/yahmed/NYWWS/Cumulative_Results/Pangolin/ \
# /home/yahmed/NYWWS/Cumulative_Results/Sequences/ \
# /home/yahmed/NYWWS/Cumulative_Results/Summary/ \
# ~/GCS/su_nywws_test_bucket/SU_Results/Cumulative_Results/"$(date +"%d-%m-%Y")"
