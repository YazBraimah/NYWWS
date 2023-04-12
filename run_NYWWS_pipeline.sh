#!/usr/bin/env bash


# Parameters
# ----------

# Name of the Google Cloud Storage bucket for rclone
RCLONE_BUCKET_NAME="gcs:su_nywws_test_bucket"

# This file lists the directories in the GCS to download BAM files from
GCS_DIRS_FILE="config/remote_dirs.txt"

# Folder to download BAM files to
BAM_FOLDER="data/raw_bam"

# Name of main conda environment of the pipeline
CONDA_ENV="nywws"

# File with the pipeline parameters
PIPELINE_CONFIG="config/pipeline_parameters.yml"

# How many jobs to run in the pipeline at once?
PIPELINE_JOBS=100

# Upload results after running the pipeline?
UPLOAD_RESULTS=true

# Name of results folder on Google Bucket
RESULTS_FOLDER=$(date +"%d-%m-%Y")


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
      -j ${PIPELINE_JOBS} \
      --use-conda \
      --configfile ${PIPELINE_CONFIG}


# Upload results
# --------------

if [ ${UPLOAD_RESULTS} = true ] ; then
    SOURCE=results
    DEST=${RCLONE_BUCKET_NAME}/SU_Results/Cumulative_Results/${RESULTS_FOLDER}
    rclone --progress --gcs-bucket-policy-only copy ${SOURCE} ${DEST}
fi
