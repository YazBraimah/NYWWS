#!/usr/bin/env bash


# Parameters
# ----------

# Name of the rclone remotes
RCLONE_BUCKET_NAME="gcs:su_nywws_test_bucket"
# Hardcoded below because I can't figure out how to deal with spaces in the dir name.
#RCLONE_ONEDRIVE_NAME="onedrive:CDC wastewater data"

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


# Update data repository
# ----------------------

cd ../NYS-WWS-Data
git pull
cp nys-wws-sars2-concentration.csv ../NYWWS/data/sample_metadata
cd ../NYWWS


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
    # Upload to GitHub
    cd ../NYS-WWS-Data
    git pull
    DEST=SU_results/$(date +"%Y%m%d")
    mkdir -p ${DEST}
    cp ../NYWWS/results/Summary/sample_info.tsv ${DEST}
    cp ../NYWWS/results/Freyja/Aggregate/freyja_parse.csv ${DEST}
    git add ${DEST}
    git commit -m "SU results for $(date +"%d %B %Y")"
    git push
    cd ../NYWWS

    # Upload to OneDrive
    SOURCE=output/BAM_processed
    rclone --progress --copy-links sync ${SOURCE} onedrive:'CDC wastewater data'/BAM-Files
    rclone --progress copyto results/Summary/SRA_table.csv onedrive:'CDC wastewater data'/BAM-Files/SRA_table.csv

    # Upload to GCS Bucket
    SOURCE=results
    DEST=${RCLONE_BUCKET_NAME}/SU_Results/Cumulative_Results/${RESULTS_FOLDER}
    rclone --progress --gcs-bucket-policy-only copy ${SOURCE} ${DEST}
fi
