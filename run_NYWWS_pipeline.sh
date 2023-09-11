#!/usr/bin/env bash

set -e

# Parameters
# ----------

# Upload results after running the pipeline?
UPLOAD_RESULTS=true

# Name of the rclone remotes
RCLONE_BUCKET_NAME="gcs:su_nywws_test_bucket"
# The name of the OneDrive remote in rclone is hardcoded below,
# because I can't figure out how to deal with spaces in the dir name.
#RCLONE_ONEDRIVE_NAME="onedrive:CDC wastewater data"

# This file lists the directories in the GCS to download BAM files from
GCS_DIRS_FILE="config/remote_dirs.txt"

# Folder to download BAM files to
BAM_FOLDER="data/raw_bam"

# Name of main conda environment of the pipeline
CONDA_ENV="nywws"

# File with the pipeline parameters
PIPELINE_CONFIG="config/pipeline_parameters.yml"

# Download BAM files
# ------------------

cat ${GCS_DIRS_FILE} | while read -r DIR
do
    SOURCE="${RCLONE_BUCKET_NAME}/${DIR}"
    DEST="${BAM_FOLDER}/${DIR}"
    mkdir -p ${DEST}
    rclone sync \
	   --progress \
	   ${SOURCE} ${DEST} \
	   --include "/202*.ptrim.bam"
done


# Update data repository
# ----------------------

cd ../NYS-WWS-Data
git pull
cp sars2-concentration.csv ../20230403_freyja-pipeline/data/sample_metadata/nys-wws-sars2-concentration.csv
cp metadata/lineage-map.csv ../20230403_freyja-pipeline/data/sample_metadata/lineage_info.csv
cp metadata/variants-of-concern.csv ../20230403_freyja-pipeline/data/sample_metadata/variants_of_concern.csv
cp nys-wws-sewersheds.csv ../20230403_freyja-pipeline/data/sample_metadata/sewershed_metadata.csv
cd ../20230403_freyja-pipeline

cp ../20230504_sample-tracking/output/sars2-sequencing-manifest.csv data/sample_metadata


# Run the pipeline
# ----------------
# Unclear to me why it needs to run twice. It updates some files the first
# time, and then re-runs all the Freyja results the second time.

# Allows activating conda env
source /home/iavascon/miniconda3/bin/activate
conda activate nywws

echo ""
echo "BAM check"
echo ""

snakemake \
    --snakefile 01_bam-check.smk \
    -c1 \
    --use-conda \
    --configfile ${PIPELINE_CONFIG}

echo ""
echo "Quality control"
echo ""

snakemake \
    --snakefile 02_quality-control.smk \
    -c20 \
    --use-conda \
    --configfile ${PIPELINE_CONFIG}

echo ""
echo "Freyja"
echo ""

snakemake \
    --snakefile 03_freyja.smk \
    --forcerun freyja_update \
    -c20 \
    --use-conda \
    --configfile ${PIPELINE_CONFIG}

echo ""
echo "Aggregate results"
echo ""

snakemake \
    --snakefile 04_aggregate.smk \
    -c20 \
    --use-conda \
    --configfile ${PIPELINE_CONFIG}

echo ""
echo "Track cryptic variants"
echo ""

snakemake \
    --snakefile 05_cryptic-variants.smk \
    -c20 \
    --use-conda \
    --configfile ${PIPELINE_CONFIG}


# Upload results
# --------------

if [ ${UPLOAD_RESULTS} = true ] ; then

    # Upload to GitHub
    cd ../NYS-WWS-Data
    git pull
    DEST=genetic-sequencing-history/$(date +"%Y%m%d")
    mkdir -p ${DEST}
    cp ../20230403_freyja-pipeline/output/results/sample_info.tsv ${DEST}
    cp ../20230403_freyja-pipeline/output/results/freyja_parse_barcode.csv ${DEST}
    cp ../20230403_freyja-pipeline/output/results/comprehensive_results_table.txt ${DEST}
    cp ../20230403_freyja-pipeline/output/results/freyja_parse.csv ./sars2-genetic-sequencing.csv
    git add .
    git commit -m "Genetic sequencing update for $(date +"%d %B %Y")"
    git push
    cd ../20230403_freyja-pipeline

    # Upload to OneDrive
    SOURCE=output/covid-filtered-BAMs
    rclone --progress --copy-links sync ${SOURCE} onedrive:'CDC wastewater data'/BAM-Files
    rclone --progress copyto output/results/SRA_table.csv onedrive:'CDC wastewater data'/BAM-Files/SRA_table.csv

    # Upload to Amazon S3
    rclone --progress copyto output/results/var.data_summary.rds s3:nystatewws/var.data_summary.rds
    rclone --progress copyto output/results/sample-id-report.tsv s3:nystatewws/covid-sample-id-report.tsv
fi
