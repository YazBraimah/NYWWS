# NYWWS (New York Waste Water Surveillance)

This is a Snakemake pipeline to identify and quantify SARS-CoV2 variants in waste water samples from across New York State using Freyja. 


## Installation instructions

### 1. GitHub command line

You need the GitHub CLI, `gh`, installed. Login to your GitHub account by running `gh login`; this allows the user running the pipeline to interact with the private GitHub repo of the project.

### 2. Folder structure

The most complete version of the pipeline should have all of these folders under the same directory:

- `20230403_freyja-pipeline/`: Contains the main pipeline scripts
- `20230504_sample-tracking/`: Script that tracks upload times of BAM files, designed to run every 10 minutes
- `20240220_bam-renames/`: Scripts used to rename BAM files that have the wrong IDs
- `NYS-WWS-Data/`: Private GitHub repo of the project

The first three should have been given to you. The GitHub repo needs to be pulled, and will require authentication because it's private. (This should be taken care of by the `gh` utility.)

It's important that those folders are all under the same directory, as different scripts will `cd` back and forth between them.

### 3. rclone

Rclone is the software used by the pipeline to interact with external cloud storage. You'll need to install it, then run `rclone setup` to configure three different cloud services that the pipeline talks to. The configuration for each will involve using a browser to log in to your account for each provider; rclone uses that account information to authenticate your access to read and write cloud files.

- Follow [these instructions](https://rclone.org/googlecloudstorage/) to set up access Google Cloud Storage, where the BAM files are uploaded to. Name the remote `gcs`.
- Follow [these instructions](https://rclone.org/s3/#configuration) to set up access to Amazon S3, where some of the output files go. Name the remote `s3`.
- Follow [these instructions](https://rclone.org/onedrive/) to set up access to Microsoft OneDrive, where some of the output files go as well. Name the remote `onedrive`.

### 4. conda

Software required for the pipeline is managed with conda environments, so you will need conda installed (I recommend either the miniconda or micromamba distributions). Once you have it installed, you can use it to install the main virtual environment, which needs to be named `nywws`. The following will install the environment:

```bash
cd 20240403_freyja-pipeline
conda env create -f envs/nywws.yml -n nywws
```

### 5. Google Cloud CLI

You need the Google Cloud CLI installed to have access to the `gcloud` command used by the sample tracking script.


## Configuration

### Pipeline parameters

Basic configuration parameters are at the top of the `20230403_freyja-pipeline/run_NYWWS_pipeline.sh` script, under the heading of "Parameters". Of those, the only mandatory change is the `CONDA_ACTIVATE` parameter, which needs to be changed with the location of your own conda installation.

Other parameters that the pipeline uses are in `20230403_freyja-pipeline/config/pipeline_parameters.yml`. These can be left mostly alone, as they refer to the internal folder structure of the pipeline. The most interesting parameters here are the QC thresholds for sample coverage; if you want to change the threshold we've been using (at least 50% of the genome at 20x or more), this is the place to do it.

Finally, `20230403_freyja-pipeline/config/remote_dirs.txt` has a list of the folders in the GCS bucket that the pipeline pulls BAM files from. This file can be changed to add or remove sequencing labs. Since New York City uses a different folder structure than the other labs, it doesn't appear in this file. Instead, NYC is treated in a special way in the `20230403_freyja-pipeline/run_NYWWS_pipeline.sh` script. 

### Sample tracking

The sample tracking script is `20230504_sample-tracking/run.sh`. It needs to be edited manually to reflect your own installation of the Google Cloud SDK and conda, as well as the path of the folder in your system. (The places to edit are indicated in the script.)

This script was designed to run every 10 minutes as a cron job. All you need to add to the crontab is the full path to the `run.sh` script.

### Other modifications

These other changes are all optional and require manual editing of various scripts that assume hard-coded values.

**Folder structure:** Both the main pipeline script and the sample tracking script assume the folder names described in the 
"Folder structure" section above. Edit them if your folder names are different.

**Rclone remote names:** The pipeline, sample tracking, and BAM rename scripts all assume that rclone remotes are named `gcs`, `s3`, and `onedrive`; you can pick different names when setting up rclone, all of those scripts need to be updated.

**Conda environment name:** The pipeline and sample tracking scripts assume that the conda environment used to run the pipeline is called `nywws`. To change that, look for the lines of code that activate the environment (start with `conda activate`).


## Execution

### Running the pipeline

The main pipeline just needs to be run inside its own folder.

```bash
cd 20230403_freyja-pipeline
bash run_NYWWS_pipeline.sh
```

The pipeline script relies on the workflow manager Snakemake to make sure that all necessary steps are run for newly downloaded BAM files, and no redundant steps are re-run for previous files. For safety, the script will crash as soon as there is any error anywhere. So if there's a crash, you need to re-run the script after fixing the cause. Again, Snakemake makes sure that any steps that ran before successfully before the crash don't need to re-run.

The various outputs produced by the pipeline will be in `20230403_freyja-pipeline/output`, while the main results files that are used downstream will be in `20230403_freyja-pipeline/output/results`.

With the parameter `UPLOAD_RESULTS=true` in the `run_NYWWS_pipeline.sh` script, the pipeline uploads results to GitHub and the cloud after it's done executing. If you want to do a local test run that does not upload, just set the parameter to `UPLOAD_RESULTS=false`.

### Sample tracking

Once the script is set as a cron job, you don't need to do anything. It will run automatically and upload its results to the Amazon S3 storage. You can run the script manually as well if you like:

```bash
cd 20230504_sample-tracking/
bash run.sh
```

If you have GNUplot installed, you can also run `bash plot.sh` to create plain-text plots of file uploads by time in `output/plots.txt`. This is useful for sanity checks, but these plots are not uploaded anywhere or used otherwise.

### Renaming BAM files

BAM files with the wrong sample ID in the Google Cloud Storage bucket need to be renamed. I've done this manually, as the process is irreversible. The sub-folders under `20240220_bam-renames/` are records of past renames for reference, and can be used as a base for a new batch of renames.

To perform BAM renames, first create a new subfolder, e.g. `20240220_bam-renames/2024XXYY_renames/`. Then, move the provided table of renames inside it (these are usually called `bam-name-fixes-2024-XX-YY.txt`). Next, run the following:

```bash
cut -f 1 bam-name-fixes-2024-XX-YY.txt | sed '1d' | awk 'BEGIN { OFS = "" } { print "**/", $1, ".ptrim.bam" }' > 01_bams-to-search.txt
rclone lsf -R --absolute --files-only gcs:su_nywws_test_bucket --include-from 01_bams-to-search.txt > 02_bams-in-bucket.txt
```

(You'll need to change the rclone command if your remote is not named `gcs`.)

Next, copy the script `20240220_bam-renames/03_make-commands.py` into your working folder, then edit it with the current filename of the table of BAM renames. After the script is edited, you need to activate the conda environment for the pipeline and run the script:

```bash
conda activate nywws
python 03_make-commands.sh
```

The script produces a final script, `04_commands.sh`. Running this script is the irreversible step that executes the renames in the bucket. Make sure that the produced commands are correct, then run them with `bash 04_commands.sh`.
