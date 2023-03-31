# NYWWS (New York Waste Water Surveillance)

#### Snakemake pipeline to identify and quantify SARS-CoV2 variants in waste water samples from across New York State using Freyja. 

#### DAG of pipeline:
![alt text](https://github.com/YazBraimah/NYWWS/blob/main/workflow.png?raw=true)


## Setup and execution

### Set up rclone

First, download and install rclone, then run `rclone setup` and follow [these instructions](https://rclone.org/googlecloudstorage/) to set it up to access Google Cloud Storage.

### Set up folders and dependencies

After cloning the repository to your local machine, first set up the pipeline's main conda environment.

```
conda env create -f envs/auto.yml -n nywws
```

Here, we gave it a name of `nywws`.

Create a directory with all the BAM files to be analyzed. BAM files are expected to have the extension `*.ptrim.bam`.

```
mkdir ./BAM_DIR/ALL_BAMS
```

### Run pipeline

Run the `scripts/make_bam_json_samples.py` program on the BAM files within the BAM directory (note the escape character before the *):

```
scripts/make_bam_json_samples.py /absolute/path/to/BAM_DIR/ALL_BAMS/\*ptrim.bam
```
This will create a file called `samples_bam.json`. This file path should be specifid in the `config.yml` file under `BAM_SAMPLES_FILE`.

Ensure that the correct paths to required inputs are specified in the `config.yml` file.

Activate the nywws environment with `conda activate nywws` and run the pipeline with this minimal `snakemake` command:

```
snakemake --use-conda -j 50 --cores 4
```

Which will run up to 50 jobs simultaneously, depending on the computing architecture.


## Running the pipeline automatically

- Make sure to set the parameters in the run script properly.
