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

Download the samples CSV file from GitHub.

### Run pipeline

Edit `config/pipeline_parameters.yml` with the correct pipeline parameters and file paths. Then, run the pipeline with this command (assuming the conda environment is called `nywws`):

```
conda run -n nywws snakemake --use-conda -j 50 --cores 4
```

Which will run up to 50 jobs simultaneously, depending on the computing architecture.


## Running the pipeline automatically

For automatically downloading and pushing files from GitHub, you will need to have the [GitHub command-line interface](https://github.com/cli/cli) installed. After doing so, you will need to create a [Personal Access Token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token). Once your token is created, store it in the environment variable `GH_TOKEN`. Note that storing credentials in a variable makes them vulnerable if the computer is compromised.

Set the parameters in the run script properly, then run it: `bash run...`.

Set up a cron job to run the script.
