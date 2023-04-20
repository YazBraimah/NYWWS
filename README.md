# NYWWS (New York Waste Water Surveillance)

This is a Snakemake pipeline to identify and quantify SARS-CoV2 variants in waste water samples from across New York State using Freyja. 

## Setup and execution

### Set up rclone

First, download and install rclone, then run `rclone setup` and follow [these instructions](https://rclone.org/googlecloudstorage/) to set it up to access Google Cloud Storage.

### Set up dependencies

You need two repositories cloned to your local machine, this one and the data repo.

After cloning the repositories, first set up the pipeline's main conda environment.

```
cd NYWWS
conda env create -f envs/nywws.yml -n nywws
```

Here, we gave it a name of `nywws`.

Finally, you need the GitHub CLI, `gh`, installed as well.

### Run pipeline

Edit `config/pipeline_parameters.yml` with the correct pipeline parameters and file paths. Then, edit `run_NYWWS_pipeline.sh` with the run parameters. Login to your GitHub account by running `gh login`. Finally, run the pipeline script:

```
bash run_NYWWS_pipeline.sh
```
