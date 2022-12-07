# NYWWS (New York Waste Water Surveillance)

#### Snakemake pipeline to identify and quantify SARS-CoV2 variants in waste water samples from across New York State. The pipline is adapted from C-WAP , a Nextflow workflow developed by the United States Food and Drug Administration, Center for Food Safety and Applied Nutrition.

#### DAG of pipeline:
![alt text](https://github.com/YazBraimah/NYWWS/blob/main/workflow.png?raw=true)


### Setup and execution

1. After cloning the repository to your local machine, first set up the base conda environment, which is in the `envs/nywws_env.yml`:

```
conda env create -f envs/nywws_env.yml
```

2. Create a directory with all the BAM files to be analyzed. BAM files are expected to have the extension `*.ptrim.bam`.

```
mkdir ./BAM_DIR/
```

3. Run the `scripts/make_bam_json_samples.py` program on the BAM files within the BAM directory (note the escape character before the *):

```
scripts/make_bam_json_samples.py /absolute/path/to/BAM_DIR/\*ptrim.bam
```
This will create a file called `samples_bam.json`. This file path should be specifid in the `config.yml` file under `BAM_SAMPLES_FILE`.

4. Ensure that the paths to the remaining dependancies (e.g., the C-WAP and LCS directories) are specified in the `config.yml` file.

5. Activate the nywws environment with `conda activate nywws` and run the pipeline with this minimal `snakemake` command:

```
snakemake --use-conda -j 10
```
Which will run 10 jobs simultaneously, depending on the computing architecture.
