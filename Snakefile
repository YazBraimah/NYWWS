from pathlib import Path

bam_folder = Path(config["bam_folder"])


rule all:
    input: "results/samples_bam.json"


rule make_samples_list:
    input: bam_folder
    output: "results/samples_bam.json"
    script: "scripts/make_bam_json_samples.py"
