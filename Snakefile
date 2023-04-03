from pathlib import Path

bam_folder = Path(config["bam_folder"])
sample_metadata = Path(config["sample_metadata"])


rule all:
    input: "output/sample_info/file_existence.tsv"


rule check_file_existence:
    input:
        manifesto = "output/sample_info/samples_bam.json",
        metadata = sample_metadata
    output: "output/sample_info/file_existence.tsv"
    script: "scripts/check_file_existence.py"


rule make_samples_list:
    input: bam_folder
    output: "output/sample_info/samples_bam.json"
    script: "scripts/make_bam_json_samples.py"
