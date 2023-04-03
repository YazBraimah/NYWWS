from pathlib import Path

BAM_FOLDER = Path(config["bam_folder"])
SAMPLE_METADATA = Path(config["sample_metadata"])


rule all:
    input: "output/sample_info/file_existence.tsv"


rule check_file_existence:
    input:
        bam_folder = "output/samples",
        metadata = SAMPLE_METADATA
    output: "output/sample_info/file_existence.tsv"
    message: "Matching downloaded BAMs with known sample IDs."
    script: "scripts/check_file_existence.py"


def get_raw_bams(wildcards):
    return list(BAM_FOLDER.glob("**/*.ptrim.bam"))

rule symlink_raw_bams:
    input: get_raw_bams
    output: directory("output/samples")
    message: "Creating symbolic links to raw BAM files."
    script: "scripts/symlink_raw_bams.py"
