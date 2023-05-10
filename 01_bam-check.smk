BAM_FOLDER = Path(config["bam_folder"])
SAMPLE_METADATA = Path(config["sample_metadata"])

rule check_file_existence:
    input:
        bam_folder = BAM_FOLDER,
        metadata = SAMPLE_METADATA
    output: "output/sample_info/file_existence.tsv"
    message: "Matching downloaded BAMs with known sample IDs."
    script: "scripts/check_file_existence.py"
