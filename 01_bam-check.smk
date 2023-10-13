BAM_FOLDER = Path(config["bam_folder"])
SAMPLE_METADATA = Path(config["sample_metadata"])

rule all:
    input: "output/sample_info/file_existence.tsv"


rule check_file_existence:
    input:
        bam_folder = BAM_FOLDER,
        metadata = SAMPLE_METADATA,
        corrupt_files = "output/sample_info/corrupt-files.txt"
    output: "output/sample_info/file_existence.tsv"
    message: "Matching downloaded BAMs with known sample IDs."
    script: "scripts/check_file_existence.py"


rule check_corrupt_files:
    input: BAM_FOLDER
    output: "output/sample_info/corrupt-files.txt"
    message: "Checking for corrupt BAM files..."
    conda: "envs/bioinfo.yml"
    shell:
        "set +e ; "
        "samtools quickcheck -qv {input}/*/*.ptrim.bam > {output} ; "
        "RESULT=$?" # last command must return non-0 exit code
