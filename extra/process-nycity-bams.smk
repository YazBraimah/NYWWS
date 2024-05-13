
from pathlib import Path

BAM_FOLDER = Path("data/nycity_before_processing")

def processed_bam_files(wildcards):
    samples = BAM_FOLDER.glob("**/*.bam")
    return [f"data/raw_bam/nycity/{bam.stem}.ptrim.bam" for bam in samples]

rule all:
    input: processed_bam_files


def find_file(wildcards):
    search = list(BAM_FOLDER.glob(f"**/{wildcards.sample_id}.bam"))
    assert len(search) == 1
    return search[0]

rule process_file:
    input: find_file
    output: "data/raw_bam/nycity/{sample_id}.ptrim.bam"
    conda: "../envs/bioinfo.yml"
    shell: "samtools reheader -c \'sed -e \"s/SN:NC_045512.2/SN:2019-nCoV/\"\' {input} > {output}"
