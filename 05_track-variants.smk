
from pathlib import Path

import pandas as pd

COVERAGE_STATUS = "output/qc/coverage/sample_coverage_status.tsv"
MUT_GROUP_FOLDER = Path(config["mutation_group_folder"])
VARIANTS_TO_TRACK = config["variants_to_track"]


rule all:
    input:
        expand(
            "output/results/variant-tracking-reports/{variant}/{result}",
            variant=VARIANTS_TO_TRACK,
            result=["report.tsv", "mutation-group-key.txt"]
        )


def mutation_group_files(wildcards):
    folder = MUT_GROUP_FOLDER / wildcards.variant
    files = list(folder.glob("*.csv"))
    return files


rule mutation_group_key:
    input: mutation_group_files
    output: "output/results/variant-tracking-reports/{variant}/mutation-group-key.txt"
    run:
        with open(output[0], "w") as f:
            for input_file in input:
                group_name = Path(input_file).stem
                with open(input_file) as g:
                    contents = g.read()
                f.write(f"=> {group_name}\n{contents}\n")
    

def read_counts(wildcards):
    valid_samples = (
        pd.read_table(COVERAGE_STATUS)
        .query("enough_coverage == 'yes'")
        .sample_id
    )
    input_files = dict()
    mutation_groups = [f.stem for f in mutation_group_files(wildcards)]
    for group in mutation_groups:
        input_files["counts:" + group] = expand(
            f"output/variant-read-counts/counts_{wildcards.variant}/{group}/{{sample}}.txt",
            sample=valid_samples
        )
    input_files["counts:total"] = expand(
        "output/variant-read-counts/counts_total/{sample}.txt",
        sample=valid_samples
    )
    return input_files


rule lineage_report:
    input:
        unpack(read_counts),
        concentration = "data/sample_metadata/nys-wws-sars2-concentration.csv",
        sewershed = "data/sample_metadata/sewershed_metadata.csv",
        freyja = "output/results/freyja_parse.csv",
        coverage = "output/qc/coverage/coverage_report.tsv"
    output: "output/results/variant-tracking-reports/{variant}/report.tsv"
    params:
        min_reads = 5,
        variant = lambda wildcards: wildcards.variant
    script: "scripts/variant-tracking-report.py"


rule count_reads_total:
    input: "output/covid-filtered/{sample}/{sample}.bam"
    output: "output/variant-read-counts/counts_total/{sample}.txt"
    conda: "envs/freyja.yml"
    shell:
        "samtools view {input} | wc -l > {output}"


rule count_reads:
    input: "output/variant-read-counts/freyja-extract_{variant}/{group}/{sample}.bam"
    output: "output/variant-read-counts/counts_{variant}/{group}/{sample}.txt"
    conda: "envs/freyja.yml"
    shell:
        "samtools view {input} | wc -l > {output}"


rule freyja_extract:
    input:
        bam = "output/covid-filtered/{sample}/{sample}.bam",
        mutations = str(MUT_GROUP_FOLDER / "{variant}" / "{group}.csv")
    output: "output/variant-read-counts/freyja-extract_{variant}/{group}/{sample}.bam"
    conda: "envs/freyja.yml"
    shell:
        "freyja extract --refname 2019-nCoV {input.mutations} {input.bam} --output {output} --same_read"
