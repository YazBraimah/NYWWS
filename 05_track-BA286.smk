
from pathlib import Path

import pandas as pd

COVERAGE_STATUS = "output/qc/coverage/sample_coverage_status.tsv"
MUT_GROUP_FOLDER = Path(config["mutation_group_folder"])


rule all:
    input:
        "output/results/BA.2.86-report/mutation-group-key.txt",
        "output/results/BA.2.86-report/BA.2.86-report.tsv"


def mutation_group_files(wildcards):
    files = list(MUT_GROUP_FOLDER.glob("*.csv"))
    return files


rule mutation_group_key:
    input: mutation_group_files
    output: "output/results/BA.2.86-report/mutation-group-key.txt"
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
    for group in mutation_groups + ["total"]:
        input_files["counts:" + group] = expand(
            f"output/BA.2.86-read-counts/counts/{group}/{{sample}}.txt",
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
    output: "output/results/BA.2.86-report/BA.2.86-report.tsv"
    params:
        min_reads = 5,
        variant = "BA.2.86"
    script: "scripts/cryptic-lineage-report.py"


def bam_for_counting_reads(wildcards):
    if wildcards.group == "total":
        return f"output/covid-filtered/{wildcards.sample}/{wildcards.sample}.bam"
    else:
        return f"output/BA.2.86-read-counts/freyja-extract/{wildcards.group}/{wildcards.sample}.bam"


rule count_reads:
    input: bam_for_counting_reads
    output: "output/BA.2.86-read-counts/counts/{group}/{sample}.txt"
    conda: "envs/freyja.yml"
    shell:
        "samtools view {input} | wc -l > {output}"


rule freyja_extract:
    input:
        bam = "output/covid-filtered/{sample}/{sample}.bam",
        mutations = str(MUT_GROUP_FOLDER / "{group}.csv")
    output: "output/BA.2.86-read-counts/freyja-extract/{group}/{sample}.bam"
    conda: "envs/freyja.yml"
    shell:
        "freyja extract --refname 2019-nCoV {input.mutations} {input.bam} --output {output} --same_read"
