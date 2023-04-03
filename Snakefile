from pathlib import Path

import pandas as pd


BAM_FOLDER = Path(config["bam_folder"])
SAMPLE_METADATA = Path(config["sample_metadata"])
BAM_HEADER = Path(config["bam_header"])
COVID_REFSEQ = Path(config["covid_refseq"])


rule all:
    input: "output/done"


def samples_with_matching_ids(wildcards):
    existence_path = checkpoints.check_file_existence.get(**wildcards).output[0]
    existence_df = pd.read_table(existence_path)
    present_samples = existence_df.sample_id[existence_df.sample_present == "ok"]
    return [
        f"results/fastQC/{sample}_fastqc.html"
        for sample in present_samples
    ]

rule all_filter_covid_reads:
    input: samples_with_matching_ids
    output: "output/done"
    message: "Done running fastqc on all samples."
    shell: "touch {output}"


rule qualimap:
    input:
        bam = "results/FastQs/{sample}/{sample}.bam",
        gtf = COVID_REFSEQ
    output: "results/BAM/{sample}/{sample}.qualimap/qualimapReport.html"
    log: "results/BAM/{sample}/qualimap.log"
    threads: 8
    resources:
        mem_mb=32000
    message: "{wildcards.sample}: Evaluating BAM mapping quality with QualiMap."
    run:
        shell('qualimap bamqc'
                ' -bam {input.bam}'
                ' --java-mem-size=32G'
                # ' -gff {input.gtf}'
                ' -outdir ' + join(OUT_DIR, 'BAM', '{wildcards.sample}', '{wildcards.sample}.qualimap') +
                ' > {log} 2>&1')


rule fastqc:
    input: "results/FastQs/{sample}/{sample}.fq.gz"
    output: "results/fastQC/{sample}_fastqc.html"
    log: "results/fastQC/{sample}.fastQC_se.log"
    conda: "envs/fastqc_env.yml"
    message: "{wildcards.sample}: Checking read quality with fastqc."
    shell: "fastqc -o results/fastQC {input} > {log} 2>&1"


rule bam_to_fastq:
    input: "results/FastQs/{sample}/{sample}.sam"
    output: "results/FastQs/{sample}/{sample}.fq.gz"
    message: "{wildcards.sample}: Converting SAM to FastQ."
    script: "scripts/sam2fastq.py"


rule filter_covid_reads:
    input:
        bam = "output/samples/{sample}.bam",
        header = BAM_HEADER
    output:
        bam = "results/FastQs/{sample}/{sample}.bam",
        sam = "results/FastQs/{sample}/{sample}.sam",
        stats_original = "results/BAM/{sample}/{sample}.cssorted.stats",
        stats_filtered = "results/BAM/{sample}/{sample}.resorted.stats"
    threads: 2
    conda: "envs/bioinfo.yml"
    message: "Filtering {wildcards.sample} to only COVID reads."
    shell:
        "samtools index {input.bam} ; "
        "samtools view -b --threads {threads} {input} 2019-nCoV "
        "| samtools reheader {input.header} - "
        "| samtools sort -o {output.bam} -@ {threads} ; " 
        "samtools view {output.bam} -o {output.sam} ; "
        "samtools stats {input.bam} | grep ^SN | cut -f 2- > {output.stats_original} ; "
        "samtools stats {output.bam} | grep ^SN | cut -f 2- > {output.stats_filtered} ; "


checkpoint check_file_existence:
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
