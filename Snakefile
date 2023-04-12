from pathlib import Path

import pandas as pd


BAM_FOLDER = Path(config["bam_folder"])
SAMPLE_METADATA = Path(config["sample_metadata"])
BAM_HEADER = Path(config["bam_header"])
FASTA_REFERENCE = Path(config["fasta_ref"])


rule all:
    input:
        "results/Summary/sample_info.tsv"


rule sample_info_report:
    input:
        existence = "output/sample_info/file_existence.tsv",
        coverage = "results/Coverage/sample_coverage_status.tsv",
        freyja = "output/sample_info/freyja_quality.tsv"
    output: "results/Summary/sample_info.tsv"
    script: "scripts/sample_status.py"


rule Freyja_parse:
    input: "results/Freyja/Aggregate/aggregated_results.tsv"
    output:
        longform = "results/Freyja/Aggregate/freyja_parse.csv",
        report = "output/sample_info/freyja_quality.tsv"
    threads: 2
    resources:
        mem_mb=4000
    message: "Parsing Freyija output."
    script: "scripts/parse_freyja_demix.py"


def freyja_demix_samples(wildcards):
    coverage_path = checkpoints.coverage_summary.get(**wildcards).output["cov_per_sample"]
    coverage_df = pd.read_table(coverage_path)
    samples_with_coverage = coverage_df.sample_id[coverage_df.enough_coverage == "yes"]
    return expand(
        "results/Freyja/Demix/Results/{sample}_freyja.demix",
        sample=samples_with_coverage
    )

rule Freyja_aggregate:
    input: freyja_demix_samples
    output:
        bt2_agg = "results/Freyja/Aggregate/aggregated_results.tsv",
        bt2_png = "results/Freyja/Aggregate/freyja_stacked_barplots.png"
    threads: 8
    resources:
        mem_mb=32000
    conda: "envs/freyja.yml"
    message: "Running Freyja Aggregate and outputting barplots."
    shell:
        "freyja aggregate"
        "  --output {output.bt2_agg}"
        "  results/Freyja/Demix/Results/"
        " ; "
        "freyja plot"
        "  {output.bt2_agg}"
        "  --output {output.bt2_png}"


rule Freyja_demix:
    input:
        bt2_variants = "results/Freyja/Variants/Results/{sample}.freyja.variants.tsv",
        bt2_depths = "results/Freyja/Variants/Results/{sample}.freyja.depths.tsv"
    output: "results/Freyja/Demix/Results/{sample}_freyja.demix"
    threads: 8
    resources:
        mem_mb=32000
    conda: "envs/freyja.yml"
    message: "{wildcards.sample}: Running Freyja Demix."
    shell:
        "freyja demix"
        "  {input.bt2_variants}"
        "  {input.bt2_depths}"
        "  --output {output}"
        # "  --confirmedonly"


rule Freyja_variants:
    input:
        updated = "results/Freyja/Update/update.ok",
        bt2_bam = "results/FastQs/{sample}/{sample}.bam",
        dna = FASTA_REFERENCE
    output:
        bt2_variants = "results/Freyja/Variants/Results/{sample}.freyja.variants.tsv",
        bt2_depths = "results/Freyja/Variants/Results/{sample}.freyja.depths.tsv"
    threads: 8
    resources:
        mem_mb=32000
    conda: "envs/freyja.yml"
    message: "{wildcards.sample}: Running Freyja Variants."
    shell:
        "freyja variants"
        "  {input.bt2_bam}"
        "  --variants {output.bt2_variants}"
        "  --depths {output.bt2_depths}"
        "  --ref {input.dna}"


rule Freyja_update:
    output: "results/Freyja/Update/update.ok"
    threads: 8
    resources:
        mem_mb=16000
    conda: "envs/freyja.yml"
    message: "Updating Freyja database."
    shell:
        "freyja update ; touch {output}"


def qualimap_of_valid_samples(wildcards):
    existence_path = checkpoints.check_file_existence.get(**wildcards).output[0]
    existence_df = pd.read_table(existence_path)
    present_samples = existence_df.sample_id[existence_df.sample_present == "ok"]
    return {
        "cov_across_ref": [
            f"results/BAM/{sample}/{sample}.qualimap/raw_data_qualimapReport/coverage_across_reference.txt"
            for sample in present_samples
        ],
        "genome_fraction_cov": [
            f"results/BAM/{sample}/{sample}.qualimap/raw_data_qualimapReport/genome_fraction_coverage.txt"
            for sample in present_samples
        ]
    }

checkpoint coverage_summary:
    input: unpack(qualimap_of_valid_samples)
    output:
        long_form_report = "results/Coverage/coverageReport.tsv",
        cov_per_sample = "results/Coverage/sample_coverage_status.tsv"
    message: "Generating coverage report of samples."
    script: "scripts/coverage_summary.py"


rule qualimap:
    input:
        bam = "results/FastQs/{sample}/{sample}.bam"
    output:
        cov_across_ref = "results/BAM/{sample}/{sample}.qualimap/raw_data_qualimapReport/coverage_across_reference.txt",
        genome_fraction_cov = "results/BAM/{sample}/{sample}.qualimap/raw_data_qualimapReport/genome_fraction_coverage.txt"
    log: "results/BAM/{sample}/qualimap.log"
    threads: 8
    resources:
        mem_mb=32000
    conda: "envs/bioinfo.yml"
    message: "{wildcards.sample}: Evaluating BAM mapping quality with QualiMap."
    shell:
        "qualimap bamqc "
        "  -bam {input.bam} "
        "  --java-mem-size=32G "
        "  -nt {threads} "
        "  -outdir results/BAM/{wildcards.sample}/{wildcards.sample}.qualimap"
        "  > {log} 2>&1"


rule fastqc:
    input: "results/FastQs/{sample}/{sample}.fq.gz"
    output: "results/fastQC/{sample}_fastqc.html"
    log: "results/fastQC/{sample}.fastQC_se.log"
    conda: "envs/bioinfo.yml"
    message: "{wildcards.sample}: Checking read quality with fastqc."
    shell: "fastqc -o results/fastQC {input} > {log} 2>&1"


rule sam_to_fastq:
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
