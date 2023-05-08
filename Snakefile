from pathlib import Path

import pandas as pd


BAM_FOLDER = Path(config["bam_folder"])
SAMPLE_METADATA = Path(config["sample_metadata"])
BAM_HEADER = Path(config["bam_header"])
FASTA_REFERENCE = Path(config["fasta_ref"])


rule all:
    input:
        "results/Summary/sample_info.tsv",
        "results/Summary/comprehensive_results_table.txt",
        "results/Summary/SRA_table.csv",
        "results/Freyja/Aggregate/freyja_parse_barcode.csv",
        "results/Summary/var.data_summary.rds"


rule prep_freyja_results:
    input:
        sewershed = "data/sample_metadata/sewershed_metadata.csv",
        variants_of_concern = "data/sample_metadata/variants_of_concern.csv",
        lineage_map = "data/sample_metadata/lineage_info.csv",
        freyja = "results/Freyja/Aggregate/freyja_parse.csv"
    output:
        rds_data = "results/Summary/var.data_summary.rds"
    conda: "envs/tidyverse.yml"
    script: "scripts/genetic-sequencing-data-prep.R"


rule summary_plots:
    input:
        cov = "results/Coverage/coverageReport.tsv",
        parse = "results/Freyja/Aggregate/freyja_parse.csv",
        samInfo = "results/Summary/sample_info.tsv"
    output:
        coveragePointPlot = "results/Summary/coveragePointPlot_by_date.pdf",
        varFreqBarPlotsMin20X = "results/Summary/Variant_frequency_barplots_by_county_min_20X.pdf",
        compResultsTable = "results/Summary/comprehensive_results_table.txt"
    params:
        samLoc = "data/sample_metadata/sampling_locations.csv",
        linInfo = "data/sample_metadata/lineage_info.csv",
        r_script = "scripts/results_summary_outputs.R"
    threads: 8
    resources:
        mem_mb=32000
    message: "Outputting summary tables and plots."
    conda: "envs/tidyverse.yml"
    shell: "Rscript {params.r_script} {input.parse} {input.cov} {params.samLoc} {input.samInfo} {params.linInfo} {output.coveragePointPlot} {output.varFreqBarPlotsMin20X} {output.compResultsTable}"


def processed_bam_links(wildcards):
    existence_path = checkpoints.check_file_existence.get(**wildcards).output[0]
    existence_df = pd.read_table(existence_path)
    present_samples = existence_df.sample_id[existence_df.sample_present == "ok"]
    return expand("output/BAM_processed/{sample}.bam", sample=present_samples)

rule SRA_table:
    input: processed_bam_links
    output: "results/Summary/SRA_table.csv"
    script: "scripts/make_sra_table.py"


rule symlink_processed_bam:
    input: "results/FastQs/{sample}/{sample}.bam"
    output: "output/BAM_processed/{sample}.bam"
    message: "{wildcards.sample}: Creating symbolic link to processed BAM file."
    shell: "ln -s $(pwd)/{input} {output}"


rule sample_info_report:
    input:
        existence = "output/sample_info/file_existence.tsv",
        coverage = "results/Coverage/sample_coverage_status.tsv",
        freyja = "output/sample_info/freyja_quality.tsv"
    output: "results/Summary/sample_info.tsv"
    script: "scripts/sample_status.py"


rule Freyja_parse_barcode:
    input: "results/Freyja/Aggregate/freyja_parse.csv"
    output:
        barcode = "results/Freyja/Aggregate/barcode.txt",
        csv = "results/Freyja/Aggregate/freyja_parse_barcode.csv"
    conda: "envs/freyja.yml"
    shell:
        "echo \"#freyja_barcode_version\" $(freyja demix --version | sed '2q;d') > {output.barcode} ; "
        "cat {output.barcode} {input} > {output.csv}"


rule Freyja_parse:
    input: "results/Freyja/Aggregate/aggregated_results.tsv"
    output:
        longform = "results/Freyja/Aggregate/freyja_parse.csv",
        report = "output/sample_info/freyja_quality.tsv"
    threads: 2
    resources:
        mem_mb=4000
    message: "Parsing Freyja output."
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
    threads: 8
    resources:
        mem_mb=32000
    conda: "envs/freyja.yml"
    message: "Running Freyja Aggregate."
    shell:
        "freyja aggregate"
        "  --output {output.bt2_agg}"
        "  results/Freyja/Demix/Results/"


rule Freyja_demix:
    input:
        updated = "results/Freyja/Update/update.ok",
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
        "  --confirmedonly"


rule Freyja_variants:
    input:
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


def find_raw_bam(wildcards):
    sample = wildcards.sample
    search = list(BAM_FOLDER.glob(f"**/{sample}.ptrim.bam"))
    assert len(search) == 1
    return str(search[0])

rule filter_covid_reads:
    input:
        bam = find_raw_bam,
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
        "samtools view -b --threads {threads} {input.bam} 2019-nCoV "
        "| samtools reheader {input.header} - "
        "| samtools sort -o {output.bam} -@ {threads} ; " 
        "samtools view {output.bam} -o {output.sam} ; "
        "samtools stats {input.bam} | grep ^SN | cut -f 2- > {output.stats_original} ; "
        "samtools stats {output.bam} | grep ^SN | cut -f 2- > {output.stats_filtered} ; "


checkpoint check_file_existence:
    input:
        bam_folder = BAM_FOLDER,
        metadata = SAMPLE_METADATA
    output: "output/sample_info/file_existence.tsv"
    message: "Matching downloaded BAMs with known sample IDs."
    script: "scripts/check_file_existence.py"
