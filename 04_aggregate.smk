
import pandas as pd

SAMPLE_EXISTENCE = Path("output/sample_info/file_existence.tsv")
FREYJA_DEMIX_REPORT = Path("output/freyja/demix_report.tsv")


rule all:
    input:
        freyja_parse = "output/results/freyja_parse.csv",
        freyja_parse_barcode = "output/results/freyja_parse_barcode.csv",
        sample_info = "output/results/sample_info.tsv",
        sra_table = "output/results/SRA_table.csv",
        comprehensive_table = "output/results/comprehensive_results_table.tsv",
        dashboard_data = "output/results/var.data_summary.rds",
        id_tracking = "output/results/sample-id-report.tsv"


rule id_tracking_report:
    input:
        sample_info = SAMPLE_EXISTENCE,
        concentration = "data/sample_metadata/nys-wws-sars2-concentration.csv",
        sample_tracking = "data/sample_metadata/sars2-sequencing-manifest.csv"
    output: "output/results/sample-id-report.tsv"
    message: "Creating sample ID report."
    script: "scripts/sample-id-report.py"


rule dashboard_results:
    input:
        sewershed = "data/sample_metadata/sewershed_metadata.csv",
        variants_of_concern = "data/sample_metadata/variants_of_concern.csv",
        lineage_map = "data/sample_metadata/lineage_info.csv",
        freyja = "output/results/freyja_parse.csv",
        concentration = "data/sample_metadata/nys-wws-sars2-concentration.csv"
    output:
        rds_data = "output/results/var.data_summary.rds"
    message: "Creating results table for the dashboard."
    conda: "envs/tidyverse.yml"
    script: "scripts/genetic-sequencing-data-prep.R"


rule summary_table:
    input:
        freyja = "output/results/freyja_parse.csv",
        lineages = "data/sample_metadata/lineage_info.csv",
        concentration= "data/sample_metadata/nys-wws-sars2-concentration.csv",
        sewersheds = "data/sample_metadata/sewershed_metadata.csv",
        samples = "output/results/sample_info.tsv"
    output:
        table = "output/results/comprehensive_results_table.tsv"
    script: "scripts/comprehensive-summary-table.py"


def bam_links_of_covid_filtered(wildcards):
    valid_samples = (
        pd.read_table(SAMPLE_EXISTENCE)
        .query("sample_present == 'ok'")
        .sample_id
    )
    return expand(
        "output/covid-filtered-BAMs/{sample}.bam",
        sample=valid_samples
    )


rule SRA_table:
    input: bam_links_of_covid_filtered
    output: "output/results/SRA_table.csv"
    script: "scripts/make_sra_table.py"


rule symlink_processed_bam:
    input: "output/covid-filtered/{sample}/{sample}.bam"
    output: "output/covid-filtered-BAMs/{sample}.bam"
    message: "{wildcards.sample}: Creating symbolic link to processed BAM file."
    shell: "ln -s $(pwd)/{input} {output}"


rule sample_info_report:
    input:
        existence = "output/sample_info/file_existence.tsv",
        coverage = "output/qc/coverage/sample_coverage_status.tsv",
        freyja_demix = "output/freyja/demix_report.tsv",
        freyja_quality = "output/freyja/freyja_quality_report.tsv"
    output: "output/results/sample_info.tsv"
    script: "scripts/sample_status.py"


rule freyja_parse_barcode:
    input:
        report = "output/results/freyja_parse.csv",
        barcode_timestamp = "output/freyja/barcodes/last_barcode_update.txt"
    output: "output/results/freyja_parse_barcode.csv"
    shell:
        "echo \"#freyja_barcode_version\" $(cat {input.barcode_timestamp}) > $TMPDIR/freyja_barcode_header ; "
        "cat $TMPDIR/freyja_barcode_header {input.report} > {output}"


rule freyja_parse:
    input: "output/freyja/aggregated_results.tsv"
    output:
        longform = "output/results/freyja_parse.csv",
        report = "output/freyja/freyja_quality_report.tsv"
    message: "Parsing Freyja output."
    script: "scripts/parse_freyja_demix.py"


def freyja_demixed_samples(wildcards):
    valid_samples = (
        pd.read_table(FREYJA_DEMIX_REPORT)
        .query("freyja_demix_status != 'solver_error'")
        .sample_id
    )
    return expand(
        "output/freyja/demix/{sample}.demix",
        sample=valid_samples
    )


rule freyja_aggregate:
    input: freyja_demixed_samples
    output: "output/freyja/aggregated_results.tsv",
    conda: "envs/freyja.yml"
    message: "Running Freyja aggregate."
    shell:
        "freyja aggregate"
        "  --output {output}"
        "  output/freyja/demix/"
