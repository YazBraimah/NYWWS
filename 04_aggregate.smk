
import pandas as pd

SAMPLE_EXISTENCE = Path("output/sample_info/file_existence.tsv")


rule all:
    input:
        freyja_parse = "output/results/freyja_parse.csv",
        freyja_parse_barcode = "output/results/freyja_parse_barcode.csv",
        sample_info = "output/results/sample_info.tsv",
        sra_table = "output/results/SRA_table.csv",
        comprehensive_table = "output/results/comprehensive_results_table.txt",
        dashboard_data = "output/results/var.data_summary.rds"


rule dashboard_results:
    input:
        sewershed = "data/sample_metadata/sewershed_metadata.csv",
        variants_of_concern = "data/sample_metadata/variants_of_concern.csv",
        lineage_map = "data/sample_metadata/lineage_info.csv",
        freyja = "output/results/freyja_parse.csv"
    output:
        rds_data = "output/results/var.data_summary.rds"
    message: "Creating results table for the dashboard."
    conda: "envs/tidyverse.yml"
    script: "scripts/genetic-sequencing-data-prep.R"


rule summary_plots:
    input:
        cov = "output/qc/coverage/coverage_report.tsv",
        parse = "output/results/freyja_parse.csv",
        samInfo = "output/results/sample_info.tsv",
        samLoc = "data/sample_metadata/sampling_locations.csv",
        linInfo = "data/sample_metadata/lineage_info.csv"
    output:
        coveragePointPlot = "output/plots/coveragePointPlot_by_date.pdf",
        varFreqBarPlotsMin20X = "output/plots/Variant_frequency_barplots_by_county_min_20X.pdf",
        compResultsTable = "output/results/comprehensive_results_table.txt"
    message: "Outputting summary tables and plots."
    conda: "envs/tidyverse.yml"
    script: "scripts/results_summary_outputs.R"


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
        "cat $TMPDIR/freyja_barcode_header {input} > {output}"


rule freyja_parse:
    input: "output/freyja/aggregated_results.tsv"
    output:
        longform = "output/results/freyja_parse.csv",
        report = "output/freyja/freyja_quality_report.tsv"
    message: "Parsing Freyja output."
    script: "scripts/parse_freyja_demix.py"


rule freyja_aggregate:
    input: freyja_demixed_samples
    output: "output/freyja/aggregated_results.tsv",
    threads: 8
    resources:
        mem_mb=32000
    conda: "envs/freyja.yml"
    message: "Running Freyja aggregate."
    shell:
        "freyja aggregate"
        "  --output {output}"
        "  output/freyja/demix/"
