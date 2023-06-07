
import pandas as pd

FASTA_REFERENCE = Path(config["fasta_ref"])
COVERAGE_REPORT = "output/qc/coverage/sample_coverage_status.tsv"


def freyja_demix_of_qc_passed_samples(wildcards):
    valid_samples = (
        pd.read_table(COVERAGE_REPORT)
        .query("enough_coverage == 'yes'")
        .sample_id
    )
    return expand(
        "output/freyja/demix/{sample}.log",
        sample=valid_samples
    )


rule freyja_demix_report:
    input:
        logs = freyja_demix_of_qc_passed_samples
    output: "output/freyja/demix_report.tsv"
    message: "Generating report of Freyja demix results."
    script: "scripts/freyja_demix_report.py"


rule freyja_demix:
    input:
        barcodes = "output/freyja/barcodes/usher_barcodes.csv",
        bt2_variants = "output/freyja/variants/{sample}_variants.tsv",
        bt2_depths = "output/freyja/variants/{sample}_depths.tsv"
    params:
        demix = "output/freyja/demix/{sample}.demix",
    log: "output/freyja/demix/{sample}.log"
    threads: 8
    resources:
        mem_mb=32000
    conda: "envs/freyja.yml"
    message: "{wildcards.sample}: Running Freyja demix."
    shell:
        "set +e ; "
        "freyja demix"
        "  {input.bt2_variants}"
        "  {input.bt2_depths}"
        "  --barcodes {input.barcodes}"
        "  --output {params.demix}"
        "  --confirmedonly"
        "&> {log} ; "
        "RESULT=$?"


rule freyja_variants:
    input:
        bt2_bam = "output/covid-filtered/{sample}/{sample}.bam",
        dna = FASTA_REFERENCE
    output:
        bt2_variants = "output/freyja/variants/{sample}_variants.tsv",
        bt2_depths = "output/freyja/variants/{sample}_depths.tsv"
    threads: 8
    resources:
        mem_mb=32000
    conda: "envs/freyja.yml"
    message: "{wildcards.sample}: Running Freyja variants."
    shell:
        "freyja variants"
        "  {input.bt2_bam}"
        "  --variants {output.bt2_variants}"
        "  --depths {output.bt2_depths}"
        "  --ref {input.dna}"


rule freyja_update:
    output:
        barcodes = "output/freyja/barcodes/usher_barcodes.csv",
        timestamp = "output/freyja/barcodes/last_barcode_update.txt"
    conda: "envs/freyja.yml"
    message: "Updating Freyja barcodes."
    shell:
        "freyja update --outdir output/freyja/barcodes"
