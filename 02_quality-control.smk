
import pandas as pd

BAM_FOLDER = Path(config["bam_folder"])
BAM_HEADER = Path(config["bam_header"])
SAMPLE_EXISTENCE = Path("output/sample_info/file_existence.tsv")


def qualimap_of_bams_with_valid_ids(wildcards):
    valid_samples = (
        pd.read_table(SAMPLE_EXISTENCE)
        .query("sample_present == 'ok'")
        .sample_id
    )
    return {
        "cov_across_ref": [
            f"output/qc/qualimap/{sample}/raw_data_qualimapReport/coverage_across_reference.txt"
            for sample in valid_samples
        ],
        "genome_fraction_cov": [
            f"output/qc/qualimap/{sample}/raw_data_qualimapReport/genome_fraction_coverage.txt"
            for sample in valid_samples
        ]
    }


rule coverage_summary:
    input: unpack(qualimap_of_bams_with_valid_ids)
    output:
        long_form_report = "output/qc/coverage/coverage_report.tsv",
        cov_per_sample = "output/qc/coverage/sample_coverage_status.tsv"
    message: "Generating coverage report of samples."
    script: "scripts/coverage_summary.py"


rule qualimap:
    input:
        bam = "output/covid-filtered/{sample}/{sample}.bam"
    output:
        cov_across_ref = "output/qc/qualimap/{sample}/raw_data_qualimapReport/coverage_across_reference.txt",
        genome_fraction_cov = "output/qc/qualimap/{sample}/raw_data_qualimapReport/genome_fraction_coverage.txt"
    params:
        empty_across = "data/qualimap_empty_results/coverage_across_reference.txt",
        empty_fraction = "data/qualimap_empty_results/genome_fraction_coverage.txt",
    log: "output/qc/qualimap/{sample}/qualimap.log"
    threads: 4
    resources:
        mem_mb=32000
    conda: "envs/bioinfo.yml"
    message: "{wildcards.sample}: Evaluating BAM mapping quality with QualiMap."
    shell:
        "set +e ; "
        "qualimap bamqc "
        "  -bam {input.bam} "
        "  --java-mem-size=32G "
        "  -nt {threads} "
        "  -outdir output/qc/qualimap/{wildcards.sample}"
        "  > {log} 2>&1 ; "
        # If qualiMap failed, copy its log to the output files.
        # It can fail because the BAM file is empty.
        "RESULT=$? ; "
        "if [ $RESULT -ne 0 ]; then "
        "  cp {params.empty_across} {output.cov_across_ref} ; "
        "  cp {params.empty_fraction} {output.genome_fraction_cov} ; "
        "fi"


rule fastqc:
    input: "output/covid-filtered/{sample}/{sample}.fq.gz"
    output: "output/qc/fastqc/{sample}_fastqc.html"
    log: "output/qc/fastqc/{sample}.fastQC_se.log"
    conda: "envs/bioinfo.yml"
    message: "{wildcards.sample}: Checking read quality with fastqc."
    shell: "fastqc -o output/qc/fastqc {input} > {log} 2>&1"


rule sam_to_fastq:
    input: "output/covid-filtered/{sample}/{sample}.sam"
    output: "output/covid-filtered/{sample}/{sample}.fq.gz"
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
        bam = "output/covid-filtered/{sample}/{sample}.bam",
        sam = "output/covid-filtered/{sample}/{sample}.sam",
        stats_original = "output/covid-filtered/{sample}/{sample}_pre-filtered-stats.txt",
        stats_filtered = "output/covid-filtered/{sample}/{sample}_post-filtered-stats.txt"
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
