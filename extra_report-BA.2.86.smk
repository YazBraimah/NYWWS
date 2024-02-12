
# This Snakefile would create an e-mail report of BA.2.86 appearances

#     snakemake \
# 	--snakefile 05_report-BA.2.86.smk \
# 	-c1 \
# 	--use-conda \
# 	--configfile ${PIPELINE_CONFIG} \
# 	--rerun-triggers mtime \
# 	--config date=$(date +"%Y%m%d")

from datetime import date
from pathlib import Path

import yaml


TODAY = config["date"]
REPORT_THRESHOLD = config["BA286_report_threshold"]

with open("data/email-config/recipients.yml") as f:
    recipients = yaml.load(f, Loader=yaml.Loader)

CARBON_COPIES = recipients["email_cc"]
RECIPIENTS = recipients["email_to"]


rule all:
    input: f"output/email-reports/sent/{TODAY}.txt"


rule send_email:
    input:
        muttrc = "data/email-config/muttrc",
        msmtprc = "data/email-config/.msmtprc",
        text = "output/email-reports/text/{date}.txt",
        table = "output/results/variant-tracking-reports/BA.2.86/{date}_BA.2.86_freyja.tsv"
    output: "output/email-reports/sent/{date}.txt"
    params:
        carbon_copies = lambda wc: " ".join(["-c " + addr for addr in CARBON_COPIES]),
        recipients = lambda wc: " ".join([addr for addr in RECIPIENTS])
    shell:
        "cat {input.text} | mutt"
        " -F {input.muttrc}"
        " -s \"Newest BA.2.86 appearances\""
        " -a {input.table}"
        " {params.carbon_copies}"
        " -- {params.recipients}"
        " ; "
        "cp {input.text} {output} ; "


def previous_report(wc):
    reports = (
        Path("output/results/variant-tracking-reports/BA.2.86/")
        .glob("[!.]*_BA.2.86_freyja.tsv")
    )
    before = [
        filename
        for filename in reports
        if filename.stem.split("_")[0] < wc["date"]
    ]
    return sorted(before)[-1]

rule compose_email:
    input:
        newest = "output/results/variant-tracking-reports/BA.2.86/{date}_BA.2.86_freyja.tsv",
        previous = previous_report
    params:
        threshold = REPORT_THRESHOLD
    output: "output/email-reports/text/{date}.txt"
    script: "scripts/compose-email-report.py"


rule make_table:
    output: "output/results/variant-tracking-reports/BA.2.86/{date}_BA.2.86_freyja.tsv",
    params:
        # This is not an input because it changes every time the pipeline runs
        # and we don't want those changes to trigger an override of past BA.2.86
        # reports
        table = "output/results/comprehensive_results_table.tsv"
    shell:
        "grep \"BA.2.86\" {params.table} | sort -k 1 | sort -r -k 5 > temp_table ; "
        "head -n 1 {params.table} > temp_header ; "
        "cat temp_header temp_table > {output} ; "
        "rm temp_header temp_table"
