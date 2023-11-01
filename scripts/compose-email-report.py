
from datetime import date, datetime
from pathlib import Path

import pandas as pd


newest_table = snakemake.input["newest"]
previous_table = snakemake.input["previous"]
callout_group_threshold = snakemake.params["threshold"]
outfile = snakemake.output[0]


newest = pd.read_table(newest_table)
previous = pd.read_table(previous_table)


# What is different?

new_records = newest.loc[
    [sid not in previous.sample_id.values for sid in newest.sample_id]
]


# Are there any new samples at or above the threshold for the callout group?

frequencies = (
    new_records
    .groupby("sample_id")
    .variant_pct
    .sum()
    .sort_values(ascending=False)
)

reported_ids = frequencies.loc[frequencies > callout_group_threshold].index.tolist()


# Compose e-mail report.

def compose_empty_report(threshold, previous_date):
    return f"Since the last report on {previous_date}, there have been no new samples analyzed with frequencies of BA.2.86 above {threshold:.1%}."

def compose_samples_report(reported_ids, df, threshold, previous_date):
    report = f"Since the last report on {previous_date}, the following samples have been analyzed and found to have frequencies of the BA.2.86 group above {threshold:.1%}:\n\n"
    sample_reports = list()
    for sid in reported_ids:
        variant_info = df.loc[df.sample_id == sid]
        sample_info = variant_info.iloc[0]
        collection_date = date.fromisoformat(sample_info.sample_collect_date).strftime("%B %d, %Y")
        total = variant_info.variant_pct.sum()
        variants = [
            f"{row.variant} ({row.variant_pct:.1%})"
            for ix, row in variant_info.iterrows()
        ]
        lines = [
            f"{sample_info.county} county, {sample_info.sewershed}",
            collection_date,
            f"Abundance: {total:.1%}",
            "Lineages: " + ", ".join(variants)
        ]
        sample_report = "\n".join(lines)
        sample_reports.append(sample_report)
    report += "\n\n".join(sample_reports) + "\n\n"
    report += f"All other new samples with BA.2.86 variants had frequencies below {threshold:.1%}."
    return report      

template = (
    "This e-mail is a report of the newest appearances of COVID variants belonging to the BA.2.86 group "
    "in the wastewater genetic surveillance sequences, as analyzed by Freyja. "
    "It was generated automatically upon completion of the pipeline on {today}."
    "\n"
    "\n"
    "{report}"
    "\n"
    "\n"
    "The complete unfiltered table of variant frequencies is attached."
    "\n"
    "\n"
)

prev_date_str = date.fromisoformat(
    Path(previous_table).stem.split("_")[0]
).strftime("%B %d, %Y")

if len(reported_ids) == 0:
    report = compose_empty_report(callout_group_threshold, prev_date_str)
else:
    report = compose_samples_report(reported_ids, newest, callout_group_threshold, prev_date_str)

message = template.format(
    today=datetime.today().strftime("%B %d, %Y at %I:%M %p"),
    report=report
)

with open(outfile, "w") as f:
    f.write(message)
