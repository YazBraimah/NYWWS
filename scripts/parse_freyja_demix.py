#!/usr/bin/env python
# coding: utf-8

from itertools import zip_longest

import pandas as pd


freyja_output = snakemake.input[0]
longform_file = snakemake.output["longform"]
report_file = snakemake.output["report"]


freyja_raw = (
    pd.read_table(freyja_output)
    .rename({"Unnamed: 0": "sample_id"}, axis="columns")
    .assign(sample_id=lambda df: df.sample_id.str.strip("_variants.tsv"))
    [["sample_id", "lineages", "abundances"]]
)
records = []
for _ix, row in freyja_raw.iterrows():
    lineages = row.lineages.split()
    abundances = [float(x) for x in row.abundances.split()]
    records.extend(zip_longest([row.sample_id], lineages, abundances, fillvalue=row.sample_id))
longform = pd.DataFrame.from_records(records, columns=["sample_id", "variant", "variant_pct"])


# Report number of variants found per sample and how many unique variant percentages
# Because bad samples get the same percentage for all variants
quality_report = (
    longform
    .groupby("sample_id")
    .nunique()
    .rename({"variant": "freyja_num_variants", "variant_pct": "freyja_num_unique_freqs"}, axis="columns")
)
quality_report.to_csv(report_file, sep="\t")


# Add an "Unidentified" variant so all samples sum to 1.0
freqs_required = 1 - longform.groupby("sample_id").variant_pct.sum()
unidentified_records = [
    (sample_id, "Unidentified", freq)
    for sample_id, freq in freqs_required.items()
    if freq > 0
]
unidentified_df = pd.DataFrame.from_records(unidentified_records, columns=["sample_id", "variant", "variant_pct"])
longform_all = (
    pd.concat([longform, unidentified_df])
    .sort_values(by=["sample_id", "variant_pct"], ascending=[True, False])
    .reset_index(drop=True)
)
longform_all.to_csv(longform_file, index=False)
