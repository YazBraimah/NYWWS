#!/usr/bin/env python
# coding: utf-8

from pathlib import Path

import pandas as pd


cov_across_ref = snakemake.input["cov_across_ref"]
genome_fraction_cov = snakemake.input["genome_fraction_cov"]

outfile_across_ref_summary = snakemake.output["long_form_report"]
outfile_samples_cov = snakemake.output["cov_per_sample"]

minimum_cov = snakemake.config["minimum_coverage"]
minimum_fraction = snakemake.config["minimum_fraction"]


across_ref_summary = (
    pd.concat({
        Path(f).parts[-3]: pd.read_table(f)
        for f in cov_across_ref
    })
    .reset_index(level=0, names="sample_id")
    .reset_index(drop=True)
)
across_ref_summary.to_csv(outfile_across_ref_summary, sep="\t", index=False)


records = []
for f in genome_fraction_cov:
    sample_id = Path(f).parts[-3]
    fraction_covs = pd.read_table(f).set_index("#Coverage (X)")
    try:
        fraction = fraction_covs.loc[minimum_cov].Coverage/100
    except:
        fraction = 0
    if fraction >= minimum_fraction:
        enough_cov = "yes"
    else:
        enough_cov = "no"
    records.append((sample_id, fraction, enough_cov))

samples_cov = pd.DataFrame.from_records(
    records,
    columns=["sample_id", f"genome_fraction_{minimum_cov}x", "enough_coverage"]
)
samples_cov.to_csv(outfile_samples_cov, sep="\t", index=False)
