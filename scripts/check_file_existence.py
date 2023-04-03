#!/usr/bin/env python
# coding: utf-8

import json
import datetime

import pandas as pd


manifesto_path = snakemake.input["manifesto"]
sample_metadata_path = snakemake.input["metadata"]
outfile = snakemake.output[0]


with open(manifesto_path) as f:
    manifesto = json.load(f)

sample_metadata = (
    pd.read_csv(sample_metadata_path)
    .assign(sample_collect_date = lambda df: pd.to_datetime(df.sample_collect_date))
)

downloaded_sample_ids = set(manifesto.keys())
post2023_sampleids = set(sample_metadata.sample_id[sample_metadata.sample_collect_date >= pd.Timestamp(2022, 12, 28)])

status = {
    "no_metadata": downloaded_sample_ids - post2023_sampleids,
    "no_file": post2023_sampleids - downloaded_sample_ids,
    "file_present": downloaded_sample_ids.intersection(post2023_sampleids)
}

status_df = pd.concat([
    pd.DataFrame({"existence": status_str}, index=list(sample_ids))
    for status_str, sample_ids in status.items()
])
status_df.index.name = "sample_id"

status_df.to_csv(outfile, sep="\t")
