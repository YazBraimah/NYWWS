#!/usr/bin/env python
# coding: utf-8

import json
import datetime
from pathlib import Path

import pandas as pd


infolder = Path(snakemake.input["bam_folder"])
sample_metadata_path = snakemake.input["metadata"]
outfile = snakemake.output[0]

bam_paths = list(infolder.glob("**/*.ptrim.bam"))
sample_ids = [bam.stem.split(".")[0] for bam in bam_paths]

sample_metadata = (
    pd.read_csv(sample_metadata_path)
    .assign(sample_collect_date = lambda df: pd.to_datetime(df.sample_collect_date))
)

downloaded_sample_ids = set(sample_ids)
post2023_sampleids = set(sample_metadata.sample_id[sample_metadata.sample_collect_date >= pd.Timestamp(2022, 12, 28)])

status = {
    "no_metadata": downloaded_sample_ids - post2023_sampleids,
    "no_file": post2023_sampleids - downloaded_sample_ids,
    "ok": downloaded_sample_ids.intersection(post2023_sampleids)
}

status_df = pd.concat([
    pd.DataFrame({"sample_present": status_str}, index=list(sample_ids))
    for status_str, sample_ids in status.items()
])
status_df.index.name = "sample_id"

# Get sequencing sites from the folders

SEQSITE_IDS = {
    "buffalo": "suny_buffalo",
    "wadsworth": "wadsworth",
    "suny_upstate": "suny_upstate",
    "rochester": "u_of_rochester",
    "nymc": "nymc",
}

def get_seqsite_id(folder):
    return SEQSITE_IDS[folder]

seq_sites = pd.DataFrame({
    "sample_id": sample_ids,
    "seq_lab_id": [get_seqsite_id(bam.parts[2]) for bam in bam_paths]
}).set_index("sample_id")

# Temporary fix to remove duplicate uploads
seq_sites = seq_sites.loc[~seq_sites.index.duplicated()]

status_df = pd.concat([status_df, seq_sites], axis="columns")

# Save results

status_df.to_csv(outfile, sep="\t")
