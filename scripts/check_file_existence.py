#!/usr/bin/env python
# coding: utf-8

import re
import json
import datetime
from pathlib import Path

import pandas as pd


infolder = Path(snakemake.input["bam_folder"])
sample_metadata_path = snakemake.input["metadata"]
corrupt_files_path = Path(snakemake.input["corrupt_files"])
outfile = snakemake.output[0]

bam_paths = list(infolder.glob("**/*.ptrim.bam"))

sample_metadata = (
    pd.read_csv(sample_metadata_path)
    .assign(sample_collect_date = lambda df: pd.to_datetime(df.sample_collect_date))
)
concentration_sampleids = set(sample_metadata.sample_id[sample_metadata.sample_collect_date >= pd.Timestamp(2022, 12, 28)])

with open(corrupt_files_path, "r") as f:
    corrupt_files = set([Path(line.strip()).stem.split(".")[0] for line in f])

SEQSITE_IDS = {
    "buffalo": "suny_buffalo",
    "wadsworth": "wadsworth",
    "suny_upstate": "suny_upstate",
    "rochester": "u_of_rochester",
    "nymc": "nymc",
}

id_format = re.compile(r"202\d[01]\d[0123]\dNY\d\d\d\d\d\d\d\d\d[A-Z]")

records = []

for bam_path in bam_paths:
    sample_id = bam_path.stem.split(".")[0]
    seqsite = SEQSITE_IDS.get(bam_path.parts[2])
    if sample_id in corrupt_files:
        status = "corrupt_bam"
    elif id_format.fullmatch(sample_id) is None:
        status = "unexpected_bam_name"
    elif int(sample_id[:8]) < 20221228:
        status = "old_sample_bam"
    elif sample_id in concentration_sampleids:
        status = "ok"
    else:
        status = "bam_without_concentration"
    records.append((sample_id, status, seqsite))
    
bam_sampleids = set([bam.stem.split(".")[0] for bam in bam_paths])
for sample_id in concentration_sampleids - bam_sampleids:
    records.append((sample_id, "concentration_without_bam"))
    
status_df = pd.DataFrame.from_records(records, columns=["sample_id", "sample_present", "seq_lab_id"])
status_df.loc[status_df.sample_id.duplicated(keep=False), "sample_present"] = "duplicate_bam"
status_df.sort_values(by=["seq_lab_id", "sample_present"]).to_csv(outfile, sep="\t", index=False)