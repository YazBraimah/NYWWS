#!/usr/bin/env python
# coding: utf-8

import pandas as pd


sample_info_path = snakemake.input["sample_info"]
concentration_path = snakemake.input["concentration"]
sample_tracking_path = snakemake.input["sample_tracking"]

outfile = snakemake.output[0]


def make_id_status(series):
    key = {
        "ok": "BAM and concentration match",
        "no_metadata": "BAM not in concentration",
        "no_file": "ID in concentration but no BAM"
    }
    return [key[s] for s in series]


sample_info = (
    pd.read_table(sample_info_path)
    .assign(
        id_status = lambda df: make_id_status(df.sample_present)
    )
    [["sample_id", "id_status", "seq_lab_id"]]
    .set_index("sample_id", verify_integrity=True)
    .rename({"seq_lab_id": "seq_lab"}, axis="columns")
)

# TEMPORARILY, we remove duplicated sample IDs from this one.
# As of June 5, 2023, the following IDs are replicated after
# Dec. 28, 2022 (start of project):
# '20221229NY071026522A', '20221229NY071026310A', '20221229NY071031518A'
pcr_lab = (
    pd.read_csv(concentration_path)
    .assign(
        sample_collect_date = lambda df: pd.to_datetime(df.sample_collect_date)
    )
    .loc[lambda df: df.sample_collect_date >= "2022-12-28"]
    .loc[lambda df: ~df.sample_id.duplicated()] 
    [["sample_id", "pcr_lab"]]
    .set_index("sample_id", verify_integrity=True)
)

date = (
    pd.read_csv(sample_tracking_path)
    .assign(
        sample_id = lambda df: [s.split(".")[0] for s in df.filename]
    )
    [["sample_id", "created_date"]]
    .set_index("sample_id", verify_integrity=True)
    .rename({"created_date": "bam_upload_date"}, axis="columns")
)


id_report = (
    sample_info.join([pcr_lab, date], how="left")
    .sort_index()
    [["id_status", "pcr_lab", "seq_lab", "bam_upload_date"]]
)

id_report.to_csv(outfile, index=True, sep="\t")
