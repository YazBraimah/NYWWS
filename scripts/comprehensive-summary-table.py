
from pathlib import Path

import pandas as pd

freyja_file = snakemake.input["freyja"]
concentration_file = snakemake.input["concentration"]
lineage_file = snakemake.input["lineages"]
sewersheds_file = snakemake.input["sewersheds"]
samples_file = snakemake.input["samples"]
outfile = snakemake.output["table"]

freyja = (
    pd.read_csv(freyja_file)
    .set_index("sample_id")
)

concentration = (
    pd.read_csv(concentration_file)
    .set_index("sample_id")
    [["sample_collect_date", "sw_id"]]
)

sewersheds = (
    pd.read_csv(sewersheds_file)
    .set_index("sw_id")
    [["county", "wwtp_name"]]
)

sample_info = (
    pd.read_table(samples_file)
    .set_index("sample_id")
    [["seq_lab_id"]]
)

lineage_map = dict()
lineage_colors = dict()
for ix, row in pd.read_csv(lineage_file).iterrows():
    lineage_map[row.lineage] = row.callout_group
    lineage_colors[row.lineage] = row.hex_code

df = (
    freyja
    .join(concentration, how="left")
    .join(sample_info, how="left")
    .join(sewersheds, on="sw_id", how="left")
    .assign(
        callout_group = lambda df: [lineage_map.get(var) for var in df.variant],
        hex_color = lambda df: [lineage_colors.get(var) for var in df.variant],
    )
    .drop(columns="sw_id")
    .rename(columns={
        "seq_lab_id": "seq_lab",
        "wwtp_name": "sewershed",
    })
    .reset_index()
    [["sample_id", "variant", "callout_group", "hex_color", "variant_pct",
      "sample_collect_date", "county", "sewershed", "seq_lab"]]
    .sort_values(
        by=["sample_id", "callout_group", "variant_pct"],
        ascending=[True, True, False]
    )
)

df.to_csv(outfile, sep="\t", index=False)
