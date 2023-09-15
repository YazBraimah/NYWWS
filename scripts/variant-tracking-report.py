from pathlib import Path

import numpy as np
import pandas as pd


CONC_DATA = snakemake.input["concentration"]
SEWER_DATA = snakemake.input["sewershed"]
FREYJA_DATA = snakemake.input["freyja"]
COVERAGE = snakemake.input["coverage"]

VARIANT = snakemake.params["variant"]
MIN_READS = snakemake.params["min_reads"]

OUTFILE = snakemake.output[0]

EVIDENCE_NAME = "logsum"


def get_record(read_count_file):
    sample_id = Path(read_count_file).stem
    with open(read_count_file, "r") as f:
        read_count = int(f.readline().strip())
    return (sample_id, read_count)


def get_records(file_list):
    records = [get_record(read_count_file) for read_count_file in file_list]
    return records


def get_read_counts(file_list, read_count_col, min_counts=None):
    records = get_records(file_list)
    df = (
        pd.DataFrame
        .from_records(records, columns=["sample_id", read_count_col])
        .set_index("sample_id")
    )
    if min_counts is not None:
        df = df.query(f"{read_count_col} >= {min_counts}")
    return df


def measure_evidence(df):
    counts = df[[col for col in df.columns if col.startswith("reads_")]]
    evidence = np.log1p(counts).sum(axis=1)
    return evidence


def get_read_counts_df():
    mut_groups = [
        s
        for s in snakemake.input.keys()
        if (s.startswith("counts:")) and (s != "counts:total")
    ]
    df = pd.concat(
        [
            get_read_counts(
                snakemake.input[group],
                "reads_" + group.split(":")[1],
                min_counts=MIN_READS
            )
            for group in sorted(mut_groups)
        ],
        axis="columns",
        join="outer"
    )
    evidence_metric = measure_evidence(df)
    df[EVIDENCE_NAME] = evidence_metric
    total_counts = get_read_counts(snakemake.input["counts:total"], "reads_total")
    df = df.join(total_counts, how="left")
    return df


if __name__ == "__main__":
    concentration_cols = ["sample_collect_date", "sw_id"]
    sewershed_cols = ["sw_id", "wwtp_name", "county", "zipcode"]

    concentration = pd.read_csv(CONC_DATA).set_index("sample_id")[concentration_cols]
    sewersheds = pd.read_csv(SEWER_DATA)[sewershed_cols].set_index("sw_id")

    coverage = (
        pd.read_table(COVERAGE)[["sample_id", "Coverage"]]
        .groupby("sample_id")
        .mean()
        .rename({"Coverage": "avg_coverage"}, axis="columns")
    )

    freyja_results = (
        pd.read_csv(FREYJA_DATA)
        .set_index("sample_id")
    )
    freyja_results = (
        freyja_results.loc[
            (freyja_results.variant == VARIANT) | (freyja_results.variant.str.startswith(VARIANT + ".")),
            ["variant", "variant_pct"]
        ]
        .rename({"variant": "freyja_variant", "variant_pct": f"freyja_abundance"}, axis="columns")
    )

    df = (
        get_read_counts_df()
        .join(coverage, how="left")
        .join(freyja_results, how="left")
        .join(concentration, how="left")
        .join(sewersheds, on="sw_id", how="left")
        .sort_values(by=EVIDENCE_NAME, ascending=False)
    )

    df.to_csv(OUTFILE, sep="\t", index=True)
