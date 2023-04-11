import pandas as pd

existence_path = snakemake.input["existence"]
coverage_path = snakemake.input["coverage"]
freyja_path = snakemake.input["freyja"]

outfile = snakemake.output[0]

def get_sample_status(row):
    if row.sample_present != "ok":
        return "PROBLEM:id_mismatch"
    elif row.enough_coverage != "yes":
        return "PROBLEM:low_coverage"
    elif (row.freyja_num_variants > 1) and (row.freyja_num_unique_freqs == 1):
        return "PROBLEM:freyja_no_freqs"
    else:
        return "ok"

result = (
    pd.concat(
        [
            pd.read_table(path).set_index("sample_id")
            for path in (existence_path, coverage_path, freyja_path)
        ],
        axis="columns",
        join="outer"
    )
    .assign(sample_status=lambda df: [get_sample_status(df.loc[ix]) for ix in df.index])
)

result.to_csv(outfile, sep="\t")
