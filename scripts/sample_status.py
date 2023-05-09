from pathlib import Path
import pandas as pd

existence_path = snakemake.input["existence"]
coverage_path = snakemake.input["coverage"]
freyja_path = snakemake.input["freyja"]
freyja_status_paths = snakemake.input["freyja_status"]

outfile = snakemake.output[0]

freyja_failures = []
for status_path in freyja_status_paths:
    sample_id = Path(status_path).stem
    with open(status_path) as f:
        if f.readline().strip() != "ok":
            freyja_failures.append(sample_id)

def get_sample_status(row):
    if row.sample_present != "ok":
        return row.sample_present
    elif row.enough_coverage != "yes":
        return "low_coverage"
    elif (row.freyja_num_variants > 1) and (row.freyja_num_unique_freqs == 1):
        return "freyja_only_one_freq"
    elif row.name in freyja_failures:
        return "freyja_demix_failed"
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
