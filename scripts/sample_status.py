from pathlib import Path
import pandas as pd

existence_path = snakemake.input["existence"]
coverage_path = snakemake.input["coverage"]
freyja_demix_path = snakemake.input["freyja_demix"]
freyja_quality_path = snakemake.input["freyja_quality"]

outfile = snakemake.output[0]

def get_sample_status(row):
    if row.sample_present != "ok":
        return row.sample_present
    elif row.enough_coverage != "yes":
        return "low_coverage"
    elif row.freyja_demix_status == "solver_error":
        return "freyja_demix_error"
    elif (row.freyja_num_variants > 1) and (row.freyja_num_unique_freqs == 1):
        return "freyja_only_one_freq"
    else:
        return "ok"

result = (
    pd.concat(
        [
            pd.read_table(path).set_index("sample_id")
            for path in (existence_path, coverage_path, freyja_demix_path, freyja_quality_path)
        ],
        axis="columns",
        join="outer"
    )
    .assign(sample_status=lambda df: [get_sample_status(df.loc[ix]) for ix in df.index])
)

result.to_csv(outfile, sep="\t")
