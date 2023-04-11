
import pandas as pd

infile = snakemake.input[0]
outfile = snakemake.output[0]

df = (
    pd.read_csv(infile)
    .groupby("sample_id")
    .nunique()
    .rename({"variant": "freyja_num_variants", "frequency": "freyja_num_unique_freqs"}, axis="columns")
)

df.to_csv(outfile, sep="\t")
