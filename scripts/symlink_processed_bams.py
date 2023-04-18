from pathlib import Path
import os

infiles = [Path(p) for p in snakemake.input]
outfolder = Path(snakemake.output[0])

outfolder.mkdir(exist_ok=True)

for infile in infiles:
    sample_id = infile.stem.replace(".bam", "")
    os.symlink(infile.resolve(), outfolder/f"{sample_id}.bam")
