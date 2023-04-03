from pathlib import Path
import os

infolder = Path(snakemake.config["bam_folder"])
outfolder = Path(snakemake.output[0])

outfolder.mkdir(exist_ok=True)

for raw_bam in Path(infolder).glob("**/*.ptrim.bam"):
    sample_id = raw_bam.stem.replace(".ptrim", "")
    os.symlink(raw_bam.resolve(), outfolder/f"{sample_id}.bam")
