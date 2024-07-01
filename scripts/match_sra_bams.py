from pathlib import Path

import pandas as pd

sra_path = snakemake.input[0]
bam_path = snakemake.params["bam_path"]

sra = pd.read_csv(sra_path)
sra_samples = set(sra.sample_name)
bam_samples = set(p.stem for p in Path(bam_path).glob("*.bam"))
obsolete_bams = bam_samples - sra_samples

for bam in obsolete_bams:
    path = Path(bam_path) / f"{bam}.bam"
    path.unlink(missing_ok=False)