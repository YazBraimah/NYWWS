from pathlib import Path
import json

infolder = snakemake.input[0]
outfile = snakemake.output[0]

def path_to_sampleid(path):
    return path.stem.split(".")[0]

bam_files = Path(infolder).glob("**/*.ptrim.bam")

manifesto = {
    path_to_sampleid(path): {"trim": str(path.resolve())}
    for path in bam_files
}

with open(outfile, "w") as f:
    json.dump(manifesto, f, indent=2)
