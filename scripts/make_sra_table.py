
from pathlib import Path
import pandas as pd

bams = [Path(p).stem for p in snakemake.input]
outfile = snakemake.output[0]

df = pd.DataFrame({
    "sample_name": bams,
    "library_ID": bams,
    "title": bams,
    "library_strategy": "AMPLICON",
    "library_source": "VIRAL RNA",
    "library_selection": "RT-PCR",
    "library_layout": "single",
    "platform": "ION_TORRENT",
    "instrument_model": "Ion Torrent Genexus",
    "design_description": "methods included in custom attributes",
    "filetype": "bam",
    "filename": [f"{sample}.bam" for sample in bams],
    "filename2": None,
    "filename3": None,
    "filename4": None,
    "assembly": "NC_045512.2",
    "fasta_file": None,
    "enrichment_kit": None,
    "amplicon_PCR_primer_scheme": None,
    "library_preparation_kit": None,
    "amplicon_size": None,
    "quality_control_method": None,
    "quality_control_method_version": None,
    "quality_control_determination": None,
    "quality_control_issues": None,
    "quality_control_details": None,
    "dehosting_method": None,
    "sequence_submitter_contact_email": None,
    "raw_sequence_data_processing_method": None,
})

df.to_csv(outfile, index=False)
