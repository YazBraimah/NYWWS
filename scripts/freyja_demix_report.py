
from pathlib import Path
import pandas as pd

freyja_logs = snakemake.input["logs"]
outfile = snakemake.output[0]

records = []

for log in freyja_logs:
    sample_id = Path(log).stem
    status = None
    with open(log) as f:
        for line in f:
            if "UserWarning: Solution may be inaccurate" in line:
                status = "solver_warning"
                break
            elif "SolverError" in line:
                status = "solver_error"
                break
        if status is None:
            status = "ok"
    records.append((sample_id, status))
    
report = pd.DataFrame.from_records(records, columns=["sample_id", "freyja_demix_status"])

report.to_csv(outfile, sep="\t", index=False)
