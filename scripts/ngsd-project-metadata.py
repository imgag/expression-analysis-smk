#!/usr/bin/env python3

import sys
from pathlib import Path
import subprocess
from io import StringIO
import pandas as pd

if len(sys.argv) == 2:
    project_name = Path(sys.argv[1]).resolve().name
else:
    project_name = Path.cwd().name

cmd = ["NGSDExportSamples", "-project", project_name, "-add_path", "SAMPLE_FOLDER"]
process = subprocess.run(cmd, capture_output=True, check=True)
csv = StringIO(process.stdout.decode())
df = pd.read_csv(csv, sep="\t")

df.rename(columns={"#name": "sample"}, inplace=True)
cols = ["sample", "name_external", "system_name", "run_name"]
df2 = df[cols].copy()
df2["file"] = df["path"] + df["sample"] + "_counts_raw.tsv"
df2["group"] = "control"
df2.sort_values(["sample"], inplace=True)

print(df2.to_csv(sep="\t", index=False))
