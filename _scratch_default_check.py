# -*- coding: utf-8 -*-
import sys
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import xafs_fit
import xafs_pipeline
import xafs_report

lep = xafs_pipeline.process_sample(
    [r"data/练手数据/02_Fe_standards/Fe_lepidocrocite.000"],
    absorber="Fe", verbose=False)
out = xafs_fit.run(lep, identify=True, fit=True,
                   absorber_symbol="Fe", verbose=False)
fit = out["fit"]["fit"]
print("r_factor", repr(fit.r_factor))
for p in fit.paths:
    print("path", p.label, repr(p.degeneracy), repr(p.r), repr(p.sigma2))
rows = xafs_report.fit_table_rows([("lepido", [p.label for p in fit.paths], fit)])
print("table")
for row in rows:
    print(row)
