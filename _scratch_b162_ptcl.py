# -*- coding: utf-8 -*-
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import xafs_fit
import xafs_pipeline

PT_CL = r"data/bmm_standards/Data/Pt/Pt-L3-PlatinumChloride.xdi"
res = xafs_pipeline.process_sample(
    [PT_CL], absorber="Pt", strict=False, verbose=False)
out = xafs_fit.run(res, identify=True, fit=True,
                   absorber_symbol="Pt", verbose=True)
sid = out["structure"]
print("status", sid["status"], "uncovered", sid["uncovered_peaks"])
for row in sid["ranked"][:3]:
    print("rank", row["id"], row["score"], row.get("de0"))
fit = out["fit"]["fit"]
print("r_factor", fit.r_factor, "nvars", fit.n_variables,
      "nidp", fit.n_independent, "de0", fit.delta_e0)
for p in fit.paths:
    print("path", p.label, "N", p.degeneracy, "+-", p.n_err,
          "R", p.r, "+-", p.delta_r_err, "s2", p.sigma2)
print("fit mode", out["fit"].get("mode"), "motif", out["fit"].get("motif"))
