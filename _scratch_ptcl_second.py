# -*- coding: utf-8 -*-
import sys
import tempfile
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import artemis
import xafs_pipeline

PT_CL = r"data/bmm_standards/Data/Pt/Pt-L3-PlatinumChloride.xdi"
res = xafs_pipeline.process_sample(
    [PT_CL], absorber="Pt", strict=False, verbose=False)


def find_ss(files, target, tol=0.06):
    best = None
    for f in files:
        p = artemis.load_feff_path(f)
        if int(getattr(p, "nleg", 2)) != 2:
            continue
        if best is None or abs(float(p.reff) - target) < abs(
                float(best[1].reff) - target):
            best = (f, p)
    if best is None or abs(float(best[1].reff) - target) > tol:
        raise RuntimeError(f"no SS at {target}: {[f.name for f in files]}")
    return best[0]


with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    for el in ("Pt", "Cl"):
        for d in (3.35, 3.45, 3.55, 3.65):
            sites = [(2.31, 0.0, 0.0, "Cl"), (-2.31, 0.0, 0.0, "Cl"),
                     (0.0, 2.31, 0.0, "Cl"), (0.0, -2.31, 0.0, "Cl"),
                     (0.0, 0.0, d, el)]
            inp = artemis.build_feff_input(
                sites, absorber="Pt", rmax=d + 1.0,
                title=f"ptcl_{el}_{d}", edge="L3")
            files = artemis.run_feff6l(inp, tmp / f"{el}_{d}")
            f1 = find_ss(files, 2.31)
            f2 = find_ss(files, d)
            fit = artemis.fit_r_space(
                res.chi_k, res.ft, [str(f1), str(f2)],
                rmin=1.7, rmax=3.8, kmin=3.0, kmax=12.0, kweight=2,
                s02=0.8075, vary_n=True,
                path_labels=["Pt-Cl", f"Pt-{el}"])
            p1, p2 = fit.paths
            print(f"el={el} d={d}: Rf={fit.r_factor:.4f}  "
                  f"N1={p1.degeneracy:.2f} R1={p1.r:.4f} "
                  f"s1={p1.sigma2:.4f} | "
                  f"N2={p2.degeneracy:.2f} R2={p2.r:.4f} "
                  f"s2={p2.sigma2:.4f}")
