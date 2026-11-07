# -*- coding: utf-8 -*-
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import artemis
import athena
import structure_id
import structure_lib
import xafs_fit
import xafs_pipeline
from absorber import get_absorber
from structure_lib import Motif

PT_CL = r"data/bmm_standards/Data/Pt/Pt-L3-PlatinumChloride.xdi"

# 4 Cl first shell + 3 Cl second shell at 3.43 A (tetrahedral directions)
s = 3.43 / np.sqrt(3.0)
sites = [[2.31, 0.0, 0.0, "Cl"], [-2.31, 0.0, 0.0, "Cl"],
         [0.0, 2.31, 0.0, "Cl"], [0.0, -2.31, 0.0, "Cl"],
         [s, s, s, "Cl"], [s, -s, -s, "Cl"], [-s, s, -s, "Cl"]]
motif = Motif(id="ptcl4_square_proto", absorber="Pt", kind="cluster",
              nuclearity="m1", source={"type": "constructed"},
              sites=sites)
motif.shell_expect = structure_lib.geometric_shells(sites)
print("shells", motif.shell_expect)
cache = structure_lib.ensure_feff(motif, edge="L3")
print("cache", cache)
for f in sorted(cache.glob("feff*.dat")):
    p = artemis.load_feff_path(f)
    if int(getattr(p, "nleg", 2)) == 2:
        print("ss", f.name, p.reff, p.degen)

res = xafs_pipeline.process_sample(
    [PT_CL], absorber="Pt", strict=False, verbose=False)
absorber = get_absorber("Pt")
k = np.asarray(res.chi_k.k)
chi = structure_id.theory_chi(motif, k, absorber, sigma2=0.005)
th = athena.ft_to_r_space(k, chi, None, res.ft.kmin, res.ft.kmax,
                          res.ft.kweight)
mask = (res.ft.r >= 1.0) & (res.ft.r <= 4.0)
print("exp peaks", structure_id._peaks(res.ft.r, res.ft.magnitude, mask))
print("theory peaks", structure_id._peaks(
    res.ft.r, np.interp(res.ft.r, th.r, th.magnitude), mask))

plan = xafs_fit.select_shell_paths(motif, edge="L3", kmin=3.0, kmax=12.0)
print("plan rmin/rmax", plan["rmin"], plan["rmax"],
      [(sh["label"], sh["reff"], sh["degen"]) for sh in plan["shells"]],
      "trunc", plan["truncated"])
fit = artemis.fit_r_space(
    res.chi_k, res.ft, [str(f) for f in plan["files"]],
    plan["rmin"], plan["rmax"], kmin=3.0, kmax=12.0, kweight=2,
    s02=0.8075, vary_n=True, path_labels=[sh["label"] for sh in plan["shells"]])
print("fit Rf", fit.r_factor, "nvars", fit.n_variables,
      "nidp", fit.n_independent)
for p in fit.paths:
    print("path", p.label, "N", p.degeneracy, "R", p.r,
          "s2", p.sigma2)
