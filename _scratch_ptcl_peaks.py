# -*- coding: utf-8 -*-
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import artemis
import athena
import structure_id
import structure_lib
import xafs_pipeline
from absorber import get_absorber

PT_CL = r"data/bmm_standards/Data/Pt/Pt-L3-PlatinumChloride.xdi"
res = xafs_pipeline.process_sample(
    [PT_CL], absorber="Pt", strict=False, verbose=False)
absorber = get_absorber("Pt")
motif = [m for m in structure_lib.load_library("Pt")
         if m.id == "ptcl4_square"][0]
ft = res.ft
mask = (ft.r >= 1.0) & (ft.r <= 4.0)
exp_peaks = structure_id._peaks(ft.r, ft.magnitude, mask)
print("exp peaks", exp_peaks)
k = np.asarray(res.chi_k.k)
for s2 in (0.003, 0.005, 0.008):
    chi = structure_id.theory_chi(motif, k, absorber, sigma2=s2)
    th = athena.ft_to_r_space(k, chi, None, ft.kmin, ft.kmax, ft.kweight)
    mag = np.interp(ft.r, th.r, th.magnitude)
    print("theory s2", s2, "peaks", structure_id._peaks(ft.r, mag, mask))
print("feff paths:")
cache = structure_lib.ensure_feff(motif, edge="L3")
for f in sorted(cache.glob("feff*.dat")):
    p = artemis.load_feff_path(f)
    print(f.name, "nleg", getattr(p, "nleg", None), "reff", p.reff,
          "degen", getattr(p, "degen", None))
