# -*- coding: utf-8 -*-
"""Scan second-shell candidates for PtCl2 against structure_id's score grid."""
import sys
import tempfile
from pathlib import Path

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import artemis
import athena
import structure_id
import xafs_pipeline
from absorber import get_absorber
from larch.xafs import feffpath

PT_CL = r"data/bmm_standards/Data/Pt/Pt-L3-PlatinumChloride.xdi"
res = xafs_pipeline.process_sample(
    [PT_CL], absorber="Pt", strict=False, verbose=False)
absorber = get_absorber("Pt")
k = np.asarray(res.chi_k.k)
ft = res.ft
mask = (ft.r >= 1.0) & (ft.r <= 4.0)
exp_peaks = structure_id._peaks(ft.r, ft.magnitude, mask)
print("exp peaks", exp_peaks)


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
        raise RuntimeError(f"no SS at {target}")
    return best[0]


def chi_from_path(path, sigma2):
    p = feffpath(str(path), s02=absorber.s02, sigma2=sigma2,
                 e0=0.0, deltar=0.0)
    p._calc_chi()
    return artemis.exafs_model_k([p], k)


with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    first_sites = [(2.31, 0.0, 0.0, "Cl"), (-2.31, 0.0, 0.0, "Cl"),
                   (0.0, 2.31, 0.0, "Cl"), (0.0, -2.31, 0.0, "Cl")]
    inp = artemis.build_feff_input(first_sites, absorber="Pt", rmax=3.3,
                                   title="first", edge="L3")
    files = artemis.run_feff6l(inp, tmp / "first")
    first_path = find_ss(files, 2.31)

    results = []
    for el in ("Cl", "Pt"):
        for d in (3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65):
            inp = artemis.build_feff_input(
                [(d, 0.0, 0.0, el)], absorber="Pt", rmax=d + 0.8,
                title=f"{el}{d}", edge="L3")
            files = artemis.run_feff6l(inp, tmp / f"{el}_{d}")
            second_path = find_ss(files, d)
            for s2 in structure_id.SIGMA2_GRID:
                chi1 = chi_from_path(first_path, s2)
                chi2 = chi_from_path(second_path, s2)
                for n2 in (6.0,):
                    chi0 = chi1 + n2 * chi2
                    best = None
                    best_mag = None
                    for de0 in structure_id.E0_GRID:
                        chi = structure_id.shift_e0(k, chi0, de0)
                        th = athena.ft_to_r_space(
                            k, chi, None, ft.kmin, ft.kmax, ft.kweight)
                        th_mag = np.interp(ft.r, th.r, th.magnitude)
                        th_c = (np.interp(ft.r, th.r, th.chi_r.real)
                                + 1j * np.interp(ft.r, th.r,
                                                 th.chi_r.imag))
                        s = structure_id.complex_score(
                            ft.chi_r, th_c, mask)
                        sm = structure_id.shape_score(
                            ft.magnitude, th_mag, mask)
                        if best is None or s < best[0]:
                            best = (s, de0, s2)
                        if best_mag is None or sm < best_mag[0]:
                            best_mag = (sm, de0, s2, th_mag)
                    peaks = structure_id._peaks(
                        ft.r, best_mag[3], mask)
                    covered = all(
                        np.any(np.abs(peaks - p) <= structure_id.PEAK_MATCH)
                        for p in exp_peaks)
                    results.append((best[0], el, d, n2, s2, best[1],
                                    tuple(np.round(peaks, 3)), covered))
    results.sort(key=lambda r: r[0])
    for row in results[:40]:
        print(f"score={row[0]:.4f} {row[1]} d={row[2]:.2f} n2={row[3]:.0f} "
              f"grid_s2={row[4]:.3f} de0={row[5]:+d} peaks={row[6]} "
              f"covered={row[7]}")
