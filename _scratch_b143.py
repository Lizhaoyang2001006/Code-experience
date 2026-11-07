# -*- coding: utf-8 -*-
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import tempfile
from pathlib import Path

import numpy as np
from larch.xafs import feffpath

import artemis
import athena
from xafs_core import ChiK

MU_FE_O = (55.845, 15.999)
TRUE_THETA = 420.0
TRUE_STATIC = 0.0012
S02 = 0.803


def find_ss(files, target, tol=0.08):
    best = None
    for f in files:
        p = artemis.load_feff_path(f)
        if int(getattr(p, "nleg", 2)) != 2:
            continue
        if best is None or abs(float(p.reff) - target) < abs(
                float(best[1].reff) - target):
            best = (f, p)
    if best is None or abs(float(best[1].reff) - target) > tol:
        raise RuntimeError(f"没有 {target} Å 附近的 SS 路径")
    return best[0]


with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    sites = []
    for axis in ((2.0, 0, 0), (-2.0, 0, 0), (0, 2.0, 0),
                 (0, -2.0, 0), (0, 0, 2.0), (0, 0, -2.0)):
        sites.append((*axis, "O"))
    for axis in ((3.68, 0, 0), (-3.68, 0, 0), (0, 3.68, 0),
                 (0, -3.68, 0), (0, 0, 3.68), (0, 0, -3.68)):
        sites.append((*axis, "O"))
    inp = artemis.build_feff_input(sites, absorber="Fe", rmax=4.5,
                                   title="debye_test")
    files = artemis.run_feff6l(inp, tmp / "feff")
    near = find_ss(files, 2.0)
    far = find_ss(files, 3.68)
    print("near/far:", near.name, far.name)
    print("reff:", artemis.load_feff_path(near).reff,
          artemis.load_feff_path(far).reff)

    k = np.arange(3.0, 12.0, 0.05)
    datasets = []
    for T in (80.0, 350.0):
        plist = []
        s2s = []
        for f in (near, far):
            reff = float(artemis.load_feff_path(f).reff)
            s2 = artemis.correlated_debye_sigma2(
                T, TRUE_THETA, reff, *MU_FE_O,
                sigma2_static=TRUE_STATIC)
            s2s.append(s2)
            p = feffpath(str(f), s02=S02, sigma2=s2, e0=0.0, deltar=0.0)
            p.degen = 6.0
            p._calc_chi()
            plist.append(p)
        print(f"T={T} sigma2 near/far:", s2s)
        chi = artemis.exafs_model_k(plist, k)
        ck = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                  mu_norm=np.ones_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
        datasets.append({
            "name": f"{T:.0f}K", "chi_k": ck, "ft": ft,
            "datfiles": [str(near), str(far)],
            "rmin": 1.4, "rmax": 4.0, "kmin": 3.0, "kmax": 12.0,
            "kweight": 2, "s02": S02, "vary_n": False,
            "path_degen": [6.0, 6.0],
            "path_labels": ["Fe-O_near", "Fe-O_far"],
            "temperature_K": T, "edge_group": "K"})

    out = artemis.fit_multi_datasets(
        datasets, share_e0="auto",
        debye_models=[{
            "paths": [(0, 0), (0, 1), (1, 0), (1, 1)],
            "mass_absorber": MU_FE_O[0],
            "mass_scatterer": MU_FE_O[1],
            "theta_d_init": 400.0,
            "k_debye": None,
        }])
    print("nvars:", out["n_variables_total"])
    print("global Rf:", out["global_r_factor"])
    print("names:", out["parameter_names"])
    print("debye:", out["debye_results"])
    for fit in out["datasets"]:
        print(fit.delta_e0, [(p.label, p.sigma2, p.sigma2_err)
                             for p in fit.paths])
