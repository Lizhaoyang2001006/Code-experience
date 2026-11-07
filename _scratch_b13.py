# -*- coding: utf-8 -*-
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from types import SimpleNamespace
import numpy as np
from larch.xafs import feffpath

import artemis
import athena
import structure_lib
import xafs_fit
from xafs_core import ChiK

S02 = 0.803


def make_res(paths, degens, noise, kmax=12.0):
    k = np.arange(3.0, kmax, 0.05)
    plist = []
    for f, d in zip(paths, degens):
        p = feffpath(str(f), s02=S02, sigma2=0.004, e0=2.0, deltar=0.0)
        p.degen = float(d)
        p._calc_chi()
        plist.append(p)
    chi = artemis.exafs_model_k(plist, k)
    if noise:
        rng = np.random.default_rng(13)
        chi = chi + rng.normal(0.0, noise, size=chi.shape)
    chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                 mu_norm=np.zeros_like(k), rbkg=1.0)
    ft = athena.ft_to_r_space(k, chi, None, 3.0, kmax, 2)
    return SimpleNamespace(chi_k=chi_k, ft=ft, absorber_symbol="Fe",
                           structure=None, merged=None, norm=None,
                           self_absorption=None)


sac = structure_lib.parse_catalytic_cluster_spec(
    "Fe", "sac:N4@2.00,axial=O@2.10,C=none")
for kmax in (12.0, 15.0):
    splan = xafs_fit.select_shell_paths(sac, edge="K", rmin=1.0, rmax=2.6,
                                        kmin=3.0, kmax=kmax)
    deg = list(splan["degens"])
    for i, sh in enumerate(splan["shells"]):
        if sh["scatterer"] == "O":
            deg[i] = 1.0
    for shared in (None, "all"):
        res_ax = make_res(splan["files"], deg, 1e-3, kmax=kmax)
        out3 = xafs_fit.fit_identified(
            res_ax, absorber_symbol="Fe",
            catalytic_model="sac:N4@2.00,axial=O@2.10,C=none",
            rmin=1.0, rmax=2.6, shared_sigma2=shared, verbose=False)
        print(f"AX kmax={kmax} shared={shared}: Rf={out3['fit'].r_factor:.5f}",
              [(p.label, round(p.degeneracy, 3), round(p.n_err, 3),
                round(p.sigma2, 4)) for p in out3["fit"].paths])
