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
import structure_lib
import xafs_fit
from xafs_core import ChiK

CIF = """data_test_hetero
_cell_length_a 4.10
_cell_length_b 4.10
_cell_length_c 4.10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'Pm-3m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Fe1 Fe 0.0 0.0 0.0 1.0
Ti1 Ti 0.5 0.5 0.5 0.7
Ru1 Ru 0.5 0.5 0.5 0.3
"""

with tempfile.TemporaryDirectory() as tmp:
    cif = Path(tmp) / "hetero.cif"
    cif.write_text(CIF, encoding="utf-8")
    motif = structure_lib.import_cif(cif, "test_hetero", "Fe")
    motif.cif_path = str(cif)
    motif.sites = structure_lib.expand_crystal(motif, rmax=5.0)
    motif.shell_expect = structure_lib.geometric_shells(motif.sites)
    print("sites:", motif.shell_expect[:4])
    plan = xafs_fit.select_shell_paths(
        motif, edge="K", rmin=2.8, rmax=4.3, kmin=3.0, kmax=12.0,
        hetero_cocupancy=True)
    print("hetero shells:", plan["hetero_shells"])
    print("files:", [Path(f).name for f in plan["files"]])
    print("labels:", plan["labels"], "degens:", plan["degens"])
    print("links:", plan["path_links"], "n_idp:", plan["n_independent"])

    s02 = 0.803
    k = np.arange(3.0, 12.0, 0.05)
    truth = {"N": [5.6, 2.4], "dr": -0.02, "sigma2": 0.004, "e0": 2.0}
    plist = []
    for f, d in zip(plan["files"], plan["degens"]):
        p = feffpath(str(f), s02=s02, sigma2=truth["sigma2"],
                     e0=truth["e0"], deltar=truth["dr"])
        p.degen = float(d)
        p._calc_chi()
        plist.append(p)
    chi = artemis.exafs_model_k(plist, k)
    rng = np.random.default_rng(7)
    chi = chi + rng.normal(0.0, 1e-3, size=chi.shape)
    ck = ChiK(k=k, chi=chi, background=np.zeros_like(k),
              mu_norm=np.zeros_like(k), rbkg=1.0)
    ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)

    fit = artemis.fit_r_space(
        ck, ft, plan["files"], rmin=plan["rmin"], rmax=plan["rmax"],
        kmin=3.0, kmax=12.0, kweight=2, s02=s02,
        vary_delta_e0=True, vary_n=True,
        path_degen=plan["degens"], path_labels=plan["labels"],
        path_links=plan["path_links"])
    print("fit nvars:", fit.n_variables, "Rf:", fit.r_factor)
    for p in fit.paths:
        print(" ", p.label, "N=", round(p.degeneracy, 4), "err=", p.n_err,
              "R=", p.r, "dR=", round(p.delta_r, 6),
              "s2=", round(p.sigma2, 6), "s2err=", p.sigma2_err)
    print("parameter names:", fit.parameter_names)

    # single-element approximation: only Ti path
    ti = [i for i, s in enumerate(plan["shells"]) if s["scatterer"] == "Ti"][0]
    fit1 = artemis.fit_r_space(
        ck, ft, [plan["files"][ti]], rmin=plan["rmin"], rmax=plan["rmax"],
        kmin=3.0, kmax=12.0, kweight=2, s02=s02,
        vary_delta_e0=True, vary_n=True,
        path_degen=[plan["degens"][ti]], path_labels=[plan["labels"][ti]])
    print("single Ti nvars:", fit1.n_variables, "Rf:", fit1.r_factor,
          "N:", fit1.paths[0].degeneracy)
