# -*- coding: utf-8 -*-
import sys
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import structure_lib
from structure_lib import Motif, expand, geometric_shells

# Cu2O cuprite
m = Motif(id="cu2o_cuprite", absorber="Cu", kind="crystal",
          nuclearity="crystal",
          source={"type": "manual", "ref": "Pn-3m, a=4.2696"},
          cell={"a": 4.2696}, spacegroup="P n -3 m",
          asym_sites=[
              {"el": "Cu", "x": 0.25, "y": 0.25, "z": 0.25, "wyckoff": "4b"},
              {"el": "O", "x": 0.0, "y": 0.0, "z": 0.0, "wyckoff": "2a"}],
          absorber_fract=[0.25, 0.25, 0.25])
try:
    sites = expand(m, rmax=4.0)
    print("cu2o sites", len(sites), geometric_shells(sites)[:4])
except Exception as exc:
    print("cu2o FAILED:", type(exc).__name__, exc)

# planar clusters
for el, n, c in (("Co", 1.96, 2.03647), ("Ni", 1.92, 2.01522),
                 ("Cu", 1.95, 2.03647)):
    sites = [(n, 0.0, 0.0, "N"), (-n, 0.0, 0.0, "N"),
             (0.0, n, 0.0, "N"), (0.0, -n, 0.0, "N")]
    for sx in (1, -1):
        for sy in (1, -1):
            sites.append((sx * c, sy * c, 0.0, "C"))
    print(el, geometric_shells(sites))

sites = [(2.31, 0.0, 0.0, "Cl"), (-2.31, 0.0, 0.0, "Cl"),
         (0.0, 2.31, 0.0, "Cl"), (0.0, -2.31, 0.0, "Cl")]
print("PtCl4", geometric_shells(sites))
