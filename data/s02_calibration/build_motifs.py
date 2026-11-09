# -*- coding: utf-8 -*-
"""生成 T5 标定专用 motif JSON（**不写入生产 _structure_lib**）
================================================================
这些结构只服务于 ``data/s02_calibration`` 的标定产物：
  - 不进入 ``_structure_lib``；
  - 不改变 43 条生产 motif 的覆盖率/溯源断言（#119）；
  - 每个 JSON 都带**展开后的 sites**，``calibrate_s02 --motif-file``
    直接交给 FEFF，不会调用 ``_cif_path`` 往生产库写生成 CIF。

用法:
    python data/s02_calibration/build_motifs.py
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import structure_lib  # noqa: E402
from structure_lib import Motif, expand_full, geometric_shells  # noqa: E402

OUT = Path(__file__).resolve().parent / "motifs"


def _crystal(mid, absorber, cell, spacegroup, asym, absorber_fract,
             note):
    m = Motif(id=mid, absorber=absorber, kind="crystal",
              nuclearity="crystal",
              source={"type": "manual", "ref": note},
              cell=cell, spacegroup=spacegroup, asym_sites=asym,
              absorber_fract=absorber_fract,
              notes=f"T5 calibration-only: {note}")
    sites, occs = expand_full(m, 6.0)
    m.sites = [[round(float(v), 6) if isinstance(v, float) else v
                for v in s] for s in sites]
    m.shell_expect = geometric_shells(m.sites, occs=occs)
    return m


def _cluster(mid, absorber, sites, note):
    m = Motif(id=mid, absorber=absorber, kind="cluster", nuclearity="m1",
              source={"type": "manual", "ref": note},
              sites=[[float(x), float(y), float(z), el]
                     for x, y, z, el in sites],
              notes=f"T5 calibration-only: {note}")
    m.shell_expect = geometric_shells(m.sites)
    return m


def build():
    made = []

    made.append(_crystal(
        "cr2o3_corundum", "Cr",
        {"a": 4.958, "b": 4.958, "c": 13.594,
         "alpha": 90, "beta": 90, "gamma": 120},
        "R -3 c",
        [{"el": "Cr", "x": 0.0, "y": 0.0, "z": 0.3475, "wyckoff": "12c"},
         {"el": "O", "x": 0.305, "y": 0.0, "z": 0.25, "wyckoff": "18e"}],
        [0.0, 0.0, 0.3475],
        "Cr2O3 corundum a=4.958 c=13.594 A"))

    made.append(_crystal(
        "mno2_rutile", "Mn",
        {"a": 4.404, "b": 4.404, "c": 2.876,
         "alpha": 90, "beta": 90, "gamma": 90},
        "P 42/m n m",
        [{"el": "Mn", "x": 0.0, "y": 0.0, "z": 0.0, "wyckoff": "2a"},
         {"el": "O", "x": 0.305, "y": 0.305, "z": 0.0, "wyckoff": "4f"}],
        [0.0, 0.0, 0.0],
        "MnO2 pyrolusite/rutile a=4.404 c=2.876 A"))

    made.append(_crystal(
        "zno_wurtzite", "Zn",
        {"a": 3.2495, "b": 3.2495, "c": 5.2069,
         "alpha": 90, "beta": 90, "gamma": 120},
        "P 63 m c",
        [{"el": "Zn", "x": 1.0 / 3.0, "y": 2.0 / 3.0, "z": 0.0,
          "wyckoff": "2b"},
         {"el": "O", "x": 1.0 / 3.0, "y": 2.0 / 3.0, "z": 0.382,
          "wyckoff": "2b"}],
        [1.0 / 3.0, 2.0 / 3.0, 0.0],
        "ZnO wurtzite a=3.2495 c=5.2069 A"))

    made.append(_crystal(
        "ruo2_rutile", "Ru",
        {"a": 4.4919, "b": 4.4919, "c": 3.1066,
         "alpha": 90, "beta": 90, "gamma": 90},
        "P 42/m n m",
        [{"el": "Ru", "x": 0.0, "y": 0.0, "z": 0.0, "wyckoff": "2a"},
         {"el": "O", "x": 0.3053, "y": 0.3053, "z": 0.0, "wyckoff": "4f"}],
        [0.0, 0.0, 0.0],
        "RuO2 rutile a=4.4919 c=3.1066 A"))

    # ---- Mo 失败尝试结构：保留以便复核“为何未标定” ----
    made.append(_crystal(
        "moo2_rutile", "Mo",
        {"a": 4.859, "b": 4.859, "c": 2.799,
         "alpha": 90, "beta": 90, "gamma": 90},
        "P 42/m n m",
        [{"el": "Mo", "x": 0.0, "y": 0.0, "z": 0.0, "wyckoff": "2a"},
         {"el": "O", "x": 0.302, "y": 0.302, "z": 0.0, "wyckoff": "4f"}],
        [0.0, 0.0, 0.0],
        "MoO2 rutile a=4.859 c=2.799 A"))

    made.append(_crystal(
        "mos2_2h", "Mo",
        {"a": 3.161, "b": 3.161, "c": 12.299,
         "alpha": 90, "beta": 90, "gamma": 120},
        "P 63/m m c",
        [{"el": "Mo", "x": 1.0 / 3.0, "y": 2.0 / 3.0, "z": 0.25,
          "wyckoff": "2c"},
         {"el": "S", "x": 1.0 / 3.0, "y": 2.0 / 3.0, "z": 0.621,
          "wyckoff": "4f"}],
        [1.0 / 3.0, 2.0 / 3.0, 0.25],
        "2H-MoS2 a=3.161 c=12.299 A"))

    r = 1.76
    made.append(_cluster(
        "moo4_tetra", "Mo",
        [(r, 0.0, 0.0, "O"), (-r, 0.0, 0.0, "O"),
         (0.0, r, 0.0, "O"), (0.0, 0.0, r, "O")],
        "MoO4 tetrahedron Mo-O 1.76 A x4"))

    # UO2 fluorite（P5-C COD 9009049）：U 的第一壳是 8×U–O 2.3678 Å，
    # 不是 α-U 金属的 U–U 壳；标定专用，不进生产 _structure_lib。
    made.append(_crystal(
        "uo2_fluorite", "U",
        {"a": 5.4682, "b": 5.4682, "c": 5.4682,
         "alpha": 90, "beta": 90, "gamma": 90},
        "F m -3 m",
        [{"el": "U", "x": 0.0, "y": 0.0, "z": 0.0, "wyckoff": "4a"},
         {"el": "O", "x": 0.25, "y": 0.25, "z": 0.25, "wyckoff": "8c"}],
        [0.0, 0.0, 0.0],
        "UO2 fluorite COD 9009049 a=5.4682 A"))

    OUT.mkdir(parents=True, exist_ok=True)
    for m in made:
        payload = m.to_json()
        payload["calibration_only"] = True
        payload["production_library"] = False
        (OUT / f"{m.id}.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=1) + "\n",
            encoding="utf-8")
        print(f"wrote {m.id}.json  sites={len(m.sites)}  "
              f"first_shell={m.shell_expect[0] if m.shell_expect else None}")


if __name__ == "__main__":
    # _cif_path 只在 expand_full 内使用；把 LIB_ROOT 临时指向系统临时
    # 目录，避免生成 CIF 落进生产 _structure_lib。
    original = structure_lib.LIB_ROOT
    with tempfile.TemporaryDirectory() as tmp:
        structure_lib.LIB_ROOT = Path(tmp)
        try:
            build()
        finally:
            structure_lib.LIB_ROOT = original
