# -*- coding: utf-8 -*-
import sys
import tempfile
from pathlib import Path

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import xafs_io
from xafs_core import XAFSConfigError

p = Path(r"data/练手数据/02_Fe_standards/Fe_standards.prj")
out = xafs_io.read_prj(p, as_scans=True)
print("meta", out["meta"])
print("labels", [s.meta["label"] for s in out["scans"]])
print("n include", len(xafs_io.read_prj(
    p, as_scans=True, include_derived=True)["scans"]))
raw = xafs_io.read_scan(
    r"data/练手数据/02_Fe_standards/Fe_lepidocrocite.100")
s = [x for x in out["scans"]
     if x.meta["label"] == "Fe_lepidocrocite.100"][0]
print("len", len(s), len(raw),
      "dE", float(np.max(np.abs(s.energy - raw.energy))),
      "corr", float(np.corrcoef(s.mu, raw.mu)[0, 1]))
with tempfile.TemporaryDirectory() as tmp:
    f = Path(tmp) / "x.prj"
    f.write_text("$VAR1 = { 'x' => 1 };", encoding="utf-8")
    try:
        xafs_io.read_prj(f)
        print("BAD: perl fragment accepted")
    except XAFSConfigError as exc:
        print("rejected OK:", str(exc)[:80])
