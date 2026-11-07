# -*- coding: utf-8 -*-
import sys, tempfile
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from pathlib import Path
import xafs_pipeline, xafs_fit, xafs_diagnostics
ROOT = Path('.')
lep = next(ROOT.glob('data/*/02_Fe_standards/Fe_lepidocrocite.000'))
res = xafs_pipeline.process_sample([str(lep)], absorber='Fe', verbose=False)
out = xafs_fit.run(res, identify=True, fit=True, absorber_symbol='Fe', verbose=False)
with tempfile.TemporaryDirectory() as td:
    pc = xafs_diagnostics.write_phase_corrected_diagnostics(res, out, td)
    cc = xafs_diagnostics.write_chi2_contour_diagnostics(res, out, td)
    print('PC', pc['r_peak_raw'], pc['r_peak_corrected'], pc['r_fit'])
    print('files', pc['csv'].stat().st_size, pc['png'].stat().st_size)
    print('CC', cc['min_chi2'], cc['chi2_grid'].shape, cc['csv'].stat().st_size, cc['png'].stat().st_size)
