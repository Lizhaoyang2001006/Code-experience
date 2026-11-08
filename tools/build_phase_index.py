# -*- coding: utf-8 -*-
"""build_phase_index.py — 生成 phase index 与覆盖报告（P1a）

用法::

    python tools/build_phase_index.py \\
        --pack _structure_lib/_phase_pack \\
        --out _structure_lib/phase_index.json \\
        --coverage _structure_lib/phase_coverage.json

P1a 的索引只来自 phase pack；``--cache`` 为 P1b 预留，当前只做提示。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import phase_pack                       # noqa: E402


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pack", default=str(phase_pack.PACK_ROOT))
    parser.add_argument("--cache", default=None,
                        help="P1b 预留：phase cache 根目录（P1a 忽略）")
    parser.add_argument("--out", default=str(phase_pack.INDEX_PATH))
    parser.add_argument("--coverage", default=str(phase_pack.COVERAGE_PATH))
    parser.add_argument("--motif-coverage",
                        default=str(phase_pack.MOTIF_COVERAGE_PATH))
    args = parser.parse_args(argv)

    pack_root = Path(args.pack)
    try:
        index = phase_pack.build_index(pack_root, out_path=args.out)
        coverage = phase_pack.build_coverage(
            pack_root, index=index, out_path=args.coverage,
            motif_coverage_path=args.motif_coverage)
    except phase_pack.PhasePackError as exc:
        print(f"[reject] {exc}", file=sys.stderr)
        return 1

    if args.cache:
        print("note: --cache 仅为 P1b 预留；P1a 索引只来自 phase pack")
    summary = index["summary"]
    print(f"entries={summary['entries']} "
          f"verified_cif={summary['verified_cif']} "
          f"elements_covered={len(summary['elements_covered'])}/"
          f"{summary['elements_total']}")
    print("verified_elements=" + ",".join(
        coverage["summary"]["verified_elements"]))
    print("not_covered_elements=" + ",".join(
        coverage["summary"]["not_covered_elements"]))
    print(f"index={Path(args.out)}")
    print(f"coverage={Path(args.coverage)}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
    raise SystemExit(main())
