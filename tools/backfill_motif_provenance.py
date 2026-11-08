# -*- coding: utf-8 -*-
"""backfill_motif_provenance.py — 给 motif 库补齐 P0b 来源字段

用法::

    python tools/backfill_motif_provenance.py
    python tools/backfill_motif_provenance.py --check

行为：

- 遍历 ``_structure_lib/<元素>/*.json``（跳过 ``_`` 开头的缓存目录）；
- 用 ``motif_provenance.classify_motif`` 分类；
- 只写 P0b 字段，不重算 sites/shell_expect；
- 写出 ``_structure_lib/provenance_report.json`` 与
  ``_structure_lib/motif_coverage.json``；
- 幂等：第二次运行不再改写任何文件（``--check`` 退出码 0）。

本脚本不会联网，也不会把现场生成的 CIF 升级为真实来源。
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import absorber as absorber_mod          # noqa: E402
import structure_lib                     # noqa: E402
from motif_provenance import (            # noqa: E402
    PERSISTED_FIELDS,
    classify_motif,
)


def _iter_motif_files(lib_root: Path):
    for element_dir in sorted(lib_root.iterdir()):
        if not element_dir.is_dir() or element_dir.name.startswith("_"):
            continue
        for path in sorted(element_dir.glob("*.json")):
            yield element_dir.name, path


def _pending_updates(data: dict, fields: dict) -> dict:
    updates = {}
    for key in PERSISTED_FIELDS:
        new = fields.get(key)
        if key == "edits" and new is None:
            new = []
        if new is None:
            # 不写 null：与 Motif.to_json() 的“跳过 None”保持一致。
            continue
        if key not in data or data[key] != new:
            updates[key] = new
    return updates


def _canonical(data: dict) -> dict:
    """按 Motif dataclass 字段顺序重排，保证 register 与 backfill 输出一致。"""
    ordered = {}
    for key in structure_lib.Motif.__dataclass_fields__:
        if key in data:
            ordered[key] = data[key]
    for key, value in data.items():
        if key not in ordered:
            ordered[key] = value
    return ordered


def _json_text(data) -> str:
    return json.dumps(data, ensure_ascii=False, indent=2) + "\n"


def _write_if_changed(path: Path, data) -> bool:
    text = _json_text(data)
    if path.exists() and path.read_text(encoding="utf-8") == text:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return True


def build_reports(entries: list[dict]) -> tuple[dict, dict]:
    status_counter = Counter(e["provenance_status"] for e in entries)
    level_counter = Counter(e["model_level"] for e in entries)
    element_counter = Counter(e["element"] for e in entries)

    report = {
        "schema_version": 1,
        "generated_by": "tools/backfill_motif_provenance.py",
        "library_root": "_structure_lib",
        "summary": {
            "motifs": len(entries),
            "by_status": dict(sorted(status_counter.items())),
            "by_model_level": dict(sorted(level_counter.items())),
            "by_element": dict(sorted(element_counter.items())),
            "unclassified": sum(
                1 for e in entries
                if e["provenance_status"] not in (
                    "verified_source_cif", "source_ref_only",
                    "generated_only", "constructed", "unverified")),
            "verified_source_cif": status_counter.get(
                "verified_source_cif", 0),
            "real_source_cifs": sum(
                1 for e in entries if e.get("source_cif_sha256")),
        },
        "entries": entries,
    }

    absorber_elements = sorted(absorber_mod.ABSORBERS.keys())
    by_element: dict[str, dict] = {}
    for element in sorted(set(absorber_elements) | set(element_counter)):
        rows = [e for e in entries if e["element"] == element]
        statuses = Counter(e["provenance_status"] for e in rows)
        by_element[element] = {
            "motifs": len(rows),
            "motif_ids": [e["id"] for e in rows],
            "verified_source_cif": statuses.get("verified_source_cif", 0),
            "source_ref_only": statuses.get("source_ref_only", 0),
            "generated_only": statuses.get("generated_only", 0),
            "constructed": statuses.get("constructed", 0),
            "unverified": statuses.get("unverified", 0),
            "real_source_cifs": sum(1 for e in rows
                                    if e.get("source_cif_sha256")),
            "formal_cn_ready": statuses.get("verified_source_cif", 0) > 0,
        }

    coverage = {
        "schema_version": 1,
        "generated_by": "tools/backfill_motif_provenance.py",
        "note": ("motif 覆盖只用于候选索引与审计；正式 CN 由真实 CIF / "
                 "phase pack 决定，不因 motif 稀疏而放行 prototype。"),
        "absorber_elements": absorber_elements,
        "covered_elements": [e for e in absorber_elements
                             if element_counter.get(e, 0) > 0],
        "missing_elements": [e for e in absorber_elements
                             if element_counter.get(e, 0) == 0],
        "library_extra_elements": sorted(set(element_counter)
                                         - set(absorber_elements)),
        "summary": {
            "elements_total": len(absorber_elements),
            "elements_covered": sum(
                1 for e in absorber_elements if element_counter.get(e, 0) > 0),
            "elements_missing": sum(
                1 for e in absorber_elements if element_counter.get(e, 0) == 0),
            "motifs": len(entries),
            "formal_cn_ready_elements": [
                e for e in absorber_elements
                if by_element.get(e, {}).get("formal_cn_ready")
            ],
        },
        "elements": by_element,
    }
    return report, coverage


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lib", default=str(ROOT / "_structure_lib"),
                        help="motif 库根目录（默认 <repo>/_structure_lib）")
    parser.add_argument("--report",
                        default=str(ROOT / "_structure_lib" /
                                    "provenance_report.json"))
    parser.add_argument("--coverage",
                        default=str(ROOT / "_structure_lib" /
                                    "motif_coverage.json"))
    parser.add_argument("--check", action="store_true",
                        help="只检查是否需要更新；需要更新时退出码 1")
    args = parser.parse_args(argv)

    lib_root = Path(args.lib)
    if not lib_root.is_dir():
        print(f"motif 库不存在: {lib_root}", file=sys.stderr)
        return 2

    entries = []
    pending_files = []
    for element, path in _iter_motif_files(lib_root):
        data = json.loads(path.read_text(encoding="utf-8"))
        motif = structure_lib.Motif.from_json(data)
        fields = classify_motif(motif)
        updates = _pending_updates(data, fields)
        data.update(updates)
        canonical = _canonical(data)
        new_text = json.dumps(canonical, ensure_ascii=False, indent=2)
        if path.read_text(encoding="utf-8") != new_text:
            pending_files.append((path, updates or {"__order__": True}))
            if not args.check:
                path.write_text(new_text, encoding="utf-8")
        entries.append({
            "element": element,
            "id": motif.id,
            "kind": motif.kind,
            "nuclearity": motif.nuclearity,
            "path": path.relative_to(ROOT).as_posix(),
            "source_type": (motif.source or {}).get("type"),
            "source_ref": fields["source_ref"],
            "source_cif_sha256": fields["source_cif_sha256"],
            "absorber_site": fields["absorber_site"],
            "occupancy_source": fields["occupancy_source"],
            "model_level": fields["model_level"],
            "provenance_status": fields["provenance_status"],
            "has_generated_cif": fields["_has_generated_cif"],
            "edits": fields["edits"],
        })

    entries.sort(key=lambda e: (e["element"], e["id"]))
    report, coverage = build_reports(entries)
    report_changed = False
    if not args.check:
        report_changed |= _write_if_changed(Path(args.report), report)
        report_changed |= _write_if_changed(Path(args.coverage), coverage)
    else:
        for target, payload in ((Path(args.report), report),
                                (Path(args.coverage), coverage)):
            if not target.exists() or \
                    target.read_text(encoding="utf-8") != _json_text(payload):
                report_changed = True

    print(f"motifs={len(entries)}")
    print("by_status=" + json.dumps(report["summary"]["by_status"],
                                    ensure_ascii=False, sort_keys=True))
    print(f"motif_json_updates={len(pending_files)}")
    print(f"reports_changed={int(report_changed)}")
    if pending_files:
        for path, updates in pending_files:
            print(f"  update {path.relative_to(ROOT).as_posix()}: "
                  f"{','.join(sorted(updates))}")
    if args.check and (pending_files or report_changed):
        return 1
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    raise SystemExit(main())
