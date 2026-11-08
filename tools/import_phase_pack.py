# -*- coding: utf-8 -*-
"""import_phase_pack.py — 把本地真实 CIF 导入 phase pack（P1a）

用法::

    python tools/import_phase_pack.py --src D:\\cifs\\seed \\
        --pack _structure_lib/_phase_pack --source COD --tags metal

- 只读本地目录，不下载、不联网；``--offline`` 仅为显式声明。
- 目录内每个 ``*.cif`` 导入一条；同目录 ``<stem>.meta.json`` 可选。
- ``--id`` / ``--meta`` 只用于单文件导入。
- 同一 ``(source, id)`` 不同 sha256 必须 ``--replace``。
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import phase_pack                       # noqa: E402


def _gather(src: Path) -> list[Path]:
    if src.is_dir():
        return sorted(p for p in src.glob("*.cif") if p.is_file())
    if src.is_file():
        return [src]
    raise SystemExit(f"源路径不存在: {src}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", required=True,
                        help="单个 CIF 或包含 *.cif 的目录")
    parser.add_argument("--pack", default=str(phase_pack.PACK_ROOT),
                        help="phase pack 根目录")
    parser.add_argument("--source", default=None,
                        help="来源标签（COD / AMCSD / user）；meta 优先")
    parser.add_argument("--id", default=None, help="单文件导入时的条目 ID")
    parser.add_argument("--meta", default=None,
                        help="单文件导入时的 meta.json 路径")
    parser.add_argument("--model-level", default=None,
                        choices=phase_pack.VERIFIED_LEVELS)
    parser.add_argument("--tags", nargs="*", default=None,
                        help="覆盖/补充 tags，例如 metal oxide")
    parser.add_argument("--notes", default=None)
    parser.add_argument("--replace", action="store_true",
                        help="同 ID 不同 sha256 时显式覆盖旧文件")
    parser.add_argument("--offline", action="store_true",
                        help="声明离线；本脚本本来就不联网")
    args = parser.parse_args(argv)

    src = Path(args.src)
    files = _gather(src)
    if not files:
        print(f"没有找到 CIF: {src}", file=sys.stderr)
        return 2
    if len(files) > 1 and (args.id or args.meta):
        print("--id/--meta 只能用于单文件导入", file=sys.stderr)
        return 2

    pack_root = Path(args.pack)
    meta = None
    if args.meta:
        meta = json.loads(Path(args.meta).read_text(encoding="utf-8"))

    imported = 0
    for path in files:
        try:
            entry = phase_pack.import_cif(
                path, pack_root,
                source=args.source,
                id=args.id,
                meta=meta,
                model_level=args.model_level,
                tags=args.tags,
                notes=args.notes,
                replace=args.replace,
            )
        except phase_pack.PhasePackError as exc:
            print(f"[reject] {path.name}: {exc}", file=sys.stderr)
            return 1
        imported += 1
        print(f"[ok] {entry['source']}:{entry['id']} "
              f"sha256={entry['cif_sha256']} "
              f"elements={','.join(entry['absorber_elements'])} "
              f"sites={','.join(entry['absorber_sites'])}")

    manifest = phase_pack.load_manifest(pack_root)
    print(f"imported={imported} manifest_entries={len(manifest['entries'])}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
    raise SystemExit(main())
