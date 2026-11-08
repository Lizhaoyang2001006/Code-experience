# -*- coding: utf-8 -*-
"""phase_pack.py — 正式 CN 的本地真实 CIF 包（P1a）
================================================================
phase pack 是**离线、可复现、可增量扩展**的结构来源，不是 motif 库，
也不是 COD/ICSD 全库镜像：

    _structure_lib/_phase_pack/
      manifest.json
      COD/1100108/<sha256>.cif
      COD/1100108/<sha256>.meta.json
      user/my_sample_Fe1/<sha256>.cif
      user/my_sample_Fe1/<sha256>.meta.json

设计约束：

- 只接收真实实验 CIF；``structure_lib._cif_path`` 现场生成的小 CIF 拒绝入库；
- 每条必须有 ``source/id/formula/spacegroup/cell/absorber_sites/
  absorber_elements/cif_sha256/experimental/model_level``；
- 同一 ``(source, id)`` 不同 sha256 默认冲突，显式 ``replace=True`` 才覆盖；
- 导入时只读本地文件，不下载、不联网；
- 索引与覆盖报告是确定性 JSON（不含时间戳/绝对路径），可直接提交。

本模块不做相解析（P1b 的 ``phase_resolver``），也不计算 N₀。
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

import absorber as absorber_mod
from motif_provenance import is_generated_cif, sha256_file
from xafs_core import XAFSConfigError

ROOT = Path(__file__).resolve().parent
PACK_ROOT = ROOT / "_structure_lib" / "_phase_pack"
MANIFEST_NAME = "manifest.json"
INDEX_PATH = ROOT / "_structure_lib" / "phase_index.json"
COVERAGE_PATH = ROOT / "_structure_lib" / "phase_coverage.json"
MOTIF_COVERAGE_PATH = ROOT / "_structure_lib" / "motif_coverage.json"

# 17 个吸收元素；顺序沿用 absorber.ABSORBERS，报告更稳定。
ABSORBER_ELEMENTS = tuple(absorber_mod.ABSORBERS)

ENTRY_REQUIRED = (
    "source", "id", "formula", "spacegroup", "cell", "absorber_sites",
    "absorber_elements", "cif_path", "cif_sha256", "experimental",
    "model_level",
)
VERIFIED_LEVELS = ("crystal_cif", "crystal_edited")

# 常见相清单：只用于覆盖报告的“还缺什么”，不参与自动选相。
COMMON_PHASES = {
    "Fe": ["metal", "FeO", "Fe3O4", "Fe2O3", "FeOOH"],
    "Co": ["metal", "CoO", "Co3O4"],
    "Ni": ["metal", "NiO", "Ni(OH)2"],
    "Cu": ["metal", "Cu2O", "CuO"],
    "Pt": ["metal", "PtO2"],
    "Cr": ["metal", "Cr2O3"],
    "Mn": ["metal", "MnO", "MnO2", "Mn3O4"],
    "Zn": ["metal", "ZnO"],
    "Mo": ["metal", "MoO2", "MoO3"],
    "Ru": ["metal", "RuO2"],
    "Rh": ["metal", "Rh2O3"],
    "Pd": ["metal", "PdO"],
    "Ag": ["metal", "Ag2O"],
    "Ir": ["metal", "IrO2"],
    "Au": ["metal", "Au2O3"],
    "Ce": ["CeO2"],
    "U": ["UO2", "U3O8"],
}

_SHA_RE = re.compile(r"^[0-9a-f]{64}$")
_ELEMENT_RE = re.compile(r"[A-Za-z]+")
_FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")


class PhasePackError(XAFSConfigError):
    """phase pack 的拒绝错误；``code`` 便于测试与报告引用。"""

    def __init__(self, code: str, message: str):
        self.code = str(code)
        self.detail = str(message)
        super().__init__("phase_pack", f"{self.code}: {self.detail}")


# ---------------------------------------------------------------
# 基础工具
# ---------------------------------------------------------------
def normalize_element(symbol) -> str:
    """把 ``Fe`` / ``Fe3+`` / ``Pt0`` 归一为 ``Fe`` / ``Fe`` / ``Pt``。"""
    m = _ELEMENT_RE.match(str(symbol or "").strip())
    return m.group(0).capitalize() if m else ""


def _norm_spacegroup(value) -> str:
    return "".join(str(value or "").split()).lower()


def _cell_dict(cell) -> dict:
    return {
        "a": round(float(cell.a), 6),
        "b": round(float(cell.b), 6),
        "c": round(float(cell.c), 6),
        "alpha": round(float(cell.alpha), 6),
        "beta": round(float(cell.beta), 6),
        "gamma": round(float(cell.gamma), 6),
    }


def _formula_elements(formula) -> set[str]:
    return {el for el, _ in _FORMULA_RE.findall(str(formula or ""))}


def _derive_formula(sites: list[dict]) -> str:
    counts: dict[str, int] = {}
    for site in sites:
        el = site["element"]
        counts[el] = counts.get(el, 0) + 1
    return "".join(f"{el}{count if count > 1 else ''}"
                   for el, count in sorted(counts.items()))


def _derive_id(stem: str) -> str:
    m = re.search(r"(\d{5,})", stem)
    if m:
        return m.group(1)
    return re.sub(r"[^0-9A-Za-z._-]+", "_", stem).strip("_") or "unknown"


def _safe_component(value, field: str) -> str:
    text = str(value or "").strip()
    if not text or text in (".", "..") or any(
            sep in text for sep in ("/", "\\", "\x00")):
        raise PhasePackError("invalid_path_component",
                             f"{field} 不能为空或含路径分隔符: {value!r}")
    return text


def _json_text(data) -> str:
    return json.dumps(data, ensure_ascii=False, indent=2) + "\n"


def _write_json_if_changed(path: Path, data) -> bool:
    text = _json_text(data)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_text(encoding="utf-8") == text:
        return False
    path.write_text(text, encoding="utf-8")
    return True


def _is_within(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except (OSError, ValueError):
        return False


def _entry_key(entry: dict) -> tuple[str, str]:
    return (str(entry.get("source", "")).lower(),
            str(entry.get("id", "")))


# ---------------------------------------------------------------
# CIF 解析
# ---------------------------------------------------------------
def parse_cif(path: str | Path) -> dict:
    """用 gemmi 解析 CIF，返回 cell / 空间群 / 位点 / 吸收元素。"""
    path = Path(path)
    if not path.is_file():
        raise PhasePackError("cif_missing", f"CIF 不存在: {path}")
    try:
        import gemmi
    except ImportError as exc:  # pragma: no cover - gemmi 是硬依赖
        raise PhasePackError("gemmi_missing", str(exc)) from exc

    try:
        doc = gemmi.cif.read(str(path))
        if len(doc) == 0:
            raise ValueError("CIF 没有数据块")
        block = doc.sole_block() if len(doc) == 1 else doc[0]
        st = gemmi.make_small_structure_from_block(block)
    except Exception as exc:
        raise PhasePackError("cif_parse_failed",
                             f"{path.name}: {exc}") from exc

    sites = []
    for index, site in enumerate(st.sites, 1):
        element = normalize_element(site.type_symbol)
        if not element:
            element = normalize_element(site.label)
        if not element:
            continue
        label = str(site.label or f"{element}{index}").strip()
        occ = float(site.occ)
        if not math.isfinite(occ):
            occ = 1.0
        sites.append({
            "label": label,
            "element": element,
            "fract": [round(float(site.fract.x), 6),
                      round(float(site.fract.y), 6),
                      round(float(site.fract.z), 6)],
            "occ": round(occ, 6),
        })
    if not sites:
        raise PhasePackError("cif_no_sites", f"{path.name}: 没有可解析的原子位点")

    formula = ""
    try:
        formula = str(block.find_value("_chemical_formula_sum") or "").strip()
    except Exception:
        formula = ""
    if not formula:
        formula = _derive_formula(sites)

    spacegroup = " ".join(str(st.spacegroup_hm or "").split())
    if not spacegroup:
        raise PhasePackError("cif_no_spacegroup",
                             f"{path.name}: 缺空间群")

    known = set(ABSORBER_ELEMENTS)
    default_sites = [s["label"] for s in sites if s["element"] in known]
    return {
        "path": str(path),
        "formula": formula,
        "formula_elements": sorted(_formula_elements(formula)
                                   | {s["element"] for s in sites}),
        "spacegroup": spacegroup,
        "cell": _cell_dict(st.cell),
        "sites": sites,
        "default_absorber_sites": default_sites,
        "absorber_elements": sorted({s["element"] for s in sites
                                     if s["element"] in known}),
    }


# ---------------------------------------------------------------
# 条目校验
# ---------------------------------------------------------------
def _normalize_absorber_sites(parsed: dict, values) -> tuple[list[str], list[str]]:
    labels = {s["label"]: s for s in parsed["sites"]}
    lowered = {k.lower(): v for k, v in labels.items()}
    selected = []
    for raw in values:
        label = str(raw).strip()
        site = labels.get(label) or lowered.get(label.lower())
        if site is None:
            raise PhasePackError(
                "absorber_site_missing",
                f"吸收位 {label!r} 不在 CIF 位点中: "
                f"{sorted(labels)}")
        if site["element"] not in ABSORBER_ELEMENTS:
            raise PhasePackError(
                "not_absorber_element",
                f"位点 {label!r} 的元素 {site['element']} 不在 17 个吸收元素中")
        if site["label"] not in selected:
            selected.append(site["label"])
    if not selected:
        raise PhasePackError("no_absorber_site",
                             "没有指定任何吸收位")
    elements = sorted({labels[label]["element"] for label in selected})
    return selected, elements


def _check_cell(parsed_cell: dict, claimed: dict) -> None:
    for key in ("a", "b", "c", "alpha", "beta", "gamma"):
        try:
            value = float(claimed[key])
        except (KeyError, TypeError, ValueError):
            raise PhasePackError("invalid_cell",
                                 f"cell.{key} 不是数字: {claimed.get(key)!r}")
        if not math.isfinite(value) or value <= 0:
            raise PhasePackError("invalid_cell",
                                 f"cell.{key} 必须是正有限数: {value!r}")
        ref = float(parsed_cell[key])
        if abs(value - ref) > max(1e-3, 1e-4 * abs(ref)):
            raise PhasePackError(
                "cell_mismatch",
                f"cell.{key} 声明 {value} 与 CIF {ref} 不一致")


def _validate_entry(pack_root: Path, entry: dict) -> None:
    if not isinstance(entry, dict):
        raise PhasePackError("invalid_entry", "manifest 条目必须是对象")
    missing = [key for key in ENTRY_REQUIRED if key not in entry]
    if missing:
        raise PhasePackError("missing_field",
                             f"manifest 条目缺字段: {missing}")

    for key in ("source", "id", "formula", "spacegroup", "cif_path"):
        if not str(entry[key] or "").strip():
            raise PhasePackError("empty_field", f"字段 {key} 不能为空")
    if not str(entry["source"]).strip() or not str(entry["id"]).strip():
        raise PhasePackError("empty_field", "source/id 不能为空")

    cif_rel = Path(str(entry["cif_path"]))
    if cif_rel.is_absolute() or ".." in cif_rel.parts:
        raise PhasePackError("unsafe_cif_path",
                             f"cif_path 必须是 pack 内相对路径: {cif_rel}")
    cif_path = pack_root / cif_rel
    if not cif_path.is_file():
        raise PhasePackError("cif_missing", f"CIF 不存在: {cif_path}")

    sha = str(entry["cif_sha256"]).lower()
    if not _SHA_RE.match(sha):
        raise PhasePackError("invalid_sha256",
                             f"cif_sha256 格式非法: {entry['cif_sha256']!r}")
    actual = sha256_file(cif_path)
    if actual != sha:
        raise PhasePackError(
            "sha256_mismatch",
            f"{cif_rel}: manifest {sha} != 文件 {actual}")

    if is_generated_cif(cif_path):
        raise PhasePackError("generated_cif_not_allowed",
                             f"{cif_rel}: 现场生成 CIF 不得进入 phase pack")

    parsed = parse_cif(cif_path)
    raw_sites = entry["absorber_sites"]
    if not isinstance(raw_sites, list) or not raw_sites:
        raise PhasePackError("no_absorber_site",
                             "absorber_sites 必须是非空数组")
    selected, elements = _normalize_absorber_sites(parsed, raw_sites)
    if selected != [str(s) for s in raw_sites]:
        raise PhasePackError(
            "absorber_site_alias",
            f"absorber_sites 必须使用 CIF 原始标签: {selected}")

    claimed_elements = sorted(str(e).capitalize()
                              for e in entry["absorber_elements"])
    if claimed_elements != elements:
        raise PhasePackError(
            "absorber_elements_mismatch",
            f"absorber_elements 声明 {claimed_elements} 与位点 {elements} 不一致")

    formula_elements = set(parsed["formula_elements"])
    lacking = [el for el in elements if el not in formula_elements]
    if lacking:
        raise PhasePackError(
            "formula_missing_element",
            f"formula {entry['formula']!r} 缺吸收元素 {lacking}")

    if _norm_spacegroup(entry["spacegroup"]) != _norm_spacegroup(
            parsed["spacegroup"]):
        raise PhasePackError(
            "spacegroup_mismatch",
            f"空间群声明 {entry['spacegroup']!r} 与 CIF "
            f"{parsed['spacegroup']!r} 不一致")
    _check_cell(parsed["cell"], entry["cell"])

    if entry["experimental"] is not True:
        raise PhasePackError("non_experimental_cif",
                             "phase pack 只接收 experimental=true 的 CIF")
    if entry["model_level"] not in VERIFIED_LEVELS:
        raise PhasePackError(
            "invalid_model_level",
            f"真实 CIF 的 model_level 必须是 {VERIFIED_LEVELS}，"
            f"收到 {entry['model_level']!r}")


def validate_entry(pack_root: str | Path, entry: dict) -> dict:
    """返回 ``{ok, code, reason}``，不抛异常。"""
    try:
        _validate_entry(Path(pack_root), entry)
        return {"ok": True, "code": "ok", "reason": ""}
    except PhasePackError as exc:
        return {"ok": False, "code": exc.code, "reason": exc.detail}


# ---------------------------------------------------------------
# 导入
# ---------------------------------------------------------------
def _load_companion_meta(src: Path, meta: dict | None) -> dict:
    if meta is not None:
        return dict(meta)
    for candidate in (src.with_suffix(".meta.json"),
                      src.with_name(src.name + ".meta.json")):
        if candidate.is_file():
            try:
                data = json.loads(candidate.read_text(encoding="utf-8"))
            except json.JSONDecodeError as exc:
                raise PhasePackError("invalid_meta",
                                     f"{candidate.name}: {exc}") from exc
            if not isinstance(data, dict):
                raise PhasePackError("invalid_meta",
                                     f"{candidate.name}: 根节点必须是对象")
            return data
    return {}


def import_cif(src_cif: str | Path, pack_root: str | Path, *,
               source: str | None = None, id: str | None = None,
               meta: dict | None = None,
               absorber_sites=None, model_level: str | None = None,
               experimental: bool | None = None,
               tags=None, notes: str | None = None,
               added_by: str = "tools/import_phase_pack.py",
               replace: bool = False) -> dict:
    """把一条本地真实 CIF 导入 phase pack，返回 manifest 条目。"""
    src = Path(src_cif)
    pack_root = Path(pack_root)
    if is_generated_cif(src):
        raise PhasePackError("generated_cif_not_allowed",
                             f"{src.name}: 现场生成 CIF 不得进入 phase pack")

    data = _load_companion_meta(src, meta)
    src_label = _safe_component(
        data.get("source") or source or "user", "source")
    motif_id = _safe_component(
        data.get("id") or id or _derive_id(src.stem), "id")

    parsed = parse_cif(src)
    sha = sha256_file(src)
    claimed_sha = data.get("cif_sha256")
    if claimed_sha and str(claimed_sha).lower() != sha:
        raise PhasePackError(
            "sha256_mismatch",
            f"meta 声明 {claimed_sha}，实际 {sha}")

    raw_sites = (data.get("absorber_sites") or absorber_sites
                 or parsed["default_absorber_sites"])
    if isinstance(raw_sites, str):
        raw_sites = [s for s in re.split(r"[,;\s]+", raw_sites) if s]
    selected, elements = _normalize_absorber_sites(parsed, raw_sites)

    formula = str(data.get("formula") or parsed["formula"]).strip()
    formula_elements = set(parsed["formula_elements"]) | _formula_elements(formula)
    lacking = [el for el in elements if el not in formula_elements]
    if lacking:
        raise PhasePackError(
            "formula_missing_element",
            f"formula {formula!r} 缺吸收元素 {lacking}")

    cell = parsed["cell"]
    if data.get("cell") is not None:
        _check_cell(cell, data["cell"])
    if data.get("spacegroup") and (
            _norm_spacegroup(data["spacegroup"]) !=
            _norm_spacegroup(parsed["spacegroup"])):
        raise PhasePackError(
            "spacegroup_mismatch",
            f"meta 空间群 {data['spacegroup']!r} 与 CIF "
            f"{parsed['spacegroup']!r} 不一致")

    is_experimental = (data.get("experimental")
                       if data.get("experimental") is not None
                       else (True if experimental is None else experimental))
    if is_experimental is not True:
        raise PhasePackError("non_experimental_cif",
                             "phase pack 只接收 experimental=true 的 CIF")
    level = str(data.get("model_level") or model_level or "crystal_cif")
    if level not in VERIFIED_LEVELS:
        raise PhasePackError(
            "invalid_model_level",
            f"真实 CIF 的 model_level 必须是 {VERIFIED_LEVELS}，收到 {level!r}")

    entry = {
        "source": src_label,
        "id": motif_id,
        "formula": formula,
        "spacegroup": parsed["spacegroup"],
        "cell": cell,
        "absorber_sites": selected,
        "absorber_elements": elements,
        "cif_path": "",  # 写入前按实际 sha 填
        "cif_sha256": sha,
        "experimental": True,
        "model_level": level,
        "tags": [str(t) for t in (data.get("tags") or tags or [])],
        "added_by": str(data.get("added_by") or added_by),
        "original_filename": src.name,
    }
    if notes is not None:
        entry["notes"] = str(notes)
    elif data.get("notes"):
        entry["notes"] = str(data["notes"])

    target_dir = pack_root / src_label / motif_id
    target_dir.mkdir(parents=True, exist_ok=True)
    target_cif = target_dir / f"{sha}.cif"
    target_meta = target_dir / f"{sha}.meta.json"
    entry["cif_path"] = target_cif.relative_to(pack_root).as_posix()

    old_files: list[Path] = []
    for existing in sorted(target_dir.glob("*.meta.json")):
        try:
            old = json.loads(existing.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        old_sha = str(old.get("cif_sha256") or "").lower()
        if old_sha and old_sha != sha:
            if not replace:
                raise PhasePackError(
                    "id_sha_conflict",
                    f"{src_label}:{motif_id} 已存在 sha256={old_sha}，"
                    f"新文件 sha256={sha}；需要显式 replace=True")
            old_files.extend([existing,
                              existing.with_name(
                                  existing.name[:-len(".meta.json")] + ".cif")])

    if target_cif.exists() and sha256_file(target_cif) != sha:
        raise PhasePackError("target_sha_mismatch",
                             f"目标文件已存在但哈希不同: {target_cif}")
    if not target_cif.exists():
        target_cif.write_bytes(src.read_bytes())
    _write_json_if_changed(target_meta, entry)
    for old in old_files:
        if old != target_cif and old != target_meta and old.exists():
            old.unlink()

    rebuild_manifest(pack_root)
    return entry


# ---------------------------------------------------------------
# manifest / index / coverage
# ---------------------------------------------------------------
def _scan_entries(pack_root: Path) -> list[dict]:
    entries = []
    for meta_path in sorted(pack_root.glob("*/*/*.meta.json")):
        parts = meta_path.relative_to(pack_root).parts
        if any(part.startswith("_") for part in parts[:-1]):
            continue
        try:
            entry = json.loads(meta_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise PhasePackError("invalid_meta",
                                 f"{meta_path}: {exc}") from exc
        if not isinstance(entry, dict):
            raise PhasePackError("invalid_meta",
                                 f"{meta_path}: 根节点必须是对象")
        actual_cif = meta_path.with_name(
            meta_path.name[:-len(".meta.json")] + ".cif")
        if not actual_cif.is_file():
            raise PhasePackError("cif_missing",
                                 f"{meta_path}: 缺对应 CIF {actual_cif.name}")
        entry["cif_path"] = actual_cif.relative_to(pack_root).as_posix()
        entries.append(entry)
    entries.sort(key=lambda e: (str(e.get("source", "")).lower(),
                                str(e.get("id", ""))))
    return entries


def rebuild_manifest(pack_root: str | Path, *, pack_id: str | None = None,
                     notes: str | None = None) -> dict:
    """扫描 pack 目录，校验全部条目并重写 manifest.json。"""
    pack_root = Path(pack_root)
    entries = _scan_entries(pack_root)
    for entry in entries:
        result = validate_entry(pack_root, entry)
        if not result["ok"]:
            raise PhasePackError(
                result["code"],
                f"{entry.get('source')}:{entry.get('id')} — {result['reason']}")

    manifest_path = pack_root / MANIFEST_NAME
    existing = {}
    if manifest_path.is_file():
        try:
            existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            existing = {}
    manifest = {
        "schema_version": 1,
        "pack_id": (pack_id or existing.get("pack_id")
                    or "formal-cn-phase-pack"),
        "notes": (notes if notes is not None
                  else existing.get("notes")
                  or "本地真实 CIF 包；不镜像 COD/ICSD，运行期离线。"),
        "entries": entries,
    }
    _write_json_if_changed(manifest_path, manifest)
    return manifest


def load_manifest(pack_root: str | Path) -> dict:
    path = Path(pack_root) / MANIFEST_NAME
    if not path.is_file():
        raise PhasePackError("manifest_missing", f"manifest 不存在: {path}")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise PhasePackError("invalid_manifest", f"{path}: {exc}") from exc
    if not isinstance(data, dict) or not isinstance(data.get("entries"), list):
        raise PhasePackError("invalid_manifest",
                             f"{path}: 缺 entries 数组")
    return data


def _repo_relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except (OSError, ValueError):
        return path.resolve().as_posix()


def build_index(pack_root: str | Path, *, out_path: str | Path | None = None,
                manifest: dict | None = None) -> dict:
    """生成 phase_index.json；所有条目重新校验 sha256 与 CIF 内容。"""
    pack_root = Path(pack_root)
    manifest = manifest or load_manifest(pack_root)
    index_entries = []
    for entry in manifest["entries"]:
        result = validate_entry(pack_root, entry)
        if not result["ok"]:
            raise PhasePackError(
                result["code"],
                f"{entry.get('source')}:{entry.get('id')} — {result['reason']}")
        cif = pack_root / entry["cif_path"]
        row = dict(entry)
        row["pack_path"] = entry["cif_path"]
        row["cif_path"] = _repo_relative(cif)
        index_entries.append(row)
    index_entries.sort(key=lambda e: (str(e["source"]).lower(), str(e["id"])))

    elements = sorted({el for row in index_entries
                       for el in row["absorber_elements"]})
    index = {
        "schema_version": 1,
        "generated_by": "tools/build_phase_index.py",
        "pack_id": manifest.get("pack_id"),
        "pack_root": _repo_relative(pack_root),
        "manifest_sha256": hashlib.sha256(
            _json_text(manifest).encode("utf-8")).hexdigest(),
        "summary": {
            "entries": len(index_entries),
            "verified_cif": len(index_entries),
            "elements_covered": elements,
            "elements_total": len(ABSORBER_ELEMENTS),
        },
        "entries": index_entries,
    }
    if out_path is not None:
        _write_json_if_changed(Path(out_path), index)
    return index


def _phase_tag_match(entry: dict, phase: str) -> bool:
    phase_low = phase.lower()
    tags = [str(t).lower() for t in entry.get("tags") or []]
    return any(tag and tag in phase_low for tag in tags)


def build_coverage(pack_root: str | Path, *,
                   index: dict | None = None,
                   out_path: str | Path | None = None,
                   motif_coverage_path: str | Path | None = None) -> dict:
    """生成 phase_coverage.json，覆盖全部 17 个吸收元素。"""
    pack_root = Path(pack_root)
    index = index or build_index(pack_root)
    motif_path = Path(motif_coverage_path or MOTIF_COVERAGE_PATH)
    motif = {}
    if motif_path.is_file():
        try:
            motif = json.loads(motif_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            motif = {}
    motif_elements = dict(motif.get("elements") or {})

    rows_by_element: dict[str, list[dict]] = {
        el: [] for el in ABSORBER_ELEMENTS}
    for row in index.get("entries", []):
        for el in row.get("absorber_elements", []):
            rows_by_element.setdefault(el, []).append(row)

    elements = {}
    for el in ABSORBER_ELEMENTS:
        rows = rows_by_element.get(el, [])
        motif_row = motif_elements.get(el) or {}
        verified = len(rows)
        motif_total = int(motif_row.get("motifs") or 0)
        missing = [phase for phase in COMMON_PHASES.get(el, [])
                   if not any(_phase_tag_match(row, phase) for row in rows)]
        if verified:
            status = "verified"
        elif motif_total:
            status = "motif_only"
        else:
            status = "not_covered"
        elements[el] = {
            "element": el,
            "verified_cif": verified,
            "prototype": (int(motif_row.get("generated_only") or 0)
                          + int(motif_row.get("unverified") or 0)
                          + int(motif_row.get("source_ref_only") or 0)),
            "cluster": int(motif_row.get("constructed") or 0),
            "motif_total": motif_total,
            "motif_source_ref_only": int(motif_row.get("source_ref_only")
                                          or 0),
            "formal_cn_ready": verified > 0,
            "coverage_status": status,
            "phase_ids": sorted(f"{row['source']}:{row['id']}" for row in rows),
            "common_phases": list(COMMON_PHASES.get(el, [])),
            "missing_common_phases": missing,
        }

    verified_elements = [el for el in ABSORBER_ELEMENTS
                         if elements[el]["formal_cn_ready"]]
    not_covered = [el for el in ABSORBER_ELEMENTS
                   if elements[el]["coverage_status"] == "not_covered"]
    motif_only = [el for el in ABSORBER_ELEMENTS
                  if elements[el]["coverage_status"] == "motif_only"]
    with_motif = [el for el in ABSORBER_ELEMENTS
                  if elements[el]["motif_total"] > 0]
    without_motif = [el for el in ABSORBER_ELEMENTS
                     if elements[el]["motif_total"] == 0]
    coverage = {
        "schema_version": 1,
        "generated_by": "tools/build_phase_index.py",
        "note": ("覆盖只描述本地来源；没有 verified CIF 的元素正式 CN "
                 "只能走 --phase-cif，不得回退 motif。"),
        "absorber_elements": list(ABSORBER_ELEMENTS),
        "summary": {
            "elements_total": len(ABSORBER_ELEMENTS),
            "elements_verified": len(verified_elements),
            "verified_elements": verified_elements,
            "elements_motif_only": len(motif_only),
            "motif_only_elements": motif_only,
            "elements_not_covered": len(not_covered),
            "not_covered_elements": not_covered,
            "elements_with_motif": len(with_motif),
            "with_motif_elements": with_motif,
            "elements_without_motif": len(without_motif),
            "without_motif_elements": without_motif,
            "verified_cif": len(index.get("entries", [])),
        },
        "elements": elements,
    }
    if out_path is not None:
        _write_json_if_changed(Path(out_path), coverage)
    return coverage


def verified_entries(index: dict, element: str) -> list[dict]:
    """返回某吸收元素的 verified 条目；只认 experimental + crystal_*。"""
    el = normalize_element(element)
    return [row for row in index.get("entries", [])
            if el in row.get("absorber_elements", [])
            and row.get("experimental") is True
            and row.get("model_level") in VERIFIED_LEVELS]
