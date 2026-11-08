# -*- coding: utf-8 -*-
"""motif_provenance.py — motif 来源分层（P0b）

正式 CN 只接受 ``verified_source_cif``。本模块只做**分类**，不计算 N₀，
也不改变任何诊断拟合路径。

分层规则（与 REVIEW_AND_REVISED_PLAN.md §12.6 一致）：

===============  ==========================  ==================
条件              provenance_status           model_level
===============  ==========================  ==================
真实 CIF+sha+位点  verified_source_cif         crystal_cif/edited
只有 COD/AMCSD 号  source_ref_only             prototype
现场生成 CIF       generated_only              prototype
cluster/mixture    constructed                 cluster
其他               unverified                  prototype/cluster
===============  ==========================  ==================

关键约束：``_structure_lib/<元素>/_cif/<id>.cif`` 是 ``_cif_path()``
现场拼出来的最小 CIF，**永远不能**作为 ``source_cif_sha256`` 的来源。
只有通过 ``Motif.cif_path`` 显式引用的外部真实 CIF 才可能升级。
"""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

PROVENANCE_STATUSES = (
    "verified_source_cif",
    "source_ref_only",
    "generated_only",
    "constructed",
    "unverified",
)

MODEL_LEVELS = (
    "crystal_cif",
    "crystal_edited",
    "prototype",
    "cluster",
)

OCCUPANCY_SOURCES = (
    "explicit",
    "implicit_1.0",
    "generated",
    "unknown",
)

# 这些字段会写回 motif JSON；``_source_cif`` 只在本模块内部使用。
PERSISTED_FIELDS = (
    "source_ref",
    "source_cif_sha256",
    "absorber_site",
    "occupancy_source",
    "model_level",
    "provenance_status",
    "edits",
)

CRYSTAL_SOURCE_TYPES = {
    "user_cif", "database_cif", "cod", "amcsd", "cif", "experimental_cif",
}
_CONSTRUCTED_TYPES = {
    "constructed", "substitution", "composite", "alloy", "sac", "dac",
    "cluster", "manual_model",
}
_DB_REF_RE = re.compile(r"\b(COD|AMCSD)\s*[-_:# ]?\s*(\d{4,})\b", re.I)
_GENERATED_MARKER = "generated_by: structure_lib._cif_path"


def sha256_file(path: str | Path) -> str:
    """流式计算文件 SHA256，避免把大 CIF 一次性读入内存。"""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_generated_cif(path: str | Path | None) -> bool:
    """判断 CIF 是否由 ``structure_lib._cif_path`` 现场生成。

    两条独立证据，任一命中即可：

    1. 路径位于 ``_structure_lib/<元素>/_cif/``；
    2. 文件头含 ``generated_by: structure_lib._cif_path`` 标记。

    路径判断覆盖库内全部生成文件；标记判断只用于人工标注过、
    但被复制到别处的生成文件。正式门禁不依赖 CIF 内容猜测。
    """
    if not path:
        return False
    p = Path(path)
    try:
        parts = [str(s).lower() for s in p.resolve().parts]
    except OSError:
        parts = [str(s).lower() for s in p.parts]
    if "_structure_lib" in parts:
        i = parts.index("_structure_lib")
        if i + 2 < len(parts) and parts[i + 2] == "_cif":
            return True
    try:
        head = p.read_text(encoding="utf-8", errors="ignore")[:512]
    except OSError:
        return False
    return _GENERATED_MARKER in head


def _fmt(value) -> str:
    try:
        return f"{float(value):.6g}"
    except (TypeError, ValueError):
        return str(value)


def _source_ref(motif) -> str | None:
    """提取稳定来源编号；只认 COD/AMCSD 或显式 source_ref。"""
    explicit = getattr(motif, "source_ref", None)
    if explicit:
        return str(explicit)
    src = dict(getattr(motif, "source", None) or {})
    src_type = str(src.get("type") or "").lower()
    for value in (src.get("ref"), src.get("id")):
        if value is None:
            continue
        m = _DB_REF_RE.search(str(value))
        if m:
            return f"{m.group(1).upper()}:{m.group(2)}"
    if src_type in ("cod", "amcsd") and src.get("id") is not None:
        return f"{src_type.upper()}:{src['id']}"
    return None


def absorber_site_label(motif) -> str | None:
    """返回吸收位标签；优先显式字段，其次分数坐标，最后唯一吸收位。"""
    explicit = getattr(motif, "absorber_site", None)
    if explicit:
        return str(explicit)
    frac = getattr(motif, "absorber_fract", None)
    if frac is not None and len(frac) == 3:
        return "frac:(" + ",".join(_fmt(v) for v in frac) + ")"
    if str(getattr(motif, "kind", "")) != "crystal":
        return None
    absorber = str(getattr(motif, "absorber", "") or "").capitalize()
    hits = [
        s for s in (getattr(motif, "asym_sites", None) or [])
        if str(s.get("el", "")).capitalize() == absorber
    ]
    if len(hits) == 1:
        s = hits[0]
        xyz = (s.get("x"), s.get("y"), s.get("z"))
        if all(v is not None for v in xyz):
            label = s.get("label") or s.get("wyckoff") or absorber
            return f"{label}@({_fmt(xyz[0])},{_fmt(xyz[1])},{_fmt(xyz[2])})"
    return None


def _cif_has_occupancy_tag(path: Path) -> bool:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return False
    return "_atom_site_occupancy" in text


def _occupancy_source(motif, real_cif: Path | None) -> str:
    explicit = getattr(motif, "occupancy_source", None)
    if explicit:
        return str(explicit)
    if real_cif is not None:
        return ("explicit" if _cif_has_occupancy_tag(real_cif)
                else "implicit_1.0")
    occs = [
        float(s.get("occ", 1.0))
        for s in (getattr(motif, "asym_sites", None) or [])
        if isinstance(s, dict) and s.get("occ") is not None
    ]
    if any(o < 0.999 for o in occs):
        return "explicit"
    return "generated"


def classify_motif(motif) -> dict:
    """分类一条 motif；返回可直接写回 JSON 的字段 + 内部 ``_source_cif``。

    该函数是纯函数，不写文件、不改 motif；重复调用结果一致。
    """
    src = dict(getattr(motif, "source", None) or {})
    src_type = str(src.get("type") or "").lower()
    kind = str(getattr(motif, "kind", "") or "")
    cif_path = getattr(motif, "cif_path", None)

    real_cif: Path | None = None
    has_generated_cif = False
    if cif_path:
        candidate = Path(cif_path)
        if candidate.exists():
            if is_generated_cif(candidate):
                has_generated_cif = True
            else:
                real_cif = candidate

    sha = None
    if real_cif is not None:
        sha = (getattr(motif, "source_cif_sha256", None)
               or sha256_file(real_cif))

    site = absorber_site_label(motif)
    source_ref = _source_ref(motif)
    occupancy_source = _occupancy_source(motif, real_cif)
    edits = list(getattr(motif, "edits", None) or [])

    if real_cif is not None and sha and site:
        status = "verified_source_cif"
        current_level = getattr(motif, "model_level", None)
        if current_level in ("crystal_cif", "crystal_edited"):
            level = current_level
        elif edits:
            level = "crystal_edited"
        else:
            level = "crystal_cif"
    elif kind in ("cluster", "mixture") or src_type in _CONSTRUCTED_TYPES:
        status = "constructed"
        level = "cluster"
    elif real_cif is not None:
        # 有真实 CIF，但缺吸收位/哈希，不能自动升级。
        status = "unverified"
        level = "prototype"
    elif source_ref:
        status = "source_ref_only"
        level = "prototype"
    elif kind == "crystal":
        status = "generated_only"
        level = "prototype"
    else:
        status = "unverified"
        level = "prototype"

    return {
        "_source_cif": str(real_cif) if real_cif is not None else None,
        "_has_generated_cif": has_generated_cif,
        "source_ref": source_ref,
        "source_cif_sha256": sha,
        "absorber_site": site,
        "occupancy_source": occupancy_source,
        "model_level": level,
        "provenance_status": status,
        "edits": edits,
    }
