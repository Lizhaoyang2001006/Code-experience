# -*- coding: utf-8 -*-
"""
phase_resolver.py — P1b 相解析与缓存（正式 CN 结构来源解析器）
================================================================

把用户指定的相（CIF 文件、精确数据库编号、SPEC 字符串、本地索引）
解析为可用于正式 CN 的 ``StructureModel``（verified_source_cif），并
将实际用到的条目缓存到本地，供后续复用。

解析优先级（严格按此顺序，不自动 Top-1）：

  1. 用户 ``--phase-cif``（可重复）—— 最高优先级，用户 CIF 始终优先；
  2. 精确 ``--phase-id SOURCE:ID`` —— 本地缓存/phase pack 优先；
     未命中且非离线才允许 fetcher（CI 默认 fetcher=None）；
  3. 本地 ``phase_index.json`` 唯一命中；
  4. ``--phase SPEC`` 唯一命中；
  5. 多候选：返回候选表，cn_table 0 行 + cif_ambiguous，不自动选；
  6. 离线未命中：offline_miss，不退回 motif。

安全约束：

  - ``--phase SPEC`` 只用 ``split(";")`` / ``split("=", 1)`` 解析；
    禁止 ``eval``/``exec``；
  - ``cell`` 必须解析为 6 个有限浮点数，非法抛 ``XAFSConfigError``；
  - 缓存复用前校验 sha256；不匹配则拒绝并重新获取；
  - 网络 fetcher 必须依赖注入；CI 默认 ``fetcher=None + --offline``，
    测试用 stub，不发真实网络请求；
  - 不做 COD/ICSD 全库镜像；只缓存实际用到的条目；
  - ``structure_lib._cif_path`` 现场生成 CIF 拒绝入库。
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path

import absorber as absorber_mod
import phase_pack
from motif_provenance import is_generated_cif, sha256_file
from structure_model import StructureModel
from xafs_core import XAFSConfigError

ROOT = Path(__file__).resolve().parent
DEFAULT_CACHE_ROOT = ROOT / "_structure_lib" / "_phase_cache"
DEFAULT_PACK_ROOT = ROOT / "_structure_lib" / "_phase_pack"
DEFAULT_INDEX_PATH = ROOT / "_structure_lib" / "phase_index.json"

ALLOWED_SPEC_KEYS = {"name", "formula", "spacegroup", "cell", "source", "id"}
NUMERIC_CELL_KEYS = ("a", "b", "c", "alpha", "beta", "gamma")


def _normalize_element(symbol: str) -> str:
    """把 'Fe3+' / 'FE' / 'fe' 归一为 'Fe'；未知则返回 capitalized。"""
    s = str(symbol).strip()
    m = re.match(r"([A-Z][a-z]?)", s)
    return m.group(1) if m else s.capitalize()

ROOT = Path(__file__).resolve().parent
DEFAULT_CACHE_ROOT = ROOT / "_structure_lib" / "_phase_cache"
DEFAULT_PACK_ROOT = ROOT / "_structure_lib" / "_phase_pack"
DEFAULT_INDEX_PATH = ROOT / "_structure_lib" / "phase_index.json"

ALLOWED_SPEC_KEYS = {"name", "formula", "spacegroup", "cell", "source", "id"}
NUMERIC_CELL_KEYS = ("a", "b", "c", "alpha", "beta", "gamma")


class PhaseResolutionError(XAFSConfigError):
    """相解析拒绝错误；携带 code 供测试与报告引用。"""

    def __init__(self, code: str, message: str):
        self.code = str(code)
        self.detail = str(message)
        super().__init__("phase_resolver", f"{self.code}: {self.detail}")


class PhaseAmbiguousError(PhaseResolutionError):
    """多候选，不自动选；调用方应把 candidates 填入 cn_table 候选表。"""


class PhaseOfflineMiss(PhaseResolutionError):
    """离线未命中，不退回 motif。"""


# ---------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------
@dataclass
class PhaseCandidate:
    """一个相候选：源信息 + CIF 路径 + 是否可走正式 CN。"""

    source: str
    id: str
    cif_path: str          # 绝对路径，供 validate_formal 直接使用
    cif_sha256: str
    absorber_elements: list[str]
    absorber_sites: list[str]
    formula: str
    spacegroup: str
    cell: dict
    experimental: bool
    model_level: str
    from_cache: bool = False
    original_filename: str | None = None
    sites: list[dict] = field(default_factory=list)   # 完整位点
    score: float | None = None


@dataclass
class PhaseResolution:
    """解析结果。"""

    absorber: str
    absorber_site: str | None
    selected: PhaseCandidate | None
    candidates: list[PhaseCandidate] = field(default_factory=list)
    status: str = "no_match"   # resolved / ambiguous / offline_miss / no_match
    structure_model: StructureModel | None = None
    reason: str = ""


# ---------------------------------------------------------------
# SPEC 解析（仅字符串拆分，禁止 eval/exec）
# ---------------------------------------------------------------
def parse_phase_spec(spec: str) -> dict:
    """解析 ``--phase SPEC``，返回 ``{name, formula, spacegroup, cell,
    source, id}``（仅包含出现的键）。

    格式示例::

        name=fcc metal;formula=Fe;spacegroup=Im-3m;
        cell=a=2.868,b=2.868,c=2.868,alpha=90,beta=90,gamma=90;
        source=COD;id=1100108

    只做字符串拆分（split(";") / split("=", 1)），绝不 eval/exec。
    cell 必须解析为 6 个正有限浮点数。
    """
    result: dict = {}
    for pair in str(spec).split(";"):
        pair = pair.strip()
        if not pair:
            continue
        if "=" not in pair:
            raise PhaseResolutionError(
                "invalid_spec", f"SPEC 片段缺少 '=': {pair!r}")
        key, value = pair.split("=", 1)
        key = key.strip().lower()
        value = value.strip()
        if key not in ALLOWED_SPEC_KEYS:
            raise PhaseResolutionError(
                "invalid_spec_key", f"SPEC 不允许的键: {key!r}")
        if key in result:
            raise PhaseResolutionError(
                "invalid_spec", f"SPEC 键重复: {key!r}")
        result[key] = value
    if "cell" in result:
        result["cell"] = _parse_cell(result["cell"])
    return result


def _parse_cell(cell_text: str) -> dict:
    text = str(cell_text).strip()
    parts = [p.strip() for p in text.split(",")]
    if len(parts) != 6:
        raise PhaseResolutionError(
            "invalid_cell",
            f"cell 需要 6 个数字 (a,b,c,alpha,beta,gamma)，收到 {len(parts)} 个: {cell_text!r}")
    values = []
    for p in parts:
        # 支持 "a=2.868" 或 "2.868" 两种写法
        if "=" in p:
            p = p.split("=", 1)[1].strip()
        try:
            values.append(float(p))
        except ValueError as exc:
            raise PhaseResolutionError(
                "invalid_cell", f"cell 包含非数字: {cell_text!r}") from exc
    for v in values:
        if not math.isfinite(v) or v <= 0:
            raise PhaseResolutionError(
                "invalid_cell", f"cell 必须是正有限数: {values}")
    return {k: v for k, v in zip(NUMERIC_CELL_KEYS, values)}


# ---------------------------------------------------------------
# 候选构建
# ---------------------------------------------------------------
def _entry_to_candidate(base_dir: Path, entry: dict) -> PhaseCandidate:
    """从 manifest entry + 基目录构建 PhaseCandidate（绝对 cif_path）。

    entry["cif_path"] 可能是 pack-relative（如 COD/.../sha.cif，来自
    manifest.json）或 repo-relative（如 _structure_lib/_phase_pack/...，
    来自 phase_index.json），两者皆能正确解析。
    """
    cif_rel = Path(str(entry["cif_path"]))
    if cif_rel.is_absolute():
        cif_path = cif_rel.resolve()
    else:
        # 先按 base_dir 解析；若不存在，再按 repo root 解析（兼容
        # phase_index.json 里的 repo-relative 路径）。
        cif_path = (base_dir / cif_rel).resolve()
        if not cif_path.is_file():
            cif_path = (ROOT / cif_rel).resolve()
    parsed = phase_pack.parse_cif(cif_path)
    return PhaseCandidate(
        source=str(entry["source"]),
        id=str(entry["id"]),
        cif_path=str(cif_path),
        cif_sha256=str(entry["cif_sha256"]).lower(),
        absorber_elements=list(entry["absorber_elements"]),
        absorber_sites=list(entry["absorber_sites"]),
        formula=str(entry["formula"]),
        spacegroup=str(entry["spacegroup"]),
        cell=dict(entry["cell"]),
        experimental=bool(entry.get("experimental", False)),
        model_level=str(entry["model_level"]),
        from_cache=False,
        original_filename=entry.get("original_filename"),
        sites=parsed["sites"],
    )


def _cif_to_candidate(cif_path: Path, source: str, motif_id: str) -> PhaseCandidate:
    """从单个用户 CIF 文件构建候选（source="user"）。"""
    if is_generated_cif(cif_path):
        raise PhaseResolutionError(
            "generated_cif_not_allowed",
            f"{cif_path.name}: 现场生成 CIF 不得作为正式来源")
    parsed = phase_pack.parse_cif(cif_path)
    sha = sha256_file(cif_path)
    return PhaseCandidate(
        source=source,
        id=motif_id,
        cif_path=str(cif_path.resolve()),
        cif_sha256=sha,
        absorber_elements=parsed["absorber_elements"],
        absorber_sites=list(parsed["default_absorber_sites"]),
        formula=parsed["formula"],
        spacegroup=parsed["spacegroup"],
        cell=parsed["cell"],
        experimental=True,
        model_level="crystal_cif",
        from_cache=False,
        original_filename=cif_path.name,
        sites=parsed["sites"],
    )


# ---------------------------------------------------------------
# 缓存（_structure_lib/_phase_cache/<source>/<id>/<sha256>.cif + .meta.json）
# ---------------------------------------------------------------
def _cache_dir(cache_root: Path, source: str, motif_id: str) -> Path:
    return cache_root / source / motif_id


def _cache_meta_path(cache_root: Path, source: str, motif_id: str,
                     sha256: str) -> Path:
    return _cache_dir(cache_root, source, motif_id) / f"{sha256}.meta.json"


def _cache_cif_path(cache_root: Path, source: str, motif_id: str,
                    sha256: str) -> Path:
    return _cache_dir(cache_root, source, motif_id) / f"{sha256}.cif"


def _cache_list(cache_root: Path, source: str, motif_id: str) -> list[tuple]:
    """列出缓存目录下全部可用条目：(sha256, entry, cif_path)。"""
    out = []
    d = _cache_dir(cache_root, source, motif_id)
    if not d.is_dir():
        return out
    for meta_path in sorted(d.glob("*.meta.json")):
        try:
            entry = json.loads(meta_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        cif_path = meta_path.with_name(
            meta_path.name[: -len(".meta.json")] + ".cif")
        if not cif_path.is_file():
            continue
        try:
            actual_sha = sha256_file(cif_path)
        except OSError:
            continue
        entry_sha = str(entry.get("cif_sha256", "")).lower()
        if entry_sha and actual_sha != entry_sha:
            continue
        out.append((actual_sha, entry, cif_path.resolve()))
    return out


def _store_to_cache(cache_root: Path, source: str, motif_id: str,
                    src_cif: Path, sha256: str, entry: dict) -> None:
    """把一条已验证的 CIF 存入缓存（幂等）。"""
    d = _cache_dir(cache_root, source, motif_id)
    d.mkdir(parents=True, exist_ok=True)
    target_cif = _cache_cif_path(cache_root, source, motif_id, sha256)
    target_meta = _cache_meta_path(cache_root, source, motif_id, sha256)
    if not target_cif.exists():
        shutil.copy2(src_cif, target_cif)
    cache_entry = dict(entry)
    cache_entry["cif_path"] = target_cif.relative_to(
        cache_root).as_posix() if cache_root in target_cif.parents \
        else str(target_cif)
    cache_entry["cif_sha256"] = sha256
    cache_entry["added_by"] = "phase_resolver cache"
    phase_pack._write_json_if_changed(target_meta, cache_entry)


# ---------------------------------------------------------------
# 结构模型构建
# ---------------------------------------------------------------
def _candidate_to_structure_model(
        candidate: PhaseCandidate, absorber: str,
        absorber_site: str | None) -> StructureModel:
    """把一个 verified 候选构造为 StructureModel（走严格链）。

    单吸收位点可以自动使用；多吸收位点必须显式 ``--absorber-site``，
    不自动 Top-1，避免把 Fe3O4 的四面体位/八面体位混为一谈。
    """
    sites = list(candidate.absorber_sites or [])
    if absorber_site:
        canonical = {s.lower(): s for s in sites}.get(
            str(absorber_site).strip().lower())
        if canonical is None:
            raise PhaseResolutionError(
                "absorber_site_missing",
                f"{candidate.source}:{candidate.id} 没有吸收位 "
                f"{absorber_site!r}；可用位点: {sites}")
        site = canonical
    elif len(sites) == 1:
        site = sites[0]
    elif len(sites) > 1:
        raise PhaseResolutionError(
            "absorber_site_ambiguous",
            f"{candidate.source}:{candidate.id} 有多个吸收位点"
            f"（{len(sites)} 个）{sites}，不自动 Top-1；"
            f"请用 --absorber-site 显式指定")
    else:
        raise PhaseResolutionError(
            "absorber_site_missing",
            f"{candidate.source}:{candidate.id} 缺少吸收位；"
            f"请用 --absorber-site 显式指定")

    source_type = "user_cif" if candidate.source == "user" else "database_cif"
    source_ref = f"{candidate.source}:{candidate.id}"

    # 完整 asym_sites（从缓存的 CIF 位点；不含占据度等细节时留空）
    asym_sites = []
    for s in candidate.sites:
        asym_sites.append({
            "label": s["label"],
            "el": s["element"],
            "x": s["fract"][0],
            "y": s["fract"][1],
            "z": s["fract"][2],
            "occ": s["occ"],
        })

    model = StructureModel(
        model_id=f"{candidate.source}:{candidate.id}",
        absorber=absorber,
        absorber_site=site,
        model_level=candidate.model_level,
        source_type=source_type,
        source_ref=source_ref,
        source_cif=candidate.cif_path,
        source_cif_sha256=candidate.cif_sha256,
        cell=candidate.cell,
        spacegroup=candidate.spacegroup,
        asym_sites=asym_sites or None,
        cif_path=candidate.cif_path,
        occupancy_policy="explicit",
        edits=[],
        assumptions={
            "p1b": "phase_resolver resolved; verified_source_cif; "
                   "formal CN gate applies",
        },
        n0_shells=[],
        warnings=[],
        provenance_status="verified_source_cif",
    )
    return model


# ---------------------------------------------------------------
# 候选评分（SPEC 匹配）
# ---------------------------------------------------------------
def _spec_match_score(candidate: PhaseCandidate, spec: dict) -> float:
    """按 SPEC 字段给候选打分；越高越匹配。不匹配返回 -1。"""
    score = 0.0
    matched = 0
    for key, value in spec.items():
        if key == "source":
            if str(candidate.source).lower() != str(value).lower():
                return -1.0
            matched += 1
            score += 1.0
        elif key == "id":
            if str(candidate.id) != str(value):
                return -1.0
            matched += 1
            score += 1.0
        elif key == "formula":
            if str(candidate.formula).lower() != str(value).lower():
                return -1.0
            matched += 1
            score += 1.0
        elif key == "spacegroup":
            if "".join(str(candidate.spacegroup).split()).lower() != \
                    "".join(str(value).split()).lower():
                return -1.0
            matched += 1
            score += 1.0
        elif key == "cell":
            try:
                spec_cell = _parse_cell(str(value))
            except PhaseResolutionError:
                return -1.0
            ok = True
            for k in NUMERIC_CELL_KEYS:
                cv = float(candidate.cell.get(k, 0) or 0)
                sv = float(spec_cell.get(k, 0) or 0)
                if sv <= 0 or abs(cv - sv) > max(1e-3, 1e-4 * abs(sv)):
                    ok = False
                    break
            if not ok:
                return -1.0
            matched += 1
            score += 1.0
    if matched == 0:
        return -1.0
    return score / len(spec)


def _filter_by_absorber(
        candidates: list[PhaseCandidate], absorber: str) -> list[PhaseCandidate]:
    """只保留吸收元素包含 absorber 的候选。"""
    return [c for c in candidates if absorber in c.absorber_elements]


# ---------------------------------------------------------------
# 主解析入口
# ---------------------------------------------------------------
def resolve_phase(
    absorber: str,
    *,
    phase_cif: list[str] | None = None,
    phase_id: str | None = None,
    phase_spec: str | None = None,
    phase_index_path: str | Path | None = None,
    phase_pack_path: str | Path | None = None,
    absorber_site: str | None = None,
    offline: bool = False,
    fetcher=None,
) -> PhaseResolution:
    """解析相来源，返回 PhaseResolution。

    ``fetcher`` 是依赖注入的回调：``fetcher(source, id) -> Path | None``，
    返回下载/获取的 CIF 本地路径，None 表示未命中。CI 默认传 None。
    """
    absorber = _normalize_element(absorber)
    pack_root = Path(phase_pack_path or DEFAULT_PACK_ROOT)
    cache_root = DEFAULT_CACHE_ROOT
    index_path = Path(phase_index_path or DEFAULT_INDEX_PATH)

    candidates: list[PhaseCandidate] = []

    # ---- 1~4 按优先级处理，一旦某阶段产生候选即停止（不自动回退）----
    if phase_cif:
        # Stage 1: user CIF (highest priority)
        for cif_str in phase_cif:
            cif_path = Path(cif_str)
            if not cif_path.is_file():
                raise PhaseResolutionError(
                    "cif_missing", f"--phase-cif file not found: {cif_path}")
            candidates.append(_cif_to_candidate(cif_path, "user",
                                                cif_path.stem))

    elif phase_id:
        # Stage 2: exact ID SOURCE:ID
        m = re.match(r"^([^:]+):(.+)$", str(phase_id))
        if not m:
            raise PhaseResolutionError(
                "invalid_phase_id",
                f"--phase-id format should be SOURCE:ID, got {phase_id!r}")
        src, mid = m.group(1), m.group(2)

        # 2a. local cache first
        cached = _cache_list(cache_root, src, mid)
        if cached:
            sha, entry, cif_path = cached[0]
            candidates.append(_entry_to_candidate(
                Path(str(cif_path.parent)),
                {**entry, "cif_path": cif_path.name}))
        else:
            # 2b. then phase pack
            if pack_root.is_dir():
                try:
                    manifest = phase_pack.load_manifest(pack_root)
                except PhaseResolutionError:
                    manifest = {"entries": []}
                pack_entries = [
                    e for e in manifest.get("entries", [])
                    if str(e.get("source", "")).lower() == src.lower()
                    and str(e.get("id", "")) == mid]
                if pack_entries:
                    entry = pack_entries[0]
                    result = phase_pack.validate_entry(pack_root, entry)
                    if result["ok"]:
                        cif_path = (pack_root / entry["cif_path"]).resolve()
                        sha = str(entry["cif_sha256"]).lower()
                        _store_to_cache(cache_root, src, mid,
                                        cif_path, sha, entry)
                        candidates.append(_entry_to_candidate(
                            pack_root, entry))
            # 2c. still unmatched: allow fetcher (non-offline only)
            if not any(c.source == src and c.id == mid for c in candidates):
                if not offline and fetcher is not None:
                    try:
                        fetched = fetcher(src, mid)
                    except Exception:
                        fetched = None
                    if fetched:
                        fetched_path = Path(fetched)
                        if is_generated_cif(fetched_path):
                            raise PhaseResolutionError(
                                "generated_cif_not_allowed",
                                f"fetcher returned generated CIF: {fetched_path}")
                        parsed = phase_pack.parse_cif(fetched_path)
                        sha = sha256_file(fetched_path)
                        entry = {
                            "source": src, "id": mid,
                            "formula": parsed["formula"],
                            "spacegroup": parsed["spacegroup"],
                            "cell": parsed["cell"],
                            "absorber_sites": parsed["default_absorber_sites"],
                            "absorber_elements": parsed["absorber_elements"],
                            "cif_path": str(fetched_path),
                            "cif_sha256": sha,
                            "experimental": True,
                            "model_level": "crystal_cif",
                        }
                        _store_to_cache(cache_root, src, mid,
                                        fetched_path, sha, entry)
                        candidates.append(_cif_to_candidate(
                            fetched_path, src, mid))

    elif index_path.is_file() or phase_spec is not None:
        # Stage 3: local phase_index.json unique match;
        # Stage 4: --phase SPEC unique match.
        try:
            index_data = json.loads(
                index_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            if phase_spec is not None and not index_path.is_file():
                index_data = {}
            else:
                raise PhaseResolutionError(
                    "invalid_phase_index",
                    f"{index_path}: {exc}") from exc
        index_entries = index_data.get("entries", []) if isinstance(
            index_data, dict) else []
        pack_for_index = Path(
            str(index_data.get("pack_root", pack_root))) if isinstance(
                index_data, dict) else pack_root

        if phase_spec is not None:
            # Stage 4: SPEC match (score, pick highest; tie => ambiguous)
            spec = parse_phase_spec(phase_spec)
            best: list[tuple[float, PhaseCandidate]] = []
            for entry in index_entries:
                if absorber not in entry.get("absorber_elements", []):
                    continue
                cand = _entry_to_candidate(pack_for_index, entry)
                score = _spec_match_score(cand, spec)
                if score >= 0:
                    best.append((score, cand))
            if best:
                best.sort(key=lambda x: x[0], reverse=True)
                top_score = best[0][0]
                candidates = [c for s, c in best if s == top_score]
        else:
            # Stage 3: match all entries by absorber element
            for entry in index_entries:
                if absorber in entry.get("absorber_elements", []):
                    candidates.append(_entry_to_candidate(
                        pack_for_index, entry))

    # ---- 去重（按 source:id+sha256）-------------------------------
    seen = set()
    uniq: list[PhaseCandidate] = []
    for c in candidates:
        key = (c.source, c.id, c.cif_sha256)
        if key not in seen:
            seen.add(key)
            uniq.append(c)
    candidates = uniq

    # ---- 只保留含该吸收元素的候选 --------------------------------
    candidates = _filter_by_absorber(candidates, absorber)

    # ---- 决定状态 -------------------------------------------------
    resolution = PhaseResolution(
        absorber=absorber,
        absorber_site=absorber_site,
        selected=None,
        candidates=candidates,
        status="no_match",
        reason="",
        structure_model=None,
    )

    if not candidates:
        resolution.status = "offline_miss" if offline else "no_match"
        resolution.reason = (
            f"未找到吸收元素 {absorber} 的相来源"
            + ("（离线模式）" if offline else ""))
        return resolution

    if len(candidates) > 1:
        resolution.status = "ambiguous"
        resolution.reason = (
            f"吸收元素 {absorber} 有 {len(candidates)} 个相候选，"
            f"不自动 Top-1；请用 --phase-id 精确指定或 --absorber-site 消歧")
        return resolution

    # ---- 唯一命中 --------------------------------------------------
    cand = candidates[0]
    resolution.selected = cand
    resolution.status = "resolved"
    resolution.reason = "唯一命中"

    # 只有 verified（experimental + crystal_*）才能构造 StructureModel
    if (cand.experimental is True
            and cand.model_level in ("crystal_cif", "crystal_edited")):
        try:
            resolution.structure_model = _candidate_to_structure_model(
                cand, absorber, absorber_site)
        except PhaseResolutionError as exc:
            # 缺 absorber_site / 多位点未指定：保留 resolved，但
            # structure_model 为 None；调用方通过 formal_cn 门禁会得到
            # 明确原因，而不是静默使用第一个位点。
            resolution.structure_model = None
            resolution.reason = exc.detail
    else:
        resolution.structure_model = None
        resolution.reason = (
            f"候选不可作为正式来源（experimental={cand.experimental}, "
            f"model_level={cand.model_level}）")

    return resolution


# ---------------------------------------------------------------
# 向后兼容的别名（供 tests / CLI 调用）
# ---------------------------------------------------------------
def resolve_structure_source(*args, **kwargs) -> PhaseResolution:
    """resolve_phase 的别名，供外部调用。"""
    return resolve_phase(*args, **kwargs)
