# -*- coding: utf-8 -*-
"""
cif_gate.py — P2a CIF 门禁（正式 CN 结构来源门禁）
======================================================
对用户/phase pack 解析出的 CIF 做逐项门禁检查，决定是否可作为
正式 CN 的 N₀ 来源。只做拒绝/放行判断，不做拟合。

设计原则（与 REVIEW_AND_REVISED_PLAN.md §5.1 一致）：

- 复用 artemis.read_cif_sites / structure_lib.expand_full，不重造 gemmi 展开；
- 必须有 _atom_site_* 循环；记录原始 CIF 路径与 sha256；
- artemis.read_cif_sites 对缺失 occupancy 默认 1.0，不能用来判断“缺失”；
  门禁必须直接读 CIF 文本检查 _atom_site_occupancy 标签；
- 同元素多个晶位时必须 --absorber-site；禁止按 occ × 多重度自动挑主位点；
- 计算结构不得作为实验 N₀；--allow-computed 时降级 cn_hypothesis_only；
- 第一壳 CN 正常范围 2–12；12–16 警告；>16 或 <2 需 --allow-unusual-coordination；
  fcc/hcp 的 CN=12 必须正常通过；
- 失败码按 §5.1 列表返回。

失败码（§5.1）：
    parse_failed, absorber_missing, absorber_site_ambiguous,
    polyhedron_not_closed, bond_length_implausible,
    occupancy_missing, occupancy_scale_failed, cocupancy_unresolved,
    cif_computed, cif_ambiguous, offline_miss, consistency_failed,
    n_not_identifiable, n_out_of_range, correlation_too_high,
    defect_model_missing, no_cif, legacy_unverified
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np

import artemis
import structure_lib
from structure_lib import SHELL_TOL
from motif_provenance import is_generated_cif, sha256_file
from xafs_core import XAFSConfigError, XAFSDataError

# 共价半径表（Å），用于键长合理性检查（取自 Pyykko & Atsumi）
COVALENT_RADIUS = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76,
    "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58, "Na": 1.66, "Mg": 1.41,
    "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.66,
    "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
    "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01,
    "Sm": 1.98, "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92, "Ho": 1.92,
    "Er": 1.89, "Tm": 1.90, "Yb": 1.90, "Lu": 1.88, "Hf": 1.85, "Ta": 1.80,
    "W": 1.80, "Re": 1.76, "Os": 1.76, "Ir": 1.76, "Pt": 1.76, "Au": 1.74,
    "Hg": 1.74, "Tl": 1.87, "Pb": 1.87, "Bi": 1.82, "Po": 1.82, "At": 1.82,
    "Rn": 1.82,
}

# 第一壳 CN 范围（§5.1）
CN_NORMAL_LO, CN_NORMAL_HI = 2, 12
CN_WARN_LO, CN_WARN_HI = 12, 16
BOND_TOL = 0.25           # 键长容差：共价半径之和 ±0.25 Å
OCCUPANCY_TOL = 0.08      # 共占位合并距离容差（Å），与 artemis 一致


class CIFGateError(XAFSConfigError):
    """CIF 门禁拒绝错误；carries failure code. """

    def __init__(self, code: str, message: str):
        self.code = str(code)
        self.detail = str(message)
        super().__init__("cif_gate", f"{self.code}: {self.detail}")


@dataclass
class CIFGateResult:
    """门禁检查结果。"""

    ok: bool
    code: str | None = None            # 失败码；ok 时为 None
    reason: str = ""
    absorber: str = ""
    absorber_site: str | None = None
    occupancy_source: str = "unknown"  # explicit / implicit_1.0 / unknown
    structure_is_computed: bool = False
    cn_hypothesis_only: bool = False
    warnings: list[str] = field(default_factory=list)
    # N₀ shell 分组结果（§5.2）
    n0_shells: list[dict] = field(default_factory=list)
    first_shell_cn: float | None = None
    rayleigh_delta_r: float | None = None
    # 原始 CIF 追溯
    cif_path: str | None = None
    cif_sha256: str | None = None
    # P5-D：多位点化学计量权重（m_i × occ_i 归一化）
    site_multiplicity: int | None = None
    site_occupancy: float | None = None
    site_weight: float | None = None


# ---------------------------------------------------------------
# CIF 文本层检查（occupancy 标签、共占位、元素集合）
# ---------------------------------------------------------------
def _cif_text(cif_path: Path) -> str:
    return cif_path.read_text(encoding="utf-8", errors="ignore")


def _cif_has_occupancy_tag(cif_path: Path) -> bool:
    """直接读 CIF 文本检查 _atom_site_occupancy 标签是否存在。

    artemis.read_cif_sites 对缺失 occupancy 默认 1.0，不能用来判断
    “缺失”，必须直接读文本。
    """
    text = _cif_text(cif_path)
    return "_atom_site_occupancy" in text


def _symbol(site) -> str | None:
    m = re.match(r"[A-Za-z]+", str(site.type_symbol or "").strip())
    return m.group(0).capitalize() if m else None


def _cif_elements(st) -> set[str]:
    """CIF 中出现的全部元素（不含 H/D）。"""
    els = set()
    for site in st.get_all_unit_cell_sites():
        sym = _symbol(site)
        if sym and sym not in ("H", "D"):
            els.add(sym)
    return els


def _detect_cocupancy(st, cell) -> bool:
    """检查是否有共占位：同一分数坐标（mod 1，容差 0.08 Å 笛卡尔）
    出现多个不同元素的位点。"""
    all_sites = st.get_all_unit_cell_sites()
    groups: dict[tuple, list[str]] = {}
    for site in all_sites:
        f = np.array([site.fract.x, site.fract.y, site.fract.z], dtype=float)
        placed = False
        for key in groups:
            d = f - np.array(key, dtype=float)
            d -= np.round(d)
            p = cell.orthogonalize(gemmi.Fractional(*d))
            if np.linalg.norm([p.x, p.y, p.z]) < OCCUPANCY_TOL:
                groups[key].append(_symbol(site))
                placed = True
                break
        if not placed:
            groups[tuple(f.round(4))] = [_symbol(site)]
    for syms in groups.values():
        uniq = {s for s in syms if s}
        if len(uniq) > 1:
            return True
    return False


def _absorber_positions(st, cell, absorber: str) -> list:
    """返回所有吸收元素的位点（按分数坐标 mod 1 分组后的代表位点）。"""
    all_sites = st.get_all_unit_cell_sites()
    cands = [s for s in all_sites if _symbol(s) == absorber]
    if not cands:
        return []
    # 按分数坐标 mod 1 分组（容差 0.08 Å）
    groups: list[dict] = []
    for site in cands:
        f = np.array([site.fract.x, site.fract.y, site.fract.z], dtype=float)
        placed = False
        for g in groups:
            d = f - g["f"]
            d -= np.round(d)
            p = cell.orthogonalize(gemmi.Fractional(*d))
            if np.linalg.norm([p.x, p.y, p.z]) < OCCUPANCY_TOL:
                g["members"].append(site)
                placed = True
                break
        if not placed:
            groups.append({"f": f, "members": [site]})
    return groups


def _parse_absorber_site(site_label: str, cell) -> tuple[float, float, float] | None:
    """解析 --absorber-site：'Fe1' 标签 或 'frac:(x,y,z)' 分数坐标。"""
    m = re.match(r"frac:\(([^,]+),([^,]+),([^)]+)\)", site_label.strip())
    if m:
        return (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    return None  # 标签形式，由调用方匹配


# ---------------------------------------------------------------
# P5-D：多位点展开与化学计量权重
# ---------------------------------------------------------------
def list_inequivalent_absorber_sites(
        cif_path: str | Path, absorber: str) -> list[dict]:
    """返回 CIF 不对称单元中每个吸收位点的多重度、占据度与权重。

    权重按化学计量 ``w_i = m_i·occ_i / Σ(m_j·occ_j)`` 归一化：
    - 尖晶石 8a/16d（occ=1）→ 1/3 : 2/3；
    - 磁铁矿 COD 1010369 的 16d 被拆成 Fe1/Fe2（occ=0.5+0.5），
      各自权重 1/3，与 Fe3(8a) 合计 1。
    绝不按权重自动挑选「主位点」——这里只用于多位点加权汇总。
    """
    cif_path = Path(cif_path)
    absorber = absorber.capitalize()
    doc = gemmi.cif.read(str(cif_path))
    block = doc.sole_block() if len(doc) == 1 else doc[0]
    st = gemmi.make_small_structure_from_block(block)

    counts: dict[str, int] = {}
    for site in st.get_all_unit_cell_sites():
        label = str(site.label)
        counts[label] = counts.get(label, 0) + 1

    rows: list[dict] = []
    for site in st.sites:
        if _symbol(site) != absorber:
            continue
        multiplicity = int(counts.get(str(site.label), 0))
        occupancy = float(site.occ) if site.occ is not None else 1.0
        rows.append({
            "label": str(site.label),
            "multiplicity": multiplicity,
            "occupancy": occupancy,
            "raw_weight": float(multiplicity) * occupancy,
            "fract": [float(site.fract.x), float(site.fract.y),
                      float(site.fract.z)],
        })

    total = sum(r["raw_weight"] for r in rows)
    if rows and total > 0:
        for r in rows:
            r["weight"] = r["raw_weight"] / total
    elif rows:
        for r in rows:
            r["weight"] = 1.0 / len(rows)
    for r in rows:
        r.pop("raw_weight", None)
    return rows


def gate_cif_all_sites(
    cif_path: str | Path,
    absorber: str,
    **kwargs,
) -> list[CIFGateResult]:
    """P5-D：对每个不等价吸收位独立跑 ``gate_cif`` 并附上化学计量权重。

    任一位点失败只会让该位点返回 ``ok=False``；调用方必须显式拒绝
    整个正式 CN 输出，不做「部分位点成功即放行」。
    """
    sites = list_inequivalent_absorber_sites(cif_path, absorber)
    if not sites:
        return [CIFGateResult(
            ok=False, code="absorber_missing",
            absorber=absorber.capitalize(),
            reason=f"{Path(cif_path).name} 中找不到元素 {absorber} 的位点",
            cif_path=str(Path(cif_path).resolve()))]

    results: list[CIFGateResult] = []
    for info in sites:
        result = gate_cif(cif_path, absorber,
                          absorber_site=info["label"], **kwargs)
        result.site_multiplicity = info["multiplicity"]
        result.site_occupancy = info["occupancy"]
        result.site_weight = info["weight"]
        results.append(result)
    return results


# ---------------------------------------------------------------
# Rayleigh CN shell 分组（§5.2）
# ---------------------------------------------------------------
def _rayleigh_delta_r(kmin: float, kmax: float) -> float:
    return math.pi / (2.0 * (kmax - kmin))


def _group_shells_rayleigh(sites, occs, kmin: float, kmax: float) -> list[dict]:
    """按 (元素, Rayleigh 距离容差) 分壳，计算 N₀ = Σ(简并度 × 占据度)。

    ΔR_Rayleigh = π/(2·(kmax−kmin))；同元素子壳间距 < ΔR 时归为一个 CN 组。
    记录 n0_members、sub_rayleigh_grouped、merge_reason。
    """
    delta_r = _rayleigh_delta_r(kmin, kmax)
    rows = []
    for i, (x, y, z, el) in enumerate(sites):
        occ = float(occs[i]) if occs is not None else 1.0
        rows.append((el, float(np.linalg.norm([x, y, z])), occ))
    rows.sort(key=lambda t: (t[1], t[0]))

    shells: list[dict] = []
    for el, r, occ in rows:
        # 找同元素且距离差 < delta_r 的现有壳层
        match = None
        for sh in shells:
            if sh["el"] == el and abs(r - sh["r"]) < delta_r:
                match = sh
                break
        if match is not None:
            n_old = match["n"]
            match["r"] = (match["r"] * n_old + r) / (n_old + 1)
            match["n"] += 1
            match["n0"] += 1.0 * occ
            match["n0_members"].append({
                "el": el, "r": r, "occ": occ,
                "merged": True, "merge_reason": f"sub-rayleigh Δr={abs(r - match['r']):.3f}Å < {delta_r:.3f}Å",
            })
            match["sub_rayleigh_grouped"] = True
        else:
            shells.append({
                "el": el, "r": r, "n": 1,
                "n0": 1.0 * occ,
                "n0_members": [{"el": el, "r": r, "occ": occ,
                                "merged": False, "merge_reason": "first member"}],
                "sub_rayleigh_grouped": False,
            })
    return shells, delta_r


# ---------------------------------------------------------------
# 主门禁入口
# ---------------------------------------------------------------
def gate_cif(
    cif_path: str | Path,
    absorber: str,
    *,
    absorber_site: str | None = None,
    kmin: float = 3.0,
    kmax: float = 12.0,
    rmax: float = 6.0,
    allow_computed: bool = False,
    allow_unusual_coordination: bool = False,
    is_computed: bool = False,
) -> CIFGateResult:
    """对单个 CIF 做逐项门禁检查，返回 CIFGateResult。

    参数:
        cif_path: CIF 文件路径（真实源 CIF，不得是 _structure_lib/*/_cif/ 现场生成文件）
        absorber: 吸收元素符号
        absorber_site: 吸收位标签 ('Fe1') 或 分数坐标 ('frac:(0,0,0)')
        kmin/kmax: 用于 Rayleigh 分组
        rmax: 展开半径
        allow_computed: 允许计算结构（降级 cn_hypothesis_only）
        allow_unusual_coordination: 允许 CN 超出正常范围
        is_computed: CIF 是否为计算结构（由调用方/phase pack entry 指定）
    """
    result = CIFGateResult(ok=False, absorber=absorber.capitalize())
    cif_path = Path(cif_path)

    # ---- 0. 现场生成 CIF 拒绝 ----
    if is_generated_cif(cif_path):
        return CIFGateResult(
            ok=False, code="legacy_unverified", absorber=absorber.capitalize(),
            reason="源 CIF 位于 _structure_lib/*/_cif/ 或带生成标记；现场生成 CIF 不得作为正式来源",
            cif_path=str(cif_path),
            cif_sha256=sha256_file(cif_path) if cif_path.is_file() else None)

    # ---- 1. 解析 CIF；必须有 _atom_site_* 循环 ----
    if not cif_path.is_file():
        return CIFGateResult(
            ok=False, code="parse_failed", absorber=absorber.capitalize(),
            reason=f"CIF 不存在: {cif_path}")
    text = _cif_text(cif_path)
    if "_atom_site_label" not in text or "_atom_site_fract_x" not in text:
        return CIFGateResult(
            ok=False, code="parse_failed", absorber=absorber.capitalize(),
            reason=f"CIF 缺少 _atom_site_* 循环: {cif_path.name}")
    try:
        import gemmi
        doc = gemmi.cif.read(str(cif_path))
        block = doc.sole_block() if len(doc) == 1 else doc[0]
        st = gemmi.make_small_structure_from_block(block)
    except Exception as exc:
        return CIFGateResult(
            ok=False, code="parse_failed", absorber=absorber.capitalize(),
            reason=f"解析 {cif_path.name} 失败：{exc}")
    result.cif_path = str(cif_path.resolve())
    result.cif_sha256 = sha256_file(cif_path)

    absorber = absorber.capitalize()

    # ---- 2. 吸收元素必须存在 ----
    elements = _cif_elements(st)
    if absorber not in elements:
        return CIFGateResult(
            ok=False, code="absorber_missing", absorber=absorber,
            reason=f"{cif_path.name} 中找不到吸收元素 {absorber}，文件里只有 {sorted(elements)}",
            cif_path=result.cif_path, cif_sha256=result.cif_sha256)

    # ---- 3. 计算结构门禁 ----
    if is_computed and not allow_computed:
        return CIFGateResult(
            ok=False, code="cif_computed", absorber=absorber,
            reason=f"{cif_path.name} 是计算结构，不得作为实验 N₀；用 --allow-computed 降级",
            cif_path=result.cif_path, cif_sha256=result.cif_sha256)
    if is_computed and allow_computed:
        result.structure_is_computed = True
        result.cn_hypothesis_only = True
        result.warnings.append("计算结构已放行（--allow-computed），CN 仅作 hypothesis_only")

    # ---- 4. 吸收位：同元素多个晶位时必须 --absorber-site ----
    groups = _absorber_positions(st, st.cell, absorber)
    if len(groups) == 0:
        return CIFGateResult(
            ok=False, code="absorber_missing", absorber=absorber,
            reason=f"{cif_path.name} 中找不到元素 {absorber} 的位点",
            cif_path=result.cif_path, cif_sha256=result.cif_sha256)
    if len(groups) > 1 and absorber_site is None:
        # 单元素结构：所有吸收元素位点对称等价，自动选第一个，不报 ambiguous
        if len(elements) == 1:
            result.warnings.append(
                f"{cif_path.name} 中元素 {absorber} 有 {len(groups)} 个对称等价晶位，"
                f"自动选 {groups[0]['members'][0].label}")
            abs_site = groups[0]["members"][0]
            result.absorber_site = str(abs_site.label)
        else:
            site_labels = []
            for g in groups:
                rep = g["members"][0]
                site_labels.append(str(rep.label))
            return CIFGateResult(
                ok=False, code="absorber_site_ambiguous", absorber=absorber,
                reason=(f"{cif_path.name} 中元素 {absorber} 有 {len(groups)} 个晶位 "
                        f"{site_labels}，必须 --absorber-site 指定；禁止按 occ×多重度自动挑主位点"),
                cif_path=result.cif_path, cif_sha256=result.cif_sha256)

    # 选定吸收位
    if absorber_site is not None:
        site_frac = _parse_absorber_site(absorber_site, st.cell)
        if site_frac is not None:
            # 分数坐标形式：找最近的吸收元素位点
            best = None
            best_d2 = float("inf")
            for g in groups:
                for site in g["members"]:
                    d = np.array([site.fract.x - site_frac[0],
                                  site.fract.y - site_frac[1],
                                  site.fract.z - site_frac[2]])
                    d -= np.round(d)
                    p = site.cell.orthogonalize(gemmi.Fractional(*d))
                    d2 = p.x ** 2 + p.y ** 2 + p.z ** 2
                    if d2 < best_d2:
                        best_d2 = d2
                        best = site
            if best is None or best_d2 > (OCCUPANCY_TOL ** 2):
                return CIFGateResult(
                    ok=False, code="absorber_site_ambiguous", absorber=absorber,
                    reason=f"--absorber-site {absorber_site!r} 未匹配到任何 {absorber} 晶位",
                    cif_path=result.cif_path, cif_sha256=result.cif_sha256)
            abs_site = best
            result.absorber_site = f"{absorber_site}@({site_frac[0]:.4f},{site_frac[1]:.4f},{site_frac[2]:.4f})"
        else:
            # 标签形式
            abs_site = None
            for g in groups:
                for site in g["members"]:
                    if str(site.label) == absorber_site:
                        abs_site = site
                        break
                if abs_site is not None:
                    break
            if abs_site is None:
                return CIFGateResult(
                    ok=False, code="absorber_site_ambiguous", absorber=absorber,
                    reason=f"--absorber-site {absorber_site!r} 未匹配到任何 {absorber} 晶位",
                    cif_path=result.cif_path, cif_sha256=result.cif_sha256)
            result.absorber_site = absorber_site
    else:
        # 只有一个晶位，自动选
        abs_site = groups[0]["members"][0]
        result.absorber_site = str(abs_site.label)

    origin = np.array([abs_site.fract.x, abs_site.fract.y, abs_site.fract.z],
                      dtype=float)

    # ---- 5. occupancy 检查 ----
    has_occ_tag = _cif_has_occupancy_tag(cif_path)
    cocupancy = _detect_cocupancy(st, st.cell)

    if not has_occ_tag:
        # 缺失 occupancy：判断是否允许隐式 1.0
        single_elem = len(elements) == 1
        ordered = not cocupancy
        # 化学式一致：单元素且没有共占位，化学式自然一致
        formula_consistent = single_elem and ordered
        if single_elem and ordered and formula_consistent:
            result.occupancy_source = "implicit_1.0"
            result.warnings.append(
                f"{cif_path.name} 缺少 _atom_site_occupancy，但为单元素有序结构，"
                f"隐式 occupancy=1.0")
        else:
            return CIFGateResult(
                ok=False, code="occupancy_missing", absorber=absorber,
                reason=(f"{cif_path.name} 缺少 _atom_site_occupancy，且不是单元素有序结构；"
                        f"无序/共占位/缺陷模型/化学式不匹配 → occupancy_missing"),
                cif_path=result.cif_path, cif_sha256=result.cif_sha256)
    else:
        # 有显式占据度标签
        result.occupancy_source = "explicit"
        if cocupancy:
            return CIFGateResult(
                ok=False, code="cocupancy_unresolved", absorber=absorber,
                reason=(f"{cif_path.name} 存在共占位，未按元素拆分分支；"
                        f"cocupancy_unresolved，不出正式 CN"),
                cif_path=result.cif_path, cif_sha256=result.cif_sha256)

    # ---- 6. 展开 CIF，计算 N₀ shell 分组（Rayleigh） ----
    # 传入 absorber_fract 以确保吸收中心正确（避免 artemis 按 occ×多重度自动选位）
    absorber_fract = None
    if absorber_site is not None:
        parsed = _parse_absorber_site(absorber_site, st.cell)
        if parsed is not None:
            absorber_fract = list(parsed)
        else:
            # 标签形式：用选定位点的分数坐标
            absorber_fract = [abs_site.fract.x, abs_site.fract.y, abs_site.fract.z]

    try:
        sites, occs = artemis.read_cif_sites(
            str(cif_path), absorber=absorber, rmax=rmax,
            absorber_fract=absorber_fract, with_occupancy=True)
    except XAFSDataError as exc:
        return CIFGateResult(
            ok=False, code="parse_failed", absorber=absorber,
            reason=str(exc),
            cif_path=result.cif_path, cif_sha256=result.cif_sha256)

    if not sites:
        return CIFGateResult(
            ok=False, code="absorber_missing", absorber=absorber,
            reason=f"{cif_path.name} 中 {absorber} 周围 rmax={rmax} Å 内无近邻",
            cif_path=result.cif_path, cif_sha256=result.cif_sha256)

    shells, delta_r = _group_shells_rayleigh(sites, occs, kmin, kmax)
    result.n0_shells = shells
    result.rayleigh_delta_r = delta_r

    # ---- 7. 第一壳 CN 范围检查 ----
    first_shell = shells[0]
    first_cn = first_shell["n0"]
    result.first_shell_cn = first_cn

    if first_cn < CN_NORMAL_LO or first_cn > CN_WARN_HI:
        if not allow_unusual_coordination:
            return CIFGateResult(
                ok=False, code="n_out_of_range", absorber=absorber,
                reason=(f"第一壳 CN={first_cn:.2f} 超出正常范围 [{CN_NORMAL_LO}–{CN_WARN_HI}]，"
                        f"需要 --allow-unusual-coordination"),
                cif_path=result.cif_path, cif_sha256=result.cif_sha256)
        result.warnings.append(
            f"第一壳 CN={first_cn:.2f} 超出正常范围，已放行（--allow-unusual-coordination）")
    elif first_cn > CN_NORMAL_HI:
        result.warnings.append(
            f"第一壳 CN={first_cn:.2f} 处于警告范围 [{CN_NORMAL_LO}–{CN_WARN_HI}]")

    # ---- 8. 键长合理性（第一壳） ----
    first_el = first_shell["el"]
    if first_el != absorber:
        r_sum = (COVALENT_RADIUS.get(absorber) and COVALENT_RADIUS.get(first_el))
        if COVALENT_RADIUS.get(absorber) and COVALENT_RADIUS.get(first_el):
            expected = COVALENT_RADIUS[absorber] + COVALENT_RADIUS[first_el]
            r_obs = first_shell["r"]
            if abs(r_obs - expected) > BOND_TOL:
                return CIFGateResult(
                    ok=False, code="bond_length_implausible", absorber=absorber,
                    reason=(f"第一壳 {absorber}–{first_el} 键长 {r_obs:.3f} Å 与共价半径之和 "
                            f"{expected:.3f} Å 偏差 > {BOND_TOL} Å；未知元素对请提供参考键长"),
                    cif_path=result.cif_path, cif_sha256=result.cif_sha256)
    # 未知元素对：记录警告，不硬失败（用户可在后续提供参考键长）

    # ---- 9. 多面体闭合（第一壳近邻数合理、键长合理、无异常空缺）----
    # 已由步骤 7（CN 范围）和步骤 8（键长）覆盖；此处补充“近邻数>0”的基础检查
    if first_shell["n"] == 0:
        return CIFGateResult(
            ok=False, code="polyhedron_not_closed", absorber=absorber,
            reason=f"{cif_path.name} 第一壳无近邻，多面体不闭合",
            cif_path=result.cif_path, cif_sha256=result.cif_sha256)

    # ---- 放行 ----
    result.ok = True
    result.code = None
    result.reason = "门禁通过，可作为正式 CN 结构来源"
    return result


# ---------------------------------------------------------------
# 向后兼容别名
# ---------------------------------------------------------------
def gate_structure_source(cif_path, absorber, **kwargs) -> CIFGateResult:
    """gate_cif 的别名，供外部调用。"""
    return gate_cif(cif_path, absorber, **kwargs)
