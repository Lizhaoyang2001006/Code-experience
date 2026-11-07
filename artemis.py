# -*- coding: utf-8 -*-
"""
artemis.py — 阶段 5~7：由 FEFF 路径到结构参数
================================================================
对应 Artemis 的功能范围。函数名与 EXAFS 方程里的每一项对应：

    FEFF 输入       build_feff_input
    CIF 读取        read_cif_sites（gemmi 对称展开 → 以吸收原子为原点的近邻坐标）
    路径计算        run_feff6l → load_feff_path
    设参数          set_path_parameters （ΔR、σ²、S₀²、ΔE₀；N 直接改 path.degen）
    理论谱          exafs_model_k （标准 EXAFS 方程求和）
    R 空间拟合      fit_r_space （复数 χ(R) 的实部与虚部同时最小二乘，可选拟合 N）
    统计量          independent_points / r_factor

依赖方向：core ← io ← athena ← **artemis**（拟合需要 athena 的 FT）。
本模块不 import larch_wt。
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

import athena
from xafs_core import (FitResult, PathFit, XAFSConfigError, XAFSDataError)

# ================================================================
# 1. FEFF 输入生成
# ================================================================
#
# 格式必须以 Artemis 导出的 feff.inp 为基准，以下是硬约束：
#   - 段名是 POTENTIALS（复数），三列表格：ipot  Z  element
#   - ATOMS 行 6 字段：x  y  z  ipot  tag  distance
#   - CONTROL / PRINT 各 4 个参数（不是 6 个）
#   - 无 EDGE 关键字，吸收边由 HOLE 卡指定；吸收原子 ipot=0
#     HOLE 1=K, 2=L1, 3=L2, 4=L3（FEFF8 约定）
# 曾经写成 POTENTIAL（单数）+ 两行 "26 Fe 2.0"，FEFF 把 26 当成势编号，
# 报 "Unique potential 26 not allowed"，连续 6 种写法全被拒。

_ELEMENTS = (
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)
ATOMIC_NUMBER = {symbol: z for z, symbol in enumerate(_ELEMENTS, start=1)}
HOLE_INDEX = {"K": 1, "L1": 2, "L2": 3, "L3": 4}

# 共价半径 (Å)，Cordero et al., Dalton Trans. 2008。只收 XAFS 常见元素；
# 缺失的键返回 None，调用方必须显式降级而不是猜。
# 用途：L1 第一壳的**峰值选择**（不是最终键长）——预期 M–L 键长
# = r(absorber)+r(ligand)，再减相位偏移，用来在多个显著峰中选最像
# 第一壳的那个。实测：Fe-O 1.98（hematite 1.35 峰而非 Fe-Fe 2.55）、
# Pt-Cl 2.38（1.90 峰而非低 R 白线伪峰 1.44）。
COVALENT_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Se": 1.20, "Br": 1.20,
    "I": 1.39,
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32,
    "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Tc": 1.47, "Ru": 1.46,
    "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "La": 2.07, "Ce": 2.04, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51,
    "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
    "Tl": 1.45, "Pb": 1.46, "Bi": 1.48, "U": 1.96,
}

# BVS 参考键长 R₀ (Å)，Brown & Altermatt, Acta Cryst. B41 (1985) 244。
# 键写法同时接受 "Fe-O" 与 ("Fe", "O")；b=0.37 Å 是同一文献的通用值。
BVS_R0 = {
    "Fe-O": 1.759, "Fe-N": 1.90, "Fe-S": 2.20, "Fe-Cl": 2.29,
    "Co-O": 1.692, "Co-N": 1.85, "Co-S": 2.16,
    "Ni-O": 1.654, "Ni-N": 1.83, "Ni-S": 2.13,
    "Cu-O": 1.679, "Cu-N": 1.85, "Cu-S": 2.15, "Cu-Cl": 2.24,
    "Mn-O": 1.790, "Mn-N": 1.95, "Cr-O": 1.720, "Zn-O": 1.704,
    "Mo-O": 1.860, "Mo-S": 2.32, "Ru-O": 1.834, "Rh-O": 1.820,
    "Pd-O": 1.792, "Ag-O": 1.842, "Pt-Cl": 2.190, "Pt-O": 1.820,
    "Ir-O": 1.820, "Au-Cl": 2.220, "Ce-O": 2.340, "U-O": 2.300,
}


def expected_bond_length(absorber: str, ligand: str) -> float | None:
    """共价半径和作为 M–L 键长初估；任一元素缺半径返回 None。"""
    ra = COVALENT_RADII.get(absorber)
    rl = COVALENT_RADII.get(ligand)
    if ra is None or rl is None:
        return None
    return ra + rl


def build_feff_input(sites, absorber: str | None = None, rmax: float = 6.0,
                     title: str = "auto", edge: str = "K") -> str:
    """
    由中心原子 + 近邻坐标生成 feff.inp 文本。

    参数:
        sites    : [(x, y, z, 元素符号), ...] 笛卡尔坐标（Å），中心固定在原点
        absorber : 中心原子元素符号
        rmax     : 团簇半径（Å），FEFF 据此扩展周期
        edge     : 吸收边 "K" / "L1" / "L2" / "L3"（决定 HOLE 卡）

    返回: feff.inp 的完整文本
    """
    if not absorber:
        raise XAFSConfigError("FEFF", "未指定吸收元素（absorber）——生成 FEFF 输入必需")
    if absorber not in ATOMIC_NUMBER:
        raise XAFSConfigError("FEFF", f"未知元素 {absorber}，请补进 ATOMIC_NUMBER")
    if edge not in HOLE_INDEX:
        raise XAFSConfigError(
            "FEFF", f"未知吸收边 {edge!r}；可选 {sorted(HOLE_INDEX)}")
    if not sites:
        raise XAFSConfigError("FEFF", "sites 为空")

    scatterers = []
    for site in sites:
        symbol = site[3] if len(site) > 3 else absorber
        if symbol not in scatterers:
            scatterers.append(symbol)

    potential_lines = [f"     0    {ATOMIC_NUMBER[absorber]:3d}   {absorber}"]
    for index, symbol in enumerate(scatterers, start=1):
        if symbol not in ATOMIC_NUMBER:
            raise XAFSConfigError("FEFF", f"未知元素 {symbol}")
        potential_lines.append(f"    {index:3d}    {ATOMIC_NUMBER[symbol]:3d}   {symbol}")

    atom_lines = [f"   0.00000    0.00000    0.00000  0    {absorber:<4s}   0.00000"]
    for site in sorted(sites, key=lambda s: sum(v * v for v in s[:3])):
        x, y, z = site[0], site[1], site[2]
        distance = float(np.sqrt(x * x + y * y + z * z))
        if distance < 0.1 or distance > rmax:
            continue
        symbol = site[3] if len(site) > 3 else absorber
        ipot = scatterers.index(symbol) + 1
        atom_lines.append(
            f"  {x:9.5f}  {y:9.5f}  {z:9.5f}  {ipot}    {symbol:<4s}  {distance:9.5f}")

    return ("* feff.inp generated by artemis.py\n"
            f"TITLE {title}\n\n"
            f"HOLE {HOLE_INDEX[edge]}   1.0\n"
            "*         mphase,mpath,mfeff,mchi\n"
            "CONTROL   1      1     1     1\n"
            "PRINT     1      0     0     0\n\n"
            f"RMAX  {rmax:.1f}\n"
            "CORRECTIONS  0.0   0.0\n\n"
            "POTENTIALS\n"
            "*    ipot   Z  element\n"
            + "\n".join(potential_lines) + "\n\n"
            "ATOMS\n"
            "*   x          y          z      ipot  tag      distance\n"
            + "\n".join(atom_lines) + "\n")


def bcc_lattice_sites(a: float = 2.866, rmax: float = 6.0, symbol: str = "Fe"):
    """体心立方金属（Fe 箔等）的近邻坐标。a=2.866 Å 时第一壳 2.482 Å ×8。"""
    basis = [(0, 0, 0), (0.5, 0.5, 0.5)]
    n = int(rmax / a) + 2
    sites = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                for b in basis:
                    p = ((i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a)
                    d = float(np.sqrt(sum(v * v for v in p)))
                    if 0.1 < d <= rmax:
                        sites.append((*p, symbol))
    return sites


def octahedron_sites(bond_length: float = 2.00, symbol: str = "O"):
    """单个 MO6 正八面体，用于无晶体结构时生成第一壳层路径。"""
    sites = []
    for axis in range(3):
        for sign in (+1, -1):
            p = [0.0, 0.0, 0.0]
            p[axis] = sign * bond_length
            sites.append((p[0], p[1], p[2], symbol))
    return sites


def read_cif_sites(cif_path, absorber: str | None = None, rmax: float = 6.0,
                   ignore=("H", "D"), absorber_fract=None,
                   with_occupancy: bool = False,
                   return_occupancy: bool = False,
                   merge_cocupancy: bool = True) -> list:
    """
    读 CIF 晶体结构，展开对称等价位置，返回以吸收原子为原点的近邻坐标。

    数学/流程:
        不对称单元 → 空间群操作展开 → 全晶胞原子
        → 分数坐标 × 晶胞矩阵 → 笛卡尔坐标（Å）
        → 吸收原子平移到原点，沿 a/b/c 周期平移补全 rmax 内的像

    参数:
        cif_path : CIF 文件路径
        absorber : 吸收原子元素符号（取第一个匹配位点作原点）
        rmax     : 团簇半径（Å），传给 build_feff_input 时保持一致
        ignore   : 忽略的元素。默认跳过 H/D——氢的背散射振幅可忽略，
                   且粉晶 CIF 里 H 的坐标常不可靠
        with_occupancy : True 时返回 (sites, occs)，occs 与 sites 平行
        return_occupancy : True 时返回 5 元组 (x, y, z, symbol, occ)；
                   默认 False 保持历史 4 元组，100% 向后兼容

    部分占据/共占据（B7-2）:
        - 读取 gemmi 的 ``site.occ`` 并钳制到 [0, 1]；
        - 展开后距离 < 0.08 Å 的散射体视为同一晶位共占据，合并为
          一个位点：元素取占据度较大者，有效占据度取
          ``min(1, Σocc)``——防止 FEFF 因两个原子坐标重合崩溃；
        - 吸收中心有多个同元素位点且未给 absorber_fract 时，按
          ``occ × 多重度`` 选择主位点（多重度 = 展开后的等价位点数）。

    限制: 默认 4 元组模式只保证几何正确，N 仍按满占据解读；
        占据敏感的定量必须用 with_occupancy/return_occupancy 通道。
        ``merge_cocupancy=False`` 时保留同一晶位的不同元素分支（B11-1
        异核共占位），不把 Ti/Ru 合并成单一元素近似。
    """
    import re

    import gemmi

    if not absorber:
        raise XAFSConfigError("CIF", "未指定吸收元素（absorber）——选吸收位必需")
    absorber = absorber.capitalize()
    ignored = {s.capitalize() for s in ignore}

    try:
        doc = gemmi.cif.read(str(cif_path))
        block = doc.sole_block() if len(doc) == 1 else doc[0]
        st = gemmi.make_small_structure_from_block(block)
        cell = st.cell
        all_sites = st.get_all_unit_cell_sites()
    except (RuntimeError, ValueError, IndexError) as exc:
        raise XAFSDataError(
            "CIF", f"解析 {Path(cif_path).name} 失败：{exc}。"
                   f"检查晶胞参数 / 空间群 / _atom_site_* 循环是否完整") from exc

    def _symbol(site):
        m = re.match(r"[A-Za-z]+", site.type_symbol.strip())
        return m.group(0).capitalize() if m else None

    def _occ(site):
        try:
            value = float(getattr(site, "occ", 1.0))
        except (TypeError, ValueError):
            return 1.0
        if not np.isfinite(value):
            return 1.0
        return float(min(1.0, max(0.0, value)))

    def _cart_site(site):
        p = cell.orthogonalize(site.fract)
        return np.array([p.x, p.y, p.z])

    # ---- 吸收中心：无显式分数坐标时按 occ × 多重度选主位点 ----
    abs_site = None
    if absorber_fract is None:
        cands = [s for s in all_sites if _symbol(s) == absorber]
        if cands:
            # 以分数坐标 mod 1 分组（容差换算成笛卡尔 0.08 Å）
            groups = []
            for s in cands:
                f = np.array([s.fract.x, s.fract.y, s.fract.z], dtype=float)
                placed = False
                for g in groups:
                    d = f - g["f"]
                    d -= np.round(d)
                    p = cell.orthogonalize(gemmi.Fractional(*d))
                    if np.linalg.norm([p.x, p.y, p.z]) < 0.08:
                        g["members"].append(s)
                        placed = True
                        break
                if not placed:
                    groups.append({"f": f, "members": [s]})
            if groups:
                def _group_weight(g):
                    mult = len(g["members"])
                    occ_sum = sum(_occ(s) for s in g["members"]) / mult
                    return occ_sum * mult
                abs_site = max(groups, key=_group_weight)["members"][0]
    else:
        fx, fy, fz = (float(v) for v in absorber_fract)
        cands = [s for s in all_sites if _symbol(s) == absorber]

        def _d2(site):
            d = np.array([site.fract.x - fx, site.fract.y - fy,
                          site.fract.z - fz])
            d -= np.round(d)          # 最近像约定
            p = cell.orthogonalize(gemmi.Fractional(*d))
            return p.x ** 2 + p.y ** 2 + p.z ** 2

        abs_site = min(cands, key=_d2) if cands else None
    if abs_site is None:
        found = sorted({sym for sym in map(_symbol, all_sites) if sym})
        raise XAFSDataError(
            "CIF", f"{Path(cif_path).name} 中找不到元素 {absorber}，"
                   f"文件里只有 {found}")

    origin = _cart_site(abs_site)

    def _lattice_vec(fx, fy, fz):
        p = cell.orthogonalize(gemmi.Fractional(fx, fy, fz))
        return np.array([p.x, p.y, p.z])

    a_vec = _lattice_vec(1, 0, 0)
    b_vec = _lattice_vec(0, 1, 0)
    c_vec = _lattice_vec(0, 0, 1)

    # 周期像范围：rmax / 轴长 + 1 圈保险
    ranges = []
    for vec in (a_vec, b_vec, c_vec):
        n = int(np.ceil(rmax / max(np.linalg.norm(vec), 1e-6))) + 1
        ranges.append(range(-n, n + 1))

    # ---- 展开 + 共占位合并（先按距离筛选，再合并，避免 O(N²) 爆炸）----
    # 按占据度降序排列，保证同一晶位保留占据度较大的元素符号。
    raw = []
    for site in all_sites:
        sym = _symbol(site)
        if sym is None or sym in ignored:
            continue
        base = _cart_site(site) - origin
        occ = _occ(site)
        for i, j, k in ((i, j, k) for i in ranges[0] for j in ranges[1]
                        for k in ranges[2]):
            pos = base + i * a_vec + j * b_vec + k * c_vec
            d = float(np.linalg.norm(pos))
            if 0.1 < d <= rmax:
                raw.append((pos, sym, occ, d))
    raw.sort(key=lambda item: (-item[2], item[1]))

    # 空间哈希：网格边长 0.08 Å，只需比较相邻 3³ 个桶。
    cell_size = 0.08
    buckets = {}
    merged = []      # [{"pos": ndarray, "sym": str, "occ": float, "d": float}]
    for pos, sym, occ, d in raw:
        key = tuple(int(v) for v in np.floor(pos / cell_size))
        match = None
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for m in buckets.get((key[0] + dx, key[1] + dy,
                                          key[2] + dz), ()):
                        if float(np.linalg.norm(pos - m["pos"])) < 0.08:
                            if merge_cocupancy or m["sym"] == sym:
                                match = m
                                break
                    if match is not None:
                        break
                if match is not None:
                    break
            if match is not None:
                break
        if match is None:
            item = {"pos": pos, "sym": sym, "occ": occ, "d": d}
            merged.append(item)
            buckets.setdefault(key, []).append(item)
        else:
            if merge_cocupancy:
                match["occ"] = float(min(1.0, match["occ"] + occ))
            else:
                # 分支模式：同一元素同一位置的周期像不重复计占据。
                match["occ"] = float(max(match["occ"], occ))
    merged.sort(key=lambda m: (m["d"], m["sym"]))

    sites = []
    occs = []
    for m in merged:
        sites.append((float(m["pos"][0]), float(m["pos"][1]),
                      float(m["pos"][2]), m["sym"]))
        occs.append(float(m["occ"]))
    if return_occupancy:
        return [(x, y, z, sym, occ)
                for (x, y, z, sym), occ in zip(sites, occs)]
    return (sites, occs) if with_occupancy else sites


def cif_hetero_cocupancy_branches(cif_path, absorber: str, rmax: float = 6.0,
                                  delta_z: int = 3,
                                  occ_min: float = 0.10) -> dict:
    """B11-1 异核共占位的双分支（Ti/Ru、Fe/Ni 等）。

    普通 ``read_cif_sites`` 在 <0.08 Å 的同一晶位上只保留占据度较大的
    元素，等效于“单元素近似”。真实催化剂的共占位位点必须拆成两支：
    每个元素有自己的散射振幅和占据度，但共享同一个 R 与 σ²（在拟合层
    用 ``shared_bonds`` 绑定）。

    返回::

        {"single_element_approx": bool,
         "branches": [{"distance", "element", "occ", "degen",
                       "shell", "n_geom"}, ...],
         "shells": [{"distance", "n_geom", "branches": [...]}, ...]}

    简并度约定：同一晶位的几何多重度 N_geom 按占据度分给各元素，
    ``Σ degen = N_geom``。只保留 |Z_a−Z_b|≥delta_z 且 occ≥occ_min 的
    异核分支；同元素共占位按单元素处理。
    """
    sites5 = read_cif_sites(cif_path, absorber=absorber, rmax=rmax,
                            return_occupancy=True, merge_cocupancy=False)
    if not sites5:
        return {"single_element_approx": True, "branches": [], "shells": []}

    # 先按空间位置聚类：同一晶位 = 距离相同且坐标差 <0.08 Å。
    sites = []
    for (x, y, z, sym, occ) in sites5:
        pos = np.array([x, y, z], dtype=float)
        d = float(np.linalg.norm(pos))
        if d <= 0.1 or occ < float(occ_min):
            continue
        sites.append({"pos": pos, "element": str(sym), "occ": float(occ),
                      "distance": d})
    if not sites:
        return {"single_element_approx": True, "branches": [], "shells": []}
    sites.sort(key=lambda s: s["distance"])

    shells = []
    for s in sites:
        placed = False
        for sh in shells:
            # 粉末 EXAFS 的“壳层”按 Reff 归并，与方向无关；同一距离
            # 上的不同晶位（不同元素分支）必须留在同一壳层里，
            # 否则简并度会被拆成 N 个 1，BVS/占位统计全错。
            if abs(sh["distance"] - s["distance"]) < 0.12:
                sh["members"].append(s)
                placed = True
                break
        if not placed:
            shells.append({"distance": s["distance"], "ref_pos": s["pos"],
                           "members": [s]})

    branches = []
    out_shells = []
    any_hetero = False
    for si, sh in enumerate(shells):
        # 几何简并度 = **不同晶位**数，而不是展开后的原子条目数。
        # 同一晶位列了 Ti(0.7)+Ru(0.3) 时，对称展开会得到 8 个 Ti 镜像
        # 和 8 个 Ru 镜像；按条目计 n_geom=16 会把占据度重复乘一遍。
        unique_pos: list[np.ndarray] = []
        for m in sh["members"]:
            if not any(np.linalg.norm(m["pos"] - u) < 0.08
                       for u in unique_pos):
                unique_pos.append(m["pos"])
        n_geom = len(unique_pos)
        by_el: dict[str, dict] = {}
        for m in sh["members"]:
            el = m["element"]
            by_el.setdefault(el, {"occ_sum": 0.0, "count": 0})
            by_el[el]["occ_sum"] += m["occ"]
            by_el[el]["count"] += 1
        hetero = False
        if len(by_el) >= 2:
            zs = sorted(ATOMIC_NUMBER.get(el, 0) for el in by_el)
            hetero = (zs[-1] - zs[0]) >= int(delta_z)
        any_hetero = any_hetero or hetero
        shell_branches = []
        for el, agg in sorted(by_el.items()):
            occ_mean = agg["occ_sum"] / max(agg["count"], 1)
            branch = {
                "distance": float(sh["distance"]),
                "element": el,
                "occ": float(occ_mean),
                "degen": float(n_geom * occ_mean),
                "shell": int(si),
                "n_geom": int(n_geom),
            }
            shell_branches.append(branch)
            branches.append(branch)
        out_shells.append({
            "distance": float(sh["distance"]),
            "n_geom": int(n_geom),
            "single_element_approx": bool(not hetero),
            "branches": shell_branches,
        })

    return {"single_element_approx": bool(not any_hetero),
            "branches": branches, "shells": out_shells}


def path_total_phase(datfile, k):
    """B11-1 返回一条 FEFF 路径的总相移 φ_tot(k)=2δ_c+φ_s。

    优先使用 FEFF 文件自带的 ``pha`` 表（FeffDatFile.pha），按请求的
    k 网格线性插值；老版本/裁剪文件没有 ``pha`` 时，退化为对 ``chi``
    做 Hilbert 解析信号提取相位，并显式记录 ``source`` 字段于返回值
    的附加属性（见 ``path_total_phase.last_source``）。

    注意：相位只定义到一个 2π 常数；FT 只关心 e^{−iφ}，常数相位不改变
    峰位，因此不需要绝对相位零点。
    """
    from larch.xafs import feffpath

    k = np.asarray(k, dtype=float)
    p = feffpath(str(datfile))
    p.create_path_params()
    fd = getattr(p, "_feffdat", None)
    if fd is not None and hasattr(fd, "pha") and hasattr(fd, "k"):
        phi = np.interp(k, np.asarray(fd.k), np.asarray(fd.pha))
        path_total_phase.last_source = "feff_pha"
        return phi
    p._calc_chi()
    from scipy.signal import hilbert
    chi = np.asarray(p.chi, dtype=float)
    kp = np.asarray(p.k, dtype=float)
    z = hilbert(chi)
    arg = np.unwrap(np.angle(z))
    phi = np.interp(k, kp, arg - 2.0 * kp * float(p.reff))
    path_total_phase.last_source = "hilbert_fallback"
    return phi


path_total_phase.last_source = None


def cif_shell_occupancy_scale(cif_path, absorber: str,
                              rmax: float = 6.0) -> list:
    """返回 CIF 展开后每个近邻的 ``(distance, element, occ)``。

    供 ``xafs_fit._ss_rows`` 按距离/元素匹配壳层平均占据度，把 FEFF
    几何简并度缩放为 ``N_cryst = N_geom·mean(occ)``。距离是几何
    距离（与 FEFF 的 Reff 同一口径，不含相移）。
    """
    sites5 = read_cif_sites(cif_path, absorber=absorber, rmax=rmax,
                            return_occupancy=True)
    return [(float(np.linalg.norm([x, y, z])), str(sym), float(occ))
            for (x, y, z, sym, occ) in sites5]


# ================================================================
# 2. FEFF 计算与路径读取
# ================================================================

def run_feff6l(feffinp, workdir) -> list:
    """
    调用 feff6l.exe，返回生成的 feffNNNN.dat 列表。

    feffinp 可以是文本内容，也可以是 .inp 文件路径。
    """
    from larch.xafs import feff6l

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    target = workdir / "feff.inp"
    text = feffinp if "\n" in str(feffinp) else Path(feffinp).read_text()
    target.write_text(text, encoding="utf-8")

    feff6l(feffinp=str(target.resolve()), folder=str(workdir.resolve()), verbose=False)
    files = sorted(workdir.glob("feff0*.dat"))
    if not files:
        raise XAFSDataError("FEFF", f"{workdir} 未生成任何 feffNNNN.dat，"
                                     f"检查 feff.inp 的 POTENTIALS / ATOMS 段")
    return files


def load_feff_path(datfile):
    """
    读入一条散射路径。

    返回 FeffPathGroup，其 params 需要用 set_path_parameters 显式写入。
    """
    from larch.xafs import feffpath

    path = feffpath(str(datfile))
    path.create_path_params()
    return path


def describe_paths(datfiles, max_show: int = 6) -> list:
    """列出路径的 (文件, Reff, 简并度, 散射腿数)，用于挑选壳层。"""
    from larch.xafs import feffpath

    summary = []
    for f in datfiles[:max_show]:
        p = feffpath(str(f))
        summary.append({"file": Path(f).name,
                        "reff": round(float(p.reff), 4),
                        "degeneracy": float(p.degen),
                        "legs": int(p.nleg)})
    return summary


# ================================================================
# 3. 参数设置与理论谱
# ================================================================

# C₃ 的物理量级：1e-4–1e-3 Å³ 是常见非谐；±0.01 Å³ 以上属于
# 相位被 C₃ 吞噬的失配补偿，不能再解释为真实分布偏斜。
C3_LIMIT = 0.01

# C₄ 的物理量级：非高斯展宽的四阶累积量。Larch 参数名 "fourth"，
# 振幅项 exp(+2/3·C₄·k⁴)。|C₄| ≥ 2e-4 Å⁴ 时 k⁴ 项在高 k 端会
# 指数爆炸，只能作为失配补偿的边界信号，不能当真实分布峰度。
C4_LIMIT = 2e-4


def set_path_parameters(path, s02=None, delta_r=None, sigma2=None,
                        delta_e0=None, third=None, fourth=None):
    """
    把拟合参数写进 path.params 并重算 χ(k)。

    **关键（踩过的坑）**:
        larch 的 __path_params 从 self.params[parname] 取值，
        而 params 一旦创建就不再跟随 self.deltar / self.sigma2 属性。
        所以必须走 pathpar_name → params[...].value → _calc_chi()，
        直接改属性或用 calc_chi_from_params() 都不会生效（实测参数逐点不变）。

    third 是三阶累积量 C₃：近邻距离分布偏斜（非高斯，如纳米颗粒表面
    张力导致的 Fe-Fe 分布不对称）时的修正项，进入相位 4/3·C₃·k³。
    fourth 是四阶累积量 C₄：振幅项 exp(+2/3·C₄·k⁴)，描述非高斯展宽。
    """
    mapping = {"s02": s02, "deltar": delta_r, "sigma2": sigma2,
               "e0": delta_e0, "third": third, "fourth": fourth}
    for name, value in mapping.items():
        if value is None:
            continue
        key = path.pathpar_name(name)
        if key in path.params:
            path.params[key].value = float(value)
    path._calc_chi()


def exafs_model_k(paths, k) -> np.ndarray:
    """
    多条路径叠加得到理论 χ(k)。

    数学（标准 EXAFS 方程）:
        χ(k) = Σ_j  S₀²·N_j / (k·R_j²) · F_j(k)
                    · exp(−2R_j/λ(k)) · exp(−2σ_j²k²)
                    · sin(2kR_j + φ_j(k))

    其中 F_j、φ_j、λ 由 FEFF 算出，R_j = Reff_j + ΔR_j。
    """
    k = np.asarray(k)
    total = np.zeros_like(k)
    for path in paths:
        total += np.interp(k, np.asarray(path.k), np.asarray(path.chi))
    return total


# ================================================================
# 4. 统计量
# ================================================================

def independent_points(kmin: float, kmax: float, rmin: float, rmax: float) -> float:
    """
    独立数据点数（Nyquist 判据）。

    数学:
        N_idp = 2·Δk·ΔR / π

    注意：**不加 2**。部分教材写成 2ΔkΔR/π + 2，属保守修正；
    larch / Artemis / IFeffit 均不加，加了会高估自由度、让过拟合更容易通过。
    """
    return float(2.0 * (kmax - kmin) * (rmax - rmin) / np.pi)


def _bond_key(absorber: str, ligand: str):
    return f"{str(absorber).capitalize()}-{str(ligand).capitalize()}"


def _path_ligand(path, absorber: str | None = None):
    """从 PathFit / dict / 标签中提取散射元素，取不到返回 None。"""
    import re

    if isinstance(path, dict):
        for key in ("ligand", "element", "scatterer", "symbol"):
            if path.get(key):
                return str(path[key]).capitalize()
        label = str(path.get("label", ""))
    else:
        label = str(getattr(path, "label", ""))
    symbols = [s.capitalize() for s in re.findall(r"[A-Z][a-z]?", label)]
    symbols = [s for s in symbols if s in ATOMIC_NUMBER]
    if absorber:
        symbols = [s for s in symbols if s != str(absorber).capitalize()]
    return symbols[0] if symbols else None


def _path_n_r(path):
    if isinstance(path, dict):
        n = path.get("N", path.get("degeneracy", path.get("n")))
        r = path.get("R", path.get("r"))
        if r is None:
            r = float(path.get("reff", 0.0)) + float(path.get("delta_r", 0.0))
    else:
        n = getattr(path, "degeneracy", None)
        r = getattr(path, "r", None)
        if r is None:
            r = float(getattr(path, "reff", 0.0)) + float(
                getattr(path, "delta_r", 0.0))
    return float(n), float(r)


def bond_valence_sum(paths, r0_table=None, b: float = 0.37,
                     target: float | None = None,
                     sigma: float | None = None,
                     absorber: str | None = None) -> float:
    """B11-2 键价和（BVS）:

        BVS = Σ_i N_i · exp[(R₀,i − R_i)/b]

    ``paths`` 每项可为 ``PathFit`` 或 dict（N/degeneracy、R/r、
    ligand/element）。R₀ 从内置 ``BVS_R0`` 取，可用 ``r0_table`` 覆盖
    （键写成 "Fe-O" 或 ("Fe","O")）。缺参考键长时抛
    ``XAFSConfigError``——不静默丢掉一条键。

    Fe–O N=6, R=2.0 → 6·exp((1.759−2.0)/0.37) ≈ 3.13（Fe³⁺ 合理）。
    """
    table = dict(BVS_R0)
    if r0_table:
        for key, val in r0_table.items():
            if isinstance(key, (tuple, list)) and len(key) == 2:
                table[_bond_key(key[0], key[1])] = float(val)
            else:
                table[str(key)] = float(val)
    if absorber is None:
        absorber = "Fe"
    if not np.isfinite(b) or b <= 0:
        raise XAFSConfigError("BVS", f"b={b!r} 必须为正")
    total = 0.0
    for path in paths:
        n, r = _path_n_r(path)
        ligand = _path_ligand(path, absorber)
        if ligand is None:
            raise XAFSConfigError(
                "BVS", f"无法确定路径 {path!r} 的散射元素；请给 "
                       f"ligand/element 字段或可解析的标签")
        key = _bond_key(absorber, ligand)
        if key not in table:
            raise XAFSConfigError(
                "BVS", f"缺少参考键长 {key}；在 r0_table 中显式提供，"
                       f"不猜测未标定的 R₀")
        r0 = float(table[key])
        if not np.isfinite(n) or n <= 0 or not np.isfinite(r) or r <= 0:
            raise XAFSConfigError(
                "BVS", f"路径 N={n!r}, R={r!r} 非物理")
        total += n * np.exp((r0 - r) / float(b))
    return float(total)


def correlated_debye_sigma2(temperature_K: float, theta_d: float, reff: float,
                            mass_absorber: float, mass_scatterer: float,
                            k_debye: float | None = None,
                            sigma2_static: float = 0.0) -> float:
    """B11-2 关联 Debye 模型的键长无序 σ²(T, R)。

    非关联项（Debye 固体）:

        σ²_un = (3ħ²/(2 μ k_B Θ_D)) · [1/4
                + (T/Θ_D)² ∫₀^{Θ_D/T} x/(eˣ−1) dx]

    关联因子:

        C(R) = 1 − sin(2·k_D·R)/(2·k_D·R)
        σ²   = σ²_static + C(R)·σ²_un

    近邻（R 小）C→0，运动强相关、σ² 被压低；远壳 C→1，回到非关联值。
    ``k_debye`` 缺省由最近邻密度估计 ``k_D ≈ (36π/R³)^(1/3)`` 并在
    返回值的 ``k_debye_source`` 中显式记录（见 ``last_details``）。
    """
    from scipy.integrate import quad

    T = float(temperature_K)
    theta = float(theta_d)
    r = float(reff)
    m1 = float(mass_absorber)
    m2 = float(mass_scatterer)
    static = float(sigma2_static)
    if not np.isfinite(T) or T <= 0:
        raise XAFSConfigError("DEBYE", f"温度 T={temperature_K!r} 必须 > 0 K")
    if not np.isfinite(theta) or theta <= 0:
        raise XAFSConfigError("DEBYE", f"Θ_D={theta_d!r} 必须 > 0 K")
    if not np.isfinite(r) or r <= 0:
        raise XAFSConfigError("DEBYE", f"Reff={reff!r} 必须 > 0 Å")
    if not (np.isfinite(m1) and m1 > 0 and np.isfinite(m2) and m2 > 0):
        raise XAFSConfigError(
            "DEBYE", f"质量必须为正: m1={mass_absorber!r}, m2={mass_scatterer!r}")
    if not np.isfinite(static) or static < 0:
        raise XAFSConfigError("DEBYE", f"σ²_static={sigma2_static!r} 不能为负")
    mu = m1 * m2 / (m1 + m2)

    upper = theta / T
    integral, _ = quad(lambda x: x / np.expm1(x) if x > 0 else 1.0,
                       0.0, upper, limit=200)
    bracket = 0.25 + (T / theta) ** 2 * float(integral)
    sigma2_un = (3.0 * HBAR2_OVER_2KB / (mu * theta)) * bracket

    if k_debye is None:
        # 保守默认：k_D=0.35 Å⁻¹ 对应低密度固体（最近邻密度量级
        # ~3×10⁻³ Å⁻³）。关键是它**不随 reff 缩放**——若用
        # k_D∝1/reff，2k_D·R 变成常数，近壳与远壳的关联因子相同，
        # 远壳 σ² 不会按物理预期增大。需要精确值时显式传 k_debye。
        kd = 0.35
        source = "conservative_default_0.35A^-1"
    else:
        kd = float(k_debye)
        source = "user"
        if not np.isfinite(kd) or kd <= 0:
            raise XAFSConfigError("DEBYE", f"k_D={k_debye!r} 必须 > 0 Å⁻¹")
    x = 2.0 * kd * r
    corr = 1.0 - (np.sin(x) / x if abs(x) > 1e-12 else 1.0)
    sigma2 = static + corr * sigma2_un
    correlated_debye_sigma2.last_details = {
        "reduced_mass_amu": float(mu),
        "sigma2_un_correlated": float(sigma2_un),
        "correlation_factor": float(corr),
        "k_debye": float(kd),
        "k_debye_source": source,
        "integral": float(integral),
    }
    return float(sigma2)


correlated_debye_sigma2.last_details = None


def sac_cluster_sites(bond_m_ligand: float = 2.00, ligand: str = "N",
                      n_coord: int = 4, axial: str | None = None,
                      axial_bond: float = 2.10, carbon_shell: bool = True,
                      ring_radius: float = 2.90) -> list:
    """B12-2 单原子催化剂（SAC）M–N₄ 平面团簇生成器。

    返回 ``[(x, y, z, element), ...]``，中心吸收原子固定在原点（不包含
    在列表里），可直接送 ``build_feff_input``。

    ``n_coord=4`` 生成正方形 M–N₄；``axial`` 给 "O"/"Cl" 时沿 ±z 加
    轴向配位；``carbon_shell=True`` 时在 R≈ring_radius 的平面内加 8 个
    C（模拟石墨烯/碳氮载体第一壳）。
    """
    d = float(bond_m_ligand)
    n = int(n_coord)
    if n not in (2, 4, 6):
        raise XAFSConfigError(
            "SAC", f"n_coord={n_coord} 只支持 2/4/6（线性/平面/八面体）")
    if not np.isfinite(d) or d <= 0:
        raise XAFSConfigError("SAC", f"bond_m_ligand={bond_m_ligand!r} 必须 > 0")
    lig = str(ligand).capitalize()
    if lig not in ATOMIC_NUMBER:
        raise XAFSConfigError("SAC", f"未知配体元素 {ligand!r}")

    sites = []
    if n == 2:
        sites += [(d, 0.0, 0.0, lig), (-d, 0.0, 0.0, lig)]
    else:
        half = d / np.sqrt(2.0)
        for sx in (-1.0, 1.0):
            for sy in (-1.0, 1.0):
                sites.append((sx * half, sy * half, 0.0, lig))
        if n == 6:
            for sz in (-1.0, 1.0):
                sites.append((0.0, 0.0, sz * d, lig))

    if axial:
        ax = str(axial).capitalize()
        if ax not in ATOMIC_NUMBER:
            raise XAFSConfigError("SAC", f"未知轴向配体元素 {axial!r}")
        if not np.isfinite(axial_bond) or axial_bond <= 0:
            raise XAFSConfigError("SAC", "axial_bond 必须 > 0 Å")
        sites.append((0.0, 0.0, float(axial_bond), ax))
        sites.append((0.0, 0.0, -float(axial_bond), ax))

    if carbon_shell:
        if not np.isfinite(ring_radius) or ring_radius <= 0:
            raise XAFSConfigError("SAC", "ring_radius 必须 > 0 Å")
        for k in range(8):
            ang = 2.0 * np.pi * k / 8.0 + np.pi / 8.0
            sites.append((float(ring_radius * np.cos(ang)),
                          float(ring_radius * np.sin(ang)), 0.0, "C"))
    return sites


def diatomic_m1m2_sites(m1: str = "Fe", m2: str = "Co",
                        bond_length: float = 2.50, ligand: str = "N",
                        n_coord: int = 4, ligand_bond: float = 2.00) -> list:
    """B12-2 双原子催化剂（DAC）M₁–M₂ 位点生成器。

    M₁ 在原点（吸收中心），M₂ 沿 +z 距离 ``bond_length``；两者各带
    ``n_coord`` 个平面配体（默认 N₄）。返回包含 M₂ 与所有配体的
    ``[(x, y, z, element), ...]``，M₁ 不重复列出。
    """
    m1 = str(m1).capitalize()
    m2 = str(m2).capitalize()
    lig = str(ligand).capitalize()
    for sym in (m1, m2, lig):
        if sym not in ATOMIC_NUMBER:
            raise XAFSConfigError("DAC", f"未知元素 {sym!r}")
    if not (np.isfinite(bond_length) and bond_length > 0):
        raise XAFSConfigError("DAC", "bond_length 必须 > 0 Å")
    if not (np.isfinite(ligand_bond) and ligand_bond > 0):
        raise XAFSConfigError("DAC", "ligand_bond 必须 > 0 Å")
    n = int(n_coord)
    if n not in (2, 4, 6):
        raise XAFSConfigError("DAC", "n_coord 只支持 2/4/6")

    z2 = float(bond_length)
    sites = [(0.0, 0.0, z2, m2)]
    half = float(ligand_bond) / np.sqrt(2.0)

    def _add_ligands(z_center):
        if n == 2:
            pts = [(ligand_bond, 0.0), (-ligand_bond, 0.0)]
        else:
            pts = [(sx * half, sy * half)
                   for sx in (-1.0, 1.0) for sy in (-1.0, 1.0)]
            if n == 6:
                sites.extend([(0.0, 0.0, z_center + ligand_bond, lig),
                              (0.0, 0.0, z_center - ligand_bond, lig)])
        for (x, y) in pts:
            sites.append((float(x), float(y), float(z_center), lig))

    _add_ligands(0.0)
    _add_ligands(z2)
    return sites


def reconstruct_pair_distribution(
        r_grid: np.ndarray,
        n: float,
        r0: float,
        sigma2: float,
        c3: float = 0.0,
        c4: float = 0.0,
) -> dict:
    """由累积量展开重构非高斯近邻径向分布函数 g(R)（纯数学，无 FEFF）。

    数学:
        z = (R − r0)/σ,  σ = sqrt(max(σ², 1e-6))
        γ1 = C3/σ³,  γ2 = C4/σ⁴
        g_raw(R) = N·φ(z)·max(0, 1 + γ1/6·H3(z) + γ2/24·H4(z))
        H3 = z³ − 3z,  H4 = z⁴ − 6z² + 3
    再按梯形积分把总面积归一回 N（非负截断会改变面积，必须重归一）。

    物理:
        C3 > 0 表示右偏（长尾在高 R 侧）——正偏态分布的峰位
        r_mode 低于平均键长 r_mean。C4 描述超额峰度（展宽/尖峰）。
        该函数只回答“给定累积量时 g(R) 长什么样”，不替代拟合；
        C3/C4 必须来自数据约束，不能用它反推不可辨识的高阶矩。

    返回:
        {"r", "g_r", "r_mode", "r_mean", "skewness", "area"}
    """
    r = np.asarray(r_grid, dtype=float)
    if r.ndim != 1 or r.size < 3:
        raise XAFSConfigError(
            "RDF", f"r_grid 必须是一维且至少 3 点，收到形状 {r.shape}")
    if np.any(np.diff(r) <= 0):
        raise XAFSConfigError("RDF", "r_grid 必须严格单调递增")
    if not np.isfinite(n) or n <= 0:
        raise XAFSConfigError("RDF", f"配位数 N={n!r} 必须为正")
    sigma = float(np.sqrt(max(float(sigma2), 1e-6)))
    if not np.isfinite(r0):
        raise XAFSConfigError("RDF", f"r0={r0!r} 非有限")
    gamma1 = float(c3) / sigma ** 3
    gamma2 = float(c4) / sigma ** 4
    z = (r - float(r0)) / sigma
    phi = np.exp(-0.5 * z ** 2) / (np.sqrt(2.0 * np.pi) * sigma)
    h3 = z ** 3 - 3.0 * z
    h4 = z ** 4 - 6.0 * z ** 2 + 3.0
    shape = np.maximum(0.0, 1.0 + gamma1 / 6.0 * h3
                       + gamma2 / 24.0 * h4)
    g_raw = float(n) * phi * shape
    trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    area_raw = float(trapz(g_raw, r))
    if area_raw <= 0 or not np.isfinite(area_raw):
        raise XAFSConfigError(
            "RDF", f"累积量给不出非负 g(R)（截断后面积={area_raw!r}）；"
                   f"C3/C4 相对 σ 过大")
    g_norm = g_raw * (float(n) / area_raw)
    area = float(trapz(g_norm, r))
    return {
        "r": r,
        "g_r": g_norm,
        "r_mode": float(r[int(np.argmax(g_norm))]),
        "r_mean": float(r0),
        "skewness": gamma1,
        "area": area,
    }


def estimate_r_noise(chi_k, kmin: float, kmax: float, kweight: int = 2,
                     dk: float = 1.0,
                     r_noise_range=(15.0, 25.0)) -> float:
    """高 R 区噪声估计 ε_R（对标 IFEFFIT 的 epsilon_k / epsilon_r）。

    物理依据：EXAFS 的物理信号在 R>15 Å 早已衰减到噪声量级；该区间
    的 |χ(R)| 只剩下数据噪声与 FT 泄漏。对实部、虚部同时取均方根：

        ε_R = sqrt( 1/(2N) · Σ_{R∈[15,25]} (Re² + Im²) )

    为什么必须用 k^w χ(k) 做 FT：IFEFFIT 的 ε 与拟合时使用同一 k 加权；
    用未加权 χ(k) 会把噪声换算错一个 k^{2w} 的尺度，χ² 随之系统偏移。
    """
    lo, hi = float(r_noise_range[0]), float(r_noise_range[1])
    if hi <= lo:
        raise XAFSConfigError(
            "FIT", f"r_noise_range 上界 {hi} 必须大于下界 {lo}")
    # larch 的 xftf 默认 R 网格只到 ~kmax（本环境实测 9.97 Å），
    # 必须显式要求扩展到噪声区上界；否则 [15,25] 内一个点都没有。
    ft_full = athena.ft_to_r_space(
        np.asarray(chi_k.k), np.asarray(chi_k.chi), None,
        kmin=float(kmin), kmax=float(kmax), kweight=int(kweight),
        dk=float(dk), rmax_out=hi)
    r = np.asarray(ft_full.r, dtype=float)
    chi_r = np.asarray(ft_full.chi_r, dtype=complex)
    mask = (r >= lo) & (r <= hi)
    n = int(mask.sum())
    if n < 8:
        raise XAFSConfigError(
            "FIT", f"高 R 噪声区 [{lo},{hi}] Å 内仅 {n} 个 FT 点"
                   f"（FT R 网格上限 {float(r.max()):.2f} Å）；"
                   f"无法估计 ε_R")
    re = chi_r[mask].real
    im = chi_r[mask].imag
    eps = float(np.sqrt(np.mean(0.5 * (re ** 2 + im ** 2))))
    if not np.isfinite(eps) or eps < 0:
        raise XAFSConfigError("FIT", f"ε_R 计算得到非法值 {eps!r}")
    return eps


def shell_f_test(fit_simple, fit_complex) -> dict:
    """增减壳层的 Hamilton / F-test（对标 Artemis 的 F-test）。

    少壳层模型（p1 个变量）与多壳层模型（p2 > p1）必须使用**同一
    R 窗与同一数据**。IFEFFIT 定义：

        F = [(χ²₁ − χ²₂)/(p₂ − p₁)] / [χ²₂/(N_idp − p₂)]

    置信度 = F 分布的 CDF(F; dfn=p₂−p₁, dfd=N_idp−p₂)。默认
    ``confidence ≥ 0.95`` 且 ``dfd ≥ 1`` 才判显著；否则新增壳层
    只是拟合噪声。
    """
    from scipy.stats import f as f_dist

    p1 = int(getattr(fit_simple, "n_variables", 0))
    p2 = int(getattr(fit_complex, "n_variables", 0))
    n_idp = float(getattr(fit_complex, "n_independent",
                          getattr(fit_simple, "n_independent", 0.0)))
    dfn = p2 - p1
    dfd = float(n_idp - p2)

    out = {
        "F": 0.0, "dfn": int(dfn), "dfd": float(dfd),
        "confidence": 0.0, "significant": False,
        "p1": p1, "p2": p2, "n_independent": float(n_idp),
        "delta_chi_square": None, "method": None,
    }
    if dfn <= 0:
        out["note"] = "多壳层模型变量数未增加，F-test 不适用"
        return out
    if dfd < 1.0:
        out["note"] = f"dfd = N_idp − p2 = {dfd:.2f} < 1，无自由度"
        return out

    chi1 = getattr(fit_simple, "chi_square", None)
    chi2 = getattr(fit_complex, "chi_square", None)
    if chi1 is not None and chi2 is not None and chi2 > 0:
        delta = float(chi1 - chi2)
        F = (delta / dfn) / (float(chi2) / dfd)
        out["method"] = "chi_square"
        out["delta_chi_square"] = delta
    else:
        # 回退：R-factor 与 χ² 只差同一常数因子，比值形式仍成立。
        r1 = float(getattr(fit_simple, "r_factor", np.nan))
        r2 = float(getattr(fit_complex, "r_factor", np.nan))
        if not (np.isfinite(r1) and np.isfinite(r2) and r2 > 0):
            out["note"] = "两个模型都缺少 χ²，R-factor 也非法"
            return out
        F = ((r1 - r2) / dfn) / (r2 / dfd)
        out["method"] = "r_factor"

    if not np.isfinite(F):
        out["note"] = "F 非有限值"
        return out
    out["F"] = float(F)
    if F > 0:
        conf = float(f_dist.cdf(F, dfn, dfd))
        out["confidence"] = conf
        out["significant"] = bool(conf >= 0.95 and dfd >= 1.0)
    else:
        out["note"] = "多壳层模型未降低残差（F ≤ 0）"
    return out


def r_factor(data_chi_r, model_chi_r, mask) -> float:
    """
    绝对拟合优度。

    数学:
        R-factor = Σ|χ_data(R) − χ_th(R)|² / Σ|χ_data(R)|²

    判据: < 0.02 优秀，< 0.05 良好，> 0.10 模型不可用。
    """
    num = float(np.sum(np.abs(data_chi_r[mask] - model_chi_r[mask]) ** 2))
    den = float(np.sum(np.abs(data_chi_r[mask]) ** 2))
    return num / den if den > 0 else np.inf


# ================================================================
# 5. R 空间最小二乘拟合
# ================================================================

def _parse_kweights(kweight) -> tuple:
    """把单权 int 或序列统一成去重后保序的正整数元组。"""
    if isinstance(kweight, (int, np.integer)):
        vals = [int(kweight)]
    else:
        vals = [int(x) for x in kweight]
    if not vals:
        raise XAFSConfigError("FIT", "kweight 不能为空")
    if any(v < 0 for v in vals):
        raise XAFSConfigError("FIT", f"kweight 不能为负: {vals}")
    out = []
    for v in vals:
        if v not in out:
            out.append(v)
    return tuple(out)


def normalize_constraints(constraints, path_labels) -> list:
    """B14-1 催化配位化学计量约束（Mini-GDS）解析与校验。

    支持两种等价输入（**绝不使用 eval/exec**）：

    1. 结构化字典列表：
       - ``{"kind": "n_sum", "target": "Fe-Cl", "source": "Fe-O",
           "total": 6.0, "coeff": 1.0}``
         ⟹ N_target = total − coeff·N_source
       - ``{"kind": "n_ratio", "target": "Fe-Cl", "source": "Fe-O",
           "ratio": 1.5}`` ⟹ N_target = ratio·N_source
       - ``{"kind": "deltar_offset", "target": "Fe-Cl", "source": "Fe-O",
           "offset": 0.02}`` ⟹ ΔR_target = ΔR_source + offset
       - ``{"kind": "bounds", "path": "Fe-O", "param": "sigma2",
           "min": 0.002, "max": 0.010}``
    2. 字符串表达式（分号或逗号分隔）：
       - ``"N(Fe-O)+N(Fe-Cl)=6.0"``
       - ``"N(Fe-Cl)=1.5*N(Fe-O)"``
       - ``"dR(Fe-Cl)=dR(Fe-O)+0.02"``

    返回带 ``target_index`` / ``source_index`` 的规范化列表；未知路径
    标签、重复目标、非法数值或成环依赖一律抛 ``XAFSConfigError``。
    """
    import re

    labels = [str(x) for x in (path_labels or [])]
    if constraints is None or constraints == "" or constraints == []:
        return []
    if not labels:
        raise XAFSConfigError(
            "FIT", "constraints 需要 path_labels 才能按壳层名绑定")
    index_of: dict[str, list[int]] = {}
    for i, lab in enumerate(labels):
        index_of.setdefault(lab, []).append(i)

    def _resolve(label, what):
        key = str(label).strip()
        hits = index_of.get(key) or []
        if not hits:
            raise XAFSConfigError(
                "FIT", f"约束 {what} 引用了未知路径标签 {key!r}；"
                       f"可选 {labels}")
        if len(hits) > 1:
            raise XAFSConfigError(
                "FIT", f"约束 {what} 的标签 {key!r} 在路径中不唯一")
        return hits[0]

    if isinstance(constraints, str):
        raw = [c.strip() for c in re.split(r"[;,]", constraints)
               if c.strip()]
    elif isinstance(constraints, (list, tuple)):
        raw = list(constraints)
    else:
        raise XAFSConfigError(
            "FIT", f"constraints 必须是字符串或列表，收到 {type(constraints).__name__}")

    out: list[dict] = []
    for item in raw:
        if isinstance(item, str):
            spec = _parse_constraint_expression(item)
        elif isinstance(item, dict):
            spec = dict(item)
        else:
            raise XAFSConfigError(
                "FIT", f"约束项必须是字符串或字典，收到 {item!r}")
        kind = str(spec.get("kind", "")).strip().lower()
        if kind == "bounds":
            try:
                path = str(spec["path"]).strip()
                param = str(spec.get("param", "sigma2")).strip().lower()
                lo = float(spec["min"])
                hi = float(spec["max"])
            except (KeyError, TypeError, ValueError) as exc:
                raise XAFSConfigError(
                    "FIT", f"bounds 约束需要 path/param/min/max: {spec!r}") from exc
            idx = _resolve(path, "bounds.path")
            if param not in ("n", "deltar", "sigma2"):
                raise XAFSConfigError(
                    "FIT", f"bounds.param={param!r} 只支持 n / deltar / sigma2")
            if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
                raise XAFSConfigError(
                    "FIT", f"bounds 区间非法: min={lo!r}, max={hi!r}")
            out.append({"kind": "bounds", "path": path, "path_index": idx,
                        "param": param, "min": lo, "max": hi})
            continue
        if kind not in ("n_sum", "n_ratio", "deltar_offset"):
            raise XAFSConfigError(
                "FIT", f"未知约束 kind={kind!r}；只支持 n_sum / n_ratio / "
                       f"deltar_offset / bounds")
        try:
            target = str(spec["target"]).strip()
            source = str(spec["source"]).strip()
        except (KeyError, TypeError) as exc:
            raise XAFSConfigError(
                "FIT", f"{kind} 约束需要 target 与 source: {spec!r}") from exc
        ti = _resolve(target, f"{kind}.target")
        si = _resolve(source, f"{kind}.source")
        if ti == si:
            raise XAFSConfigError(
                "FIT", f"约束 target 与 source 同为 {target!r}，自约束无意义")
        rec = {"kind": kind, "target": target, "source": source,
               "target_index": ti, "source_index": si}
        try:
            if kind == "n_sum":
                total = float(spec["total"])
                coeff = float(spec.get("coeff", 1.0))
                if not np.isfinite(total) or total <= 0:
                    raise XAFSConfigError(
                        "FIT", f"n_sum.total={total!r} 必须 > 0")
                if not np.isfinite(coeff) or coeff <= 0:
                    raise XAFSConfigError(
                        "FIT", f"n_sum.coeff={coeff!r} 必须 > 0")
                rec.update(total=total, coeff=coeff)
            elif kind == "n_ratio":
                ratio = float(spec["ratio"])
                if not np.isfinite(ratio) or ratio < 0:
                    raise XAFSConfigError(
                        "FIT", f"n_ratio.ratio={ratio!r} 必须 ≥ 0")
                rec.update(ratio=ratio)
            else:
                offset = float(spec.get("offset", 0.0))
                if not np.isfinite(offset):
                    raise XAFSConfigError(
                        "FIT", f"deltar_offset.offset={offset!r} 非有限")
                rec.update(offset=offset)
        except KeyError as exc:
            raise XAFSConfigError(
                "FIT", f"{kind} 约束缺少字段 {exc.args[0]!r}") from exc
        out.append(rec)

    # ---- 同一参数被多条约束指向同一 target → 过约束，拒绝 ----
    seen_n: dict[int, int] = {}
    seen_dr: dict[int, int] = {}
    for rec in out:
        if rec["kind"] in ("n_sum", "n_ratio"):
            if rec["target_index"] in seen_n:
                raise XAFSConfigError(
                    "FIT", f"路径 {rec['target']!r} 的 N 被多条约束重复绑定")
            seen_n[rec["target_index"]] = rec["source_index"]
        elif rec["kind"] == "deltar_offset":
            if rec["target_index"] in seen_dr:
                raise XAFSConfigError(
                    "FIT", f"路径 {rec['target']!r} 的 ΔR 被多条约束重复绑定")
            seen_dr[rec["target_index"]] = rec["source_index"]
    for mapping, what in ((seen_n, "N"), (seen_dr, "ΔR")):
        for start in mapping:
            seen = set()
            k = start
            while k in mapping:
                if k in seen:
                    raise XAFSConfigError(
                        "FIT", f"{what} 约束成环（存在自我依赖），拒绝执行")
                seen.add(k)
                k = mapping[k]
    return out


def _parse_constraint_expression(text: str) -> dict:
    """解析一条约束字符串表达式（正则，不使用 eval）。"""
    import re

    s = str(text).strip()
    m = re.fullmatch(
        r"N\(\s*(?P<a>[^()]+?)\s*\)\s*\+\s*N\(\s*(?P<b>[^()]+?)\s*\)"
        r"\s*=\s*(?P<c>[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)", s)
    if m:
        return {"kind": "n_sum", "source": m.group("a"),
                "target": m.group("b"), "total": float(m.group("c")),
                "coeff": 1.0}
    m = re.fullmatch(
        r"N\(\s*(?P<t>[^()]+?)\s*\)\s*=\s*"
        r"(?P<r>[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*\*\s*"
        r"N\(\s*(?P<s>[^()]+?)\s*\)", s)
    if m:
        return {"kind": "n_ratio", "target": m.group("t"),
                "source": m.group("s"), "ratio": float(m.group("r"))}
    m = re.fullmatch(
        r"dR\(\s*(?P<t>[^()]+?)\s*\)\s*=\s*dR\(\s*(?P<s>[^()]+?)\s*\)"
        r"(?:\s*(?P<sign>[+-])\s*(?P<off>(?:\d+\.?\d*|\.\d+)"
        r"(?:[eE][-+]?\d+)?))?", s, re.I)
    if m:
        off = float(m.group("off") or 0.0)
        if m.group("sign") == "-":
            off = -off
        return {"kind": "deltar_offset", "target": m.group("t"),
                "source": m.group("s"), "offset": off}
    raise XAFSConfigError(
        "FIT", f"无法解析约束表达式 {text!r}；支持 "
               f"'N(A)+N(B)=C' / 'N(B)=r*N(A)' / 'dR(B)=dR(A)+off'")


def fit_r_space(chi_k, ft, datfiles, rmin: float, rmax: float,
                kmin: float = 3.0, kmax: float = 12.0, kweight=2,
                s02: float | None = None, vary_delta_e0: bool = True,
                vary_n=False, vary_c3=False, vary_s02=False,
                vary_sigma2=True, sigma2_init: float = 0.005,
                deltar_lim: tuple = (-0.3, 0.3),
                path_degen: list | None = None,
                path_labels: list | None = None,
                path_groups: list | None = None,
                path_group_names: list | None = None,
                vary_weights: bool = False,
                weight_init: list | None = None,
                sigma2_groups: list | None = None,
                sigma2_group_names: list | None = None,
                ms_links: list | None = None,
                isotropic_deltar: bool = False,
                vary_c4=False,
                bvs_restraint: dict | None = None,
                path_links: list | None = None,
                constraints=None,
                fit_space: str = "r",
                wt_sigma: float = 1.5,
                wt_rmin: float | None = None,
                wt_rmax: float | None = None,
                wt_kmin: float | None = None,
                wt_kmax: float | None = None):
    """
    在 R 空间对 χ(R) 的实部与虚部同时最小二乘。

    数学:
        min Σ_j [ (Re χ_data(R_j) − Re χ_th(R_j))²
                + (Im χ_data(R_j) − Im χ_th(R_j))² ]

    拟合变量: 每条路径的 ΔR、σ²；全局共享的 ΔE₀；
    可选每条路径的 N（配位数）。
    **只拟合 |χ(R)| 会丢相位**：实测用幅度拟合时 ΔE₀ 不可识别，
    Jacobian 直接奇异（Singular matrix）。

    vary_n:
        False / None   N 用 FEFF 晶体学简并度（固定）
        True           所有路径的 N 都作为变量
        [0, 2, ...]    只把指定下标的路径 N 设为变量（如只放开第一壳层）
    vary_c3: 同 vary_n 的开关语义，放开各路径三阶累积量 C₃（Å³）。
        仅用于近邻分布显著偏斜的体系（纳米颗粒表面张力、强非谐）；
        C₃ 与 ΔR 部分相关，常规体系不要开（会吃掉 R 的精度）。
    vary_sigma2: True / False / 路径下标列表。False 时 σ² 全部固定在
        sigma2_init——标定（calibrate_s02.py）用它把问题降成
        (S₀², ΔE₀) 两变量近线性拟合，否则 4 变量在窄窗内过参数化
        （实测 Fe 箔 kw=1/3 直接塌到参数边界、Rf 0.63，Kw 敏感）。
    sigma2_groups: 可选，列表的列表，如 [[0, 1], [2]] 表示路径 0、1
        共享一个 σ² 变量，路径 2 独立。用于“多壳层共享 σ²”的简化模型：
        减少自由参数、避免 σ² 与 ΔR 过度相关。组内路径输出相同的
        σ² 和 σ²_err。与 vary_sigma2 互斥（提供其一即可）。
    sigma2_init: 各路径 σ² 初值（默认 0.005）；σ² 被固定时即锁定值。
    deltar_lim: ΔR 的拟合边界 (下, 上)，默认 ±0.3。标定时收紧到
        ±0.02（RT 金属箔热膨胀位移量级），防止 ΔR 与 ΔE₀ 相 stray。
    vary_s02: **标定模式专用**（calibrate_s02.py）——N 固定为晶体学
        简并度、S₀² 作为唯一振幅变量放下。生产拟合严禁用它
        （S₀²·N 简并，会同时抬高两者）；s02=None 且 vary_s02 时
        用 1.0 作初值（标定不需要先验标定值）。
    path_degen: 可选，与 datfiles 等长；拟合前覆盖各路径的初始简并度。
        用于“合并壳层”的代表路径：把同壳多条 SS 路径的简并度加总到
        代表路径上，拟合出的 N 才是壳层总配位数，而不是单条路径的 N。
    path_labels: 可选，与 datfiles 等长；写入 PathFit.label 的壳层名。
    path_groups / path_group_names / vary_weights / weight_init:
        **mixture 相分数拟合**。path_groups 是每条路径所属组分下标
        （0..K-1，连续），vary_weights=True 时用 softmax 参数化放开
        K 个组分权重（恒正、和=1），振幅按 ``N_crystal·w_k`` 缩放。
        此时禁止同时 vary_n（N 与 w 只以乘积出现，严格简并）。
        输出 FitResult.mixture_weights = [{group, weight, weight_err}]；
        PathFit.degeneracy 是**有效简并度** N_crystal·w（不是配位数），
        n_err 由权重误差传播得到。权重误差来自协方差 + softmax 雅可比。
        weight_init 为各组分初始权重（默认均匀）。
    S₀² 与 N 在振幅上只以乘积 S₀²·N 出现，强关联不可同时放开——
    这里固定 S₀²（由标样定）、放开 N，是 Artemis 的标准做法。
    N 的实现：larch 的 _calc_chi 直接读 path.degen 属性乘进振幅
    （与 ΔR/σ² 走 params 字典不同），故每轮迭代直接改属性。

    ``isotropic_deltar``：True 时所有单散射路径共用唯一全局晶格膨胀
    系数 α（槽位 ``(-5, "alpha")``），ΔR_i = α·Reff_i，边界 ±0.08
    （±8% 晶格应变）。多壳层立方/高对称体系可把 M 个 ΔR 变量压缩成
    1 个，避免远壳层因 N_idp 不足被截断。
    ``vary_c4``：同 vary_c3 的开关语义，放开各路径四阶累积量 C₄
    （槽位 ``(i, "c4")``，边界 ±C4_LIMIT = ±2e-4 Å⁴）。

    ``bvs_restraint``（B11-2）：``{"target", "sigma", "r0_table"}``，
    把 ``(BVS−target)/sigma`` 追加到残差；只增加一行约束、不增加自由
    变量。输出 FitResult.bvs/bvs_target/bvs_chi2。
    ``fit_space``（B11-3）：``"r"``（默认，历史逐位一致）/ ``"q"``
    （R 窗反傅里叶）/ ``"wt"``（二维复数小波）。不同空间的 R-factor
    口径不同，不可横向比较。

    返回 FitResult。
    """
    from scipy.optimize import least_squares

    fit_space = str(fit_space or "r").strip().lower()
    if fit_space not in ("r", "q", "wt"):
        raise XAFSConfigError(
            "FIT", f"fit_space={fit_space!r} 只支持 'r' / 'q' / 'wt'")

    if s02 is None and not vary_s02:
        raise XAFSConfigError(
            "FIT", "s02 未提供——S₀² 必须显式传入（absorber.py 标定值，"
                   "Fe=0.803 箔标定）；无标定元素的 CN 定量一律拒绝")
    kw_tuple = _parse_kweights(kweight)
    multi_kw = len(kw_tuple) > 1
    # 主权重：生产约定 kw=2；不在元组中时取中位/第一个。
    primary_kw = 2 if 2 in kw_tuple else kw_tuple[len(kw_tuple) // 2]
    paths = [load_feff_path(f) for f in datfiles]
    n_path = len(paths)
    if path_degen is not None:
        if len(path_degen) != n_path:
            raise XAFSConfigError(
                "FIT", f"path_degen 长度 {len(path_degen)} != 路径数 {n_path}")
        for path, degen0 in zip(paths, path_degen):
            path.degen = float(degen0)
    if path_labels is not None and len(path_labels) != n_path:
        raise XAFSConfigError(
            "FIT", f"path_labels 长度 {len(path_labels)} != 路径数 {n_path}")
    # 晶体学简并度快照：权重模式下每轮迭代按 N_crystal·w 重新写 degen。
    crystal_degen = np.asarray([float(p.degen) for p in paths], dtype=float)

    # ---- 受约束多散射（MS）：路径参数由母单散射壳层按几何比例驱动 ----
    # 关键：MS 路径**不建立任何自由变量槽位**，n_var 不变，因此不会因
    # 加 MS 而触发 N_idp 过拟合截断。
    ms_by_path: dict[int, tuple] = {}
    if ms_links:
        for link in ms_links:
            idx = int(link["ms_path_index"])
            parent = int(link["parent_shell_index"])
            if not (0 <= idx < n_path and 0 <= parent < n_path):
                raise XAFSConfigError(
                    "FIT", f"ms_links 下标越界: ms={idx}, parent={parent}")
            if idx == parent or idx in ms_by_path:
                raise XAFSConfigError(
                    "FIT", f"ms_links 路径 {idx} 重复或指向自身")
            if parent in ms_by_path:
                raise XAFSConfigError(
                    "FIT", f"MS 母壳层 {parent} 不能又是受约束 MS 路径")
            ms_by_path[idx] = (parent, float(link["reff_ratio"]),
                               float(link["degen_ratio"]))
    ms_path_set = set(ms_by_path)

    # ---- σ² 共享组：多条路径共享一个 σ² 变量，减少自由参数 ----
    sigma2_group_of_path = {}   # path_idx -> group_idx
    n_s2_group = 0
    if sigma2_groups and ms_path_set:
        # MS 路径的 σ² 由母壳层驱动，不参与共享组变量。名称必须与
        # 保留下来的组同步过滤——按原下标截断会在中间组被删空时错位。
        filtered = []
        filtered_names = []
        for gi, group in enumerate(sigma2_groups):
            keep = [int(pi) for pi in group if int(pi) not in ms_path_set]
            if keep:
                filtered.append(keep)
                if sigma2_group_names is not None \
                        and gi < len(sigma2_group_names):
                    filtered_names.append(sigma2_group_names[gi])
        sigma2_groups = filtered
        if sigma2_group_names is not None:
            sigma2_group_names = filtered_names
    if sigma2_groups:
        n_s2_group = len(sigma2_groups)
        if sigma2_group_names is not None \
                and len(sigma2_group_names) != n_s2_group:
            raise XAFSConfigError(
                "FIT", f"sigma2_group_names 长度 {len(sigma2_group_names)} "
                       f"!= 组数 {n_s2_group}")
        for gi, group in enumerate(sigma2_groups):
            for pi in group:
                pi = int(pi)
                if pi < 0 or pi >= n_path:
                    raise XAFSConfigError(
                        "FIT", f"sigma2_groups 路径下标 {pi} 超出范围 0–{n_path-1}")
                if pi in sigma2_group_of_path:
                    raise XAFSConfigError(
                        "FIT", f"路径 {pi} 出现在多个 σ² 共享组中，"
                               f"每条路径只能属于一个组")
                sigma2_group_of_path[pi] = gi
    sigma2_group_names = (list(sigma2_group_names)
                          if sigma2_group_names is not None
                          else [f"s2grp{g}" for g in range(n_s2_group)])

    # ---- B13-2 单数据集内路径绑定（异核共占位双分支）----
    link_deltar: dict[int, int] = {}
    link_sigma2: dict[int, int] = {}
    link_n: dict[int, tuple[int, float]] = {}
    for li, link in enumerate(path_links or []):
        if not isinstance(link, dict):
            raise XAFSConfigError("FIT", f"path_links[{li}] 必须是字典")
        try:
            tgt = int(link["target"])
            src = int(link["source"])
        except (KeyError, TypeError, ValueError) as exc:
            raise XAFSConfigError(
                "FIT", f"path_links[{li}] 需要整数 target/source") from exc
        for val, what in ((tgt, "target"), (src, "source")):
            if val < 0 or val >= n_path:
                raise XAFSConfigError(
                    "FIT", f"path_links[{li}].{what}={val} 超出路径范围 "
                           f"0–{n_path - 1}")
        if tgt == src:
            raise XAFSConfigError(
                "FIT", f"path_links[{li}] 的 target==source，自绑定无约束")
        if tgt in ms_path_set:
            raise XAFSConfigError(
                "FIT", f"path_links[{li}] 的 target 是受约束 MS 路径")
        if link.get("share_deltar", True):
            if tgt in link_deltar:
                raise XAFSConfigError(
                    "FIT", f"路径 {tgt} 被多条 share_deltar 绑定")
            link_deltar[tgt] = src
        if link.get("share_sigma2", True):
            if tgt in link_sigma2:
                raise XAFSConfigError(
                    "FIT", f"路径 {tgt} 被多条 share_sigma2 绑定")
            link_sigma2[tgt] = src
        ratio = link.get("degen_ratio")
        if ratio is not None:
            ratio = float(ratio)
            if not np.isfinite(ratio) or ratio <= 0:
                raise XAFSConfigError(
                    "FIT", f"path_links[{li}].degen_ratio={ratio!r} 必须 > 0")
            if tgt in link_n:
                raise XAFSConfigError(
                    "FIT", f"路径 {tgt} 被多条 degen_ratio 绑定")
            link_n[tgt] = (src, ratio)
    if link_deltar and isotropic_deltar:
        raise XAFSConfigError(
            "FIT", "path_links 的 share_deltar 与 isotropic_deltar 互斥："
                   "全局 α 已经把所有 ΔR 绑在一起")
    if link_sigma2 and sigma2_groups:
        for tgt in link_sigma2:
            if tgt in sigma2_group_of_path:
                raise XAFSConfigError(
                    "FIT", f"路径 {tgt} 同时被 share_sigma2 与 σ² 共享组约束")

    # ---- B14-1 催化化学计量约束（Mini-GDS）----
    fit_labels = [str(path_labels[i]) if path_labels is not None
                  else Path(datfiles[i]).name for i in range(n_path)]
    con_list = normalize_constraints(constraints, fit_labels)
    con_n: dict[int, dict] = {}
    con_dr: dict[int, dict] = {}
    for rec in con_list:
        if rec["kind"] in ("n_sum", "n_ratio"):
            con_n[int(rec["target_index"])] = rec
        elif rec["kind"] == "deltar_offset":
            con_dr[int(rec["target_index"])] = rec
    if con_dr and isotropic_deltar:
        raise XAFSConfigError(
            "FIT", "deltar_offset 约束与 isotropic_deltar 互斥")
    for tgt in set(con_n) & set(link_n):
        raise XAFSConfigError("FIT", f"路径 {tgt} 的 N 同时被 path_links 与 constraints 绑定")
    for tgt in set(con_dr) & set(link_deltar):
        raise XAFSConfigError("FIT", f"路径 {tgt} 的 ΔR 同时被 path_links 与 constraints 绑定")

    for rec in con_list:
        if rec["kind"] != "bounds" and int(rec["target_index"]) in ms_path_set:
            raise XAFSConfigError(
                "FIT", f"约束目标 {rec['target']!r} 是受约束 MS 路径，"
                       f"其 N/ΔR 由母壳层驱动，不能再被化学计量约束")

    # ---- 哪些路径的 N 作为变量 ----
    if vary_n is True:
        n_paths = set(range(n_path))
    elif vary_n:
        n_paths = {int(i) for i in vary_n}
        bad = n_paths - set(range(n_path))
        if bad:
            raise XAFSConfigError(
                "FIT", f"vary_n 索引 {sorted(bad)} 超出路径范围 0–{n_path - 1}")
    else:
        n_paths = set()
    # MS 路径的 N 由母壳层驱动；即使用户传 True 也不给它独立变量。
    n_paths -= ms_path_set
    # B13-2/B14-1：被比例/守恒约束绑定的 N 不是自由变量。
    n_paths -= set(link_n) | set(con_n)
    c3_paths = set(range(n_path)) if vary_c3 is True else \
        ({int(i) for i in vary_c3} if vary_c3 else set())
    c3_paths -= ms_path_set
    c4_paths = set(range(n_path)) if vary_c4 is True else \
        ({int(i) for i in vary_c4} if vary_c4 else set())
    c4_paths -= ms_path_set
    bad = c4_paths - set(range(n_path))
    if bad:
        raise XAFSConfigError(
            "FIT", f"vary_c4 索引 {sorted(bad)} 超出路径范围 0–{n_path - 1}")

    # ---- mixture 权重（vary_weights）----
    weight_slot_ids = []
    weight_groups = None
    weight_names = None
    n_group = 0
    w0 = None
    if vary_weights:
        if n_paths:
            raise XAFSConfigError(
                "FIT", "vary_weights 与 vary_n 不能同时放开：振幅只以 "
                       "N·w 的乘积出现，两者同时自由会严格简并")
        if path_groups is None:
            raise XAFSConfigError(
                "FIT", "vary_weights 需要 path_groups（每条路径的组分下标）")
        if len(path_groups) != n_path:
            raise XAFSConfigError(
                "FIT", f"path_groups 长度 {len(path_groups)} != 路径数 {n_path}")
        weight_groups = np.asarray([int(g) for g in path_groups], dtype=int)
        if weight_groups.min() < 0:
            raise XAFSConfigError("FIT", "path_groups 不能含负数")
        group_ids = sorted(set(weight_groups.tolist()))
        if group_ids != list(range(len(group_ids))):
            raise XAFSConfigError(
                "FIT", f"path_groups 必须是 0..K-1 连续下标，收到 {group_ids}")
        n_group = len(group_ids)
        if n_group < 2:
            raise XAFSConfigError(
                "FIT", "vary_weights 至少需要 2 个组分（单组分权重恒为 1）")
        if weight_init is None:
            w0 = np.full(n_group, 1.0 / n_group)
        else:
            w0 = np.asarray(weight_init, dtype=float)
            if w0.shape != (n_group,):
                raise XAFSConfigError(
                    "FIT", f"weight_init 长度 {len(w0)} != 组分数 {n_group}")
            if np.any(w0 <= 0) or abs(float(w0.sum()) - 1.0) > 1e-6:
                raise XAFSConfigError("FIT", "weight_init 必须为正且和为 1")
        weight_names = (list(path_group_names) if path_group_names is not None
                        else [f"group{k}" for k in range(n_group)])
        if len(weight_names) != n_group:
            raise XAFSConfigError(
                "FIT", f"path_group_names 长度 {len(weight_names)} != 组分数 {n_group}")

    # ---- 变量槽位: 每条路径 (ΔR[, σ²][, N][, C₃])，末尾可选 ΔE₀ ----
    slots = []                       # [(path_idx, 'deltar'|'sigma2'|'n'|'c3'|'e0'), ...]
    x0, lower, upper = [], [], []
    if vary_sigma2 is True:
        s2_paths = set(range(n_path))
    else:
        s2_paths = ({int(i) for i in vary_sigma2} if vary_sigma2 else set())
        bad = s2_paths - set(range(n_path))
        if bad:
            raise XAFSConfigError(
                "FIT", f"vary_sigma2 索引 {sorted(bad)} 超出路径范围 0–{n_path - 1}")
    # 与 σ² 共享组互斥
    if sigma2_groups and vary_sigma2 is not True:
        raise XAFSConfigError(
            "FIT", "sigma2_groups 与 vary_sigma2 列表模式不可同时使用；"
                   "请只提供 sigma2_groups，或只提供 vary_sigma2")

    # 1) 每条 SS 路径的 ΔR（MS 路径由母壳层比例驱动，不建槽位）
    # B7-1：isotropic_deltar 时不再逐路径建 ΔR，只建唯一全局 α。
    use_iso_deltar = bool(isotropic_deltar) and n_path > 1
    if use_iso_deltar:
        slots.append((-5, "alpha"))
        x0.append(0.0)
        lower.append(-0.08); upper.append(0.08)
    else:
        for i, path in enumerate(paths):
            if i in ms_path_set or i in link_deltar or i in con_dr:
                continue
            slots.append((i, "deltar"))
            x0.append(0.0); lower.append(deltar_lim[0]); upper.append(deltar_lim[1])

    # 2) σ²：共享组（每组一个变量）或独立
    s2_independent = set()           # 拥有独立 σ² 变量的路径下标
    if sigma2_groups:
        for gi, group in enumerate(sigma2_groups):
            slots.append((-4, gi))   # -4 标记 σ² 共享组，gi=组索引
            x0.append(float(sigma2_init)); lower.append(1e-4); upper.append(0.04)
    else:
        for i in sorted(s2_paths - ms_path_set - set(link_sigma2)):
            slots.append((i, "sigma2"))
            x0.append(float(sigma2_init)); lower.append(1e-4); upper.append(0.04)
            s2_independent.add(i)

    # 3) N
    for i in sorted(n_paths):
        # 必须按**该路径自己的**晶体学简并度设初值与边界；旧版误用
        # 循环残留变量（上一条路径的 degen），多壳层时会把近壳 N 的
        # 上界错误压到 3×远壳简并度，直接截断真实配位数。
        degen = float(crystal_degen[i])
        slots.append((i, "n"))
        x0.append(degen)
        lower.append(max(0.5, 0.1 * degen))
        upper.append(3.0 * degen)
    # 4) C3
    for i in sorted(c3_paths):
        slots.append((i, "c3"))
        x0.append(0.0); lower.append(-C3_LIMIT); upper.append(C3_LIMIT)
    # 4.5) C4（B7-3）
    for i in sorted(c4_paths):
        slots.append((i, "c4"))
        x0.append(0.0); lower.append(-C4_LIMIT); upper.append(C4_LIMIT)
    # 5) S0²（标定模式）
    if vary_s02:
        slots.append((-2, "s02"))
        x0.append(1.0 if s02 is None else float(s02))
        lower.append(0.3); upper.append(1.3)
    # 6) mixture 权重 softmax logit
    if vary_weights:
        for k in range(n_group - 1):
            p0 = np.log(max(float(w0[k]), 1e-6) / max(float(w0[-1]), 1e-6))
            slots.append((-3, k))
            x0.append(float(np.clip(p0, -8.0, 8.0)))
            lower.append(-8.0); upper.append(8.0)
            weight_slot_ids.append(len(slots) - 1)
    # 7) 全局 ΔE₀
    if vary_delta_e0:
        slots.append((-1, "e0"))
        x0.append(0.0); lower.append(-15.0); upper.append(15.0)
    n_var = len(slots)
    x0 = np.asarray(x0, dtype=float)
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    slot_index = {s: idx for idx, s in enumerate(slots)}

    # ---- B14-1 约束边界：bounds 覆盖 + n_sum 的物理非负性 ----
    for rec in con_list:
        if rec["kind"] != "bounds":
            continue
        key = (int(rec["path_index"]), rec["param"])
        if key not in slot_index:
            raise XAFSConfigError(
                "FIT", f"bounds 目标 {rec['path']}:{rec['param']} 不是自由变量"
                       f"（可能被约束/固定），无法施加边界")
        j = slot_index[key]
        lower[j] = max(float(lower[j]), float(rec["min"]))
        upper[j] = min(float(upper[j]), float(rec["max"]))
        if lower[j] >= upper[j]:
            raise XAFSConfigError(
                "FIT", f"bounds 与现有边界不相容: {rec['path']}:{rec['param']} "
                       f"[{lower[j]}, {upper[j]}]")
    for rec in con_list:
        if rec["kind"] != "n_sum":
            continue
        key = (int(rec["source_index"]), "n")
        if key not in slot_index:
            continue
        j = slot_index[key]
        ub = float(rec["total"]) / float(rec["coeff"])
        if ub <= float(lower[j]):
            raise XAFSConfigError(
                "FIT", f"n_sum 约束不可行：total/coeff={ub:.3f} ≤ "
                       f"{rec['source']} 的 N 下界 {lower[j]:.3f}")
        upper[j] = min(float(upper[j]), ub)

    # 自由变量名：与 covariance/correlation 行列一一对应，报告层据此解释
    # 强相关（N–σ²–ΔE0）而不是把形式误差当独立精度。
    _kind_names = {"deltar": "ΔR", "sigma2": "σ²", "n": "N",
                   "c3": "C3", "c4": "C4", "s02": "S0²", "e0": "ΔE0",
                   "alpha": "lattice_alpha"}
    parameter_names = []
    for (pi, pkind) in slots:
        if pi == -4:
            parameter_names.append(f"s2grp:{sigma2_group_names[pkind]}")
        elif pi >= 0:
            plabel = (path_labels[pi] if path_labels is not None
                      else Path(datfiles[pi]).name)
            parameter_names.append(f"{plabel}:{_kind_names[pkind]}")
        elif pi == -3:
            parameter_names.append(f"w:{weight_names[pkind]}")
        else:
            parameter_names.append(_kind_names[pkind])

    def _link_root(mapping, key, what):
        """沿 target→source 链走到根；成环时显式报错。"""
        seen = set()
        k = int(key)
        while k in mapping:
            if k in seen:
                raise XAFSConfigError("FIT", f"{what} 绑定成环: 路径 {k}")
            seen.add(k)
            k = int(mapping[k])
        return k

    def _n_ratio_chain(key):
        """N_key = prod·N_root，返回 (root, prod)。"""
        prod = 1.0
        k = int(key)
        seen = set()
        while k in link_n:
            if k in seen:
                raise XAFSConfigError("FIT", f"degen_ratio 绑定成环: 路径 {k}")
            seen.add(k)
            k, ratio = link_n[k]
            k = int(k)
            prod *= float(ratio)
        return k, prod

    def _n_constraint_chain(key):
        """返回 (root, coeff, const)，满足 N_key = const + coeff·N_root。"""
        if key in con_n:
            rec = con_n[key]
            root, coeff, const = _n_constraint_chain(int(rec["source_index"]))
            if rec["kind"] == "n_sum":
                return root, -float(rec["coeff"]) * coeff, \
                    float(rec["total"]) - float(rec["coeff"]) * const
            return root, float(rec["ratio"]) * coeff, float(rec["ratio"]) * const
        return int(key), 1.0, 0.0

    def _dr_constraint_chain(key):
        """返回 (root, offset)，满足 ΔR_key = ΔR_root + offset。"""
        off = 0.0
        k = int(key)
        seen = set()
        while k in con_dr:
            if k in seen:
                raise XAFSConfigError("FIT", f"deltar_offset 约束成环: 路径 {k}")
            seen.add(k)
            rec = con_dr[k]
            off += float(rec["offset"])
            k = int(rec["source_index"])
        return k, off

    def softmax_weights(x):
        """当前 logit → 组分权重（最后一个组分为参考，logit=0）。"""
        p = np.asarray([x[slot_index[(-3, k)]] for k in range(n_group - 1)]
                       + [0.0], dtype=float)
        p = p - p.max()
        e = np.exp(p)
        return e / e.sum()

    # 模型必须在**数据的 k 网格**上计算，否则 FT 后相位对不齐。
    ft_k = np.asarray(chi_k.k)
    data_r = ft.chi_r
    r_grid = ft.r
    mask = (r_grid >= rmin) & (r_grid <= rmax)
    if mask.sum() < 8:
        raise XAFSConfigError("FIT", f"R 区间 [{rmin},{rmax}] 内仅 {mask.sum()} 点")

    # ---- 多 kweight 联合拟合：逐权重准备实验谱与尺度因子 ----
    # 单权模式直接用传入 ft 的复数谱，保证与历史结果逐位一致。
    # 多权模式对同一 χ(k) 用相同 (kmin,kmax,dk) 重新 FT；分支残差必须
    # 除以该分支实验谱的 RMS 范数 s_w，否则 kw=3 的幅度（~k² 倍）会
    # 完全统治最小二乘，退化成单权 kw=3 拟合。
    data_r_by_kw: dict[int, np.ndarray] = {}
    scale_by_kw: dict[int, float] = {}
    for w in kw_tuple:
        if not multi_kw and w == int(getattr(ft, "kweight", primary_kw)):
            data_w = data_r
        else:
            data_w = athena.ft_to_r_space(
                ft_k, np.asarray(chi_k.chi), None, kmin, kmax, w, ft.dk).chi_r
        data_r_by_kw[w] = data_w
        if multi_kw:
            scale_by_kw[w] = float(np.sqrt(np.mean(np.abs(data_w[mask]) ** 2)))
            if not np.isfinite(scale_by_kw[w]) or scale_by_kw[w] <= 0:
                scale_by_kw[w] = 1.0
        else:
            scale_by_kw[w] = 1.0

    # ---- B11-3 多空间：q（反 FT）或 wt（二维复数小波）的数据侧 ----
    data_space_by_kw: dict = {}
    space_mask_by_kw: dict = {}
    space_scale_by_kw: dict = {}
    wt_r_window = None
    if fit_space == "q":
        for w in kw_tuple:
            bt = athena.back_ft_to_q_space(
                r_grid, data_r_by_kw[w], rmin=float(rmin), rmax=float(rmax),
                dr=float(getattr(ft, "dk", 0.1)))
            q = np.asarray(bt.q, dtype=float)
            mq = (q >= float(kmin)) & (q <= float(kmax))
            if int(mq.sum()) < 8:
                raise XAFSConfigError(
                    "FIT", f"q 空间 [{kmin},{kmax}] 内仅 {int(mq.sum())} 点")
            data_space_by_kw[w] = (q, np.asarray(bt.chi_q, dtype=complex))
            space_mask_by_kw[w] = mq
            s = float(np.sqrt(np.mean(np.abs(bt.chi_q[mq]) ** 2)))
            space_scale_by_kw[w] = s if np.isfinite(s) and s > 0 else 1.0
    elif fit_space == "wt":
        import larch_wt

        rlo = float(wt_rmin if wt_rmin is not None else rmin)
        rhi = float(wt_rmax if wt_rmax is not None else rmax)
        if rhi <= rlo:
            raise XAFSConfigError("FIT", f"wt R 窗无效: [{rlo},{rhi}]")
        wt_r_window = (rlo, rhi)
        wt_rmax_eff = max(rhi, float(rmax))
        for w in kw_tuple:
            wm = larch_wt.continuous_wavelet_transform(
                chi_k, sigma=float(wt_sigma), kmin=float(kmin),
                kmax=min(float(kmax), float(np.max(chi_k.k))),
                kweight=int(w), rmax=wt_rmax_eff, normalize=False)
            kc = np.asarray(wm.k_centers, dtype=float)
            klo = float(wt_kmin if wt_kmin is not None else kc.min())
            khi = float(wt_kmax if wt_kmax is not None else kc.max())
            mk = (kc >= klo) & (kc <= khi)
            mr = (np.asarray(wm.r) >= rlo) & (np.asarray(wm.r) <= rhi)
            if int(mk.sum()) < 3 or int(mr.sum()) < 3:
                raise XAFSConfigError(
                    "FIT", f"WT 窗口过小: k={int(mk.sum())}, R={int(mr.sum())}")
            data_space_by_kw[w] = (kc, np.asarray(wm.r)[mr],
                                   np.asarray(wm.complex_map)[np.ix_(mk, mr)])
            space_mask_by_kw[w] = (mk, mr)
            block = data_space_by_kw[w][2]
            s = float(np.sqrt(np.mean(np.abs(block) ** 2)))
            space_scale_by_kw[w] = s if np.isfinite(s) and s > 0 else 1.0

    def _model_in_space(model_k, w):
        """把模型 χ(k) 变到当前 fit_space 的数据网格（复数）。"""
        if fit_space == "r":
            return _ft_model(model_k, w)
        if fit_space == "q":
            model_r = _ft_model(model_k, w)
            bt = athena.back_ft_to_q_space(
                r_grid, model_r, rmin=float(rmin), rmax=float(rmax),
                dr=float(getattr(ft, "dk", 0.1)))
            q_data = data_space_by_kw[w][0]
            return (np.interp(q_data, bt.q, np.asarray(bt.chi_q).real)
                    + 1j * np.interp(q_data, bt.q, np.asarray(bt.chi_q).imag))
        import larch_wt
        wm = larch_wt.continuous_wavelet_transform(
            _chi_k_like(model_k), sigma=float(wt_sigma), kmin=float(kmin),
            kmax=min(float(kmax), float(np.max(_chi_k_like(model_k).k))),
            kweight=int(w), rmax=max(wt_r_window[1], float(rmax)),
            normalize=False)
        return np.asarray(wm.complex_map)

    def _chi_k_like(model_k):
        from types import SimpleNamespace
        return SimpleNamespace(k=np.asarray(ft_k, dtype=float),
                               chi=np.asarray(model_k, dtype=float))

    # ---- B13-2/B14-1 派生参数求值：沿 link/constraint 链回代到自由变量 ----
    def _dr_of(x, key, seen=None):
        seen = set() if seen is None else seen
        if key in seen:
            raise XAFSConfigError("FIT", f"ΔR 绑定成环: 路径 {key}")
        seen.add(key)
        if key in link_deltar:
            return _dr_of(x, int(link_deltar[key]), seen)
        if key in con_dr:
            rec = con_dr[key]
            return _dr_of(x, int(rec["source_index"]), seen) + float(rec["offset"])
        sk = (int(key), "deltar")
        return float(x[slot_index[sk]]) if sk in slot_index else 0.0

    def _s2_of(x, key, seen=None):
        seen = set() if seen is None else seen
        if key in seen:
            raise XAFSConfigError("FIT", f"σ² 绑定成环: 路径 {key}")
        seen.add(key)
        if key in link_sigma2:
            return _s2_of(x, int(link_sigma2[key]), seen)
        if key in sigma2_group_of_path:
            gi = sigma2_group_of_path[key]
            sk = (-4, gi)
            return float(x[slot_index[sk]]) if sk in slot_index \
                else float(sigma2_init)
        sk = (int(key), "sigma2")
        return float(x[slot_index[sk]]) if sk in slot_index \
            else float(sigma2_init)

    def _n_of(x, key, seen=None):
        seen = set() if seen is None else seen
        if key in seen:
            raise XAFSConfigError("FIT", f"N 绑定成环: 路径 {key}")
        seen.add(key)
        if key in link_n:
            src, ratio = link_n[key]
            return float(ratio) * _n_of(x, int(src), seen)
        if key in con_n:
            rec = con_n[key]
            src_val = _n_of(x, int(rec["source_index"]), seen)
            if rec["kind"] == "n_sum":
                return float(rec["total"]) - float(rec["coeff"]) * src_val
            return float(rec["ratio"]) * src_val
        sk = (int(key), "n")
        return (float(x[slot_index[sk]]) if sk in slot_index
                else float(crystal_degen[key]))

    def apply_params(x):
        """把变量向量写回路径：ΔR/σ²/E0 走 params 通道，N 改 degen 属性。"""
        by_path = {}
        for (i, kind), v in zip(slots, x):
            v = float(v)
            if i == -4:
                # σ² 共享组：kind 是组索引，写到组内所有路径
                for pi in sigma2_groups[kind]:
                    by_path.setdefault(pi, {})["sigma2"] = v
            elif i >= 0:
                by_path.setdefault(i, {})[kind] = v
            elif i == -2:
                by_path.setdefault(-2, {})["s02"] = v
            elif i == -1:
                by_path.setdefault(-1, {})["e0"] = v
        # 派生参数注入 by_path，后续统一按“自由变量”写回路径。
        for tgt in set(link_deltar) | set(con_dr):
            by_path.setdefault(tgt, {})["deltar"] = _dr_of(x, tgt)
        for tgt in link_sigma2:
            by_path.setdefault(tgt, {})["sigma2"] = _s2_of(x, tgt)
        for tgt in set(link_n) | set(con_n):
            by_path.setdefault(tgt, {})["n"] = _n_of(x, tgt)
        s02_val = by_path.get(-2, {}).get("s02", s02)
        alpha_val = (float(x[slot_index[(-5, "alpha")]])
                     if use_iso_deltar else None)
        w = softmax_weights(x) if vary_weights else None
        for i, path in enumerate(paths):
            v = by_path.get(i, {})
            if w is not None:
                # 相分数拟合：有效简并度 = 晶体学 N × 组分权重。
                path.degen = float(crystal_degen[i]
                                   * w[int(weight_groups[i])])
            elif "n" in v:
                path.degen = v["n"]
            if alpha_val is not None and i not in ms_path_set:
                dr_value = alpha_val * float(path.reff)
            else:
                dr_value = v.get("deltar", 0.0)
            set_path_parameters(
                path, s02=s02_val,
                delta_r=dr_value,
                sigma2=v.get("sigma2", sigma2_init),
                delta_e0=by_path.get(-1, {}).get("e0", 0.0),
                third=v.get("c3"), fourth=v.get("c4"))
        # 受约束 MS：ΔR/σ² 按路径长度比例、N 按母壳层有效简并度缩放。
        for idx, (parent, reff_ratio, degen_ratio) in ms_by_path.items():
            if use_iso_deltar:
                dr_p = float(alpha_val) * float(paths[parent].reff)
            else:
                dr_p = float(x[slot_index[(parent, "deltar")]])
            if parent in sigma2_group_of_path:
                s2_key = (-4, sigma2_group_of_path[parent])
                s2_p = (float(x[slot_index[s2_key]])
                        if s2_key in slot_index else float(sigma2_init))
            elif (parent, "sigma2") in slot_index:
                s2_p = float(x[slot_index[(parent, "sigma2")]])
            else:
                s2_p = float(sigma2_init)
            if vary_weights:
                n_parent = float(crystal_degen[parent]
                                 * w[int(weight_groups[parent])])
            elif parent in n_paths:
                n_parent = float(x[slot_index[(parent, "n")]])
            else:
                n_parent = float(crystal_degen[parent])
            paths[idx].degen = float(degen_ratio * n_parent)
            set_path_parameters(
                paths[idx], s02=s02_val, delta_r=dr_p * reff_ratio,
                sigma2=s2_p * reff_ratio,
                delta_e0=by_path.get(-1, {}).get("e0", 0.0), third=None)

    def _ft_model(model_k, w):
        """把模型 χ(k) 用权重 w 变到 R 网格并插值到数据 R 轴。"""
        model_ft = athena.ft_to_r_space(ft_k, model_k, None, kmin, kmax,
                                        w, ft.dk)
        return np.interp(r_grid, model_ft.r, model_ft.chi_r.real) \
            + 1j * np.interp(r_grid, model_ft.r, model_ft.chi_r.imag)

    def model_chi_r(w=None):
        """单权兼容入口：返回主权重的模型 χ(R)。"""
        return _ft_model(exafs_model_k(paths, ft_k),
                         primary_kw if w is None else w)

    # ---- B11-2 BVS 软约束配置 ----
    bvs_cfg = None
    if bvs_restraint is not None:
        if not isinstance(bvs_restraint, dict):
            raise XAFSConfigError("FIT", "bvs_restraint 必须是字典")
        try:
            bvs_target_val = float(bvs_restraint["target"])
            bvs_sigma_val = float(bvs_restraint.get("sigma", 0.1))
        except (KeyError, TypeError, ValueError) as exc:
            raise XAFSConfigError(
                "FIT", "bvs_restraint 需要 target 与 sigma（正数）") from exc
        if not np.isfinite(bvs_target_val):
            raise XAFSConfigError("FIT", "bvs_restraint.target 非有限")
        if not np.isfinite(bvs_sigma_val) or bvs_sigma_val <= 0:
            raise XAFSConfigError("FIT", "bvs_restraint.sigma 必须 > 0")
        bvs_cfg = {
            "target": bvs_target_val,
            "sigma": bvs_sigma_val,
            "r0_table": bvs_restraint.get("r0_table"),
            "absorber": bvs_restraint.get("absorber"),
            "b": float(bvs_restraint.get("b", 0.37)),
        }

    def _current_bvs():
        rows = []
        for i, path in enumerate(paths):
            label = (path_labels[i] if path_labels is not None
                     else Path(datfiles[i]).name)
            # ΔR 必须读 params 里的当前值：larch 的 set_path_parameters
            # 只改 params、不同步 path.deltar 属性（实测恒为 0）。
            # 若用属性，BVS 会在 Reff 处计算，约束与最终 R 脱节。
            try:
                dr_now = float(
                    path.params[path.pathpar_name("deltar")].value)
            except (AttributeError, KeyError):
                dr_now = float(getattr(path, "deltar", 0.0) or 0.0)
            rows.append({
                "label": label,
                "N": float(path.degen),
                "R": float(path.reff) + dr_now,
                "ligand": _path_ligand(label, bvs_cfg["absorber"]),
            })
        return bond_valence_sum(
            rows, r0_table=bvs_cfg["r0_table"], b=bvs_cfg["b"],
            absorber=bvs_cfg["absorber"])

    def residual(x):
        apply_params(x)
        model_k = exafs_model_k(paths, ft_k)
        parts = []
        for w in kw_tuple:
            if fit_space == "r":
                model_w = _ft_model(model_k, w)
                data_w = data_r_by_kw[w]
                dr = (model_w.real - data_w.real)[mask]
                di = (model_w.imag - data_w.imag)[mask]
                if multi_kw:
                    dr = dr / scale_by_kw[w]
                    di = di / scale_by_kw[w]
            elif fit_space == "q":
                model_w = _model_in_space(model_k, w)
                data_w = data_space_by_kw[w][1]
                mq = space_mask_by_kw[w]
                dr = (model_w.real - data_w.real)[mq] / space_scale_by_kw[w]
                di = (model_w.imag - data_w.imag)[mq] / space_scale_by_kw[w]
            else:
                mk, mr = space_mask_by_kw[w]
                model_full = _model_in_space(model_k, w)
                model_w = model_full[np.ix_(mk, mr)]
                data_w = data_space_by_kw[w][2]
                dr = ((model_w.real - data_w.real)
                      / space_scale_by_kw[w]).ravel()
                di = ((model_w.imag - data_w.imag)
                      / space_scale_by_kw[w]).ravel()
            parts.append(dr)
            parts.append(di)
        if bvs_cfg is not None:
            bvs_now = _current_bvs()
            parts.append(np.asarray(
                [(bvs_now - bvs_cfg["target"]) / bvs_cfg["sigma"]]))
        return np.concatenate(parts)

    solution = least_squares(residual, x0, bounds=(lower, upper))

    dof = max(len(solution.fun) - n_var, 1)
    covariance = np.linalg.pinv(solution.jac.T @ solution.jac) * (2 * solution.cost / dof)
    errors = np.sqrt(np.abs(np.diag(covariance)))
    # 相关矩阵：cov_ij / (σ_i σ_j)。变量数可到几十，JSON 里直接存完整
    # 矩阵（报告层只摘最强相关对）；零误差槽位不产生 NaN。
    denom = np.outer(errors, errors)
    correlation = np.divide(covariance, denom,
                            out=np.zeros_like(covariance), where=denom > 0)
    np.fill_diagonal(correlation, 1.0)

    # ---- mixture 权重与误差（delta method，softmax 雅可比）----
    mixture_weights = None
    w_final = None
    w_err_final = None
    if vary_weights:
        w_final = softmax_weights(solution.x)
        cov_pp = covariance[np.ix_(weight_slot_ids, weight_slot_ids)]
        # G[k, j] = ∂w_k/∂p_j；参考组分 logit 固定为 0（第 K 行无自由列）
        G = np.zeros((n_group, n_group - 1), dtype=float)
        for k in range(n_group):
            for j in range(n_group - 1):
                G[k, j] = w_final[k] * ((1.0 if k == j else 0.0) - w_final[j])
        cov_w = G @ cov_pp @ G.T
        w_err_final = np.sqrt(np.abs(np.diag(cov_w)))
        mixture_weights = [
            {"group": weight_names[k], "weight": float(w_final[k]),
             "weight_err": float(w_err_final[k])}
            for k in range(n_group)]

    # ---- 预计算 σ² 共享组的最优值与误差 ----
    s2_group_values = {}
    s2_group_errors = {}
    if sigma2_groups:
        for gi in range(n_s2_group):
            slot_key = (-4, gi)
            if slot_key in slot_index:
                s2_group_values[gi] = float(
                    solution.x[slot_index[slot_key]])
                s2_group_errors[gi] = float(
                    errors[slot_index[slot_key]])

    # 先写回最优解，MS 路径的缩放值也才是最终值。
    apply_params(solution.x)
    alpha_val = (float(solution.x[slot_index[(-5, "alpha")]])
                 if use_iso_deltar else None)
    alpha_err = (float(errors[slot_index[(-5, "alpha")]])
                 if use_iso_deltar else None)
    # ---- B13-2/B14-1 输出层：派生参数误差沿同一链传播 ----
    def _dr_err_of(key, seen=None):
        seen = set() if seen is None else seen
        if key in seen:
            raise XAFSConfigError("FIT", f"ΔR 误差链成环: 路径 {key}")
        seen.add(key)
        if key in link_deltar:
            return _dr_err_of(int(link_deltar[key]), seen)
        if key in con_dr:
            return _dr_err_of(int(con_dr[key]["source_index"]), seen)
        sk = (int(key), "deltar")
        return float(errors[slot_index[sk]]) if sk in slot_index else None

    def _s2_err_of(key, seen=None):
        seen = set() if seen is None else seen
        if key in seen:
            raise XAFSConfigError("FIT", f"σ² 误差链成环: 路径 {key}")
        seen.add(key)
        if key in link_sigma2:
            return _s2_err_of(int(link_sigma2[key]), seen)
        if key in sigma2_group_of_path:
            gi = sigma2_group_of_path[key]
            return s2_group_errors.get(gi)
        sk = (int(key), "sigma2")
        return float(errors[slot_index[sk]]) if sk in slot_index else None

    def _n_err_of(key, seen=None):
        seen = set() if seen is None else seen
        if key in seen:
            raise XAFSConfigError("FIT", f"N 误差链成环: 路径 {key}")
        seen.add(key)
        if vary_weights:
            return float(crystal_degen[key]
                         * w_err_final[int(weight_groups[key])])
        if key in link_n:
            src, ratio = link_n[key]
            e = _n_err_of(int(src), seen)
            return None if e is None else float(ratio) * e
        if key in con_n:
            rec = con_n[key]
            e = _n_err_of(int(rec["source_index"]), seen)
            if e is None:
                return None
            coeff = (float(rec["coeff"]) if rec["kind"] == "n_sum"
                     else float(rec["ratio"]))
            return abs(coeff) * e
        sk = (int(key), "n")
        return float(errors[slot_index[sk]]) if sk in slot_index else None

    fit_paths_list = []
    for i, path in enumerate(paths):
        if i in ms_by_path:
            parent, reff_ratio, degen_ratio = ms_by_path[i]
            if use_iso_deltar:
                dr_p = float(alpha_val) * float(paths[parent].reff)
                dr_err_p = float(alpha_err) * float(paths[parent].reff)
            else:
                dr_p = float(solution.x[slot_index[(parent, "deltar")]])
                dr_err_p = float(errors[slot_index[(parent, "deltar")]])
            if parent in sigma2_group_of_path:
                s2_key = (-4, sigma2_group_of_path[parent])
                s2_p = float(solution.x[slot_index[s2_key]])
                s2_err_p = float(errors[slot_index[s2_key]])
            elif (parent, "sigma2") in slot_index:
                s2_p = float(solution.x[slot_index[(parent, "sigma2")]])
                s2_err_p = float(errors[slot_index[(parent, "sigma2")]])
            else:
                s2_p, s2_err_p = float(sigma2_init), None
            if vary_weights:
                n_err_p = float(crystal_degen[parent]
                                * w_err_final[int(weight_groups[parent])])
            elif parent in n_paths:
                n_err_p = float(errors[slot_index[(parent, "n")]])
            else:
                n_err_p = None
            fit_paths_list.append(PathFit(
                label=(path_labels[i] if path_labels is not None
                       else Path(datfiles[i]).name),
                reff=float(path.reff),
                degeneracy=float(path.degen),
                delta_r=float(dr_p * reff_ratio),
                delta_r_err=float(dr_err_p * reff_ratio),
                sigma2=float(s2_p * reff_ratio),
                sigma2_err=(None if s2_err_p is None
                            else float(s2_err_p * reff_ratio)),
                n_err=(None if n_err_p is None
                       else float(n_err_p * degen_ratio)),
                sigma2_at_bound=bool(s2_p * reff_ratio <= 5e-4),
                nleg=int(getattr(path, "nleg", 2)),
                ms_parent=(path_labels[parent]
                           if path_labels is not None else str(parent))))
            continue
        if i in sigma2_group_of_path or i in link_sigma2 \
                or (i, "sigma2") in slot_index:
            s2_val = _s2_of(solution.x, i)
            s2_err = _s2_err_of(i)
            s2_at_bound = s2_val <= 5e-4
        else:
            s2_val = float(sigma2_init)
            s2_err = None
            s2_at_bound = False
        c3_val = (float(solution.x[slot_index[(i, "c3")]])
                  if i in c3_paths else None)
        c3_err_val = (float(errors[slot_index[(i, "c3")]])
                      if i in c3_paths else None)
        c4_val = (float(solution.x[slot_index[(i, "c4")]])
                  if i in c4_paths else None)
        c4_err_val = (float(errors[slot_index[(i, "c4")]])
                      if i in c4_paths else None)
        if use_iso_deltar:
            dr_val = float(alpha_val) * float(paths[i].reff)
            dr_err_val = float(alpha_err) * float(paths[i].reff)
        else:
            dr_val = _dr_of(solution.x, i)
            dr_err_val = _dr_err_of(i)
        fit_paths_list.append(PathFit(
            label=(path_labels[i] if path_labels is not None
                   else Path(datfiles[i]).name),
            reff=float(path.reff),
            # 权重模式：有效简并度 N_crystal·w（不是配位数）；普通模式：拟合 N
            degeneracy=(float(crystal_degen[i] * w_final[int(weight_groups[i])])
                        if vary_weights else float(path.degen)),
            delta_r=dr_val,
            delta_r_err=dr_err_val,
            sigma2=float(s2_val),
            sigma2_err=s2_err,
            n_err=((float(crystal_degen[i]
                          * w_err_final[int(weight_groups[i])]))
                   if vary_weights else _n_err_of(i)),
            c3=c3_val,
            c3_err=c3_err_val,
            c3_at_bound=bool(c3_val is not None
                             and abs(c3_val) >= 0.8 * C3_LIMIT),
            c4=c4_val,
            c4_err=c4_err_val,
            c4_at_bound=bool(c4_val is not None
                             and abs(c4_val) >= 0.8 * C4_LIMIT),
            sigma2_at_bound=s2_at_bound))

    # R-factor 在最优解处重算：必须给**所有**路径都设上最优参数。
    # 旧版只设了 paths[0]，多路径时其余路径仍是残差函数最后一次迭代的值，
    # 导致双壳层拟合的 R-factor 不可信。
    apply_params(solution.x)
    model_k_final = exafs_model_k(paths, ft_k)

    def _space_rfactor(w):
        if fit_space == "r":
            return float(r_factor(data_r_by_kw[w],
                                  _ft_model(model_k_final, w), mask))
        if fit_space == "q":
            mq = space_mask_by_kw[w]
            data_w = data_space_by_kw[w][1][mq]
            model_w = _model_in_space(model_k_final, w)[mq]
            ones = np.ones(mq.sum(), dtype=bool)
            return float(r_factor(data_w, model_w, ones))
        mk, mr = space_mask_by_kw[w]
        data_w = data_space_by_kw[w][2].ravel()
        model_w = _model_in_space(model_k_final, w)[np.ix_(mk, mr)].ravel()
        ones = np.ones(data_w.size, dtype=bool)
        return float(r_factor(data_w, model_w, ones))

    r_factors_by_kweight = {
        int(w): _space_rfactor(w) for w in kw_tuple}

    # ---- B11-2 BVS 最终值与约束 χ² ----
    bvs_value = None
    bvs_target_out = None
    bvs_chi2_out = None
    if bvs_cfg is not None:
        bvs_value = float(_current_bvs())
        bvs_target_out = float(bvs_cfg["target"])
        bvs_chi2_out = float(
            ((bvs_value - bvs_target_out) / bvs_cfg["sigma"]) ** 2)

    # ---- B12-1 AIC / BIC（以数据残差 RSS 计，不含 BVS 约束行） ----
    data_fun = (solution.fun[:-1] if bvs_cfg is not None
                else solution.fun)
    rss = float(np.sum(data_fun ** 2))
    n_obs = int(data_fun.size)
    aic = bic = None
    if n_obs > 0 and rss > 0:
        aic = float(n_obs * np.log(rss / n_obs) + 2.0 * n_var)
        bic = float(n_obs * np.log(rss / n_obs) + n_var * np.log(n_obs))

    # ---- B4-1 统计量：高 R 噪声 ε_R、χ²、χ²_ν ----
    # 用主 kweight 分支计算，与 r_factor 的汇报口径一致；多权时各分支
    # 共用 N_idp（权重不增加独立信息量，见 B1-1 契约）。
    n_idp = independent_points(kmin, kmax, rmin, rmax)
    epsilon_r = None
    chi_square = None
    reduced_chi_square = None
    try:
        epsilon_r = estimate_r_noise(
            chi_k, kmin, kmax, kweight=int(primary_kw), dk=float(ft.dk))
    except XAFSConfigError:
        epsilon_r = None
    if epsilon_r is not None and epsilon_r > 0:
        data_main = data_r_by_kw[int(primary_kw)]
        model_main = _ft_model(model_k_final, int(primary_kw))
        n_grid = int(mask.sum())
        diff2 = float(np.sum(np.abs(data_main[mask] - model_main[mask]) ** 2))
        chi_square = (n_idp * diff2
                      / (2.0 * max(n_grid, 1) * epsilon_r ** 2))
        reduced_chi_square = chi_square / max(n_idp - n_var, 1.0)

    return FitResult(
        paths=fit_paths_list,
        s02=(float(solution.x[slot_index[(-2, "s02")]])
             if vary_s02 else float(s02)),
        s02_err=(float(errors[slot_index[(-2, "s02")]]) if vary_s02 else None),
        delta_e0=float(solution.x[slot_index[(-1, "e0")]]) if vary_delta_e0 else 0.0,
        delta_e0_err=float(errors[slot_index[(-1, "e0")]]) if vary_delta_e0 else 0.0,
        r_factor=r_factors_by_kweight[int(primary_kw)],
        n_independent=float(n_idp),
        n_variables=n_var,
        epsilon_r=epsilon_r,
        chi_square=chi_square,
        reduced_chi_square=reduced_chi_square,
        lattice_strain=alpha_val,
        lattice_strain_err=alpha_err,
        kmin=float(kmin), kmax=float(kmax),
        rmin=float(rmin), rmax=float(rmax),
        mixture_weights=mixture_weights,
        parameter_names=parameter_names,
        correlation=correlation,
        sigma2_groups=sigma2_groups,
        sigma2_group_names=sigma2_group_names,
        kweights=kw_tuple,
        r_factors_by_kweight=r_factors_by_kweight,
        model_r=np.asarray(r_grid, dtype=float),
        model_chi_r=np.asarray(
            _ft_model(model_k_final, primary_kw), dtype=complex),
        model_k=np.asarray(model_k_final, dtype=float),
        model_k_grid=np.asarray(ft_k, dtype=float),
        aic=aic, bic=bic,
        bvs=bvs_value, bvs_target=bvs_target_out,
        bvs_chi2=bvs_chi2_out,
        fit_space=fit_space)


# ================================================================
# 5b. 多数据集全局受约束联合拟合（B8-1）
# ================================================================

# ħ²/(2 k_B) 换算到 amu·Å²·K：
#   ħ = 1.054571817e-34 J·s, k_B = 1.380649e-23 J/K
#   1 amu = 1.66053906660e-27 kg, 1 Å² = 1e-20 m²
#   → 24.254316 amu·Å²·K
HBAR2_OVER_2KB = 24.254316


def einstein_sigma2(temperature_K: float, theta_e: float,
                    reduced_mass_amu: float,
                    sigma2_static: float = 0.0) -> float:
    """爱因斯坦模型的热无序 σ²(T)（量子谐振子精确式）。

    数学:
        σ²(T) = σ²_static + (ħ² / (2 μ k_B Θ_E)) · coth(Θ_E / (2T))
        ħ²/(2 k_B) = 24.254316 amu·Å²·K

    物理:
        T → 0 时 coth → 1，保留零点振动；T ≫ Θ_E 时
        σ² ≈ σ²_static + 2·(24.254316/(μΘ_E))·(T/Θ_E)，线性于 T。
        高温近似不能替代本式——变温联合拟合的全部信息来自
        coth 的非线性温度依赖。

    参数:
        temperature_K   : 样品温度 (K)，必须 > 0
        theta_e         : 爱因斯坦温度 (K)，必须 > 0
        reduced_mass_amu: 键折合质量 (amu)，必须 > 0
        sigma2_static   : 静态无序 (Å²)，必须 ≥ 0

    返回:
        σ²(T) (Å²)
    """
    T = float(temperature_K)
    theta = float(theta_e)
    mu = float(reduced_mass_amu)
    static = float(sigma2_static)
    if not np.isfinite(T) or T <= 0:
        raise XAFSConfigError("FIT", f"温度 T={temperature_K!r} 必须 > 0 K")
    if not np.isfinite(theta) or theta <= 0:
        raise XAFSConfigError(
            "FIT", f"爱因斯坦温度 Θ_E={theta_e!r} 必须 > 0 K")
    if not np.isfinite(mu) or mu <= 0:
        raise XAFSConfigError(
            "FIT", f"折合质量 μ={reduced_mass_amu!r} 必须 > 0 amu")
    if not np.isfinite(static) or static < 0:
        raise XAFSConfigError(
            "FIT", f"σ²_static={sigma2_static!r} 不能为负")
    x = theta / (2.0 * T)
    thermal = (HBAR2_OVER_2KB / (mu * theta)) / np.tanh(x)
    return float(static + thermal)


def profile_chi2_contour(fit_func, x_name: str, y_name: str,
                         x_values, y_values, fixed: dict | None = None,
                         sigma_levels=(2.30, 6.18, 11.83)) -> dict:
    """B12-1 二维 χ² 剖面扫描（供诊断图等高线）。

    对 ``x_name × y_name`` 的每个网格点调用 ``fit_func(**params)``，
    其余参数由 ``fixed`` 固定。``fit_func`` 可返回标量 χ²、``FitResult``
    或 ``(chi2, ...)`` 元组。``levels`` 是相对最小值的 Δχ² 阈值
    （二维 1σ/2σ/3σ 分别 ≈2.30/6.18/11.83）。
    """
    xv = np.asarray(x_values, dtype=float)
    yv = np.asarray(y_values, dtype=float)
    if xv.ndim != 1 or yv.ndim != 1 or xv.size < 2 or yv.size < 2:
        raise XAFSConfigError("PROFILE", "x_values/y_values 必须是一维且 ≥2 点")
    fixed = dict(fixed or {})
    grid = np.full((yv.size, xv.size), np.nan, dtype=float)
    for iy, yv_i in enumerate(yv):
        for ix, xv_i in enumerate(xv):
            params = dict(fixed)
            params[x_name] = float(xv_i)
            params[y_name] = float(yv_i)
            try:
                res = fit_func(**params)
            except Exception:  # noqa: BLE001 - 失败网格保留 NaN，不伪造
                continue
            if hasattr(res, "chi_square"):
                val = res.chi_square
            elif isinstance(res, (tuple, list)):
                val = res[0]
            else:
                val = res
            try:
                grid[iy, ix] = float(val)
            except (TypeError, ValueError):
                grid[iy, ix] = np.nan
    finite = grid[np.isfinite(grid)]
    levels = ([float(np.min(finite) + d) for d in sigma_levels]
              if finite.size else [])
    return {"x": xv, "y": yv, "chi2_grid": grid,
            "min_chi2": float(np.min(finite)) if finite.size else None,
            "levels": levels, "delta_chi2": list(sigma_levels)}


def fit_multi_datasets(
    datasets: list[dict],
    shared_bonds: list[dict] | None = None,
    einstein_models: list[dict] | None = None,
    share_e0="auto",
    debye_models: list[dict] | None = None,
) -> dict:
    """多数据集全局受约束联合拟合（Co-refinement）。

    支持三类痛点:
      1. 双金属双吸收边（如 Pt L3 + Co K）：用 ``shared_bonds`` 强制
         异核键满足 R_A-B = R_B-A、σ²_A-B = σ²_B-A；
      2. 变温多谱（如 80 K + 350 K）：用 ``einstein_models`` 以
         σ²_static 与 Θ_E 两个全局参数驱动每个温度下的 σ²(T)；
      3. 关联 Debye 模型（B14-3）：用 ``debye_models`` 以
         (σ²_static, Θ_D) 驱动同一温度下**随 R_eff 变化**的 σ²(T,R)，
         近壳运动相关压低、远壳回到非关联值。

    数学:
        每个数据集 d 在自身 mask_d 内的实/虚部残差除以该数据集实验谱
        的 RMS 范数 s_d = sqrt(mean(|χ_data,d(R)|²))，再拼接成全局残差：
            min Σ_d ||(χ_model,d − χ_data,d)/s_d||²
        这样 80 K 与 350 K 的振幅差异不会让某一个数据集统治最小二乘。

    参数:
        datasets : 每项
            {"name", "chi_k", "ft", "datfiles", "rmin", "rmax",
             "kmin", "kmax", "kweight", "s02", "vary_n",
             "path_degen", "path_labels", "temperature_K"}
        shared_bonds : 每项
            {"target": (ds_b, path_b), "source": (ds_a, path_a),
             "share_deltar": bool, "share_sigma2": bool,
             "degen_ratio": float | None}
            目标路径不建独立 ΔR/σ²/N 槽位，复用源路径；N 按
            degen_ratio 缩放（None 表示不共享 N）。
        einstein_models : 每项
            {"paths": [(ds, path), ...], "reduced_mass_amu": float,
             "theta_e_init": float = 350.0}
            组内路径的 σ² 全部由 (σ²_static, Θ_E) 两个全局变量按
            各自数据集温度计算。
        debye_models : 每项
            {"paths": [(ds, path), ...], "mass_absorber": float,
             "mass_scatterer": float, "theta_d_init": float = 400.0,
             "k_debye": float | None}
            组内路径的 σ² 由 (σ²_static, Θ_D) 按
            ``correlated_debye_sigma2(T, Θ_D, R_eff, m_a, m_s)`` 实时
            计算；与 einstein/share_sigma2 对同一路径互斥。

    返回:
        {"datasets": [FitResult, ...], "n_independent_total": float,
         "n_variables_total": int, "global_r_factor": float,
         "einstein_results": [{"theta_e", "theta_e_err",
                               "sigma2_static", "sigma2_static_err"}, ...],
         "debye_results": [{"theta_d", "theta_d_err",
                            "sigma2_static", "sigma2_static_err",
                            "paths", ...}, ...],
         "parameter_names": list, "correlation": np.ndarray}

    注意:
        ``share_e0`` 控制 ΔE₀ 的分组（B12-1）：
        True 所有数据集共享一个 ΔE₀（历史行为）；False 每个数据集独立；
        "auto"（默认）按 ``datasets[di]["edge_group"]`` 分组，缺省用
        ``edge`` 字符串，再缺省全部共享。同组共享一个 ΔE₀。
        双吸收边场景应让不同边落在不同 edge_group，而不是硬凑一个
        无物理意义的共同值。
    """
    from scipy.optimize import least_squares

    if not datasets:
        raise XAFSConfigError("FIT", "datasets 不能为空")
    n_ds = len(datasets)
    required = ("chi_k", "ft", "datfiles", "rmin", "rmax",
                "kmin", "kmax", "kweight")
    for di, d in enumerate(datasets):
        for key in required:
            if d.get(key) is None:
                raise XAFSConfigError(
                    "FIT", f"数据集 {di} 缺少字段 {key!r}")
        if not d["datfiles"]:
            raise XAFSConfigError("FIT", f"数据集 {di} 的 datfiles 为空")
        if d.get("s02") is None:
            raise XAFSConfigError(
                "FIT", f"数据集 {di} 未提供 s02；S₀² 必须显式给定")

    # ---- 载入路径并冻结晶体学简并度 ----
    ds_paths: list[list] = []
    crystal_degen: list[list[float]] = []
    for di, d in enumerate(datasets):
        paths = [load_feff_path(f) for f in d["datfiles"]]
        degs = d.get("path_degen")
        if degs is not None:
            if len(degs) != len(paths):
                raise XAFSConfigError(
                    "FIT", f"数据集 {di} 的 path_degen 长度 {len(degs)} "
                           f"!= 路径数 {len(paths)}")
            for p, deg in zip(paths, degs):
                p.degen = float(deg)
        labels = d.get("path_labels")
        if labels is not None and len(labels) != len(paths):
            raise XAFSConfigError(
                "FIT", f"数据集 {di} 的 path_labels 长度 {len(labels)} "
                       f"!= 路径数 {len(paths)}")
        ds_paths.append(paths)
        crystal_degen.append([float(p.degen) for p in paths])

    def _check_path(di, pi, what):
        if di < 0 or di >= n_ds:
            raise XAFSConfigError(
                "FIT", f"{what} 的数据集下标 {di} 超出 0–{n_ds - 1}")
        if pi < 0 or pi >= len(ds_paths[di]):
            raise XAFSConfigError(
                "FIT", f"{what} 的路径下标 {pi} 超出数据集 {di} 的 "
                       f"0–{len(ds_paths[di]) - 1}")

    # ---- 跨数据集路径绑定 ----
    bound_deltar: dict[tuple, tuple] = {}
    bound_sigma2: dict[tuple, tuple] = {}
    bound_n: dict[tuple, tuple] = {}
    for bi, bond in enumerate(shared_bonds or []):
        try:
            tgt = tuple(int(v) for v in bond["target"])
            src = tuple(int(v) for v in bond["source"])
        except (KeyError, TypeError, ValueError) as exc:
            raise XAFSConfigError(
                "FIT", f"shared_bonds[{bi}] 的 target/source 必须是 "
                       f"(数据集下标, 路径下标)") from exc
        if len(tgt) != 2 or len(src) != 2:
            raise XAFSConfigError(
                "FIT", f"shared_bonds[{bi}] 的 target/source 必须各含 2 个整数")
        _check_path(tgt[0], tgt[1], f"shared_bonds[{bi}].target")
        _check_path(src[0], src[1], f"shared_bonds[{bi}].source")
        if tgt == src:
            raise XAFSConfigError(
                "FIT", f"shared_bonds[{bi}] 的 target 与 source 相同；"
                       f"自绑定不产生约束")
        if bond.get("share_deltar"):
            if tgt in bound_deltar:
                raise XAFSConfigError(
                    "FIT", f"路径 {tgt} 被多个 share_deltar 绑定")
            bound_deltar[tgt] = src
        if bond.get("share_sigma2"):
            if tgt in bound_sigma2:
                raise XAFSConfigError(
                    "FIT", f"路径 {tgt} 被多个 share_sigma2 绑定")
            bound_sigma2[tgt] = src
        ratio = bond.get("degen_ratio")
        if ratio is not None:
            ratio = float(ratio)
            if not np.isfinite(ratio) or ratio <= 0:
                raise XAFSConfigError(
                    "FIT", f"shared_bonds[{bi}].degen_ratio={ratio!r} "
                           f"必须为正")
            if tgt in bound_n:
                raise XAFSConfigError(
                    "FIT", f"路径 {tgt} 被多个 degen_ratio 绑定")
            bound_n[tgt] = (src, ratio)

    # ---- 爱因斯坦变温组 ----
    einstein_groups = list(einstein_models or [])
    einstein_of: dict[tuple, int] = {}
    for gi, model in enumerate(einstein_groups):
        try:
            mu = float(model["reduced_mass_amu"])
            plist = model["paths"]
        except (KeyError, TypeError, ValueError) as exc:
            raise XAFSConfigError(
                "FIT", f"einstein_models[{gi}] 缺少 reduced_mass_amu/paths") from exc
        if not np.isfinite(mu) or mu <= 0:
            raise XAFSConfigError(
                "FIT", f"einstein_models[{gi}].reduced_mass_amu={mu!r} "
                       f"必须 > 0")
        if not plist:
            raise XAFSConfigError(
                "FIT", f"einstein_models[{gi}].paths 不能为空")
        for item in plist:
            if len(item) != 2:
                raise XAFSConfigError(
                    "FIT", f"einstein_models[{gi}].paths 每项必须是 "
                           f"(数据集下标, 路径下标)")
            di, pi = int(item[0]), int(item[1])
            _check_path(di, pi, f"einstein_models[{gi}].paths")
            key = (di, pi)
            if key in einstein_of:
                raise XAFSConfigError(
                    "FIT", f"路径 {key} 出现在多个爱因斯坦组中")
            if key in bound_sigma2:
                raise XAFSConfigError(
                    "FIT", f"路径 {key} 同时被 share_sigma2 与爱因斯坦模型约束")
            temp = datasets[di].get("temperature_K")
            if temp is None:
                raise XAFSConfigError(
                    "FIT", f"数据集 {di} 参与爱因斯坦模型但未给 temperature_K")
            if float(temp) <= 0:
                raise XAFSConfigError(
                    "FIT", f"数据集 {di} 的 temperature_K={temp!r} 必须 > 0")
            einstein_of[key] = gi

    # ---- 关联 Debye 变温组（B14-3）----
    # 与 einstein_models 的区别：爱因斯坦模型只按温度给出 σ²(T)，
    # 不区分键长；Debye 模型用关联因子 C(R)=1−sin(2k_D R)/(2k_D R)
    # 让近壳运动强相关（σ² 小）、远壳回到非关联值（σ² 大）。六方
    # 密堆/岩盐等固体中同一 Θ_D 必须同时约束近壳与远壳。
    debye_groups = list(debye_models or [])
    debye_of: dict[tuple, int] = {}
    for gi, model in enumerate(debye_groups):
        try:
            m_abs = float(model["mass_absorber"])
            m_sc = float(model["mass_scatterer"])
            plist = model["paths"]
        except (KeyError, TypeError, ValueError) as exc:
            raise XAFSConfigError(
                "FIT", f"debye_models[{gi}] 缺少 "
                       f"mass_absorber/mass_scatterer/paths") from exc
        if not (np.isfinite(m_abs) and m_abs > 0):
            raise XAFSConfigError(
                "FIT", f"debye_models[{gi}].mass_absorber={m_abs!r} 必须 > 0")
        if not (np.isfinite(m_sc) and m_sc > 0):
            raise XAFSConfigError(
                "FIT", f"debye_models[{gi}].mass_scatterer={m_sc!r} 必须 > 0")
        if not plist:
            raise XAFSConfigError(
                "FIT", f"debye_models[{gi}].paths 不能为空")
        for item in plist:
            if len(item) != 2:
                raise XAFSConfigError(
                    "FIT", f"debye_models[{gi}].paths 每项必须是 "
                           f"(数据集下标, 路径下标)")
            di, pi = int(item[0]), int(item[1])
            _check_path(di, pi, f"debye_models[{gi}].paths")
            key = (di, pi)
            if key in debye_of:
                raise XAFSConfigError(
                    "FIT", f"路径 {key} 出现在多个 Debye 组中")
            if key in einstein_of:
                raise XAFSConfigError(
                    "FIT", f"路径 {key} 同时被 Debye 与爱因斯坦模型约束")
            if key in bound_sigma2:
                raise XAFSConfigError(
                    "FIT", f"路径 {key} 同时被 share_sigma2 与 Debye 模型约束")
            temp = datasets[di].get("temperature_K")
            if temp is None:
                raise XAFSConfigError(
                    "FIT", f"数据集 {di} 参与 Debye 模型但未给 temperature_K")
            if float(temp) <= 0:
                raise XAFSConfigError(
                    "FIT", f"数据集 {di} 的 temperature_K={temp!r} 必须 > 0")
            debye_of[key] = gi

    # ---- 每数据集的 N 变量集合 ----
    ds_n_paths: list[set[int]] = []
    for di, d in enumerate(datasets):
        vn = d.get("vary_n", False)
        n_path = len(ds_paths[di])
        if vn is True:
            nset = set(range(n_path))
        elif vn:
            nset = {int(i) for i in vn}
            bad = nset - set(range(n_path))
            if bad:
                raise XAFSConfigError(
                    "FIT", f"数据集 {di} 的 vary_n 索引 {sorted(bad)} "
                           f"超出 0–{n_path - 1}")
        else:
            nset = set()
        ds_n_paths.append(nset)

    # ---- B12-1 ΔE₀ 分组（share_e0）----
    if share_e0 is True:
        e0_raw_groups = ["global"] * n_ds
    elif share_e0 is False:
        e0_raw_groups = [f"ds{di}" for di in range(n_ds)]
    elif str(share_e0).lower() == "auto":
        e0_raw_groups = []
        for di, d in enumerate(datasets):
            key = d.get("edge_group")
            if key is None:
                key = d.get("edge")
            e0_raw_groups.append("global" if key is None else str(key))
        if all(k == "global" for k in e0_raw_groups):
            e0_raw_groups = ["global"] * n_ds
    else:
        raise XAFSConfigError(
            "FIT", f"share_e0={share_e0!r} 只支持 True / False / 'auto'")
    e0_group_order = []
    for key in e0_raw_groups:
        if key not in e0_group_order:
            e0_group_order.append(key)
    e0_group_index = {key: gi for gi, key in enumerate(e0_group_order)}
    e0_group_of_ds = [e0_group_index[key] for key in e0_raw_groups]
    n_e0_groups = len(e0_group_order)

    def e0_slot_key(gi):
        # 单组时沿用历史键，保证旧回归的参数名/数值逐位不变。
        return ("global", "e0") if n_e0_groups == 1 else ("e0grp", gi)

    # ---- 全局变量槽位 ----
    slots: list[tuple] = []
    x0: list[float] = []
    lower: list[float] = []
    upper: list[float] = []

    def add_slot(key, init, lo, hi):
        slots.append(key)
        x0.append(float(init))
        lower.append(float(lo))
        upper.append(float(hi))

    for di in range(n_ds):
        for pi in range(len(ds_paths[di])):
            key = (di, pi)
            if key not in bound_deltar:
                add_slot((di, pi, "deltar"), 0.0, -0.3, 0.3)
            if key not in bound_sigma2 and key not in einstein_of \
                    and key not in debye_of:
                add_slot((di, pi, "sigma2"), 0.005, 1e-4, 0.04)
            if pi in ds_n_paths[di] and key not in bound_n:
                deg = crystal_degen[di][pi]
                add_slot((di, pi, "n"), deg,
                         max(0.5, 0.1 * deg), 3.0 * deg)
    for gi, model in enumerate(einstein_groups):
        add_slot(("einstein", gi, "sigma2_static"), 0.001, 0.0, 0.015)
        add_slot(("einstein", gi, "theta_e"),
                 float(model.get("theta_e_init", 350.0)), 80.0, 1200.0)
    for gi, model in enumerate(debye_groups):
        add_slot(("debye", gi, "sigma2_static"), 0.001, 0.0, 0.015)
        add_slot(("debye", gi, "theta_d"),
                 float(model.get("theta_d_init", 400.0)), 80.0, 1200.0)
    for gi in range(n_e0_groups):
        add_slot(e0_slot_key(gi), 0.0, -15.0, 15.0)

    x0_a = np.asarray(x0, dtype=float)
    lower_a = np.asarray(lower, dtype=float)
    upper_a = np.asarray(upper, dtype=float)
    slot_index = {s: i for i, s in enumerate(slots)}
    n_var = len(slots)

    # ---- 槽位解析（支持绑定链） ----
    def _slot_value(x, key):
        return float(x[slot_index[key]])

    def resolved_dr(x, key):
        seen = set()
        k = key
        while k in bound_deltar:
            if k in seen:
                raise XAFSConfigError(
                    "FIT", f"share_deltar 绑定成环: {k}")
            seen.add(k)
            k = bound_deltar[k]
        return _slot_value(x, (k[0], k[1], "deltar"))

    def resolved_sigma2(x, key):
        seen = set()
        k = key
        while True:
            if k in einstein_of:
                gi = einstein_of[k]
                T = float(datasets[k[0]]["temperature_K"])
                mu = float(einstein_groups[gi]["reduced_mass_amu"])
                return einstein_sigma2(
                    T, _slot_value(x, ("einstein", gi, "theta_e")), mu,
                    _slot_value(x, ("einstein", gi, "sigma2_static")))
            if k in debye_of:
                gi = debye_of[k]
                model = debye_groups[gi]
                T = float(datasets[k[0]]["temperature_K"])
                reff = float(ds_paths[k[0]][k[1]].reff)
                return correlated_debye_sigma2(
                    T, _slot_value(x, ("debye", gi, "theta_d")), reff,
                    float(model["mass_absorber"]),
                    float(model["mass_scatterer"]),
                    k_debye=model.get("k_debye"),
                    sigma2_static=_slot_value(
                        x, ("debye", gi, "sigma2_static")))
            if k in bound_sigma2:
                if k in seen:
                    raise XAFSConfigError(
                        "FIT", f"share_sigma2 绑定成环: {k}")
                seen.add(k)
                k = bound_sigma2[k]
                continue
            return _slot_value(x, (k[0], k[1], "sigma2"))

    def resolved_n(x, key):
        seen = set()
        k = key
        factor = 1.0
        while k in bound_n:
            if k in seen:
                raise XAFSConfigError("FIT", f"degen_ratio 绑定成环: {k}")
            seen.add(k)
            k, ratio = bound_n[k][0], bound_n[k][1]
            factor *= ratio
        slot = (k[0], k[1], "n")
        base = (_slot_value(x, slot) if slot in slot_index
                else crystal_degen[k[0]][k[1]])
        return factor * base

    def apply_params(x):
        for di, paths in enumerate(ds_paths):
            e0 = _slot_value(x, e0_slot_key(e0_group_of_ds[di]))
            s02 = float(datasets[di]["s02"])
            for pi, path in enumerate(paths):
                key = (di, pi)
                path.degen = resolved_n(x, key)
                set_path_parameters(
                    path, s02=s02,
                    delta_r=resolved_dr(x, key),
                    sigma2=resolved_sigma2(x, key),
                    delta_e0=e0, third=None)

    # ---- 各数据集的实验谱、mask 与 RMS 尺度 ----
    ds_ft_k = [np.asarray(d["chi_k"].k, dtype=float) for d in datasets]
    ds_data_r = [np.asarray(d["ft"].chi_r, dtype=complex) for d in datasets]
    ds_r_grid = [np.asarray(d["ft"].r, dtype=float) for d in datasets]
    ds_mask = []
    ds_scale = []
    for di, d in enumerate(datasets):
        m = (ds_r_grid[di] >= float(d["rmin"])) & \
            (ds_r_grid[di] <= float(d["rmax"]))
        if m.sum() < 8:
            raise XAFSConfigError(
                "FIT", f"数据集 {di} 的 R 区间 "
                       f"[{d['rmin']},{d['rmax']}] 内仅 {m.sum()} 点")
        s = float(np.sqrt(np.mean(np.abs(ds_data_r[di][m]) ** 2)))
        ds_mask.append(m)
        ds_scale.append(s if np.isfinite(s) and s > 0 else 1.0)

    def model_chi_r(di):
        d = datasets[di]
        model_k = exafs_model_k(ds_paths[di], ds_ft_k[di])
        model_ft = athena.ft_to_r_space(
            ds_ft_k[di], model_k, None,
            float(d["kmin"]), float(d["kmax"]),
            int(d["kweight"]), float(getattr(d["ft"], "dk", 1.0)))
        return (np.interp(ds_r_grid[di], model_ft.r, model_ft.chi_r.real)
                + 1j * np.interp(ds_r_grid[di], model_ft.r,
                                 model_ft.chi_r.imag))

    def residual(x):
        apply_params(x)
        parts = []
        for di in range(n_ds):
            model = model_chi_r(di)
            data = ds_data_r[di]
            m = ds_mask[di]
            s = ds_scale[di]
            parts.append((model.real - data.real)[m] / s)
            parts.append((model.imag - data.imag)[m] / s)
        return np.concatenate(parts)

    solution = least_squares(residual, x0_a, bounds=(lower_a, upper_a))

    dof = max(len(solution.fun) - n_var, 1)
    covariance = (np.linalg.pinv(solution.jac.T @ solution.jac)
                  * (2 * solution.cost / dof))
    errors = np.sqrt(np.abs(np.diag(covariance)))
    denom = np.outer(errors, errors)
    correlation = np.divide(
        covariance, denom, out=np.zeros_like(covariance), where=denom > 0)
    np.fill_diagonal(correlation, 1.0)

    # ---- 最优参数写回并计算各数据集 R-factor ----
    apply_params(solution.x)
    ds_rf = []
    for di in range(n_ds):
        model = model_chi_r(di)
        ds_rf.append(float(r_factor(ds_data_r[di], model, ds_mask[di])))
    scaled_res = residual(solution.x)
    den = 0.0
    for di in range(n_ds):
        m = ds_mask[di]
        den += float(np.sum(np.abs(ds_data_r[di][m]) ** 2)
                     / ds_scale[di] ** 2)
    global_rf = float(np.sum(scaled_res ** 2) / den) if den > 0 else np.inf

    # ---- 误差传播 ----
    def _slot_err(key):
        idx = slot_index.get(key)
        return None if idx is None else float(errors[idx])

    def resolved_dr_err(key):
        seen = set()
        k = key
        while k in bound_deltar:
            if k in seen:
                raise XAFSConfigError(
                    "FIT", f"share_deltar 绑定成环: {k}")
            seen.add(k)
            k = bound_deltar[k]
        return _slot_err((k[0], k[1], "deltar"))

    def resolved_sigma2_err(key):
        seen = set()
        k = key
        while True:
            if k in einstein_of:
                gi = einstein_of[k]
                i_static = slot_index[("einstein", gi, "sigma2_static")]
                i_theta = slot_index[("einstein", gi, "theta_e")]
                T = float(datasets[k[0]]["temperature_K"])
                mu = float(einstein_groups[gi]["reduced_mass_amu"])
                theta = float(solution.x[i_theta])
                xval = theta / (2.0 * T)
                coth = 1.0 / np.tanh(xval)
                csch2 = (1.0 / np.sinh(xval) ** 2
                         if xval < 350.0 else 0.0)
                a = HBAR2_OVER_2KB / mu
                d_theta = -a / theta ** 2 * (coth + xval * csch2)
                var = (errors[i_static] ** 2
                       + (d_theta * errors[i_theta]) ** 2
                       + 2.0 * d_theta * covariance[i_static, i_theta])
                return float(np.sqrt(abs(var)))
            if k in debye_of:
                gi = debye_of[k]
                model = debye_groups[gi]
                i_static = slot_index[("debye", gi, "sigma2_static")]
                i_theta = slot_index[("debye", gi, "theta_d")]
                T = float(datasets[k[0]]["temperature_K"])
                reff = float(ds_paths[k[0]][k[1]].reff)
                m_abs = float(model["mass_absorber"])
                m_sc = float(model["mass_scatterer"])
                kd = model.get("k_debye")
                theta = float(solution.x[i_theta])
                static = float(solution.x[i_static])
                f0 = correlated_debye_sigma2(
                    T, theta, reff, m_abs, m_sc, k_debye=kd,
                    sigma2_static=static)
                # 规范要求的有限差分步长：σ²_static 1e-4 Å²、Θ_D 0.5 K。
                d_s = 1e-4
                f_s = correlated_debye_sigma2(
                    T, theta, reff, m_abs, m_sc, k_debye=kd,
                    sigma2_static=static + d_s)
                d_static = (f_s - f0) / d_s
                d_t = 0.5
                f_t = correlated_debye_sigma2(
                    T, theta + d_t, reff, m_abs, m_sc, k_debye=kd,
                    sigma2_static=static)
                d_theta = (f_t - f0) / d_t
                var = (d_static ** 2 * errors[i_static] ** 2
                       + d_theta ** 2 * errors[i_theta] ** 2
                       + 2.0 * d_static * d_theta
                       * covariance[i_static, i_theta])
                return float(np.sqrt(abs(var)))
            if k in bound_sigma2:
                if k in seen:
                    raise XAFSConfigError(
                        "FIT", f"share_sigma2 绑定成环: {k}")
                seen.add(k)
                k = bound_sigma2[k]
                continue
            return _slot_err((k[0], k[1], "sigma2"))

    def resolved_n_err(key):
        seen = set()
        k = key
        factor = 1.0
        while k in bound_n:
            if k in seen:
                raise XAFSConfigError(
                    "FIT", f"degen_ratio 绑定成环: {k}")
            seen.add(k)
            src, ratio = bound_n[k]
            k = src
            factor *= ratio
        e = _slot_err((k[0], k[1], "n"))
        return None if e is None else factor * e

    # ---- 全局变量名 ----
    _kind_names = {"deltar": "ΔR", "sigma2": "σ²", "n": "N"}
    parameter_names = []
    for slot in slots:
        if slot[0] == "global":
            parameter_names.append("ΔE0")
        elif slot[0] == "e0grp":
            parameter_names.append(f"ΔE0:{e0_group_order[slot[1]]}")
        elif slot[0] == "einstein":
            parameter_names.append(
                f"einstein{slot[1]}:"
                + ("sigma2_static" if slot[2] == "sigma2_static"
                   else "theta_e"))
        elif slot[0] == "debye":
            parameter_names.append(
                f"debye{slot[1]}:"
                + ("sigma2_static" if slot[2] == "sigma2_static"
                   else "theta_d"))
        else:
            di, pi, kind = slot
            labels = datasets[di].get("path_labels")
            label = (labels[pi] if labels is not None
                     else Path(datasets[di]["datfiles"][pi]).name)
            name = datasets[di].get("name", f"ds{di}")
            parameter_names.append(f"{name}:{label}:{_kind_names[kind]}")

    def variables_for_dataset(di):
        """该数据集实际受影响的全局变量数（含被绑定源变量）。"""
        keys = set()

        def walk(key):
            if key in bound_deltar:
                walk(bound_deltar[key])
            else:
                keys.add((key[0], key[1], "deltar"))
            k = key
            seen = set()
            while k in bound_sigma2:
                if k in seen:
                    break
                seen.add(k)
                k = bound_sigma2[k]
            if k in einstein_of:
                gi = einstein_of[k]
                keys.add(("einstein", gi, "sigma2_static"))
                keys.add(("einstein", gi, "theta_e"))
            elif k in debye_of:
                gi = debye_of[k]
                keys.add(("debye", gi, "sigma2_static"))
                keys.add(("debye", gi, "theta_d"))
            else:
                keys.add((k[0], k[1], "sigma2"))
            k = key
            seen = set()
            while k in bound_n:
                if k in seen:
                    break
                seen.add(k)
                k = bound_n[k][0]
            if (k[0], k[1], "n") in slot_index:
                keys.add((k[0], k[1], "n"))

        for pi in range(len(ds_paths[di])):
            walk((di, pi))
        e0_key = e0_slot_key(e0_group_of_ds[di])
        if e0_key in slot_index:
            keys.add(e0_key)
        return len(keys)

    # ---- B12-1 AIC / BIC ----
    # 必须在 variables_for_dataset 定义之后计算：每个数据集的 k 是
    # 该数据集实际受影响的全局变量数（含被 shared_bonds 绑定的源变量），
    # 不能用全局 n_var 冒充。
    def _aic_bic(rss, n_obs, k):
        if n_obs <= 0 or rss <= 0:
            return None, None
        base = n_obs * np.log(rss / n_obs)
        return float(base + 2.0 * k), float(base + k * np.log(n_obs))

    ds_aic, ds_bic = [], []
    for di in range(n_ds):
        model_d = model_chi_r(di)
        m = ds_mask[di]
        rss_d = float(np.sum(np.abs(ds_data_r[di][m] - model_d[m]) ** 2))
        n_d = 2 * int(m.sum())
        k_d = int(variables_for_dataset(di))
        a_d, b_d = _aic_bic(rss_d, n_d, k_d)
        ds_aic.append(a_d)
        ds_bic.append(b_d)
    global_aic, global_bic = _aic_bic(
        float(np.sum(scaled_res ** 2)), int(scaled_res.size), int(n_var))

    # ---- 每个数据集的 FitResult ----
    n_idp_total = 0.0
    fit_results = []
    for di, d in enumerate(datasets):
        paths = ds_paths[di]
        labels = d.get("path_labels")
        model_r = model_chi_r(di)
        model_k = exafs_model_k(paths, ds_ft_k[di])
        fit_paths = []
        for pi, path in enumerate(paths):
            key = (di, pi)
            dr = resolved_dr(solution.x, key)
            dr_err = resolved_dr_err(key)
            s2 = resolved_sigma2(solution.x, key)
            s2_err = resolved_sigma2_err(key)
            n_val = resolved_n(solution.x, key)
            n_err = resolved_n_err(key)
            fit_paths.append(PathFit(
                label=(labels[pi] if labels is not None
                       else Path(d["datfiles"][pi]).name),
                reff=float(path.reff),
                degeneracy=float(n_val),
                delta_r=float(dr),
                delta_r_err=float(dr_err if dr_err is not None else 0.0),
                sigma2=float(s2),
                sigma2_err=(None if s2_err is None else float(s2_err)),
                n_err=(None if n_err is None else float(n_err)),
                sigma2_at_bound=bool(s2 <= 5e-4),
                nleg=int(getattr(path, "nleg", 2))))
        n_idp = independent_points(
            float(d["kmin"]), float(d["kmax"]),
            float(d["rmin"]), float(d["rmax"]))
        n_idp_total += n_idp
        fit_results.append(FitResult(
            paths=fit_paths,
            s02=float(d["s02"]),
            delta_e0=float(solution.x[slot_index[e0_slot_key(
                e0_group_of_ds[di])]]),
            delta_e0_err=float(errors[slot_index[e0_slot_key(
                e0_group_of_ds[di])]]),
            r_factor=float(ds_rf[di]),
            n_independent=float(n_idp),
            n_variables=int(variables_for_dataset(di)),
            kmin=float(d["kmin"]), kmax=float(d["kmax"]),
            rmin=float(d["rmin"]), rmax=float(d["rmax"]),
            parameter_names=list(parameter_names),
            correlation=correlation,
            model_r=np.asarray(ds_r_grid[di], dtype=float),
            model_chi_r=np.asarray(model_r, dtype=complex),
            model_k=np.asarray(model_k, dtype=float),
            model_k_grid=np.asarray(ds_ft_k[di], dtype=float),
            aic=ds_aic[di], bic=ds_bic[di]))

    # ---- 爱因斯坦结果 ----
    einstein_results = []
    for gi in range(len(einstein_groups)):
        i_static = slot_index[("einstein", gi, "sigma2_static")]
        i_theta = slot_index[("einstein", gi, "theta_e")]
        einstein_results.append({
            "theta_e": float(solution.x[i_theta]),
            "theta_e_err": float(errors[i_theta]),
            "sigma2_static": float(solution.x[i_static]),
            "sigma2_static_err": float(errors[i_static])})

    # ---- 关联 Debye 结果 ----
    debye_results = []
    for gi, model in enumerate(debye_groups):
        i_static = slot_index[("debye", gi, "sigma2_static")]
        i_theta = slot_index[("debye", gi, "theta_d")]
        debye_results.append({
            "theta_d": float(solution.x[i_theta]),
            "theta_d_err": float(errors[i_theta]),
            "sigma2_static": float(solution.x[i_static]),
            "sigma2_static_err": float(errors[i_static]),
            "paths": [(int(item[0]), int(item[1]))
                      for item in model["paths"]],
            "mass_absorber": float(model["mass_absorber"]),
            "mass_scatterer": float(model["mass_scatterer"]),
            "k_debye": model.get("k_debye"),
        })

    return {
        "datasets": fit_results,
        "n_independent_total": float(n_idp_total),
        "n_variables_total": int(n_var),
        "global_r_factor": global_rf,
        "einstein_results": einstein_results,
        "debye_results": debye_results,
        "parameter_names": parameter_names,
        "correlation": correlation,
        "aic": global_aic,
        "bic": global_bic,
        "global_aic": global_aic,
        "global_bic": global_bic,
        "e0_groups": [e0_group_order[gi] for gi in range(n_e0_groups)],
        "e0_group_of_dataset": list(e0_group_of_ds),
    }


# ================================================================
# 6. 自测
# ================================================================

if __name__ == "__main__":
    import sys
    import tempfile
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    # ---- CIF → sites 与内置 bcc 生成器交叉验证 ----
    bcc_cif = "\n".join([
        "data_fe_bcc",
        "_cell_length_a 2.866",
        "_cell_length_b 2.866",
        "_cell_length_c 2.866",
        "_cell_angle_alpha 90",
        "_cell_angle_beta 90",
        "_cell_angle_gamma 90",
        "_symmetry_space_group_name_H-M 'I m -3 m'",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "Fe1 Fe 0 0 0",
    ])
    with tempfile.TemporaryDirectory() as td:
        cif = Path(td) / "fe_bcc.cif"
        cif.write_text(bcc_cif, encoding="utf-8")
        sites_cif = read_cif_sites(cif, absorber="Fe", rmax=6.0)
    d_cif = sorted(round(float(np.sqrt(s[0]**2 + s[1]**2 + s[2]**2)), 4)
                   for s in sites_cif)
    d_ref = sorted(round(float(np.sqrt(s[0]**2 + s[1]**2 + s[2]**2)), 4)
                   for s in bcc_lattice_sites(2.866, 6.0))
    assert d_cif == d_ref, f"CIF 展开与内置 bcc 不一致: {d_cif[:8]} vs {d_ref[:8]}"
    print(f"CIF→sites 校验: {len(sites_cif)} 个近邻，与 bcc_lattice_sites 一致")

    workdir = Path("./_artemis_selftest")
    print("生成 Fe 箔（bcc）的 feff.inp ...")
    inp = build_feff_input(bcc_lattice_sites(2.866, 6.0), absorber="Fe", rmax=6.0,
                           title="Fe foil bcc")
    files = run_feff6l(inp, workdir)
    print(f"  路径数: {len(files)}")
    for row in describe_paths(files, 2):
        print(f"  {row['file']}: Reff={row['reff']:.3f} Å  "
              f"简并度={row['degeneracy']:.0f}  腿数={row['legs']}")
