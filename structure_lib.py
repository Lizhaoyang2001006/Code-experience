# -*- coding: utf-8 -*-
"""
structure_lib.py — motif 注册表（结构假设库）
================================================================
每条 motif 是"一种可能的局域结构"，双源统一：

    kind=crystal  晶体骨架：{cell, spacegroup, asym_sites}
                  —— 即"空间群 + 晶胞 + 特征点(Wyckoff 坐标)"，
                     gemmi 展开成全晶胞（与 read_cif_sites 同一引擎）
    kind=cluster  团簇/分子：绝对直角坐标 [(x,y,z,el), ...]
                  —— SAC、双核、DFT 模型等无空间群对象

注册（一次性）流程：
    骨架/坐标 → 展开 → 以吸收原子为原点的近邻坐标表
    → 几何壳层指纹 shell_expect（自动计算，预筛用）
    → motif JSON 落盘 _structure_lib/<absorber>/<id>.json
    → 查重（指纹比对，容差 TOL_DUP）

库是人工精管的——在线检索（COD API 等）只是往注册表里
添条目的手段，失败不影响库内既有条目（去除网页抓取脆弱性）。
motif JSON 带 version 字段，meta.json 记录库版本以便复现。

自测: python structure_lib.py  → 加载全部 motif 并打印指纹。
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import artemis
from xafs_core import XAFSConfigError, XAFSDataError

LIB_ROOT = Path(__file__).resolve().parent / "_structure_lib"
FEFF_CACHE = Path(__file__).resolve().parent / "_feff_cache"
# 理论小波标尺缓存（reference_element_k_peaks）：与 motif FEFF 缓存分开，
# 避免不同用途的理论路径互相覆盖。
CACHE_ROOT = LIB_ROOT / "_wt_ref_cache"
# B13-1 动态催化 motif（SAC/DAC）缓存：按结构参数哈希分目录，可复现、不进主库。
DYNAMIC_CACHE = LIB_ROOT / "_dynamic_cache"
# B13-2 异核共占位分支的单散射 FEFF 缓存。
HETERO_CACHE = LIB_ROOT / "_hetero_cache"

SHELL_TOL = 0.08        # 几何壳层合并容差 (Å)：小于此间距并为一层
TOL_DUP = 0.05          # 查重容差 (Å)：指纹层位差小于此视为重复
LIB_VERSION = 1         # 注册表格式版本，变更时 +1 并写迁移说明


@dataclass
class Motif:
    """一条结构假设。sites/shell_expect 在注册时生成并随 JSON 持久化。"""
    id: str
    absorber: str
    kind: str                          # "crystal" | "cluster" | "mixture"
    nuclearity: str                    # m1/m2/m3_6/cluster/np/crystal
    source: dict = field(default_factory=dict)
    cell: dict | None = None           # crystal: a,b,c,alpha,beta,gamma
    spacegroup: str | None = None
    asym_sites: list | None = None     # crystal: [{el,x,y,z,occ?,wyckoff?}]
    absorber_fract: list | None = None # 多位点吸收原子选位（分数坐标）
    cif_path: str | None = None        # 可选：直接引用外部 CIF（含 _atom_site_occupancy）
    cache_dir: str | None = None       # B13-1：动态 motif 的 FEFF 缓存根目录
                                       # （None → 全局 FEFF_CACHE，历史行为）
    sites: list | None = None          # 展开后近邻坐标 [[x,y,z,el],...]
    shell_expect: list | None = None   # 几何指纹 [{el,r,n},...]
    components: list | None = None     # kind=mixture: [{id, weight}, ...]
    size_model: dict | None = None     # 纳米颗粒模型的壳层 N 缩放先验
    notes: str = ""
    version: int = LIB_VERSION

    # ---------- 序列化 ----------
    def to_json(self) -> dict:
        return {k: v for k, v in self.__dict__.items() if v is not None}

    @staticmethod
    def from_json(d: dict) -> "Motif":
        known = {f for f in Motif.__dataclass_fields__}
        return Motif(**{k: v for k, v in d.items() if k in known})


def geometric_shells(sites, absorber_tol: float = SHELL_TOL,
                     occs=None) -> list:
    """
    由近邻坐标表计算几何壳层指纹：[{el, r, n}]，按 r 升序。

    距离差 < absorber_tol 的同元素近邻并为一层（n=计数，r=均值）。
    occs 提供时另计 n_occ（占据率加权和）——固溶体/无序体系的
    N 先验；n 仍是几何计数，两者并存不互相覆盖。
    这是**免 FEFF 的预筛指纹**——只含几何信息，不含散射相移，
    因此用于"排除"而非"确认"。
    """
    rows = []
    for i, (x, y, z, el) in enumerate(sites):
        occ = float(occs[i]) if occs is not None else 1.0
        rows.append((el, float(np.linalg.norm([x, y, z])), occ))
    rows.sort(key=lambda t: (t[1], t[0]))
    shells = []
    for el, r, occ in rows:
        if shells and shells[-1]["el"] == el \
                and abs(r - shells[-1]["r"]) < absorber_tol:
            sh = shells[-1]
            sh["r"] = (sh["r"] * sh["n"] + r) / (sh["n"] + 1)
            sh["n"] += 1
            sh["n_occ"] += occ
        else:
            shells.append({"el": el, "r": r, "n": 1, "n_occ": occ})
    if occs is None:
        for sh in shells:
            del sh["n_occ"]
    return shells


def _cif_path(motif: Motif) -> Path:
    """由晶体骨架现场拼最小 CIF 并落盘（留档可复查），返回路径。

    若 motif.cif_path 指向外部 CIF（例如含部分占据/共占据的实验 CIF），
    直接使用它，不再由 asym_sites 重建。
    """
    if motif.cif_path:
        external = Path(motif.cif_path)
        if external.exists():
            return external
    c = motif.cell or {}
    sg = motif.spacegroup
    if not c or not sg or not motif.asym_sites:
        raise XAFSDataError("motif", f"{motif.id}: 晶体骨架缺 cell/spacegroup/asym_sites")
    lines = [
        f"data_{motif.id}",
        f"_cell_length_a {c['a']}",
        f"_cell_length_b {c.get('b', c['a'])}",
        f"_cell_length_c {c.get('c', c['a'])}",
        f"_cell_angle_alpha {c.get('alpha', 90)}",
        f"_cell_angle_beta {c.get('beta', 90)}",
        f"_cell_angle_gamma {c.get('gamma', 90)}",
        f"_symmetry_space_group_name_H-M '{sg}'",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_occupancy",
    ]
    for i, s in enumerate(motif.asym_sites, 1):
        lines.append(
            f"{s['el']}{i} {s['el']} {s['x']} {s['y']} {s['z']} "
            f"{s.get('occ', 1.0)}")
    cif_dir = LIB_ROOT / motif.absorber / "_cif"
    cif_dir.mkdir(parents=True, exist_ok=True)
    cif_path = cif_dir / f"{motif.id}.cif"
    cif_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return cif_path


def expand_crystal(motif: Motif, rmax: float = 6.0) -> list:
    """
    晶体 motif → 近邻坐标表。

    做法：由骨架现场拼一个最小 CIF 文本，复用 artemis.read_cif_sites
    （gemmi 展开）——**对称性逻辑全库只有一份**，不在这里重造。
    """
    return artemis.read_cif_sites(_cif_path(motif), motif.absorber, rmax=rmax,
                                  absorber_fract=motif.absorber_fract)


def expand_full(motif: Motif, rmax: float = 6.0):
    """expand 的占据率保留版：返回 (sites, occs|None)。cluster 无占据信息。"""
    if motif.kind == "crystal":
        return artemis.read_cif_sites(_cif_path(motif), motif.absorber,
                                      rmax=rmax,
                                      absorber_fract=motif.absorber_fract,
                                      with_occupancy=True)
    return expand(motif, rmax), None


def expand(motif: Motif, rmax: float = 6.0) -> list:
    """双源统一入口：任何 motif → 近邻坐标表。"""
    if motif.kind == "mixture":
        raise XAFSDataError(
            "motif", f"{motif.id}: mixture 没有独立坐标——用 "
                     f"components 引用的单一 motif 展开（χ 为加权求和）")
    if motif.kind == "crystal":
        return expand_crystal(motif, rmax)
    if motif.kind == "cluster":
        if not motif.sites:
            raise XAFSDataError("motif", f"{motif.id}: cluster 缺 sites")
        # 约定：吸收原子在原点（不入列），sites 只含近邻
        return [tuple(s) for s in motif.sites
                if float(np.linalg.norm(s[:3])) > 0.1]
    raise XAFSDataError("motif", f"{motif.id}: 未知 kind={motif.kind!r}")


def register(motif: Motif, rmax: float = 6.0, overwrite: bool = False) -> Path:
    """
    展开 motif、计算指纹、查重、落盘。返回 JSON 路径。

    查重：同库内已有 motif 的 shell_expect 与本次计算的层位差
    全部 < TOL_DUP 且元素序列一致 → 判重复，拒绝注册（防止同一结构
    以不同来源/设置重复入库，后续识别时无法区分该算谁）。

    kind=mixture（多位点吸收）：不展开坐标、不算几何指纹，改为校验
    components——引用的单一 motif 必须同库已存在（因此先注册位点，再
    注册混合条目）、权重为正且总和为 1、且不得嵌套 mixture。
    """
    lib_dir = LIB_ROOT / motif.absorber
    out = lib_dir / f"{motif.id}.json"
    if motif.kind == "mixture":
        _validate_components(motif, load_library(motif.absorber))
        motif.sites = None
        motif.shell_expect = None
    else:
        motif.sites, occs = expand_full(motif, rmax)
        motif.sites = [list(map(lambda v: round(v, 6) if isinstance(v, float) else v, s))
                       for s in motif.sites]
        motif.shell_expect = [
            {k: (round(v, 4) if isinstance(v, float) else v)
             for k, v in s.items() if k != "n_occ" or occs is not None}
            for s in geometric_shells(motif.sites, occs=occs)]
    if out.exists() and not overwrite:
        raise XAFSConfigError("motif", f"{motif.id} 已存在（overwrite=False）")
    for other in load_library(motif.absorber):
        if other.id != motif.id and _same_fingerprint(other, motif):
            raise XAFSConfigError(
                "motif", f"{motif.id} 与已有 {other.id} 指纹重复，拒绝注册")
    lib_dir.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(motif.to_json(), ensure_ascii=False, indent=2),
                   encoding="utf-8")
    return out


def _validate_components(motif: Motif, library: list) -> None:
    """mixture 条目的 components 校验（register 调用）。"""
    comps = motif.components
    if not comps:
        raise XAFSConfigError("motif", f"{motif.id}: mixture 缺 components")
    ids, total = [], 0.0
    for comp in comps:
        try:
            cid, w = comp["id"], float(comp["weight"])
        except (KeyError, TypeError, ValueError):
            raise XAFSConfigError(
                "motif", f"{motif.id}: components 项须为 {{\"id\", \"weight\"}}，"
                         f"得到 {comp!r}") from None
        if w <= 0:
            raise XAFSConfigError("motif", f"{motif.id}: 权重须 >0，得到 {w}")
        total += w
        ids.append(cid)
        found = [m for m in library if m.id == cid]
        if not found:
            raise XAFSConfigError(
                "motif", f"{motif.id}: 引用的 {cid} 不在 _structure_lib/"
                         f"{motif.absorber}——先注册各位点，再注册混合条目")
        if found[0].components:
            raise XAFSConfigError("motif", f"{motif.id}: 不允许嵌套 mixture")
    if abs(total - 1.0) > 1e-6:
        raise XAFSConfigError(
            "motif", f"{motif.id}: 权重总和={total}，应为 1")
    if len(set(ids)) != len(ids):
        raise XAFSConfigError("motif", f"{motif.id}: components 引用重复")


def _same_fingerprint(a: Motif, b: Motif) -> bool:
    sa, sb = a.shell_expect or [], b.shell_expect or []
    if len(sa) != len(sb):
        return False
    return all(x["el"] == y["el"] and abs(x["r"] - y["r"]) < TOL_DUP
               and x["n"] == y["n"] for x, y in zip(sa, sb))


def load_library(absorber: str) -> list:
    """加载某元素的全部 motif（按 id 排序）。"""
    lib_dir = LIB_ROOT / absorber
    if not lib_dir.is_dir():
        return []
    return [Motif.from_json(json.loads(p.read_text(encoding="utf-8")))
            for p in sorted(lib_dir.glob("*.json"))]


def ensure_feff(motif: Motif, rmax: float = 6.0, edge: str = "K") -> Path:
    """
    保证 motif 的 FEFF 缓存存在，返回缓存目录。

    缓存键 = motif id + 近邻坐标哈希前 8 位（非 K 边再加边名）——
    坐标/吸收边变了缓存自动失效，不会出现"库改了但 FEFF 是旧的"
    或"L3 motif 复用了 K 边路径"这类静默错误。

    mixture 条目没有自己的 FEFF：调用方应分别对 components 的各位点
    取理论路径（structure_id.theory_chi / xafs_fit._ss_rows 已处理）。
    """
    import hashlib
    if motif.kind == "mixture":
        raise XAFSDataError(
            "motif", f"{motif.id}: mixture 自己不跑 FEFF——"
                     f"请对 components 引用的各位点 ensure_feff")
    key = hashlib.md5(
        json.dumps(motif.sites, sort_keys=True).encode()).hexdigest()[:8]
    suffix = "" if edge == "K" else f"_{edge.lower()}"
    # 动态 motif（B13-1）把 FEFF 输出固定在自己的哈希目录，避免与主库
    # 同名缓存互串、也保证临时构建的位点在进程重启后仍可复现。
    cache_root = Path(motif.cache_dir) if motif.cache_dir else FEFF_CACHE
    cache = cache_root / f"lib_{motif.id}_{key}{suffix}"
    if not list(cache.glob("feff*.dat")):
        inp = artemis.build_feff_input(motif.sites, absorber=motif.absorber,
                                       rmax=rmax, title=motif.id, edge=edge)
        artemis.run_feff6l(inp, cache)
    return cache


def import_cif(cif_path, motif_id: str, absorber: str,
               nuclearity: str = "crystal", source: dict | None = None,
               absorber_fract=None) -> Motif:
    """
    CIF 导入器：gemmi 抽取骨架（cell + 空间群 + 非对称单元）→ motif。

    这是"网上抓到 CIF 后入库"的通道：只保留骨架，其余丢弃，
    等价于"只取特征点"，但接受任意来源的完整小 CIF。
    """
    import gemmi
    doc = gemmi.cif.read(str(cif_path))
    block = doc.sole_block() if len(doc) == 1 else doc[0]
    st = gemmi.make_small_structure_from_block(block)
    cell = st.cell
    asym = []
    for site in st.sites:
        el = site.type_symbol.strip()
        import re
        m = re.match(r"[A-Za-z]+", el)
        if not m:
            continue
        asym.append({"el": m.group(0).capitalize(),
                     "x": round(float(site.fract.x), 5),
                     "y": round(float(site.fract.y), 5),
                     "z": round(float(site.fract.z), 5),
                     "occ": round(float(site.occ), 4),
                     "label": site.label})
    motif = Motif(
        id=motif_id, absorber=absorber, kind="crystal", nuclearity=nuclearity,
        source=source or {"type": "cif", "file": str(cif_path)},
        cell={"a": round(cell.a, 5), "b": round(cell.b, 5),
              "c": round(cell.c, 5), "alpha": round(cell.alpha, 3),
              "beta": round(cell.beta, 3), "gamma": round(cell.gamma, 3)},
        spacegroup=st.spacegroup_hm,
        asym_sites=asym, absorber_fract=absorber_fract)
    return motif


# ================================================================
# 替位编辑（掺杂体系候选名册的机械生成）
# ================================================================

# 常见离子对的参考键长 (Å)：Shannon 半径量级的手管表。
# 查不到就必须显式给 relax_to——不猜。
BOND_REF = {
    ("Fe", "O", "tet"): 1.86, ("Fe", "O", "oct"): 2.00,
    ("Co", "O", "tet"): 1.83, ("Co", "O", "oct"): 1.93,
    ("Ni", "O", "oct"): 2.05, ("Cu", "O", "oct"): 1.96,
    ("Si", "O", "tet"): 1.61, ("Al", "O", "tet"): 1.75,
    ("Fe", "N", None): 1.93, ("Co", "N", None): 1.90,
    ("Fe", "C", None): 1.99,
}


def substitute(host: Motif, host_el: str, dopant: str, site_fract,
               motif_id: str | None = None, relax_to: float | None = None,
               coord: str | None = None, rmax: float = 6.0) -> Motif:
    """
    把宿主晶体指定位点替换为掺杂元素，返回 cluster motif（掺杂位在原点）。

    局域弛豫：第一壳（d ≤ 最近邻×1.25）整体径向缩放 f = relax_to/d_min，
    **保持层内劈裂**（如八面体 Jahn-Teller 的长短键不被抹平）；
    远程几何保持宿主——"局域畸变、长程受限"的标准近似，
    深度弛豫属 DFT 范畴（算完用坐标重注册即可）。

    relax_to=None 时查 BOND_REF[(dopant, 配体, coord)]，查不到报错不猜。
    """
    host_el = host_el.capitalize()
    dopant = dopant.capitalize()
    sites = artemis.read_cif_sites(_cif_path(host), absorber=host_el,
                                   absorber_fract=site_fract, rmax=rmax)
    if not sites:
        raise XAFSDataError("substitute",
                            f"{host.id} 的 {host_el}@{site_fract} 展开为空")
    dists = [float(np.linalg.norm(s[:3])) for s in sites]
    d_min = min(dists)
    cutoff = d_min * 1.25
    ligands = {s[3] for s, d in zip(sites, dists) if d <= cutoff}
    if relax_to is None:
        if len(ligands) != 1:
            raise XAFSConfigError(
                "substitute", f"第一壳含多种配体 {sorted(ligands)}，"
                              f"须显式给 relax_to")
        lig = next(iter(ligands))
        relax_to = BOND_REF.get((dopant, lig, coord)) \
            or BOND_REF.get((dopant, lig, None))
        if relax_to is None:
            raise XAFSConfigError(
                "substitute", f"BOND_REF 无 ({dopant},{lig},{coord})，"
                              f"须显式给 relax_to")
    f = relax_to / d_min
    new_sites = [[x * f, y * f, z * f, el] if d <= cutoff else [x, y, z, el]
                 for (x, y, z, el), d in zip(sites, dists)]
    mid = motif_id or f"{host.id}__{dopant.lower()}_site"
    return Motif(
        id=mid, absorber=dopant, kind="cluster", nuclearity="m1",
        source={"type": "substitution", "host": host.id, "host_el": host_el,
                "site": [float(v) for v in site_fract],
                "relax_to": relax_to},
        sites=new_sites,
        notes=f"{host.id} 的 {host_el} 位替位为 {dopant}；第一壳缩放 "
              f"{f:.3f} 至 {relax_to:.2f} Å（层内劈裂保留），远程保持宿主。")


def substitute_series(host: Motif, host_el: str, dopant: str,
                      relax_to: float | None = None, coord: str | None = None,
                      rmax: float = 6.0) -> list:
    """
    对宿主中 host_el 的**每个结晶学不等价位点**逐一替位，返回候选名册。
    等价位点在展开后指纹相同，注册阶段查重会自动去重（register_series）。
    """
    targets = [s for s in (host.asym_sites or [])
               if s["el"] == host_el.capitalize()]
    if not targets:
        raise XAFSDataError("substitute",
                            f"{host.id} 的 asym_sites 里没有 {host_el}")
    out = []
    for s in targets:
        label = str(s.get("wyckoff") or s.get("label")
                    or len(out) + 1).replace(" ", "")
        out.append(substitute(host, host_el, dopant,
                              [s["x"], s["y"], s["z"]],
                              motif_id=f"{host.id}__{dopant.lower()}_{label}",
                              relax_to=relax_to, coord=coord, rmax=rmax))
    return out


def register_series(motifs, rmax: float = 6.0):
    """批量注册；指纹重复（等价位点）跳过并回报，不中断。"""
    made, skipped = [], []
    for m in motifs:
        try:
            made.append(register(m, rmax=rmax))
        except XAFSConfigError as exc:
            skipped.append(f"{m.id}: {exc.message}")
    return made, skipped


# ================================================================
# 纳米颗粒几何模型（表面-核心平均的 N 先验）
# ================================================================

def nanoparticle_model(host: Motif, size_nm: float,
                       rmax_shell: float = 4.0) -> dict:
    """
    球形颗粒的几何配位统计：对球内每个吸收原子数其壳层近邻，
    得各壳层平均 N 相对体相的缩放因子（拟合 N 的先验区间）。

    **只改 N 先验，不改散射物理**——FEFF 路径仍用体相 motif 的；
    纯几何硬球截断，不含表面重构/弛豫；<1 nm 统计噪声大，结果只作先验。
    """
    if not host.shell_expect:
        raise XAFSDataError("nanoparticle",
                            f"{host.id} 缺 shell_expect，先注册")
    R = size_nm * 10.0 / 2.0
    sites = expand(host, rmax=R + rmax_shell + 1.0)
    origin = np.zeros(3)
    abs_el = host.absorber
    atoms = [(np.array(s[:3]), s[3]) for s in sites
             if np.linalg.norm(s[:3]) <= R]
    centers = [origin] + [p for p, el in atoms if el == abs_el]
    atoms_plus = atoms + [(origin, abs_el)]      # 中心原子也是别人的近邻
    bulk = [sh for sh in host.shell_expect if sh["r"] <= rmax_shell]
    counts = np.zeros(len(bulk))
    for c in centers:
        for p, el in atoms_plus:
            d = float(np.linalg.norm(p - c))
            if d < 0.1 or d > rmax_shell:
                continue
            for j, sh in enumerate(bulk):
                if sh["el"] == el and abs(d - sh["r"]) < 0.15:
                    counts[j] += 1.0
                    break
    shells = []
    for sh, cnt in zip(bulk, counts):
        n_avg = float(cnt) / max(len(centers), 1)
        shells.append({"el": sh["el"], "r": sh["r"], "n_bulk": sh["n"],
                       "n_avg": round(n_avg, 3),
                       "scale": round(n_avg / sh["n"], 4) if sh["n"] else None})
    return {"size_nm": size_nm, "radius_A": R, "n_absorber": len(centers),
            "shells": shells}


def nanoparticle_motif(host: Motif, size_nm: float,
                       motif_id: str | None = None) -> Motif:
    """把尺寸模型挂进 motif 副本（size_model 字段）。调用方负责注册。"""
    model = nanoparticle_model(host, size_nm)
    data = host.to_json()
    data.update({"id": motif_id or f"{host.id}_np{size_nm:g}nm",
                 "nuclearity": "np", "size_model": model,
                 "notes": (host.notes or "")
                          + f"｜{size_nm} nm 球形颗粒几何模型，"
                            f"壳层 N 先验见 size_model"})
    return Motif.from_json(data)


# ================================================================
# B13-1 动态催化 motif（SAC / DAC）：构建、缓存与字符串规格解析
# ================================================================

# 常见配体元素：DAC 字符串里 "Co@2.50" 与 "N4@2.00" 靠它区分金属伙伴/配体。
_LIGAND_ELEMENTS = {"N", "O", "C", "S", "P", "F", "Cl", "Br", "I", "B", "Se"}
_CATALYTIC_DICT_KEYS = {
    "sac": {"type", "ligand", "n_coord", "bond_m_ligand", "axial",
            "axial_bond", "carbon_shell", "ring_radius"},
    "dac": {"type", "partner", "bond_length", "ligand", "n_coord",
            "ligand_bond"},
    "alloy": {"type", "partner", "x_partner", "r_aa", "r_ab",
              "n_total", "lattice"},
}

# 第一配位壳方向（单位向量）。fcc/hcp 用立方八面体 12 方向，
# bcc 用立方体 8 个顶点方向——与体相第一壳配位数一致。
_FCC12_DIRS = [
    (1.0, 1.0, 0.0), (1.0, -1.0, 0.0), (-1.0, 1.0, 0.0), (-1.0, -1.0, 0.0),
    (1.0, 0.0, 1.0), (1.0, 0.0, -1.0), (-1.0, 0.0, 1.0), (-1.0, 0.0, -1.0),
    (0.0, 1.0, 1.0), (0.0, 1.0, -1.0), (0.0, -1.0, 1.0), (0.0, -1.0, -1.0),
]
_BCC8_DIRS = [
    (1.0, 1.0, 1.0), (1.0, 1.0, -1.0), (1.0, -1.0, 1.0), (1.0, -1.0, -1.0),
    (-1.0, 1.0, 1.0), (-1.0, 1.0, -1.0), (-1.0, -1.0, 1.0),
    (-1.0, -1.0, -1.0),
]


def _dynamic_cache_dir(payload: dict) -> Path:
    """按结构参数哈希给出动态 motif 的固定缓存目录（可复现、跨进程稳定）。"""
    import hashlib
    key = hashlib.md5(
        json.dumps(payload, sort_keys=True, default=str).encode()
    ).hexdigest()[:12]
    out = DYNAMIC_CACHE / f"{payload.get('type', 'motif')}_{key}"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _finalize_dynamic_motif(motif: Motif, payload: dict) -> Motif:
    """计算几何指纹、落盘 motif.json、绑定 FEFF 缓存目录。"""
    motif.sites = [[round(float(v), 6) if isinstance(v, float) else v
                    for v in s] for s in motif.sites]
    motif.shell_expect = geometric_shells(motif.sites)
    cache_dir = _dynamic_cache_dir(payload)
    motif.cache_dir = str(cache_dir)
    (cache_dir / "motif.json").write_text(
        json.dumps(motif.to_json(), ensure_ascii=False, indent=2),
        encoding="utf-8")
    return motif


def build_sac_motif(
    absorber: str,
    ligand: str = "N",
    n_coord: int = 4,
    bond_m_ligand: float = 2.00,
    axial: str | None = None,
    axial_bond: float = 2.10,
    carbon_shell: bool = True,
    ring_radius: float = 2.90,
) -> Motif:
    """B13-1 构建单原子催化剂（SAC）动态 motif。

    几何由 ``artemis.sac_cluster_sites`` 生成（吸收原子在原点），
    FEFF 缓存落在 ``_structure_lib/_dynamic_cache/<hash>/``——不写进
    ``_structure_lib/<absorber>/*.json``，因此不会改变默认结构识别排序。
    """
    sym = str(absorber).capitalize()
    if sym not in artemis.ATOMIC_NUMBER:
        raise XAFSConfigError("SAC", f"未知吸收元素 {absorber!r}")
    lig = str(ligand).capitalize()
    sites = artemis.sac_cluster_sites(
        bond_m_ligand=float(bond_m_ligand), ligand=lig,
        n_coord=int(n_coord), axial=axial, axial_bond=float(axial_bond),
        carbon_shell=bool(carbon_shell), ring_radius=float(ring_radius))
    payload = {"type": "sac", "absorber": sym, "ligand": lig,
               "n_coord": int(n_coord), "bond_m_ligand": float(bond_m_ligand),
               "axial": (None if not axial else str(axial).capitalize()),
               "axial_bond": float(axial_bond),
               "carbon_shell": bool(carbon_shell),
               "ring_radius": float(ring_radius)}
    motif_id = f"sac_{sym}_{lig}{int(n_coord)}_{float(bond_m_ligand):.2f}"
    if axial:
        motif_id += f"_ax{payload['axial']}{float(axial_bond):.2f}"
    motif_id += (f"_C{float(ring_radius):.2f}" if carbon_shell else "_noC")
    motif = Motif(
        id=motif_id, absorber=sym, kind="cluster", nuclearity="mononuclear",
        source={"type": "constructed", "module": "B13-1",
                "ref": "SAC M-N4 平面 + 可选轴向配体 / 碳壳层"},
        sites=[list(s) for s in sites],
        notes=(f"动态 SAC 假设：{sym}-{lig}{int(n_coord)}"
               f"@{float(bond_m_ligand):.2f} Å"
               + (f"，轴向 {payload['axial']}@{float(axial_bond):.2f} Å"
                  if axial else "，无轴向配体")
               + ("，含碳载体壳层" if carbon_shell else "，无碳壳层")
               + "；几何为构造值，定量前建议用 DFT 精化。"))
    return _finalize_dynamic_motif(motif, payload)


def build_dac_motif(
    absorber: str,
    partner: str,
    bond_length: float = 2.50,
    ligand: str = "N",
    n_coord: int = 4,
    ligand_bond: float = 2.00,
) -> Motif:
    """B13-1 构建双原子催化剂（DAC）动态 motif（含 M₁–M₂ 金属键路径）。"""
    sym = str(absorber).capitalize()
    par = str(partner).capitalize()
    for el in (sym, par):
        if el not in artemis.ATOMIC_NUMBER:
            raise XAFSConfigError("DAC", f"未知元素 {el!r}")
    if sym == par:
        raise XAFSConfigError("DAC", "absorber 与 partner 不能相同（同核双原子请用 SAC 系列）")
    lig = str(ligand).capitalize()
    sites = artemis.diatomic_m1m2_sites(
        m1=sym, m2=par, bond_length=float(bond_length), ligand=lig,
        n_coord=int(n_coord), ligand_bond=float(ligand_bond))
    payload = {"type": "dac", "absorber": sym, "partner": par,
               "bond_length": float(bond_length), "ligand": lig,
               "n_coord": int(n_coord), "ligand_bond": float(ligand_bond)}
    motif_id = (f"dac_{sym}_{par}_{float(bond_length):.2f}"
                f"_{lig}{int(n_coord)}_{float(ligand_bond):.2f}")
    motif = Motif(
        id=motif_id, absorber=sym, kind="cluster", nuclearity="binuclear",
        source={"type": "constructed", "module": "B13-1",
                "ref": "DAC M1-M2 直接金属键 + 双方平面配体"},
        sites=[list(s) for s in sites],
        notes=(f"动态 DAC 假设：{sym}–{par} {float(bond_length):.2f} Å"
               f" + {sym}-{lig}{int(n_coord)}@{float(ligand_bond):.2f} Å；"
               "M–M 路径是否存在必须用 Hamilton F-test 判定，不能只看 Rf。"))
    return _finalize_dynamic_motif(motif, payload)


def build_alloy_motif(
    absorber: str,
    partner: str,
    x_partner: float = 0.5,
    r_aa: float | None = None,
    r_ab: float | None = None,
    n_total: int = 12,
    lattice: str = "fcc",
) -> Motif:
    """B17-2 构建双金属纳米合金动态 motif（免 CIF）。

    第一配位壳按 fcc 立方八面体 12 方向（bcc 用 8 顶点方向）放置：
    前 ``n_ab`` 个方向放 partner、其余放 absorber，形成同核 A–A 与
    异核 A–B 两条可拟合路径。缓存落 ``_dynamic_cache/alloy_<hash>``，
    不进入主结构库，默认识别排序零影响。
    """
    sym = str(absorber).capitalize()
    par = str(partner).capitalize()
    lat = str(lattice).strip().lower()
    if sym not in artemis.ATOMIC_NUMBER:
        raise XAFSConfigError("ALLOY", f"未知吸收元素 {absorber!r}")
    if par not in artemis.ATOMIC_NUMBER:
        raise XAFSConfigError("ALLOY", f"未知 partner 元素 {partner!r}")
    if sym == par:
        raise XAFSConfigError(
            "ALLOY", "absorber 与 partner 不能相同（同核请用晶体 motif）")
    if lat not in ("fcc", "hcp", "bcc"):
        raise XAFSConfigError(
            "ALLOY", f"lattice={lattice!r} 只支持 'fcc' / 'hcp' / 'bcc'")
    xb = float(x_partner)
    if not (np.isfinite(xb) and 0.0 < xb < 1.0):
        raise XAFSConfigError(
            "ALLOY", f"x_partner={x_partner!r} 必须在 (0, 1) 内")
    n_tot = int(n_total)
    if n_tot < 2:
        raise XAFSConfigError("ALLOY", f"n_total={n_total!r} 必须 ≥ 2")

    dirs = _BCC8_DIRS if lat == "bcc" else _FCC12_DIRS
    if n_tot > len(dirs):
        raise XAFSConfigError(
            "ALLOY", f"n_total={n_tot} 超过 {lat} 第一壳方向数 {len(dirs)}")
    dirs = [np.asarray(d, dtype=float) / np.linalg.norm(d) for d in dirs[:n_tot]]

    r_cov_a = artemis.COVALENT_RADII.get(sym)
    r_cov_b = artemis.COVALENT_RADII.get(par)
    if r_cov_a is None or r_cov_b is None:
        raise XAFSConfigError(
            "ALLOY", f"缺少共价半径：{sym}={r_cov_a!r}, {par}={r_cov_b!r}")
    raa = 2.0 * float(r_cov_a) if r_aa is None else float(r_aa)
    rab = (float(r_cov_a) + float(r_cov_b)
           if r_ab is None else float(r_ab))
    for name, val in (("r_aa", raa), ("r_ab", rab)):
        if not (np.isfinite(val) and val > 0.0):
            raise XAFSConfigError(
                "ALLOY", f"{name}={val!r} 必须是有限正数")
    # 两距离小于 0.01 Å 时视为同一壳层：用平均值强制精确相等。
    if abs(raa - rab) < 0.01:
        same = 0.5 * (raa + rab)
        raa = rab = same

    n_ab = max(1, min(n_tot - 1, int(round(n_tot * xb))))
    n_aa = n_tot - n_ab
    sites = []
    for i, u in enumerate(dirs):
        if i < n_ab:
            sites.append([float(u[0] * rab), float(u[1] * rab),
                          float(u[2] * rab), par])
        else:
            sites.append([float(u[0] * raa), float(u[1] * raa),
                          float(u[2] * raa), sym])

    import hashlib

    payload = {"type": "alloy", "absorber": sym, "partner": par,
               "x_partner": xb, "r_aa": raa, "r_ab": rab,
               "n_total": n_tot, "lattice": lat}
    key = hashlib.md5(
        json.dumps(payload, sort_keys=True, default=str).encode()
    ).hexdigest()
    motif = Motif(
        id=f"alloy_{sym.lower()}_{par.lower()}_{key[:8]}",
        absorber=sym, kind="cluster", nuclearity="alloy",
        source={"type": "dynamic_alloy", "partner": par,
                "x_partner": xb, "r_aa": raa, "r_ab": rab,
                "n_total": n_tot, "lattice": lat},
        sites=sites,
        notes=(f"动态双金属合金假设：{sym}–{sym} {raa:.3f} Å ×{n_aa}，"
               f"{sym}–{par} {rab:.3f} Å ×{n_ab}（x_{par}={xb:.3f}，"
               f"{lat} 第一壳 {n_tot} 方向）；"
               "无序/有序/偏析由 Cowley 参数判定。"))
    return _finalize_dynamic_motif(motif, payload)


def parse_catalytic_cluster_spec(absorber: str, spec) -> Motif:
    """B13-1 把 ``str`` / ``dict`` / 现成 ``Motif`` 统一成动态催化 motif。

    字符串语法（逗号分隔，大小写不敏感）::

        SAC: "sac:N4@2.00,axial=O@2.10,C@2.90"
             "sac:N4@2.00,axial=O@2.10,C=none"
        DAC: "dac:Co@2.50,N4@2.00"
        合金: "alloy:Co@2.64,Pt@2.76,x=0.25,n=12"
             "alloy:partner=Co,r_ab=2.64,r_aa=2.76,x=0.25"

    字典语法见 ``build_sac_motif`` / ``build_dac_motif`` /
    ``build_alloy_motif`` 的参数名。
    """
    if isinstance(spec, Motif):
        if str(spec.absorber).capitalize() != str(absorber).capitalize():
            raise XAFSConfigError(
                "motif", f"催化 motif 的吸收元素 {spec.absorber!r} 与 "
                         f"--absorber {absorber!r} 不一致")
        return spec
    if isinstance(spec, dict):
        return _catalytic_from_dict(absorber, spec)
    if isinstance(spec, str):
        return _catalytic_from_dict(
            absorber, _parse_catalytic_string(spec, absorber))
    raise XAFSConfigError(
        "motif", f"无法解析催化模型 {spec!r}：需 str / dict / Motif")


def _catalytic_from_dict(absorber: str, spec: dict) -> Motif:
    if not isinstance(spec, dict):
        raise XAFSConfigError("motif", f"催化模型必须是字典，收到 {type(spec).__name__}")
    kind = str(spec.get("type", "sac")).strip().lower()
    if kind not in ("sac", "dac", "alloy"):
        raise XAFSConfigError(
            "motif", f"催化模型 type={kind!r} 只支持 sac / dac / alloy")
    unknown = set(spec) - _CATALYTIC_DICT_KEYS[kind]
    if unknown:
        raise XAFSConfigError(
            "motif", f"{kind} 规格含未知键 {sorted(unknown)}；"
                     f"允许 {sorted(_CATALYTIC_DICT_KEYS[kind])}")
    if kind == "alloy":
        if spec.get("partner") is None:
            raise XAFSConfigError(
                "motif", "alloy 规格必须给 partner（第二金属元素）")
        return build_alloy_motif(
            absorber, str(spec["partner"]),
            x_partner=float(spec.get("x_partner", 0.5)),
            r_aa=(None if spec.get("r_aa") is None else float(spec["r_aa"])),
            r_ab=(None if spec.get("r_ab") is None else float(spec["r_ab"])),
            n_total=int(spec.get("n_total", 12)),
            lattice=str(spec.get("lattice", "fcc")))
    if kind == "sac":
        axial = spec.get("axial")
        return build_sac_motif(
            absorber, ligand=spec.get("ligand", "N"),
            n_coord=int(spec.get("n_coord", 4)),
            bond_m_ligand=float(spec.get("bond_m_ligand", 2.00)),
            axial=(None if axial in (None, "", "none", "None") else str(axial)),
            axial_bond=float(spec.get("axial_bond", 2.10)),
            carbon_shell=bool(spec.get("carbon_shell", True)),
            ring_radius=float(spec.get("ring_radius", 2.90)))
    if spec.get("partner") is None:
        raise XAFSConfigError("motif", "dac 规格必须给 partner（M₂ 元素）")
    return build_dac_motif(
        absorber, str(spec["partner"]),
        bond_length=float(spec.get("bond_length", 2.50)),
        ligand=spec.get("ligand", "N"),
        n_coord=int(spec.get("n_coord", 4)),
        ligand_bond=float(spec.get("ligand_bond", 2.00)))


def _parse_catalytic_string(spec: str, absorber: str | None = None) -> dict:
    """解析 SAC/DAC/合金精简字符串（不使用 eval/exec）。"""
    import re

    text = str(spec).strip()
    low = text.lower()
    if low.startswith("sac:"):
        kind, body = "sac", text[4:]
    elif low.startswith("dac:"):
        kind, body = "dac", text[4:]
    elif low.startswith("alloy:"):
        return _parse_alloy_string(absorber, text[6:], spec)
    else:
        raise XAFSConfigError(
            "motif", "催化模型字符串必须以 'sac:' / 'dac:' / 'alloy:' "
                     f"开头，收到 {spec!r}")
    out: dict = {"type": kind}
    tokens = [t.strip() for t in body.split(",") if t.strip()]
    if not tokens:
        raise XAFSConfigError("motif", f"催化模型字符串缺少几何参数: {spec!r}")

    for tok in tokens:
        m = re.fullmatch(r"axial\s*=\s*([A-Za-z]{1,2})\s*@\s*([0-9.]+)", tok, re.I)
        if m:
            out["axial"] = m.group(1).capitalize()
            out["axial_bond"] = float(m.group(2))
            continue
        m = re.fullmatch(r"C\s*=\s*(none|false|off|no)", tok, re.I)
        if m:
            out["carbon_shell"] = False
            continue
        m = re.fullmatch(r"C\s*@\s*([0-9.]+)", tok, re.I)
        if m:
            out["carbon_shell"] = True
            out["ring_radius"] = float(m.group(1))
            continue
        m = re.fullmatch(r"([A-Za-z]{1,2})(\d*)\s*@\s*([0-9.]+)", tok)
        if not m:
            raise XAFSConfigError(
                "motif", f"无法解析催化模型片段 {tok!r}（示例 N4@2.00 / "
                         f"axial=O@2.10 / C@2.90 / C=none）")
        el = m.group(1).capitalize()
        n_coord = int(m.group(2)) if m.group(2) else None
        dist = float(m.group(3))
        if el not in artemis.ATOMIC_NUMBER:
            raise XAFSConfigError("motif", f"未知元素 {el!r}（片段 {tok!r}）")
        if kind == "sac":
            if el in _LIGAND_ELEMENTS and "ligand" not in out:
                out["ligand"] = el
                if n_coord is not None:
                    out["n_coord"] = n_coord
                out["bond_m_ligand"] = dist
            elif el == "C":
                out["carbon_shell"] = True
                out["ring_radius"] = dist
            else:
                raise XAFSConfigError(
                    "motif", f"SAC 字符串无法归类片段 {tok!r}：配体须为 "
                             f"{sorted(_LIGAND_ELEMENTS)}，碳壳层写 C@R")
        else:
            if el in _LIGAND_ELEMENTS:
                out["ligand"] = el
                if n_coord is not None:
                    out["n_coord"] = n_coord
                out["ligand_bond"] = dist
            else:
                if "partner" in out and out["partner"] != el:
                    raise XAFSConfigError(
                        "motif", f"DAC 只能有一个金属 partner，收到 "
                                 f"{out['partner']} 与 {el}")
                out["partner"] = el
                out["bond_length"] = dist
    if kind == "dac" and "partner" not in out:
        raise XAFSConfigError(
            "motif", f"DAC 字符串缺少金属 partner（示例 'dac:Co@2.50,N4@2.00'）: {spec!r}")
    return out


def _parse_alloy_string(absorber: str | None, body: str, original) -> dict:
    """解析 ``alloy:Co@2.64,Pt@2.76,x=0.25,n=12`` 形式（无 eval/exec）。"""
    import re

    out: dict = {"type": "alloy"}
    sym = (str(absorber).capitalize() if absorber is not None else None)
    tokens = [t.strip() for t in str(body).split(",") if t.strip()]
    if not tokens:
        raise XAFSConfigError(
            "motif", f"合金字符串缺少几何参数: {original!r}")
    for tok in tokens:
        m = re.fullmatch(r"(?:x|x_partner)\s*=\s*([0-9.]+)", tok, re.I)
        if m:
            out["x_partner"] = float(m.group(1))
            continue
        m = re.fullmatch(r"(?:n|n_total)\s*=\s*(\d+)", tok, re.I)
        if m:
            out["n_total"] = int(m.group(1))
            continue
        m = re.fullmatch(r"partner\s*=\s*([A-Za-z]{1,2})", tok, re.I)
        if m:
            out["partner"] = m.group(1).capitalize()
            continue
        m = re.fullmatch(r"(r_aa|r_ab)\s*=\s*([0-9.]+)", tok, re.I)
        if m:
            out[m.group(1).lower()] = float(m.group(2))
            continue
        m = re.fullmatch(r"lattice\s*=\s*(fcc|hcp|bcc)", tok, re.I)
        if m:
            out["lattice"] = m.group(1).lower()
            continue
        m = re.fullmatch(r"([A-Za-z]{1,2})\s*@\s*([0-9.]+)", tok)
        if not m:
            raise XAFSConfigError(
                "motif", f"无法解析合金模型片段 {tok!r}（示例 "
                         f"'alloy:Co@2.64,Pt@2.76,x=0.25,n=12'）")
        el = m.group(1).capitalize()
        dist = float(m.group(2))
        if el not in artemis.ATOMIC_NUMBER:
            raise XAFSConfigError("motif", f"未知元素 {el!r}（片段 {tok!r}）")
        if sym is not None and el == sym:
            out["r_aa"] = dist
            continue
        if "partner" in out and out["partner"] != el:
            raise XAFSConfigError(
                "motif", f"合金模型只能有一个 partner，收到 "
                         f"{out['partner']} 与 {el}")
        out["partner"] = el
        out["r_ab"] = dist
    if "partner" not in out:
        raise XAFSConfigError(
            "motif", f"合金字符串缺少 partner（示例 "
                     f"'alloy:Co@2.64,Pt@2.76,x=0.25,n=12'）: {original!r}")
    return out


if __name__ == "__main__":
    for el_dir in sorted(LIB_ROOT.glob("*")):
        if not el_dir.is_dir() or el_dir.name.startswith("_"):
            continue
        for m in load_library(el_dir.name):
            fp = "  ".join(f"{s['el']}@{s['r']:.2f}x{s['n']}"
                           for s in (m.shell_expect or [])[:6])
            print(f"[{m.absorber}] {m.id:<22} {m.kind:<7} {m.nuclearity:<9} "
                  f"{len(m.sites or [])} sites | {fp}")

    # ---- 工具自测（用 Fe 种子库，不产生新注册）----
    print("\n== 工具自测 ==")
    _lib = {m.id: m for m in load_library("Fe")}

    # 替位：hematite 的 Fe 位换成 Co（BOND_REF Co-O oct=1.93），
    # 层内劈裂必须保留（1.94→1.93，2.11→2.099）
    sub = substitute(_lib["hematite"], "Fe", "Co", [0.0, 0.0, 0.3553],
                     coord="oct")
    sh = geometric_shells(sub.sites)
    assert abs(sh[0]["r"] - 1.93) < 0.03 and abs(sh[1]["r"] - 2.099) < 0.03, sh[:2]
    print(f"  替位: hematite→Co 第一壳 {sh[0]['r']:.3f}/{sh[1]['r']:.3f} Å"
          f"（期望 1.930/2.099，劈裂保留）")

    # 纳米模型：bcc 2 nm → 各壳层 scale ∈ (0,1]
    model = nanoparticle_model(_lib["fe_bcc"], 2.0)
    s1 = model["shells"][0]["scale"]
    assert 0.3 < s1 < 1.0, model["shells"][:2]
    print(f"  纳米模型: bcc 2 nm 吸收原子 {model['n_absorber']} 个，"
          f"壳层缩放 {[s['scale'] for s in model['shells'][:3]]}")

    # occupancy：lepidocrocite 全占据 → n_occ == n
    _sites, _occs = expand_full(_lib["lepidocrocite"])
    sho = geometric_shells(_sites, occs=_occs)
    assert all(abs(x["n_occ"] - x["n"]) < 1e-6 for x in sho)
    print(f"  occupancy: 全占据样品 n_occ==n 通过（{len(_occs)} 近邻）")
    print("自测通过")
