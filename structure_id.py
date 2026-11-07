# -*- coding: utf-8 -*-
"""
structure_id.py — S6.5 结构识别（指纹匹配 + 覆盖清单）
================================================================
在 S5/S6 之后、CN 拟合之前回答："实验 χ(R) 长得像库里的谁？"

方法（复数 χ(R) 匹配，不是经验相位偏移）:
    每个 motif → FEFF 理论 χ(k)（自然简并度、s02、nleg≤4）
    → σ² 网格 × ΔE0 网格 → 与实验同窗口做 xftf
    → 复数谱残差：score = min_α ‖χe(R) − α·χt(R)‖² / ‖χe(R)‖²
      (R∈[rmin,rmax]，α 为复数尺度，同时吸收振幅与全局相位)
    score 越小越像；Δscore < INDIST_TOL 的候选**并列报出**——
    该分辨窗口内不可分辨时不装确定。

覆盖清单（机制 2 的雏形）:
    实验 |χ(R)| 里每个显著峰（≥15% 窗口最大值、prominence ≥2%、
    半高宽 ≥0.08 Å）必须在最优候选的理论谱里有对应峰
    （|ΔR| ≤ 0.15 Å），否则列入 uncovered_peaks——未被解释的特征
    意味着假设空间不完备，应扩库而不是忽略。prominence/宽度门槛
    用于排除低 R 背景鼓包与单点噪声尖刺（Fe 箔 1.595 Å 峰高 15.6%
    但 prominence 仅 1.6%；两扫描合并的 Co 合成数据在 1.626 Å 出现
    FWHM≈0.05 Å 的假峰）。
    mixture（多位点吸收）例外：理论 χ 是各位点 χ 的加权和，复合谱
    中分量峰受干涉位移/压制（磁铁矿 tet 位 Fe-Fe@3.48 峰在复合谱中
    降到 14.7% 被滤掉、1.35 Å 只余斜坡）——覆盖判定复合谱 ∪ 各位
    点分量谱，任一位点解释即算覆盖。

边界（诚实声明）:
    - 指纹是"候选排序"，不是"结构确证"；最终 CN 仍须逐候选全拟合
    - `partial` 表示最优候选未解释全部显著峰（多位点/混合相/库缺条目），
      此时不做自动多壳层 CN，只做 L1 第一壳
    - mixture motif（多位点吸收）是各位点 χ 的加权求和，共用
      S₀²/σ²/ΔE₀——位点间 σ² 差异不被指纹分辨，定量时是系统误差来源
    - cluster motif 若含共线构型（如平面 FeN₄），MS 路径计入——
      这里取 nleg ≤ 4 的全部路径，不做 nleg==2 截断
    - 弱结晶/非晶（如水铁矿）无可靠 motif，识别预期落空 → 只报第一壳

自测: python structure_id.py
      用 lepidocrocite.000 跑全流程，lepidocrocite 必须排名第一。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import artemis
import athena
import structure_lib
from absorber import get_absorber
from xafs_core import XAFSConfigError, XAFSDataError, XAFSError

INDIST_TOL = 0.05      # Δscore 小于此 → 不可分辨，并列报告
PEAK_H = 0.15          # 覆盖清单的显著峰阈值（占窗口峰高比例）
PEAK_MATCH = 0.15      # 峰位匹配容差 (Å)
PEAK_WIDTH_MIN = 0.08  # 显著峰最小半高宽 (Å)：窄于此的是噪声尖刺，不是壳层
                       # （EXAFS 的 R 分辨率约 π/Δk≈0.26 Å，真实壳层不会这么窄）
PEAK_PROM = 0.02       # 最小 prominence（占窗口峰高比例）：去除骑在
                       # 第一壳肩上的低 R 背景鼓包（Fe 箔 1.595 Å 仅 1.6%）
PEAK_EDGE = 0.10       # 峰位距窗口边缘的最小距离 (Å)：窗口截断处的
                       # “峰”是 FT 边界效应，位置和宽度都不可靠
FIT_R = (1.0, 4.0)     # 指纹比较的 R 窗口 (Å)
E0_GRID = tuple(range(-8, 10, 2))   # 每候选扫描的 ΔE0 网格 (eV)
SIGMA2_GRID = (0.003, 0.006, 0.009)  # 指纹比较的 σ² 网格 (Å²)
CONFIDENT_SCORE = 0.16  # 复数谱相对残差低于此才可能判为 confident
CONFIDENT_MARGIN = 0.03  # 最优与次优的最小差距；过小判 ambiguous
                         # （0.02 会让 Hedenbergite 这类非库矿物误判 confident；
                         #  0.03 下真实磁铁矿/针铁矿/方铁矿仍保留 margin）
L1_TRIGGER = CONFIDENT_SCORE  # 兼容旧接口；新逻辑见 status


C_K2E = 0.2625         # k²[Å⁻²] = C · (E−E0)[eV]

S02_MIN = 0.4          # B18-2：运行时 S₀² 覆盖的合法区间
S02_MAX = 1.3


def validate_s02_override(s02) -> float:
    """校验运行时 S₀² 覆盖值（0.4–1.3），绝不改 absorber.py 静态表。"""
    try:
        val = float(s02)
    except (TypeError, ValueError) as exc:
        raise XAFSConfigError(
            "structure_id", f"s02={s02!r} 必须是数值") from exc
    if not np.isfinite(val) or not (S02_MIN <= val <= S02_MAX):
        raise XAFSConfigError(
            "structure_id",
            f"s02={s02!r} 超出合法区间 [{S02_MIN}, {S02_MAX}]")
    return val


def absorber_with_s02(symbol: str, s02) -> object:
    """返回带运行时 S₀² 的 AbsorberConfig 副本（dataclasses.replace）。"""
    import dataclasses
    return dataclasses.replace(
        get_absorber(symbol), s02=validate_s02_override(s02))


def theory_chi(motif, k, absorber, sigma2=0.003):
    """motif → FEFF 理论 χ(k)（含 MS，nleg≤4，自然简并度，求和一次）。

    kind=mixture（多位点吸收）：χ = Σ 权重·χ_位点。物理上这是正确的
    线性组合（不同吸收位点的 χ 独立产生后按占据加权求和），前提是
    各位点共用 S₀²/σ²/ΔE₀——对 N 简并的指纹比较足够，定量时各位点
    σ² 的差异是系统误差来源之一。

    注意：larch 的 FEFF path 对象用 ``nleg`` 表示散射腿数；早期版本误写
    ``legs``，``getattr`` 的默认值让所有路径都通过过滤，MS 截断实际失效。
    这里同时兼容只暴露 ``legs`` 的旧对象，但优先读 ``nleg``。
    """
    from larch.xafs import feffpath
    comps = getattr(motif, "components", None)
    if comps:
        lib = {m.id: m for m in structure_lib.load_library(motif.absorber)}
        total = None
        for comp in comps:
            chi_c = theory_chi(lib[comp["id"]], k, absorber, sigma2=sigma2)
            w = float(comp["weight"])
            total = w * chi_c if total is None else total + w * chi_c
        return total
    cache = structure_lib.ensure_feff(motif, edge=absorber.edge)
    paths = []
    for f in sorted(cache.glob("feff*.dat")):
        p = feffpath(str(f), s02=absorber.s02 or 1.0, sigma2=sigma2,
                     e0=0.0, deltar=0.0)
        nleg = getattr(p, "nleg", getattr(p, "legs", 2))
        if nleg <= 4 and getattr(p, "degen", 1) >= 1 \
                and getattr(p, "reff", 0) <= 5.5:
            p._calc_chi()
            paths.append(p)
    if not paths:
        raise XAFSDataError("structure_id", f"{motif.id}: 缓存内无有效路径")
    return artemis.exafs_model_k(paths, k)


def shift_e0(k, chi, de0):
    """
    用 k 轴重采样模拟 ΔE0 平移：k'² = k² + C_K2E·ΔE0。

    为什么不用 larch 路径参数：裸 feffpath 的 e0 属性在部分 larch
    版本里不被 _calc_chi() 消费（实测 χ 完全不变，见开发记录），
    而 ΔE0 的物理效果就是 k 网格的平方根平移——直接在求和后的
    χ(k) 上重采样，既不依赖 larch 内部行为，也省掉逐路径重算。
    """
    ks = np.sqrt(np.clip(k * k + C_K2E * de0, 0.0, None))
    return np.interp(ks, k, chi, left=0.0, right=0.0)


def _peaks(r, mag, mask, height_frac=PEAK_H,
           width_min: float = PEAK_WIDTH_MIN,
           prom_frac: float = PEAK_PROM,
           edge_margin: float = PEAK_EDGE):
    """窗口内显著峰位：高度、prominence、宽度三重门槛。

    - 高度 ≥ height_frac·max：排除小峰；
    - prominence ≥ prom_frac·max：排除骑在强峰肩上的低 R 背景鼓包
      （Fe 箔 1.595 Å 峰高 15.6% 但 prominence 仅 1.6%，不是壳层）；
    - 峰位离窗口边缘 ≥ edge_margin：窗口截断处的“峰”是 FT 边界效应，
      位置/宽度都不可靠（Vivianite 3.958 Å 就在 4.0 Å 边缘）；
    - 半高宽 ≥ width_min：排除单点噪声尖刺（两扫描合并的 Co 合成
      数据在 1.626 Å 出现 FWHM≈0.05 Å 的尖刺，单扫描时不存在）。
    EXAFS 的径向分辨率由 k 窗口决定，真实壳层峰不可能窄于 ~0.1 Å。
    """
    from scipy.signal import find_peaks, peak_prominences, peak_widths
    m = mag.copy()
    m[~mask] = 0.0
    if m.max() <= 0:
        return np.array([])
    peak_max = float(m.max())
    idx, _ = find_peaks(m, height=height_frac * peak_max)
    if not len(idx):
        return np.array([])
    prom = peak_prominences(mag, idx)[0]
    idx = idx[prom >= prom_frac * peak_max]
    if not len(idx):
        return np.array([])
    r_lo = float(np.min(r[mask])) + edge_margin
    r_hi = float(np.max(r[mask])) - edge_margin
    idx = idx[(r[idx] >= r_lo) & (r[idx] <= r_hi)]
    if not len(idx):
        return np.array([])
    if len(r) > 1:
        dr = float(np.median(np.diff(r)))
        # 宽度必须在**未置零**的原始 |χ(R)| 上算：把窗口外置 0 会
        # 改变峰两侧的局部基线，窗口边缘的窄尖刺会被算成 0.6 Å 宽。
        widths = peak_widths(mag, idx, rel_height=0.5)[0] * dr
        idx = idx[widths >= width_min]
    return r[idx]


def shape_score(exp_mag, th_mag, mask):
    """最优幅度缩放下的相对残差平方——只比形状，不比绝对强度。"""
    e, t = exp_mag[mask], th_mag[mask]
    alpha = float(np.dot(e, t) / np.dot(t, t)) if np.dot(t, t) > 0 else 0.0
    den = float(np.dot(e, e))
    return float(np.sum((e - alpha * t) ** 2)) / den if den > 0 else np.inf


def complex_score(exp_chi_r, th_chi_r, mask):
    """复数 χ(R) 的相对残差——同时比较幅度与相位。

    允许一个复数尺度 α（吸收整体振幅与全局相位），但**不允许**逐点
    调整相位；因此比只看 |χ(R)| 更能区分径向峰位相近、但相位/配位
    不同的候选。真实数据上这是识别磁铁矿与绿锈类结构的关键。
    """
    e, t = exp_chi_r[mask], th_chi_r[mask]
    den = float(np.vdot(e, e).real)
    if den <= 0:
        return np.inf
    t2 = float(np.vdot(t, t).real)
    if t2 <= 0:
        return np.inf
    alpha = np.vdot(t, e) / t2
    residual = e - alpha * t
    return float(np.vdot(residual, residual).real / den)


def identify(ft, k, absorber_symbol: str | None = None, kmin: float = 3.0,
             kmax: float = 12.0, kweight: int = 2, fit_r=FIT_R,
             verbose: bool = True, s02: float | None = None) -> dict:
    """
    对实验 ft（FourierTransform）跑全库指纹匹配。

    k 必须与 ft 同源（同一 χ(k) 网格）——理论谱在此网格上计算，
    保证比较严格同窗口、同网格。

    ``s02``（B18-2）：运行时覆盖 S₀²（0.4–1.3），用于未标定元素的
    显式解锁；不修改全局 ``absorber.ABSORBERS`` 静态表。

    返回 dict:
        ranked            按 score 升序的候选表
        indistinguishable 与最优 Δscore < INDIST_TOL 的并列候选
        uncovered_peaks   最优候选解释不了的实验峰位 (Å)。mixture 候选
                          按"复合谱 ∪ 各分量谱分别解释"判定——χ 是
                          各位点 χ 的加权和，复合峰位/相对峰高受干涉
                          位移，与单一 motif 的判定语义不同
    """
    if not absorber_symbol:
        raise XAFSConfigError(
            "structure_id", "未指定吸收元素——S6.5 需要它选 motif 库与 S₀²；"
                            "用 process_sample(absorber=...) 或 --absorber 提供")
    absorber = (get_absorber(absorber_symbol) if s02 is None
                else absorber_with_s02(absorber_symbol, s02))
    if not absorber.calibrated:
        raise XAFSDataError(
            "structure_id", f"{absorber.symbol} 未标定 S₀²，"
                            f"识别结果不能用于 CN 定量，拒绝运行")
    library = structure_lib.load_library(absorber.symbol)
    if not library:
        raise XAFSDataError("structure_id",
                            f"_structure_lib/{absorber.symbol} 为空，先建库")
    mask = (ft.r >= fit_r[0]) & (ft.r <= fit_r[1])
    exp_peaks = _peaks(ft.r, ft.magnitude, mask)

    ranked = []
    theory_cache = {}
    comp_cache = {}     # mixture：best_combo 下各位点自己的理论 |χ(R)|
    for motif in library:
        comps = getattr(motif, "components", None)
        if comps:
            lib = {m.id: m for m in structure_lib.load_library(motif.absorber)}
            comp_objs = [lib[c["id"]] for c in comps]
            weights = [float(c["weight"]) for c in comps]
        else:
            comp_objs, weights = None, None
        best = None       # 复数谱最优：用于排序
        best_mag = None   # |χ(R)| 最优：用于峰覆盖清单
        for sigma2 in SIGMA2_GRID:
            chi_comps = ([theory_chi(m, k, absorber, sigma2=sigma2)
                          for m in comp_objs] if comp_objs else [])
            chi0 = (sum(w * chi for w, chi in zip(weights, chi_comps))
                    if comp_objs else theory_chi(motif, k, absorber,
                                                 sigma2=sigma2))
            for de0 in E0_GRID:
                chi = shift_e0(k, chi0, de0)
                th = athena.ft_to_r_space(k, chi, None, kmin, kmax, kweight)
                th_mag = np.interp(ft.r, th.r, th.magnitude)
                th_complex = (
                    np.interp(ft.r, th.r, th.chi_r.real)
                    + 1j * np.interp(ft.r, th.r, th.chi_r.imag))
                s = complex_score(ft.chi_r, th_complex, mask)
                s_mag = shape_score(ft.magnitude, th_mag, mask)
                if best is None or s < best["score"]:
                    best = {"score": s, "de0": de0, "sigma2": sigma2}
                if best_mag is None or s_mag < best_mag["score"]:
                    # 各分量谱仅在成为新 |χ(R)| 最优时才算（同 de0/σ²）：
                    # 27 组合×每分量全存纯属浪费，且每分量只需一次 FT。
                    cmags = None
                    if comp_objs:
                        cmags = []
                        for chi0_c in chi_comps:
                            th_c = athena.ft_to_r_space(
                                k, shift_e0(k, chi0_c, de0), None,
                                kmin, kmax, kweight)
                            cmags.append(np.interp(ft.r, th_c.r,
                                                   th_c.magnitude))
                    best_mag = {"score": s_mag, "mag": th_mag, "comp": cmags}
        # 峰覆盖是幅度谱的性质：用 magnitude 最优的理论谱检查，
        # 避免“复数最优但幅度峰位略偏”造成假 uncovered。
        theory_cache[motif.id] = best_mag["mag"]
        comp_cache[motif.id] = best_mag.get("comp") or []
        ranked.append({"id": motif.id, "kind": motif.kind,
                       "nuclearity": motif.nuclearity,
                       "score": round(best["score"], 4),
                       "score_mag": round(best_mag["score"], 4),
                       "sigma2": best["sigma2"], "de0": best["de0"]})
    ranked.sort(key=lambda d: d["score"])

    best = ranked[0]
    # 最优候选自己不参与“并列”：旧版把 best 自身也放进 indist，
    # 导致 indistinguishable 永远非空，下游无法判断是否真有歧义。
    indist = [d["id"] for d in ranked[1:]
              if d["score"] - best["score"] < INDIST_TOL]
    second = ranked[1]["score"] if len(ranked) > 1 else None
    # 只有一条候选时 margin 无定义：写 None 而不是 inf，
    # 否则 meta.json 会出现非标准 JSON 的 Infinity。
    margin = (float(second - best["score"]) if second is not None else None)

    th_peaks = _peaks(ft.r, theory_cache[best["id"]], mask)
    # mixture（多位点吸收）：χ = Σ w·χ_位点，复合谱里分量峰互相干涉——
    # 实测磁铁矿 tet 位 Fe-Fe@3.48 的峰在复合谱中相对峰高 14.7%<15% 被
    # 过滤、1.35 Å 处只余上升斜坡无局部极大。这是求和干涉的产物，
    # 不是模型缺壳层，故覆盖按"复合谱 ∪ 各分量分别解释"判定：
    # 任一位点在 PEAK_MATCH 内有显著峰，该实验峰即被模型解释。
    for cmag in (comp_cache.get(best["id"]) or []):
        th_peaks = np.unique(np.concatenate(
            [th_peaks, _peaks(ft.r, cmag, mask)]))
    uncovered = [round(float(p), 3) for p in exp_peaks
                 if not np.any(np.abs(th_peaks - p) <= PEAK_MATCH)]
    coverage_ok = not uncovered

    # 四态判定：低分 + 足够 margin 只是必要条件；最优候选还必须
    # 解释全部显著峰。覆盖失败时给 partial——它可能仍是正确物相，
    # 但模型不完整（多位点/混合相/库缺条目），不足以直接做 CN。
    if best["score"] > CONFIDENT_SCORE:
        status = "no_match"
    elif margin is not None and margin < CONFIDENT_MARGIN:
        status = "ambiguous"
    elif not coverage_ok:
        status = "partial"
    else:
        status = "confident"
    result = {"ranked": ranked, "indistinguishable": indist,
              "uncovered_peaks": uncovered, "coverage_ok": coverage_ok,
              "status": status,
              "confidence": {"best_score": best["score"],
                             "second_score": second,
                             "margin": margin,
                             "score_limit": CONFIDENT_SCORE,
                             "margin_limit": CONFIDENT_MARGIN},
              "meta": {"absorber": absorber.symbol,
                       "lib_version": structure_lib.LIB_VERSION,
                       "n_motifs": len(library),
                       "fit_r": list(fit_r), "kmin": kmin, "kmax": kmax,
                       "score": "complex chi(R), sigma2 grid",
                       "sigma2_grid": list(SIGMA2_GRID)}}
    if not coverage_ok:
        result["warning"] = (
            f"最优候选未解释实验峰 {uncovered} Å；模型可能不完整"
            f"（多位点/混合相/库缺条目），confident 不代表结构确证")
    if verbose:
        print(f"  [S6.5] 结构识别（库 {len(library)} 条，"
              f"R∈[{fit_r[0]},{fit_r[1]}] Å）")
        for d in ranked:
            tag = " ★" if d["id"] == best["id"] else (
                " ~" if d["id"] in indist else "")
            print(f"         {d['id']:<18} {d['nuclearity']:<8} "
                  f"score={d['score']:.3f}  ΔE0={d['de0']:+.0f} eV"
                  f"  σ²={d['sigma2']:.3f}{tag}")
        margin_txt = "n/a（仅 1 条候选）" if margin is None else f"{margin:.3f}"
        print(f"         判定: {status}（best={best['score']:.3f}, "
              f"margin={margin_txt}）")
        if indist:
            print(f"         ⚠ 不可分辨并列: {indist}")
        if uncovered:
            print(f"         ⚠ 未解释实验峰 {uncovered} Å → 假设空间不完备，应扩库")
    return result


def identify_result(res, absorber_symbol: str | None = None, **kw):
    """PipelineResult 入口：元素优先取参数，缺省回落 res.absorber_symbol。

    识别结果带 ``status``:
        confident   best score 低、margin 足够，且解释全部显著峰
                    （mixture 候选按"复合谱∪各分量分别解释"判覆盖）
        partial     best 分数/margin 通过，但仍有未解释峰——模型不完整
                    （多位点残缺/混合相/库缺条目），不能直接做 CN
        ambiguous   best 与次优过于接近，不能唯一确定
        no_match    best score 高于阈值，库内很可能没有正确结构

    partial / ambiguous / no_match 会自动尝试 L1 第一壳拟合；L1 只给
    第一壳，长程结构保持未知。L1 的配体/几何可用 ``l1_ligand`` /
    ``l1_geometry`` 覆盖（默认 O / oct，适用于多数氧化物，但不是通用答案）。
    """
    l1_ligand = kw.pop("l1_ligand", "O")
    l1_geometry = kw.pop("l1_geometry", "oct")
    vary_c3 = kw.pop("vary_c3", False)
    fit_space = kw.pop("fit_space", "r")
    s02 = kw.pop("s02", None)
    absorber_symbol = absorber_symbol or getattr(res, "absorber_symbol", None)
    out = identify(res.ft, res.chi_k.k, absorber_symbol=absorber_symbol,
                   kmin=res.ft.kmin, kmax=res.ft.kmax,
                   kweight=res.ft.kweight, s02=s02, **kw)
    # L1 降级：库内无相似候选或候选不可分辨时，合成团簇只拟第一壳——
    # 弱结晶/非晶/库外体系唯一诚实的定量路径（长程不可知）。
    if out["status"] != "confident":
        try:
            out["l1_first_shell"] = l1_first_shell_fit(
                res, absorber_symbol=absorber_symbol,
                ligand=l1_ligand, geometry=l1_geometry,
                vary_c3=vary_c3, fit_space=fit_space, s02=s02,
                verbose=kw.get("verbose", True))
            margin = out["confidence"]["margin"]
            margin_txt = "n/a" if margin is None else f"{margin:.3f}"
            out["l1_note"] = (
                f"识别状态 {out['status']}（best={out['ranked'][0]['score']:.3f}, "
                f"margin={margin_txt}）：仅第一壳可定量，长程结构未知")
        except XAFSError as exc:
            # 识别本身仍应返回候选表；L1 失败作为可读诊断，不吞掉。
            out["l1_error"] = str(exc)
    return out


# ================================================================
# L1 降级：无结构第一壳拟合（合成团簇）
# ================================================================

_GEOM = {"oct": [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
                 (0, 0, 1), (0, 0, -1)],
         "tet": [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]}


def _first_shell_peak(ft, rmin: float = 1.0, rmax: float = 2.6,
                      height_frac: float = 0.30,
                      absorber: str | None = None,
                      ligand: str | None = None,
                      phase_shift: float = 0.45) -> float:
    """返回 R 窗口内最可能的第一壳峰位（相位未修正的 |χ(R)|）。

    选择规则（按优先级）：
    1. **配体感知**：由共价半径和估预期键长，减去相位偏移得到预期
       表观峰位，在显著峰（高度 ≥ height_frac·max）里选最近者。
       这解决了两个真实反例：Fe 氧化物里 Fe–Fe 第二壳（2.55）比
       Fe–O 第一壳（1.35）更强，argmax 会选错；Pt L3 白线在低 R
       （1.2–1.4）产生显著伪峰，"最低显著峰"会选错。
    2. 半径缺失/无显著峰 → 退回**最低显著峰**（重金属第二壳不会
       被误当第一壳）。
    3. 仍无峰 → 窗口内 argmax。
    """
    from scipy.signal import find_peaks
    m = (ft.r >= rmin) & (ft.r <= rmax)
    if not m.any():
        return float("nan")
    rr, mag = ft.r[m], ft.magnitude[m]
    idx, _ = find_peaks(mag, height=height_frac * float(mag.max()))
    if len(idx) and absorber and ligand:
        r_exp = artemis.expected_bond_length(absorber, ligand)
        if r_exp is not None:
            app_exp = r_exp - phase_shift
            return float(rr[idx[int(np.argmin(np.abs(rr[idx] - app_exp)))]])
    if len(idx):
        return float(rr[idx[0]])
    return float(rr[int(np.argmax(mag))])


def synthetic_shell_motif(absorber_symbol: str, ligand: str = "O",
                          geometry: str = "oct", r0: float = 2.0):
    """合成第一壳团簇：吸收原子在原点 + n 个配体（oct=6 / tet=4）。"""
    import math
    if geometry not in _GEOM:
        raise XAFSConfigError(
            "structure_id", f"未知 L1 几何 {geometry!r}；可选 {sorted(_GEOM)}")
    dirs = _GEOM[geometry]
    norm = math.sqrt(sum(v * v for v in dirs[0]))
    sites = [[r0 * d[0] / norm, r0 * d[1] / norm, r0 * d[2] / norm, ligand]
             for d in dirs]
    return structure_lib.Motif(
        id=f"_l1_{absorber_symbol}{ligand}_{geometry}_{r0:.2f}",
        absorber=absorber_symbol, kind="cluster", nuclearity="m1",
        source={"type": "synthetic", "purpose": "L1 first-shell fallback"},
        sites=sites)


def l1_first_shell_fit(res, absorber_symbol: str | None = None,
                       ligand: str = "O", geometry: str = "oct",
                       r0: float | None = None,
                       vary_c3: bool | list | str = False,
                       fit_space: str = "r",
                       s02: float | None = None,
                       verbose: bool = True) -> dict:
    """
    无结构第一壳拟合：合成团簇的 SS 路径拟 N/R/σ²/ΔE0。

    适用：弱结晶/非晶/库外体系——此时只有第一壳可信，
    长程拟合是强行给数据安结构。r0 默认取**配体感知选出的第一壳峰**
    +0.45 Å（相位粗估）：由共价半径和估预期键长，在显著峰里选最接近
    者（避免把 Fe 氧化物的 Fe–Fe 第二壳或 Pt L3 的低 R 白线伪峰当
    第一壳）；拟合会精修 ΔR，初值偏差不进最终结果。R 窗上限随表观
    峰扩展（≥2.3，≤3.0 Å），长键金属不会被截在窗外。

    ``ligand`` / ``geometry`` 是显式结构假设，默认 O/oct 只适用于多数
    氧化物；硫化物、氯化物、四面体体系必须由调用方覆盖。输出里会保留
    这两个假设，避免把 L1 的 N 误解成“无模型 CN”。
    """
    absorber_symbol = absorber_symbol or getattr(res, "absorber_symbol", None)
    if not absorber_symbol:
        raise XAFSConfigError("structure_id",
                              "L1 第一壳拟合需要吸收元素（--absorber）以取 S₀²")
    absorber = (get_absorber(absorber_symbol) if s02 is None
                else absorber_with_s02(absorber_symbol, s02))
    ligand = ligand.strip().capitalize()
    if ligand not in artemis.ATOMIC_NUMBER:
        raise XAFSConfigError(
            "structure_id", f"L1 配体 {ligand!r} 未收录于 artemis.ATOMIC_NUMBER")
    if geometry not in _GEOM:
        raise XAFSConfigError(
            "structure_id", f"未知 L1 几何 {geometry!r}；可选 {sorted(_GEOM)}")
    if r0 is None:
        peak = _first_shell_peak(res.ft, absorber=absorber.symbol,
                                 ligand=ligand)
        if not np.isfinite(peak):
            peak = float(res.r1)
        r0 = float(peak) + 0.45
    motif = synthetic_shell_motif(absorber.symbol, ligand, geometry, r0)
    cache = structure_lib.ensure_feff(motif, edge=absorber.edge)
    rows = artemis.describe_paths(sorted(cache.glob("feff*.dat")), 40)
    ss = [r for r in rows if r["legs"] == 2 and abs(r["reff"] - r0) <= 0.35]
    if not ss:
        raise XAFSDataError("structure_id",
                            f"L1 团簇在 {r0:.2f}±0.35 Å 内无 SS 路径")
    datfiles = [cache / r["file"] for r in ss]
    # R 窗上限必须随第一壳表观峰走：固定 2.3 Å 会把长键（Pt–Pt 2.77 /
    # Fe–Fe 2.48）的模型峰截在窗外，曾给出 Rf=1.0 的伪失败。取表观峰
    # +0.15 Å（主瓣右肩），下限 2.3，上限 3.0（不放进第二壳）。
    rmax_fit = min(3.0, max(2.3, (r0 - 0.45) + 0.15))
    if isinstance(vary_c3, str):
        key = vary_c3.strip().lower()
        c3_arg = ([0] if key == "first" else
                  (True if key == "all" else None))
    elif vary_c3 is True:
        c3_arg = True
    elif isinstance(vary_c3, (list, tuple, set)):
        c3_arg = sorted({int(i) for i in vary_c3}) or None
    else:
        c3_arg = None

    def _run_fit(c3):
        return artemis.fit_r_space(
            res.chi_k, res.ft, datfiles,
            rmin=1.0, rmax=rmax_fit,
            kmin=res.ft.kmin, kmax=res.ft.kmax,
            kweight=res.ft.kweight, s02=absorber.s02,
            vary_delta_e0=True, vary_n=True,
            vary_c3=(c3 if c3 is not None else False),
            fit_space=fit_space)

    fit = _run_fit(c3_arg)
    c3_dropped_by_nidp = False
    if c3_arg and fit.n_variables >= fit.n_independent:
        c3_dropped_by_nidp = True
        if c3_arg is True or tuple(c3_arg) != (0,):
            fit = _run_fit([0])
        if fit.n_variables >= fit.n_independent:
            fit = _run_fit(None)

    # 同一第一壳可能被 FEFF 拆成多条 SS 路径（例如 oct 畸变）。
    # N 必须是壳层总和；旧版只读 paths[0] 会系统性低估配位数。
    n_total = float(sum(p.degeneracy for p in fit.paths))
    n_err = (float(np.sqrt(sum((p.n_err or 0.0) ** 2 for p in fit.paths)))
             if any(p.n_err is not None for p in fit.paths) else None)
    p_main = max(fit.paths, key=lambda p: p.degeneracy)
    out = {"geometry": geometry, "ligand": ligand, "r0": round(r0, 4),
           "N": round(n_total, 2),
           "N_err": round(n_err, 2) if n_err is not None else None,
           "R": round(p_main.r, 4), "R_err": round(p_main.delta_r_err, 4),
           "sigma2": round(p_main.sigma2, 4), "rf": round(fit.r_factor, 4),
           "delta_e0": round(fit.delta_e0, 2),
           "c3": (None if p_main.c3 is None else float(p_main.c3)),
           "c3_err": (None if p_main.c3_err is None
                      else float(p_main.c3_err)),
           "c3_at_bound": bool(any(p.c3_at_bound for p in fit.paths)),
           "c3_dropped_by_nidp": bool(c3_dropped_by_nidp),
           "n_paths": len(fit.paths), "n_independent": fit.n_independent,
           "n_variables": fit.n_variables,
           "warning": (f"仅第一壳可信；长程结构未知（L1 降级，"
                       f"假设 {absorber.symbol}–{ligand} {geometry}）")}
    if abs(out["delta_e0"]) >= 14.5:
        out["warning"] += ("；ΔE0 已接近 ±15 eV 边界，"
                           "第一壳模型/配体假设可能不充分")
    if verbose:
        print(f"  [L1] 结构不确定 → 第一壳拟合（{geometry} {ligand} 合成团簇）: "
              f"N={out['N']}±{out['N_err']}, R={out['R']} Å, "
              f"σ²={out['sigma2']}, Rf={out['rf']}")
        print(f"       ⚠ {out['warning']}")
    return out


def reference_element_k_peaks(
        absorber_symbol: str,
        elements: tuple[str, ...] = ("O", "Cl", "Fe"),
        r_bond: float = 2.10,
        sigma2: float = 0.005,
        kmin: float = 3.0,
        kmax: float = 12.0,
        kweight: int = 2,
        wt_sigma: float = 1.5,
        edge: str | None = None,
) -> dict[str, float]:
    """候选散射元素的理论小波包络峰位标尺（含 σ² 调制）。

    为什么需要它：小波 k_peak 不是元素指纹常数。包络峰位受
    ``k^w·exp(-2σ²k²)`` 与窗口截断共同调制——高无序的 Fe–Fe 壳
    （σ²≈0.015 Å²）理论峰可被压到 5–6 Å⁻¹，与低无序 Fe–O/Fe–Cl
    重合。盲判配体前必须用**同一组 (kmin,kmax,kweight,sigma,σ²)**
    的理论路径画出标尺，而不是背"3d 金属峰在 7 Å⁻¹"。

    实现：每个候选元素生成 MO6 八面体单壳团簇 → FEFF6 → 第一条
    SS 路径（s02=1, ΔR=0, 指定 σ²）→ 同一 WT 参数 → 在
    ``[r_bond-1.0, r_bond+0.6] Å`` 内定位 2D 峰并取一维切片 k_peak。
    FEFF 结果缓存在 ``structure_lib.CACHE_ROOT``，同一
    (absorber, edge, element, r_bond) 只跑一次。

    依赖铁律：`larch_wt` 在函数体内局部导入，模块顶层仍只有
    artemis/athena/structure_lib，无循环依赖。

    返回: ``{element: k_peak (Å⁻¹)}``。
    """
    from larch.xafs import feffpath

    import larch_wt
    from xafs_core import ChiK

    absorber = get_absorber(absorber_symbol)
    edge = edge or absorber.edge
    if edge not in artemis.HOLE_INDEX:
        raise XAFSConfigError(
            "WT", f"未知吸收边 {edge!r}；可选 {sorted(artemis.HOLE_INDEX)}")
    if not (r_bond > 0):
        raise XAFSConfigError("WT", f"r_bond={r_bond} 必须 > 0 Å")
    if not (sigma2 >= 0):
        raise XAFSConfigError("WT", f"sigma2={sigma2} 不能为负")
    if not elements:
        return {}

    cache_root = structure_lib.CACHE_ROOT
    rmax_wt = max(4.0, float(r_bond) + 1.0)
    k_grid = np.arange(0.0, float(kmax) + 2.0 + 1e-9, 0.05)
    out: dict[str, float] = {}
    for element in elements:
        el = str(element).strip().capitalize()
        if el not in artemis.ATOMIC_NUMBER:
            raise XAFSConfigError(
                "WT", f"未知散射元素 {element!r}，未收录于 artemis.ATOMIC_NUMBER")
        cache = cache_root / (
            f"wt_ref_{absorber.symbol}_{edge}_{el}_"
            f"{int(round(float(r_bond) * 100))}")
        if not list(cache.glob("feff*.dat")):
            sites = artemis.octahedron_sites(float(r_bond), symbol=el)
            inp = artemis.build_feff_input(
                sites, absorber=absorber.symbol, rmax=rmax_wt,
                title=f"wt_ref_{absorber.symbol}_{el}", edge=edge)
            artemis.run_feff6l(inp, cache)

        # 取 Reff 最接近目标键长的那条单散射路径
        candidates = []
        for f in sorted(cache.glob("feff*.dat")):
            path = feffpath(str(f))
            if int(getattr(path, "nleg", 2)) != 2:
                continue
            candidates.append((abs(float(path.reff) - float(r_bond)), f))
        if not candidates:
            raise XAFSDataError(
                "WT", f"{el}@{r_bond:.2f} Å 的 FEFF 缓存中没有 nleg=2 路径")
        _, datfile = min(candidates, key=lambda item: item[0])
        path = artemis.load_feff_path(datfile)
        artemis.set_path_parameters(
            path, s02=1.0, delta_r=0.0, sigma2=float(sigma2),
            delta_e0=0.0)
        chi = np.interp(k_grid, np.asarray(path.k), np.asarray(path.chi))
        ck = ChiK(k=k_grid, chi=chi,
                  background=np.zeros_like(k_grid),
                  mu_norm=np.zeros_like(k_grid), rbkg=1.0)
        wt = larch_wt.continuous_wavelet_transform(
            ck, sigma=float(wt_sigma), kmin=float(kmin), kmax=float(kmax),
            kweight=int(kweight), rmax=rmax_wt, normalize=False)
        # 理论上表观 R = r_bond − 相位偏移；窗口给足 ±1 Å 容差。
        _, r_peak, _ = larch_wt.locate_shell(
            wt, max(0.5, float(r_bond) - 1.0), float(r_bond) + 0.6)
        profile = larch_wt.extract_k_profile(wt, r_peak, half_width=0.30)
        out[el] = float(profile["k_peak"])
    return out


if __name__ == "__main__":
    import xafs_pipeline   # 仅自测需要；模块级不依赖归约链（保持依赖线性）
    root = Path(__file__).resolve().parent
    data = root / "data" / "练手数据" / "02_Fe_standards" / "Fe_lepidocrocite.000"
    res = xafs_pipeline.process_sample([data], absorber="Fe", verbose=True)
    out = identify_result(res)
    top = out["ranked"][0]["id"]
    print(f"\n自测: top1 = {top}（期望 lepidocrocite）")
    if top != "lepidocrocite":
        sys.exit("自测失败：标准样未命中自身 motif")

    # L1 机制直测（已知样品）：N 应 ≈6、R 应 ≈2.0
    l1 = l1_first_shell_fit(res, geometry="oct", verbose=True)
    assert 4.5 < l1["N"] < 7.5, l1
    assert 1.90 < l1["R"] < 2.15, l1
    print("自测通过")
