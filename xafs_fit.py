# -*- coding: utf-8 -*-
"""
xafs_fit.py — 拟合链编排（S6.5 结构识别 + S7 自动/指定结构拟合）
================================================================
归约链（xafs_pipeline，S1–S6）与拟合链在这里分界：

    PipelineResult 进 → S6.5 结构识别
                     → S7 拟合（识别到 confident 结构，或显式给 motif）

关键决策：
  - 识别状态不是 confident（ambiguous / partial / no_match）时，默认只做
    L1 第一壳拟合；长程结构未知。``allow_ambiguous=True`` 才允许用 Top-1
    强制拟合，并在结果里标记。
  - 自动拟合的路径按 Reff + 散射元素分壳层；同一壳层的多条 SS 路径
    合并为一条代表路径，简并度加总，拟合出的 N 是**壳层总配位数**。
    这是 Artemis 里“合并不可分辨壳层”的自动化版本，代表路径的近似
    仍会在结果里保留 warning。
  - S₀² 只从 absorber.py 的标定值取；未标定元素拒绝 CN 定量。
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np

import artemis
import cn_table
import structure_id
import structure_lib
from absorber import get_absorber
from structure_model import CNError, StructureModel, validate_formal
from xafs_core import XAFSConfigError, XAFSDataError, XAFSError

# 形式误差 ≠ 绝对精度。N 的误差棒只来自最小二乘协方差，
# 不包含 S₀² 传递、代表路径/壳层合并、N–σ²–ΔE₀ 相关和归一化偏差。
# 按 Fe 基准（Fe 箔第二壳 9.45 vs 晶体学 6；lepidocrocite Fe–Fe 7.6 vs 6）
# 实测，CN 的绝对系统误差约 ±20–40%（1σ，量级估计），跨体系比较时
# 不应把 ±0.3 当成“配位数已精确到小数第二位”。
CN_SYSTEMATIC_BUDGET = {
    "relative_1sigma": 0.30,
    "components": {
        "S02_transfer": 0.15,
        "representative_path_and_shell_merge": 0.25,
        "N_sigma2_e0_correlation": 0.15,
        "normalization": 0.10,
    },
    "note": ("n_err 是形式误差；CN 绝对值的系统误差约 ±20–40%，"
             "报告时应注明预算，不能把形式误差当绝对精度。"),
}

# 荧光自吸收：没有样品组成/厚度就无法做可验证的定量校正。
# 这里只做保守筛查：荧光模式下边跳变越大，表观振幅被压低越严重。
# 阈值是**风险门禁**，不是校正因子；高风险时拒绝 CN/相分数结论。
SELF_ABSORPTION_HIGH_EDGE_STEP = 1.0   # Δμ ≥ 1.0：浓样/厚样，必须校正后才定量
SELF_ABSORPTION_WARN_EDGE_STEP = 0.3   # 0.3–1.0：可能受影响，只警告


def _self_absorption_diag(res) -> dict:
    """荧光模式自吸收风险（不校正，只决定能否定量）。

    荧光信号 ∝ 吸收量，但出射荧光会被样品自身再吸收；浓度/厚度越大，
    EXAFS 振幅被压低越多，表现为 N 偏小。可靠校正需要样品组成、密度、
    厚度与几何，本管线不假设这些量，因此在荧光模式下：
      - Δμ ≥ 1.0 → risk=high，cn_gate=True（拒绝 CN/相分数定量）
      - 0.3 ≤ Δμ < 1.0 → risk=possible，只警告
      - Δμ < 0.3 → risk=low（稀样，自吸收通常可忽略）
    透射/文件自带 xmu 不受此门禁影响。
    """
    merged = getattr(res, "merged", None)
    mode = ((getattr(merged, "meta", None) or {}).get("mode")
            if merged is not None else None) or "unknown"
    norm = getattr(res, "norm", None)
    edge_step = float(getattr(norm, "edge_step", float("nan")))
    sa = getattr(res, "self_absorption", None)
    if sa and sa.get("corrected"):
        gate = not bool(sa.get("reliable", False))
        return {
            "mode": mode, "edge_step": edge_step, "corrected": True,
            "reliable": bool(sa.get("reliable", False)),
            "risk": ("corrected_unreliable" if gate else "corrected"),
            "cn_gate": gate,
            "method": sa.get("method"),
            "formula": sa.get("formula"),
            "thickness_um": sa.get("thickness_um"),
            "density": sa.get("density"),
            "anginp": sa.get("anginp"), "angout": sa.get("angout"),
            "line": sa.get("line"),
            "correction_factor_median": sa.get("correction_factor_median"),
            "note": ("已做 Booth 自吸收校正；可靠性取决于总化学式、厚度、"
                     "有效密度与角度的准确性（" + str(sa.get("note", "")) + "）"
                     + ("；" + str(sa["warning"]) if sa.get("warning") else "")),
        }
    diag = {"mode": mode, "edge_step": edge_step, "corrected": False,
            "risk": "none", "cn_gate": False,
            "note": "非荧光模式，不涉及自吸收"}
    if mode != "fluor":
        return diag
    if not math.isfinite(edge_step):
        diag.update(risk="unknown", cn_gate=True,
                    note="荧光模式但边跳变缺失，无法判断自吸收，拒绝定量")
    elif edge_step >= SELF_ABSORPTION_HIGH_EDGE_STEP:
        diag.update(risk="high", cn_gate=True,
                    note=(f"荧光 Δμ={edge_step:.3f} ≥ "
                          f"{SELF_ABSORPTION_HIGH_EDGE_STEP:.1f}："
                          "浓样/厚样自吸收显著，未校正时振幅与 N 被压低，"
                          "需透射谱或组成/厚度已知的外部校正"))
    elif edge_step >= SELF_ABSORPTION_WARN_EDGE_STEP:
        diag.update(risk="possible",
                    note=(f"荧光 Δμ={edge_step:.3f} 处于可能受自吸收影响的"
                          "区间；若 N 系统性偏低，优先怀疑此项"))
    else:
        diag.update(risk="low",
                    note="荧光稀样，自吸收通常可忽略")
    return diag


def _max_offdiag_correlation(fit) -> tuple[float, str | None, str | None]:
    """返回 (|ρ|max, 变量名 i, 变量名 j)；无相关矩阵时返回 (0, None, None)。"""
    corr = getattr(fit, "correlation", None)
    names = getattr(fit, "parameter_names", None)
    if corr is None or names is None:
        return 0.0, None, None
    c = np.array(corr, dtype=float, copy=True)
    if c.ndim != 2 or c.size == 0:
        return 0.0, None, None
    np.fill_diagonal(c, 0.0)
    ij = np.unravel_index(np.argmax(np.abs(c)), c.shape)
    return abs(float(c[ij])), names[ij[0]], names[ij[1]]


def systematic_budget(fit, flags: dict | None = None,
                      self_abs: dict | None = None,
                      truncated_shells: int = 0,
                      fit_weights: bool = False) -> dict:
    """按**本次拟合的实际条件**生成 CN 系统误差预算（1σ，相对值）。

    基础预算 0.30 是历史实测锚点（见 CN_SYSTEMATIC_BUDGET）。这里不重算
    基础项，只把**额外可归因的系统误差**与基础值按平方和合成：

      - Booth 自吸收已校正且 reliable：+0.10（化学式/厚度/角度误差残留）
      - 荧光自吸收未校正但处于 possible 区间：+0.15（振幅压低方向未知）
      - 最强参数相关 |ρ| ≥ 0.9：+0.10；0.8–0.9：+0.05
      - 因 N_idp 截断远壳层：+0.05
      - σ² 撞下界（形状失配）：+0.10

    宁可将不可归因项显式列出，也不把 0.30 当成放之四海皆准的常数。
    """
    flags = flags or {}
    self_abs = self_abs or {}
    base = float(CN_SYSTEMATIC_BUDGET["relative_1sigma"])
    extras: dict[str, float] = {}
    notes: list[str] = []

    if self_abs.get("corrected") and self_abs.get("reliable"):
        extras["self_absorption_residual"] = 0.10
        notes.append("Booth 校正后残留（化学式/厚度/角度误差）")
    elif (not self_abs.get("corrected")
          and self_abs.get("risk") == "possible"):
        extras["self_absorption_uncorrected"] = 0.15
        notes.append("荧光自吸收可能压低振幅且未校正")

    rho, ni, nj = _max_offdiag_correlation(fit)
    if rho >= 0.9:
        extras["parameter_correlation"] = 0.10
        notes.append(f"参数强相关 |ρ|={rho:.3f}"
                     + (f"（{ni} vs {nj}）" if ni else ""))
    elif rho >= 0.8:
        extras["parameter_correlation"] = 0.05
        notes.append(f"参数相关 |ρ|={rho:.3f}"
                     + (f"（{ni} vs {nj}）" if ni else ""))

    if truncated_shells > 0:
        extras["truncated_shells"] = 0.05
        notes.append(f"{truncated_shells} 条远壳层因 N_idp 截断")
    if flags.get("sigma2_at_bound"):
        extras["sigma2_at_bound"] = 0.10
        notes.append("σ² 撞下界，模型形状失配")

    total = math.hypot(base, *extras.values())
    return {
        "relative_1sigma": round(total, 4),
        "base_relative_1sigma": base,
        "extra_components": extras,
        "components": dict(CN_SYSTEMATIC_BUDGET["components"]),
        "max_parameter_correlation": round(rho, 4),
        "note": ("基础预算 "
                 f"±{100 * base:.0f}% 与本次额外系统项按平方和合成"
                 + ("：" + "；".join(notes) if notes else "（无额外项）")
                 + "。n_err 仍是形式误差，不能当绝对精度。"),
    }


def _first_shell_cn(fit_out) -> dict:
    """从 motif / L1 拟合产物中取第一壳 CN 摘要（只读，不做推断）。"""
    out = {"n_fit": None, "n0": None, "label": None, "reff": None,
           "cn_reliable": None, "mode": None}
    if not isinstance(fit_out, dict):
        return out
    out["mode"] = fit_out.get("mode")
    out["cn_reliable"] = fit_out.get("cn_reliable")
    if fit_out.get("mode") == "motif_weights":
        return out
    l1 = fit_out.get("l1")
    if fit_out.get("mode") == "l1_first_shell" and isinstance(l1, dict):
        out["n_fit"] = (float(l1["N"]) if l1.get("N") is not None else None)
        out["reff"] = (float(l1["R"]) if l1.get("R") is not None else None)
        out["label"] = "L1 first shell"
        return out
    fit_obj = fit_out.get("fit")
    paths = list(getattr(fit_obj, "paths", None) or [])
    if paths:
        p = paths[0]
        degen = getattr(p, "degeneracy", None)
        out["n_fit"] = float(degen) if degen is not None else None
        out["reff"] = float(getattr(p, "r", None)) \
            if getattr(p, "r", None) is not None else None
    shells = fit_out.get("shells") or []
    if shells:
        sh = shells[0]
        out["n0"] = (float(sh["degen"]) if sh.get("degen") is not None
                     else None)
        out["label"] = sh.get("label")
        if out["reff"] is None and sh.get("reff") is not None:
            out["reff"] = float(sh["reff"])
    return out


def xanes_exafs_cross_check(res, fit_out=None, *,
                            white_line=None,
                            white_line_window=(-6.0, 40.0),
                            thresholds: dict | None = None) -> dict:
    """XANES–EXAFS 联合软校验：白线/边步/归一化 vs 第一壳 CN。

    **只产生警告，不改变 ``cn_reliable``、不拒绝 CN、不写正式表。**
    物理上白线强度受氧化态、配位几何、自吸收与多重散射共同影响，
    不能由 CN 单独解释；因此本函数只做**保守极端值筛查**，并把
    ``heuristic=True`` / ``affects_cn_reliable=False`` 写进产物。

    检查项：
      1. 归一化平台是否在 0.85–1.15；
      2. 边步是否非正/过弱/过强（原始 μ 跳变，不是归一化后平台）；
      3. XANES 导数双峰 ``e0_jump_risk``；
      4. 白线拟合是否可靠、白线位置相对 E0 是否越界；
      5. K 边极端对比：高 CN + 极弱白线，或低 CN + 极强白线；
      6. 荧光自吸收风险（若 fit_out 带该诊断）。

    ``white_line`` 可注入预计算结果（测试/批量复用），避免重复拟合。
    """
    th = {
        "plateau_range": (0.85, 1.15),
        "edge_step_min": 0.05,
        "edge_step_max": 5.0,
        "white_line_area_low": 0.5,
        "white_line_area_high": 12.0,
        "high_cn": 8.0,
        "low_cn": 1.5,
        "white_line_offset_range": (-8.0, 40.0),
    }
    if thresholds:
        th.update(thresholds)

    norm = getattr(res, "norm", None)
    xf = dict(getattr(res, "xanes_features", None) or {})
    fit_out = fit_out if fit_out is not None else getattr(res, "fit", None)
    absorber = getattr(res, "absorber_symbol", None)
    edge = "K"
    if absorber:
        try:
            edge = get_absorber(absorber).edge
        except (XAFSConfigError, KeyError):
            edge = str(xf.get("edge") or "K")
    else:
        edge = str(xf.get("edge") or "K")

    wl = white_line
    if wl is None:
        if norm is None:
            wl = {"available": False, "reliable": False,
                  "note": "无 NormalizedScan，无法拟合白线"}
        else:
            try:
                import athena
                wl = athena.fit_whiteline_area(
                    norm, window=white_line_window)
                wl = dict(wl)
                wl["available"] = True
            except Exception as exc:  # noqa: BLE001 - 软校验不得中断主链
                wl = {"available": False, "reliable": False,
                      "note": f"白线拟合失败: {type(exc).__name__}: {exc}"}

    first = _first_shell_cn(fit_out)
    warnings: list[dict] = []
    checks: dict[str, dict] = {}

    def warn(code, message):
        warnings.append({"code": code, "message": message})

    # 1. 归一化平台
    plateau = float(getattr(norm, "plateau", float("nan"))) \
        if norm is not None else float("nan")
    p_lo, p_hi = th["plateau_range"]
    if not math.isfinite(plateau):
        checks["normalization_plateau"] = {
            "status": "unavailable", "value": None}
    elif p_lo <= plateau <= p_hi:
        checks["normalization_plateau"] = {
            "status": "ok", "value": round(plateau, 4)}
    else:
        checks["normalization_plateau"] = {
            "status": "warn", "value": round(plateau, 4)}
        warn("normalization_plateau_out_of_range",
             f"边后平台={plateau:.3f} 不在 {p_lo:.2f}–{p_hi:.2f}；"
             f"EXAFS 振幅/CN 与 XANES 白线均可能被归一化偏差污染")

    # 2. 边步
    estep = float(getattr(norm, "edge_step", float("nan"))) \
        if norm is not None else float("nan")
    if not math.isfinite(estep) or estep <= 0:
        checks["edge_step"] = {"status": "unavailable", "value": None}
    elif estep < th["edge_step_min"]:
        checks["edge_step"] = {"status": "warn", "value": round(estep, 4)}
        warn("edge_step_weak",
             f"原始边步={estep:.3g} 过弱（阈值 {th['edge_step_min']}）："
             f"样品过稀/过薄，EXAFS 信噪比与 CN 可靠性下降")
    elif estep > th["edge_step_max"]:
        checks["edge_step"] = {"status": "warn", "value": round(estep, 4)}
        warn("edge_step_strong",
             f"原始边步={estep:.3g} 过强（阈值 {th['edge_step_max']}）："
             f"厚样/浓样自吸收可能压低 EXAFS 振幅")
    else:
        checks["edge_step"] = {"status": "ok", "value": round(estep, 4)}

    # 3. E0 导数双峰
    if xf.get("e0_jump_risk") is True:
        checks["e0_jump_risk"] = {"status": "warn"}
        warn("e0_jump_risk",
             "XANES 上升沿导数双峰：E0 可能在两个台阶间跳变 3–6 eV，"
             "相位/振幅解释需复核")
    elif "e0_jump_risk" in xf:
        checks["e0_jump_risk"] = {"status": "ok"}
    else:
        checks["e0_jump_risk"] = {"status": "unavailable"}

    # 4. 白线可靠性与位置
    wl_reliable = bool(wl.get("reliable")) if isinstance(wl, dict) else False
    wl_area = (float(wl.get("white_line_area_norm") or 0.0)
               if isinstance(wl, dict) else None)
    wl_energy = (wl.get("white_line_energy")
                 if isinstance(wl, dict) else None)
    if not wl_reliable:
        checks["white_line_fit"] = {
            "status": "unavailable",
            "note": (wl or {}).get("note", "白线拟合不可靠")}
    else:
        checks["white_line_fit"] = {
            "status": "ok", "area_norm": round(wl_area, 4)}
        e0 = float(getattr(norm, "e0", float("nan"))) \
            if norm is not None else float("nan")
        if wl_energy is not None and math.isfinite(e0):
            offset = float(wl_energy) - e0
            lo, hi = th["white_line_offset_range"]
            if not (lo <= offset <= hi):
                checks["white_line_position"] = {
                    "status": "warn", "offset_eV": round(offset, 2)}
                warn("white_line_position_outlier",
                     f"白线位置 E_wl−E0={offset:+.1f} eV 超出 "
                     f"[{lo:+.0f}, {hi:+.0f}] eV；检查 E0 选择/边前背景/"
                     f"价态归属")
            else:
                checks["white_line_position"] = {
                    "status": "ok", "offset_eV": round(offset, 2)}

    # 5. K 边极端 CN–白线对比（L 边白线由 2p→5d 主导，不做此判据）
    n_fit = first.get("n_fit")
    if (str(edge).upper().startswith("K") and wl_reliable
            and n_fit is not None and first.get("cn_reliable") is not False):
        if n_fit >= th["high_cn"] and wl_area <= th["white_line_area_low"]:
            checks["cn_white_line_contrast"] = {
                "status": "warn", "n_fit": round(n_fit, 2),
                "white_line_area_norm": round(wl_area, 4)}
            warn("high_cn_low_white_line",
                 f"第一壳 N≈{n_fit:.1f} 但白线净面积仅 {wl_area:.2f} eV："
                 f"白线极弱与高 CN 同时出现，需排查归一化/自吸收/多相混合")
        elif n_fit <= th["low_cn"] and wl_area >= th["white_line_area_high"]:
            checks["cn_white_line_contrast"] = {
                "status": "warn", "n_fit": round(n_fit, 2),
                "white_line_area_norm": round(wl_area, 4)}
            warn("low_cn_high_white_line",
                 f"第一壳 N≈{n_fit:.1f} 但白线净面积达 {wl_area:.2f} eV："
                 f"强白线通常含氧化态/多重散射贡献，CN 可能被低估")
        else:
            checks["cn_white_line_contrast"] = {
                "status": "ok", "n_fit": round(n_fit, 2),
                "white_line_area_norm": round(wl_area, 4)}
    else:
        checks["cn_white_line_contrast"] = {"status": "unavailable"}

    # 6. 自吸收风险（fit_out 已给出诊断时）
    self_abs = (fit_out.get("self_absorption")
                if isinstance(fit_out, dict) else None)
    if isinstance(self_abs, dict):
        risk = self_abs.get("risk")
        if risk in ("high", "unknown", "corrected_unreliable"):
            checks["self_absorption"] = {"status": "warn", "risk": risk}
            warn("self_absorption_risk",
                 self_abs.get("note") or f"荧光自吸收风险={risk}")
        else:
            checks["self_absorption"] = {"status": "ok", "risk": risk}
    else:
        checks["self_absorption"] = {"status": "unavailable"}

    # 7. CN 可靠性信息（不改它，只记录）
    if isinstance(fit_out, dict) and first.get("cn_reliable") is False:
        checks["cn_reliability"] = {
            "status": "info", "cn_reliable": False}
    elif isinstance(fit_out, dict):
        checks["cn_reliability"] = {
            "status": "ok", "cn_reliable": first.get("cn_reliable")}
    else:
        checks["cn_reliability"] = {"status": "unavailable"}

    status = ("warning" if warnings
              else ("ok" if any(c["status"] == "ok"
                                for c in checks.values()) else "unavailable"))
    return {
        "enabled": True,
        "heuristic": True,
        "affects_cn_reliable": False,
        "status": status,
        "absorber": absorber,
        "edge": edge,
        "first_shell": first,
        "white_line": wl,
        "checks": checks,
        "warnings": warnings,
        "n_warnings": len(warnings),
        "note": ("XANES–EXAFS 软校验只做保守极端值筛查；白线强度不能"
                 "单独证明或否定 CN，本结果不改变 cn_reliable/正式 CN 门禁"),
    }


def _find_motif(absorber: str, motif_id: str) -> structure_lib.Motif:
    for motif in structure_lib.load_library(absorber):
        if motif.id == motif_id:
            return motif
    raise XAFSConfigError(
        "FIT", f"结构库 {absorber} 中不存在 motif {motif_id!r}")


def _ss_rows(motif: structure_lib.Motif, edge: str) -> list:
    """返回 motif 的全部单散射路径摘要（已按 Reff 排序）。"""
    cache = structure_lib.ensure_feff(motif, edge=edge)
    rows = []
    for f in sorted(cache.glob("feff*.dat")):
        path = artemis.load_feff_path(f)
        if int(getattr(path, "nleg", 2)) != 2:
            continue
        geom = getattr(path, "geom", None) or []
        scatterer = geom[-1][0] if geom else "?"
        rows.append({"file": f, "reff": float(path.reff),
                     "degen": float(path.degen), "scatterer": scatterer})
    rows.sort(key=lambda r: r["reff"])
    # ---- B7-2 部分占据：按距离/元素匹配壳层平均占据度缩放简并度 ----
    # 只对确有 occ<0.999 或直接引用外部 CIF 的晶体 motif 做展开，
    # 避免全占据体系每次识别都重复展开晶胞。
    needs_occ = (motif.kind == "crystal"
                 and (motif.cif_path is not None
                      or any(float(s.get("occ", 1.0)) < 0.999
                             for s in (motif.asym_sites or []))))
    if needs_occ:
        try:
            scale = artemis.cif_shell_occupancy_scale(
                structure_lib._cif_path(motif), motif.absorber, rmax=6.0)
        except XAFSError as exc:
            raise CNError(
                "occupancy_scale_failed",
                f"{motif.id}: 无法展开 CIF 占据度；拒绝按满占据继续"
                f"（{getattr(exc, 'message', str(exc))}）") from exc
        if not scale:
            raise CNError(
                "occupancy_scale_failed",
                f"{motif.id}: CIF 占据度展开为空；拒绝按满占据继续")
        for r in rows:
            matches = [occ for (dist, el, occ) in scale
                       if el == r["scatterer"]
                       and abs(dist - float(r["reff"])) <= 0.05]
            if matches:
                occ_mean = float(np.mean(matches))
                if abs(occ_mean - 1.0) > 1e-12:
                    r["degen"] = float(r["degen"]) * occ_mean
                    r["occupancy_scale"] = occ_mean
    return rows


def _collinear_ms_rows(motif: structure_lib.Motif, edge: str,
                       rmax: float, min_abs_cos: float = 0.85,
                       min_reff_ratio: float = 1.5) -> list:
    """筛选近共线多散射（MS）路径：nleg∈{3,4}、reff≤rmax、所有内角≈0/180°。

    为什么必须几何筛选而不是把所有 nleg>2 都塞进去：FEFF 对每个团簇
    输出几十条 MS 路径（NiO rmax=4.5 时有 81 条），多数是 90° 弱路径；
    只有共线/近共线路径（如 Ni–O–Ni、Ni–O–Ni–O，reff≈2×2.0885）振幅
    足以影响 3–4 Å 拟合。

    判据用 ``|cos|`` 而不是带符号的 ``cos``：FEFF 的 ``geom`` 是路径
    原子序列，最后一条腿由末原子回到吸收原子是**隐式**的。对共线路径，
    中间原子处“指向上一原子”和“指向下一原子”的两个向量要么同向
    （继续前行，0° 转折）要么反向（180° 回散射），二者都满足 |cos|≈1；
    90° 路径给出 |cos|≈0。旧实现只看 ``max(cos)`` 且漏掉回吸收原子的
    闭合腿，会把 NiO 的 feff0007/0009 全部误删。
    """
    cache = structure_lib.ensure_feff(motif, edge=edge)
    rows = []
    for f in sorted(cache.glob("feff*.dat")):
        path = artemis.load_feff_path(f)
        nleg = int(getattr(path, "nleg", 2))
        if nleg not in (3, 4):
            continue
        reff = float(path.reff)
        if reff > rmax:
            continue
        geom = getattr(path, "geom", None) or []
        if len(geom) < nleg:
            continue
        try:
            pts = np.asarray([[float(g[4]), float(g[5]), float(g[6])]
                              for g in geom[:nleg]], dtype=float)
        except (TypeError, ValueError, IndexError):
            continue
        # 闭合路径：末原子 → 吸收原子（GEOM 不重复列吸收原子）。
        closed = np.vstack([pts, pts[0]])
        cosines = []
        for i in range(1, len(closed) - 1):
            v1 = closed[i - 1] - closed[i]
            v2 = closed[i + 1] - closed[i]
            n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
            if n1 <= 1e-8 or n2 <= 1e-8:
                cosines = []
                break
            cosines.append(float(np.dot(v1, v2) / (n1 * n2)))
        if not cosines or min(abs(c) for c in cosines) < min_abs_cos:
            continue
        first_dist = float(np.linalg.norm(pts[1] - pts[0]))
        # 聚焦判据：Reff ≈ 2·R_first（nleg=3 的三体 180° 路径与 nleg=4
        # 的往返路径）。首腿本身等于 Reff 的“长腿”路径（如 Ni–Ni–O）
        # 不是聚焦路径，不纳入受约束 MS。
        if first_dist <= 0 or reff < min_reff_ratio * first_dist:
            continue
        try:
            path_id = int(f.stem.replace("feff", ""))
        except ValueError:
            path_id = 0
        rows.append({
            "file": f, "reff": reff, "degen": float(path.degen),
            "nleg": nleg, "scatterer": str(geom[-1][0]),
            "first_el": str(geom[1][0]), "first_dist": first_dist,
            "sequence": "-".join(str(g[0]) for g in geom[:nleg]),
            "reff_ratio": float(reff / first_dist),
            "path_id": path_id,
        })
    rows.sort(key=lambda r: (r["reff"], -r["degen"]))
    return rows


def _mixture_rows(motif: structure_lib.Motif, edge: str, merge_tol: float,
                  apply_weights: bool = True):
    """展开 SS 路径；mixture 则按权重缩放简并度并收紧合并容差。

    返回 (rows, merge_tol_effective)。普通 motif 原样透传；
    mixture 的 χ 是各位点 χ 的加权线性组合，对应到路径层就是
    把每个分量的 degen 乘以该位点权重，再做**同一**窗口拟合。

    apply_weights=False 用于**权重拟合**：保留晶体学简并度、给每条
    路径打上 ``group``（组分下标）与 ``group_name``，由 fit_r_space
    在迭代中按 N_crystal·w 缩放——否则权重会被预先乘死，无法再拟合。
    """
    if motif.kind != "mixture":
        return _ss_rows(motif, edge), merge_tol
    lib = {m.id: m for m in structure_lib.load_library(motif.absorber)}
    rows = []
    for ci, comp in enumerate(motif.components):
        w = float(comp["weight"])
        for r in _ss_rows(lib[comp["id"]], edge):
            d = dict(r)
            d["group"] = ci
            d["group_name"] = comp["id"]
            if apply_weights:
                d["degen"] = float(r["degen"]) * w
            rows.append(d)
    rows.sort(key=lambda r: r["reff"])
    # 两吸收位点的近距壳（磁铁矿 Fe-O 1.886 / 2.060）只差 0.17 Å，
    # 单位点容差 0.20 会把它们错误合并；物理分辨 π/(2Δk) 允许用
    # 更紧的容差分开子壳层。
    return rows, min(merge_tol, 0.10)


def _hetero_branch_row(motif: structure_lib.Motif, edge: str,
                       element: str, distance: float) -> dict:
    """为异核共占位分支生成一条单散射 FEFF 路径（缓存到 _hetero_cache）。

    做法：把共占位位点纯化为该分支元素的单个近邻原子跑一次 FEFF，
    取该距离处的 nleg==2 路径；简并度由占位数在调用方覆盖。
    """
    element = str(element).capitalize()
    distance = float(distance)
    cache = (structure_lib.HETERO_CACHE
             / f"{motif.id}_{element}_{distance:.3f}")
    if not list(cache.glob("feff*.dat")):
        inp = artemis.build_feff_input(
            [(distance, 0.0, 0.0, element)], absorber=motif.absorber,
            rmax=distance + 0.8,
            title=f"{motif.id}_{element}_{distance:.2f}", edge=edge)
        artemis.run_feff6l(inp, cache)
    best = None
    for f in sorted(cache.glob("feff*.dat")):
        p = artemis.load_feff_path(f)
        if int(getattr(p, "nleg", 2)) != 2:
            continue
        if abs(float(p.reff) - distance) > 0.25:
            continue
        if best is None or (abs(float(p.reff) - distance)
                            < abs(float(best[0].reff) - distance)):
            best = (p, f)
    if best is None:
        raise XAFSDataError(
            "FIT", f"异核分支 {motif.absorber}-{element}@{distance:.2f} Å "
                   f"在 FEFF 缓存 {cache} 中无 SS 路径")
    return {"file": best[1], "reff": float(best[0].reff),
            "degen": float(best[0].degen), "scatterer": element}


def _hetero_shells_in_window(motif: structure_lib.Motif,
                             rmin: float, rmax: float) -> list:
    """返回窗口内需要拆分的异核共占位壳层（single_element_approx=False）。"""
    if not motif.cif_path:
        raise XAFSConfigError(
            "FIT", f"{motif.id}: hetero_cocupancy 需要 motif.cif_path "
                   f"（含 _atom_site_occupancy 的共占位 CIF）")
    hetero = artemis.cif_hetero_cocupancy_branches(
        motif.cif_path, motif.absorber, rmax=max(4.0, rmax + 1.5))
    out = []
    for si, sh in enumerate(hetero.get("shells") or []):
        if sh.get("single_element_approx"):
            continue
        if not (rmin - 0.5 <= float(sh["distance"]) <= rmax + 0.5):
            continue
        item = dict(sh)
        item["shell_index"] = si
        out.append(item)
    return out


def _replace_hetero_rows(motif: structure_lib.Motif, edge: str,
                         rows: list, hetero_shells: list) -> list:
    """用每分支 FEFF 路径替换共占位壳层的单元素近似行。"""
    out = list(rows)
    for sh in hetero_shells:
        distance = float(sh["distance"])
        els = {str(b["element"]).capitalize() for b in sh["branches"]}
        out = [r for r in out
               if not (str(r["scatterer"]).capitalize() in els
                       and abs(float(r["reff"]) - distance) <= 0.25)]
        for bi, branch in enumerate(sh["branches"]):
            row = _hetero_branch_row(
                motif, edge, str(branch["element"]).capitalize(), distance)
            row["degen"] = float(branch["degen"])
            row["hetero_group"] = int(sh["shell_index"])
            row["hetero_branch"] = int(bi)
            row["hetero_n_geom"] = int(sh["n_geom"])
            out.append(row)
    return out


def select_shell_paths(motif: structure_lib.Motif, edge: str = "K",
                       rmin: float | None = None, rmax: float | None = None,
                       merge_tol: float = 0.20, max_paths: int = 12,
                       kmin: float = 3.0, kmax: float = 12.0,
                       group_by_component: bool = False,
                       include_ms: bool = False, max_ms: int = 3,
                       isotropic_deltar: bool = False,
                       hetero_cocupancy: bool = False) -> dict:
    """按 (散射元素, Reff) 分壳层，返回代表路径 + 总简并度。

    ``merge_tol`` 默认 0.20 Å，仍小于 π/(2Δk)（Δk=6 时约 0.26 Å）——
    比 0.15 Å 更符合“不可分辨壳层”的物理判据：针铁矿 Fe-O
    1.95×3 + 2.10×3（ΔR=0.15 Å）应合并成一层 N=6，否则两条 N
    强相关、误差爆炸。代表路径取该壳层中简并度最大的一条，N 初值
    取壳层总简并度。

    自动选壳层还必须服从信息量上限：每条路径放开 N/ΔR/σ² 是 3 个
    变量，加 1 个全局 ΔE0，故 ``3n+1`` 必须小于独立点数 N_idp。
    超出时从**最远壳层**开始截断（近壳信噪比高、相关弱），并在
    ``truncated`` 里返回被丢弃的壳层，避免静默过拟合。

    kind=mixture（多位点吸收）：各分量 SS 路径按权重缩放简并度后
    统一分组，并把合并容差收紧到 ≤0.10 Å——磁铁矿 8a/16d 第一壳
    1.886 与 2.060 Å 仅差 0.17 Å，默认 0.20 容差会错误合并成一层。
    """
    if group_by_component and motif.kind != "mixture":
        raise XAFSConfigError(
            "FIT", f"{motif.id} 不是 mixture，group_by_component 仅用于多位点混合相")
    # 权重拟合必须保留**未乘权重**的晶体学简并度，否则 N·w 被预先乘死。
    rows, merge_tol = _mixture_rows(
        motif, edge, merge_tol, apply_weights=not group_by_component)
    if not rows:
        raise XAFSDataError("FIT", f"{motif.id}: 没有可用单散射路径")
    if rmin is None:
        rmin = max(1.0, min(r["reff"] for r in rows) - 0.6)
    if rmax is None:
        # 默认止步 3.5 Å：覆盖第一壳 + 常见第二壳（Fe 箔 2.87、
        # Fe–Fe 3.06），但不把 3.6–4.4 Å 的弱远壳全部放进自动拟合。
        # 需要第三壳时显式 --fit-rmax 3.9，并自行承担过拟合风险。
        rmax = min(3.5, max(r["reff"] for r in rows) + 0.3)

    # ---- B13-2 异核共占位：把单元素近似拆成每元素一条分支路径 ----
    hetero_shells = []
    if hetero_cocupancy:
        hetero_shells = _hetero_shells_in_window(
            motif, rmin=float(rmin), rmax=float(rmax))
        if hetero_shells:
            rows = _replace_hetero_rows(motif, edge, rows, hetero_shells)

    kept = [r for r in rows if rmin - 0.5 <= r["reff"] <= rmax + 0.5]
    if not kept:
        raise XAFSConfigError(
            "FIT", f"{motif.id}: R∈[{rmin:.2f},{rmax:.2f}] 内无 SS 路径")
    groups = []
    for row in kept:
        row_key = row.get("group") if group_by_component else None
        # 同一 (散射元素, 组分) 内按距离容差合并；不同组分永不合并，
        # 否则磁铁矿 tet/oct 的权重会被并进同一条路径而不可辨识。
        match = None
        for g in groups:
            if g["scatterer"] != row["scatterer"]:
                continue
            if group_by_component and g["group"] != row_key:
                continue
            if abs(row["reff"] - g["reff_mean"]) <= merge_tol:
                if (match is None
                        or abs(row["reff"] - g["reff_mean"])
                        < abs(row["reff"] - match["reff_mean"])):
                    match = g
        if match is not None:
            g = match
            old_n = g["degen"]
            new_n = old_n + row["degen"]
            g["reff_mean"] = ((g["reff_mean"] * old_n
                               + row["reff"] * row["degen"]) / new_n)
            g["degen"] = new_n
            g["members"].append(row)
            if row["degen"] > g["rep"]["degen"]:
                g["rep"] = row
            continue
        groups.append({"scatterer": row["scatterer"], "reff_mean": row["reff"],
                       "degen": row["degen"], "rep": row,
                       "members": [row],
                       "group": row_key,
                       "group_name": (row.get("group_name")
                                      if group_by_component else None)})
    # 逐行插入可能打乱距离顺序；截断策略依赖"从最远壳开始丢"。
    groups.sort(key=lambda g: g["reff_mean"])

    # 只保留中心落在拟合窗口内的壳层；窗口外的弱路径不参与自动拟合。
    groups = [g for g in groups
              if rmin - 0.20 <= g["reff_mean"] <= rmax + 0.15]
    if not groups:
        raise XAFSConfigError(
            "FIT", f"{motif.id}: R∈[{rmin:.2f},{rmax:.2f}] 内无有效壳层")

    n_idp = artemis.independent_points(kmin, kmax, rmin, rmax)
    if isotropic_deltar:
        # ΔR 已压缩为唯一全局 α：每壳只剩 N、σ² 两个变量，加 α、ΔE₀
        # 两个全局变量，信息量上限为 2n + 2 < N_idp。
        max_by_idp = max(1, int((n_idp - 2.0) // 2.0))
    else:
        max_by_idp = max(1, int((n_idp - 1.0) // 3.0))
    limit = min(max_paths, max_by_idp)
    truncated = groups[limit:]
    groups = groups[:limit]

    files, degens, labels, shells = [], [], [], []
    for g in groups:
        el = g["scatterer"]
        # 用代表距离做唯一标签，避免 "Fe-O2" 让人误以为是第二壳；
        # 多壳层表里真正重要的是 Reff/N/R，而不是序号。
        base = f"{motif.absorber}-{el}@{g['reff_mean']:.2f}"
        label = (f"{g['group_name']}:{base}" if group_by_component else base)
        files.append(g["rep"]["file"])
        degens.append(g["degen"])
        labels.append(label)
        shells.append({"label": label, "scatterer": el,
                       "reff": round(g["rep"]["reff"], 4),
                       "reff_mean": round(g["reff_mean"], 4),
                       "degen": g["degen"], "n_members": len(g["members"]),
                       "group": (g["group_name"] if group_by_component
                                 else None)})

    # ---- 权重拟合所需的组分映射 ----
    # fit_r_space 要求 path_groups 是 0..K-1 连续下标，而组分可能因
    # R 窗/信息量截断缺失，故重映射并同时给出原始组分下标。
    path_groups = None
    group_names = None
    component_ids = None
    weight_degenerate_pairs = []
    if group_by_component:
        used = sorted({int(g["group"]) for g in groups})
        remap = {old: i for i, old in enumerate(used)}
        path_groups = [remap[int(g["group"])] for g in groups]
        component_ids = used
        group_names = []
        for old in used:
            name = next(g["group_name"] for g in groups
                        if int(g["group"]) == old)
            group_names.append(name)
        # 不同组分在同一 (散射元素, 距离) 上出现时，只有振幅之和可辨识，
        # 单组分权重严格简并——必须显式标记，不能给出"相分数"。
        for a in range(len(groups)):
            for b in range(a + 1, len(groups)):
                ga, gb = groups[a], groups[b]
                if (int(ga["group"]) != int(gb["group"])
                        and ga["scatterer"] == gb["scatterer"]
                        and abs(ga["reff_mean"] - gb["reff_mean"])
                        <= merge_tol):
                    weight_degenerate_pairs.append(
                        [f"{ga['group_name']}@{ga['reff_mean']:.2f}",
                         f"{gb['group_name']}@{gb['reff_mean']:.2f}"])

    # ---- 受约束多散射：追加近共线 nleg>2 路径，不新增自由变量 ----
    ms_links = []
    if include_ms:
        for row in _collinear_ms_rows(motif, edge, rmax + 0.15):
            parent = None
            for gi, g in enumerate(groups):
                if g["scatterer"] != row["first_el"]:
                    continue
                if abs(g["reff_mean"] - row["first_dist"]) <= 0.05:
                    parent = gi
                    break
            if parent is None:
                continue
            idx = len(files)
            # 标签必须带完整原子序列：同一 Reff 下 Ni–O–O 与
            # Ni–O–Ni–O 是不同路径，只写末原子会产生重名。
            label = (f"{row.get('sequence', motif.absorber + '-' + row['scatterer'])}"
                     f"@{row['reff']:.2f}(MS{row['nleg']}"
                     f"#{row.get('path_id', 0)})")
            files.append(row["file"])
            degens.append(row["degen"])
            labels.append(label)
            shells.append({"label": label, "scatterer": row["scatterer"],
                           "reff": round(row["reff"], 4),
                           "reff_mean": round(row["reff"], 4),
                           "degen": row["degen"], "n_members": 1,
                           "group": None, "kind": "ms",
                           "parent": labels[parent]})
            ms_links.append({
                "ms_path_index": idx,
                "parent_shell_index": parent,
                "reff_ratio": float(row["reff"]
                                    / groups[parent]["reff_mean"]),
                "degen_ratio": float(row["degen"]
                                     / groups[parent]["degen"]),
                "nleg": int(row["nleg"]),
            })
            if len(ms_links) >= max_ms:
                break

    # ---- B13-2 双分支绑定：同壳层不同元素分支共享 ΔR/σ²，各自放开 N ----
    path_links = []
    if hetero_cocupancy and hetero_shells:
        by_group: dict[int, list[int]] = {}
        for gi, g in enumerate(groups):
            hg = g["rep"].get("hetero_group")
            if hg is None:
                continue
            by_group.setdefault(int(hg), []).append(gi)
        for _hg, idxs in sorted(by_group.items()):
            idxs.sort(key=lambda gi: int(
                groups[gi]["rep"].get("hetero_branch", 0)))
            src = idxs[0]
            for tgt in idxs[1:]:
                path_links.append({"target": int(tgt), "source": int(src),
                                   "share_deltar": True,
                                   "share_sigma2": True,
                                   "degen_ratio": None})

    return {"files": files, "degens": degens, "labels": labels,
            "shells": shells, "rmin": float(rmin), "rmax": float(rmax),
            "merge_tol": merge_tol, "n_independent": round(n_idp, 2),
            "max_paths_by_idp": max_by_idp, "ms_links": ms_links,
            "groups": path_groups, "group_names": group_names,
            "component_ids": component_ids,
            "weight_degenerate_pairs": weight_degenerate_pairs,
            "hetero_shells": [{"distance": float(sh["distance"]),
                               "n_geom": int(sh["n_geom"]),
                               "branches": [dict(b) for b in sh["branches"]]}
                              for sh in hetero_shells],
            "path_links": path_links,
            "truncated": [{"label": (f"{g['group_name']}:{motif.absorber}-"
                                     f"{g['scatterer']}@{g['reff_mean']:.2f}"
                                     if group_by_component else
                                     f"{motif.absorber}-{g['scatterer']}"
                                     f"@{g['reff_mean']:.2f}"),
                           "reff_mean": round(g["reff_mean"], 4),
                           "degen": g["degen"]} for g in truncated]}


def _normalize_c3_paths(vary_c3, n_path: int):
    """把 CLI/API 的 vary_c3 归一成 artemis 接受的 None / True / [idx...]。"""
    if vary_c3 is None or vary_c3 is False:
        return None
    if vary_c3 is True:
        return True
    if isinstance(vary_c3, str):
        key = vary_c3.strip().lower()
        if key in ("", "none", "false", "off"):
            return None
        if key == "first":
            return [0]
        if key == "all":
            return True
        raise XAFSConfigError(
            "FIT", f"vary_c3 字符串只支持 'first'/'all'，收到 {vary_c3!r}")
    if isinstance(vary_c3, (list, tuple, set)):
        idx = sorted({int(i) for i in vary_c3})
        for i in idx:
            if i < 0 or i >= n_path:
                raise XAFSConfigError(
                    "FIT", f"vary_c3 路径下标 {i} 超出范围 0–{n_path - 1}")
        return idx or None
    raise XAFSConfigError(
        "FIT", f"vary_c3 必须是 bool 或路径下标列表，收到 {vary_c3!r}")


def rdf_table(fit) -> dict:
    """B9-1：把拟合出的 C3/C4 累积量重构为非高斯 g(R)。

    只对真正放开了 C3 或 C4 的路径生成——常规拟合没有非高斯信息，
    不能凭空画一条不对称分布。返回:
        {"rows": [{"label", "r", "g_r"}, ...],
         "summary": [{"label", "r_mode", "r_mean", "skewness", "area",
                      "c3", "c4"}, ...]}
    """
    rows, summary = [], []
    for p in fit.paths:
        if p.c3 is None and p.c4 is None:
            continue
        r_grid = np.linspace(max(0.5, p.r - 0.8), p.r + 1.0, 401)
        rdf = artemis.reconstruct_pair_distribution(
            r_grid, n=float(p.degeneracy), r0=float(p.r),
            sigma2=float(p.sigma2), c3=float(p.c3 or 0.0),
            c4=float(p.c4 or 0.0))
        rows.append({"label": p.label, "r": rdf["r"], "g_r": rdf["g_r"]})
        summary.append({"label": p.label, "r_mode": rdf["r_mode"],
                        "r_mean": rdf["r_mean"],
                        "skewness": rdf["skewness"],
                        "area": rdf["area"], "c3": p.c3, "c4": p.c4})
    return {"rows": rows, "summary": summary}


def write_rdf_csv(path, rdf: dict):
    """把 rdf_table 的结果写成 long-form CSV（每路径一段）。"""
    import csv
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "r_A", "g_r", "r_mode", "r_mean",
                    "skewness", "area"])
        meta = {m["label"]: m for m in rdf.get("summary", [])}
        for block in rdf.get("rows", []):
            m = meta.get(block["label"], {})
            for riv, g in zip(block["r"], block["g_r"]):
                w.writerow([block["label"], f"{riv:.6f}", f"{g:.6e}",
                            m.get("r_mode"), m.get("r_mean"),
                            m.get("skewness"), m.get("area")])
    return path


def _s02_calibrated_flag(absorber, s02_source: str | None) -> bool:
    """正式表 S₀² 是否来自“已标定”来源。

    静态表已标定（absorber.calibrated=True）或 T5 独立标定产物
    （``s02_calibration:`` 前缀）都算已标定；普通 ``--s02`` 手填值
    仍为 False。这不改变 #60 的静态表约束。
    """
    return (bool(getattr(absorber, "calibrated", False))
            or str(s02_source or "").startswith("s02_calibration:"))


def _formal_cn_meta(formal_cn: bool, structure_model: StructureModel | None,
                    status_override: str | None = None,
                    reason_override: str | None = None) -> dict:
    """给拟合结果附加诊断/正式双轨标记。

    P0a 只负责“拒绝或挂起”，绝不产生 ``n_fit`` 点值：
    - ``formal_cn=False``：仅标记 diagnostic_only；
    - 无模型：``no_cif``；
    - L1：``l1_diagnostic_only``；
    - 模型来源合格：``pending_consistency``；P3 条件 CN 成功后由
      ``build_formal_cn_block`` 覆盖为 ``ok``/``consistency_failed``；
    - 模型来源不合格：返回 ``validate_formal`` 的拒绝码。
    """
    if not formal_cn:
        return {"diagnostic_only": True, "cn_status": "diagnostic_only"}
    if status_override:
        status = status_override
        reason = reason_override
    elif structure_model is None:
        status = "no_cif"
        reason = "未提供 StructureModel（无合格结构来源）"
    else:
        verdict = validate_formal(structure_model)
        if verdict["ok"]:
            status = "pending_consistency"
            reason = "P0a 尚未实现固定 N 一致性检验；正式 CN 从 P3 开始"
        else:
            status = verdict["code"] or "unverified_source"
            reason = verdict["reason"]
    return {
        "diagnostic_only": True,
        "cn_status": status,
        "formal_cn": {
            "status": status,
            "rows": [],
            "reject_reason": reason,
            "model": (structure_model.to_json()
                      if structure_model is not None else None),
        },
    }


def fit_identified(res, motif_id: str | None = None,
                   absorber_symbol: str | None = None,
                   rmin: float | None = None, rmax: float | None = None,
                   allow_ambiguous: bool = False,
                   l1_ligand: str = "O", l1_geometry: str = "oct",
                   fit_weights: bool = False,
                   shared_sigma2: str | list | None = None,
                   kweights: tuple | list | None = None,
                   include_ms: bool = False,
                   vary_c3: bool | list | str = False,
                   vary_c4: bool | list | str = False,
                   isotropic_deltar: bool = False,
                   fit_space: str = "r",
                   bvs_restraint: dict | None = None,
                   catalytic_model=None,
                   hetero_cocupancy: bool = False,
                   constraints=None,
                   s02: float | None = None,
                   s02_source: str | None = None,
                   verbose: bool = True,
                   formal_cn: bool = False,
                   structure_model: StructureModel | None = None,
                   gate_result=None) -> dict:
    """对识别结果做 S7 拟合；未 confident 时返回 L1 第一壳结果。

    ``fit_weights=True`` 仅对 mixture（多位点吸收）生效：固定每条路径的
    晶体学简并度，放开各组分权重（softmax，和=1）与每壳 ΔR/σ²/ΔE₀，
    给出相分数。此时 N 列是**有效简并度 N_crystal·w**，不是配位数。

    ``shared_sigma2`` 对标 Artemis 的路径约束：多条路径共享一个 σ²
    变量，减少自由参数、避免 σ² 与 ΔR 过度相关。
        None        每条路径独立 σ²（默认）
        "all"       所有路径共享一个 σ²
        [[0,1],[2]] 显式分组：路径 0、1 共享 σ²，路径 2 独立
    与 fit_weights 互斥（权重模式下 σ² 已按组分缩放，不再共享）。

    ``kweights``：多 k 权重联合拟合（对标 Artemis），如 ``(1, 2, 3)``。
    None 时沿用单权 ``res.ft.kweight``（默认 2），行为与历史完全一致。

    ``fit_space``（B11-3）：``"r"``（默认）/ ``"q"`` / ``"wt"``，透传给
    ``artemis.fit_r_space``。不同空间的 Rf 口径不同，不可横向比较。
    ``bvs_restraint``（B11-2）：``{"target", "sigma", "r0_table"?}`` 软约束
    键价和；L1 合成团簇不提供 R₀ 语义，遇到 L1 降级时显式拒绝而不是
    静默忽略。
    """
    absorber_symbol = absorber_symbol or getattr(res, "absorber_symbol", None)
    if not absorber_symbol:
        raise XAFSConfigError(
            "FIT", "缺少吸收元素：给 absorber_symbol 或 --absorber")
    absorber = (get_absorber(absorber_symbol) if s02 is None
                else structure_id.absorber_with_s02(absorber_symbol, s02))
    if absorber.s02 is None:
        raise XAFSConfigError(
            "FIT", f"{absorber.symbol} 未标定 S₀²，拒绝 CN 定量；"
                   f"先在 absorber.py 用标准样标定")
    s02_label = (s02_source if s02_source is not None else
                 ("--s02 runtime override" if s02 is not None
                  else "absorber.py"))
    if str(s02_label).startswith("s02_calibration:") and s02 is None:
        raise XAFSConfigError(
            "FIT", "T5 S₀² 来源已给出但 s02=None；拒绝在没有数值的情况下"
                   "记录标定溯源")
    self_abs = _self_absorption_diag(res)

    # ---- B13-1 显式催化模型：绕过晶相库搜索，直接进 S7 ----
    catalytic_motif = None
    if catalytic_model is not None:
        catalytic_motif = structure_lib.parse_catalytic_cluster_spec(
            absorber.symbol, catalytic_model)
        if str(catalytic_motif.absorber).capitalize() != absorber.symbol.capitalize():
            raise XAFSConfigError(
                "FIT", f"催化模型的吸收元素 {catalytic_motif.absorber!r} 与 "
                       f"--absorber {absorber.symbol!r} 不一致")

    sid = getattr(res, "structure", None)
    if catalytic_motif is not None:
        if sid is None:
            sid = {
                "status": "catalytic_model",
                "absorber": absorber.symbol,
                "ranked": [{"id": catalytic_motif.id,
                            "score": None, "de0": 0.0,
                            "kind": catalytic_motif.kind,
                            "nuclearity": catalytic_motif.nuclearity}],
                "indistinguishable": [],
                "uncovered_peaks": [],
                "coverage_ok": True,
                "confidence": {
                    "margin": None,
                    "note": "用户显式给定催化模型，未做结构库搜索"},
                "warning": ("使用 --catalytic-model 显式模型；识别状态"
                            "不代表库内排序结果，长程结构由该模型强制假设"),
            }
        if motif_id is None:
            motif_id = catalytic_motif.id
    else:
        if sid is None:
            sid = structure_id.identify_result(
                res, absorber_symbol=absorber_symbol,
                l1_ligand=l1_ligand, l1_geometry=l1_geometry,
                vary_c3=vary_c3, fit_space=fit_space, s02=s02,
                verbose=verbose)

    if catalytic_model is None and motif_id is None:
        status = sid.get("status")
        if status == "confident" or allow_ambiguous:
            motif_id = sid["ranked"][0]["id"]
        else:
            l1 = sid.get("l1_first_shell")
            if l1 is None:
                detail = sid.get("l1_error", "L1 拟合未生成")
                raise XAFSDataError(
                    "FIT", f"识别状态 {status}，不能自动全拟合：{detail}")
            l1_c3_bound = bool(l1.get("c3_at_bound", False))
            if bvs_restraint is not None:
                raise XAFSConfigError(
                    "FIT", "bvs_restraint 需要结构模型给出的配体与参考键长；"
                           "识别状态不是 confident，L1 合成团簇的 BVS 无物理"
                           "意义，拒绝静默忽略该约束")
            l1_reliable = (l1.get("rf", 1.0) <= 0.08
                           and abs(l1.get("delta_e0", 0.0)) < 15.0
                           and not self_abs["cn_gate"]
                           and not l1_c3_bound)
            l1_flags = {"sigma2_at_bound":
                        float(l1.get("sigma2", 1.0)) <= 5e-4,
                        "c3_at_bound": l1_c3_bound}
            out = {"mode": "l1_first_shell", "status": status,
                   "coverage_ok": sid.get("coverage_ok"),
                   "cn_reliable": l1_reliable,
                   "s02": absorber.s02, "s02_source": s02_label,
                   "flags": {"c3_dropped_by_nidp":
                             bool(l1.get("c3_dropped_by_nidp", False)),
                             "c3_at_bound": l1_c3_bound},
                   "structure": sid, "l1": l1,
                   "self_absorption": self_abs,
                   "cn_systematics": systematic_budget(
                       None, flags=l1_flags, self_abs=self_abs),
                   "warning": l1.get(
                       "warning", "仅第一壳可信；长程结构未知")}
            if self_abs["risk"] in ("high", "unknown",
                                    "corrected_unreliable"):
                prefix = ("荧光自吸收校正不可靠" if self_abs["corrected"]
                          else "荧光自吸收未校正")
                out["warning"] += ("；" + prefix + "（"
                                   + self_abs["note"] + "）")
            elif self_abs["risk"] == "possible":
                out["warning"] += ("；荧光自吸收风险（"
                                   + self_abs["note"] + "）")
            if fit_weights:
                out["warning"] += ("；识别状态不是 confident，"
                                   "--fit-weights 未执行（返回 L1 第一壳）")
            out.update(_formal_cn_meta(
                formal_cn, structure_model,
                status_override=("l1_diagnostic_only" if formal_cn else None),
                reason_override=("L1 合成团簇只做诊断，不进入正式 CN 表"
                                 if formal_cn else None)))
            return out

    motif = (catalytic_motif if catalytic_motif is not None
             else _find_motif(absorber.symbol, motif_id))
    if fit_weights and motif.kind != "mixture":
        raise XAFSConfigError(
            "FIT", f"{motif.id} 不是混合相（kind={motif.kind}），"
                   f"--fit-weights 只用于多位点吸收的 mixture")
    plan = select_shell_paths(motif, edge=absorber.edge,
                              rmin=rmin, rmax=rmax,
                              kmin=res.ft.kmin, kmax=res.ft.kmax,
                              group_by_component=fit_weights,
                              include_ms=include_ms,
                              isotropic_deltar=isotropic_deltar,
                              hetero_cocupancy=hetero_cocupancy)

    # ---- P3 条件 CN：固定 N → 自由 N（首壳单一比例因子）→ Hamilton F-test ----
    conditional_cn = None
    conditional_cn_error = ""
    n0_mapping = None
    if formal_cn and gate_result is not None and gate_result.ok:
        try:
            n0_mapping = map_n0_to_plan(plan, gate_result.n0_shells,
                                        kmin=res.ft.kmin, kmax=res.ft.kmax)
            if n0_mapping["mapping_ok"]:
                conditional_cn = run_conditional_cn(
                    res, plan, n0_mapping,
                    absorber_symbol=absorber.symbol, s02=absorber.s02,
                    kmin=res.ft.kmin, kmax=res.ft.kmax,
                    rmin=plan["rmin"], rmax=plan["rmax"],
                    verbose=verbose)
            else:
                conditional_cn_error = (
                    "N₀ 映射失败: " + n0_mapping["mapping_reason"])
            if verbose and conditional_cn is not None:
                print(f"  [P3] 条件 CN: N/N₀={conditional_cn['n_ratio']} "
                      f"({conditional_cn['n_status']}), "
                      f"consistency_passed={conditional_cn['consistency_passed']}")
        except (XAFSConfigError, XAFSDataError) as exc:
            conditional_cn_error = (
                "条件 CN 跳过: " + getattr(exc, "message", str(exc)))
            if verbose:
                print(f"  [P3] {conditional_cn_error}")
    elif formal_cn and structure_model is not None:
        conditional_cn_error = "缺少 CIF gate 结果，未运行条件 CN"
        if verbose:
            print(f"  [P3] {conditional_cn_error}")

    if verbose:
        print(f"  [S7] motif={motif.id}  R∈[{plan['rmin']:.2f},{plan['rmax']:.2f}] Å"
              f"  S₀²={absorber.s02}  N_idp={plan['n_independent']:.1f}")
        for sh in plan["shells"]:
            print(f"       {sh['label']:<10} Reff={sh['reff']:.3f} Å"
                  f"  N0={sh['degen']:.2f}  ({sh['n_members']} 条合并)")
        if plan["truncated"]:
            names = [t["label"] for t in plan["truncated"]]
            print(f"       ⚠ 信息量上限：最多 {plan['max_paths_by_idp']} 条壳层"
                  f"（3n+1 < N_idp），已截断 {names}")

    weight_init = None
    if fit_weights:
        all_ids = list(range(len(motif.components)))
        missing = [i for i in all_ids if i not in plan["component_ids"]]
        if missing:
            names = [motif.components[i]["id"] for i in missing]
            raise XAFSConfigError(
                "FIT", "权重拟合缺少组分路径：" + ", ".join(names)
                       + "；这些组分在 R 窗内没有可用壳层或被 N_idp 截断，"
                         "相分数不可辨识。请扩大 --fit-rmax 或减小壳层数")
        # 初始权重按识别阶段给出的组分比例；仅用于优化起点。
        weight_init = [float(motif.components[i]["weight"])
                       for i in plan["component_ids"]]
        total_w = sum(weight_init)
        if total_w <= 0:
            raise XAFSConfigError("FIT", "mixture 初始权重之和必须为正")
        weight_init = [w / total_w for w in weight_init]

    # ---- σ² 共享组（对标 Artemis 路径约束）----
    sigma2_groups = None
    sigma2_group_names = None
    if shared_sigma2 is not None and not fit_weights:
        n_path = len(plan["files"])
        if shared_sigma2 == "all":
            sigma2_groups = [list(range(n_path))]
            sigma2_group_names = ["all_shells"]
        elif isinstance(shared_sigma2, list):
            # 显式分组：[[0, 1], [2]] 验证下标合法
            for g in shared_sigma2:
                for pi in g:
                    if int(pi) < 0 or int(pi) >= n_path:
                        raise XAFSConfigError(
                            "FIT", f"shared_sigma2 路径下标 {pi} 超出范围 0–{n_path-1}")
            sigma2_groups = [list(map(int, g)) for g in shared_sigma2]
            sigma2_group_names = [f"s2grp{g}" for g in range(len(sigma2_groups))]
        else:
            raise XAFSConfigError(
                "FIT", f"shared_sigma2 必须是 'all' 或路径下标列表，收到 {shared_sigma2!r}")

    requested_c3 = _normalize_c3_paths(vary_c3, len(plan["files"]))
    requested_c4 = _normalize_c3_paths(vary_c4, len(plan["files"]))
    if fit_weights and requested_c3:
        raise XAFSConfigError(
            "FIT", "C3 与 mixture 相分数拟合不能同时开启："
                   "权重模式的有效简并度 N_crystal·w 不是配位数，"
                   "C3 无法解释为近邻分布偏斜")
    if isotropic_deltar and fit_weights:
        raise XAFSConfigError(
            "FIT", "isotropic_deltar 与 mixture 相分数拟合不能同时开启："
                   "各组分晶格膨胀不一致时全局 α 无物理意义")

    def _fit_with_cumulants(c3_arg, c4_arg):
        return artemis.fit_r_space(
            res.chi_k, res.ft, plan["files"],
            rmin=plan["rmin"], rmax=plan["rmax"],
            kmin=res.ft.kmin, kmax=res.ft.kmax,
            kweight=(tuple(kweights) if kweights else res.ft.kweight),
            s02=absorber.s02, vary_delta_e0=True, vary_n=not fit_weights,
            vary_weights=fit_weights,
            vary_c3=(c3_arg if c3_arg is not None else False),
            vary_c4=(c4_arg if c4_arg is not None else False),
            isotropic_deltar=isotropic_deltar,
            path_degen=plan["degens"], path_labels=plan["labels"],
            path_groups=(plan["groups"] if fit_weights else None),
            path_group_names=(plan["group_names"] if fit_weights else None),
            weight_init=weight_init,
            sigma2_groups=sigma2_groups,
            sigma2_group_names=sigma2_group_names,
            ms_links=(plan["ms_links"] or None),
            path_links=(plan.get("path_links") or None),
            constraints=constraints,
            fit_space=fit_space,
            bvs_restraint=bvs_restraint)

    # 自动 N_idp 防护（B9-1：C3 与 C4 完全对称）：
    # 先按用户请求放开；若 n_var ≥ N_idp，降级到仅第一壳；仍超限则
    # 彻底关闭，并分别记录 c3_dropped_by_nidp / c4_dropped_by_nidp。
    # C4 是更高阶项，优先降 C4 以保留用户更常需要的 C3。
    c3_effective = requested_c3
    c4_effective = requested_c4
    fit = _fit_with_cumulants(c3_effective, c4_effective)
    c3_dropped_by_nidp = False
    c4_dropped_by_nidp = False

    if c4_effective and fit.n_variables >= fit.n_independent:
        c4_dropped_by_nidp = True
        if c4_effective is True or tuple(c4_effective) != (0,):
            c4_effective = [0]
            fit = _fit_with_cumulants(c3_effective, c4_effective)
        if fit.n_variables >= fit.n_independent:
            c4_effective = None
            fit = _fit_with_cumulants(c3_effective, c4_effective)

    if c3_effective and fit.n_variables >= fit.n_independent:
        c3_dropped_by_nidp = True
        if c3_effective is True or tuple(c3_effective) != (0,):
            c3_effective = [0]
            fit = _fit_with_cumulants(c3_effective, c4_effective)
        if fit.n_variables >= fit.n_independent:
            c3_effective = None
            fit = _fit_with_cumulants(c3_effective, c4_effective)

    # ---- B13-1 DAC 金属键存在性：含/不含 M–M 的 Hamilton F-test ----
    # 只看 Rf 无法回答“到底有没有 M–M 键”：多一条路径总会降低残差。
    # 必须在同一 k/R 窗内比较两个嵌套模型，并把 N_MM 与其形式误差
    # 一起纳入判据（N 不显著大于误差时不算检出金属键）。
    mm_bond_ftest = None
    if (catalytic_motif is not None
            and str(catalytic_motif.nuclearity) == "binuclear"
            and not fit_weights):
        partner = None
        for sh in plan["shells"]:
            el = str(sh.get("scatterer") or "")
            if (el and el != absorber.symbol
                    and el not in structure_lib._LIGAND_ELEMENTS):
                partner = el
                break
        mm_idx = [i for i, sh in enumerate(plan["shells"])
                  if partner is not None and sh["scatterer"] == partner]
        mm_set = set(mm_idx)
        keep = [i for i in range(len(plan["files"])) if i not in mm_set]
        if partner is not None and mm_idx and keep:
            try:
                simple_fit = artemis.fit_r_space(
                    res.chi_k, res.ft,
                    [plan["files"][i] for i in keep],
                    rmin=plan["rmin"], rmax=plan["rmax"],
                    kmin=res.ft.kmin, kmax=res.ft.kmax,
                    kweight=(tuple(kweights) if kweights
                             else res.ft.kweight),
                    s02=absorber.s02, vary_delta_e0=True, vary_n=True,
                    path_degen=[plan["degens"][i] for i in keep],
                    path_labels=[plan["labels"][i] for i in keep],
                    fit_space=fit_space)
                ftest = artemis.shell_f_test(simple_fit, fit)
                mm_path = fit.paths[mm_idx[0]]
                conf = ftest.get("confidence")
                p_value = (None if conf is None
                           else float(max(0.0, 1.0 - float(conf))))
                n_err = mm_path.n_err
                significant = bool(
                    p_value is not None and p_value < 0.05
                    and n_err is not None
                    and float(mm_path.degeneracy) > float(n_err))
                mm_bond_ftest = {
                    "partner": partner,
                    "f_stat": (None if ftest.get("F") is None
                               else float(ftest["F"])),
                    "p_value": p_value,
                    "significant": significant,
                    "delta_r_factor": float(simple_fit.r_factor
                                            - fit.r_factor),
                    "n_mm": float(mm_path.degeneracy),
                    "n_mm_err": (None if n_err is None else float(n_err)),
                    "r_mm": float(mm_path.r),
                    "dfn": int(ftest.get("dfn", 0)),
                    "dfd": float(ftest.get("dfd", 0.0)),
                    "method": ftest.get("method"),
                    "simple_r_factor": float(simple_fit.r_factor),
                }
            except XAFSError as exc:
                mm_bond_ftest = {
                    "partner": partner, "f_stat": None, "p_value": None,
                    "significant": False, "delta_r_factor": None,
                    "error": getattr(exc, "message", str(exc))}

    # 形式误差 >1 或 >50% 都说明 N 没被数据约束住（Fe 箔第二壳
    # N=9.45±2.40 是典型：Rf 很好但 N 明显偏高，不能当可靠 CN）。
    n_err_bad = any(
        p.n_err is not None and (p.n_err > 1.0
                                 or p.n_err > 0.5 * p.degeneracy)
        for p in fit.paths)
    # σ² 撞下界不是"振动很小"，而是模型形状与数据失配时 σ² 被反向
    # 压缩的证据（goethite 中距离壳 0.0001–0.0007 Å²）。此时 N 由
    # 振幅补偿撑起来，形式误差可能仍好看，必须独立判不可靠。
    sigma2_bound_bad = any(p.sigma2_at_bound for p in fit.paths)
    c3_bound_bad = any(p.c3_at_bound for p in fit.paths)
    c4_bound_bad = any(p.c4_at_bound for p in fit.paths)
    c3_active_shells = [p.label for p in fit.paths if p.c3 is not None]
    coverage_ok = bool(sid.get("coverage_ok", True))
    base_reliable = not (fit.overfit or fit.r_factor > 0.08
                         or abs(fit.delta_e0) >= 15.0 or n_err_bad
                         or sigma2_bound_bad or c3_bound_bad or c4_bound_bad
                         or not coverage_ok
                         or self_abs["cn_gate"])

    # mixture 权重拟合的可靠性是**相分数**的可靠性，和配位数不是一回事：
    # 权重撞 softmax 边界、误差过大、或两组分在同一壳层不可分辨，都会
    # 让相分数失去物理含义。
    weight_at_bound = False
    weight_err_bad = False
    if fit_weights and fit.mixture_weights:
        weight_at_bound = any(w["weight"] <= 0.01 or w["weight"] >= 0.99
                              for w in fit.mixture_weights)
        weight_err_bad = any(
            (not math.isfinite(w["weight_err"]))
            or w["weight_err"] > 0.15
            or w["weight_err"] > 0.5 * w["weight"]
            for w in fit.mixture_weights)
    weight_degenerate = bool(plan["weight_degenerate_pairs"])
    weights_reliable = None
    if fit_weights:
        weights_reliable = (base_reliable and not weight_at_bound
                            and not weight_err_bad and not weight_degenerate)
    # 权重模式下 N 列是有效简并度 N_crystal·w，不是配位数，故 CN 不成立。
    cn_reliable = base_reliable and not fit_weights
    flags = {
        "overfit": bool(fit.overfit),
        "r_factor_bad": bool(fit.r_factor > 0.08),
        "delta_e0_bad": bool(abs(fit.delta_e0) >= 15.0),
        "n_err_bad": bool(n_err_bad),
        "sigma2_at_bound": bool(sigma2_bound_bad),
        "coverage_bad": bool(not coverage_ok),
        # 撞界的壳层标签，报告层可直接列出"哪些壳不可信"
        "sigma2_bound_shells": [p.label for p in fit.paths
                                if p.sigma2_at_bound],
        "weight_at_bound": bool(weight_at_bound),
        "weight_err_bad": bool(weight_err_bad),
        "weight_degenerate": weight_degenerate,
        "weight_degenerate_pairs": plan["weight_degenerate_pairs"],
        "self_absorption_gate": bool(self_abs["cn_gate"]),
        "self_absorption_risk": self_abs["risk"],
        "c3_requested": bool(requested_c3),
        "c3_dropped_by_nidp": bool(c3_dropped_by_nidp),
        "c3_at_bound": bool(c3_bound_bad),
        "c3_active_shells": c3_active_shells,
        "c3_bound_shells": [p.label for p in fit.paths if p.c3_at_bound],
        "c4_requested": bool(requested_c4),
        "c4_dropped_by_nidp": bool(c4_dropped_by_nidp),
        "c4_at_bound": bool(c4_bound_bad),
        "c4_active_shells": [p.label for p in fit.paths if p.c4 is not None],
        "c4_bound_shells": [p.label for p in fit.paths if p.c4_at_bound],
    }

    # ---- B17 催化剂精细结构反演：纳米粒径 / 合金短程有序 ----
    nanoparticle_size = None
    alloy_sro = None
    metal_shells = [sh for sh in plan["shells"] if sh.get("scatterer")]
    if (not fit_weights and metal_shells
            and len(metal_shells) == len(plan["shells"])
            and all(sh["scatterer"] == absorber.symbol
                    for sh in metal_shells)):
        lat = ("bcc" if motif.id in ("fe_bcc", "mo_bcc", "cr_bcc")
               else "hcp" if motif.id in ("co_hcp", "zn_hcp") else "fcc")
        n1_bulk = {"fcc": 12.0, "hcp": 12.0, "bcc": 8.0}[lat]
        if float(fit.paths[0].degeneracy) < 0.96 * n1_bulk:
            n2_val = None
            if (len(fit.paths) > 1
                    and plan["shells"][1]["scatterer"] == absorber.symbol):
                n2_val = float(fit.paths[1].degeneracy)
            est = artemis.estimate_nanoparticle_size(
                n1=float(fit.paths[0].degeneracy),
                r1=float(fit.paths[0].r),
                n1_err=(None if fit.paths[0].n_err is None
                        else float(fit.paths[0].n_err)),
                n2=n2_val, lattice=lat)
            if est.get("regime") != "bulk_limit":
                nanoparticle_size = est
    if (str(getattr(motif, "nuclearity", "")) == "alloy"
            and not fit_weights):
        partner = str((motif.source or {}).get("partner") or "")
        aa_paths, ab_paths = [], []
        for sh, p in zip(plan["shells"], fit.paths):
            if sh["scatterer"] == absorber.symbol:
                aa_paths.append(p)
            elif partner and sh["scatterer"] == partner:
                ab_paths.append(p)
        if partner and (aa_paths or ab_paths):
            n_aa = float(sum(p.degeneracy for p in aa_paths))
            n_ab = float(sum(p.degeneracy for p in ab_paths))

            def _quad_err(paths):
                errs = [p.n_err for p in paths if p.n_err is not None]
                if len(errs) != len(paths):
                    return None
                return float(np.sqrt(sum(float(e) ** 2 for e in errs)))

            alloy_sro = artemis.compute_cowley_sro(
                n_aa=n_aa, n_ab=n_ab,
                x_partner=float((motif.source or {}).get("x_partner", 0.5)),
                n_aa_err=_quad_err(aa_paths),
                n_ab_err=_quad_err(ab_paths))
            n_total = n_aa + n_ab
            if n_total < 11.5:
                r_vals = [(float(p.degeneracy), float(p.r))
                          for p in aa_paths + ab_paths]
                r_mean = (sum(n * r for n, r in r_vals) / n_total
                          if n_total > 0 else None)
                if r_mean is not None:
                    est = artemis.estimate_nanoparticle_size(
                        n1=n_total, r1=r_mean,
                        n1_err=_quad_err(aa_paths + ab_paths),
                        lattice=str((motif.source or {}).get(
                            "lattice", "fcc")))
                    if est.get("regime") != "bulk_limit":
                        nanoparticle_size = est

    # ---- B18-1 Rayleigh 分辨极限与轻元素 Z 简并门禁 ----
    res_shells = []
    for pi, p in enumerate(fit.paths):
        sh = plan["shells"][pi] if pi < len(plan["shells"]) else {}
        res_shells.append({"label": p.label,
                           "scatterer": sh.get("scatterer", ""),
                           "r": float(p.r)})
    resolution_check = artemis.check_shell_resolution_and_z_degeneracy(
        res_shells, kmin=res.ft.kmin, kmax=res.ft.kmax)
    flags["sub_rayleigh_z_degeneracy"] = bool(
        not resolution_check["resolvable"])

    # ---- B4-1 增减壳层 F-test：同一 R 窗，仅保留第 1 壳的对照拟合 ----
    # 必须用**同一数据窗口**，否则 χ² 求和的 R 区间不同，F 不成立。
    # 权重拟合与相分数不是配位数问题，不纳入壳层 F-test。
    shell_f_test = None
    if len(plan["shells"]) >= 2 and not fit_weights:
        try:
            simple_fit = artemis.fit_r_space(
                res.chi_k, res.ft, [plan["files"][0]],
                rmin=plan["rmin"], rmax=plan["rmax"],
                kmin=res.ft.kmin, kmax=res.ft.kmax,
                kweight=(tuple(kweights) if kweights else res.ft.kweight),
                s02=absorber.s02, vary_delta_e0=True, vary_n=True,
                path_degen=[plan["degens"][0]],
                path_labels=[plan["labels"][0]],
                fit_space=fit_space)
            shell_f_test = artemis.shell_f_test(simple_fit, fit)
        except XAFSError as exc:
            shell_f_test = {"F": None, "confidence": None,
                            "significant": False,
                            "error": getattr(exc, "message", str(exc))}

    out = {"mode": ("motif_weights" if fit_weights else "motif"),
           "status": sid.get("status"), "motif": motif.id,
           "coverage_ok": sid.get("coverage_ok"), "cn_reliable": cn_reliable,
           "weights_reliable": weights_reliable,
           "mixture_weights": (fit.mixture_weights if fit_weights else None),
           "flags": flags,
           "structure": sid, "rmin": plan["rmin"], "rmax": plan["rmax"],
           "s02": absorber.s02, "s02_source": s02_label,
           "shells": plan["shells"], "fit": fit,
           "fit_space": fit.fit_space,
           "aic": fit.aic, "bic": fit.bic,
           "bvs": fit.bvs, "bvs_target": fit.bvs_target,
           "bvs_chi2": fit.bvs_chi2,
           "truncated_shells": plan["truncated"],
           "ms_links": plan["ms_links"],
           "path_files": [str(f) for f in plan["files"]],
           "path_links": plan.get("path_links") or [],
           "hetero_cocupancy": bool(hetero_cocupancy),
           "hetero_shells": plan.get("hetero_shells") or [],
           "nanoparticle_size": nanoparticle_size,
           "alloy_sro": alloy_sro,
           "resolution_check": resolution_check,
           "mm_bond_ftest": mm_bond_ftest,
           "shell_f_test": shell_f_test,
           "c3_dropped_by_nidp": bool(c3_dropped_by_nidp),
           "c3_at_bound": bool(c3_bound_bad),
           "c4_at_bound": bool(c4_bound_bad),
           "self_absorption": self_abs,
           "cn_systematics": systematic_budget(
               fit, flags=flags, self_abs=self_abs,
               truncated_shells=len(plan["truncated"]),
               fit_weights=fit_weights)}
    if conditional_cn is not None:
        out["conditional_cn"] = conditional_cn
    if fit_weights and fit.mixture_weights:
        w_max_err = max(float(w["weight_err"])
                        for w in fit.mixture_weights)
        out["weights_systematics"] = {
            "relative_1sigma": round(math.hypot(0.20, w_max_err), 4),
            "components": {"component_model_binding": 0.20,
                           "fitted_weight_error": round(w_max_err, 4)},
            "note": ("相分数误差 = 拟合形式误差与结构模型绑定误差"
                     "（0.20）平方和；组分数由模型选择决定，"
                     "不是独立物相定量。"),
        }
    if sigma2_groups:
        out["sigma2_groups"] = sigma2_groups
        out["sigma2_group_names"] = sigma2_group_names
        out["sigma2_shared"] = True
    else:
        out["sigma2_shared"] = False
    if fit_weights:
        out["warning"] = ("相分数是各组分 χ 的线性组合权重（有效简并度 "
                          "N_crystal·w），N 列不是配位数；组分数与结构模型"
                          "绑定，不是独立物相定量")
    else:
        out["warning"] = ("同壳层路径以代表路径 + 总简并度近似；"
                          "代表路径换一条会改变 F(k)/φ(k)，N 有系统误差")
    if plan["truncated"]:
        out["warning"] += ("；受 N_idp 限制，自动拟合截断了远壳层 "
                           + ", ".join(t["label"] for t in plan["truncated"]))
    if (flags["sub_rayleigh_z_degeneracy"]
            and shared_sigma2 is None and constraints is None):
        out["warning"] += ("；壳层间距低于 Rayleigh 极限且散射元素近简并，"
                           "建议绑定 shared_sigma2/constraints 或用 BVS 消歧")
    if c3_dropped_by_nidp:
        out["warning"] += ("；开启 C3 后自由变量数 ≥ N_idp，"
                           "已按 c3_dropped_by_nidp 自动降级/关闭 C3")
    if c4_bound_bad:
        out["warning"] += ("；C4 达到 ±80% 边界（±1.6e-4 Å⁴），"
                           "非高斯展宽项在吞噬失配，CN 不可靠")
    if isotropic_deltar and fit.r_factor > 0.08:
        out["warning"] += (f"；各向同性晶格约束 Rf={fit.r_factor:.4f} > 0.08："
                           "各壳层位移不服从同一 α·Reff 缩放，"
                           "不要用该约束下的 N 作配位数结论")
    if not sid.get("coverage_ok", True):
        out["warning"] += ("；最优候选未解释部分实验峰，"
                           "结构模型可能不完整，CN 解释需谨慎")
    if self_abs["risk"] in ("high", "unknown", "corrected_unreliable"):
        prefix = ("荧光自吸收校正不可靠" if self_abs["corrected"]
                  else "荧光自吸收未校正")
        out["warning"] += "；" + prefix + "（" + self_abs["note"] + "）"
    elif self_abs["risk"] == "possible":
        out["warning"] += "；荧光自吸收风险（" + self_abs["note"] + "）"
    reasons = []
    if fit.overfit:
        reasons.append("变量数 ≥ 独立点数")
    if fit.r_factor > 0.08:
        reasons.append(f"Rf={fit.r_factor:.3f}>0.08")
    if abs(fit.delta_e0) >= 15.0:
        reasons.append(f"|ΔE0|={abs(fit.delta_e0):.1f} eV≥15")
    if n_err_bad:
        reasons.append("N 形式误差过大")
    if sigma2_bound_bad:
        reasons.append("σ² 撞下界（形状失配证据"
                       + (f"：{', '.join(flags['sigma2_bound_shells'])}"
                          if flags["sigma2_bound_shells"] else "")
                       + "）")
    if c3_bound_bad:
        reasons.append("C3 撞 ±0.01 Å³ 边界（非谐项在吞噬相位失配"
                       + (f"：{', '.join(flags['c3_bound_shells'])}"
                          if flags["c3_bound_shells"] else "") + "）")
    if not coverage_ok:
        reasons.append("覆盖清单未通过/模型不完整")
    if self_abs["cn_gate"]:
        if self_abs["corrected"]:
            reasons.append("Booth 自吸收校正不可靠（截断/非有限值）")
        else:
            reasons.append(f"荧光自吸收未校正（Δμ="
                           f"{self_abs['edge_step']:.3g}，"
                           f"risk={self_abs['risk']}）")
    if fit_weights:
        if weight_at_bound:
            reasons.append("组分权重撞 softmax 边界")
        if weight_err_bad:
            reasons.append("组分权重误差过大")
        if weight_degenerate:
            pairs = "; ".join(" vs ".join(p)
                              for p in plan["weight_degenerate_pairs"])
            reasons.append(f"不同组分在同一壳层不可分辨（{pairs}）")
        if not weights_reliable:
            out["warning"] += ("；相分数不可靠（" + "、".join(reasons)
                               + "），不可作为物相定量结论")
    elif not cn_reliable:
        out["warning"] += ("；CN 不可靠（" + "、".join(reasons)
                           + "），不可直接用于配位数结论")
    if verbose:
        print(f"       ΔE0={fit.delta_e0:+.2f}±{fit.delta_e0_err:.2f} eV  "
              f"Rf={fit.r_factor:.4f}  vars={fit.n_variables}/"
              f"N_idp={fit.n_independent:.1f}"
              + ("" if fit.lattice_strain is None else
                 f"  α={fit.lattice_strain:+.5f}"
                 f"±{fit.lattice_strain_err:.5f}"))
        if fit_weights and fit.mixture_weights:
            print("       相分数（χ 线性组合权重，和=1）：")
            for w in fit.mixture_weights:
                print(f"         {w['group']:<18} w={w['weight']:.4f}"
                      f"±{w['weight_err']:.4f}")
        ntag = "N_eff" if fit_weights else "N"
        for p in fit.paths:
            print(f"       {p.label:<10} {ntag}={p.degeneracy:.2f}"
                  f"±{p.n_err:.2f}" if p.n_err is not None else
                  f"       {p.label:<10} {ntag}={p.degeneracy:.2f}", end="")
            print(f"  R={p.r:.3f} Å  σ²={p.sigma2:.4f}"
                  + (f"  C3={p.c3:+.2e} Å³" if p.c3 is not None else "")
                  + (f"  C4={p.c4:+.2e} Å⁴" if p.c4 is not None else "")
                  + ("  ⚠σ²撞下界" if p.sigma2_at_bound else "")
                  + ("  ⚠C3撞边界" if p.c3_at_bound else "")
                  + ("  ⚠C4撞边界" if p.c4_at_bound else ""))
        if fit_weights and not weights_reliable:
            print(f"       ⚠ 相分数不可靠：{'、'.join(reasons)}")
        elif not cn_reliable:
            print(f"       ⚠ CN 不可靠：{'、'.join(reasons)}")
    out.update(_formal_cn_meta(formal_cn, structure_model))
    if formal_cn and structure_model is not None:
        formal_block = out.get("formal_cn") or {}
        if formal_block.get("status") == "pending_consistency":
            if conditional_cn is not None and n0_mapping is not None:
                block = build_formal_cn_block(
                    conditional_cn=conditional_cn, plan=plan,
                    n0_mapping=n0_mapping,
                    structure_model=structure_model,
                    absorber=absorber,
                    cn_systematics=out.get("cn_systematics"),
                    s02_source=s02_label,
                    data_sha256=(getattr(res, "data_sha256", "") or ""))
            else:
                block = {
                    "status": "consistency_failed",
                    "rows": [],
                    "reject_reason": (conditional_cn_error
                                      or "条件 CN 未运行"),
                    "model": structure_model.to_json(),
                }
            out["formal_cn"] = block
            out["cn_status"] = block["status"]
    if formal_cn and catalytic_model is not None:
        out["formal_cn"] = {
            "status": "hypothesis_model_only",
            "rows": [],
            "reject_reason": (
                "catalytic_model 属于无实验 CIF 的催化位点假设模型，"
                "仅输出诊断 CN，不产出正式 CN"
            ),
            "model": None,
        }
        out["cn_status"] = "hypothesis_model_only"
    return out


def run(res, identify: bool = True, fit: bool = False,
        absorber_symbol: str | None = None, fit_motif: str | None = None,
        fit_rmin: float | None = None, fit_rmax: float | None = None,
        allow_ambiguous: bool = False, l1_ligand: str = "O",
        l1_geometry: str = "oct", fit_weights: bool = False,
        shared_sigma2: str | list | None = None,
        kweights: tuple | list | None = None,
        include_ms: bool = False,
        vary_c3: bool | list | str = False,
        vary_c4: bool | list | str = False,
        isotropic_deltar: bool = False,
        fit_space: str = "r",
        bvs_restraint: dict | None = None,
        catalytic_model=None,
        hetero_cocupancy: bool = False,
        constraints=None,
        s02: float | None = None,
        s02_source: str | None = None,
        verbose: bool = True,
        formal_cn: bool = False,
        structure_model: StructureModel | None = None,
        gate_result=None,
        xanes_cross_check: bool = False) -> dict:
    """S6.5（可选 S7）的唯一入口；结果回填 res.structure / res.fit。

    ``s02``（B18-2）为运行时 S₀² 覆盖：仅构造局部 AbsorberConfig 副本，
    不修改 ``absorber.ABSORBERS``；未标定元素需显式给出该值。
    ``xanes_cross_check``（T8）只追加 XANES–EXAFS 软校验警告，
    不改变拟合结果与 ``cn_reliable``。
    """
    out: dict = {}
    if identify or fit:
        # 显式催化模型不需要先跑库搜索：fit_identified 会自建状态。
        if catalytic_model is None:
            out["structure"] = structure_id.identify_result(
                res, absorber_symbol=absorber_symbol,
                l1_ligand=l1_ligand, l1_geometry=l1_geometry,
                vary_c3=vary_c3, fit_space=fit_space, s02=s02,
                verbose=verbose)
            res.structure = out["structure"]
    if fit:
        fit_out = fit_identified(
            res, motif_id=fit_motif, absorber_symbol=absorber_symbol,
            rmin=fit_rmin, rmax=fit_rmax,
            allow_ambiguous=allow_ambiguous,
            l1_ligand=l1_ligand, l1_geometry=l1_geometry,
            fit_weights=fit_weights,
            shared_sigma2=shared_sigma2, kweights=kweights,
            include_ms=include_ms, vary_c3=vary_c3, vary_c4=vary_c4,
            isotropic_deltar=isotropic_deltar,
            fit_space=fit_space, bvs_restraint=bvs_restraint,
            catalytic_model=catalytic_model,
            hetero_cocupancy=hetero_cocupancy,
            constraints=constraints,
            s02=s02,
            s02_source=s02_source,
            verbose=verbose,
            formal_cn=formal_cn,
            structure_model=structure_model,
            gate_result=gate_result)
        out["fit"] = fit_out
        res.fit = fit_out
        out["diagnostic_only"] = True
        out["cn_status"] = fit_out.get("cn_status", "diagnostic_only")
        if formal_cn:
            out["formal_cn"] = fit_out.get("formal_cn")
    elif formal_cn:
        out["formal_cn"] = {
            "status": "no_cif",
            "rows": [],
            "reject_reason": "未运行 S7 拟合（fit=False）",
            "model": (structure_model.to_json()
                      if structure_model is not None else None),
        }
        out["diagnostic_only"] = True
        out["cn_status"] = "no_cif"
    if xanes_cross_check:
        cc = xanes_exafs_cross_check(res, fit_out=out.get("fit"))
        out["xanes_cross_check"] = cc
        res.xanes_cross_check = cc
        if isinstance(out.get("fit"), dict):
            out["fit"]["xanes_cross_check"] = cc
        if verbose and cc["warnings"]:
            for w in cc["warnings"]:
                print(f"  [T8] ⚠ {w['code']}: {w['message']}")
    return out


def run_operando_exafs_series(
    results: list,
    absorber_symbol: str,
    motif: str | structure_lib.Motif | None = None,
    catalytic_model=None,
    rmin: float | None = None,
    rmax: float | None = None,
    kmin: float = 3.0,
    kmax: float = 12.0,
    kweight: int = 2,
    s02: float | None = None,
    share_sigma2_across_series: bool = False,
    names: list[str] | None = None,
) -> dict:
    """B19-1：原位序列 EXAFS 联合精修与 N(t)/R(t)/σ²(t) 轨迹。

    全序列共享**单一全局 ΔE₀**（``share_e0=True``），消除逐点 ΔE₀
    漂移经 ΔE₀–ΔR 相关转化成的虚假键长抖动；每个时间点独立放开
    各壳层 N/ΔR。``share_sigma2_across_series=True`` 时同一壳层在
    全序列共享 σ²，显著压低单点 N(t) 的形式误差。
    """
    if not results or len(results) < 3:
        raise XAFSConfigError(
            "OPERANDO", f"run_operando_exafs_series 至少需要 3 个时间点，"
                        f"收到 {0 if not results else len(results)}")
    if names is not None and len(names) != len(results):
        raise XAFSConfigError(
            "OPERANDO", f"names 长度 {len(names)} != 时间点数 {len(results)}")

    absorber = (get_absorber(absorber_symbol) if s02 is None
                else structure_id.absorber_with_s02(absorber_symbol, s02))
    if absorber.s02 is None:
        raise XAFSConfigError(
            "OPERANDO",
            f"{absorber.symbol} 未标定 S₀²；原位 EXAFS 序列联合精修需要 "
            f"显式 s02（0.4–1.3）")

    if catalytic_model is not None:
        target = structure_lib.parse_catalytic_cluster_spec(
            absorber.symbol, catalytic_model)
    elif motif is not None:
        if isinstance(motif, structure_lib.Motif):
            target = motif
            if str(target.absorber).capitalize() != absorber.symbol.capitalize():
                raise XAFSConfigError(
                    "OPERANDO", f"motif 吸收元素 {target.absorber!r} 与 "
                                f"absorber_symbol {absorber.symbol!r} 不一致")
        else:
            target = _find_motif(absorber.symbol, str(motif))
    else:
        sid = structure_id.identify_result(
            results[0], absorber_symbol=absorber.symbol, s02=s02,
            verbose=False)
        target = _find_motif(absorber.symbol, sid["ranked"][0]["id"])

    plan = select_shell_paths(
        target, edge=absorber.edge, rmin=rmin, rmax=rmax,
        kmin=kmin, kmax=kmax)
    labels = list(plan["labels"])
    n_points = len(results)
    names_list = (list(names) if names is not None
                  else [f"point{i}" for i in range(n_points)])

    datasets = []
    for di, res in enumerate(results):
        datasets.append({
            "name": names_list[di],
            "chi_k": res.chi_k,
            "ft": res.ft,
            "datfiles": list(plan["files"]),
            "rmin": float(plan["rmin"]),
            "rmax": float(plan["rmax"]),
            "kmin": float(kmin),
            "kmax": float(kmax),
            "kweight": int(kweight),
            "s02": float(absorber.s02),
            "vary_n": True,
            "path_degen": list(plan["degens"]),
            "path_labels": labels,
        })

    shared_bonds = None
    if share_sigma2_across_series:
        shared_bonds = []
        for pi in range(len(plan["files"])):
            for di in range(1, n_points):
                shared_bonds.append({
                    "target": (di, pi),
                    "source": (0, pi),
                    "share_deltar": False,
                    "share_sigma2": True,
                    "degen_ratio": None,
                })

    multi_out = artemis.fit_multi_datasets(
        datasets, shared_bonds=shared_bonds, share_e0=True)

    trajectories = {}
    for pi, label in enumerate(labels):
        trajectories[label] = {
            "reff": float(plan["shells"][pi]["reff"]),
            "N": [float(f.paths[pi].degeneracy)
                  for f in multi_out["datasets"]],
            "N_err": [f.paths[pi].n_err for f in multi_out["datasets"]],
            "R": [float(f.paths[pi].r) for f in multi_out["datasets"]],
            "R_err": [f.paths[pi].delta_r_err
                      for f in multi_out["datasets"]],
            "sigma2": [float(f.paths[pi].sigma2)
                       for f in multi_out["datasets"]],
            "sigma2_err": [f.paths[pi].sigma2_err
                           for f in multi_out["datasets"]],
        }

    return {
        "names": names_list,
        "motif": target.id,
        "shells": labels,
        "global_r_factor": float(multi_out["global_r_factor"]),
        "shared_delta_e0": float(multi_out["datasets"][0].delta_e0),
        "shared_delta_e0_err": multi_out["datasets"][0].delta_e0_err,
        "r_factors": [float(f.r_factor) for f in multi_out["datasets"]],
        "trajectories": trajectories,
        "multi_fit": multi_out,
    }


# ================================================================
# P2b：固定 N 一致性 —— N₀ shell 到 FEFF 路径的映射
# ================================================================
def map_n0_to_plan(plan: dict, n0_shells: list[dict],
                   kmin: float = 3.0, kmax: float = 12.0) -> dict:
    """把 CIFGateResult 的 n0_shells 映射到 plan["files"] 的每条 FEFF 路径。

    匹配规则：用 (scatterer, reff) 匹配 n0_shells 中的壳层；容差 0.5 Å。
    一个 CN 组可对应多条路径（子壳），每条路径的 path_degen 用该路径
    对应壳层的 N₀（不取组总和）。

    返回 dict：
        mapping_ok: bool          是否成功映射全部路径
        mapping_reason: str       失败原因（ok 时为空）
        n0_per_path: list[float]  每条路径的 N₀（用于 path_degen）
        path_links: list[dict]    同组多条路径的 degen_ratio 绑定
        n_var_fixed_n: int        固定 N 变量数 = 2n + 1
        n_var_free_n_first: int   自由 N（首壳）变量数 = 2n + 2
        n_paths: int              路径总数
    """
    n_paths = len(plan["files"])
    plan_shells = plan.get("shells", [])
    if len(plan_shells) != n_paths:
        return {
            "mapping_ok": False,
            "mapping_reason": (f"plan shells({len(plan_shells)}) != "
                               f"files({n_paths}), 无法一一映射"),
        }

    n0_per_path = [0.0] * n_paths
    shell_idx_per_path = [-1] * n_paths
    used_shell_idx: set[int] = set()

    # ---- 为每条路径匹配 n0_shell ----------------------------------
    for pi, sh in enumerate(plan_shells):
        scatterer = str(sh.get("scatterer", "")).capitalize()
        reff = float(sh.get("reff", sh.get("reff_mean", 0.0)))
        best_si = None
        best_dist = float("inf")
        for si, n0_sh in enumerate(n0_shells):
            if si in used_shell_idx:
                continue
            if str(n0_sh.get("el", "")).capitalize() != scatterer:
                continue
            dist = abs(float(n0_sh.get("r", 0.0)) - reff)
            if dist < best_dist and dist < 0.5:
                best_dist = dist
                best_si = si
        if best_si is None:
            return {
                "mapping_ok": False,
                "mapping_reason": (f"路径 {pi} ({scatterer}@{reff:.3f} Å) "
                                   f"无匹配 N₀ 壳层"),
            }
        used_shell_idx.add(best_si)
        shell_idx_per_path[pi] = best_si
        n0_per_path[pi] = float(n0_shells[best_si]["n0"])

    # ---- path_links：只绑定同一个 N₀ 壳层内的子壳 ----
    # 不同 N₀ 壳层是不同配位壳，绝不能互相绑定；否则远壳 N 会跟着
    # 首壳 N 漂移，违反“首壳放开、远壳固定 N₀”的条件 CN 约定。
    groups_by_shell: dict[int, list[int]] = {}
    for pi, si in enumerate(shell_idx_per_path):
        groups_by_shell.setdefault(si, []).append(pi)

    path_links: list[dict] = []
    for si, idxs in groups_by_shell.items():
        if len(idxs) <= 1:
            continue
        idxs.sort(key=lambda pi: float(plan_shells[pi].get("reff", 0.0)))
        ref = idxs[0]
        ref_n0 = n0_per_path[ref]
        for tgt in idxs[1:]:
            tgt_n0 = n0_per_path[tgt]
            ratio = (tgt_n0 / ref_n0) if ref_n0 > 0 else 1.0
            path_links.append({
                "target": int(tgt),
                "source": int(ref),
                "degen_ratio": round(ratio, 6),
                "share_deltar": False,
                "share_sigma2": False,
            })

    n_var_fixed_n = 2 * n_paths + 1      # 每壳 ΔR + σ²，加全局 ΔE₀
    n_var_free_n_first = 2 * n_paths + 2 # 首壳多放开一个 N

    return {
        "mapping_ok": True,
        "mapping_reason": "",
        "n0_per_path": [round(v, 6) for v in n0_per_path],
        "path_links": path_links,
        "n_var_fixed_n": n_var_fixed_n,
        "n_var_free_n_first": n_var_free_n_first,
        "n_paths": n_paths,
    }


# ================================================================
# P2b：select_shell_paths 信息量模式扩展
# ================================================================
_orig_select_shell_paths = select_shell_paths


def select_shell_paths(motif, edge="K", rmin=None, rmax=None,
                       merge_tol=0.20, max_paths=12, kmin=3.0, kmax=12.0,
                       group_by_component=False, include_ms=False,
                       max_ms=3, isotropic_deltar=False,
                       hetero_cocupancy=False,
                       n_mode: str = "free"):
    """select_shell_paths 的包装，支持 P2b 固定 N 信息量模式。

    n_mode:
        "free"      诊断/自由 N（默认）：3n+1 < N_idp
        "fixed_n"   固定 N：2n+1 < N_idp（每壳 ΔR+σ²，加 ΔE₀）
        "conditional" 条件 CN：2n+2 < N_idp（首壳多放开一个 N）
    """
    if n_mode not in ("free", "fixed_n", "conditional"):
        raise XAFSConfigError(
            "FIT", f"n_mode 必须是 free/fixed_n/conditional，收到 {n_mode!r}")
    plan = _orig_select_shell_paths(
        motif, edge=edge, rmin=rmin, rmax=rmax, merge_tol=merge_tol,
        max_paths=max_paths, kmin=kmin, kmax=kmax,
        group_by_component=group_by_component, include_ms=include_ms,
        max_ms=max_ms, isotropic_deltar=isotropic_deltar,
        hetero_cocupancy=hetero_cocupancy)
    if n_mode == "free":
        return plan
    # 重算信息量上限
    n_idp = plan["n_independent"]
    n_paths = len(plan["files"])
    if n_mode == "fixed_n":
        max_by_idp = max(1, int((n_idp - 1.0) // 2.0))
    else:  # conditional
        max_by_idp = max(1, int((n_idp - 2.0) // 2.0))
    limit = min(max_paths, max_by_idp)
    plan["n_paths_fixed_n"] = n_paths
    plan["n_var_fixed_n"] = 2 * n_paths + 1
    plan["n_var_conditional"] = 2 * n_paths + 2
    # 截截断：按新上限重新截断
    groups_count = len(plan["files"])
    if groups_count > limit:
        plan["truncated"] = plan.get("truncated", []) + [
            {"label": f"truncated_by_{n_mode}", "note":
             f"信息量上限 {n_mode}: {groups_count} > {limit}"}]
        plan["files"] = plan["files"][:limit]
        plan["degens"] = plan["degens"][:limit]
        plan["labels"] = plan["labels"][:limit]
        plan["shells"] = plan["shells"][:limit]
    return plan


# ================================================================
# P3：条件 CN —— 固定 N → 自由 N（首壳单一比例因子）
# ================================================================
CN_RATIO_OK_LO, CN_RATIO_OK_HI = 0.7, 1.3
CN_RATIO_WARN_LO, CN_RATIO_WARN_HI = 0.6, 1.4


def build_conditional_fit_args(
        plan: dict, n0_mapping: dict, ref_idx: int | None = None) -> dict:
    """构建条件 CN 拟合的 artemis.fit_r_space 参数。

    第一步（固定 N）：vary_n=[], path_degen=n0_per_path，所有壳层 N 固定为 N₀。
    第二步（条件 CN）：vary_n=[ref_idx]，只放开首壳参考路径 N；
    其他首壳路径用 path_links.degen_ratio = N0_j/N0_ref 绑定到参考路径；
    远壳 N 保持固定。

    返回 dict：
        fit_args: dict  透传给 artemis.fit_r_space 的参数
        ref_idx: int    首壳参考路径下标（N₀ 最大的首壳路径）
        n0_total: float 首壳 N₀ 总和
    """
    n0_per_path = n0_mapping["n0_per_path"]
    path_links = list(n0_mapping.get("path_links") or [])
    n_paths = len(n0_per_path)
    plan_shells = plan.get("shells", [])

    # 找出首壳路径（reff 最小的若干条，scatterer 可能相同或不同）
    first_shell_idxs = []
    if plan_shells:
        min_r = min(float(s.get("reff", s.get("reff_mean", 0.0)))
                    for s in plan_shells)
        for i, s in enumerate(plan_shells):
            r = float(s.get("reff", s.get("reff_mean", 0.0)))
            if abs(r - min_r) < 1e-6:
                first_shell_idxs.append(i)

    if not first_shell_idxs:
        first_shell_idxs = [0]

    # 参考路径：首壳中 N₀ 最大的
    ref_idx = ref_idx if ref_idx is not None else \
        max(first_shell_idxs, key=lambda i: n0_per_path[i])
    ref_idx = int(ref_idx)
    if ref_idx < 0 or ref_idx >= n_paths:
        raise XAFSConfigError(
            "FIT", f"ref_idx {ref_idx} 超出路径范围 0–{n_paths-1}")

    n0_total = sum(n0_per_path[i] for i in first_shell_idxs)

    # path_links：非参考首壳路径绑定到参考路径
    # 先清除已有的首壳间链接，重新按 N₀ 比例建立
    links = {}
    for link in path_links:
        links[link["target"]] = link

    for i in first_shell_idxs:
        if i == ref_idx:
            continue
        ratio = n0_per_path[i] / n0_per_path[ref_idx] if n0_per_path[ref_idx] > 0 else 1.0
        links[i] = {
            "target": i, "source": ref_idx,
            "degen_ratio": round(ratio, 6),
            "share_deltar": False,
            "share_sigma2": False,
        }

    fit_args = {
        "vary_n": [ref_idx],           # 只放开首壳参考路径 N
        "path_degen": [round(v, 6) for v in n0_per_path],  # 所有路径 N₀ 初值
        "path_links": list(links.values()),
    }
    return {
        "fit_args": fit_args,
        "ref_idx": ref_idx,
        "n0_total": round(n0_total, 6),
        "first_shell_idxs": first_shell_idxs,
    }


def compute_n_ratio(fit_result: dict, n0_mapping: dict,
                    ref_idx: int | None = None) -> dict:
    """从拟合结果计算 N/N₀。

    N_fit_total = Σ 首壳拟合路径的 N
    N0_total   = Σ 首壳 N₀
    N/N0       = N_fit_total / N0_total

    返回 dict：
        n_ratio: float | None
        n_fit_total: float
        n0_total: float
        n_status: str  "ok" / "warning" / "out_of_range" / "undef"
    """
    paths = fit_result.get("paths", [])
    n0_per_path = n0_mapping.get("n0_per_path", [])
    plan_shells = fit_result.get("shells", []) or \
        n0_mapping.get("plan_shells", [])

    if not paths or not n0_per_path:
        return {"n_ratio": None, "n_fit_total": None, "n0_total": None,
                "n_status": "undef"}

    # 首壳路径：reff 最小的
    first_idxs = []
    if plan_shells:
        min_r = min(float(s.get("reff", s.get("reff_mean", 0.0)))
                    for s in plan_shells)
        for i, s in enumerate(plan_shells):
            r = float(s.get("reff", s.get("reff_mean", 0.0)))
            if abs(r - min_r) < 1e-6:
                first_idxs.append(i)
    if not first_idxs:
        first_idxs = [0] if len(paths) > 0 else []

    n_fit_total = 0.0
    n0_total = 0.0
    for i in first_idxs:
        if i < len(paths):
            n_fit_total += float(paths[i].degen or 0.0)
        if i < len(n0_per_path):
            n0_total += float(n0_per_path[i])

    if n0_total <= 0:
        return {"n_ratio": None, "n_fit_total": round(n_fit_total, 6),
                "n0_total": round(n0_total, 6), "n_status": "undef"}

    n_ratio = n_fit_total / n0_total
    status = check_n_ratio(n_ratio)

    return {
        "n_ratio": round(n_ratio, 6),
        "n_fit_total": round(n_fit_total, 6),
        "n0_total": round(n0_total, 6),
        "n_status": status,
    }


def check_n_ratio(n_ratio: float) -> str:
    """N/N₀ 范围检查（P3 门禁）。

    0.7–1.3: ok（通过）
    0.6–0.7 / 1.3–1.4: warning（警告）
    <0.6 或 >1.4: out_of_range（需要定量缺陷模型，否则拒绝）
    """
    if n_ratio is None:
        return "undef"
    if n_ratio < CN_RATIO_WARN_LO or n_ratio > CN_RATIO_WARN_HI:
        return "out_of_range"
    if n_ratio < CN_RATIO_OK_LO or n_ratio > CN_RATIO_OK_HI:
        return "warning"
    return "ok"


# ================================================================
# P3：缺陷模型解析
# ================================================================
def parse_defect_model(spec: str | dict) -> dict:
    """解析缺陷模型规范，返回标准化 dict。

    支持 JSON 字符串或 dict：
        {"type": "nanoparticle", "diameter_nm": 2.5, "shape": "spherical", "lattice": "fcc"}
        {"type": "vacancy", "site": "O1", "fraction": 0.1}
        {"type": "custom", "expected_n_ratio": 0.75, "uncertainty": 0.10, "justification": "..."}

    返回 dict 包含：type, expected_n_ratio, uncertainty, model_dependent (bool)
    """
    if isinstance(spec, str):
        import json as _json
        try:
            spec = _json.loads(spec)
        except _json.JSONDecodeError as exc:
            raise XAFSConfigError(
                "FIT", f"缺陷模型 JSON 解析失败: {exc}") from exc
    if not isinstance(spec, dict):
        raise XAFSConfigError(
            "FIT", "缺陷模型必须是 JSON 对象或 dict")

    mtype = str(spec.get("type", "")).lower()
    if mtype not in ("nanoparticle", "vacancy", "custom"):
        raise XAFSConfigError(
            "FIT", f"未知缺陷模型类型: {mtype!r}（支持 nanoparticle/vacancy/custom）")

    result = {"type": mtype, "model_dependent": True}

    if mtype == "nanoparticle":
        dia = spec.get("diameter_nm")
        lattice = str(spec.get("lattice", "")).lower()
        if dia is None:
            raise XAFSConfigError(
                "FIT", "nanoparticle 缺陷模型需要 diameter_nm")
        if lattice not in ("fcc", "hcp", "bcc"):
            raise XAFSConfigError(
                "FIT", f"nanoparticle 仅对 fcc/hcp/bcc 开放，收到 {lattice!r}")
        u = 1.0 / float(dia)   # u = r1/d，简化：取 r1=1Å
        expected = 1.0 - 1.5 * u + 0.5 * u ** 3
        result["expected_n_ratio"] = round(max(0.0, expected), 6)
        result["uncertainty"] = float(spec.get("uncertainty", 0.10))
        result["diameter_nm"] = float(dia)
        result["lattice"] = lattice

    elif mtype == "vacancy":
        site = spec.get("site")
        fraction = spec.get("fraction")
        if site is None or fraction is None:
            raise XAFSConfigError(
                "FIT", "vacancy 缺陷模型需要 site 与 fraction")
        result["expected_n_ratio"] = round(1.0 - float(fraction), 6)
        result["uncertainty"] = float(spec.get("uncertainty", 0.10))
        result["site"] = str(site)
        result["fraction"] = float(fraction)

    else:  # custom
        ratio = spec.get("expected_n_ratio")
        uncertainty = spec.get("uncertainty")
        if ratio is None or uncertainty is None:
            raise XAFSConfigError(
                "FIT", "custom 缺陷模型需要 expected_n_ratio 与 uncertainty")
        result["expected_n_ratio"] = float(ratio)
        result["uncertainty"] = float(uncertainty)
        result["justification"] = str(spec.get("justification", ""))

    return result


def check_defect_model(n_ratio: float, defect_model: dict) -> dict:
    """检查 N/N₀ 是否在缺陷模型预期的不确定度内。

    返回 dict：
        ok: bool
        status: str  "defect_ok" / "defect_model_mismatch" / "defect_model_missing"
        detail: str
    """
    if defect_model is None:
        return {"ok": False, "status": "defect_model_missing",
                "detail": "N/N₀ 超出范围但未提供定量缺陷模型"}
    expected = defect_model.get("expected_n_ratio")
    uncertainty = defect_model.get("uncertainty", 0.10)
    if expected is None:
        return {"ok": False, "status": "defect_model_missing",
                "detail": "缺陷模型缺少 expected_n_ratio"}
    if abs(n_ratio - expected) <= uncertainty:
        return {"ok": True, "status": "defect_ok",
                "detail": f"N/N₀={n_ratio:.3f} 在模型预期 {expected:.3f}±{uncertainty:.2f} 范围内"}
    return {"ok": False, "status": "defect_model_mismatch",
            "detail": f"N/N₀={n_ratio:.3f} 偏离模型预期 {expected:.3f}±{uncertainty:.2f}"}




# ================================================================
# P3：条件 CN 端到端 —— 固定 N → 自由 N（首壳单一比例因子）→ Hamilton F-test
# ================================================================
def run_conditional_cn(
    res, plan: dict, n0_mapping: dict,
    absorber_symbol: str, s02: float,
    ref_idx: int | None = None,
    kmin: float = 3.0, kmax: float = 12.0,
    rmin: float | None = None, rmax: float | None = None,
    verbose: bool = True
) -> dict:
    """P3 条件 CN 端到端：固定 N 拟合 → 条件 CN（首壳单一比例因子）
    → 自由 N 对照 → Hamilton F-test 一致性判据。

    参数:
        res: PipelineResult（S1–S6 产物，含 chi_k, ft）
        plan: select_shell_paths 的输出（files, degens, labels, shells, rmin, rmax）
        n0_mapping: map_n0_to_plan 的输出（n0_per_path, path_links, ...）
        absorber_symbol: 吸收元素符号
        s02: 标定的 S₀² 值
        ref_idx: 首壳参考路径下标（None 时自动取 N₀ 最大的首壳路径）
        kmin/kmax: FT 窗口
        rmin/rmax: 拟合窗口（None 时用 plan 中的值）
        verbose: 是否打印进度

    返回 dict:
        conditional_cn: {
            n_ratio, n_fit_total, n0_total, n_status,   # N/N₀ 及状态
            ref_idx, first_shell_idxs,                  # 参考路径与首壳路径
            fixed_n_fit, conditional_n_fit, free_n_fit, # 三种拟合的摘要
            hamilton_f_test,                            # F-test 结果
            consistency_passed,                         # 一致性通过标志
            fixed_n_paths, conditional_n_paths          # 路径级 N 值
        }
    """
    n0_per_path = n0_mapping["n0_per_path"]
    path_links = list(n0_mapping.get("path_links") or [])
    n_paths = len(n0_per_path)
    plan_shells = plan.get("shells", [])

    if rmin is None:
        rmin = float(plan.get("rmin", 1.0))
    if rmax is None:
        rmax = float(plan.get("rmax", 4.0))

    # ---- 找出首壳路径（reff 最小）----
    first_shell_idxs = []
    if plan_shells:
        min_r = min(float(s.get("reff", s.get("reff_mean", 0.0)))
                    for s in plan_shells)
        for i, s in enumerate(plan_shells):
            r = float(s.get("reff", s.get("reff_mean", 0.0)))
            if abs(r - min_r) < 1e-6:
                first_shell_idxs.append(i)
    if not first_shell_idxs:
        first_shell_idxs = [0]

    # 参考路径：首壳中 N₀ 最大的
    if ref_idx is None:
        ref_idx = max(first_shell_idxs, key=lambda i: n0_per_path[i])
    ref_idx = int(ref_idx)
    if ref_idx < 0 or ref_idx >= n_paths:
        raise XAFSConfigError(
            "FIT", f"ref_idx {ref_idx} 超出路径范围 0–{n_paths-1}")

    n0_total = sum(float(n0_per_path[i]) for i in first_shell_idxs)
    if n0_total <= 0:
        raise XAFSConfigError(
            "FIT", f"首壳 N₀ 总和为 {n0_total}，无法计算 N/N₀")

    # ---- 第一步：固定 N 拟合（vary_n=[], path_degen=n0_per_path）----
    if verbose:
        print(f"  [P3-1] 固定 N 拟合：n_paths={n_paths}, "
              f"vary_n=[], N 固定为 N₀")
    fit_fixed = artemis.fit_r_space(
        res.chi_k, res.ft, plan["files"],
        rmin=rmin, rmax=rmax,
        kmin=res.ft.kmin, kmax=res.ft.kmax,
        kweight=res.ft.kweight,
        s02=s02, vary_delta_e0=True, vary_n=[],
        path_degen=[round(v, 6) for v in n0_per_path],
        path_labels=plan["labels"],
        fit_space="r")

    # ---- 第二步：条件 CN（vary_n=[ref_idx]，首壳非参考路径用 path_links 绑定）----
    links = {}
    for link in path_links:
        links[link["target"]] = link

    for i in first_shell_idxs:
        if i == ref_idx:
            continue
        ratio = n0_per_path[i] / n0_per_path[ref_idx] if n0_per_path[ref_idx] > 0 else 1.0
        links[i] = {
            "target": i, "source": ref_idx,
            "degen_ratio": round(ratio, 6),
            "share_deltar": False,
            "share_sigma2": False,
        }

    if verbose:
        print(f"  [P3-2] 条件 CN 拟合：vary_n=[{ref_idx}], "
              f"path_links={len(links)} 条（首壳非参考路径绑定到参考路径）")
    fit_conditional = artemis.fit_r_space(
        res.chi_k, res.ft, plan["files"],
        rmin=rmin, rmax=rmax,
        kmin=res.ft.kmin, kmax=res.ft.kmax,
        kweight=res.ft.kweight,
        s02=s02, vary_delta_e0=True, vary_n=[ref_idx],
        path_degen=[round(v, 6) for v in n0_per_path],
        path_labels=plan["labels"],
        path_links=list(links.values()),
        fit_space="r")

    # ---- 第三步：自由 N 对照（所有 N 可变，path_degen=n0_per_path 作为初值）----
    if verbose:
        print(f"  [P3-3] 自由 N 对照拟合：vary_n=True，所有 N 自由")
    fit_free = artemis.fit_r_space(
        res.chi_k, res.ft, plan["files"],
        rmin=rmin, rmax=rmax,
        kmin=res.ft.kmin, kmax=res.ft.kmax,
        kweight=res.ft.kweight,
        s02=s02, vary_delta_e0=True, vary_n=True,
        path_degen=[round(v, 6) for v in n0_per_path],
        path_labels=plan["labels"],
        fit_space="r")

    # ---- 计算 N/N₀（基于条件 CN 拟合结果）----
    n_fit_total = 0.0
    for i in first_shell_idxs:
        if i < len(fit_conditional.paths):
            n_fit_total += float(fit_conditional.paths[i].degeneracy or 0.0)

    n_ratio = n_fit_total / n0_total
    n_status = check_n_ratio(n_ratio)

    # ---- Hamilton F-test：固定 N vs 条件 CN（嵌套模型比较）----
    try:
        ftest = artemis.shell_f_test(fit_fixed, fit_conditional)
    except XAFSError as exc:
        ftest = {"F": None, "confidence": None, "significant": False,
                 "error": getattr(exc, "message", str(exc))}

    # 一致性判据：
    # - N/N₀ 在 ok/warning 范围
    # - 自由 N 不显著改善（Hamilton F-test 不显著）→ 固定 N 充分
    # - 不过拟合，Rf 合理
    ftest_significant = bool(ftest.get("significant", False))
    consistency_passed = (
        n_status in ("ok", "warning")
        and not ftest_significant
        and not fit_conditional.overfit
        and float(fit_conditional.r_factor) <= 0.08
    )

    if verbose:
        print(f"  [P3-4] N/N₀ = {n_ratio:.4f} ({n_status}), "
              f"F-test significant={ftest_significant}, "
              f"consistency_passed={consistency_passed}")

    return {
        "n_ratio": round(n_ratio, 6),
        "n_fit_total": round(n_fit_total, 6),
        "n0_total": round(n0_total, 6),
        "n_status": n_status,
        "ref_idx": ref_idx,
        "first_shell_idxs": first_shell_idxs,
        "fixed_n_fit": {
            "r_factor": float(fit_fixed.r_factor),
            "n_variables": int(fit_fixed.n_variables),
            "n_independent": float(fit_fixed.n_independent),
            "delta_e0": float(fit_fixed.delta_e0),
            "overfit": bool(fit_fixed.overfit),
        },
        "conditional_n_fit": {
            "r_factor": float(fit_conditional.r_factor),
            "n_variables": int(fit_conditional.n_variables),
            "n_independent": float(fit_conditional.n_independent),
            "delta_e0": float(fit_conditional.delta_e0),
            "overfit": bool(fit_conditional.overfit),
        },
        "free_n_fit": {
            "r_factor": float(fit_free.r_factor),
            "n_variables": int(fit_free.n_variables),
            "n_independent": float(fit_free.n_independent),
            "delta_e0": float(fit_free.delta_e0),
            "overfit": bool(fit_free.overfit),
        },
        "hamilton_f_test": {
            "F": ftest.get("F"),
            "confidence": ftest.get("confidence"),
            "significant": ftest_significant,
            "p_value": (None if ftest.get("confidence") is None
                        else float(max(0.0, 1.0 - float(ftest.get("confidence"))))),
        },
        "consistency_passed": bool(consistency_passed),
        "kmin": float(res.ft.kmin),
        "kmax": float(res.ft.kmax),
        "kweight": int(res.ft.kweight),
        "fixed_n_paths": [
            {"label": p.label, "reff": float(p.r),
             "N": float(p.degeneracy),
             "N_err": (float(p.n_err) if p.n_err is not None else None),
             "delta_r": float(p.delta_r),
             "sigma2": float(p.sigma2)}
            for p in fit_fixed.paths
        ],
        "conditional_n_paths": [
            {"label": p.label, "reff": float(p.r),
             "N": float(p.degeneracy),
             "N_err": (float(p.n_err) if p.n_err is not None else None),
             "delta_r": float(p.delta_r),
             "sigma2": float(p.sigma2)}
            for p in fit_conditional.paths
        ],
    }


def build_formal_cn_block(
    *,
    conditional_cn: dict,
    plan: dict,
    n0_mapping: dict,
    structure_model: StructureModel,
    absorber,
    cn_systematics: dict | None = None,
    s02_source: str = "absorber.py",
    data_sha256: str = "",
    software_version: str = cn_table.SOFTWARE_VERSION,
) -> dict:
    """把条件 CN 结果转成正式 CN 块。

    - ``consistency_passed=True``：``status="ok"``，每个拟合路径一行；
      第一壳 ``n_status="fitted"``，远壳 ``n_status="fixed"``。
    - 否则：``status="consistency_failed"``，``rows=[]``，拒绝原因可复核。
    - 每行包含 §4.4 要求的全部 provenance 字段：CIF 来源、sha256、
      吸收位、N₀、模型级别、S₀²、数据哈希与软件版本。
    """
    model = structure_model.to_json()
    cond_fit = conditional_cn.get("conditional_n_fit") or {}
    rf = float(cond_fit.get("r_factor") or 0.0)
    overall_n_status = str(conditional_cn.get("n_status") or "")
    ftest = conditional_cn.get("hamilton_f_test") or {}
    reasons = []
    if overall_n_status not in ("ok", "warning"):
        reasons.append(
            f"N/N₀={conditional_cn.get('n_ratio')} "
            f"({overall_n_status}) 超出范围，且未提供定量缺陷模型")
    if bool(ftest.get("significant")):
        p_value = ftest.get("p_value")
        reasons.append(
            "自由 N 显著改善（Hamilton F-test"
            + (f" p={float(p_value):.4g}" if p_value is not None else "")
            + "），固定 N 不一致")
    if bool(cond_fit.get("overfit")):
        reasons.append("条件 CN 拟合变量数 ≥ N_idp")
    if rf > 0.08:
        reasons.append(f"条件 CN Rf={rf:.4f} > 0.08")
    if not conditional_cn.get("consistency_passed"):
        return {
            "status": "consistency_failed",
            "rows": [],
            "reject_reason": "；".join(reasons) or "条件 CN 一致性未通过",
            "model": model,
        }

    paths = conditional_cn.get("conditional_n_paths") or []
    first_shell = {int(i) for i in
                   (conditional_cn.get("first_shell_idxs") or [])}
    n0_per_path = list(n0_mapping.get("n0_per_path") or [])
    shells = plan.get("shells") or []
    systematics = (cn_systematics or {}).get("relative_1sigma")
    row_gate = "passed" if overall_n_status == "ok" else "warning"
    row_flags = ("consistency_passed" if overall_n_status == "ok"
                 else "consistency_passed;n_ratio_warning")
    assumptions = (
        "first_shell_single_scale_factor;"
        "higher_shells_fixed_at_N0;"
        "hamilton_f_test_nonsignificant;"
        "data_sha256_over_input_files"
    )

    rows: list[dict] = []
    for i, shell in enumerate(shells):
        if i >= len(paths):
            break
        path = paths[i]
        fitted = i in first_shell
        n_fit = float(path.get("N") or 0.0)
        n0 = (float(n0_per_path[i]) if i < len(n0_per_path)
              else float(shell.get("degen") or 0.0))
        n_ratio = ((n_fit / n0) if fitted and n0 > 0
                   else (1.0 if not fitted else ""))
        row = {column: "" for column in cn_table.CN_TABLE_COLUMNS}
        row.update({
            "phase_id": structure_model.model_id,
            "model_level": structure_model.model_level,
            "source_type": structure_model.source_type,
            "source_ref": structure_model.source_ref,
            "source_cif_sha256": structure_model.source_cif_sha256,
            "absorber_site": structure_model.absorber_site,
            "shell_group": "first_shell" if fitted else "higher_shell",
            "shell_label": shell.get("label", ""),
            "scatterer": shell.get("scatterer", ""),
            "r_cif": float(shell.get("reff")
                           or shell.get("reff_mean") or 0.0),
            "n0": n0,
            "n0_members": shell.get("n_members", ""),
            "n_fit": n_fit,
            "n_status": "fitted" if fitted else "fixed",
            "n_ratio": n_ratio,
            "n_err_formal": (path.get("N_err") if fitted else ""),
            "n_err_systematic": (systematics if fitted else ""),
            "delta_r": path.get("delta_r", ""),
            "sigma2": path.get("sigma2", ""),
            "delta_e0": cond_fit.get("delta_e0", ""),
            "rf": rf,
            "s02": absorber.s02,
            "s02_source": s02_source,
            "s02_calibrated": _s02_calibrated_flag(absorber, s02_source),
            "kmin": conditional_cn.get("kmin", ""),
            "kmax": conditional_cn.get("kmax", ""),
            "kweight": conditional_cn.get("kweight", ""),
            "rmin": plan.get("rmin", ""),
            "rmax": plan.get("rmax", ""),
            "n_idp": cond_fit.get("n_independent", ""),
            "n_variables": cond_fit.get("n_variables", ""),
            "cn_label": "cn_conditional",
            "gate_status": row_gate,
            "reject_reason": "",
            "flags": row_flags,
            "assumptions": assumptions,
            "data_sha256": data_sha256,
            "software_version": software_version,
        })
        rows.append(row)
    return {"status": "ok", "rows": rows, "reject_reason": None,
            "model": model}


# ================================================================
# P5-D：多位点正式 CN（--absorber-site all）
# ================================================================
def _site_fractional(model: StructureModel) -> list[float]:
    """从 StructureModel 取吸收位的分数坐标。"""
    label = str(model.absorber_site or "")
    m = re.match(r"frac:\(([^,]+),([^,]+),([^)]+)\)", label)
    if m:
        return [float(m.group(1)), float(m.group(2)), float(m.group(3))]
    for site in (model.asym_sites or []):
        if str(site.get("label")) == label:
            return [float(site["x"]), float(site["y"]), float(site["z"])]
    raise CNError(
        "absorber_site_missing",
        f"{model.model_id}: 找不到吸收位 {label!r} 的分数坐标")


def _site_motif_from_model(model: StructureModel) -> structure_lib.Motif:
    """由 CIF + 吸收位坐标构建该位点的团簇 motif（不写生产结构库）。"""
    if not model.source_cif:
        raise CNError("no_cif", f"{model.model_id}: 缺少 source_cif")
    sites, occs = artemis.read_cif_sites(
        model.source_cif, absorber=model.absorber, rmax=6.0,
        absorber_fract=_site_fractional(model), with_occupancy=True)
    if not sites:
        raise CNError(
            "absorber_missing",
            f"{model.model_id}@{model.absorber_site}: rmax=6.0 Å 内无近邻")
    # motif id 会参与 FEFF 缓存目录名（Windows 不允许 ':'），必须做
    # 文件系统安全化；只影响诊断缓存，不改变 phase/site 溯源字段。
    safe_id = re.sub(r"[^A-Za-z0-9_.@-]+", "_",
                     f"{model.model_id}@{model.absorber_site}")
    return structure_lib.Motif(
        id=safe_id,
        absorber=model.absorber,
        kind="cluster",
        nuclearity="m1",
        source={"type": "cif_site", "ref": model.source_ref or ""},
        sites=sites,
        shell_expect=structure_lib.geometric_shells(sites, occs=occs),
        absorber_site=model.absorber_site,
        source_ref=model.source_ref,
        source_cif_sha256=model.source_cif_sha256,
        model_level=model.model_level,
        provenance_status=model.provenance_status,
    )


def _num(value) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _first_shell_stats(rows: list[dict]) -> dict:
    """把一个位点 formal 行里的第一壳行汇总成 N₀/N_fit/误差/键长。"""
    first = [r for r in rows if r.get("shell_group") == "first_shell"]
    n0 = sum(_num(r.get("n0")) or 0.0 for r in first)
    n_fit = sum(_num(r.get("n_fit")) or 0.0 for r in first)
    errors = [e for e in (_num(r.get("n_err_formal")) for r in first)
              if e is not None]
    r_terms = [(float(r["r_cif"]), _num(r.get("n0")) or 0.0)
               for r in first if _num(r.get("r_cif")) is not None]
    sig_terms = [(_num(r.get("sigma2")), _num(r.get("n0")) or 0.0)
                 for r in first if _num(r.get("sigma2")) is not None]
    dr_terms = [(_num(r.get("delta_r")), _num(r.get("n0")) or 0.0)
                for r in first if _num(r.get("delta_r")) is not None]
    e0_terms = [(_num(r.get("delta_e0")), _num(r.get("n0")) or 0.0)
                for r in first if _num(r.get("delta_e0")) is not None]
    rf = [_num(r.get("rf")) for r in first]
    rf = [v for v in rf if v is not None]
    scatterers = sorted({str(r.get("scatterer") or "") for r in first
                         if r.get("scatterer")})
    return {
        "n0": n0,
        "n_fit": n_fit,
        "n_err": (math.sqrt(sum(e * e for e in errors))
                  if errors else None),
        "r_cif": (sum(r * n for r, n in r_terms) / n0 if n0 > 0
                  and r_terms else None),
        "sigma2": (sum(s * n for s, n in sig_terms) / n0 if n0 > 0
                   and sig_terms else None),
        "delta_r": (sum(v * n for v, n in dr_terms) / n0 if n0 > 0
                    and dr_terms else None),
        "delta_e0": (sum(v * n for v, n in e0_terms) / n0 if n0 > 0
                     and e0_terms else None),
        "rf": (max(rf) if rf else None),
        "scatterer": (scatterers[0] if len(scatterers) == 1 else "mixed"),
        "n_rows": len(first),
    }


def run_multi_site_formal_cn(
    res,
    *,
    site_models: list[StructureModel],
    gate_results,
    s02: float | None = None,
    s02_source: str = "absorber.py",
    data_sha256: str = "",
    software_version: str = cn_table.SOFTWARE_VERSION,
    verbose: bool = False,
) -> dict:
    """P5-D：对 ``--absorber-site all`` 的每个不等价位点独立做条件 CN。

    输出契约：
    - 每个位点各自的门禁、条件 CN 行（``flags`` 带 ``multi_site_member``
      与 ``site_weight``）；
    - 仅当**全部**位点 ``build_formal_cn_block`` 返回 ``status="ok"`` 时，
      追加一行 ``absorber_site="ALL_WEIGHTED"`` 的化学计量加权汇总；
    - 任一位点失败：该位点写拒绝行，且**不输出**加权汇总行。
    """
    if not site_models:
        raise XAFSConfigError("FIT", "多位点正式 CN 缺少结构模型")
    if len(site_models) != len(gate_results):
        raise XAFSConfigError(
            "FIT", f"多位点模型数({len(site_models)}) != 门禁结果数"
                   f"({len(gate_results)})")
    absorber_symbol = site_models[0].absorber
    absorber = (get_absorber(absorber_symbol) if s02 is None
                else structure_id.absorber_with_s02(absorber_symbol, s02))
    if absorber.s02 is None:
        raise XAFSConfigError(
            "FIT", f"{absorber.symbol} 未标定 S₀²，拒绝多位点 CN 定量")

    rows: list[dict] = []
    per_site: list[dict] = []
    all_ok = True
    for model, gate in zip(site_models, gate_results):
        record = {
            "absorber_site": model.absorber_site,
            "site_weight": getattr(gate, "site_weight", None),
            "site_multiplicity": getattr(gate, "site_multiplicity", None),
            "site_occupancy": getattr(gate, "site_occupancy", None),
            "status": "gate_failed",
            "rows": [],
        }
        if not getattr(gate, "ok", False):
            record["gate_code"] = getattr(gate, "code", None)
            record["gate_reason"] = getattr(gate, "reason", "")
            all_ok = False
            rows.append(cn_table.refusal_row(
                "gate_failed", getattr(gate, "reason", "") or "CIF gate 拒绝",
                phase_id=model.model_id,
                model_level=model.model_level,
                source_type=model.source_type,
                source_ref=model.source_ref,
                source_cif_sha256=model.source_cif_sha256,
                absorber_site=model.absorber_site,
                s02=absorber.s02,
                s02_source=s02_source,
                s02_calibrated=_s02_calibrated_flag(absorber, s02_source),
            ))
            per_site.append(record)
            continue

        try:
            motif = _site_motif_from_model(model)
            plan = select_shell_paths(
                motif, edge=absorber.edge,
                kmin=res.ft.kmin, kmax=res.ft.kmax)
            n0_mapping = map_n0_to_plan(
                plan, gate.n0_shells,
                kmin=res.ft.kmin, kmax=res.ft.kmax)
            if not n0_mapping["mapping_ok"]:
                raise CNError("n0_mapping_failed", n0_mapping["mapping_reason"])
            conditional = run_conditional_cn(
                res, plan, n0_mapping,
                absorber_symbol=absorber.symbol, s02=absorber.s02,
                kmin=res.ft.kmin, kmax=res.ft.kmax,
                rmin=plan["rmin"], rmax=plan["rmax"],
                verbose=False)
            block = build_formal_cn_block(
                conditional_cn=conditional, plan=plan,
                n0_mapping=n0_mapping, structure_model=model,
                absorber=absorber, s02_source=s02_source,
                data_sha256=data_sha256,
                software_version=software_version)
            if block["status"] != "ok":
                all_ok = False
                record.update({
                    "status": block["status"],
                    "reject_reason": block.get("reject_reason"),
                })
                rows.append(cn_table.refusal_row(
                    block["status"], block.get("reject_reason") or "",
                    phase_id=model.model_id,
                    model_level=model.model_level,
                    source_type=model.source_type,
                    source_ref=model.source_ref,
                    source_cif_sha256=model.source_cif_sha256,
                    absorber_site=model.absorber_site,
                    s02=absorber.s02,
                    s02_source=s02_source,
                    s02_calibrated=_s02_calibrated_flag(
                        absorber, s02_source),
                ))
                per_site.append(record)
                continue

            weight = float(getattr(gate, "site_weight", 0.0) or 0.0)
            for row in block["rows"]:
                flags = row.get("flags") or ""
                row["flags"] = (
                    flags + ";" if flags else ""
                ) + f"multi_site_member;site_weight={weight:.6f}"
            rows.extend(block["rows"])
            stats = _first_shell_stats(block["rows"])
            record.update({
                "status": "ok",
                "site_weight": weight,
                "n0": stats["n0"],
                "n_fit": stats["n_fit"],
                "n_ratio": (stats["n_fit"] / stats["n0"]
                            if stats["n0"] > 0 else None),
                "rf": stats["rf"],
                "rows": len(block["rows"]),
            })
            per_site.append(record)
        except (XAFSConfigError, XAFSDataError, CNError) as exc:
            all_ok = False
            reason = getattr(exc, "message", str(exc))
            record.update({"status": "error", "reject_reason": reason})
            rows.append(cn_table.refusal_row(
                "multi_site_error", reason,
                phase_id=model.model_id,
                model_level=model.model_level,
                source_type=model.source_type,
                source_ref=model.source_ref,
                source_cif_sha256=model.source_cif_sha256,
                absorber_site=model.absorber_site,
                s02=absorber.s02,
                s02_source=s02_source,
                s02_calibrated=_s02_calibrated_flag(absorber, s02_source),
            ))
            per_site.append(record)

    result = {
        "status": "ok" if all_ok else "multi_site_failed",
        "rows": rows,
        "per_site": per_site,
        "multi_site": True,
        "site_weights": {
            str(r["absorber_site"]): r.get("site_weight") for r in per_site},
        "reject_reason": (None if all_ok else
                          "至少一个吸收位点未通过门禁/条件 CN；"
                          "不输出化学计量加权汇总"),
        "model": site_models[0].to_json(),
        "models": [m.to_json() for m in site_models],
    }
    if not all_ok:
        return result

    # ---- 全部位点通过：化学计量加权汇总行 ----
    weights = [float(getattr(g, "site_weight", 0.0) or 0.0)
               for g in gate_results]
    site_stats = [_first_shell_stats(block_rows)
                  for block_rows in ([r for r in rows
                                      if r.get("absorber_site") == m.absorber_site]
                                     for m in site_models)]
    n0_total = sum(w * s["n0"] for w, s in zip(weights, site_stats))
    n_fit_total = sum(w * s["n_fit"] for w, s in zip(weights, site_stats))
    r_num = sum(w * s["n0"] * (s["r_cif"] or 0.0)
                for w, s in zip(weights, site_stats))
    sig_num = sum(w * s["n0"] * (s["sigma2"] or 0.0)
                  for w, s in zip(weights, site_stats))
    err_terms = [w * s["n_err"] for w, s in zip(weights, site_stats)
                 if s["n_err"] is not None]
    rf_vals = [s["rf"] for s in site_stats if s["rf"] is not None]
    first_rows = [r for r in rows if r.get("shell_group") == "first_shell"]
    scatterers = sorted({str(r.get("scatterer"))
                         for r in first_rows if r.get("scatterer")})
    dr_num = sum(w * s["n0"] * (s["delta_r"] or 0.0)
                 for w, s in zip(weights, site_stats))
    e0_num = sum(w * s["n0"] * (s["delta_e0"] or 0.0)
                 for w, s in zip(weights, site_stats))

    summary = {column: "" for column in cn_table.CN_TABLE_COLUMNS}
    summary.update({
        "phase_id": site_models[0].model_id,
        "model_level": site_models[0].model_level,
        "source_type": site_models[0].source_type,
        "source_ref": site_models[0].source_ref,
        "source_cif_sha256": site_models[0].source_cif_sha256,
        "absorber_site": "ALL_WEIGHTED",
        "shell_group": "first_shell_weighted",
        "shell_label": "+".join(
            f"{m.absorber_site}({w:.3f})"
            for m, w in zip(site_models, weights)),
        "scatterer": (scatterers[0] if len(scatterers) == 1 else "mixed"),
        "r_cif": (r_num / n0_total if n0_total > 0 else ""),
        "n0": n0_total,
        "n0_members": sum(int(s["n_rows"]) for s in site_stats),
        "n_fit": n_fit_total,
        "n_status": "weighted_average",
        "n_ratio": (n_fit_total / n0_total if n0_total > 0 else ""),
        "n_err_formal": (math.sqrt(sum(e * e for e in err_terms))
                         if err_terms else ""),
        "n_err_systematic": CN_SYSTEMATIC_BUDGET["relative_1sigma"],
        "delta_r": (dr_num / n0_total if n0_total > 0 else ""),
        "sigma2": (sig_num / n0_total if n0_total > 0 else ""),
        "delta_e0": (e0_num / n0_total if n0_total > 0 else ""),
        "rf": (max(rf_vals) if rf_vals else ""),
        "s02": absorber.s02,
        "s02_source": s02_source,
        "s02_calibrated": _s02_calibrated_flag(absorber, s02_source),
        "kmin": min(_num(r.get("kmin")) or 0.0 for r in first_rows)
                if first_rows else "",
        "kmax": max(_num(r.get("kmax")) or 0.0 for r in first_rows)
                if first_rows else "",
        "kweight": first_rows[0].get("kweight", "") if first_rows else "",
        "rmin": min(_num(r.get("rmin")) or 0.0 for r in first_rows)
                if first_rows else "",
        "rmax": max(_num(r.get("rmax")) or 0.0 for r in first_rows)
                if first_rows else "",
        "n_idp": min(_num(r.get("n_idp")) or 0.0 for r in first_rows)
                 if first_rows else "",
        "n_variables": sum(int(_num(r.get("n_variables")) or 0)
                           for r in first_rows),
        "cn_label": "cn_conditional_multi_site",
        "gate_status": ("passed"
                        if all(r.get("gate_status") == "passed"
                               for r in first_rows) else "warning"),
        "reject_reason": "",
        "flags": "multi_site_weighted_summary",
        "assumptions": (
            "site_weight=m_i*occ_i_normalized;"
            "each_site_fit_independent;"
            "first_shell_weighted_by_stoichiometry;"
            "no_partial_site_release"),
        "data_sha256": data_sha256,
        "software_version": software_version,
    })
    rows.append(summary)
    if verbose:
        print(f"  [P5-D] {len(site_models)} 个吸收位点全部通过；"
              f"加权 N₀={n0_total:.4f}  N_fit={n_fit_total:.4f}")
    return result
