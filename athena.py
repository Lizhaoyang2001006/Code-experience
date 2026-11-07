# -*- coding: utf-8 -*-
"""
athena.py — 阶段 1~4：由 μ(E) 到 |χ(R)|
================================================================
对应 Athena 的功能范围。函数名与研究流程里的数学步骤一一对应：

    μ(E)          mu_from_*                     （见 xafs_io）
    能量校准       find_e0_derivative_max / find_e0_half_height
    多扫描合并     align_by_reference / merge_scans
    归一化         fit_linear_pre_edge → compute_edge_step → normalize_absorption
    背景扣除       wave_number → extract_chi_k（内部调 AUTOBK）
    傅里叶变换     ft_to_r_space（内部调 xftf，窗口由其统一处理）

依赖：xafs_core（常数与数据结构）、xafs_io（读写）、larch（AUTOBK 算法）。
**不 import artemis / larch_wt**，可单独运行 `python athena.py` 自测。
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from larch import Group
from larch.xafs import autobk
from xraydb import (material_mu, material_mu_components, mu_elam, xray_edge,
                    xray_line)

from xafs_core import (K2, BackFourierTransform, ChiK, FourierTransform,
                       NormalizedScan, RawScan, XAFSConfigError, XAFSDataError)


# ================================================================
# 1. 多扫描合并
# ================================================================

def align_by_reference(scans, edge_energy: float, tolerance: float = 0.2):
    """
    用参考通道的边位移把各扫描平移到同一能量刻度。

    数学:
        ΔE_i = E_ref_standard − E_ref_measured,i
        E_corrected = E + ΔE_i          （整体平移，不做插值）

    约束:
        样品与参考**共用同一个 ΔE**，不各自对齐。
        参考通道必须「独立」——即其边位与样品边位相差 >3 eV。
        若两者边位几乎相同，说明该通道测的是同一样品（多为下游电离室），
        此时拒绝对齐，而不是硬凑一个平移量。

    返回: (对齐后的 scans, info 字典)
    """
    e0_refs_collected = []
    for i, scan in enumerate(scans):
        # guess=edge_energy：文件头可能无 e0_nominal（SIMPLE 格式），
        # 边窗口必须用调用方已解析的标称边能，绝不回退到任何默认元素
        info = determine_e0(scan, guess=edge_energy)
        if info.get("e0_reference") is None:
            raise XAFSDataError(
                "S2", f"扫描 {i+1} 无参考通道（Ir），无法做能量对齐")
        shift = info["e0_sample"] - info["e0_reference"]
        if abs(shift) < 3.0:
            raise XAFSDataError(
                "S2", f"扫描 {i+1} 的参考通道边位与样品仅差 {shift:+.2f} eV，"
                      f"测的是同一吸收体（多为样品下游电离室），不可用于对齐")
        e0_refs_collected.append(info["e0_reference"])

    e0_ref = np.asarray(e0_refs_collected)
    shifts = edge_energy - e0_ref

    aligned = []
    for scan, dE in zip(scans, shifts):
        moved = RawScan(
            energy=np.asarray(scan.energy) + dE,
            i0=scan.i0, it=scan.it, ifluor=scan.ifluor, iref=scan.iref,
            dtime=scan.dtime, mu=scan.mu, sigma_mu=scan.sigma_mu,
            meta=dict(scan.meta))
        moved.meta["align_shift_eV"] = float(dE)
        aligned.append(moved)

    e0_after = np.array([determine_e0(s, guess=edge_energy)["e0_sample"]
                         for s in aligned])
    residual = float(np.std(e0_after))
    return aligned, {
        "e0_reference_per_scan": np.round(e0_ref, 3).tolist(),
        "shifts_eV": np.round(shifts, 3).tolist(),
        "e0_sample_after": np.round(e0_after, 3).tolist(),
        "residual_std_eV": round(residual, 3),
        "ok": bool(residual < tolerance),
        "note": f"对齐后样品边位标准差 {residual:.3f} eV",
    }


def merge_scans(scans, reject_outliers: bool = False,
                z_threshold: float = 3.5) -> RawScan:
    """
    多次扫描求算术平均（可选异常扫描剔除 + 逐点标准差）。

    数学:
        μ̄(E) = (1/N) Σ_i μ_i(E)
        σ_μ(E) = std_i μ_i(E)（保留扫描，ddof=1）

    网格一致时直接平均；不一致时线性插值到公共网格。
    **不用三次样条**——吸收边是陡跳变，样条会过冲。

    异常扫描判据（B6-1，``reject_outliers=True`` 且至少 3 条扫描）：
        d_i = RMS_E[μ_i(E) − median_j μ_j(E)]
        d_i > median(d) + z·max(1.4826·MAD(d), 1e-4)
        且 d_i > 0.02·ptp(median spectrum)
    保底护栏：至少保留 2 条最接近中位谱的扫描，绝不全部剔光。
    默认 ``reject_outliers=False`` 时 energy/mu 与原实现逐位一致，
    仅额外填充 ``sigma_mu`` 与 ``meta["merge_info"]``。
    """
    if not scans:
        raise XAFSConfigError("S1.5", "没有可合并的扫描")

    grids_equal = all(
        len(s.energy) == len(scans[0].energy)
        and np.allclose(s.energy, scans[0].energy)
        for s in scans)

    if grids_equal:
        energy = np.asarray(scans[0].energy)
        mu_stack = np.array([s.mu for s in scans])
        interpolated = False
    else:
        lo = max(s.energy.min() for s in scans)
        hi = min(s.energy.max() for s in scans)
        energy = np.arange(lo, hi, 0.2)
        mu_stack = np.array([np.interp(energy, s.energy, s.mu) for s in scans])
        interpolated = True

    n_scan = len(scans)
    # 默认路径保持与原实现完全相同的求和顺序与数值。
    mu_mean = mu_stack.mean(axis=0)

    # ---- B6-1 相异度筛查：中位谱 RMS 偏离 + MAD 门禁 ----
    deviations = np.zeros(n_scan, dtype=float)
    mu_median = None
    if n_scan >= 2:
        mu_median = np.median(mu_stack, axis=0)
        deviations = np.sqrt(np.mean((mu_stack - mu_median) ** 2, axis=1))

    kept = np.ones(n_scan, dtype=bool)
    rejected_indices: list[int] = []
    if reject_outliers and n_scan >= 3:
        d_med = float(np.median(deviations))
        mad = float(np.median(np.abs(deviations - d_med)))
        robust_sigma = max(1.4826 * mad, 1e-4)
        threshold = d_med + float(z_threshold) * robust_sigma
        ptp = float(np.ptp(mu_median)) if mu_median is not None else 0.0
        for i, d_i in enumerate(deviations):
            if d_i > threshold and d_i > 0.02 * ptp:
                kept[i] = False
        # 保底：至少保留 2 条最接近中位谱的扫描
        if int(kept.sum()) < 2:
            order = np.argsort(deviations)[:2]
            kept[:] = False
            kept[order] = True
        rejected_indices = [int(i) for i in np.flatnonzero(~kept)]
        if rejected_indices:
            mu_mean = mu_stack[kept].mean(axis=0)

    single = n_scan == 1
    rms = [float(np.sqrt(np.mean((m - mu_mean) ** 2))) for m in mu_stack]
    if int(kept.sum()) >= 2:
        sigma_mu = np.std(mu_stack[kept], axis=0, ddof=1)
        sigma_mu = np.asarray(sigma_mu, dtype=float)
    else:
        sigma_mu = None

    meta = dict(scans[0].meta)
    meta.update({
        "n_scans": n_scan,
        "grids_identical": bool(grids_equal),
        "interpolated": interpolated,
        "scan_rms": None if single else rms,
        "repeatability_note": "单扫描，无法评估重复性" if single else None,
        "merge_info": {
            "n_input": n_scan,
            "n_kept": int(kept.sum()),
            "rejected_indices": rejected_indices,
            "scan_rms_deviations": [float(v) for v in deviations],
            "z_threshold": float(z_threshold),
            "reject_outliers": bool(reject_outliers),
        },
    })
    # 插值分支下 iref 长度与新网格不匹配，保留只会诱导后续误用，直接置 None。
    # 网格一致时保留第一条的 iref（长度天然对齐）。
    iref = scans[0].iref if grids_equal else None
    return RawScan(energy=energy, mu=mu_mean, i0=None, it=None,
                   iref=iref, sigma_mu=sigma_mu, meta=meta)


def _copy_raw_scan(scan: RawScan) -> RawScan:
    """深拷贝 RawScan 的逐点数组（不共享可变状态）。"""
    return RawScan(
        energy=np.array(scan.energy, dtype=float, copy=True),
        i0=None if scan.i0 is None else np.array(scan.i0, dtype=float, copy=True),
        it=None if scan.it is None else np.array(scan.it, dtype=float, copy=True),
        ifluor=None if scan.ifluor is None else np.array(scan.ifluor,
                                                         dtype=float, copy=True),
        iref=None if scan.iref is None else np.array(scan.iref,
                                                     dtype=float, copy=True),
        dtime=None if scan.dtime is None else np.array(scan.dtime,
                                                       dtype=float, copy=True),
        mu=None if scan.mu is None else np.array(scan.mu, dtype=float, copy=True),
        sigma_mu=(None if scan.sigma_mu is None else
                  np.array(scan.sigma_mu, dtype=float, copy=True)),
        meta=dict(scan.meta or {}),
    )


def apply_deadtime_correction(scan: RawScan) -> RawScan:
    """B2-3 荧光探测器死时间校正（逐点）。

    约定（对标 Athena 的 dtime 列）：
        dtime >= 1.0   → 倍乘因子：I_corr = I_raw · dtime
        0 < dtime < 1  → 死时间分数：I_corr = I_raw / (1 − dtime)
        dtime <= 0 / 非有限值 → 该点不校正（factor=1）

    只改荧光强度并重算 μ = I_f^corr / I_0；若没有 I_f 但已有 μ，则按
    同一 factor 缩放 μ（等价于 I_f 校正后除以同一 I_0）。透射模式不做
    死时间校正（透射探测器没有荧光死时间语义）。
    """
    out = _copy_raw_scan(scan)
    mode = (scan.meta or {}).get("mode")
    if mode not in (None, "fluor"):
        # 透射/文件自带 xmu 的 dtime 列是积分时间，不是死时间因子；
        # 直接乘/除会把整条谱系统性地缩放。
        out.meta["deadtime_corrected"] = False
        out.meta["deadtime_note"] = f"mode={mode}，非荧光模式，跳过死时间校正"
        return out
    if scan.dtime is None:
        out.meta["deadtime_corrected"] = False
        out.meta["deadtime_note"] = "无 dtime 列"
        return out
    d = np.asarray(scan.dtime, dtype=float)
    if d.ndim == 0:
        d = np.full(len(scan.energy), float(d))
    if d.shape != np.asarray(scan.energy).shape:
        if d.size == 1:
            d = np.full(len(scan.energy), float(d.ravel()[0]))
        else:
            raise XAFSConfigError(
                "S1.2", f"dtime 长度 {d.shape} 与能量 {np.shape(scan.energy)} 不一致")
    factor = np.ones_like(d, dtype=float)
    mul = np.isfinite(d) & (d >= 1.0)
    frac = np.isfinite(d) & (d > 0.0) & (d < 1.0)
    factor[mul] = d[mul]
    factor[frac] = 1.0 / (1.0 - d[frac])
    n_bad = int(np.sum(~(mul | frac)))
    n_eff = int(np.sum(mul | frac))

    if scan.ifluor is not None:
        out.ifluor = np.asarray(scan.ifluor, dtype=float) * factor
        if scan.i0 is not None:
            i0 = np.asarray(scan.i0, dtype=float)
            out.mu = out.ifluor / np.where(i0 != 0, i0, np.nan)
        elif scan.mu is not None:
            out.mu = np.asarray(scan.mu, dtype=float) * factor
        else:
            raise XAFSConfigError(
                "S1.2", "死时间校正需要 I_0（重算 μ = I_f/I_0）或已有 μ")
    elif scan.mu is not None:
        out.mu = np.asarray(scan.mu, dtype=float) * factor
    else:
        raise XAFSConfigError(
            "S1.2", "死时间校正需要荧光强度 I_f 或已算好的 μ")
    out.meta["deadtime_corrected"] = bool(n_eff > 0)
    out.meta["deadtime_note"] = (f"{n_eff}/{len(d)} 点有效校正"
                                 + (f"，{n_bad} 点无有效因子" if n_bad else ""))
    out.meta["deadtime_factor_range"] = (float(np.min(factor)),
                                         float(np.max(factor)))
    return out


def deglitch_ranges(scan: RawScan, ranges: list[tuple[float, float]]):
    """B2-3 显式能量区间去 glitch：直接删点，不插值伪造数据。

    安全门禁（对标 Athena 的 deglitch 保守原则）：
      - 删除比例 > 15% → 拒绝（多半是区间写错/整段谱被掏空）；
      - 吸收边核心 [E0−10, E0+10] eV 剩余 < 5 点 → 拒绝；
      - E0 未知（meta 无 e0_nominal）时跳过核心区检查并记录。

    返回 ``(new_scan, info)``；``info`` 含 n_total/n_removed/fraction/ranges。
    """
    if not ranges:
        raise XAFSConfigError("S1.2", "deglitch 区间列表为空")
    energy = np.asarray(scan.energy, dtype=float)
    if energy.ndim != 1 or energy.size == 0:
        raise XAFSConfigError("S1.2", "deglitch 需要一维非空能量数组")
    mask = np.zeros(energy.size, dtype=bool)
    norm_ranges = []
    for item in ranges:
        try:
            lo, hi = float(item[0]), float(item[1])
        except (TypeError, ValueError, IndexError) as exc:
            raise XAFSConfigError(
                "S1.2", f"deglitch 区间格式错误: {item!r}（应 (e_lo, e_hi)）") from exc
        if not (np.isfinite(lo) and np.isfinite(hi)):
            raise XAFSConfigError("S1.2", f"deglitch 区间含非有限值: {item!r}")
        if hi < lo:
            lo, hi = hi, lo
        norm_ranges.append((lo, hi))
        mask |= (energy >= lo) & (energy <= hi)

    n_total = int(energy.size)
    n_removed = int(mask.sum())
    fraction = float(n_removed / n_total) if n_total else 0.0
    if fraction > 0.15:
        raise XAFSConfigError(
            "S1.2", f"deglitch 将删除 {n_removed}/{n_total} = {fraction:.1%} "
                    f"数据点（>15% 安全上限）；请缩小区间或检查能量单位")

    e0 = None
    if scan.meta:
        e0 = scan.meta.get("e0_nominal")
    core_remaining = None
    if e0 is not None and np.isfinite(float(e0)):
        core = (energy >= float(e0) - 10.0) & (energy <= float(e0) + 10.0)
        core_remaining = int(np.sum(core & ~mask))
        if int(core.sum()) >= 5 and core_remaining < 5:
            raise XAFSConfigError(
                "S1.2", f"deglitch 把吸收边核心 [{float(e0)-10:.1f},"
                        f"{float(e0)+10:.1f}] eV 掏空（剩余 {core_remaining} 点 <5）；"
                        f"E0 将不可靠，拒绝执行")

    keep = ~mask
    out = _copy_raw_scan(scan)
    out.energy = energy[keep]
    for attr in ("i0", "it", "ifluor", "iref", "dtime", "mu"):
        arr = getattr(scan, attr)
        if arr is not None:
            arr = np.asarray(arr, dtype=float)
            if arr.shape == energy.shape:
                setattr(out, attr, arr[keep])
            else:
                setattr(out, attr, None)
    info = {
        "n_total": n_total,
        "n_removed": n_removed,
        "fraction": fraction,
        "ranges": norm_ranges,
        "edge_core_remaining": core_remaining,
        "e0_nominal": None if e0 is None else float(e0),
    }
    out.meta["deglitch"] = info
    return out, info


def despike_scan(scan: RawScan, kernel: int = 5,
                 threshold: float = 3.0) -> tuple[RawScan, dict]:
    """S1.2 单点去毛刺（中值滑窗 + MAD 门禁）。

    注意与 ``deglitch_ranges`` 的分工：
      - ``despike_scan`` 仅适合宽度 1–2 个采样点的孤立电子学毛刺；
      - 宽度 ≥ kernel//2 + 1（如 5–15 eV 的单色器晶体 glitch）会让滑窗中值
        本身被抬高，必须用 ``deglitch_ranges`` 显式删点，严禁指望自动滤波。
    """
    out = _copy_raw_scan(scan)
    mu = np.asarray(out.mu, dtype=float).copy()
    half = max(1, int(kernel) // 2)
    win = 2 * half + 1
    padded = np.pad(mu, half, mode="edge")
    med = np.array([np.median(padded[i:i + win]) for i in range(len(mu))])
    resid = np.abs(mu - med)
    mad = float(np.median(resid)) + 1e-12
    spike_mask = resid > float(threshold) * mad
    n_spike = int(spike_mask.sum())
    if n_spike:
        mu[spike_mask] = med[spike_mask]
        out.mu = mu
    info = {"n_spikes": n_spike, "kernel": win, "threshold_mad": float(threshold)}
    out.meta["despike"] = info
    return out, info


# ================================================================
# 2. 能量校准：E0 的三种定义
# ================================================================

def find_e0_derivative_max(energy, mu, lo, hi, smooth_eV: float = 2.0) -> float:
    """
    一阶导数最大值处的能量。

    数学:
        E0 = argmax_E  dμ/dE     ,  E ∈ [lo, hi]

    实现要点:
        ① 先重采样到 0.1 eV 再求导，避免原始不等间距导致的偏大导数
        ② 平滑用 uniform_filter1d(mode="nearest")。
           **不能用 np.convolve(mode="same")**——它在边界补零，
           会把窗口起点做成一个虚假的巨大正导数峰（实测会把 E0 判到窗口下沿）。
    """
    from scipy.ndimage import uniform_filter1d

    grid = np.arange(lo, hi, 0.1)
    resampled = np.interp(grid, energy, mu)
    size = max(int(round(smooth_eV / 0.1)), 1)
    if size > 1:
        resampled = uniform_filter1d(resampled, size=size, mode="nearest")
    derivative = np.gradient(resampled, grid)
    return float(grid[int(np.argmax(derivative))])


def find_e0_half_height(energy, mu, e0_guess: float, window_eV: float = 80.0) -> float:
    """
    归一化谱 μ_norm = 0.5 处的能量（half-height 定义）。

    数学:
        μ_norm(E0) = 0.5,   其中 μ_norm = (μ − y_pre) / Δμ

    **必须先归一化再找 0.5**。μ 的绝对值由样品厚度与稀释度决定：
    γ-FeOOH 的 μ 全程在 1.36–2.42，直接判 `mu > 0.5` 会在第一个点就成立，
    插值一路外推到 9593 eV（实测）。
    与归一化的循环依赖用两段式消除：先用导数峰做初值归一化，再取半高点。
    """
    from larch.xafs import pre_edge

    provisional = Group(energy=energy, mu=mu)
    pre_edge(provisional, e0=e0_guess, pre1=-150, pre2=-30,
             norm1=150, norm2=600, nnorm=2)

    mask = ((provisional.energy >= e0_guess - window_eV)
            & (provisional.energy <= e0_guess + window_eV))
    e = np.asarray(provisional.energy)[mask]
    y = np.asarray(provisional.norm)[mask]
    if len(e) < 4:
        return float(e0_guess)

    diff = y - 0.5
    crossings = np.where((diff[:-1] < 0) & (diff[1:] >= 0))[0]
    if len(crossings) == 0:
        return float(e0_guess)
    i = int(crossings[0])
    return float(e[i] + (0.5 - y[i]) * (e[i + 1] - e[i]) / (y[i + 1] - y[i]))


def determine_e0(scan, method: str = "derivative", manual_e0: float | None = None,
                 guess: float | None = None) -> dict:
    """
    判定样品的吸收边位置，并同时给出参考通道的边位供交叉校验。

    method:
        "derivative"   一阶导数峰（默认，最常用）
        "half_height"  归一化半高点（可把 ΔE0 压低约 3 eV）
        "manual"       手工指定，带范围校验

    返回字典含:
        e0_sample            用于 k 轴零点
        e0_reference         参考通道边位，只用于判断刻度漂移，不参与 k 轴
        e0_derivative_peak   导数法结果（half_height 模式下作为对照保留）
    """
    if scan.mu is None:
        raise XAFSDataError("S2", "缺少 μ(E)，无法判定 E0")
    energy = np.asarray(scan.energy)
    mu = np.asarray(scan.mu)

    nominal = guess or scan.meta.get("e0_nominal")
    if nominal is None:
        raise XAFSConfigError(
            "S2", "无法确定吸收边搜索窗口：文件头无 e0_nominal，也未给吸收元素"
                  "（absorber 查表）或边能量 guess——元素信息是 S2 的必需输入")
    lo, hi = nominal - 30, nominal + 50

    out = {"e0_nominal": scan.meta.get("e0_nominal"), "e0_guess": nominal}
    e0_deriv = find_e0_derivative_max(energy, mu, lo, hi)
    out["e0_derivative_peak"] = e0_deriv

    method = (method or "derivative").lower()
    if method == "manual":
        if manual_e0 is None:
            raise XAFSConfigError("S2", "method='manual' 但未给出 manual_e0")
        if not (energy.min() <= manual_e0 <= energy.max()):
            raise XAFSConfigError(
                "S2", f"manual_e0={manual_e0} 超出数据范围 "
                      f"{energy.min():.0f}–{energy.max():.0f}")
        out["e0_sample"] = float(manual_e0)
        out["e0_method"] = "manual"
    elif method == "half_height":
        out["e0_sample"] = find_e0_half_height(energy, mu, e0_deriv)
        out["e0_method"] = "half_height"
    elif method == "derivative":
        out["e0_sample"] = e0_deriv
        out["e0_method"] = "derivative"
    else:
        raise XAFSConfigError(
            "S2", f"未知 method={method}，可选 derivative / half_height / manual")

    if scan.iref is not None and scan.i0 is not None:
        mu_ref = np.log(scan.i0 / scan.iref)
        e0_ref = find_e0_derivative_max(energy, mu_ref, lo, hi)
        shift = out["e0_sample"] - e0_ref
        out["e0_reference"] = e0_ref
        out["chemical_shift"] = shift
        out["reference_independent"] = bool(abs(shift) > 3.0)
        if not out["reference_independent"]:
            out["reference_warning"] = (
                f"参考通道边位与样品仅差 {shift:+.2f} eV，测的是同一吸收体"
                f"（多为样品下游电离室），不可用于能量校准")
    return out


# ================================================================
# 3. 归一化
# ================================================================

def fit_linear_pre_edge(energy, mu, e0: float, pre_window=(-150, -30)):
    """
    边前线性基线。

    数学:
        y_pre(E) = a·E + b,   a,b 由 [E0+pre_lo, E0+pre_hi] 内最小二乘确定

    返回 (a, b, y_pre 数组)
    """
    mask = ((energy >= e0 + pre_window[0]) & (energy <= e0 + pre_window[1]))
    if mask.sum() < 4:
        raise XAFSDataError(
            "S3", f"边前仅 {int(mask.sum())} 点，不足以做线性拟合；"
                  f"请放宽 pre_window 或检查数据起点")
    slope, intercept = np.polyfit(energy[mask], mu[mask], 1)
    return float(slope), float(intercept), slope * energy + intercept


def compute_edge_step(energy, mu, pre_line, e0: float, post_window=(150, 600)) -> float:
    """
    吸收阶跃高度。

    数学:
        Δμ = median[ μ(E) − y_pre(E) ] ,  E ∈ [E0+post_lo, E0+post_hi]

    用中位数而非均值，避免残余 EXAFS 振荡拉偏。
    """
    mask = ((energy >= e0 + post_window[0]) & (energy <= e0 + post_window[1]))
    if mask.sum() < 4:
        raise XAFSDataError(
            "S3", f"边后仅 {int(mask.sum())} 点，无法确定边跳变")
    return float(np.median(mu[mask] - pre_line[mask]))


def normalize_absorption(scan, e0: float, pre_window=(-150, -30),
                         post_window=(150, 600), engine: str = "larch") -> NormalizedScan:
    """
    归一化吸收谱。

    数学:
        μ_norm(E) = (μ(E) − y_pre(E)) / Δμ

    物理:
        扣掉边前基线以消除样品厚度与基底吸收；再以 Δμ 归一，
        使不同浓度样品的边跳变都等于 1，从而可横向比较。
    """
    energy = np.asarray(scan.energy)
    mu = np.asarray(scan.mu)

    if engine == "larch":
        # 与 Athena 一致：边后区域用二次多项式拟合，外推到 E0 得 edge_step。
        # native 模式改用中位数法，边后平台更接近 1（实测 1.002 vs 0.912），
        # 但 edge_step 会小约 10%，与 Athena 的数字对不上。
        from larch.xafs import pre_edge
        g = Group(energy=energy, mu=mu)
        pre_edge(g, e0=e0, pre1=pre_window[0], pre2=pre_window[1],
                 norm1=post_window[0], norm2=post_window[1], nnorm=2)
        if float(g.edge_step) <= 0:
            raise XAFSDataError(
                "S3", f"larch pre_edge 给出 Δμ = {g.edge_step:.4g} ≤ 0，"
                      f"边跳变方向异常，检查数据或改用 engine='native'")
        plateau_mask = ((energy > e0 + 200)
                        & (energy < e0 + min(post_window[1], 600)))
        plateau = float(np.mean(np.asarray(g.norm)[plateau_mask]))
        # Larch 的 g.flat 已经扣掉二次边后背景，但外推到远区的残差会让
        # 平台偏离 1（LEPIDO 实测 0.925）。Athena 的 flattened 显示约定
        # 是平台中心≈1，这里按 [E0+100, E0+400] 均值归一。
        flat = np.asarray(getattr(g, "flat", g.norm), dtype=float).copy()
        flat = _renormalize_flat_plateau(energy, flat, e0)
        return NormalizedScan(energy=np.asarray(g.energy), mu=np.asarray(g.mu),
                              mu_norm=np.asarray(g.norm),
                              pre_edge_line=np.asarray(g.pre_edge),
                              edge_step=float(g.edge_step),
                              plateau=plateau, e0=float(e0),
                              flattened=flat)

    _, _, pre_line = fit_linear_pre_edge(energy, mu, e0, pre_window)
    edge_step = compute_edge_step(energy, mu, pre_line, e0, post_window)
    if edge_step <= 0:
        raise XAFSDataError("S3", f"Δμ = {edge_step:.4g} ≤ 0，边跳变方向异常")
    mu_norm = (mu - pre_line) / edge_step

    plateau_mask = (energy > e0 + 200) & (energy < e0 + min(post_window[1], 600))
    plateau = float(np.mean(mu_norm[plateau_mask])) if plateau_mask.any() else np.nan

    flat = _flatten_post_edge(energy, mu_norm, e0, post_window)
    return NormalizedScan(energy=energy, mu=mu, mu_norm=mu_norm,
                          pre_edge_line=pre_line, edge_step=edge_step,
                          plateau=plateau, e0=float(e0), flattened=flat)


def normalize_mback(scan, e0: float, elem: str, edge: str = "K",
                    pre_window=(-200, -30), post_window=(60, 150),
                    fit_order: int = 2) -> NormalizedScan:
    """B10-1 短边后窗口的 Cromer–Liberman（MBACK）归一化。

    传统 ``normalize_absorption`` 在边后仅 60–150 eV 可用时，二次多项式
    外推会把平台拉偏（实测可差 10% 以上）。MBACK 改用一个**理论截面**
    做振幅标尺，只拟合一个比例因子和一条低阶基线：

        μ_meas(E) ≈ a·μ_th(E) + b0 + b1·(E−E0) + …

    其中 μ_th 来自 ``xraydb.mu_elam``（Cromer–Liberman 参数化）。
    归一化:

        μ_norm(E) = (μ_meas(E) − B(E)) / (a·Δμ_th)
        Δμ_th = μ_th(E0+50) − μ_th(E0−50)

    边界条件（不伪造）:
      - 边前/边后窗口各至少 3 点，否则 ``XAFSDataError``；
      - ``a <= 0`` 或 ``a·Δμ_th <= 0``：理论与实测边跳方向矛盾，拒绝。

    ``fit_order`` 是基线多项式阶数（1=常数，2=线性，默认 2）。
    """
    energy = np.asarray(scan.energy, dtype=float)
    mu = np.asarray(scan.mu, dtype=float)
    if energy.shape != mu.shape:
        raise XAFSDataError(
            "S3", f"MBACK 输入长度不一致: E={energy.shape}, μ={mu.shape}")
    e0 = float(e0)
    elem = str(elem)
    mu_th = np.asarray(mu_elam(elem, energy), dtype=float)
    if energy.size < 8 or not np.all(np.isfinite(mu_th)):
        raise XAFSDataError("S3", f"MBACK 理论截面不可用（元素 {elem}）")

    pre = ((energy >= e0 + float(pre_window[0]))
           & (energy <= e0 + float(pre_window[1])))
    post = ((energy >= e0 + float(post_window[0]))
            & (energy <= e0 + float(post_window[1])))
    if int(pre.sum()) < 3 or int(post.sum()) < 3:
        raise XAFSDataError(
            "S3", f"MBACK 窗口点数不足: 边前 {int(pre.sum())}、"
                  f"边后 {int(post.sum())}（各需 ≥3）")

    mask = pre | post
    x = energy[mask] - e0
    order = max(1, int(fit_order))
    # 设计矩阵：[μ_th, 1, ΔE, ΔE², ...]，基线阶数 = order
    cols = [mu_th[mask]] + [x ** p for p in range(order)]
    design = np.column_stack(cols)
    coef, _, rank, _ = np.linalg.lstsq(design, mu[mask], rcond=None)
    if rank < design.shape[1]:
        raise XAFSDataError(
            "S3", f"MBACK 设计矩阵秩亏（rank={rank}/{design.shape[1]}），"
                  f"窗口内理论截面与多项式共线")
    a = float(coef[0])
    base_coef = np.asarray(coef[1:], dtype=float)
    baseline = np.zeros_like(energy)
    for p, c in enumerate(base_coef):
        baseline += c * (energy - e0) ** p

    if not np.isfinite(a) or a <= 0:
        raise XAFSDataError(
            "S3", f"MBACK 振幅因子 a={a:.4g} ≤ 0：实测谱与 {elem} 理论截面"
                  f"边跳方向不一致，拒绝归一化")

    jump_pair = np.asarray(mu_elam(elem, np.array([e0 + 50.0, e0 - 50.0])),
                           dtype=float)
    jump = float(jump_pair[0] - jump_pair[1])
    edge_step = a * jump
    if not np.isfinite(edge_step) or edge_step <= 0:
        raise XAFSDataError(
            "S3", f"MBACK 边跳 a·Δμ_th={edge_step:.4g} ≤ 0，拒绝归一化")

    mu_norm = (mu - baseline) / edge_step
    plateau_mask = ((energy >= e0 + float(post_window[0]))
                    & (energy <= e0 + float(post_window[1])))
    plateau = (float(np.mean(mu_norm[plateau_mask]))
               if plateau_mask.any() else float("nan"))
    flat = _flatten_post_edge(
        energy, mu_norm, e0,
        (float(post_window[0]), max(float(post_window[1]), 300.0)))
    mback = {
        "a": a,
        "b": float(base_coef[0]) if base_coef.size else 0.0,
        "baseline_coeffs": [float(c) for c in base_coef],
        "theory_edge_step": float(jump),
        "edge_step": float(edge_step),
        "post_window": [float(post_window[0]), float(post_window[1])],
        "pre_window": [float(pre_window[0]), float(pre_window[1])],
        "element": elem,
        "edge": str(edge),
        "fit_order": order,
    }
    return NormalizedScan(energy=energy, mu=mu, mu_norm=mu_norm,
                          pre_edge_line=baseline, edge_step=float(edge_step),
                          plateau=plateau, e0=e0, flattened=flat,
                          engine="mback", mback=mback)


def _renormalize_flat_plateau(energy, flat, e0: float):
    """把 flattened 谱的边后平台中心归一到 1.0（E<E0 不动）。"""
    energy = np.asarray(energy, dtype=float)
    out = np.asarray(flat, dtype=float).copy()
    sel = energy >= e0
    for lo, hi in ((100.0, 400.0), (30.0, 600.0)):
        mask = (energy >= e0 + lo) & (energy <= e0 + hi)
        if int(mask.sum()) >= 5:
            scale = float(np.mean(out[mask]))
            if np.isfinite(scale) and scale > 1e-6:
                out[sel] = out[sel] / scale
            break
    return out


def _flatten_post_edge(energy, mu_norm, e0: float,
                       post_window=(150, 600)):
    """native 分支配平：边后二次背景扣除 + 平台归一。

    E≥E0 时 ``flat = μ_norm − [P(E) − P(E0)]``，P 是拟合到 μ_norm 的
    二次多项式；E<E0 保持 μ_norm。之所以对 μ_norm 而非原始 μ 拟合，
    是因为原始 μ 的二次外推到 E0 会把平台抬到 1.5 左右（实测），
    与 Athena flattened 的“平台≈1”约定不符。
    """
    energy = np.asarray(energy, dtype=float)
    mu_norm = np.asarray(mu_norm, dtype=float)
    out = mu_norm.copy()
    sel = energy >= e0
    mask = ((energy >= e0 + post_window[0])
            & (energy <= e0 + post_window[1]))
    if int(mask.sum()) >= 5:
        coef = np.polyfit(energy[mask], mu_norm[mask], 2)
        post = np.polyval(coef, energy)
        base = float(np.polyval(coef, e0))
        out[sel] = mu_norm[sel] - (post[sel] - base)
    return _renormalize_flat_plateau(energy, out, e0)


def fit_whiteline_area(norm: NormalizedScan,
                       window: tuple[float, float] = (-6.0, 40.0)) -> dict:
    """B9-2 白线净面积：arctan 阶跃 + Pseudo-Voigt。

    K 边白线常与边前峰/主边台阶叠在一起；L2/L3 边没有 1s→3d 四极
    边前峰，主特征就是白线。统一模型:

        μ(E) = b0 + b1(E−E0)
               + h·[1/2 + arctan((E−E0)/w)/π]      （主边阶跃）
               + PV(E; A, E_wl, γ, η)               （白线，A=积分面积）

    拟合后 ``A`` 即扣掉阶跃后的净白线面积（μ_norm 已按边跳=1 归一，
    故单位是 eV；``white_line_area`` 乘回 ``edge_step`` 得到原始 μ·eV）。
    窗口点数 <8 或拟合失败时返回 ``reliable=False`` 与零面积，不伪造。
    """
    from scipy.optimize import least_squares

    energy = np.asarray(norm.energy, dtype=float)
    y_all = (np.asarray(norm.flattened, dtype=float)
             if norm.flattened is not None
             else np.asarray(norm.mu_norm, dtype=float))
    e0 = float(norm.e0)
    emin, emax = e0 + float(window[0]), e0 + float(window[1])
    sel = (energy >= emin) & (energy <= emax)
    E, y = energy[sel], y_all[sel]
    out = {
        "white_line_area": 0.0,
        "white_line_area_norm": 0.0,
        "whiteline_peak_height": 0.0,
        "white_line_energy": None,
        "r_factor": None,
        "reliable": False,
        "window_eV": [float(emin), float(emax)],
        "note": "",
    }
    if E.size < 8:
        out["note"] = f"白线窗口仅 {E.size} 点（<8），安全降级"
        return out

    y_max = float(np.max(y))
    y_min = float(np.min(y))
    span = max(y_max - y_min, 1e-9)
    edge0 = float(np.median(y[:max(3, E.size // 10)]))
    # 初值：阶跃高度 ≈ 远端平台 − 近端基线；白线中心取窗口内最大值。
    tail = float(np.median(y[-max(3, E.size // 10):]))
    h0 = max(tail - edge0, 0.05 * span)
    i_max = int(np.argmax(y))
    A0 = max(0.2 * span * 3.0, 1e-6)
    x0 = np.array([edge0, 0.0, h0, 1.5, A0, float(E[i_max]), 1.5, 0.5])

    def _model(p):
        b0, b1, h, w, A, ec, gam, eta = p
        step = h * (0.5 + np.arctan((E - e0) / max(w, 1e-6)) / np.pi)
        return (b0 + b1 * (E - e0) + step
                + _pseudo_voigt_area(E, A, ec, gam, eta))

    lb = [-np.inf, -np.inf, 0.0, 0.2, 0.0, e0 - 8.0, 0.3, 0.0]
    ub = [np.inf, np.inf, 3.0 * span, 15.0, 50.0 * span * max(emax - emin, 1.0),
          e0 + 20.0, 10.0, 1.0]
    try:
        sol = least_squares(lambda p: _model(p) - y, np.clip(x0, lb, ub),
                            bounds=(lb, ub), max_nfev=40000)
        p = np.asarray(sol.x, dtype=float)
        resid = _model(p) - y
    except Exception as exc:  # noqa: BLE001 - 拟合失败必须可解释地降级
        out["note"] = f"白线拟合失败: {exc}"
        return out

    b0, b1, h, w, A, ec, gam, eta = [float(v) for v in p]
    if not (np.all(np.isfinite(p)) and A >= 0 and h >= 0):
        out["note"] = "白线拟合得到非物理参数（面积/阶跃为负或非有限）"
        return out
    pv_at_center = float(_pseudo_voigt_area(
        np.array([ec]), A, ec, gam, eta)[0])
    denom = float(np.sum((y - np.mean(y)) ** 2))
    out.update({
        "white_line_area": float(A * float(norm.edge_step)),
        "white_line_area_norm": float(A),
        "whiteline_peak_height": pv_at_center,
        "white_line_energy": ec,
        "r_factor": (float(np.sum(resid ** 2) / denom)
                     if denom > 0 else None),
        "reliable": True,
        "note": "",
        "edge_step_fit": h,
        "edge_width_eV": w,
    })
    return out


def extract_xanes_features(norm: NormalizedScan,
                           e0_nominal: float | None = None,
                           edge: str = "K",
                           absorber: str | None = None) -> dict:
    """B2-1 XANES 定量特征 + 颈部肩峰/导数双峰指纹。

    输入优先用 ``norm.flattened``（边后展平谱）；缺失时退回 ``mu_norm``。
    所有能量窗口以本次样品的 E0 为锚点，输出为 JSON 可序列化的标量。

    颈部判据（对标 Athena 的 derivative doublet 预警）：
      1. 在上升沿 E∈[E0−5, E_WL] 且 μ_flat∈[0.15, 0.85] 内算平滑一阶导；
      2. 找 ≥2 个局部极大值（间距 ≥1.5 eV，次峰 ≥0.35×主峰）；
      3. 次峰 ≥0.70×主峰 → ``e0_jump_risk=True``：E0 可能在两个台阶间
         跳变 3–6 eV，后续归一化/FT 相位会随之漂移。
    """
    energy = np.asarray(norm.energy, dtype=float)
    e0 = float(norm.e0)
    mu_flat = norm.flattened if norm.flattened is not None else norm.mu_norm
    mu_flat = np.asarray(mu_flat, dtype=float)
    if energy.shape != mu_flat.shape:
        raise XAFSConfigError(
            "S3", f"XANES 特征提取长度不一致: E={energy.shape}, μ={mu_flat.shape}")

    out = {
        "e_half_height": None,
        "half_height_shift": None,
        "white_line_energy": None,
        "white_line_height": None,
        "pre_edge_peak_energy": None,
        "pre_edge_peak_height": 0.0,
        "pre_edge_area": 0.0,
        "pre_edge_centroid": None,
        "pre_edge_voigt_area": 0.0,
        "pre_edge_n_peaks": 0,
        "pre_edge_voigt_r_factor": None,
        "neck_shoulder_present": False,
        "neck_shoulder_energy": None,
        "derivative_split_eV": 0.0,
        "e0_jump_risk": False,
        "derivative_peaks_eV": [],
        "derivative_peak_heights": [],
    }

    # ---- 白线（flattened 谱）----
    wl_mask = (energy >= e0 - 2.0) & (energy <= e0 + 30.0)
    if wl_mask.any():
        idx = int(np.argmax(mu_flat[wl_mask]))
        out["white_line_energy"] = float(energy[wl_mask][idx])
        out["white_line_height"] = float(mu_flat[wl_mask][idx])

    # ---- 半高能量 E1/2 ----
    hh_mask = (energy >= e0 - 15.0) & (energy <= e0 + 20.0)
    E_h, mu_h = energy[hh_mask], mu_flat[hh_mask]
    above = np.flatnonzero(mu_h >= 0.5)
    if above.size:
        j = int(above[0])
        if j == 0:
            e_half = float(E_h[0])
        else:
            y0, y1 = float(mu_h[j - 1]), float(mu_h[j])
            if abs(y1 - y0) < 1e-12:
                e_half = float(E_h[j])
            else:
                frac = (0.5 - y0) / (y1 - y0)
                e_half = float(E_h[j - 1] + frac * (E_h[j] - E_h[j - 1]))
        out["e_half_height"] = e_half
        if e0_nominal is not None:
            out["half_height_shift"] = float(e_half - e0_nominal)

    edge_type = str(edge or "K").strip().upper()
    out["edge_type"] = edge_type
    if absorber is not None:
        out["absorber"] = str(absorber)
    is_l_edge = edge_type.startswith("L")

    if is_l_edge:
        # B9-2：L2/L3 边没有 1s→3d 四极边前峰，跳过边前峰分峰，
        # 改为拟合主边阶跃 + 白线，输出净面积。字段名与 K 边不混用。
        wl = fit_whiteline_area(norm)
        for key in [k for k in out if k.startswith("pre_edge_")]:
            out.pop(key, None)
        out.update({
            "white_line_area": float(wl.get("white_line_area") or 0.0),
            "white_line_area_norm": float(wl.get("white_line_area_norm") or 0.0),
            "whiteline_peak_height": float(wl.get("whiteline_peak_height") or 0.0),
            "white_line_fit_r_factor": wl.get("r_factor"),
            "white_line_fit_reliable": bool(wl.get("reliable")),
            "white_line_fit_note": wl.get("note", ""),
        })
        if wl.get("white_line_energy") is not None:
            out["white_line_energy"] = float(wl["white_line_energy"])
            out["white_line_height"] = float(
                wl.get("whiteline_peak_height") or out["white_line_height"] or 0.0)
        return out

    # ---- 边前峰（扣两端点连线基线）----
    pre_mask = (energy >= e0 - 12.0) & (energy <= e0 - 2.0)
    if int(pre_mask.sum()) >= 3:
        Ep, mp = energy[pre_mask], mu_flat[pre_mask]
        if Ep[-1] > Ep[0]:
            baseline = mp[0] + (mp[-1] - mp[0]) * (Ep - Ep[0]) / (Ep[-1] - Ep[0])
        else:
            baseline = np.zeros_like(mp)
        resid = mp - baseline
        k = int(np.argmax(resid))
        out["pre_edge_peak_energy"] = float(Ep[k])
        out["pre_edge_peak_height"] = float(max(0.0, resid[k]))
        positive = np.maximum(resid, 0.0)
        trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")
        out["pre_edge_area"] = float(max(0.0, trapz(positive, Ep)))

    # ---- B6-3 边前峰分峰（扣主边拖尾 + Pseudo-Voigt）----
    try:
        vor = fit_pre_edge_peaks(norm)
    except (XAFSError, ValueError, FloatingPointError):
        vor = None
    if vor and vor.get("n_peaks", 0) >= 1 \
            and vor.get("centroid_energy") is not None:
        out["pre_edge_centroid"] = float(vor["centroid_energy"])
        out["pre_edge_voigt_area"] = float(vor["total_area"])
        out["pre_edge_n_peaks"] = int(vor["n_peaks"])
        out["pre_edge_voigt_r_factor"] = vor.get("r_factor")
    else:
        # 窗口点数不足/拟合失败：安全降级到两点直线基线的峰位
        out["pre_edge_centroid"] = out["pre_edge_peak_energy"]

    # ---- 颈部肩峰：上升沿一阶导数双峰 ----
    wl_e = out["white_line_energy"]
    rise_hi = float(wl_e) if wl_e is not None else e0 + 30.0
    rise_mask = ((energy >= e0 - 5.0) & (energy <= rise_hi)
                 & (mu_flat >= 0.15) & (mu_flat <= 0.85))
    Er, mr = energy[rise_mask], mu_flat[rise_mask]
    if len(Er) >= 7:
        # 按能量排序（文件偶尔有重复/倒序点）
        order = np.argsort(Er)
        Er, mr = Er[order], mr[order]
        dE = float(np.median(np.diff(Er))) if len(Er) > 1 else 0.0
        if np.isfinite(dE) and dE > 0:
            D1 = np.gradient(mr, Er)
            from scipy.signal import find_peaks, savgol_filter
            npts = int(round(1.5 / dE))
            npts = max(5, npts)
            if npts % 2 == 0:
                npts += 1
            if npts > len(Er):
                npts = len(Er) if len(Er) % 2 == 1 else len(Er) - 1
            if npts >= 5 and npts > 3:
                D1 = savgol_filter(D1, npts, 3)
            peaks, _ = find_peaks(D1, distance=max(1, int(1.5 / dE)))
            if len(peaks):
                heights = D1[peaks]
                main_h = float(np.max(heights))
                if main_h > 0:
                    sig = [int(p) for p, h in zip(peaks, heights)
                           if float(h) >= 0.35 * main_h]
                    out["derivative_peaks_eV"] = [
                        float(Er[p]) for p in sorted(sig, key=lambda q: Er[q])]
                    out["derivative_peak_heights"] = [
                        float(D1[p]) for p in sorted(sig, key=lambda q: Er[q])]
                    if len(sig) >= 2:
                        top = sorted(sig, key=lambda p: D1[p],
                                     reverse=True)[:2]
                        e_sorted = sorted(float(Er[p]) for p in top)
                        h_sorted = sorted(float(D1[p]) for p in top)
                        split = e_sorted[1] - e_sorted[0]
                        if split >= 1.5:
                            out["neck_shoulder_present"] = True
                            out["neck_shoulder_energy"] = e_sorted[0]
                            out["derivative_split_eV"] = float(split)
                            hmax = max(h_sorted)
                            out["e0_jump_risk"] = bool(
                                hmax > 0 and min(h_sorted) >= 0.70 * hmax)
    return out


# ================================================================
# 4. 背景扣除
# ================================================================

def _pseudo_voigt_area(E, A: float, E_c: float, gamma: float,
                       eta: float) -> np.ndarray:
    """面积参数化 Pseudo-Voigt（A 为积分面积，单位同 y·eV）。"""
    E = np.asarray(E, dtype=float)
    x = E - float(E_c)
    lor = (float(gamma) / np.pi) / (x ** 2 + float(gamma) ** 2)
    gau = (np.sqrt(np.log(2.0)) / (float(gamma) * np.sqrt(np.pi))
           * np.exp(-np.log(2.0) * (x / float(gamma)) ** 2))
    return float(A) * (float(eta) * lor + (1.0 - float(eta)) * gau)


def fit_pre_edge_peaks(
        norm: NormalizedScan,
        window: tuple[float, float] = (-14.0, -2.0),
        max_peaks: int = 2,
) -> dict:
    """边前峰主边拖尾基线扣除 + 多 Pseudo-Voigt 分峰（Wilke 质心）。

    为什么两点直线基线不够：边前峰常贴近主吸收边上升沿，主边的
    Lorentzian 拖尾在窗口内是弯曲的；直线基线会过度切割边前峰，
    使 ``pre_edge_area`` 系统性偏小，也无法拆分 t2g/eg 双峰。

    基线:
        B(E) = a + b·(E−E0) + c/((E−(E0+2))² + 3²),  c ≥ 0
        仅用窗口两端各 1.5 eV 的边界带拟合，避免把边前峰本身当基线。
    峰形:
        面积参数化 Pseudo-Voigt（A=积分面积，γ=半高半宽，η=洛伦兹权重）。
    双峰判据:
        RSS2 < 0.70·RSS1 且峰间距 ≥0.8 eV 且次峰面积 ≥10% 主峰面积。

    返回:
        centroid_energy / total_area / n_peaks / peaks / r_factor /
        baseline / residual / reliable
        面积加权质心 E_c = Σ A_k E_k / Σ A_k；分母为 0 时回退到
        最大峰中心并置 ``reliable=False``。
    """
    from scipy.optimize import least_squares

    energy = np.asarray(norm.energy, dtype=float)
    y_all = (np.asarray(norm.flattened, dtype=float)
             if norm.flattened is not None
             else np.asarray(norm.mu_norm, dtype=float))
    e0 = float(norm.e0)
    emin, emax = e0 + float(window[0]), e0 + float(window[1])
    if emax <= emin:
        raise XAFSConfigError(
            "S3", f"边前窗口无效: [{window[0]}, {window[1]}] eV")
    sel = (energy >= emin) & (energy <= emax)
    E = energy[sel]
    y = y_all[sel]
    out = {
        "centroid_energy": None, "total_area": 0.0, "n_peaks": 0,
        "peaks": [], "r_factor": None, "baseline": None, "residual": None,
        "reliable": False, "window_eV": [float(emin), float(emax)],
        "note": "",
    }
    if E.size < 6:
        out["note"] = f"边前窗口仅 {E.size} 点（<6），安全降级"
        return out

    # ---- 1) 主边拖尾基线：两端边界带联合拟合 ----
    band = (((E >= emin) & (E <= min(emin + 1.5, emax)))
            | ((E >= max(emax - 1.5, emin)) & (E <= emax)))
    Eb, yb = E[band], y[band]
    e_edge_tail = e0 + 2.0
    gamma0 = 3.0

    def _baseline(p, x):
        return p[0] + p[1] * (x - e0) + p[2] / ((x - e_edge_tail) ** 2
                                                + gamma0 ** 2)

    if Eb.size >= 4:
        x0 = np.array([float(np.median(yb)), 0.0, 0.0])
        try:
            bl = least_squares(
                lambda p: _baseline(p, Eb) - yb, x0,
                bounds=([-np.inf, -np.inf, 0.0], [np.inf, np.inf, np.inf]),
                max_nfev=20000)
            base = _baseline(bl.x, E)
            if not np.all(np.isfinite(base)):
                raise ValueError("baseline non-finite")
        except Exception:
            coef = np.polyfit(Eb, yb, 1)
            base = np.polyval(coef, E)
    else:
        coef = np.polyfit(E, y, 1)
        base = np.polyval(coef, E)
    y_net = np.maximum(y - base, 0.0)
    out["baseline"] = base
    if float(np.max(y_net)) <= 1e-12:
        out["note"] = "扣基线后无边前峰信号"
        return out

    dE = float(np.median(np.diff(E))) if E.size > 1 else 0.0
    if not np.isfinite(dE) or dE <= 0:
        out["note"] = "能量网格非单调，无法分峰"
        return out

    ec_lo, ec_hi = emin + 0.5, emax - 0.5
    a_max = max(float(np.trapezoid(y_net, E)) * 5.0, 1e-6)

    def _fit_n(n_peaks, starts):
        best = None
        for x0 in starts:
            def fun(p):
                model = np.zeros_like(E)
                for k in range(n_peaks):
                    off = 4 * k
                    model = model + _pseudo_voigt_area(
                        E, p[off], p[off + 1], p[off + 2], p[off + 3])
                return model - y_net
            lb, ub = [], []
            for _ in range(n_peaks):
                lb.extend([0.0, ec_lo, 0.4, 0.0])
                ub.extend([a_max, ec_hi, 2.5, 1.0])
            try:
                sol = least_squares(fun, np.clip(x0, lb, ub),
                                    bounds=(lb, ub), max_nfev=40000)
            except Exception:
                continue
            rss = float(np.sum(sol.fun ** 2))
            if best is None or rss < best[0]:
                best = (rss, np.asarray(sol.x, dtype=float))
        return best

    # 单峰初值：面积/峰高/中心/宽度粗估
    i_max = int(np.argmax(y_net))
    area_total = float(np.trapezoid(y_net, E))
    gamma_init = 1.0
    x1 = [area_total, float(E[i_max]), gamma_init, 0.5]
    fit1 = _fit_n(1, [x1])
    if fit1 is None:
        out["note"] = "单峰拟合未收敛"
        return out
    rss1, p1 = fit1

    best_n, best_p, best_rss = 1, p1, rss1
    if max_peaks >= 2:
        # 第二峰初值：距主峰 ≥0.8 eV 的最强局部极大（无则取最强点）
        from scipy.signal import find_peaks
        peak_idx, _ = find_peaks(y_net, prominence=1e-12)
        c_main = float(p1[1])
        cand = [int(i) for i in peak_idx
                if abs(float(E[int(i)]) - c_main) >= 0.8]
        if cand:
            c_second = float(E[max(cand, key=lambda i: y_net[i])])
        else:
            far = np.abs(E - c_main) >= 0.8
            c_second = (float(E[far][int(np.argmax(y_net[far]))])
                        if np.any(far) else c_main + 1.0)
        h_main = float(np.interp(c_main, E, y_net))
        h_second = float(np.interp(c_second, E, y_net))
        frac = h_second / (h_main + h_second) if (h_main + h_second) > 0 else 0.5
        area_second = area_total * max(0.05, min(0.95, frac))
        area_main = max(area_total - area_second, 1e-9)
        starts = [
            [area_main, c_main, 1.0, 0.5, area_second, c_second, 1.0, 0.5],
            [area_second, c_second, 1.0, 0.5, area_main, c_main, 1.0, 0.5],
            [area_main * 0.7, c_main, 0.8, 0.5,
             area_second * 1.3, c_second, 1.2, 0.5],
        ]
        fit2 = _fit_n(2, starts)
        if fit2 is not None:
            rss2, p2 = fit2
            centers = [float(p2[1]), float(p2[5])]
            areas = [abs(float(p2[0])), abs(float(p2[4]))]
            sep = abs(centers[0] - centers[1])
            ratio = min(areas) / max(max(areas), 1e-12)
            if rss2 < 0.70 * rss1 and sep >= 0.8 and ratio >= 0.10:
                best_n, best_p, best_rss = 2, p2, rss2

    peaks = []
    for k in range(best_n):
        off = 4 * k
        peaks.append({"area": float(abs(best_p[off])),
                      "center": float(best_p[off + 1]),
                      "hwhm": float(best_p[off + 2]),
                      "eta": float(best_p[off + 3])})
    peaks.sort(key=lambda p: p["center"])
    model = np.zeros_like(E)
    for p in peaks:
        model = model + _pseudo_voigt_area(
            E, p["area"], p["center"], p["hwhm"], p["eta"])
    denom_rss = float(np.sum(y_net ** 2))
    total_area = float(sum(p["area"] for p in peaks))
    if total_area > 0:
        centroid = float(sum(p["area"] * p["center"] for p in peaks)
                         / total_area)
        reliable = True
        note = ""
    else:
        centroid = float(peaks[0]["center"]) if peaks else None
        reliable = False
        note = "分峰总面积为 0，质心回退到最大峰中心"
    out.update({
        "centroid_energy": centroid,
        "total_area": total_area,
        "n_peaks": int(best_n),
        "peaks": peaks,
        "r_factor": (float(best_rss / denom_rss)
                     if denom_rss > 0 else None),
        "model": model,
        "residual": y_net - model,
        "reliable": bool(reliable),
        "note": note,
    })
    return out


def wave_number(energy, e0: float) -> np.ndarray:
    """
    光电子波矢。

    数学:
        k = sqrt(2·m_e·(E − E0) / ħ²) = sqrt(K2·(E − E0)) ≈ 0.5123·sqrt(E − E0)

    单位: E 用 eV，k 得 Å⁻¹。E < E0 处返回 0。
    """
    return np.sqrt(np.maximum(K2 * (np.asarray(energy) - e0), 0.0))


def extract_chi_k(norm: NormalizedScan, rbkg: float = 1.2, kweight: int = 2,
                  chi_std=None, k_std=None, clamp_lo: float | None = None,
                  clamp_hi: float | None = None) -> ChiK:
    """
    扣除原子背景，提取 EXAFS 振荡。

    数学:
        μ_norm(k) = χ(k) + μ0(k)
        约束 μ0 在 R 空间不含 R > Rbkg 的成分，
        实现为  χ = F⁻¹{ R > Rbkg } F[μ_norm]

    Rbkg 的物理含义:
        ≈ 0.5–0.7 × 第一配位壳层峰位。Rbkg 越大背景越"软"，
        R<1 Å 的残留越少，但过大会吃掉真实配位峰。

    为什么用 larch 的 autobk 而不是自己写:
        同一份数据，自研 R 空间滤波仅改过渡带宽度（0.06 → 0.25），
        第一峰结果会在 0.86 与 1.47 Å 之间跳；标准实现稳定给出 1.47。

    B10-1 先验与钳制（全部可选，默认不改变历史行为）:
        chi_std/k_std  标准谱 χ(k) 先验，保护低 R 真实信号；
        clamp_lo/hi    抑制样条两端甩尾的权重（0–1）。
        若当前 larch 版本不支持这些关键字，显式降级为无先验 AUTOBK，
        并在 ``ChiK.chi_std_used`` 保留 False（不假装用上了）。
    """
    def _run(extra: dict):
        g = Group(energy=np.asarray(norm.energy, dtype=float),
                  mu=np.asarray(norm.mu_norm, dtype=float))
        g.norm = np.asarray(norm.mu_norm, dtype=float)
        g.e0 = float(norm.e0)
        g.edge_step = 1.0
        autobk(g, rbkg=rbkg, kweight=kweight, e0=norm.e0, **extra)
        return g

    extra: dict = {}
    priors_requested = any(v is not None for v in
                           (chi_std, k_std, clamp_lo, clamp_hi))
    if chi_std is not None or k_std is not None:
        if chi_std is None or k_std is None:
            raise XAFSConfigError(
                "S5", "chi_std 与 k_std 必须成对给出（先验谱需要自己的 k 轴）")
        chi_std_arr = np.asarray(chi_std, dtype=float)
        k_std_arr = np.asarray(k_std, dtype=float)
        if chi_std_arr.ndim != 1 or k_std_arr.ndim != 1:
            raise XAFSConfigError("S5", "chi_std/k_std 必须是一维数组")
        if chi_std_arr.size != k_std_arr.size or chi_std_arr.size < 4:
            raise XAFSConfigError(
                "S5", f"chi_std/k_std 长度不一致或过短: "
                      f"{chi_std_arr.size} vs {k_std_arr.size}")
        if not (np.all(np.isfinite(chi_std_arr))
                and np.all(np.isfinite(k_std_arr))):
            raise XAFSConfigError("S5", "chi_std/k_std 含非有限值")
        extra["chi_std"] = chi_std_arr
        extra["k_std"] = k_std_arr
    if clamp_lo is not None:
        extra["clamp_lo"] = float(clamp_lo)
    if clamp_hi is not None:
        extra["clamp_hi"] = float(clamp_hi)

    degraded = False
    try:
        g = _run(extra)
    except TypeError:
        if not priors_requested:
            raise
        # 当前 larch 版本不接受先验/钳制关键字：显式降级。
        degraded = True
        g = _run({})

    return ChiK(k=np.asarray(g.k), chi=np.asarray(g.chi),
                background=np.asarray(g.bkg), mu_norm=np.asarray(g.norm),
                rbkg=float(rbkg),
                chi_std_used=bool("chi_std" in extra and not degraded),
                clamp_lo=(float(clamp_lo) if clamp_lo is not None
                          and not degraded else None),
                clamp_hi=(float(clamp_hi) if clamp_hi is not None
                          and not degraded else None))


def high_k_snr(chi_k: ChiK, sigma_mu=None, edge_step: float | None = None,
               kweight: int = 2) -> dict:
    """B9-2 由 μ(E) 噪声估计高 k 信噪比与可信 k 上限。

    噪声换算（常数近似，逐点 σ 更细时直接传数组）:

        σ_χ(k) = σ_μ / |edge_step|

    信号包络用 χ(k) 的滚动 RMS（窗口 ≈ 0.8 Å⁻¹），避免单点噪声造成
    假交叉。``snr_k < 2`` 且连续 3 点确认，记为该谱的 ``k_snr2``。

    诚实降级：没有 ``sigma_mu`` 时返回 ``snr_k=None``、``k_snr2=None``，
    并在 ``reason`` 说明，绝不编造一条“看起来很稳”的曲线。
    """
    k = np.asarray(chi_k.k, dtype=float)
    chi = np.asarray(chi_k.chi, dtype=float)
    out = {"snr_k": None, "k_snr2": None, "reliable_kmax": None,
           "sigma_chi": None, "reason": "", "window_A": None}
    if sigma_mu is None:
        out["reason"] = ("未提供 sigma_mu：无法估计 μ(E) 噪声，"
                         "不生成 SNR 判据")
        return out
    if edge_step is None:
        edge_step = 1.0
        out["reason"] = ("未提供 edge_step：按已归一化谱 edge_step=1 解释 "
                         "sigma_mu")
    edge_step = float(edge_step)
    if not np.isfinite(edge_step) or abs(edge_step) < 1e-12:
        raise XAFSConfigError(
            "S5", f"edge_step={edge_step:.4g} 不能用于噪声换算")
    sig = np.asarray(sigma_mu, dtype=float)
    if sig.ndim == 0:
        sigma_chi = np.full_like(k, float(sig) / abs(edge_step))
    elif sig.ndim == 1 and sig.size == k.size:
        sigma_chi = sig / abs(edge_step)
    else:
        raise XAFSConfigError(
            "S5", f"sigma_mu 需为标量或与 k 等长（{k.size}），"
                  f"收到 shape={sig.shape}；不猜测能量网格重采样")
    if not np.all(np.isfinite(sigma_chi)) or np.any(sigma_chi <= 0):
        raise XAFSConfigError("S5", "sigma_mu 必须全部为有限正数")
    out["sigma_chi"] = float(np.median(sigma_chi))

    dk = float(np.median(np.diff(k))) if k.size > 1 else 0.05
    if not np.isfinite(dk) or dk <= 0:
        raise XAFSConfigError("S5", "k 网格非单调，无法计算滚动 SNR")
    nwin = int(round(0.8 / dk))
    nwin = max(3, nwin if nwin % 2 == 1 else nwin + 1)
    nwin = min(nwin, k.size if k.size % 2 == 1 else k.size - 1)
    if nwin < 3:
        raise XAFSConfigError("S5", f"k 点数 {k.size} 太少（需 ≥4）")
    out["window_A"] = float(nwin * dk)
    half = nwin // 2
    padded = np.pad(chi, half, mode="edge")
    local_ms = np.array([
        np.mean(padded[i:i + nwin] ** 2) for i in range(k.size)])
    envelope = np.sqrt(np.maximum(local_ms, 0.0))
    snr = envelope / sigma_chi
    out["snr_k"] = snr

    below = snr < 2.0
    k_snr2 = None
    for i in range(max(0, below.size - 2)):
        if below[i] and below[i + 1] and below[i + 2]:
            k_snr2 = float(k[i])
            break
    out["k_snr2"] = k_snr2
    out["reliable_kmax"] = (k_snr2 if k_snr2 is not None
                            else float(k.max()))
    if k_snr2 is None and not out["reason"]:
        out["reason"] = "全 k 区间滚动 SNR 均 ≥2"
    return out


def merge_mca_channels(channels, dtime=None, z_threshold: float = 3.5,
                       min_channels: int = 2):
    """B9-3 多通道荧光探测器（MCA/多电离室）坏道剔除与加权合并。

    步骤:
      1. 逐道死时间校正（与 ``apply_deadtime_correction`` 同一约定）；
      2. 以通道间中位谱为参考，计算每道的稳健 z 分数
         ``max|I−I_med| / (1.4826·MAD)``；Bragg spike 会在少数能量点
         把 max 拉爆，因此比整道 RMS 更敏感；
      3. 按 z 从大到小剔除，但保底 ``min_channels`` 道；
      4. 以各道边阶（用 90 分位−10 分位近似）为权重加权求和——
         边阶高的通道信噪比更好；边阶全为 0 时退化为等权。

    返回 ``(summed, info)``；``info`` 含 ``rejected``、``n_kept``、
    ``weights``、``scores``。``channels`` 形状 ``(n_energy, n_channel)``。
    """
    arr = np.asarray(channels, dtype=float)
    if arr.ndim == 1:
        arr = arr[:, None]
    if arr.ndim != 2 or arr.shape[1] < 1:
        raise XAFSConfigError(
            "S1", f"MCA channels 需为 (n_energy, n_channel)，收到 {arr.shape}")
    n_energy, n_ch = arr.shape
    if not np.all(np.isfinite(arr)):
        raise XAFSConfigError("S1", "MCA channels 含非有限值")

    corr = arr.copy()
    dt = None
    if dtime is not None:
        dt = np.asarray(dtime, dtype=float)
        if dt.ndim == 1:
            if dt.size != n_ch:
                raise XAFSConfigError(
                    "S1", f"dtime 长度 {dt.size} 与通道数 {n_ch} 不一致")
            factor = np.where(dt >= 1.0, dt,
                              np.where((dt > 0) & (dt < 1), 1.0 / (1.0 - dt), 1.0))
            corr = corr * factor[None, :]
        elif dt.shape == (n_energy, n_ch):
            factor = np.where(dt >= 1.0, dt,
                              np.where((dt > 0) & (dt < 1), 1.0 / (1.0 - dt), 1.0))
            corr = corr * factor
        else:
            raise XAFSConfigError(
                "S1", f"dtime 需为 (n_ch,) 或 (n_energy, n_ch)，收到 {dt.shape}")

    med = np.median(corr, axis=1, keepdims=True)
    diff = corr - med
    mad = 1.4826 * np.median(np.abs(diff), axis=0)
    global_scale = float(np.median(mad[mad > 0])) if np.any(mad > 0) else 0.0
    if global_scale <= 0:
        global_scale = float(np.median(np.abs(diff))) or 1.0
    # 双判据：整道稳健 RMS 超标（增益/死时间异常）或少量极大离群点
    # （Bragg spike）。单用 max|diff| 会被 400 点高斯噪声的多重比较
    # 误伤——干净道 max z 也会到 3.5–4。
    rms_z = np.sqrt(np.mean(diff ** 2, axis=0)) / global_scale
    max_z = np.max(np.abs(diff), axis=0) / global_scale
    n_out = np.sum(np.abs(diff) > float(z_threshold) * global_scale, axis=0)
    bad = (rms_z > float(z_threshold)) | (
        (max_z > 2.0 * float(z_threshold)) & (n_out >= 2))
    scores = rms_z

    order = [int(i) for i in np.argsort(scores)[::-1]]
    rejected: list[int] = []
    keep = list(range(n_ch))
    for idx in order:
        if len(keep) <= max(1, int(min_channels)):
            break
        if bad[idx]:
            rejected.append(int(idx))
            keep.remove(int(idx))
    rejected.sort()

    spans = (np.percentile(corr[:, keep], 90, axis=0)
             - np.percentile(corr[:, keep], 10, axis=0))
    spans = np.maximum(spans, 0.0)
    if float(np.sum(spans)) > 0:
        weights = spans / float(np.sum(spans))
    else:
        weights = np.full(len(keep), 1.0 / len(keep))
    summed = np.sum(corr[:, keep] * weights[None, :], axis=1)

    info = {
        "rejected": rejected,
        "n_kept": len(keep),
        "kept": keep,
        "weights": [float(w) for w in weights],
        "scores": [float(s) for s in scores],
        "max_z": [float(s) for s in max_z],
        "outlier_counts": [int(v) for v in n_out],
        "z_threshold": float(z_threshold),
        "deadtime_corrected": bool(dt is not None),
        "n_channels": int(n_ch),
    }
    return summed, info


# ================================================================
# 4.5 荧光自吸收校正（Booth & Bridges 2005）
# ================================================================

def _booth_mu_terms(formula: str, elem: str, edge: str, line: str, e0: float,
                    k, density: float):
    """Booth 公式所需的光滑吸收系数（μm⁻¹）与荧光线能量。

    返回 (mu_t, mu_a, mu_f, e_fluor, comp)：
        mu_t  总吸收（photo）
        mu_a  吸收元素贡献（photo，**已扣边前基线**）
        mu_f  荧光能量处的总吸收（total）
    全部已乘密度并换算为 μm⁻¹；d 用 μm 时 α·d 无量纲。

    能量锚点用**表值边能量**（xraydb），不是样品 E0：
        Demeter/Athena 的 Booth 实现同样用 k²/K2 + E_edge + 0.1。
        自吸收会把 E0 拉低（实测合成谱 7111.3 eV < Fe K 7112 eV），
        若用样品 E0，k=0 落在表值边以下，μ_a 扣基线后 ≤0，校正失败。
        以表值边为锚还能避免边能量差造成的 kink。

    为什么 μ_a 要扣边前基线：
        Booth 公式里的 χ 是相对振荡 (μ_a−μ̄_a)/μ̄_a，μ̄_a 是可振荡的
        边跳部分；归一化谱已经把边前吸收扣掉了。Athena/Demeter 的
        Booth 实现同样用 cross_section(E) − cross_section(edge−200 eV)。
        不扣基线会把 μ_a 高估约 15%（Fe K 边），修正因子系统偏大。
    """
    try:
        e_fluor = float(xray_line(elem, line).energy)
        e_edge = float(xray_edge(elem, edge).energy)
        energy = np.asarray(k, dtype=float) ** 2 / K2 + e_edge + 0.1
        comp = material_mu_components(formula, energy, density=density,
                                      kind="photo")
        comp_below = material_mu_components(formula, e_edge - 200.0,
                                            density=density, kind="photo")
        mu_f_cm = float(material_mu(formula, e_fluor, density=density,
                                    kind="total"))
    except Exception as exc:
        raise XAFSConfigError(
            "S4.5", f"无法计算 {formula} / {elem} {edge} {line} 的截面："
                    f"{exc}") from exc
    if elem not in comp["elements"]:
        raise XAFSConfigError(
            "S4.5", f"化学式 {formula} 不含吸收元素 {elem}；"
                    f"校正必须用样品总组成，不能只写活性相")

    mass = float(comp["mass"])
    mu_t_cm = np.zeros_like(energy, dtype=float)
    for el in comp["elements"]:
        stoich, atomic_mass, mu_mass = comp[el]
        mu_t_cm += (density * float(stoich) * float(atomic_mass)
                    * np.asarray(mu_mass, dtype=float) / mass)
    stoich_a, mass_a, mu_a_mass = comp[elem]
    mu_a_cm = (density * float(stoich_a) * float(mass_a)
               * np.asarray(mu_a_mass, dtype=float) / mass)
    stoich_b, mass_b, mu_b_mass = comp_below[elem]
    mu_a_below_cm = (density * float(stoich_b) * float(mass_b)
                     * float(mu_b_mass) / float(comp_below["mass"]))
    mu_a_cm = mu_a_cm - mu_a_below_cm
    if np.any(mu_a_cm <= 0):
        raise XAFSConfigError(
            "S4.5", f"{elem} 扣边前基线后 μ_a ≤ 0；检查 E0/边能量是否匹配")
    return mu_t_cm * 1e-4, mu_a_cm * 1e-4, mu_f_cm * 1e-4, e_fluor, comp


def booth_self_absorption(chi_k: ChiK, formula: str, elem: str, edge: str,
                          line: str, e0: float, thickness_um: float,
                          density: float, anginp: float = 45.0,
                          angout: float = 45.0) -> tuple[ChiK, dict]:
    """在 k 空间做荧光自吸收校正（Booth & Bridges, Phys. Scr. 2005）。

    物理:
        荧光强度 ∝ μ_a(E) / [μ_T(E) + g·μ_f] × [1 − exp(−(μ_T/sinφ
        + μ_f/sinθ)d)]。测量得到的 χ_exp 因此被样品自身的吸收压低，
        且压低程度随 k 变化。Booth & Bridges 把上式对 χ 反解成闭式
        （Eq. 4），得到 χ(k) = f(χ_exp(k), μ_T, μ_a, μ_f, d, φ, θ)。

    输入（缺一不可）:
        formula      样品**总**化学式（含基体/稀释剂），如 Fe3O4、
                     ((NH4)2SO4)0.1(H2O)54.8；不是只写活性相
        elem         吸收元素（Fe/Pt/...）
        edge, line   吸收边与荧光线（K→Ka，L3→La）
        thickness_um 样品厚度（μm）；厚样可取大于信息深度的值
        density      有效密度（g/cm³）；粉末压片用几何密度，不是晶体密度
        anginp/out   入射/出射角（相对样品表面，度）

    假设与边界:
        - 均匀、平表面样品；不均匀/颗粒/孔洞会让有效 d、ρ 失真
        - μ_T/μ_a 用 xraydb 的光滑截面，保留吸收边但没有 EXAFS 振荡；
          μ_a 已扣边前基线（与 Athena/Demeter 的 Booth 实现一致），
          能量锚定表值边 E_edge+0.1（不随样品 E0 漂移）
        - 只改振幅、不改相位；返回新的 ChiK，background/mu_norm 不修改
        - 薄样/厚样由 d 自动落入相应极限；d→0 时 χ_corr→χ_exp
        - 自吸收会压低/移动 E0（合成厚样实测 7111.3 vs 表值 7112 eV）；
          校正不修 E0，拟合前应用参考箔或已知边能量校准
        - 化学式/厚度/密度/角度任一不可靠，校正就不可靠；info 里给出
          correction_factor 与 reliable 标记，报告层不得静默引用
    """
    formula = (formula or "").strip()
    if not formula:
        raise XAFSConfigError("S4.5", "Booth 自吸收校正需要样品总化学式")
    if not (thickness_um > 0):
        raise XAFSConfigError("S4.5", f"厚度必须 >0 μm，收到 {thickness_um}")
    if not (density > 0):
        raise XAFSConfigError("S4.5", f"密度必须 >0 g/cm³，收到 {density}")
    for name, ang in (("anginp", anginp), ("angout", angout)):
        if not (0.0 < ang <= 90.0):
            raise XAFSConfigError(
                "S4.5", f"{name}={ang}° 必须在 (0, 90] 内（相对样品表面；"
                        f"90° = 垂直入射/出射）")

    mu_t, mu_a, mu_f, e_fluor, comp = _booth_mu_terms(
        formula, elem, edge, line, e0, chi_k.k, density)
    sin_phi = np.sin(np.deg2rad(anginp))
    sin_theta = np.sin(np.deg2rad(angout))
    g = sin_phi / sin_theta
    alpha = mu_t + g * mu_f
    if np.any(alpha <= 0):
        raise XAFSDataError("S4.5", "Booth 校正的 α=μ_T+g·μ_f 出现非正值，"
                                    "检查化学式/荧光线/能量范围")
    arg = alpha * float(thickness_um) / sin_phi
    exp_neg = np.exp(-np.clip(arg, 0.0, 700.0))
    gamma = 1.0 - exp_neg
    beta = mu_a * float(thickness_um) * alpha / sin_phi * exp_neg

    chi_exp = np.asarray(chi_k.chi, dtype=float)
    chi_corr = np.array(chi_exp, dtype=float, copy=True)
    thin = arg < 1e-6
    # arg>20 时 exp(-arg)<2e-9，厚样极限与 Eq.4 的差异远小于数据误差；
    # 提前切到厚样分支还能避开 β→0 时二次公式的灾难性相消。
    thick = arg > 20.0
    mid = ~(thin | thick)

    # 中间厚度：Eq. 4（β>0）
    if np.any(mid):
        b = beta[mid]
        a = alpha[mid]
        ma = mu_a[mid]
        gm = gamma[mid]
        ce = chi_exp[mid]
        term = gm * (a - ma * (ce + 1.0)) + b
        disc = term ** 2 + 4.0 * a * b * gm * ce
        disc = np.clip(disc, 0.0, None)
        chi_corr[mid] = (-term + np.sqrt(disc)) / (2.0 * b)

    # 厚样极限（d→∞，β→0）：χ = χ_exp / [1 − (μ_a/α)(χ_exp+1)]
    clipped = 0
    if np.any(thick):
        denom = 1.0 - (mu_a[thick] / alpha[thick]) * (chi_exp[thick] + 1.0)
        ok = denom > 1e-6
        clipped = int(np.count_nonzero(~ok))
        chi_corr[thick] = np.where(ok, chi_exp[thick] / np.where(ok, denom, 1.0),
                                   chi_exp[thick])

    finite = np.all(np.isfinite(chi_corr))
    ratio = np.where(np.abs(chi_exp) > 1e-8, chi_corr / chi_exp, np.nan)
    ratio_finite = ratio[np.isfinite(ratio)]
    e_edge_tab = float(xray_edge(elem, edge).energy)
    e0_offset = float(e0) - e_edge_tab
    e0_warning = None
    if abs(e0_offset) > 3.0:
        e0_warning = (f"E0={e0:.2f} eV 与表值 {edge} 边 {e_edge_tab:.2f} eV "
                      f"相差 {e0_offset:+.2f} eV；自吸收可能移动了 E0，"
                      f"拟合前应用参考箔/已知边能量校准")
    info = {
        "method": "Booth & Bridges 2005 Eq.4 (k-space)",
        "corrected": True,
        "formula": formula, "elem": elem, "edge": edge, "line": line,
        "thickness_um": float(thickness_um), "density": float(density),
        "anginp": float(anginp), "angout": float(angout), "g": float(g),
        "fluor_energy_eV": float(e_fluor),
        "energy_anchor_eV": e_edge_tab + 0.1,
        "e0_offset_eV": e0_offset,
        "mu_a_baseline": "subtracted (E_edge-200 eV)",
        "mu_f_cm": float(mu_f * 1e4),
        "mu_t_cm_range": [float(np.min(mu_t) * 1e4),
                          float(np.max(mu_t) * 1e4)],
        "mu_a_cm_range": [float(np.min(mu_a) * 1e4),
                          float(np.max(mu_a) * 1e4)],
        "correction_factor_median": (float(np.nanmedian(ratio_finite))
                                     if ratio_finite.size else None),
        "correction_factor_range": ([float(np.nanmin(ratio_finite)),
                                     float(np.nanmax(ratio_finite))]
                                    if ratio_finite.size else None),
        "thin_limit_points": int(np.count_nonzero(thin)),
        "thick_limit_points": int(np.count_nonzero(thick)),
        "thick_limit_clipped": clipped,
        "reliable": bool(finite and clipped == 0),
        "warning": e0_warning,
        "note": ("Booth 厚样/任意厚度近似；要求均匀样品与准确的总化学式、"
                 "厚度、有效密度、角度。只改振幅不改相位，也不修 E0"
                 "（自吸收会把 E0 拉低，拟合前需用参考箔校准）；化学式"
                 "偏差会直接进入校正因子，不能把校正后的 N 当无系统误差"
                 "结果。"),
    }
    return ChiK(k=np.asarray(chi_k.k), chi=chi_corr,
                background=chi_k.background, mu_norm=chi_k.mu_norm,
                rbkg=chi_k.rbkg), info


def fluo_xanes_self_absorption(
        norm: NormalizedScan,
        formula: str, elem: str, edge: str, line: str,
        thickness_um: float, density: float,
        anginp: float = 45.0, angout: float = 45.0,
) -> tuple[NormalizedScan, dict]:
    """XANES 能量域荧光自吸收反解（Tröger / FLUO，与 EXAFS Booth 隔离）。

    物理:
        浓样荧光测量的归一化 XANES 满足厚样压缩关系
            μ_exp(E) = (α+1)·μ_true(E) / (α + μ_true(E))
        闭式反解（Tröger 方程）
            μ_true(E) = α·μ_exp(E) / (α + 1 − μ_exp(E))
        其中 α = [μ_bg(E_ref) + g·μ_f] / μ_a(E_ref)，
        E_ref = E0 + 50 eV，g = sinθ_in / sinθ_out。

    为什么不能直接用 Larch fluo_corr 全谱：全谱反解会把高能 EXAFS
    放大数倍，污染 χ(k)。本函数**只改 NormalizedScan 的 XANES**
    （mu_norm / flattened），绝不生成或修改 ChiK；EXAFS 自吸收仍走
    S4.5 的 ``booth_self_absorption``。

    物理护栏:
      1. 仅对 E ≥ E0−15 eV 且 μ_exp > 0 的点应用；
      2. 分母 D = α+1−μ_exp 的下界为 0.25·(α+1)；触界即钳制并置
         ``reliable=False``，防止白线尖峰处分母趋零炸出无穷大；
      3. 返回新对象，不原地修改输入 norm。

    可靠性仍取决于总化学式、厚度、有效密度与角度；这些量不准确时
    反解出的白线高度与 LCF 权重不可当作无系统误差结果。
    """
    formula = (formula or "").strip()
    if not formula:
        raise XAFSConfigError("S4.5", "XANES 荧光自吸收校正需要样品总化学式")
    if not (thickness_um > 0):
        raise XAFSConfigError("S4.5", f"厚度必须 >0 μm，收到 {thickness_um}")
    if not (density > 0):
        raise XAFSConfigError("S4.5", f"密度必须 >0 g/cm³，收到 {density}")
    for name, ang in (("anginp", anginp), ("angout", angout)):
        if not (0.0 < ang <= 90.0):
            raise XAFSConfigError(
                "S4.5", f"{name}={ang}° 必须在 (0, 90] 内")

    energy = np.asarray(norm.energy, dtype=float)
    e0 = float(norm.e0)
    e_edge_tab = float(xray_edge(elem, edge).energy)
    e_ref = e0 + 50.0
    # _booth_mu_terms 的输入是 k；锚点仍是表值边 +0.1 eV（与 Booth 一致）
    k_ref = np.sqrt(max(K2 * (e_ref - (e_edge_tab + 0.1)), 0.0))
    mu_t, mu_a, mu_f, e_fluor, comp = _booth_mu_terms(
        formula, elem, edge, line, e0, np.asarray([k_ref], dtype=float),
        density)
    sin_phi = np.sin(np.deg2rad(anginp))
    sin_theta = np.sin(np.deg2rad(angout))
    g = sin_phi / sin_theta
    mu_bg = float(mu_t[0] - mu_a[0])
    alpha = (mu_bg + g * float(mu_f)) / float(mu_a[0])
    if not np.isfinite(alpha) or alpha <= 0:
        raise XAFSDataError(
            "S4.5", f"Tröger α={alpha!r} 非正/非有限；检查化学式/荧光线/边")

    floor = 0.25 * (alpha + 1.0)
    applied = (energy >= e0 - 15.0)

    def _invert(y):
        out = np.array(y, dtype=float, copy=True)
        y = np.asarray(y, dtype=float)
        mask = applied & (y > 0)
        if not np.any(mask):
            return out, True, 0
        denom = alpha + 1.0 - y[mask]
        n_clamped = int(np.count_nonzero(denom < floor))
        denom = np.maximum(denom, floor)
        out[mask] = alpha * y[mask] / denom
        return out, n_clamped == 0, n_clamped

    mu_corr, ok_mu, clamp_mu = _invert(norm.mu_norm)
    if norm.flattened is not None:
        flat_corr, ok_flat, clamp_flat = _invert(norm.flattened)
    else:
        flat_corr, ok_flat, clamp_flat = None, True, 0

    finite = bool(np.all(np.isfinite(mu_corr))
                  and (flat_corr is None or np.all(np.isfinite(flat_corr))))
    reliable = bool(ok_mu and ok_flat and finite
                    and (clamp_mu + clamp_flat) == 0)
    info = {
        "method": "Tröger/FLUO inverse (XANES energy domain)",
        "corrected": True,
        "xanes_corrected": True,
        "exafs_modified": False,
        "formula": formula, "elem": elem, "edge": edge, "line": line,
        "thickness_um": float(thickness_um), "density": float(density),
        "anginp": float(anginp), "angout": float(angout), "g": float(g),
        "energy_ref_eV": float(e_ref),
        "fluor_energy_eV": float(e_fluor),
        "alpha": float(alpha),
        "denominator_floor": float(floor),
        "clamped_points": int(clamp_mu + clamp_flat),
        "applied_from_eV": float(e0 - 15.0),
        "reliable": reliable,
        "warning": (None if reliable else
                    "Tröger 分母触底被钳制：白线尖峰过强或 α 不可靠，"
                    "校正后的白线/LCF 权重不能当无系统误差结果"),
        "note": ("仅修 XANES；EXAFS χ(k) 必须走独立的 Booth 校正。"
                 "总化学式/厚度/有效密度/角度任一不可靠，α 与反解结果"
                 "即不可靠。"),
    }
    corrected = NormalizedScan(
        energy=np.array(energy, copy=True),
        mu=np.array(norm.mu, dtype=float, copy=True),
        mu_norm=mu_corr,
        pre_edge_line=np.array(norm.pre_edge_line, dtype=float, copy=True),
        edge_step=float(norm.edge_step),
        plateau=float(norm.plateau),
        e0=e0,
        flattened=(None if flat_corr is None
                   else np.asarray(flat_corr, dtype=float)))
    return corrected, info


# ================================================================
# 5. Rbkg 寻优
# ================================================================

def measure_shell_metrics(chi_k: ChiK, kmin=3.0, kmax=12.0, kweight=2):
    """
    用一次 FT 量出「第一壳层」与「低 R 残留」两个指标。

    返回:
        r1         第一峰位置 (Å)
        amplitude  第一峰高度
        leakage    R < 0.8 Å 的最大幅值
        ratio      leakage / amplitude
    """
    ft = ft_to_r_space(chi_k.k, chi_k.chi, None, kmin, kmax, kweight)
    mag = ft.magnitude

    # 第一壳表观峰窗口：上限 2.6 覆盖金属第一壳到重金属（如 Pt–Pt 2.77 →
    # 表观 ~2.45）；元素普适性要求，不再按 Fe 校准
    shell = (ft.r > 1.0) & (ft.r < 2.6)
    r1 = float(ft.r[shell][np.argmax(mag[shell])]) if shell.any() else np.nan
    amplitude = float(mag[shell].max()) if shell.any() else 0.0

    low = ft.r < 0.8
    leakage = float(mag[low].max()) if low.any() else 0.0
    return {"r1": r1, "amplitude": amplitude, "leakage": leakage,
            "ratio": leakage / amplitude if amplitude > 0 else np.inf}


def optimize_rbkg(norm: NormalizedScan, rbkg_range=(0.7, 1.6), n_grid=10,
                  kmin=3.0, kmax=12.0, kweight=2, leak_veto=0.20):
    """
    在 Rbkg 网格上寻优。

    判据（乘积型）:
        score = (1 − 低R残留比) × (主峰强度 / 全网格最大主峰强度)
        veto  : 低R残留比 > 0.20 直接淘汰

    为什么不用「主峰强度最大化」:
        实测 Rbkg 从 0.7 扫到 1.6，主峰强度为
        1.794 / 1.755 / 1.713 / 1.693 / 1.701 / 1.810 / 1.290 / 1.339 / 1.432 / 0.434
        **非单调**。按强度最大会选中 Rbkg=1.2 并否决最干净的 1.5。

    返回: (最优结果字典, 全表)
    """
    table = []
    for rbkg in np.linspace(rbkg_range[0], rbkg_range[1], n_grid):
        chi_k = extract_chi_k(norm, rbkg=rbkg, kweight=kweight)
        metrics = measure_shell_metrics(chi_k, kmin, kmax, kweight)
        table.append({"rbkg": float(rbkg), **metrics})

    amp_max = max(row["amplitude"] for row in table) or 1.0
    for row in table:
        row["score"] = (1.0 - min(row["ratio"], 1.0)) * (row["amplitude"] / amp_max)
        row["veto"] = bool(row["ratio"] > leak_veto)

    valid = [row for row in table if not row["veto"]]
    if valid:
        s_max = max(r["score"] for r in valid)
        # 并列裁决：与最优相差 <5% 时取最小 Rbkg。
        # 单调爬向网格边缘是「样条吃壳层→泄漏比降低→得分虚高」的病理形态
        # （NiO 实测：0.84→0.96 单调爬到 1.60，N(Ni-O) 被吃掉一半），
        # 得分不可分辨时低 Rbkg 对第一壳振幅更安全。
        tied = [r for r in valid if r["score"] >= 0.95 * s_max]
        best = min(tied, key=lambda r: r["rbkg"])
    else:
        best = min(table, key=lambda r: r["ratio"])

    # 置信度：区分「寻优得到的」与「在噪声里挑出来的」
    scores = [row["score"] for row in table]
    valid_fraction = len(valid) / len(table)
    roughness = float(np.median(np.abs(np.diff(scores))) / (max(scores) or 1.0)) if len(scores) > 1 else 0.0
    top2 = sorted((r["score"] for r in valid), reverse=True)[:2]
    margin = float((top2[0] - top2[1]) / top2[0]) if len(top2) > 1 and top2[0] > 0 else 0.0

    if valid_fraction >= 0.6 and roughness < 0.10 and margin > 0.05:
        confidence = "high"
    elif valid_fraction >= 0.4 and roughness < 0.20:
        confidence = "medium"
    else:
        confidence = "low"

    best.update({"confidence": confidence,
                 "confidence_metrics": {"valid_fraction": round(valid_fraction, 2),
                                        "roughness": round(roughness, 3),
                                        "margin": round(margin, 3)}})
    return best, table


# ================================================================
# 6. 傅里叶变换
# ================================================================

def back_ft_to_q_space(r, chi_r, rmin=1.0, rmax=2.3, dr=0.1) -> BackFourierTransform:
    """R 窗反傅里叶变换（Back-FT）：χ(R) → χ(q)。

    对标 Athena 的「反傅里叶滤波」与 Artemis 的 q 空间拟合诊断：只保留
    [rmin, rmax] 内的壳层，反变换回波数空间，得到该壳层的振荡相位与
    包络。q 轴由 larch.xftr 决定（0–30 Å⁻¹，步长 0.05）。

    数学:
        χ(q) = ∫[rmin,rmax] χ(R)·W_R(R)·e^{+2iqR} dR
    与正变换 e^{−2ikR} 互为逆变换；R 窗用 Hanning（dr 为过渡宽度）。
    """
    r = np.asarray(r, dtype=float)
    chi_r = np.asarray(chi_r, dtype=complex)
    if r.shape != chi_r.shape:
        raise XAFSConfigError(
            "S5", f"back-FT 输入长度不一致: r={r.shape}, chi_r={chi_r.shape}")
    mask = (r >= rmin) & (r <= rmax)
    if int(mask.sum()) < 5:
        raise XAFSConfigError(
            "S5", f"back-FT R 区间 [{rmin},{rmax}] 内仅 {int(mask.sum())} 点")
    from larch.xafs import xftr

    g = Group(r=r, chir=chi_r)
    xftr(g, rmin=float(rmin), rmax=float(rmax), dr=float(dr),
         window="hanning")
    q = np.asarray(g.q, dtype=float)
    chi_q = (np.asarray(g.chiq_re, dtype=float)
             + 1j * np.asarray(g.chiq_im, dtype=float))
    if q.shape != chi_q.shape:
        # larch 的 xftr 在粗 R 网格上会输出更长的 q 轴（到 30 Å⁻¹），
        # 但 chiq 只保留 nfft/2 个有效点（实测 1024）。此时 q 的多余
        # 尾部没有对应的 χ(q)，截断而不是插值/猜测。
        if q.size > chi_q.size and q.ndim == 1 and chi_q.ndim == 1:
            q = q[:chi_q.size]
        else:
            raise XAFSConfigError("S5", "back-FT 输出 q 与 χ(q) 长度不一致")
    return BackFourierTransform(
        q=q, chi_q=chi_q, rmin=float(rmin), rmax=float(rmax),
        dr=float(dr), r_window=np.asarray(g.rwin, dtype=float))


def ft_to_r_space(k, chi, r=None, kmin=3.0, kmax=12.0, kweight=2, dk=1.0,
                  rmax_out: float | None = None,
                  phase_correct=None) -> FourierTransform:
    # r 参数保留是为了向后兼容，但**实际被忽略**——R 网格由 xftf 内部决定。
    """
    傅里叶变换到 R 空间。

    数学:
        χ(R) = ∫[kmin,kmax] k^n · χ(k) · W(k) · e^{−2ikR} dk

    关键（易错）:
        相位因子必须是 e^{−2ikR}，与 larch.xftf 一致。
        若误写成 e^{−2πikR}，R 轴会被压缩 π 倍——
        Fe–S 的 1.82 Å 会跑到 0.58 Å，拟合窗口随之全部错位
        （实测会把 R-factor 从 0.02 顶到 0.87）。

    B11-1 相移校正:
        ``phase_correct`` 给 FEFF 路径总相移 φ_tot(k)（长度与 k 相同）。
        变换前乘 ``exp(−iφ_tot)``，峰位从“表观 R”（约 Reff−0.3~0.5 Å）
        回到真实 Reff。默认 None，历史数值逐位不变。
        实现用线性性：xftf(re+i·im)=xftf(re)+i·xftf(im)。

    返回 FourierTransform（含复数 χ(R)，实部虚部都能取）。
    """
    k = np.asarray(k)
    chi = np.asarray(chi)

    mask = (k >= kmin) & (k <= kmax)
    if mask.sum() < 8:
        raise XAFSConfigError("S5", f"k 区间 [{kmin},{kmax}] 内仅 {mask.sum()} 点")

    # 内部调用 larch.xftf 而不是自行搭 FFT。
    # 原因：FT 的归一化因子、窗函数端点处理、相位符号有多种约定，
    # 自研版本与 larch 实测差 0.053 Å（1.420 vs 1.473），且对 dk 不敏感；
    # 而 0.05 Å 会直接传导到拟合出的 R。既然 Athena/Artemis 用的就是这套，
    # 这里保持一致，避免"看起来对、实际错位"。
    from larch import Group
    from larch.xafs import xftf

    xftf_kw = dict(kmin=kmin, kmax=kmax, kweight=kweight, dk=dk,
                   window="hanning")
    # 默认不传 rmax_out，保持历史 R 网格逐位一致（现有 FT/拟合/报告都
    # 消费这条网格）。只有高 R 噪声估计等明确需要远区的调用才扩展。
    if rmax_out is not None:
        xftf_kw["rmax_out"] = float(rmax_out)

    def _one(y):
        g = Group(k=np.asarray(k), chi=np.asarray(y))
        xftf(g, **xftf_kw)
        return (np.asarray(g.chir_re) + 1j * np.asarray(g.chir_im),
                np.asarray(g.r))

    if phase_correct is None:
        chi_r, r_native = _one(chi)
    else:
        phi = np.asarray(phase_correct, dtype=float)
        if phi.ndim != 1 or phi.size != k.size:
            raise XAFSConfigError(
                "S5", f"phase_correct 需与 k 等长（{k.size}），"
                      f"收到 shape={phi.shape}")
        if not np.all(np.isfinite(phi)):
            raise XAFSConfigError("S5", "phase_correct 含非有限值")
        rotated = chi * np.exp(-1j * phi)      # Re + i·Im
        ft_re, r_native = _one(np.real(rotated))
        ft_im, r_im = _one(np.imag(rotated))
        if r_im.shape != r_native.shape:
            raise XAFSConfigError("S5", "相移校正 FT 的 R 网格不一致")
        # FT(re + i·im) = FT(re) + i·FT(im)
        chi_r = ft_re + 1j * ft_im

    return FourierTransform(r=r_native, chi_r=chi_r, kmin=float(kmin),
                            kmax=float(kmax), kweight=int(kweight), dk=float(dk))


# ================================================================
# 7. 自测
# ================================================================

if __name__ == "__main__":
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    from xafs_io import read_scan

    target = (Path(__file__).resolve().parent / "data" / "练手数据" /
              "02_Fe_standards" / "Fe_lepidocrocite.000")
    scan = read_scan(target)

    e0_info = determine_e0(scan)
    print(f"E0 (derivative) : {e0_info['e0_sample']:.2f} eV")
    print(f"E0 (half_height): {find_e0_half_height(scan.energy, scan.mu, e0_info['e0_sample']):.2f} eV")

    norm = normalize_absorption(scan, e0_info["e0_sample"])
    print(f"edge_step       : {norm.edge_step:.4f}")
    print(f"边后平台        : {norm.plateau:.4f}")

    best, table = optimize_rbkg(norm)
    print(f"Rbkg 选中       : {best['rbkg']:.2f} Å  置信度 {best['confidence']}")

    chi_k = extract_chi_k(norm, rbkg=best["rbkg"])
    ft = ft_to_r_space(chi_k.k, chi_k.chi, np.arange(0.0, 4.5, 0.02))
    shell = (ft.r > 1.0) & (ft.r < 2.6)
    print(f"FT 第一峰       : {ft.r[shell][np.argmax(ft.magnitude[shell])]:.3f} Å")
