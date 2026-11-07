# -*- coding: utf-8 -*-
"""
larch_wt.py — 阶段 6：连续小波变换
================================================================
对应「Larch + Python」这一档功能。**larch 本身没有小波变换函数**
（`larch.xafs` 里只有 feff6l / feffit / xftf，没有任何 wavelet），
所以这一层必须自己实现。

数学:
    W(R, k) = ∫ χ(k')·k'^n·g((k' − k)/σ)·e^{−2ik'R} dk'

    g(u) = (1/√σ)·exp(−u²/2)         零载频高斯包络

**为什么不用带载频的经典 Morlet** `e^{i k₀(k'−k)/σ}`：
    它的峰位满足 R_wt = k₀/(2σ) − R₀，是镜像关系。
    实测 R₀ = 1.5 / 2.0 / 2.5 / 3.0 → 输出 3.50 / 3.00 / 2.50 / 2.00，
    正好反过来。改用零载频包络后，合成信号偏差 0.00 Å，
    R 轴与 |χ(R)| 严格对齐。

**σ 不能太小**：σ < 1.2 时相邻壳层粘连。实测 σ=0.8 时
Fe–O 与 Fe–S 挤到 1.82 / 2.06 分不开；σ=1.5 才回到 1.44 / 2.66。
代价是 k 方向定位变差，有效 kc 区间被压到 [kmin+σ, kmax−σ]。

依赖方向：core ← io ← athena ← **larch_wt**（读 ChiK）。
本模块不 import artemis。
"""
from __future__ import annotations

import numpy as np

from xafs_core import ChiK, WaveletMap, XAFSConfigError


def morlet_envelope(u, sigma: float):
    """
    零载频高斯包络（Morlet 的实包络）。

    数学:
        g(u) = (1/√σ)·exp(−u²/2),   u = (k' − k)/σ
    """
    return np.exp(-0.5 * u ** 2) / np.sqrt(sigma)


def continuous_wavelet_transform(chi_k: ChiK, sigma: float = 1.5,
                                 kmin: float = 3.0, kmax: float = 14.0,
                                 kweight: int = 2, rmax: float = 4.0,
                                 dk: float = 0.1,
                                 normalize: bool = True,
                                 phase_correct=None) -> WaveletMap:
    """
    连续小波变换，得到 |W(R,k)| 二维图。

    参数:
        chi_k  : athena.extract_chi_k 的输出
        sigma  : 包络宽度（Å⁻¹）。≥1.5 才能分开相邻壳层
        kmin   : 小波分析的 k 下限
        kmax   : k 上限
        rmax   : R 轴上限（Å）
        normalize: True 时每张图各自归一到峰值 1（默认，便于看形状）；
                   False 时保留原始幅度，供多 kweight/σ **跨图比较**
                   （kweight 改变 k 加权，σ 改变包络积分，各自归一后
                   强度不可比）。
        phase_correct: B11-1 可选 FEFF 总相移 φ_tot(k)，长度与 chi_k.k
                   相同（或等长均匀 k 数组）。变换前乘 e^{−iφ}，峰位从
                   表观 R 移到真实 Reff；默认 None 时历史数值逐位不变。

    返回 WaveletMap，magnitude 形状为 (len(k_centers), len(r))。
    """
    # kmax 不得超出数据范围：np.interp 对越界查询用边界值外推，
    # 窗口跨出数据末端会造出平台状伪影（run_xafs 调用处有同样防护）。
    kmax = min(kmax, float(np.max(chi_k.k)))

    if sigma < 1.0:
        raise XAFSConfigError(
            "WT", f"σ={sigma} 过小，相邻壳层会粘连；建议 σ ≥ 1.2")

    # 1) 重采样到均匀 k 网格；两端加 Hanning 渐变，避免边界阶跃
    #    （不加 taper 时，窗口跨越边界只剩半边，会在 R≈0 造出强假峰）
    target = np.arange(0.0, kmax + 1.5, dk)
    chi_resampled = np.zeros_like(target)
    inside = (target >= kmin) & (target <= kmax)
    chi_resampled[inside] = np.interp(target[inside], chi_k.k, chi_k.chi)

    # B11-1：相位校正。φ 先在原始 k 网格上插值到均匀 target；
    # 窗外的 taper 为 0，不会引入相位跳变。
    if phase_correct is not None:
        phi = np.asarray(phase_correct, dtype=float)
        if phi.ndim != 1:
            raise XAFSConfigError("WT", "phase_correct 必须是一维数组")
        if phi.size == np.asarray(chi_k.k).size:
            phi_t = np.interp(target, np.asarray(chi_k.k), phi)
        elif phi.size == target.size:
            phi_t = phi
        else:
            raise XAFSConfigError(
                "WT", f"phase_correct 长度 {phi.size} 与 k 网格 "
                      f"{np.asarray(chi_k.k).size} 不匹配")
        if not np.all(np.isfinite(phi_t)):
            raise XAFSConfigError("WT", "phase_correct 含非有限值")
        chi_resampled = chi_resampled * np.exp(-1j * phi_t)

    ramp = max(int(1.0 / dk), 1)
    taper = np.ones_like(target)
    taper[:ramp] = 0.5 * (1 - np.cos(np.pi * np.arange(ramp) / ramp))
    taper[-ramp:] = taper[:ramp][::-1]
    weighted = chi_resampled * target ** kweight * taper

    # 2) 滑窗中心只在窗口完整的位置取，避免半边窗
    k_centers = np.arange(kmin + sigma, kmax - sigma + 1e-9, dk)
    r_grid = np.arange(0.0, rmax + 1e-9, 0.02)

    wave_map = np.zeros((len(k_centers), len(r_grid)))
    c_map = np.zeros((len(k_centers), len(r_grid)), dtype=complex)
    for i, kc in enumerate(k_centers):
        u = (target - kc) / sigma
        envelope = morlet_envelope(u, sigma)
        transform = np.exp(-2j * np.outer(r_grid, target)) @ (weighted * envelope)
        c_row = transform * dk
        c_map[i, :] = c_row
        wave_map[i, :] = np.abs(c_row)

    peak = wave_map.max()
    if normalize and peak > 0:
        # 幅值与复数矩阵必须除以**同一个** peak，否则 |complex_map|
        # 与 magnitude 不再逐点相等，复数残差会混入人为标度误差。
        wave_map = wave_map / peak
        c_map = c_map / peak

    return WaveletMap(k_centers=k_centers, r=r_grid,
                      magnitude=wave_map, sigma=float(sigma),
                      complex_map=c_map)


def locate_shell(wave_map: WaveletMap, r_lo: float, r_hi: float):
    """
    在给定 R 窗口内找最强响应点，返回 (k, R, 相对强度)。

    注意：窗口内若同时存在两条壳层分支（例如 Fe3O4 的 Fe–O@1.4 与
    Fe–Fe@2.7），argmax 会在两分支间来回跳；跨样品比较时应改用
    ``locate_shell_near`` 锚定到 FT 已选中的峰。
    """
    columns = (wave_map.r >= r_lo) & (wave_map.r <= r_hi)
    if not columns.any():
        raise XAFSConfigError("WT", f"R 窗口 [{r_lo},{r_hi}] 内无数据点")
    sub = wave_map.magnitude[:, columns]
    i, j = np.unravel_index(int(np.argmax(sub)), sub.shape)
    return (float(wave_map.k_centers[i]),
            float(wave_map.r[columns][j]),
            float(sub[i, j]))


def extract_k_profile(wave_map: WaveletMap, r_center: float,
                      half_width: float = 0.25) -> dict:
    """提取指定壳层 R_i 处的一维小波包络切片 |W(k; R_i)|。

    为什么需要一维切片：二维伪彩图里重叠的 M–O（轻元素，包络峰
    k≈4–5 Å⁻¹）与 M–M（3d 金属，k≈6.5–8.5 Å⁻¹；Pt 等重元素 k>9 Å⁻¹）
    只表现为颜色深浅，肉眼很难分辨。沿 R_i 抽一条 |W(k)| 曲线，峰位
    直接回答“这个壳层主要由哪个 k 的光电子贡献”，双峰则提示同一 R
    窗口里叠了两种散射体。

    返回字段:
        k / ridge_profile / slice_profile : 一维曲线（同一 k 网格）
        r_track    : 每个 k 处包络最大值所在的 R（嵴线轨迹）
        k_peak     : 嵴线主峰对应的 k（Å⁻¹）
        r_peak     : 主峰处的 R（Å）
        slice_r    : 单行切片使用的实际 R（最接近 r_center 的网格）
        bimodal / k_secondary / secondary_ratio / peak_separation :
                     双峰干涉判据（次峰 ≥20% 主峰且间距 ≥1.5 Å⁻¹）
    """
    if not np.isfinite(r_center):
        raise XAFSConfigError("WT", f"r_center={r_center} 非有限值")
    if not np.isfinite(half_width) or half_width <= 0:
        raise XAFSConfigError("WT", f"half_width={half_width} 必须 > 0")

    r = np.asarray(wave_map.r, dtype=float)
    mag = np.asarray(wave_map.magnitude, dtype=float)
    k = np.asarray(wave_map.k_centers, dtype=float)
    if mag.shape != (len(k), len(r)):
        raise XAFSConfigError(
            "WT", f"magnitude 形状 {mag.shape} 与 k×R 网格 "
                  f"({len(k)}, {len(r)}) 不一致")

    lo, hi = float(r_center) - float(half_width), float(r_center) + float(half_width)
    cols = np.flatnonzero((r >= lo) & (r <= hi))
    if cols.size == 0:
        raise XAFSConfigError(
            "WT", f"R 窗 [{lo:.2f}, {hi:.2f}] Å 内无 WT 网格点")

    sub = mag[:, cols]
    ridge = sub.max(axis=1)
    ridge_cols = np.argmax(sub, axis=1)
    r_track = r[cols][ridge_cols]
    idx_center = int(np.argmin(np.abs(r - float(r_center))))

    # ---- B5-1 相位：沿嵴线与中心单行 ----
    # phase_ridge 取每个 k 处包络最大值所在的 R 的复数相位；拍频节点
    # 处相位会跳变 π，仅看模值无法与"双元素叠加"区分。
    complex_map = getattr(wave_map, "complex_map", None)
    if complex_map is not None:
        c_map = np.asarray(complex_map)
        phase_ridge = np.angle(c_map[np.arange(len(k)), cols[ridge_cols]])
        phase_slice = np.angle(c_map[:, idx_center]).copy()
    else:
        phase_ridge = None
        phase_slice = None

    i_peak = int(np.argmax(ridge))
    k_peak = float(k[i_peak])
    peak_amp = float(ridge[i_peak])

    # ---- 双峰判据：轻元素 M–O 与 M–M 在同一 R 窗口重叠时，包络会出现
    # 两个相距 ≥1.5 Å⁻¹ 的峰；单壳层的高斯型包络只有一个主峰。----
    bimodal = False
    k_secondary = None
    secondary_ratio = 0.0
    separation = 0.0
    beat_node_phase_jump = 0.0
    if peak_amp > 0 and len(ridge) >= 7:
        from scipy.signal import find_peaks
        # 轻微平滑以抑制网格离散造成的锯齿，不改变峰位（dk=0.1）。
        try:
            from scipy.ndimage import gaussian_filter1d
            smooth = gaussian_filter1d(ridge, sigma=1.0)
        except Exception:      # pragma: no cover - scipy 总是可用
            smooth = ridge
        dk = float(np.median(np.diff(k))) if len(k) > 1 else 0.1
        peaks, props = find_peaks(
            smooth, prominence=max(0.03 * peak_amp, 1e-12),
            distance=max(1, int(round(1.0 / max(dk, 1e-6)))))
        if len(peaks) >= 2:
            order = sorted(peaks, key=lambda p: smooth[p], reverse=True)
            p_main, p_second = int(order[0]), int(order[1])
            h_main, h_second = float(smooth[p_main]), float(smooth[p_second])
            if h_second >= 0.20 * h_main:
                sep = abs(float(k[p_main]) - float(k[p_second]))
                if sep >= 1.5:
                    bimodal = True
                    k_secondary = float(k[p_second])
                    secondary_ratio = (h_second / h_main) if h_main > 0 else 0.0
                    separation = sep
                    # 两峰之间的最小幅值谷：若复数相位在此处跳变 π，
                    # 说明是同元素长短亚壳层的相消拍频，而不是双元素叠加。
                    if phase_ridge is not None:
                        lo_i, hi_i = sorted((p_main, p_second))
                        valley = lo_i + int(np.argmin(ridge[lo_i:hi_i + 1]))
                        a = max(0, valley - 2)
                        b = min(len(k), valley + 3)
                        if b - a >= 2:
                            ph = np.unwrap(phase_ridge[a:b])
                            beat_node_phase_jump = float(
                                np.max(np.abs(np.diff(ph))))

    return {
        "k": k,
        "r_center": float(r_center),
        "r_lo": lo,
        "r_hi": hi,
        "ridge_profile": ridge,
        "slice_profile": mag[:, idx_center].copy(),
        "phase_ridge": phase_ridge,
        "phase_slice": phase_slice,
        "r_track": r_track,
        "slice_r": float(r[idx_center]),
        "k_peak": k_peak,
        "k_peak_slice": float(k[int(np.argmax(mag[:, idx_center]))]),
        "r_peak": float(r_track[i_peak]),
        "peak_amplitude": peak_amp,
        "bimodal": bool(bimodal),
        "k_secondary": k_secondary,
        "secondary_ratio": float(secondary_ratio),
        "peak_separation": float(separation),
        "beat_node_phase_jump": float(beat_node_phase_jump),
    }


def compute_model_wt_residual(chi_k_data: ChiK, model_k: np.ndarray,
                              sigma: float = 1.5, kmin: float = 3.0,
                              kmax: float = 12.0, kweight: int = 2,
                              rmax: float = 4.0) -> dict:
    """实验与模型在**完全相同** WT 参数下的二维幅值残差。

    契约（对标 Larch-WT / Hama 的实验-模型 WT 对比）：
      - 两侧都用 ``(sigma, kmin, kmax, kweight, rmax, normalize=False)``；
      - 统一除以 ``wt_data.magnitude.max()`` 后再相减，避免模型振幅
        整体缩放把残差图染成单一颜色；
      - 二维 WT R-factor = Σ(|W_data|−|W_model|)² / Σ|W_data|²。

    返回:
        wt_data / wt_model : WaveletMap（未归一化，共享同一物理标度）
        residual_map       : 归一化幅值差，形状 = wt_data.magnitude.shape
        r_factor           : 二维 WT R-factor
        scale              : 归一化因子（wt_data 的峰值）
    """
    model_k = np.asarray(model_k, dtype=float)
    k_data = np.asarray(chi_k_data.k, dtype=float)
    if model_k.shape != k_data.shape:
        if model_k.size == k_data.size:
            model_k = model_k.reshape(k_data.shape)
        else:
            raise XAFSConfigError(
                "WT", f"model_k 长度 {model_k.size} 与数据 k 网格 "
                      f"{k_data.size} 不一致")
    if not np.all(np.isfinite(model_k)):
        raise XAFSConfigError("WT", "model_k 含非有限值")

    model_ck = ChiK(
        k=k_data,
        chi=model_k,
        background=np.asarray(chi_k_data.background, dtype=float),
        mu_norm=np.asarray(chi_k_data.mu_norm, dtype=float),
        rbkg=float(getattr(chi_k_data, "rbkg", 0.0)),
    )

    wt_data = continuous_wavelet_transform(
        chi_k_data, sigma=float(sigma), kmin=float(kmin),
        kmax=float(kmax), kweight=int(kweight), rmax=float(rmax),
        normalize=False)
    wt_model = continuous_wavelet_transform(
        model_ck, sigma=float(sigma), kmin=float(kmin),
        kmax=float(kmax), kweight=int(kweight), rmax=float(rmax),
        normalize=False)

    scale = float(wt_data.magnitude.max())
    if not np.isfinite(scale) or scale <= 0:
        raise XAFSConfigError("WT", "实验 WT 峰值为 0，无法归一化残差")
    data_n = wt_data.magnitude / scale
    model_n = wt_model.magnitude / scale
    residual = data_n - model_n
    denom = float(np.sum(data_n ** 2))
    r_2d = (float(np.sum(residual ** 2) / denom) if denom > 0 else float("inf"))

    # ---- B5-1 复数域残差：幅值相等但相位相反的模型不再被判为 0 ----
    data_c = getattr(wt_data, "complex_map", None)
    model_c = getattr(wt_model, "complex_map", None)
    if data_c is not None and model_c is not None:
        diff_c = np.asarray(data_c) - np.asarray(model_c)
        complex_residual_map = np.abs(diff_c) / scale
        denom_c = float(np.sum(np.abs(data_c) ** 2))
        complex_r = (float(np.sum(np.abs(diff_c) ** 2) / denom_c)
                     if denom_c > 0 else float("inf"))
    else:                              # 旧版对象/自定义 WaveletMap 的安全降级
        complex_residual_map = None
        complex_r = None

    return {
        "wt_data": wt_data,
        "wt_model": wt_model,
        "residual_map": residual,
        "complex_residual_map": complex_residual_map,
        "normalized_data": data_n,
        "normalized_model": model_n,
        "scale": scale,
        "r_factor": r_2d,
        "complex_r_factor": complex_r,
        "sigma": float(sigma),
        "kmin": float(kmin),
        "kmax": float(kmax),
        "kweight": int(kweight),
        "rmax": float(rmax),
    }


def locate_shell_near(wave_map: WaveletMap, r_center: float,
                      half_width: float = 0.5, r_floor: float = 0.8):
    """
    在 FT 已选峰位 r_center ± half_width 内定位 WT 响应。

    为什么需要这个函数（实测教训）：
        宽窗 [1.0, 2.6] 会同时装进 Fe3O4 的 Fe–O（WT R≈1.42，k≈4.5）
        与 Fe–Fe（WT R≈2.72，k≈5.8）两条分支，二者强度只差约 3%。
        荧光自吸收校正把低 k 振幅抬高 7% 后，argmax 就从 2.60 跳到
        1.42——但 FT 第一峰完全没动（仍是 2.577），说明这是选峰窗口
        的伪跳变，不是校正改变了壳层。
        以 FT 峰为中心取窄窗，WT 只回答"这个 FT 峰的能量来自哪个 k"，
        R 方向不再与 FT 竞争选峰。

    参数:
        r_center   FT 第一峰位（Å）
        half_width 半窗宽（Å）；默认 0.5，覆盖 σ=1.5 的 R 分辨率
        r_floor    R 下限，避免贴 0 的残余背景伪峰
    """
    if not np.isfinite(r_center):
        raise XAFSConfigError("WT", f"r_center={r_center} 非有限值")
    if half_width <= 0:
        raise XAFSConfigError("WT", f"half_width={half_width} 必须 > 0")
    r_lo = max(float(r_floor), float(r_center) - float(half_width))
    r_hi = float(r_center) + float(half_width)
    return locate_shell(wave_map, r_lo, r_hi)


# ================================================================
# 自测
# ================================================================

if __name__ == "__main__":
    import sys
    from pathlib import Path
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    import athena
    from xafs_io import read_scan

    target = (Path(__file__).resolve().parent / "data" / "练手数据" /
              "02_Fe_standards" / "Fe_lepidocrocite.000")
    scan = read_scan(target)
    e0 = athena.determine_e0(scan)["e0_sample"]
    norm = athena.normalize_absorption(scan, e0)
    chi_k = athena.extract_chi_k(norm, rbkg=1.20)

    wave_map = continuous_wavelet_transform(chi_k, sigma=1.5, kmax=14.0)
    k_o, r_o, amp_o = locate_shell(wave_map, 1.0, 2.0)
    k_f, r_f, amp_f = locate_shell(wave_map, 2.2, 3.3)
    print(f"WT 网格      : {wave_map.magnitude.shape[0]} k × {wave_map.magnitude.shape[1]} R")
    print(f"第一壳最强   : k={k_o:.1f} Å⁻¹, R={r_o:.2f} Å, 相对强度={amp_o:.2f}")
    print(f"第二壳最强   : k={k_f:.1f} Å⁻¹, R={r_f:.2f} Å, 相对强度={amp_f:.2f}")
