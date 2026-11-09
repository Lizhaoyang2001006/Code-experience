# -*- coding: utf-8 -*-
"""
xafs_pipeline.py — S1–S6 共享编排（run_xafs 与 xafs_report 的唯一上游）
================================================================
把「读数 → 合并 → E0 → 归一化 → Rbkg → χ(k) → FT → 小波」串成一次调用。

    输入: 原始数据文件路径列表（可多次扫描）
    输出: PipelineResult（xafs_core），含各阶段 dataclass 与关键标量
    关键决策:
      - 对齐失败（参考通道非独立/缺失）→ 降级为直接合并，不中断流程
      - 边后平台 ∉ [0.85, 1.15]      → 抛 XAFSError（归一化不可信）
      - rbkg 显式给出                 → 跳过网格寻优，仅计算壳层指标

**不含任何物理计算**——所有数学都在 athena / larch_wt 里；
**不写任何文件**——落盘是入口层（run_xafs / xafs_report）的事。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import athena
import larch_wt
import xafs_io
from absorber import default_fluorescence_line, get_absorber
from xafs_core import K2, PipelineResult, XAFSDataError, XAFSError


def configure_console():
    """把 stdout/stderr 切到 UTF-8，避免 GBK 控制台遇到 Å 崩溃。

    处理链会打印 Å、Å⁻¹、σ² 等符号；Python 在中文 Windows 默认用
    GBK 编码 stdout，``Å`` (U+00C5) 不在 GBK 中。入口层和
    process_sample 都调用一次，保证脚本/管道/报告共用同一行为。
    """
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8", errors="replace")
            except (AttributeError, ValueError, OSError):
                pass


def process_sample(paths, align: bool = False, absorber: str | None = None,
                   edge_energy: float | None = None,
                   e0_method: str = "derivative", manual_e0: float | None = None,
                   rbkg: float | None = None, kmin: float = 3.0,
                   kmax: float = 12.0, kweight: int = 2,
                   post_window=(150, 600), wt_sigma: float = 1.5,
                   strict: bool = True,
                   despike: bool = False,
                   deglitch: list[tuple[float, float]] | None = None,
                   deadtime: bool = False,
                   reject_outlier_scans: bool = False,
                   mback: bool = False,
                   chi_std=None, chi_std_k=None,
                   clamp_lo: float | None = None,
                   clamp_hi: float | None = None,
                   self_absorption_formula: str | None = None,
                   self_absorption_thickness: float | None = None,
                   self_absorption_density: float | None = None,
                   self_absorption_angles: tuple = (45.0, 45.0),
                   self_absorption_line: str | None = None,
                   verbose: bool = True) -> PipelineResult:
    """跑完 S1–S6（算法元素无关）。absorber 只做三件配置事：
    对齐查参考边能量、S2 定边搜索窗口、记入 res.absorber_symbol 传给拟合链。"""
    configure_console()
    paths = [str(p) for p in paths]
    abs_cfg = get_absorber(absorber) if absorber else None

    # ---- S1 解析 ----
    scans = xafs_io.read_scans(paths)
    if verbose:
        for scan, path in zip(scans, paths):
            print(f"  [S1] {Path(path).name:<28} "
                  f"{scan.meta['format_detected']:<6} "
                  f"mode={scan.meta['mode']:<6} n={len(scan):<5} "
                  f"E={scan.energy.min():.0f}-{scan.energy.max():.0f}")

    # ---- S1.2 多通道 MCA 坏道剔除（自动；仅在文件解析出 mca 时触发）----
    for scan in scans:
        mca = (scan.meta or {}).get("mca")
        if not mca or mca.get("channels") is None:
            continue
        summed, mca_info = athena.merge_mca_channels(
            mca["channels"], dtime=mca.get("dtime"))
        scan.ifluor = summed
        if scan.i0 is not None:
            scan.mu = summed / scan.i0
        else:
            scan.mu = summed
        # 通道级死时间已在 merge 内逐道校正，清空 dtime 防止二次校正。
        scan.dtime = None
        scan.meta["mca_merge"] = mca_info
        scan.meta["deadtime_corrected"] = bool(mca_info["deadtime_corrected"])
        if verbose:
            print(f"  [S1.2] MCA: 保留 {mca_info['n_kept']}/"
                  f"{mca_info['n_channels']} 道，剔除 {mca_info['rejected']}")

    # ---- S1.2 死时间校正（可选；逐点 dtime，默认关闭）----
    if deadtime:
        corrected = []
        for scan in scans:
            before = scan.meta.get("deadtime_corrected")
            scan = athena.apply_deadtime_correction(scan)
            corrected.append(scan)
            if verbose:
                print(f"  [S1.2] 死时间校正: {scan.meta.get('deadtime_note', '')}"
                      + ("" if scan.meta.get("deadtime_corrected") else "（无有效因子，未改动）"))
        scans = corrected
        if verbose and not any(s.meta.get("deadtime_corrected") for s in scans):
            print("       ⚠ 请求了 --deadtime，但所有扫描都没有有效 dtime 因子")

    # ---- S1.2 显式能量区间 deglitch（删点，不插值）----
    if deglitch:
        cleaned = []
        for scan in scans:
            scan, info = athena.deglitch_ranges(scan, deglitch)
            cleaned.append(scan)
            if verbose:
                print(f"  [S1.2] deglitch: 删除 {info['n_removed']}/{info['n_total']}"
                      f" 点（{info['fraction']:.1%}），区间 {info['ranges']}")
        scans = cleaned

    # ---- S1.2 去毛刺（可选；中值滑窗 + MAD 门禁，默认关闭）----
    if despike:
        total_spikes = 0
        despiked_scans = []
        for scan in scans:
            scan_d, d_info = athena.despike_scan(scan, kernel=5, threshold=3.0)
            total_spikes += int(d_info["n_spikes"])
            despiked_scans.append(scan_d)
        scans = despiked_scans
        if verbose:
            print(f"  [S1.2] 去毛刺: 替换 {total_spikes} 个离群点"
                  + ("（窗口=5, 阈值=3×MAD）" if total_spikes else ""))

    # ---- S1.5 合并（可选先对齐）----
    if align and len(scans) > 1:
        if edge_energy is None:     # 参考边能量按吸收元素查表（无元素硬编码）
            if abs_cfg is None:
                raise XAFSError("S1.5", "对齐需要参考边能量："
                                        "给 absorber（自动查表）或 edge_energy")
            edge_energy = abs_cfg.e0_edge
        try:
            scans, align_info = athena.align_by_reference(scans, edge_energy)
        except XAFSDataError as exc:
            # 参考通道非独立（如 Ir 测的是同一样品）或缺失时降级：
            # 跳过对齐直接合并，不让整条流水线崩掉。
            if verbose:
                print(f"  [S1.5] 参考通道不可用，已跳过对齐直接合并（{exc.message}）")
        else:
            if verbose:
                print(f"  [S1.5] 能量对齐: ΔE={align_info['shifts_eV']}  "
                      f"{align_info['note']}")
    merged = athena.merge_scans(
        scans, reject_outliers=bool(reject_outlier_scans))
    if verbose and reject_outlier_scans:
        info = (merged.meta or {}).get("merge_info") or {}
        print(f"  [S1.5] 多扫描相异度筛查: 保留 {info.get('n_kept')}/"
              f"{info.get('n_input')} 条，剔除下标 "
              f"{info.get('rejected_indices')}")

    # ---- S2 E0 ----
    e0_info = athena.determine_e0(merged, method=e0_method, manual_e0=manual_e0,
                                  guess=abs_cfg.e0_edge if abs_cfg else None)
    if verbose:
        print(f"  [S2] E0 = {e0_info['e0_sample']:.2f} eV ({e0_info['e0_method']})"
              + (f"   参考 {(e0_info.get('e0_reference') or float('nan')):.2f} eV"
                 if e0_info.get("e0_reference") else ""))
        if e0_info.get("reference_warning"):
            print(f"       ⚠ {e0_info['reference_warning']}")

    # ---- S3 归一化 ----
    edge_kind = (abs_cfg.edge if abs_cfg is not None else "K")
    if mback:
        if abs_cfg is None:
            raise XAFSError("S3", "MBACK 归一化需要 --absorber（理论截面按元素查表）")
        norm = athena.normalize_mback(merged, e0_info["e0_sample"],
                                      elem=abs_cfg.symbol, edge=edge_kind,
                                      post_window=post_window)
    else:
        norm = athena.normalize_absorption(merged, e0_info["e0_sample"],
                                           post_window=post_window)
    xanes_features = athena.extract_xanes_features(
        norm, edge=edge_kind, absorber=(abs_cfg.symbol if abs_cfg else None))
    if verbose:
        print(f"  [S3] edge_step = {norm.edge_step:.4f}   "
              f"边后平台 = {norm.plateau:.4f}")
        if xanes_features.get("neck_shoulder_present"):
            print(f"       XANES 颈部肩峰 {xanes_features['neck_shoulder_energy']:.1f} eV"
                  f"（导数双峰间距 {xanes_features['derivative_split_eV']:.1f} eV"
                  + ("，E0 跳变风险" if xanes_features["e0_jump_risk"] else "") + "）")
    if not (0.85 <= norm.plateau <= 1.15):
        msg = (f"边后平台 = {norm.plateau:.3f}（阈值 0.85–1.15）。"
               f"稀释样品的低 Z 基体呈 E^-3 衰减，直线外推扣不掉，"
               f"此时 FT 峰位与幅度均不可引用")
        # strict=True（run_xafs 默认）：未知样品宁可中断也不产出误导结果；
        # strict=False（xafs_report）：已知标样（如厚 Fe 箔平台 0.79）降级为警告。
        if strict:
            raise XAFSError(
                "S3", f"归一化失败：{msg}"
                      "（如为已知厚样/稀释标样，可加 --no-strict 降级为警告）")
        if verbose:
            print(f"       ⚠ {msg}（strict=False，仅警告）")

    # ---- S4 Rbkg ----
    if rbkg is None:
        best, table = athena.optimize_rbkg(norm, kmin=kmin, kmax=kmax,
                                           kweight=kweight)
    else:
        probe = athena.extract_chi_k(norm, rbkg=rbkg, kweight=kweight,
                                     chi_std=chi_std, k_std=chi_std_k,
                                     clamp_lo=clamp_lo, clamp_hi=clamp_hi)
        best = {"rbkg": float(rbkg),
                **athena.measure_shell_metrics(probe, kmin, kmax, kweight),
                "confidence": "manual"}
        table = []
    if verbose:
        print(f"  [S4] Rbkg = {best['rbkg']:.2f} Å   置信度 {best['confidence']}"
              + (f"  (有效候选 {best['confidence_metrics']['valid_fraction']:.0%},"
                 f" 崎岖度 {best['confidence_metrics']['roughness']:.3f})"
                 if "confidence_metrics" in best else ""))

    # ---- S5 背景 + FT ----
    chi_k = athena.extract_chi_k(norm, rbkg=best["rbkg"], kweight=kweight,
                                 chi_std=chi_std, k_std=chi_std_k,
                                 clamp_lo=clamp_lo, clamp_hi=clamp_hi)
    # σ_μ(E) 与 χ(k) 不同网格：用 k(E)=sqrt((E−E0)/K2) 把噪声插值到 k 轴，
    # 再交给 high_k_snr 做 σ_χ=σ_μ/edge_step。绝不对 σ_μ 直接按索引使用。
    sigma_mu_k = None
    if merged.sigma_mu is not None:
        e_at_k = float(norm.e0) + np.asarray(chi_k.k, dtype=float) ** 2 / K2
        sigma_mu_k = np.interp(e_at_k, np.asarray(norm.energy, dtype=float),
                               np.asarray(merged.sigma_mu, dtype=float))
    snr = athena.high_k_snr(chi_k, sigma_mu=sigma_mu_k,
                            edge_step=norm.edge_step, kweight=kweight)
    # ---- S4.5 荧光自吸收校正（可选；只改 χ(k) 振幅，不改 μ(E)）----
    sa_info = None
    if self_absorption_formula:
        mode = (merged.meta or {}).get("mode")
        if mode != "fluor":
            raise XAFSError(
                "S4.5", f"自吸收校正只适用于荧光模式，当前 mode={mode}；"
                        f"透射谱/文件自带 xmu 不能反推荧光自吸收")
        if abs_cfg is None:
            raise XAFSError(
                "S4.5", "自吸收校正需要吸收元素（--absorber）来确定边/荧光线")
        if (self_absorption_thickness is None
                or self_absorption_density is None):
            raise XAFSError(
                "S4.5", "Booth 校正需要样品厚度（μm）与有效密度（g/cm³）；"
                        "缺一不可，不能假设")
        line = self_absorption_line or default_fluorescence_line(abs_cfg.edge)
        chi_k, sa_info = athena.booth_self_absorption(
            chi_k, formula=self_absorption_formula, elem=abs_cfg.symbol,
            edge=abs_cfg.edge, line=line, e0=e0_info["e0_sample"],
            thickness_um=float(self_absorption_thickness),
            density=float(self_absorption_density),
            anginp=float(self_absorption_angles[0]),
            angout=float(self_absorption_angles[1]))
        # B6-2：XANES 区域独立反解（Tröger/FLUO）。必须在 chi_k 已由原始
        # norm 生成之后执行，确保 EXAFS 路径 100% 不受能量域校正影响。
        norm_xanes, xa_info = athena.fluo_xanes_self_absorption(
            norm, formula=self_absorption_formula, elem=abs_cfg.symbol,
            edge=abs_cfg.edge, line=line,
            thickness_um=float(self_absorption_thickness),
            density=float(self_absorption_density),
            anginp=float(self_absorption_angles[0]),
            angout=float(self_absorption_angles[1]))
        sa_info["xanes_self_absorption"] = xa_info
        sa_info["xanes_corrected"] = True
        norm = norm_xanes
        xanes_features = athena.extract_xanes_features(
            norm, edge=edge_kind,
            absorber=(abs_cfg.symbol if abs_cfg else None))
        if verbose:
            print(f"  [S4.5] Booth 自吸收校正: {sa_info['formula']}  "
                  f"d={sa_info['thickness_um']:g} μm  "
                  f"ρ={sa_info['density']:g} g/cm³  "
                  f"φ/θ={sa_info['anginp']:g}/{sa_info['angout']:g}°  "
                  f"校正因子中位数={sa_info['correction_factor_median']:.3f}")
            print(f"  [S4.5] XANES Tröger 反解: α={xa_info['alpha']:.3f}  "
                  f"reliable={xa_info['reliable']}  "
                  f"白线高度 {xanes_features.get('white_line_height'):.3f}")
            if not sa_info["reliable"]:
                print("       ⚠ 校正不可靠："
                      + sa_info.get("note", "检查化学式/厚度/密度/角度"))
            if sa_info.get("warning"):
                print("       ⚠ " + sa_info["warning"])
            if not xa_info["reliable"]:
                print("       ⚠ " + str(xa_info.get("warning")))
    ft = athena.ft_to_r_space(chi_k.k, chi_k.chi, None, kmin, kmax, kweight)
    # 第一壳表观峰窗口（1.0, 2.6）：元素普适——覆盖轻配体到重金属
    # （Pt–Pt 2.77 → 表观 ~2.45）；此前 2.2 是 Fe 校准的残留
    shell = (ft.r > 1.0) & (ft.r < 2.6)
    r1 = float(ft.r[shell][np.argmax(ft.magnitude[shell])]) if shell.any() \
        else np.nan
    if verbose:
        print(f"  [S5] FT 第一峰 = {r1:.3f} Å   k_max = {chi_k.k.max():.1f} Å⁻¹")

    # ---- S6 小波 ----
    wave_map = larch_wt.continuous_wavelet_transform(
        chi_k, sigma=wt_sigma, kmin=kmin, kmax=min(14.0, chi_k.k.max()),
        kweight=kweight)
    # WT 的 R 窗必须以 FT 已选峰为中心：宽窗 [1.0,2.6] 在多壳体系里会
    # 同时装进两条分支，argmax 会在它们之间跳（Fe3O4 实测：S4.5 把低 k
    # 抬高 7% 就把 WT 峰从 2.60 翻到 1.42，而 FT 峰纹丝不动）。
    if np.isfinite(r1):
        k_wt, r_wt, _ = larch_wt.locate_shell_near(wave_map, r1)
    else:  # 无 FT 峰可比时退回宽窗，仍给出诊断值
        k_wt, r_wt, _ = larch_wt.locate_shell(wave_map, 1.0, 2.6)
    if verbose:
        print(f"  [S6] WT 第一峰 = {r_wt:.3f} Å (k={k_wt:.1f})   "
              f"与 FT 差 {abs(r_wt - r1):.3f} Å")

    # 本模块止步于 S6：结构识别（S6.5）与拟合（S7）属拟合链，
    # 由 xafs_fit 统一编排——归约链不得反向调用拟合侧模块。
    return PipelineResult(scans=scans, merged=merged, e0=e0_info, norm=norm,
                          rbkg=best, chi_k=chi_k, ft=ft, wt=wave_map,
                          r1=r1, r_wt=r_wt, k_wt=k_wt, rbkg_table=table,
                          self_absorption=sa_info,
                          xanes_features=xanes_features,
                          snr=snr,
                          absorber_symbol=absorber)
