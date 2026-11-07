# -*- coding: utf-8 -*-
"""
xafs_diagnostics.py — 把 S1–S6 的「黑盒内部」画出来
================================================================
每个处理阶段一张诊断图，回答一个问题：

    S1 原始扫描   各扫描一致吗？通道角色对吗？        s1_raw_scans.png
    S2 E0 判定    导数峰 / 参考边 选得合理吗？        s2_e0.png
    S3 归一化     边前/边后基线外推靠谱吗？           s3_normalization.png
    S4 Rbkg 寻优  网格得分曲线是「寻优」还是「挑噪声」？ s4_rbkg.png
    S5 χ(k)→FT    窗口切哪了？第一峰在哪？            s5_chik_chir.png
    S6 小波       k–R 二维响应长什么样？              s6_wavelet.png

    输入: PipelineResult（xafs_pipeline.process_sample 的产物）
    输出: outdir 下 6 张 PNG（Agg 后端，不弹窗）
    关键决策: 只读 PipelineResult，不重算任何物理量——
              诊断图必须与产线结果严格同源。
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import athena
import artemis
import larch_wt
from xafs_core import PipelineResult, XAFSConfigError

SCAN_COLORS = ("#1b6ca8", "#d1495b", "#2a9d8f", "#e9c46a")
WT_SWEEP_SIGMAS = (1.0, 1.5, 2.0)
WT_SWEEP_KWEIGHTS = (1, 2, 3)


def _style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "axes.linewidth": 1.0, "font.size": 10,
    })


def _save(fig, outdir: Path, name: str, made: list):
    fig.savefig(outdir / name, dpi=200, bbox_inches="tight")
    plt.close(fig)
    made.append(outdir / name)


def _s1(res: PipelineResult, outdir, made):
    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(10, 3.6))
    e0 = res.e0["e0_sample"]
    for i, (scan, c) in enumerate(zip(res.scans, SCAN_COLORS)):
        ax_full.plot(scan.energy, scan.mu, color=c, lw=0.7,
                     label=f"scan {i + 1}")
        ax_zoom.plot(scan.energy, scan.mu, color=c, lw=0.9)
    ax_full.plot(res.merged.energy, res.merged.mu, "k-", lw=1.4,
                 label="merged")
    if res.merged.sigma_mu is not None:
        sig = np.asarray(res.merged.sigma_mu, dtype=float)
        if sig.shape == np.asarray(res.merged.energy).shape:
            ax_full.fill_between(res.merged.energy,
                                 res.merged.mu - sig, res.merged.mu + sig,
                                 color="0.5", alpha=0.35, lw=0,
                                 label=r"merged $\pm\sigma_\mu$")
            ax_zoom.fill_between(res.merged.energy,
                                 res.merged.mu - sig, res.merged.mu + sig,
                                 color="0.5", alpha=0.35, lw=0)
    ax_full.axvline(e0, ls="--", color="0.4", lw=0.8)
    ax_full.set_xlabel("Energy (eV)")
    ax_full.set_ylabel(r"$\mu(E)$")
    ax_full.legend(frameon=False, fontsize=8)
    ax_full.set_title("S1  full range", fontsize=10)
    ax_zoom.set_xlim(e0 - 120, e0 + 180)
    ax_zoom.axvline(e0, ls="--", color="0.4", lw=0.8)
    ax_zoom.set_xlabel("Energy (eV)")
    ax_zoom.set_title("S1  edge zoom (scan overlap)", fontsize=10)
    _save(fig, outdir, "s1_raw_scans.png", made)


def _s2(res: PipelineResult, outdir, made):
    e, mu = res.merged.energy, res.merged.mu
    e0 = res.e0["e0_sample"]
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.plot(e, mu, "b-", lw=1.0)
    ax.axvline(e0, ls="--", color="r", lw=1.0)
    ax.text(e0 + 8, mu.min() + 0.05 * (mu.max() - mu.min()),
            f"$E_0$ = {e0:.2f} eV\n({res.e0['e0_method']})", color="r",
            fontsize=9)
    if res.e0.get("e0_reference"):
        ax.axvline(res.e0["e0_reference"], ls=":", color="g", lw=1.0,
                   label=f"reference {res.e0['e0_reference']:.2f} eV")
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlim(e0 - 150, e0 + 250)
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel(r"$\mu(E)$")
    ax.set_title("S2  E0 determination", fontsize=10)
    # 内嵌一阶导数：E0 应落在导数主峰上
    ax_in = ax.inset_axes([0.52, 0.52, 0.44, 0.42])
    dmu = np.gradient(mu, e)
    ax_in.plot(e, dmu / dmu.max(), "k-", lw=0.8)
    ax_in.axvline(e0, ls="--", color="r", lw=0.8)
    ax_in.set_xlim(e0 - 60, e0 + 80)
    ax_in.set_title("d$\\mu$/dE (norm.)", fontsize=8)
    ax_in.tick_params(labelsize=7)
    _save(fig, outdir, "s2_e0.png", made)


def _s3(res: PipelineResult, outdir, made):
    n = res.norm
    fig, (ax_raw, ax_n) = plt.subplots(1, 2, figsize=(10, 3.6))
    ax_raw.plot(n.energy, n.mu, "b-", lw=0.9, label=r"$\mu(E)$")
    ax_raw.plot(n.energy, n.pre_edge_line, "r--", lw=1.0,
                label="pre-edge line (extrap.)")
    ax_raw.axvline(n.e0, ls="--", color="0.4", lw=0.8)
    ax_raw.annotate(f"edge step = {n.edge_step:.4f}",
                    xy=(n.e0, n.mu.max()), xytext=(n.e0 + 120, n.mu.max()),
                    fontsize=9, arrowprops=dict(arrowstyle="-", lw=0.7))
    ax_raw.set_xlabel("Energy (eV)")
    ax_raw.set_ylabel(r"$\mu(E)$")
    ax_raw.legend(frameon=False, fontsize=8)
    ax_raw.set_title("S3  pre-edge / edge step", fontsize=10)
    ax_n.plot(n.energy, n.mu_norm, "b-", lw=0.9)
    ax_n.axhline(1.0, ls=":", color="0.4", lw=0.8)
    ax_n.axvline(n.e0, ls="--", color="0.4", lw=0.8)
    ax_n.text(0.03, 0.05, f"post-edge plateau = {n.plateau:.4f}"
                          "  (accept 0.85–1.15)",
              transform=ax_n.transAxes, fontsize=9)
    ax_n.set_xlabel("Energy (eV)")
    ax_n.set_ylabel(r"Normalized $\mu(E)$")
    ax_n.set_title("S3  normalized", fontsize=10)
    _save(fig, outdir, "s3_normalization.png", made)


def _s4(res: PipelineResult, outdir, made):
    ck = res.chi_k
    # mu_norm/background 在能量网格上（与 k 不等长），先换算到 k 轴再画
    k_full = athena.wave_number(res.norm.energy, res.norm.e0)
    fig, (ax_bkg, ax_grid) = plt.subplots(1, 2, figsize=(10, 3.6))
    ax_bkg.plot(k_full, ck.mu_norm, "b-", lw=0.9, label=r"$\mu(E(k))$ norm")
    ax_bkg.plot(k_full, ck.background, "r-", lw=1.1,
                label=rf"spline bkg, $R_{{bkg}}$={ck.rbkg:.2f} Å")
    ax_bkg.set_xlabel(r"$k$ ($\AA^{-1}$)")
    ax_bkg.legend(frameon=False, fontsize=8)
    ax_bkg.set_title("S4  background spline", fontsize=10)
    table = res.rbkg_table
    if table:
        xs = [row["rbkg"] for row in table]
        ys = [row["score"] for row in table]
        ax_grid.plot(xs, ys, "ko-", ms=4, lw=0.8)
        ax_grid.axvline(res.rbkg["rbkg"], ls="--", color="r", lw=1.0)
        ax_grid.set_xlabel(r"$R_{bkg}$ (Å)")
        ax_grid.set_ylabel("score (arb.)")
        ax_grid.set_title(f"S4  grid search → {res.rbkg['confidence']}",
                          fontsize=10)
    else:
        ax_grid.text(0.5, 0.5, f"manual Rbkg = {res.rbkg['rbkg']:.2f} Å\n"
                               "(grid search skipped)",
                     transform=ax_grid.transAxes, ha="center", fontsize=10)
        ax_grid.set_title("S4  grid search", fontsize=10)
    _save(fig, outdir, "s4_rbkg.png", made)


def _s5(res: PipelineResult, outdir, made):
    ck, ft = res.chi_k, res.ft
    w = ft.kweight
    fig, (ax_k, ax_r) = plt.subplots(1, 2, figsize=(10, 3.6))
    ax_k.plot(ck.k, ck.k ** w * ck.chi, "k-", lw=0.8)
    ax_k.axvspan(ft.kmin, ft.kmax, color="orange", alpha=0.15,
                 label=f"FT window {ft.kmin:.0f}–{ft.kmax:.0f} Å$^{{-1}}$")
    ax_k.set_xlabel(r"$k$ ($\AA^{-1}$)")
    ax_k.set_ylabel(rf"$k^{w}\chi(k)$")
    ax_k.legend(frameon=False, fontsize=8, loc="upper right")
    ax_k.set_title("S5  weighted χ(k) + FT window", fontsize=10)
    ax_r.plot(ft.r, ft.magnitude, "k-", lw=1.0)
    if np.isfinite(res.r1):
        ax_r.axvline(res.r1, ls="--", color="r", lw=0.9)
        ax_r.text(res.r1 + 0.05, 0.9 * ax_r.get_ylim()[1],
                  f"first peak {res.r1:.3f} Å", color="r", fontsize=9)
    ax_r.set_xlim(0, 4.5)
    ax_r.set_xlabel(r"$R$ ($\AA$)")
    ax_r.set_ylabel(r"$|\chi(R)|$")
    ax_r.set_title("S5  |χ(R)| (not phase-corrected)", fontsize=10)
    _save(fig, outdir, "s5_chik_chir.png", made)


def _s6(res: PipelineResult, outdir, made):
    wt = res.wt
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    vmax = float(wt.magnitude.max()) or 1.0
    im = ax.contourf(wt.k_centers, wt.r, wt.magnitude.T / vmax, levels=21,
                     cmap="jet", vmin=0, vmax=1)
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label(r"$|W(R,k)|$ (normalized)", fontsize=9)
    ax.set_xlabel(r"$k$ ($\AA^{-1}$)")
    ax.set_ylabel(r"$R$ ($\AA$)")
    ax.set_ylim(0, wt.r.max())
    if np.isfinite(res.r_wt):
        ax.axhline(res.r_wt, color="k", ls="--", lw=0.8)
        ax.text(0.02, res.r_wt + 0.07,
                f"first shell  R≈{res.r_wt:.2f} Å",
                transform=ax.get_yaxis_transform(), fontsize=9,
                bbox=dict(facecolor="white", alpha=0.75, edgecolor="none",
                          pad=1.2))
    ax.text(0.98, 0.03, "R not phase-corrected", transform=ax.transAxes,
            ha="right", fontsize=8,
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none",
                      pad=1.2))
    ax.set_title(rf"S6  zero-carrier CWT ($\sigma$={wt.sigma})", fontsize=10)
    _save(fig, outdir, "s6_wavelet.png", made)


def _s7_ft_wt_joint(res: PipelineResult, outdir, made):
    """S6 扩展：FT–WT 联合诊断图。

    左图：|χ(R)|，标注 FT 第一峰 r1。
    右图：WT contour，叠加水平线 R=r1 与垂直线 k=k_wt，交点标出，
    直接回答「这个 FT 峰主要来自哪个 k」。

    为什么需要这个图（实测教训）：
        宽窗 [1.0,2.6] 在多壳体系会同时装进两条分支，WT argmax 会跳壳。
        例如 Fe3O4 的 Fe–O (R≈1.42, k≈4.5) 与 Fe–Fe (R≈2.72, k≈5.8)
        强度只差 3%；Booth 自吸收校正把低 k 抬高 7% 后，WT 峰从 2.60
        跳到 1.42，而 FT 峰纹丝不动。联合图能让用户肉眼核对：
        FT 峰位 r1 在 WT 上对应的 k 是否合理，从而判定 WT 选峰是否可信。
    """
    ck, ft, wt = res.chi_k, res.ft, res.wt
    fig, (ax_r, ax_k) = plt.subplots(1, 2, figsize=(10, 3.6))

    # ---- 左图：|χ(R)|，标注 FT 第一峰 ----
    ax_r.plot(ft.r, ft.magnitude, "k-", lw=1.0)
    if np.isfinite(res.r1):
        ax_r.axvline(res.r1, ls="--", color="r", lw=1.0)
        ax_r.text(res.r1 + 0.05, 0.9 * ax_r.get_ylim()[1],
                  f"FT r1 = {res.r1:.3f} Å", color="r", fontsize=9)
    ax_r.set_xlim(0, 4.5)
    ax_r.set_xlabel(r"$R$ ($\AA$)")
    ax_r.set_ylabel(r"$|\chi(R)|$")
    ax_r.set_title("FT  |χ(R)|", fontsize=10)

    # ---- 右图：WT contour，叠加 FT 峰位与 WT 峰位 ----
    vmax = float(wt.magnitude.max()) or 1.0
    im = ax_k.contourf(wt.k_centers, wt.r, wt.magnitude.T / vmax, levels=21,
                       cmap="jet", vmin=0, vmax=1)
    fig.colorbar(im, ax=ax_k, pad=0.02, label=r"$|W(R,k)|$ (normalized)")
    ax_k.set_xlabel(r"$k$ ($\AA^{-1}$)")
    ax_k.set_ylabel(r"$R$ ($\AA$)")
    ax_k.set_ylim(0, wt.r.max())

    # 水平线：FT 峰位 r1
    if np.isfinite(res.r1):
        ax_k.axhline(res.r1, color="r", ls="--", lw=0.9)
    # 垂直线：WT 峰位 k_wt
    if np.isfinite(res.k_wt):
        ax_k.axvline(res.k_wt, color="b", ls="-.", lw=0.9)

    # 交点标注
    if np.isfinite(res.r1) and np.isfinite(res.k_wt):
        ax_k.plot(res.k_wt, res.r1, "wo", ms=6, mec="k", mew=1.0,
                  label=f"FT-WT anchor\nR={res.r1:.2f} Å, k={res.k_wt:.1f} Å$^{{-1}}$")
        ax_k.legend(frameon=False, fontsize=8, loc="upper right")

    ax_k.set_title("WT  R–k response", fontsize=10)
    _save(fig, outdir, "s7_ft_wt_joint.png", made)


def wt_sweep(res: PipelineResult, sigmas=WT_SWEEP_SIGMAS,
             kweights=WT_SWEEP_KWEIGHTS, rmax: float = 4.0):
    """在 FT 第一峰锚定下，扫描多组 (σ, kweight) 的 WT 响应。

    为什么必须固定锚点：宽窗 argmax 会在 Fe3O4 的 Fe–O/Fe–Fe 分支间跳
    （见 locate_shell_near 注释）。这里用 FT 峰 r1 ±0.5 Å 锚定，回答的是
    “同一个 FT 峰在不同 σ / kweight 下对应哪个 k”，而不是重新选壳。

    为什么按 kweight 行分别归一：kweight=3 的幅度约是 kweight=1 的 k² ≈ 36 倍
    （k≈6 Å⁻¹ 实测），用单一全局最大值做色标会把 σ=1.0/1.5 的 kweight=1、2
    整块压成深蓝（实测 s8 图 kweight=1 全黑）。因此**同一行（固定 kweight）
    共享一个最大值**，使 σ 之间可比；不同 kweight 行之间幅度不可比，
    这是 k 加权的定义决定的，不是数据处理。
    """
    r1 = float(res.r1) if np.isfinite(res.r1) else 1.8
    ft = res.ft
    maps = []
    rows = []
    for sigma in sigmas:
        for kw in kweights:
            wt = larch_wt.continuous_wavelet_transform(
                res.chi_k, sigma=float(sigma), kmin=ft.kmin, kmax=ft.kmax,
                kweight=int(kw), rmax=rmax, normalize=False)
            k_p, r_p, amp = larch_wt.locate_shell_near(wt, r1, half_width=0.5)
            maps.append((sigma, kw, wt))
            rows.append({
                "sigma": float(sigma), "kweight": int(kw),
                "k_peak": round(k_p, 3), "r_peak": round(r_p, 3),
                "amplitude": round(amp, 6),
                "offset_from_ft": round(r_p - r1, 3),
            })
    gmax = max(float(wt.magnitude.max()) for _, _, wt in maps) or 1.0
    row_max = {}
    for _, kw, wt in maps:
        row_max[int(kw)] = max(row_max.get(int(kw), 0.0),
                               float(wt.magnitude.max()))
    return {"rows": rows, "maps": maps, "global_max": gmax,
            "row_max": row_max, "anchor_r": r1}


def _render_wt_sweep(sweep: dict, png: Path):
    """把 wt_sweep 的 3×3 网格画到 png（供诊断阶段与单独导出共用）。

    色标按**行**（固定 kweight）归一：同一行内 σ 可比，跨行幅度不可比。
    """
    maps = sweep["maps"]
    gmax = sweep["global_max"]
    row_max = sweep.get("row_max") or {}
    r1 = sweep["anchor_r"]
    k_lo = min(float(wt.k_centers.min()) for _, _, wt in maps)
    k_hi = max(float(wt.k_centers.max()) for _, _, wt in maps)
    fig, axes = plt.subplots(len(WT_SWEEP_SIGMAS), len(WT_SWEEP_KWEIGHTS),
                             figsize=(10.2, 8.4), sharex=True, sharey=True)
    for ax, (sigma, kw, wt) in zip(axes.ravel(), maps):
        scale = row_max.get(int(kw), gmax) or gmax
        mesh = ax.pcolormesh(wt.k_centers, wt.r,
                             wt.magnitude.T / scale, cmap="jet",
                             vmin=0, vmax=1, shading="auto")
        row = next(r for r in sweep["rows"]
                   if r["sigma"] == sigma and r["kweight"] == kw)
        ax.axhline(r1, color="w", ls="--", lw=0.8, alpha=0.9)
        ax.plot(row["k_peak"], row["r_peak"], "x", color="w", ms=8,
                mew=1.6)
        ax.text(0.03, 0.94,
                f"k={row['k_peak']:.1f}  R={row['r_peak']:.2f}\n"
                f"ΔR={row['offset_from_ft']:+.2f} Å",
                transform=ax.transAxes, va="top", fontsize=7.5, color="w",
                bbox=dict(facecolor="black", alpha=0.45, edgecolor="none",
                          pad=1.5))
        ax.set_title(rf"$\sigma$={sigma}, kweight={kw}", fontsize=9)
        ax.set_ylim(0, wt.r.max())
        ax.set_xlim(k_lo, k_hi)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$k$ ($\AA^{-1}$)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$R$ ($\AA$)")
    fig.suptitle(
        "WT sweep anchored at FT first peak "
        rf"($R$={r1:.2f} Å; color scaled per kweight row; dashed = FT peak)",
        fontsize=10)
    cbar = fig.colorbar(mesh, ax=axes, pad=0.015, fraction=0.025)
    cbar.set_label(r"$|W(R,k)|$ / max of same kweight row", fontsize=9)
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _write_wt_sweep_csv(rows: list, path: Path):
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["sigma", "kweight", "k_peak", "r_peak",
                            "amplitude", "offset_from_ft"])
        w.writeheader()
        w.writerows(rows)


def write_back_ft_comparison(res: PipelineResult, fit_out, outdir,
                             stem: str | None = None):
    """Back-FT 诊断：R 窗 + q 空间实验/模型相位与包络。

    有 S7 拟合时用 ``FitResult.model_chi_r`` 做同一 R 窗反变换，比较
    实验与模型的 Re χ(q)；无拟合时只给实验谱（仍可用于判断壳层相位）。
    输出 ``<stem>_backft.png``（默认 stem="s5b_backft"）。
    """
    _style()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fit = (fit_out or {}).get("fit")
    if fit is not None:
        rmin, rmax = float(fit.rmin), float(fit.rmax)
    else:
        rmin = max(1.0, float(res.r1) - 0.5)
        rmax = min(3.4, float(res.r1) + 0.5)
    if not np.isfinite(rmin) or rmax <= rmin:
        rmin, rmax = 1.0, 2.6
    bft = athena.back_ft_to_q_space(res.ft.r, res.ft.chi_r, rmin, rmax)
    bft_model = None
    if fit is not None and getattr(fit, "model_chi_r", None) is not None:
        bft_model = athena.back_ft_to_q_space(
            fit.model_r, fit.model_chi_r, rmin, rmax)

    fig, (ax_r, ax_q) = plt.subplots(1, 2, figsize=(10.4, 3.9))
    rmask = (res.ft.r >= rmin) & (res.ft.r <= rmax)
    ax_r.plot(res.ft.r, res.ft.magnitude, color="#1b6ca8", lw=1.3,
              label="data")
    if fit is not None and getattr(fit, "model_chi_r", None) is not None:
        ax_r.plot(fit.model_r, np.abs(fit.model_chi_r), color="#d1495b",
                  lw=1.1, ls="--", label="model")
    ax_r.axvspan(rmin, rmax, color="#f0a028", alpha=0.12,
                 label=f"back-FT window {rmin:.2f}–{rmax:.2f} Å")
    ax_r.set_xlim(0, min(6.0, float(res.ft.r.max())))
    ax_r.set_xlabel(r"$R$ ($\AA$)")
    ax_r.set_ylabel(r"$|\chi(R)|$")
    ax_r.set_title("R space: filtered shell", fontsize=10)
    ax_r.legend(frameon=False, fontsize=8)

    ax_q.plot(bft.q, bft.real, color="#1b6ca8", lw=1.2,
              label=r"Re $\chi(q)$")
    ax_q.plot(bft.q, bft.magnitude, color="#68737f", lw=0.9, alpha=0.8,
              label=r"$|\chi(q)|$")
    ax_q.plot(bft.q, -bft.magnitude, color="#68737f", lw=0.9, alpha=0.5)
    if bft_model is not None:
        ax_q.plot(bft_model.q, bft_model.real, color="#d1495b", lw=1.0,
                  ls="--", label=r"Re $\chi_{model}(q)$")
    ax_q.set_xlim(0, min(16.0, float(bft.q.max())))
    ax_q.set_xlabel(r"$q$ ($\AA^{-1}$)")
    ax_q.set_ylabel(r"filtered $\chi$")
    ax_q.set_title("q space: phase + envelope", fontsize=10)
    ax_q.legend(frameon=False, fontsize=8, ncol=2)
    fig.suptitle(
        f"S5b  Back-FT of R∈[{rmin:.2f},{rmax:.2f}] Å "
        r"(Hanning, $e^{+2iqR}$)", fontsize=10)
    png = outdir / (f"{stem}.png" if stem else "s5b_backft.png")
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return png, bft, bft_model


def _s8_wt_sweep(res: PipelineResult, outdir, made):
    """S8：多 kweight / σ 的 WT 对比（同一锚点、按 kweight 行归一色标）。"""
    sweep = wt_sweep(res)
    png = Path(outdir) / "s8_wt_sweep.png"
    _render_wt_sweep(sweep, png)
    made.append(png)
    csv_path = Path(outdir) / "wt_sweep.csv"
    _write_wt_sweep_csv(sweep["rows"], csv_path)
    made.append(csv_path)


def write_wt_sweep(res: PipelineResult, outdir, stem: str | None = None):
    """单独输出 WT 扫描图 + CSV（供 xafs_report --wt-sweep 调用）。"""
    _style()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sweep = wt_sweep(res)
    png = outdir / (f"{stem}.png" if stem else "s8_wt_sweep.png")
    csv_path = outdir / (f"{stem}.csv" if stem else "wt_sweep.csv")
    _render_wt_sweep(sweep, png)
    _write_wt_sweep_csv(sweep["rows"], csv_path)
    return png, csv_path, sweep["rows"]


def _fit_result(fit_out):
    """兼容三种传入：FitResult 本身 / fit_identified 的 dict / run() 的外层 dict。"""
    obj = fit_out
    for _ in range(3):
        if obj is None or not isinstance(obj, dict):
            break
        obj = obj.get("fit")
    if obj is not None and hasattr(obj, "paths") and hasattr(obj, "model_k"):
        return obj
    return None


def _ft_anchor_near_ft(r, mag, r_bond: float) -> float:
    """把键长 R_bond 映射到 FT/WT 的表观 R 峰（相位位移约 −0.3～−0.5 Å）。

    在 [R_bond−0.95, R_bond−0.15] 内取 |χ(R)| 最大处；若窗口内无点，
    退回 R_bond−0.4。这样第一壳 Fe–O（2.00 → 1.47）与第二壳
    Fe–Fe（3.07 → 2.70）都能锚到自己的 FT 峰，不会串壳。
    """
    r = np.asarray(r, dtype=float)
    mag = np.asarray(mag, dtype=float)
    lo = max(0.8, float(r_bond) - 0.95)
    hi = min(float(r.max()), float(r_bond) - 0.15)
    if hi > lo:
        sel = (r >= lo) & (r <= hi)
        if int(sel.sum()) > 0:
            return float(r[sel][int(np.argmax(mag[sel]))])
    return float(max(0.8, float(r_bond) - 0.4))


def write_wt_shell_profiles(res: PipelineResult, fit_out, outdir,
                            stem: str | None = None,
                            ref_k_peaks: dict[str, float] | None = None):
    """B3-1/B5-1 壳层一维小波包络 + 实验/模型二维 WT 残差三联图。

    为什么必须做“实验 vs 模型”而不只看实验图：二维 WT 伪彩图里
    模型残差会把「哪些 R–k 区域没拟合好」直接暴露出来；单壳层的一维
    切片则回答第二壳到底是轻元素（k≈4–5 Å⁻¹）、3d 金属（6.5–8.5）
    还是重金属（>9）。两者共用同一组 WT 参数，禁止各自的归一化。

    ``ref_k_peaks``（B5-2）：候选散射元素的理论 k_peak 标尺
    （``structure_id.reference_element_k_peaks`` 的输出）。传入时在
    每个切片顶部画灰色点划线，避免把高 σ² 金属壳误判为轻配体。

    输出:
        <stem>_wt_profiles.png : N 个壳层切片 + 1 个二维残差面板
        <stem>_wt_profiles.csv : 每个壳层的 k / data / model 包络与相位
    """
    _style()
    fit = _fit_result(fit_out)
    if fit is None or getattr(fit, "model_k", None) is None:
        raise XAFSConfigError(
            "WT", "write_wt_shell_profiles 需要含 model_k 的 S7 FitResult"
                 "（先跑 S7 拟合）")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = stem or (res.stem or "wt_profiles")

    # ---- 用拟合路径的键长定位第一/第二壳，再映射到 FT 表观 R ----
    ss_paths = [p for p in sorted(fit.paths, key=lambda q: q.reff)
                if int(getattr(p, "nleg", 2)) == 2]
    shells = []
    seen = []
    for p in ss_paths:
        r_bond = float(p.r)
        if any(abs(r_bond - s) < 0.05 for s in seen):
            continue
        seen.append(r_bond)
        shells.append((p, _ft_anchor_near_ft(res.ft.r, res.ft.magnitude,
                                             r_bond)))
        if len(shells) >= 2:
            break
    if not shells:
        raise XAFSConfigError(
            "WT", "拟合结果里没有单散射路径，无法生成壳层 WT 切片")

    residual = larch_wt.compute_model_wt_residual(
        res.chi_k, fit.model_k,
        sigma=float(res.wt.sigma),
        kmin=float(fit.kmin),
        kmax=float(min(fit.kmax, np.max(res.chi_k.k))),
        kweight=int(res.ft.kweight),
        rmax=float(res.wt.r.max()))
    data_map = residual["wt_data"]
    model_map = residual["wt_model"]
    scale = residual["scale"] or 1.0

    profiles = []
    for p, r_anchor in shells:
        d = larch_wt.extract_k_profile(data_map, r_anchor, half_width=0.25)
        m = larch_wt.extract_k_profile(model_map, r_anchor, half_width=0.25)
        profiles.append({
            "label": p.label, "r_bond": float(p.r),
            "r_anchor": float(r_anchor),
            "data": d, "model": m,
            "k_peak_data": float(d["k_peak"]),
            "k_peak_model": float(m["k_peak"]),
        })

    # ---- CSV：长表，所有壳层共用 k 网格（由同一个 WT 参数决定） ----
    csv_path = outdir / f"{prefix}_wt_profiles.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["shell_index", "shell_label", "r_bond_A",
                         "r_anchor_A", "k_invA",
                         "data_ridge", "model_ridge",
                         "data_slice", "model_slice",
                         "data_phase_ridge", "model_phase_ridge"])
        for i, item in enumerate(profiles, start=1):
            d, m = item["data"], item["model"]
            d_phase = (d.get("phase_ridge") if d.get("phase_ridge") is not None
                       else np.zeros_like(d["k"]))
            m_phase = (m.get("phase_ridge") if m.get("phase_ridge") is not None
                       else np.zeros_like(m["k"]))
            for kk, dr, mr, ds, ms, dp, mp in zip(
                    d["k"], d["ridge_profile"] / scale,
                    m["ridge_profile"] / scale,
                    d["slice_profile"] / scale,
                    m["slice_profile"] / scale, d_phase, m_phase):
                writer.writerow([i, item["label"], f"{item['r_bond']:.4f}",
                                 f"{item['r_anchor']:.4f}", f"{float(kk):.4f}",
                                 f"{float(dr):.8e}", f"{float(mr):.8e}",
                                 f"{float(ds):.8e}", f"{float(ms):.8e}",
                                 f"{float(dp):.6f}", f"{float(mp):.6f}"])

    # ---- 图：N 个壳层切片 + 残差 ----
    n_shells = len(profiles)
    fig, axes = plt.subplots(
        1, n_shells + 1, figsize=(4.1 * (n_shells + 1), 3.6),
        gridspec_kw={"width_ratios": [1.0] * n_shells + [1.25]})
    axes = np.atleast_1d(axes)
    for ax, item in zip(axes[:n_shells], profiles):
        d, m = item["data"], item["model"]
        ax.plot(d["k"], d["ridge_profile"] / scale, color="#1b6ca8",
                lw=1.7, label="experiment")
        ax.plot(m["k"], m["ridge_profile"] / scale, color="#d1495b",
                lw=1.5, ls="--", label="fit model")
        ax.axvline(d["k_peak"], color="#1b6ca8", ls=":", lw=0.9)
        ax.axvline(m["k_peak"], color="#d1495b", ls=":", lw=0.9)
        ax.set_xlabel(r"$k$ ($\AA^{-1}$)")
        ax.set_ylabel(r"$|W(k; R_i)|$ (norm.)")
        ax.set_title(f"{item['label']}  R$_i$={item['r_anchor']:.2f} Å\n"
                     f"k$_{{peak}}$: data={d['k_peak']:.2f}, "
                     f"model={m['k_peak']:.2f} Å$^{{-1}}$",
                     fontsize=9)
        ax.legend(frameon=False, fontsize=8)
        ax.set_xlim(float(d["k"].min()), float(d["k"].max()))
        if ref_k_peaks:
            # B5-2 理论元素标尺：灰色点划线，标签贴在轴顶部
            ymin, ymax = ax.get_ylim()
            for el_name, k_ref in sorted(ref_k_peaks.items(),
                                         key=lambda kv: float(kv[1])):
                if not (d["k"].min() <= k_ref <= d["k"].max()):
                    continue
                ax.axvline(float(k_ref), color="#777777", ls="-.", lw=0.8,
                           alpha=0.85)
                ax.text(float(k_ref), ymax * 0.98, str(el_name),
                        rotation=90, ha="right", va="top", fontsize=7,
                        color="#555555")

    ax_res = axes[-1]
    kc = data_map.k_centers
    rr = data_map.r
    rmap = residual["residual_map"]
    vmax = float(np.max(np.abs(rmap))) or 1.0
    mesh = ax_res.pcolormesh(kc, rr, rmap.T, cmap="RdBu_r",
                             vmin=-vmax, vmax=vmax, shading="auto")
    fig.colorbar(mesh, ax=ax_res, pad=0.02,
                 label=r"$|W_{data}|-|W_{model}|$ (norm.)")
    ax_res.set_xlabel(r"$k$ ($\AA^{-1}$)")
    ax_res.set_ylabel(r"$R$ ($\AA$)")
    ax_res.set_title(f"WT residual (2D $R_{{WT}}$={residual['r_factor']:.4f}"
                     + ("" if residual.get("complex_r_factor") is None else
                        f", complex={residual['complex_r_factor']:.4f}") + ")",
                     fontsize=9)
    fig.tight_layout()
    png = outdir / f"{prefix}_wt_profiles.png"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return {
        "png": png,
        "csv": csv_path,
        "shells": [{"label": it["label"], "r_bond": it["r_bond"],
                    "r_anchor": it["r_anchor"],
                    "k_peak_data": it["k_peak_data"],
                    "k_peak_model": it["k_peak_model"]}
                   for it in profiles],
        "r_factor": float(residual["r_factor"]),
        "complex_r_factor": residual.get("complex_r_factor"),
        "residual": residual,
    }


STAGES = (_s1, _s2, _s3, _s4, _s5, _s6, _s7_ft_wt_joint, _s8_wt_sweep)


def write_diagnostics(res: PipelineResult, outdir) -> list:
    """逐阶段画诊断图，返回生成的 PNG 路径列表。"""
    _style()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    made: list = []
    for stage in STAGES:
        stage(res, outdir, made)
    return made


# ================================================================
# B14-2 审稿级诊断：相移校正真键长 FT/WT 与二维 χ² 置信椭圆
# ================================================================

def _fit_payload(fit_out) -> dict:
    """兼容 ``xafs_fit.run`` 外层输出与 ``run["fit"]`` 内层字典。"""
    if not isinstance(fit_out, dict):
        return {}
    inner = fit_out.get("fit")
    if isinstance(inner, dict) and ("mode" in inner or "fit" in inner):
        return inner
    return fit_out


def _first_shell_datfile(fit_payload: dict):
    """取第一壳 FEFF 路径：优先 S7 path_files，回退 L1 合成团簇缓存。"""
    files = fit_payload.get("path_files") or []
    if files:
        return Path(files[0])
    l1 = fit_payload.get("l1") or {}
    if not l1:
        return None
    import structure_id
    import structure_lib
    motif = structure_id.synthetic_shell_motif(
        str(fit_payload.get("absorber") or l1.get("absorber") or "Fe"),
        ligand=str(l1.get("ligand", "O")),
        geometry=str(l1.get("geometry", "oct")),
        r0=float(l1.get("r0", 2.0)))
    cache = structure_lib.ensure_feff(motif)
    rows = artemis.describe_paths(sorted(cache.glob("feff*.dat")), 40)
    ss = [r for r in rows if r["legs"] == 2]
    if not ss:
        return None
    r0 = float(l1.get("r0", ss[0]["reff"]))
    best = min(ss, key=lambda r: abs(float(r["reff"]) - r0))
    return cache / best["file"]


def _peak_in_window(r, mag, lo=0.8, hi=3.2):
    r = np.asarray(r, dtype=float)
    mag = np.asarray(mag, dtype=float)
    m = (r >= lo) & (r <= hi) & np.isfinite(mag)
    if not m.any():
        return float("nan")
    rr = r[m]
    return float(rr[int(np.argmax(mag[m]))])


def write_phase_corrected_diagnostics(res: PipelineResult, fit_out, outdir,
                                      stem: str | None = None) -> dict:
    """B14-2 相移校正真键长诊断（FT + 二维 WT）。

    用第一壳 FEFF 路径的总相移 φ_tot(k) 对 χ(k) 做 e^{−iφ} 校正，峰位
    从表观 R（Reff−0.3~0.5 Å）回到真实键长；导出 CSV 与双栏 PNG。
    """
    _style()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = stem or getattr(res, "stem", "") or "diagnostics"
    payload = _fit_payload(fit_out)
    fit = payload.get("fit")
    datfile = _first_shell_datfile(payload)
    if datfile is None or not Path(datfile).is_file():
        raise XAFSConfigError(
            "diag", "相移校正需要第一壳 FEFF 路径（S7 path_files 或 L1 缓存）")
    chi_k = res.chi_k
    k = np.asarray(chi_k.k, dtype=float)
    ft = res.ft
    kmin, kmax = float(ft.kmin), float(ft.kmax)
    kw = int(ft.kweight)
    dk = float(getattr(ft, "dk", 1.0))
    r_fit = (float(fit.paths[0].r) if fit is not None and fit.paths
             else float(payload.get("l1", {}).get("R", np.nan)))
    # FEFF 的 pha 表以 FEFF 内部 E0 为原点，而管线 χ(k) 以 e0_sample 为原点。
    # 直接把 φ_tot(k) 套在数据 k 上会残留 ΔE0 相移（Fe 实测峰位 1.90 →
    # 2.00 Å 差 0.10 Å）。改用 k_feff = sqrt(k² − 0.2625·ΔE) 在 FEFF 的
    # k 网格上取相移，ΔE = e0_sample + ΔE0_fit − e0_edge。
    delta_e0 = 0.0
    if fit is not None:
        delta_e0 = float(getattr(fit, "delta_e0", 0.0) or 0.0)
    elif payload.get("l1"):
        delta_e0 = float(payload["l1"].get("delta_e0", 0.0) or 0.0)
    e0_sample = float((res.e0 or {}).get("e0_sample", np.nan))
    e0_edge = None
    _symbol = getattr(res, "absorber_symbol", None) or payload.get("absorber")
    if _symbol:
        try:
            from absorber import get_absorber
            e0_edge = getattr(get_absorber(str(_symbol)), "e0_edge", None)
        except Exception:  # noqa: BLE001 - 缺标定表时退化为不补偿
            e0_edge = None
    delta_e = np.nan
    if e0_edge is not None and np.isfinite(e0_sample):
        delta_e = e0_sample + delta_e0 - float(e0_edge)
    if np.isfinite(delta_e):
        k_feff = np.sqrt(np.maximum(k * k - 0.2625 * delta_e, 1e-6))
        phi = artemis.path_total_phase(datfile, k_feff)
    else:
        phi = artemis.path_total_phase(datfile, k)
    raw = athena.ft_to_r_space(k, chi_k.chi, None, kmin, kmax, kw, dk)
    corr = athena.ft_to_r_space(k, chi_k.chi, None, kmin, kmax, kw, dk,
                                phase_correct=phi)
    grid = np.asarray(raw.r, dtype=float)
    corr_r = np.interp(grid, np.asarray(corr.r, dtype=float), corr.chi_r.real)
    corr_i = np.interp(grid, np.asarray(corr.r, dtype=float), corr.chi_r.imag)
    mag_raw = np.abs(raw.chi_r)
    mag_corr = np.hypot(corr_r, corr_i)
    # 相移校正后的真键长峰在 r_fit 附近；窗口随 r_fit 走，避免抓到
    # 远壳或低 R 泄漏（固定 [0.8,3.2] 在 NiO 上会选中 2.55 Å 的 Ni–Ni）。
    lo = max(0.8, r_fit - 0.7) if np.isfinite(r_fit) else 0.8
    hi = r_fit + 0.4 if np.isfinite(r_fit) else 3.2
    r_peak_raw = _peak_in_window(grid, mag_raw, lo, hi)
    r_peak_corrected = _peak_in_window(grid, mag_corr, lo, hi)

    csv_path = outdir / f"{stem}_ft_phase_corrected.csv"
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["R_A", "mag_raw", "real_raw",
                    "mag_phase_corrected", "real_phase_corrected"])
        for i in range(grid.size):
            w.writerow([f"{grid[i]:.8g}", f"{mag_raw[i]:.10g}",
                        f"{raw.chi_r.real[i]:.10g}",
                        f"{mag_corr[i]:.10g}", f"{corr_r[i]:.10g}"])

    wt = larch_wt.continuous_wavelet_transform(
        chi_k, sigma=1.5, kmin=kmin, kmax=min(12.0, kmax),
        kweight=kw, rmax=3.5, phase_correct=phi)

    png_path = outdir / f"{stem}_phase_corrected.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    ax = axes[0]
    ax.plot(grid, mag_raw, color="#9aa0a6", lw=1.2,
            label="raw |χ(R)|")
    ax.plot(grid, mag_corr, color="#1b6ca8", lw=1.5,
            label="phase-corrected |χ(R)|")
    if np.isfinite(r_fit):
        ax.axvline(r_fit, color="#d1495b", ls="--", lw=1.0,
                   label=f"fit R₁ = {r_fit:.3f} Å")
    ax.set_xlim(0.5, 3.2)
    ax.set_xlabel("R (Å)")
    ax.set_ylabel("|χ(R)|")
    ax.set_title("Phase-corrected FT: apparent → true bond length")
    ax.legend(fontsize=8, frameon=False)
    ax2 = axes[1]
    mesh = ax2.pcolormesh(wt.r, wt.k_centers, wt.magnitude,
                          shading="auto", cmap="viridis")
    ax2.set_xlim(0.5, 3.2)
    ax2.set_xlabel("R (Å)")
    ax2.set_ylabel("k (Å⁻¹)")
    ax2.set_title("Phase-corrected wavelet |W(R,k)|")
    fig.colorbar(mesh, ax=ax2, label="normalized |W|")
    fig.tight_layout()
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return {"csv": csv_path, "png": png_path,
            "r_peak_raw": r_peak_raw,
            "r_peak_corrected": r_peak_corrected,
            "r_fit": r_fit,
            "delta_e": delta_e,
            "phase_source": getattr(artemis.path_total_phase, "last_source", None)}


def write_chi2_contour_diagnostics(res: PipelineResult, fit_out, outdir,
                                   stem: str | None = None,
                                   n_grid: int = 15) -> dict | None:
    """B14-2 二维 χ² 置信椭圆（第一壳 N–σ²，固定其余最优参数）。

    Δχ² = 2.30/6.18/11.83 对应二维 1σ/2σ/3σ；导出 CSV 与 PNG。
    """
    _style()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = stem or getattr(res, "stem", "") or "diagnostics"
    payload = _fit_payload(fit_out)
    fit = payload.get("fit")
    files = payload.get("path_files") or []
    if fit is None or not fit.paths or not files:
        return None
    if len(files) != len(fit.paths):
        return None
    paths = [artemis.load_feff_path(f) for f in files]
    for p, pf in zip(paths, fit.paths):
        p.degen = float(pf.degeneracy)
        artemis.set_path_parameters(
            p, s02=float(fit.s02), delta_r=float(pf.delta_r),
            sigma2=float(pf.sigma2), delta_e0=float(fit.delta_e0))
    p0, pf0 = paths[0], fit.paths[0]

    x_name = "N" if pf0.n_err is not None else "deltar"
    y_name = "sigma2"
    x0 = (float(pf0.degeneracy) if x_name == "N" else float(pf0.delta_r))
    x_err = pf0.n_err if x_name == "N" else pf0.delta_r_err
    y0 = float(pf0.sigma2)
    y_err = pf0.sigma2_err
    x_half = max(3.0 * float(x_err or 0.0), 0.20 * max(abs(x0), 1e-3))
    y_half = max(3.0 * float(y_err or 0.0), 0.20 * max(abs(y0), 1e-4))
    xv = np.linspace(max(0.0, x0 - x_half), x0 + x_half, int(n_grid))
    yv = np.linspace(max(1e-4, y0 - y_half), y0 + y_half, int(n_grid))

    k = np.asarray(res.chi_k.k, dtype=float)
    r_data = np.asarray(res.ft.r, dtype=float)
    data = np.asarray(res.ft.chi_r, dtype=complex)
    mask = (r_data >= float(fit.rmin)) & (r_data <= float(fit.rmax))
    n_pts = int(mask.sum())
    if n_pts < 8:
        return None
    eps = fit.epsilon_r
    denom = float(np.sum(np.abs(data[mask]) ** 2)) or 1.0

    def chi2_at(**params):
        n_val = float(params.get(
            "N", pf0.degeneracy if x_name != "N" else x0))
        dr_val = float(params.get(
            "deltar", pf0.delta_r if x_name != "deltar" else x0))
        s2_val = float(params.get("sigma2", y0))
        p0.degen = n_val
        artemis.set_path_parameters(
            p0, s02=float(fit.s02), delta_r=dr_val,
            sigma2=s2_val, delta_e0=float(fit.delta_e0))
        model_k = artemis.exafs_model_k(paths, k)
        model_ft = athena.ft_to_r_space(
            k, model_k, None, float(fit.kmin), float(fit.kmax),
            int(res.ft.kweight), float(getattr(res.ft, "dk", 1.0)))
        m_re = np.interp(r_data, model_ft.r, model_ft.chi_r.real)
        m_im = np.interp(r_data, model_ft.r, model_ft.chi_r.imag)
        diff2 = float(np.sum((data[mask].real - m_re[mask]) ** 2
                             + (data[mask].imag - m_im[mask]) ** 2))
        if eps is not None and np.isfinite(eps) and eps > 0:
            return float(fit.n_independent * diff2
                         / (2.0 * n_pts * eps ** 2))
        return float(fit.n_independent * diff2 / (2.0 * n_pts * denom))

    contour = artemis.profile_chi2_contour(
        chi2_at, x_name, y_name, xv, yv)
    if contour.get("min_chi2") is None:
        return None
    csv_path = outdir / f"{stem}_chi2_contour.csv"
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow([x_name, y_name, "chi2"])
        grid = contour["chi2_grid"]
        for iy, y in enumerate(contour["y"]):
            for ix, x in enumerate(contour["x"]):
                w.writerow([f"{x:.8g}", f"{y:.8g}", f"{grid[iy, ix]:.10g}"])

    png_path = outdir / f"{stem}_chi2_contour.png"
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    cs = ax.contourf(contour["x"], contour["y"], contour["chi2_grid"],
                     levels=24, cmap="viridis")
    if contour["levels"]:
        try:
            ax.contour(contour["x"], contour["y"], contour["chi2_grid"],
                       levels=contour["levels"], colors="w",
                       linewidths=0.8, linestyles="--")
        except ValueError:
            pass
    ax.plot([x0], [y0], marker="*", color="#d1495b", ms=12,
            label="optimum")
    ax.set_xlabel("N" if x_name == "N" else r"$\Delta R$ (Å)")
    ax.set_ylabel(r"$\sigma^2$ (Å$^2$)")
    ax.set_title(r"$\chi^2$ confidence contour "
                 r"($\Delta\chi^2$=2.30/6.18/11.83)")
    ax.legend(fontsize=8, frameon=False)
    fig.colorbar(cs, ax=ax, label=r"$\chi^2$")
    fig.tight_layout()
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return {"csv": csv_path, "png": png_path,
            "x_name": x_name, "y_name": y_name,
            "min_chi2": float(contour["min_chi2"]),
            "chi2_grid": contour["chi2_grid"],
            "x": contour["x"], "y": contour["y"],
            "levels": contour["levels"]}


# ================================================================
# B15-1 原位（Operando）序列诊断：浓度/物种分数/描述符四宫格
# ================================================================

def write_operando_diagnostics(operando_out: dict, outdir,
                               stem: str = "operando") -> dict:
    """把 ``xafs_lcf.run_operando_series`` 的结果写成 CSV + 四宫格图。

    产物：
      - ``<stem>_operando_concentrations.csv``：每条谱的 MCR 组分分数、
        LCF 标样权重与 E0/半高能/白线/FT 第一峰描述符；
      - ``<stem>_operando_components.csv``：能量网格与 MCR 纯组分谱；
      - ``<stem>_operando_summary.png``：四宫格演化图。
    """
    _style()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    names = list(operando_out.get("names") or [])
    n = len(names)
    mcr = operando_out.get("mcr_als") or {}
    C = np.asarray(mcr.get("concentrations", np.zeros((0, n))),
                   dtype=float)
    if C.ndim != 2 or C.shape[1] != n:
        C = np.zeros((int(operando_out.get("n_components") or 0), n))
    k = int(C.shape[0])
    S = np.asarray(mcr.get("components", np.zeros((0, k))), dtype=float)
    xg = np.asarray(mcr.get("x", []), dtype=float)
    std_names = list(operando_out.get("accepted_standards") or [])
    lcf_series = list(operando_out.get("lcf_series") or [])
    descriptors = list(operando_out.get("descriptors") or [])

    # ---- 浓度 / 权重 / 描述符总表 ----
    conc_csv = outdir / f"{stem}_operando_concentrations.csv"
    header = (["index", "name"]
              + [f"mcr_component_{i + 1}" for i in range(k)]
              + [f"lcf_{s}" for s in std_names]
              + ["e0_eV", "e_half_height_eV",
                 "white_line_height", "ft_r1_A"])
    with open(conc_csv, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for i in range(n):
            row = [i, names[i]]
            row += [f"{C[j, i]:.10g}" if j < C.shape[0] else ""
                    for j in range(k)]
            weights = (lcf_series[i].get("weights")
                       if i < len(lcf_series) else {}) or {}
            row += [f"{float(weights.get(s, np.nan)):.10g}"
                    if weights.get(s) is not None else ""
                    for s in std_names]
            d = descriptors[i] if i < len(descriptors) else {}
            for key in ("e0", "e_half_height",
                        "white_line_height", "ft_r1"):
                val = d.get(key)
                row.append("" if val is None else f"{float(val):.10g}")
            w.writerow(row)

    # ---- MCR 纯组分谱 ----
    comp_csv = outdir / f"{stem}_operando_components.csv"
    with open(comp_csv, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["energy_eV"]
                   + [f"mcr_component_{i + 1}" for i in range(k)])
        for i in range(xg.size):
            w.writerow([f"{xg[i]:.10g}"]
                       + [f"{S[i, j]:.10g}" if i < S.shape[0] else ""
                          for j in range(k)])

    # ---- 四宫格演化图 ----
    png_path = outdir / f"{stem}_operando_summary.png"
    idx = np.arange(n)
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))
    ax = axes[0, 0]
    for j in range(k):
        ax.plot(idx, C[j, :], marker="o", lw=1.4,
                label=f"MCR component {j + 1}")
    ax.set_xlabel("spectrum index")
    ax.set_ylabel("MCR concentration")
    ax.set_title("MCR-ALS species fractions")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[0, 1]
    if std_names:
        for s in std_names:
            vals = [float(((lcf_series[i].get("weights") or {}).get(s, np.nan))
                          if i < len(lcf_series) else np.nan)
                    for i in range(n)]
            ax.plot(idx, vals, marker="s", lw=1.2, label=s)
        ax.set_ylabel("LCF weight")
        ax.legend(fontsize=8, frameon=False)
    else:
        ax.text(0.5, 0.5, "no accepted standards",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_ylabel("LCF weight")
    ax.set_xlabel("spectrum index")
    ax.set_title("LCF standard weights (SPOIL-filtered)")

    ax = axes[1, 0]
    e0 = np.asarray([np.nan if (i >= len(descriptors)
                                or descriptors[i].get("e0") is None)
                     else float(descriptors[i]["e0"])
                     for i in range(n)], dtype=float)
    ehh = np.asarray([np.nan if (i >= len(descriptors)
                                 or descriptors[i].get("e_half_height") is None)
                      else float(descriptors[i]["e_half_height"])
                      for i in range(n)], dtype=float)
    wl = np.asarray([np.nan if (i >= len(descriptors)
                                or descriptors[i].get("white_line_height") is None)
                     else float(descriptors[i]["white_line_height"])
                     for i in range(n)], dtype=float)
    ax.plot(idx, e0 - np.nanmean(e0) if np.isfinite(e0).any() else e0,
            marker="o", color="#1b6ca8", label="E0 − mean (eV)")
    ax.plot(idx, ehh - np.nanmean(ehh) if np.isfinite(ehh).any() else ehh,
            marker="^", color="#2a9d8f", label="E1/2 − mean (eV)")
    ax.set_xlabel("spectrum index")
    ax.set_ylabel("energy shift (eV)")
    ax2 = ax.twinx()
    ax2.plot(idx, wl, marker="d", color="#d1495b", lw=1.1,
             label="white-line height")
    ax2.set_ylabel("white-line height", color="#d1495b")
    ax.set_title("Valence / structure descriptors")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, frameon=False, loc="best")

    ax = axes[1, 1]
    for j in range(k):
        ax.plot(xg, S[:, j], lw=1.4, label=f"component {j + 1}")
    ax.set_xlabel("energy relative to E0 (eV)")
    ax.set_ylabel("normalized pure component")
    ax.set_title("MCR-ALS pure component spectra (energy-sorted)")
    ax.legend(fontsize=8, frameon=False)

    fig.tight_layout()
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return {"concentrations_csv": conc_csv,
            "components_csv": comp_csv,
            "png": png_path}


def write_operando_exafs_diagnostics(exafs_series_out: dict, outdir,
                                     stem: str = "operando_exafs") -> dict:
    """B19-1：原位 EXAFS 序列 N(t)/R(t)/σ²(t) 轨迹落盘。

    产物：
      - ``<stem>_exafs_evolution.csv``：每个时间点一行，含各壳层
        N/N_err/R/R_err/σ²；
      - ``<stem>_exafs_evolution.png``：上 N(t)±σ_N、下 R(t)±σ_R。
    """
    _style()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    names = list(exafs_series_out.get("names") or [])
    shells = list(exafs_series_out.get("shells") or [])
    traj = exafs_series_out.get("trajectories") or {}
    r_factors = list(exafs_series_out.get("r_factors") or [])
    n = len(names)
    delta_e0 = exafs_series_out.get("shared_delta_e0")

    csv_path = outdir / f"{stem}_exafs_evolution.csv"
    header = ["index", "name", "r_factor", "delta_e0"]
    for label in shells:
        header += [f"{label}_N", f"{label}_N_err",
                   f"{label}_R", f"{label}_R_err", f"{label}_sigma2"]
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for i in range(n):
            rf = r_factors[i] if i < len(r_factors) else None
            row = [i, names[i],
                   "" if rf is None else f"{float(rf):.10g}",
                   "" if delta_e0 is None else f"{float(delta_e0):.10g}"]
            for label in shells:
                t = traj.get(label) or {}

                def _val(key):
                    vals = t.get(key) or []
                    if i >= len(vals) or vals[i] is None:
                        return ""
                    return f"{float(vals[i]):.10g}"

                row += [_val("N"), _val("N_err"), _val("R"),
                        _val("R_err"), _val("sigma2")]
            w.writerow(row)

    png_path = outdir / f"{stem}_exafs_evolution.png"
    idx = np.arange(n)
    fig, axes = plt.subplots(2, 1, figsize=(8.6, 7.0), sharex=True)
    colors = ("#1b6ca8", "#d1495b", "#2a9d8f", "#e9c46a", "#6a4c93")
    for j, label in enumerate(shells):
        t = traj.get(label) or {}
        c = colors[j % len(colors)]
        nv = np.asarray([np.nan if v is None else float(v)
                         for v in (t.get("N") or [])], dtype=float)
        ne = np.asarray([0.0 if v is None else float(v)
                         for v in (t.get("N_err") or [])], dtype=float)
        rv = np.asarray([np.nan if v is None else float(v)
                         for v in (t.get("R") or [])], dtype=float)
        re = np.asarray([0.0 if v is None else float(v)
                         for v in (t.get("R_err") or [])], dtype=float)
        axes[0].errorbar(idx, nv, yerr=ne, marker="o", lw=1.3,
                         color=c, capsize=2.5, label=label)
        axes[1].errorbar(idx, rv, yerr=re, marker="s", lw=1.3,
                         color=c, capsize=2.5, label=label)
    axes[0].set_ylabel("N(t)")
    axes[0].set_title("Operando EXAFS: coordination number evolution")
    axes[0].legend(fontsize=8, frameon=False)
    axes[1].set_xlabel("spectrum index")
    axes[1].set_ylabel("R(t) (Å)")
    axes[1].set_title("Operando EXAFS: bond length evolution")
    axes[1].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return {"csv": csv_path, "png": png_path}

