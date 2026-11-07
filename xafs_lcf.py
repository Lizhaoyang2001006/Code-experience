# -*- coding: utf-8 -*-
"""
xafs_lcf.py — XANES / χ(k) 线性组合拟合（LCF）与 PCA 主成分数判定
====================================================================
对标 Athena 的 Linear Combination Fitting 与六标准法（LCF + PCA），
解决**多物相机械混合**样品的相比例定量（区别于单晶多位点的
``--fit-weights``：后者是同一吸收位点的组分权重，本模块是不同端元
谱的线性组合）。

依赖：``xafs_core``、``athena``。不 import 拟合链/报告层，保持单向依赖。
"""
from __future__ import annotations

import itertools

import numpy as np

import athena  # noqa: F401  （保持“依赖 athena”的架构约定）
from xafs_core import XAFSConfigError


# ================================================================
# 1. 公共网格与带约束最小二乘
# ================================================================

def _as_xy(x, y, label: str):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape:
        raise XAFSConfigError(
            "LCF", f"{label} 的 x/y 长度不一致: {x.shape} vs {y.shape}")
    if x.size < 3:
        raise XAFSConfigError("LCF", f"{label} 只有 {x.size} 个点，不足以做 LCF")
    order = np.argsort(x)
    return x[order], y[order]


def _common_grid(target_x, target_y, standards: dict, xmin: float,
                 xmax: float):
    """把目标谱与所有端元插值到公共网格（目标谱在交集内的点）。"""
    if not standards:
        raise XAFSConfigError("LCF", "至少需要一个端元谱")
    tx, ty = _as_xy(target_x, target_y, "target")
    lo, hi = float(xmin), float(xmax)
    lo = max(lo, float(tx.min()))
    hi = min(hi, float(tx.max()))
    prepared = {}
    for name, (sx, sy) in standards.items():
        sx, sy = _as_xy(sx, sy, f"standard {name!r}")
        lo = max(lo, float(sx.min()))
        hi = min(hi, float(sx.max()))
        prepared[name] = (sx, sy)
    if hi <= lo:
        raise XAFSConfigError(
            "LCF", f"公共能量/波数区间为空: [{lo:.3f}, {hi:.3f}]")
    mask = (tx >= lo) & (tx <= hi)
    if int(mask.sum()) < 5:
        raise XAFSConfigError(
            "LCF", f"公共区间 [{lo:.3f}, {hi:.3f}] 内仅 {int(mask.sum())} 个目标点")
    x = tx[mask]
    y = ty[mask]
    names = list(prepared)
    cols = [np.interp(x, prepared[n][0], prepared[n][1]) for n in names]
    return x, y, names, np.column_stack(cols)


def _fit_subset(x, y, S, names, sum_to_one: bool,
                allow_energy_shift: bool):
    """在给定端元子集上做非负（和=1）最小二乘，返回权重/误差/R-factor。"""
    from scipy.optimize import minimize

    n, k = S.shape
    if k == 0:
        raise XAFSConfigError("LCF", "端元子集为空")

    shifts = np.zeros(k, dtype=float)
    if allow_energy_shift:
        def model(params):
            w = params[:k]
            sh = params[k:]
            shifted = np.column_stack([
                np.interp(x - sh[i], x, S[:, i]) for i in range(k)])
            return shifted @ w, shifted
        x0 = np.r_[np.full(k, 1.0 / k), np.zeros(k)]
        bounds = [(0.0, 1.0)] * k + [(-2.0, 2.0)] * k
        cons = ([{"type": "eq", "fun": lambda p: float(np.sum(p[:k]) - 1.0)}]
                if sum_to_one else [])
        res = minimize(lambda p: float(np.sum((y - model(p)[0]) ** 2)),
                       x0, method="SLSQP", bounds=bounds, constraints=cons)
        w = np.asarray(res.x[:k], dtype=float)
        shifts = np.asarray(res.x[k:], dtype=float)
        y_fit, S_fit = model(res.x)
    else:
        if sum_to_one and k == 1:
            w = np.array([1.0])
        elif sum_to_one:
            def obj(wp):
                w = np.r_[wp, 1.0 - float(np.sum(wp))]
                return float(np.sum((y - S @ w) ** 2))
            x0 = np.full(k - 1, 1.0 / k)
            bounds = [(0.0, 1.0)] * (k - 1)
            cons = [{"type": "ineq",
                     "fun": lambda wp: 1.0 - float(np.sum(wp))}]
            res = minimize(obj, x0, method="SLSQP", bounds=bounds,
                           constraints=cons)
            w = np.r_[np.asarray(res.x, dtype=float),
                      1.0 - float(np.sum(res.x))]
        else:
            def obj(w):
                return float(np.sum((y - S @ w) ** 2))
            res = minimize(obj, np.full(k, 1.0 / k), method="SLSQP",
                           bounds=[(0.0, 1.0)] * k)
            w = np.asarray(res.x, dtype=float)
        y_fit = S @ w
        S_fit = S

    # ---- 形式误差：σ²·(JᵀJ)⁻¹；和=1 时用降维参数化传播 ----
    resid = y - y_fit
    dof = max(n - k, 1)
    sigma2 = float(np.sum(resid ** 2) / dof)
    if sum_to_one and k > 1:
        Jp = S_fit[:, :-1] - S_fit[:, -1:]
        cov_p = sigma2 * np.linalg.pinv(Jp.T @ Jp)
        A = np.vstack([np.eye(k - 1), -np.ones((1, k - 1))])
        cov = A @ cov_p @ A.T
    else:
        cov = sigma2 * np.linalg.pinv(S_fit.T @ S_fit)
    w_err = np.sqrt(np.maximum(np.diag(cov), 0.0))

    denom = float(np.sum(y ** 2))
    r_factor = float(np.sum(resid ** 2) / denom) if denom > 0 else np.inf
    return {
        "names": list(names),
        "weights": {n: float(v) for n, v in zip(names, w)},
        "weight_err": {n: float(v) for n, v in zip(names, w_err)},
        "weight_vector": w,
        "r_factor": r_factor,
        "x": np.asarray(x, dtype=float),
        "y": np.asarray(y, dtype=float),
        "y_fit": np.asarray(y_fit, dtype=float),
        "shift": {n: float(v) for n, v in zip(names, shifts)},
        "significant": bool(all(w[i] > 2.0 * w_err[i] for i in range(k))),
        "n_points": int(n),
    }


# ================================================================
# 2. LCF 主入口
# ================================================================

def run_lcf(target_x, target_y, standards: dict, xmin: float, xmax: float,
            sum_to_one: bool = True, allow_energy_shift: bool = False) -> dict:
    """线性组合拟合。

    参数:
        target_x/target_y : 目标谱（XANES 用相对能量 eV；EXAFS 用 k Å⁻¹）
        standards         : {名称: (x, y)}；端元谱
        xmin/xmax         : 拟合区间
        sum_to_one        : 权重非负且和=1（默认）
        allow_energy_shift: 允许每个端元在 ±2 eV 内平移（默认关闭）

    返回:
        ``best``（最优组合）与 ``combinations``（组分>3 时枚举 2/3 组分）。
        顶层同时给出 ``weights``/``weight_err``/``r_factor`` 便于报告层直接读取。
    """
    x, y, names, S = _common_grid(target_x, target_y, standards,
                                  xmin, xmax)
    k_total = len(names)
    combinations = []
    if k_total <= 3:
        best = _fit_subset(x, y, S, names, sum_to_one,
                           allow_energy_shift)
        combinations = [best]
    else:
        for size in (2, 3):
            if size > k_total:
                continue
            for idx in itertools.combinations(range(k_total), size):
                sub_names = [names[i] for i in idx]
                sub_S = S[:, list(idx)]
                fit = _fit_subset(x, y, sub_S, sub_names, sum_to_one,
                                  allow_energy_shift)
                combinations.append(fit)
        if not combinations:
            raise XAFSConfigError("LCF", "没有可拟合的 2/3 组分组合")
        # 先按“所有组分权重显著”再按 R-factor 排序
        combinations.sort(key=lambda f: (not f["significant"],
                                         f["r_factor"]))
        best = combinations[0]

    return {
        "best": best,
        "combinations": combinations,
        "weights": best["weights"],
        "weight_err": best["weight_err"],
        "weight_vector": best["weight_vector"],
        "r_factor": best["r_factor"],
        "x": best["x"],
        "y": best["y"],
        "y_fit": best["y_fit"],
        "names": best["names"],
        "sum_to_one": bool(sum_to_one),
        "allow_energy_shift": bool(allow_energy_shift),
        "x_range": (float(best["x"][0]), float(best["x"][-1])),
    }


# ================================================================
# 3. PCA / 主成分数判定
# ================================================================

def run_pca(spectra: list, xmin: float, xmax: float) -> dict:
    """对一组谱做 SVD，给出方差解释率、Malinowski IND 与 SPOIL。

    ``spectra``: [(名称, x, y), ...]。所有谱插值到第一列谱在公共区间内
    的能量/波数网格；先扣列均值再 SVD（PCA 标准做法）。
    """
    if len(spectra) < 2:
        raise XAFSConfigError("PCA", "至少需要 2 条谱")
    xs = []
    prepared = []
    for name, sx, sy in spectra:
        sx, sy = _as_xy(sx, sy, f"spectrum {name!r}")
        xs.append(sx)
        prepared.append((str(name), sx, sy))
    lo = max(float(xmin), max(float(sx.min()) for sx in xs))
    hi = min(float(xmax), min(float(sx.max()) for sx in xs))
    if hi <= lo:
        raise XAFSConfigError(
            "PCA", f"公共区间为空: [{lo:.3f}, {hi:.3f}]")
    grid = prepared[0][1]
    grid = grid[(grid >= lo) & (grid <= hi)]
    if grid.size < 5:
        raise XAFSConfigError("PCA", "公共网格点数不足 5")
    Y = np.column_stack([np.interp(grid, sx, sy)
                         for _, sx, sy in prepared])
    mean = Y.mean(axis=1, keepdims=True)
    Yc = Y - mean
    U, s, Vt = np.linalg.svd(Yc, full_matrices=False)
    evals = s ** 2
    total = float(np.sum(evals))
    evr = (evals / total) if total > 0 else np.zeros_like(evals)
    cumulative = np.cumsum(evr)
    n_comp = int(np.searchsorted(cumulative, 0.99) + 1)
    n_comp = max(1, min(n_comp, len(s)))

    # Malinowski indicator：IND_k = sqrt(Σ_{i>k} λ_i / ((n−k)(m−k)))
    n_samples = Y.shape[1]
    n_energy = Y.shape[0]
    ind = []
    for k in range(min(n_samples, n_energy)):
        num = float(np.sum(evals[k:]))
        den = float(max((n_samples - k) * (n_energy - k), 1))
        ind.append(float(np.sqrt(num / den)))
    n_ind = int(np.argmin(ind)) if ind else 0

    # SPOIL：用前 n_comp 个主成分重建，残差范数 / 原谱范数
    k = n_comp
    recon_c = U[:, :k] @ (U[:, :k].T @ Yc)
    recon = recon_c + mean
    spoils = {}
    for i, (name, _, _) in enumerate(prepared):
        denom = float(np.sum(Y[:, i] ** 2))
        spoils[name] = (float(np.sqrt(np.sum((Y[:, i] - recon[:, i]) ** 2)
                                     / denom)) if denom > 0 else float("nan"))

    return {
        "x": grid,
        "eigenvalues": evals.tolist(),
        "explained_variance_ratio": evr.tolist(),
        "cumulative_variance": cumulative.tolist(),
        "n_components": n_comp,
        "ind": ind,
        "n_components_ind": n_ind,
        "spoils": spoils,
        "mean": mean[:, 0].tolist(),
        "components": U.tolist(),
        "scores": (U * s).tolist(),
        "loadings": Vt.tolist(),
        "n_samples": n_samples,
        "n_energy": n_energy,
    }


def _grid_on_pca(pca: dict, x, y):
    """把候选谱插值到 PCA 网格并扣列均值。"""
    grid = np.asarray(pca["x"], dtype=float)
    mean = np.asarray(pca["mean"], dtype=float)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    yg = np.interp(grid, x, y)
    return yg - mean


def target_transform(pca: dict, candidates, n_components: int | None = None,
                     spoil_good: float = 3.0,
                     spoil_poor: float = 6.0) -> dict:
    """B10-2 Malinowski Target Transformation（SPOIL 判据）。

    对候选谱 c 投影到前 k 个主成分并重建:

        p = U_kᵀ c,   r = U_k p
        RET   = sqrt(Σ(c−r)²/(n−k))
        REP   = sqrt(Σr²/n)
        SPOIL = RET/REP

    ``c`` 已扣 PCA 训练集均值。SPOIL 越小越像训练集里的真实端元：
    ``< spoil_good`` good，``< spoil_poor`` fair，``< 2·spoil_poor`` poor，
    否则 unacceptable。REP=0 时 SPOIL=inf 并标记 ``degenerate``。
    """
    grid = np.asarray(pca["x"], dtype=float)
    U = np.asarray(pca["components"], dtype=float)
    n = int(grid.size)
    k = int(n_components or pca.get("n_components") or 1)
    k = max(1, min(k, U.shape[1]))
    Uk = U[:, :k]
    out = {}
    for item in candidates:
        if isinstance(item, dict):
            name = str(item.get("name", f"cand{len(out)}"))
            x, y = item["x"], item["y"]
        else:
            name, x, y = item
        c = _grid_on_pca(pca, x, y)
        p = Uk.T @ c
        r = Uk @ p
        dof = max(n - k, 1)
        ret = float(np.sqrt(np.sum((c - r) ** 2) / dof))
        rep = float(np.sqrt(np.sum(r ** 2) / n))
        degenerate = rep <= 1e-14
        spoil = float("inf") if degenerate else ret / rep
        if degenerate:
            rating = "unacceptable"
        elif spoil <= float(spoil_good):
            rating = "good"
        elif spoil <= float(spoil_poor):
            rating = "fair"
        elif spoil <= 2.0 * float(spoil_poor):
            rating = "poor"
        else:
            rating = "unacceptable"
        out[name] = {"spoils": spoil, "ret": ret, "rep": rep,
                     "rating": rating, "degenerate": bool(degenerate),
                     "n_components": k}
    return out


def mcr_als(spectra: list, n_components: int, max_iter: int = 200,
            tol: float = 1e-6, nonneg: bool = True, closure: bool = False,
            seed: int = 0) -> dict:
    """B10-2 原位序列谱的 MCR-ALS 盲解谱。

    模型 ``Y ≈ S·C``（Y: n_energy×n_samples；S: n_energy×k 组分谱；
    C: k×n_samples 浓度）。PCA/SVD 初始化后交替非负最小二乘:

        C ← argmin_{C≥0} ||Y − S C||
        S ← argmin_{S≥0} ||Y − S C||

    ``closure=True`` 时每轮把每个样品的浓度归一（和为 1）。
    返回 ``concentrations``、``components``、``r_factor``、迭代数与收敛标志。
    ``seed`` 仅为 API 兼容预留（SVD 初始化是确定性的）。
    """
    if len(spectra) < 2:
        raise XAFSConfigError("MCR", "至少需要 2 条谱")
    if int(n_components) < 1:
        raise XAFSConfigError("MCR", f"n_components={n_components} 必须 ≥1")
    xs = []
    prepared = []
    for name, sx, sy in spectra:
        sx, sy = _as_xy(sx, sy, f"spectrum {name!r}")
        xs.append(sx)
        prepared.append((str(name), sx, sy))
    lo = max(float(sx.min()) for sx in xs)
    hi = min(float(sx.max()) for sx in xs)
    if hi <= lo:
        raise XAFSConfigError("MCR", f"公共区间为空: [{lo:.3f}, {hi:.3f}]")
    grid = prepared[0][1]
    grid = grid[(grid >= lo) & (grid <= hi)]
    if grid.size < 5:
        raise XAFSConfigError("MCR", "公共网格点数不足 5")
    Y = np.column_stack([np.interp(grid, sx, sy)
                         for _, sx, sy in prepared])
    n_energy, n_samples = Y.shape
    k = max(1, min(int(n_components), n_energy, n_samples))

    def _nnls_like(A, B):
        """逐列解 min ||A·X − B||；nonneg 时用精确 NNLS。"""
        if nonneg:
            from scipy.optimize import nnls
            X = np.column_stack([nnls(A, B[:, j])[0]
                                 for j in range(B.shape[1])])
            return X
        X, *_ = np.linalg.lstsq(A, B, rcond=None)
        return X

    # ---- 初始化 ----
    if nonneg:
        # SIMPLISMA 风格“最纯行”：行纯度 = 行标准差 / 行均值，
        # 贪心选择 k 个互相尽量远的能量点，用这些行作为浓度初值。
        purity = Y.std(axis=1) / (Y.mean(axis=1) + 1e-12)
        chosen: list[int] = []
        for _ in range(k):
            if chosen:
                dist = np.min(np.abs(np.arange(n_energy)[:, None]
                                     - np.asarray(chosen)[None, :]), axis=1)
                score = purity * dist / max(n_energy - 1, 1)
            else:
                score = purity
            idx = int(np.argmax(score))
            if idx in chosen:
                idx = int(np.argmax(np.where(np.isin(np.arange(n_energy),
                                                     chosen), -np.inf, purity)))
            chosen.append(idx)
        C = np.maximum(Y[chosen, :].copy(), 0.0)
        S = np.maximum(_nnls_like(C.T, Y.T).T, 0.0)
    else:
        U, s, Vt = np.linalg.svd(Y, full_matrices=False)
        S = U[:, :k] * np.sqrt(np.maximum(s[:k], 0.0))[None, :]
        C = (np.sqrt(np.maximum(s[:k], 0.0))[:, None] * Vt[:k, :])

    converged = False
    iterations = 0
    prev_r = None
    for it in range(1, int(max_iter) + 1):
        iterations = it
        C = _nnls_like(S, Y)
        if closure:
            denom = np.sum(C, axis=0, keepdims=True)
            C = np.divide(C, denom, out=np.zeros_like(C), where=denom > 1e-14)
        S = _nnls_like(C.T, Y.T).T
        recon = S @ C
        rss = float(np.sum((Y - recon) ** 2))
        denom_y = float(np.sum(Y ** 2)) or 1.0
        rf = float(np.sqrt(rss / denom_y))
        if prev_r is not None and abs(prev_r - rf) <= float(tol):
            converged = True
            break
        prev_r = rf

    recon = S @ C
    rss = float(np.sum((Y - recon) ** 2))
    denom_y = float(np.sum(Y ** 2)) or 1.0
    # MCR 的 S/C 存在尺度不确定性；统一把每条组分谱归一到大 peak=1，
    # 浓度按同一因子缩放（峰归一化是解谱软件的标准输出约定）。
    scale = np.maximum(np.max(S, axis=0), 1e-14)
    S_out = S / scale[None, :]
    C_out = C * scale[:, None]
    return {
        "x": grid,
        "concentrations": C_out,
        "components": S_out,
        "component_scale": scale,
        "names": [name for name, _, _ in prepared],
        "r_factor": float(np.sqrt(rss / denom_y)),
        "iterations": int(iterations),
        "converged": bool(converged),
        "n_components": int(k),
    }


# ================================================================
# 4. 直接从 NormalizedScan / ChiK 做 LCF
# ================================================================

def _scan_xy(scan, mode: str):
    if mode == "xanes":
        if not hasattr(scan, "energy") or not hasattr(scan, "e0"):
            raise XAFSConfigError("LCF", "XANES 模式需要 NormalizedScan")
        x = np.asarray(scan.energy, dtype=float) - float(scan.e0)
        y = (np.asarray(scan.flattened, dtype=float)
             if getattr(scan, "flattened", None) is not None
             else np.asarray(scan.mu_norm, dtype=float))
        return x, y
    if mode == "exafs":
        if not hasattr(scan, "k") or not hasattr(scan, "chi"):
            raise XAFSConfigError("LCF", "EXAFS 模式需要 ChiK")
        return (np.asarray(scan.k, dtype=float),
                np.asarray(scan.chi, dtype=float))
    raise XAFSConfigError("LCF", f"未知 LCF 模式 {mode!r}（xanes / exafs）")


def lcf_from_scans(target, standards: dict, mode: str = "xanes",
                   xmin: float | None = None, xmax: float | None = None,
                   sum_to_one: bool = True,
                   allow_energy_shift: bool = False) -> dict:
    """从 ``NormalizedScan``（XANES）或 ``ChiK``（EXAFS）直接做 LCF。

    XANES 默认区间为相对 E0 的 [−20, +60] eV；EXAFS 默认 [3,12] Å⁻¹。
    """
    tx, ty = _scan_xy(target, mode)
    std_xy = {name: _scan_xy(scan, mode)
              for name, scan in standards.items()}
    if xmin is None:
        xmin = -20.0 if mode == "xanes" else 3.0
    if xmax is None:
        xmax = 60.0 if mode == "xanes" else 12.0
    return run_lcf(tx, ty, std_xy, float(xmin), float(xmax),
                   sum_to_one=sum_to_one,
                   allow_energy_shift=allow_energy_shift)


def lcf_summary(result: dict) -> dict:
    """把 run_lcf 结果压成可落盘/报告的摘要（不写整条谱）。"""
    best = result.get("best") or {}
    return {
        "weights": best.get("weights"),
        "weight_err": best.get("weight_err"),
        "r_factor": best.get("r_factor"),
        "components": best.get("names"),
        "n_components": len(best.get("names") or []),
        "x_range": result.get("x_range"),
        "sum_to_one": result.get("sum_to_one"),
        "combinations": [
            {"components": c.get("names"),
             "weights": c.get("weights"),
             "weight_err": c.get("weight_err"),
             "r_factor": c.get("r_factor"),
             "significant": c.get("significant")}
            for c in (result.get("combinations") or [])
        ],
    }


# ================================================================
# 4.5 原位（Operando）时间序列一键解谱（B15-1）
# ================================================================

def _operando_item(item, idx: int, mode: str):
    """归一化一条原位谱 → (name, x, y, payload)。

    payload 用于描述符轨迹：PipelineResult / NormalizedScan 保留对象，
    (name, x, y) 纯数组项为 None。
    """
    if isinstance(item, (tuple, list)):
        if len(item) == 2:
            name, payload = item[0], item[1]
            if hasattr(payload, "norm"):
                x, y = _scan_xy(payload.norm, mode)
                return str(name), x, y, payload
            if hasattr(payload, "energy") and hasattr(payload, "e0"):
                x, y = _scan_xy(payload, mode)
                return str(name), x, y, payload
            if (isinstance(payload, (tuple, list))
                    and len(payload) == 2):
                return (str(name), np.asarray(payload[0], dtype=float),
                        np.asarray(payload[1], dtype=float), None)
            raise XAFSConfigError(
                "OPERANDO", f"第 {idx} 项的二元组第二元素无法识别为 "
                            f"NormalizedScan/PipelineResult/(x, y)")
        if len(item) == 3:
            return (str(item[0]), np.asarray(item[1], dtype=float),
                    np.asarray(item[2], dtype=float), None)
    if hasattr(item, "norm"):
        x, y = _scan_xy(item.norm, mode)
        return f"scan_{idx:02d}", x, y, item
    if hasattr(item, "energy") and hasattr(item, "e0"):
        x, y = _scan_xy(item, mode)
        return f"scan_{idx:02d}", x, y, item
    raise XAFSConfigError(
        "OPERANDO", f"第 {idx} 项既不是 (name, NormalizedScan)、"
                    f"(name, x, y)，也不是 PipelineResult/NormalizedScan")


def _operando_descriptor(name: str, payload) -> dict:
    """从 PipelineResult / NormalizedScan 提取价态与结构描述符。"""
    out = {"name": name, "e0": None, "e_half_height": None,
           "white_line_height": None, "ft_r1": None}
    scan = None
    if payload is None:
        return out
    if hasattr(payload, "norm"):
        scan = payload.norm
        out["ft_r1"] = getattr(payload, "r1", None)
        e0_obj = getattr(payload, "e0", None)
        if isinstance(e0_obj, dict):
            out["e0"] = e0_obj.get("e0_sample", e0_obj.get("e0"))
        elif e0_obj is not None:
            out["e0"] = float(e0_obj)
    elif hasattr(payload, "energy") and hasattr(payload, "e0"):
        scan = payload
    if scan is None:
        return out
    if out["e0"] is None:
        out["e0"] = float(scan.e0)
    try:
        feat = athena.extract_xanes_features(scan, edge="K")
        out["e_half_height"] = feat.get("e_half_height")
        out["white_line_height"] = feat.get("white_line_height")
    except Exception:  # noqa: BLE001 - 描述符失败不应中断整条序列
        pass
    return out


def _operando_standard_xy(standards: dict | None, mode: str) -> dict:
    """候选标样 → {name: (x, y)}；支持 NormalizedScan / PipelineResult / 数组。"""
    out = {}
    for name, val in (standards or {}).items():
        if hasattr(val, "norm"):
            out[str(name)] = _scan_xy(val.norm, mode)
        elif hasattr(val, "energy") and hasattr(val, "e0"):
            out[str(name)] = _scan_xy(val, mode)
        elif isinstance(val, (tuple, list)) and len(val) == 2:
            out[str(name)] = (np.asarray(val[0], dtype=float),
                              np.asarray(val[1], dtype=float))
        else:
            raise XAFSConfigError(
                "OPERANDO", f"标样 {name!r} 不是 NormalizedScan/"
                            f"PipelineResult/(x, y)")
    return out


def run_operando_series(
    items: list,
    standards: dict | None = None,
    n_components: int | None = None,
    mode: str = "xanes",
    xmin: float = -20.0,
    xmax: float = 60.0,
    closure: bool = True,
    max_iter: int = 800,
    spoil_good: float = 3.0,
    spoil_poor: float = 6.0,
) -> dict:
    """B15-1 原位时间序列一键解谱（PCA → SPOIL → LCF → MCR-ALS）。

    参数
    ----
    items : 至少 3 条谱。元素可为
        ``(name, NormalizedScan)`` / ``(name, PipelineResult)`` /
        ``(name, x, y)`` / 裸 ``PipelineResult`` / 裸 ``NormalizedScan``。
    standards : {名称: NormalizedScan | PipelineResult | (x, y)}，可选。
        先用 SPOIL 剔除伪标样，再对每条原位谱做 LCF。
    n_components : 主成分数；None 时取
        ``max(2, min(PCA 99% 方差分量数, len(items)-1))``。

    步骤
    ----
    1. PCA（``run_pca``）定物种数与 SPOIL 参考；
    2. Target Transform（``target_transform``）过滤 ``rating ==
       "unacceptable"``（SPOIL > 2·spoil_poor）的候选标样；
    3. 对每条谱跑 ``run_lcf`` 得标样权重矩阵；
    4. ``mcr_als`` 盲解谱，并把纯组分按特征能量（窗口内主峰）升序排序，
       消除 SVD/NNLS 的组分排列随机性；
    5. 从 PipelineResult/NormalizedScan 提取 E0、半高能、白线高度与
       FT 第一峰位，构成描述符轨迹。

    返回
    ----
    ``{"names", "n_components", "pca", "target_transform",
       "accepted_standards", "rejected_standards", "lcf_series",
       "lcf_weights", "mcr_als", "descriptors", "warnings", "mode"}``
    """
    if items is None or len(items) < 3:
        raise XAFSConfigError(
            "OPERANDO", f"原位序列至少需要 3 条谱，收到 {0 if items is None else len(items)}")
    if mode not in ("xanes", "exafs"):
        raise XAFSConfigError("OPERANDO", f"未知 mode={mode!r}")
    prepared = [_operando_item(item, i, mode) for i, item in enumerate(items)]
    spectra = [(name, x, y) for name, x, y, _ in prepared]
    names = [name for name, _, _, _ in prepared]

    pca = run_pca(spectra, xmin, xmax)
    if n_components is None:
        k = max(2, min(int(pca.get("n_components") or 2), len(items) - 1))
    else:
        k = int(n_components)
    if k < 1 or k > len(items) - 1:
        raise XAFSConfigError(
            "OPERANDO", f"n_components={k} 必须在 1–{len(items) - 1}")
    k = min(k, len(pca["components"][0]))

    std_xy = _operando_standard_xy(standards, mode)
    warnings: list[str] = []
    tt_out = {}
    accepted: list[str] = []
    rejected: list[str] = []
    if std_xy:
        tt_out = target_transform(
            pca, [(n, xy[0], xy[1]) for n, xy in std_xy.items()],
            n_components=k, spoil_good=spoil_good,
            spoil_poor=spoil_poor)
        accepted = [n for n in std_xy
                    if tt_out[n]["rating"] != "unacceptable"]
        rejected = [n for n in std_xy if n not in accepted]
        if not accepted:
            warnings.append(
                "所有候选标样 SPOIL 均 unacceptable，退回全部标样并标记警告")
            accepted = list(std_xy)
            rejected = []

    lcf_series = []
    lcf_weights = np.zeros((len(accepted), len(spectra)), dtype=float)
    if accepted:
        fit_std = {n: std_xy[n] for n in accepted}
        for j, (name, sx, sy) in enumerate(spectra):
            res = run_lcf(sx, sy, fit_std, xmin, xmax)
            weights = {n: float((res.get("weights") or {}).get(n, 0.0))
                       for n in accepted}
            errs = {n: float((res.get("weight_err") or {}).get(n, 0.0))
                    for n in accepted}
            for i, n in enumerate(accepted):
                lcf_weights[i, j] = weights[n]
            lcf_series.append({"name": name, "weights": weights,
                               "weight_err": errs,
                               "r_factor": float(res.get("r_factor", np.nan))})
    else:
        for name, _, _ in spectra:
            lcf_series.append({"name": name, "weights": {},
                               "weight_err": {}, "r_factor": None})

    mcr_out = mcr_als(spectra, n_components=k, max_iter=max_iter,
                      closure=closure)

    # ---- 按特征能量（窗口内主峰）升序重排 MCR 组分 ----
    xg = np.asarray(mcr_out["x"], dtype=float)
    S = np.asarray(mcr_out["components"], dtype=float)
    C = np.asarray(mcr_out["concentrations"], dtype=float)
    region = (xg >= float(xmin)) & (xg <= float(xmax)) & np.isfinite(xg)
    if not region.any():
        region = np.ones_like(xg, dtype=bool)
    char_e = []
    for j in range(S.shape[1]):
        col = S[:, j]
        mask = region & np.isfinite(col)
        if not mask.any():
            char_e.append(np.inf)
        else:
            char_e.append(float(xg[mask][int(np.argmax(col[mask]))]))
    order = np.argsort(np.asarray(char_e, dtype=float))
    mcr_out = dict(mcr_out)
    mcr_out["components"] = S[:, order]
    mcr_out["concentrations"] = C[order, :]
    mcr_out["component_scale"] = np.asarray(
        mcr_out.get("component_scale", np.ones(S.shape[1])))[order]
    mcr_out["characteristic_energy"] = [float(char_e[i]) for i in order]
    mcr_out["component_order"] = [int(i) for i in order]

    descriptors = [_operando_descriptor(name, payload)
                   for name, _, _, payload in prepared]

    return {
        "names": names,
        "n_components": int(k),
        "pca": pca,
        "target_transform": tt_out,
        "accepted_standards": accepted,
        "rejected_standards": rejected,
        "lcf_series": lcf_series,
        "lcf_weights": lcf_weights,
        "mcr_als": mcr_out,
        "descriptors": descriptors,
        "warnings": warnings,
        "mode": mode,
        "x_range": (float(xmin), float(xmax)),
    }


# ================================================================
# 5. 自测
# ================================================================

if __name__ == "__main__":
    # 三条纯端元 + 已知权重，验证回收
    rng = np.random.default_rng(0)
    x = np.linspace(-20, 60, 401)
    s1 = 1.0 / (1.0 + np.exp(-(x + 5) / 2.0))
    s2 = 1.0 / (1.0 + np.exp(-(x - 5) / 2.0))
    s3 = np.exp(-((x - 15) / 8.0) ** 2)
    y = 0.35 * s1 + 0.65 * s2 + rng.normal(0, 1e-3, x.size)
    res = run_lcf(x, y, {"a": (x, s1), "b": (x, s2), "c": (x, s3)},
                  -20, 60)
    print("LCF weights:", res["weights"], "Rf=", res["r_factor"])
    pca = run_pca([("m1", x, 0.2 * s1 + 0.8 * s2),
                   ("m2", x, 0.5 * s1 + 0.5 * s2),
                   ("m3", x, 0.8 * s1 + 0.2 * s2),
                   ("m4", x, 0.5 * s1 + 0.5 * s3),
                   ("m5", x, 0.5 * s2 + 0.5 * s3),
                   ("m6", x, 0.3 * s1 + 0.3 * s2 + 0.4 * s3)], -20, 60)
    print("PCA n_comp:", pca["n_components"],
          "cum3:", float(np.sum(pca["explained_variance_ratio"][:3])))
