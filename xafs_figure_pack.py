# -*- coding: utf-8 -*-
"""xafs_figure_pack.py — 催化剂 XAFS 结果的出版级图包
=====================================================

本模块**只读取已有输出**，不重新拟合、不改任何物理判定。它把一个
``run_xafs.py`` / ``xafs_report.py`` 输出目录整理成一套可直接用于论文、
组会和审稿回复的图包：

- ``Figure1_master.png``：一个催化剂结构的 8 面板总图
  (a) XANES  (b) χ(k)  (c) FT  (d) 小波
  (e) 拟合壳层  (f) 结构识别  (g) 正式 CN 门禁  (h) 结构卡片
- ``Figure2_confidence.png``：置信与过程仪表盘
  (a) 门禁状态板  (b) 决策轨迹  (c) 参数来源
  (d) 门禁表  (e) 决策表  (f) 关键参数表
- ``Figure3_diagnostics.png``：拟合与误差诊断
  (a) 参数相关矩阵  (b) 路径 N–R 图  (c) CN 误差预算
  (d) Rbkg 寻优曲线  (e) 归一化分解  (f) 拟合窗口表
  (g) 稳定性检查表  (h) 扫描/数据质量表
- ``Figure4_features.png``：数据特征与输入清单
  (a) dμ/dE  (b) q 空间  (c) FT 峰识别
  (d) CIF 配位壳层  (e) 特征总表  (f) 输入文件清单
- ``Figure2g_cn_table.png`` / ``Figure2h_fit_table.png``：表格直接出图
- 每个面板同名 ``.csv``：保存该面板实际绘制的数据
- ``data/``：自包含数据包（源 JSON/CSV 副本、全路径表、全参数表、
  拟合变体、相关矩阵、CIF 壳层、结构候选、门禁表、决策表、
  汇总 JSON、面板索引与 ``README.md``）
- ``figure_manifest.json``：Figure 1(a) ↔ CSV ↔ 源文件 ↔ 图注的映射
- ``figure_pack.md``：人可读的图包索引

用法::

    python xafs_figure_pack.py out/RuO2
    python xafs_figure_pack.py out/RuO2 --stem Ru-K-RuO2 --formats png,pdf
    python xafs_figure_pack.py out/FeNC --dpi 300 --quiet

也可以在原入口直接生成::

    python run_xafs.py sample.xdi --absorber Fe --identify --fit --figure-pack
    python xafs_report.py sample.xdi --absorber Fe --figure-pack

设计边界：
- 不依赖旧 ``plot_five_panel`` / ``xafs_report`` 的绘图逻辑；旧图仍保留。
- 缺失的面板会写成 ``status=skipped``，不会伪造数据。
- CSV 保存的是**画图时实际使用的列**，不是重新计算的物理量。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle


# ---------------------------------------------------------------
# 基础样式
# ---------------------------------------------------------------
STATUS_COLORS = {
    "passed": "#2e7d32",
    "ok": "#2e7d32",
    "warning": "#f9a825",
    "warned": "#f9a825",
    "failed": "#c62828",
    "not_applicable": "#9e9e9e",
    "not_available": "#9e9e9e",
    "skipped": "#bdbdbd",
    "info": "#1565c0",
    "diagnostic": "#1565c0",
    "diagnostic_only": "#1565c0",
    "consistency_failed": "#c62828",
    "multi_site_failed": "#c62828",
    "no_cif": "#c62828",
    "legacy_unverified": "#f9a825",
    "hypothesis_model_only": "#9e9e9e",
}
SOURCE_COLORS = {
    "cli": "#1565c0",
    "config": "#2e7d32",
    "default": "#90a4ae",
    "auto": "#ef6c00",
    "calibrated": "#6a1b9a",
    "runtime_override": "#c62828",
}
BAND_COLORS = {
    "reject": "#ffcdd2",
    "warn": "#fff9c4",
    "ok": "#c8e6c9",
}
GATE_ORDER = (
    "normalization",
    "self_absorption",
    "structure_id",
    "cif",
    "fit",
    "formal_cn",
    "xanes_cross_check",
)
GATE_STATUSES = ("passed", "warning", "failed", "not_applicable", "not_available")


class FigurePackError(RuntimeError):
    """图包无法生成时抛出。"""


def _style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        # Microsoft YaHei 覆盖中文；数学上标用 mathtext 渲染（该字体缺
        # U+207B），避免表格中的中文说明变成方框。
        "font.sans-serif": ["Microsoft YaHei", "DejaVu Sans", "Arial"],
        "mathtext.fontset": "dejavusans",
        "axes.linewidth": 0.9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })


# ---------------------------------------------------------------
# 通用读取
# ---------------------------------------------------------------
def _read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        raise FigurePackError(f"JSON 读取失败: {path}: {exc}") from exc


def _read_csv_rows(path: Path) -> list[dict]:
    try:
        with path.open("r", encoding="utf-8-sig", errors="replace",
                       newline="") as fh:
            return list(csv.DictReader(fh))
    except Exception as exc:  # noqa: BLE001
        raise FigurePackError(f"CSV 读取失败: {path}: {exc}") from exc


def _read_numeric(path: Path) -> tuple[list[str] | None, np.ndarray]:
    """读取 2D 数值表；CSV 首行若含非数字则当表头。"""
    if path.suffix.lower() == ".csv":
        with path.open("r", encoding="utf-8-sig", errors="replace",
                       newline="") as fh:
            rows = [r for r in csv.reader(fh) if r]
        if not rows:
            return None, np.empty((0, 0))
        header: list[str] | None = None
        try:
            float(rows[0][0])
        except (TypeError, ValueError):
            header = [str(x).strip() for x in rows[0]]
            rows = rows[1:]
        data = []
        for row in rows:
            try:
                data.append([float(x) for x in row])
            except (TypeError, ValueError):
                continue
        return header, np.asarray(data, dtype=float)
    try:
        arr = np.loadtxt(path, comments="#", ndmin=2)
    except Exception as exc:  # noqa: BLE001
        raise FigurePackError(f"数值表读取失败: {path}: {exc}") from exc
    return None, np.asarray(arr, dtype=float)


def _read_wavelet(path: Path) -> dict | None:
    """读取 ``_wavelet.csv`` / ``_wt.dat`` → k, R, z。"""
    try:
        if path.suffix.lower() == ".csv":
            with path.open("r", encoding="utf-8-sig", errors="replace",
                           newline="") as fh:
                rows = [r for r in csv.reader(fh) if r]
            if len(rows) < 2:
                return None
            header = rows[0]
            r = np.asarray([float(x) for x in header[1:]], dtype=float)
            ks, zs = [], []
            for row in rows[1:]:
                if not row:
                    continue
                ks.append(float(row[0]))
                zs.append([float(x) for x in row[1:]])
            if not ks:
                return None
            k = np.asarray(ks, dtype=float)
            z = np.asarray(zs, dtype=float)
        else:
            lines = [ln.strip() for ln in
                     path.read_text(encoding="utf-8", errors="replace").splitlines()
                     if ln.strip()]
            # 头部形如 "# k \ R  0.00  0.02 ..."：先找带 R 的注释行，
            # 不能先把所有 '#' 行都删掉，否则第一行数据会被当成表头。
            header_idx = None
            header_parts: list[str] = []
            for idx, ln in enumerate(lines):
                stripped = ln.lstrip("# ").strip()
                parts = stripped.split()
                if not parts:
                    continue
                if parts[0].lower() == "k" and any(
                        ("\\" in p) or p.lower() == "r" for p in parts[1:3]):
                    header_idx = idx
                    header_parts = parts
                    break
            if header_idx is None:
                body = [ln for ln in lines if not ln.startswith("#")]
                if len(body) < 2:
                    return None
                header_parts = body[0].split()
                header_idx = lines.index(body[0]) - 1
                data_lines = body[1:]
            else:
                data_lines = [ln for ln in lines[header_idx + 1:]
                              if not ln.startswith("#")]
            r_vals = []
            for token in header_parts:
                try:
                    r_vals.append(float(token))
                except ValueError:
                    continue
            if not r_vals:
                return None
            r = np.asarray(r_vals, dtype=float)
            ks, zs = [], []
            for ln in data_lines:
                parts = ln.split()
                if len(parts) < 2:
                    continue
                ks.append(float(parts[0]))
                zs.append([float(x) for x in parts[1:]])
            if not ks:
                return None
            k = np.asarray(ks, dtype=float)
            z = np.asarray(zs, dtype=float)
        ncol = min(len(r), z.shape[1] if z.ndim == 2 else 0)
        if ncol == 0:
            return None
        return {"k": k, "r": r[:ncol], "z": z[:, :ncol], "source": path}
    except Exception:  # noqa: BLE001 - 单个面板失败不应中断整包
        return None


def _first_file(outdir: Path, stem: str, suffixes: tuple[str, ...]) -> Path | None:
    for suffix in suffixes:
        path = outdir / f"{stem}{suffix}"
        if path.is_file():
            return path
    return None


def discover_stem(outdir: Path, stem: str | None = None) -> str:
    """从输出目录自动识别样品 stem。"""
    outdir = Path(outdir)
    if stem:
        return str(stem)
    metas = [p for p in outdir.glob("*_meta.json")
             if p.name != "operando_meta.json"]
    if metas:
        preferred = outdir / f"{outdir.name}_meta.json"
        if preferred.is_file():
            return preferred.name[:-len("_meta.json")]
        metas.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        return metas[0].name[:-len("_meta.json")]
    for pattern, tail in (("*_chik.csv", "_chik.csv"),
                          ("*_chi.dat", "_chi.dat"),
                          ("*_xanes.csv", "_xanes.csv"),
                          ("*_mu.csv", "_mu.csv"),
                          ("*_mu.dat", "_mu.dat")):
        found = sorted(outdir.glob(pattern))
        if len(found) == 1:
            return found[0].name[:-len(tail)]
    raise FigurePackError(
        f"无法在 {outdir} 识别样品 stem；请显式传 --stem。"
        " 需要 *_meta.json / *_chik.csv / *_xanes.csv 之一。")


def _xanes_from(ctx: dict) -> dict | None:
    outdir, stem = ctx["outdir"], ctx["stem"]
    path = _first_file(outdir, stem, ("_xanes.csv", "_mu.csv", "_mu.dat"))
    if path is None:
        return None
    header, arr = _read_numeric(path)
    if arr.ndim != 2 or arr.shape[0] < 2:
        return None
    if path.name.endswith("_xanes.csv") and header:
        idx = {str(h).strip(): i for i, h in enumerate(header)}
        x = arr[:, idx.get("energy_eV", 0)]
        for name in ("mu_norm", "mu_flat", "mu"):
            if name in idx:
                y = arr[:, idx[name]]
                break
        else:
            y = arr[:, 1]
    elif arr.shape[1] >= 3:
        x, y = arr[:, 0], arr[:, 2]
    else:
        x, y = arr[:, 0], arr[:, 1]
    return {"energy": x, "mu_norm": y, "source": path}


def _chik_from(ctx: dict) -> dict | None:
    outdir, stem = ctx["outdir"], ctx["stem"]
    path = _first_file(outdir, stem, ("_chik.csv", "_chi.dat"))
    if path is None:
        return None
    header, arr = _read_numeric(path)
    if arr.ndim != 2 or arr.shape[0] < 2:
        return None
    if path.name.endswith("_chik.csv") and header:
        idx = {str(h).strip(): i for i, h in enumerate(header)}
        k = arr[:, idx.get("k_inv_A", 0)]
        chi = arr[:, idx.get("chi", 1)]
        k2chi = (arr[:, idx["k2chi"]] if "k2chi" in idx
                 else k * k * chi)
    else:
        k = arr[:, 0]
        chi = arr[:, 1]
        k2chi = arr[:, 2] if arr.shape[1] > 2 else k * k * chi
    return {"k": k, "chi": chi, "k2chi": k2chi, "source": path}


def _ft_from(ctx: dict) -> dict | None:
    outdir, stem = ctx["outdir"], ctx["stem"]
    path = _first_file(outdir, stem, ("_chir.csv", "_ft.dat"))
    if path is None:
        return None
    header, arr = _read_numeric(path)
    if arr.ndim != 2 or arr.shape[0] < 2:
        return None
    if path.name.endswith("_chir.csv") and header:
        idx = {str(h).strip(): i for i, h in enumerate(header)}
        r = arr[:, idx.get("R_A", 0)]
        re = arr[:, idx.get("chi_R_re", 1)]
        im = arr[:, idx.get("chi_R_im", 2)]
        mag = arr[:, idx["magnitude"]] if "magnitude" in idx else np.hypot(re, im)
    else:
        r = arr[:, 0]
        if arr.shape[1] >= 4:
            re, im, mag = arr[:, 1], arr[:, 2], arr[:, 3]
        elif arr.shape[1] == 3:
            re, im = arr[:, 1], arr[:, 2]
            mag = np.hypot(re, im)
        else:
            re = np.full_like(r, np.nan)
            im = np.full_like(r, np.nan)
            mag = arr[:, 1]
    return {"r": r, "re": re, "im": im, "mag": mag, "source": path}


def _chiq_from(ctx: dict) -> dict | None:
    outdir, stem = ctx["outdir"], ctx["stem"]
    path = _first_file(outdir, stem, ("_chiq.csv", "_chiq.dat"))
    if path is None:
        return None
    header, arr = _read_numeric(path)
    if arr.ndim != 2 or arr.shape[0] < 2:
        return None
    if header:
        idx = {str(h).strip(): i for i, h in enumerate(header)}
        q = arr[:, idx.get("q_inv_A", 0)]
        re = arr[:, idx["chiq_re"]] if "chiq_re" in idx else arr[:, 1]
        im = arr[:, idx["chiq_im"]] if "chiq_im" in idx else arr[:, 2]
        mag = (arr[:, idx["chiq_mag"]] if "chiq_mag" in idx
               else np.hypot(re, im))
    else:
        q = arr[:, 0]
        if arr.shape[1] >= 4:
            re, im, mag = arr[:, 1], arr[:, 2], arr[:, 3]
        else:
            re = arr[:, 1]
            im = arr[:, 2] if arr.shape[1] > 2 else np.zeros_like(q)
            mag = np.hypot(re, im)
    return {"q": q, "re": re, "im": im, "mag": mag, "source": path}


def _as_float(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _dig(payload: Any, *keys: str) -> Any:
    """安全嵌套取值：任一层缺失或为 None 时返回 None。"""
    cur = payload
    for key in keys:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(key)
    return cur


def _meta_e0(meta: dict) -> float | None:
    e0 = meta.get("e0")
    if isinstance(e0, dict):
        return _as_float(e0.get("e0_sample") or e0.get("e0_derivative_peak")
                         or e0.get("e0_nominal"))
    return _as_float(e0)


def _meta_absorber(meta: dict) -> str:
    if meta.get("absorber"):
        return str(meta["absorber"])
    params = (meta.get("params") or {}).get("effective") or {}
    return str(params.get("absorber") or "")


def _meta_edge(meta: dict, absorber: str) -> str:
    edge = meta.get("edge")
    if edge:
        return str(edge)
    return f"{absorber}-edge" if absorber else "XAFS"


def _sample_name(ctx: dict) -> str:
    meta = ctx.get("meta") or {}
    src = meta.get("source") or {}
    if isinstance(src, dict) and src.get("path"):
        return Path(str(src["path"])).stem
    return str(ctx["stem"])


def load_context(outdir: Path, stem: str | None = None) -> dict:
    """读取一个输出目录中所有可用的图包输入。"""
    outdir = Path(outdir)
    if not outdir.is_dir():
        raise FigurePackError(f"输出目录不存在: {outdir}")
    stem = discover_stem(outdir, stem)
    meta_path = outdir / f"{stem}_meta.json"
    meta = _read_json(meta_path) if meta_path.is_file() else {}
    fit_path = outdir / f"{stem}_fit.json"
    fit = _read_json(fit_path) if fit_path.is_file() else None
    sid_path = outdir / f"{stem}_structure_id.json"
    structure = _read_json(sid_path) if sid_path.is_file() else None
    if structure is None and meta.get("structure_id"):
        structure = dict(meta["structure_id"])
    if structure is not None:
        structure = dict(structure)
        if not structure.get("ranked") and structure.get("top3"):
            structure["ranked"] = list(structure["top3"])
    cn_path = (outdir / "cn_table.csv"
               if (outdir / "cn_table.csv").is_file()
               else outdir / f"{stem}_cn_table.csv")
    cn_rows = _read_csv_rows(cn_path) if cn_path.is_file() else []
    fit_table_path = _first_file(outdir, stem, ("_fit_table.csv", "_fit.csv"))
    fit_table = _read_csv_rows(fit_table_path) if fit_table_path else []
    ctx = {
        "outdir": outdir,
        "stem": stem,
        "meta": meta,
        "meta_path": meta_path if meta_path.is_file() else None,
        "fit": fit,
        "fit_path": fit_path if fit_path.is_file() else None,
        "structure": structure,
        "structure_path": sid_path if sid_path.is_file() else None,
        "cn_rows": cn_rows,
        "cn_path": cn_path if cn_path.is_file() else None,
        "fit_table": fit_table,
        "fit_table_path": fit_table_path,
        "xanes": _xanes_from({"outdir": outdir, "stem": stem}),
        "chik": _chik_from({"outdir": outdir, "stem": stem}),
        "ft": _ft_from({"outdir": outdir, "stem": stem}),
        "chiq": _chiq_from({"outdir": outdir, "stem": stem}),
        "wt": None,
    }
    wt_path = _first_file(outdir, stem, ("_wavelet.csv", "_wt.dat"))
    if wt_path is not None:
        ctx["wt"] = _read_wavelet(wt_path)
    if not ctx["meta"] and not any(ctx[k] for k in ("xanes", "chik", "ft")):
        raise FigurePackError(
            f"{outdir} 中没有可识别的 XAFS 输出；至少需要 "
            "*_meta.json / *_xanes.csv / *_chik.csv / *_ft.csv 之一。")
    return ctx


# ---------------------------------------------------------------
# 格式化 / 表格
# ---------------------------------------------------------------
def _fmt(value: Any, nd: int = 4) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "—"
        return f"{value:.{nd}g}"
    return str(value)


def _fmt_pm(value: Any, err: Any, nd: int = 3) -> str:
    v = _as_float(value)
    if v is None:
        return "—"
    e = _as_float(err)
    return f"{v:.{nd}f}±{e:.{nd}f}" if e is not None else f"{v:.{nd}f}"


def _short(value: Any, limit: int = 64) -> str:
    if value is None:
        return "—"
    if isinstance(value, (dict, list, tuple)):
        try:
            text = json.dumps(value, ensure_ascii=False, separators=(",", ":"))
        except Exception:  # noqa: BLE001
            text = str(value)
    else:
        text = str(value)
    text = text.replace("\n", " ")
    return text if len(text) <= limit else text[:limit - 1] + "…"


def _status_color(status: Any) -> str:
    return STATUS_COLORS.get(str(status).lower(), "#9e9e9e")


def _display(value: Any) -> str:
    """表格/图中显示用：把字体缺字的上/下标转成纯 ASCII。"""
    text = str(value)
    for src, dst in (("Å⁻¹", "Å^-1"), ("Å⁻²", "Å^-2"),
                     ("Å⁻³", "Å^-3"), ("⁻¹", "^-1"),
                     ("⁻²", "^-2"), ("⁻³", "^-3")):
        text = text.replace(src, dst)
    for i in range(10):
        text = text.replace(chr(0x2080 + i), str(i))
    return text


def _draw_missing(ax, reason: str) -> None:
    ax.axis("off")
    ax.text(0.5, 0.55, "Not available", ha="center", va="center",
            fontsize=11, color="#616161", transform=ax.transAxes)
    ax.text(0.5, 0.42, reason, ha="center", va="center", fontsize=7.5,
            color="#757575", wrap=True, transform=ax.transAxes)


def _render_table(ax, headers: list[str], rows: list[list[str]],
                  *, status_col: int | None = None, fontsize: float = 7.0,
                  title: str | None = None) -> None:
    ax.axis("off")
    if title:
        ax.set_title(_display(title), fontsize=9.5, pad=6)
    if not rows:
        ax.text(0.5, 0.5, "No rows", ha="center", va="center",
                fontsize=9, color="#757575", transform=ax.transAxes)
        return
    headers = [_display(h) for h in headers]
    rows = [[_display(c) for c in row] for row in rows]
    ncol = len(headers)
    widths = []
    for i, head in enumerate(headers):
        lens = [len(str(head))] + [len(str(r[i])) for r in rows[:20]]
        widths.append(max(8, min(46, max(lens))))
    total = float(sum(widths))
    col_widths = [w / total for w in widths]
    table = ax.table(
        cellText=rows, colLabels=headers, loc="upper center",
        colWidths=col_widths, cellLoc="left", colLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)
    table.scale(1.0, 1.18)
    for (row, col), cell in table.get_celld().items():
        cell.set_linewidth(0.4)
        cell.set_edgecolor("#cfd8dc")
        if row == 0:
            cell.set_facecolor("#eceff1")
            cell.set_text_props(weight="bold")
        elif status_col is not None and col == status_col:
            status = rows[row - 1][col]
            color = _status_color(status)
            cell.set_facecolor(color)
            cell.get_text().set_color("white" if status in (
                "passed", "failed", "ok") else "#212121")


def _summary_counts(gates: dict) -> dict:
    summary = gates.get("summary") or {}
    return {
        "passed": list(summary.get("passed") or []),
        "warned": list(summary.get("warned") or []),
        "failed": list(summary.get("failed") or []),
        "not_applicable": list(summary.get("not_applicable") or []),
    }


def _formal_block(ctx: dict) -> dict:
    meta = ctx.get("meta") or {}
    fit = ctx.get("fit") or {}
    return meta.get("formal_cn") or fit.get("formal_cn") or {}


def _fit_rows(ctx: dict) -> list[dict]:
    """从 cn_table / fit.json / structure_id / meta 提取拟合壳层。"""
    rows: list[dict] = []
    cn_rows = ctx.get("cn_rows") or []
    cn_by_shell: dict[str, dict] = {}
    for row in cn_rows:
        n_fit = _as_float(row.get("n_fit"))
        if n_fit is None or n_fit <= 0:
            continue
        shell = (row.get("shell_label") or row.get("shell_group")
                 or row.get("absorber_site") or "shell")
        r_cif = _as_float(row.get("r_cif"))
        delta_r = _as_float(row.get("delta_r"))
        r = (r_cif + delta_r) if (r_cif is not None and delta_r is not None) else r_cif
        item = {
            "shell": str(shell),
            "n": n_fit,
            "n_err": _as_float(row.get("n_err_formal")),
            "n0": _as_float(row.get("n0")),
            "n_ratio": _as_float(row.get("n_ratio")),
            "r": r,
            "r_err": None,
            "sigma2": _as_float(row.get("sigma2")),
            "rf": _as_float(row.get("rf")),
            "mode": "formal_cn",
            "status": row.get("gate_status") or "unknown",
            "note": row.get("reject_reason") or "",
            "source": "cn_table.csv",
        }
        rows.append(item)
        cn_by_shell[item["shell"]] = item

    fit = ctx.get("fit") or {}
    fit_obj = fit.get("fit") if isinstance(fit.get("fit"), dict) else None
    shell_degen: dict[str, float | None] = {}
    for shell in (fit.get("shells") or []):
        if isinstance(shell, dict) and shell.get("label"):
            shell_degen[str(shell["label"])] = _as_float(shell.get("degen"))
    if fit_obj:
        for path in (fit_obj.get("paths") or []):
            if not isinstance(path, dict):
                continue
            label = str(path.get("label") or "shell")
            if label in cn_by_shell:
                # 已有正式 CN 行时，用 fit 的误差/R 补强
                item = cn_by_shell[label]
                item["r"] = item["r"] if item["r"] is not None else _as_float(path.get("R"))
                item["r_err"] = _as_float(path.get("R_err"))
                item["sigma2"] = (item["sigma2"] if item["sigma2"] is not None
                                  else _as_float(path.get("sigma2")))
                continue
            n = _as_float(path.get("N"))
            if n is None:
                continue
            rows.append({
                "shell": label,
                "n": n,
                "n_err": _as_float(path.get("N_err")),
                "n0": shell_degen.get(label),
                "n_ratio": None,
                "r": _as_float(path.get("R")),
                "r_err": _as_float(path.get("R_err")),
                "sigma2": _as_float(path.get("sigma2")),
                "rf": _as_float(fit_obj.get("r_factor")),
                "mode": str(fit.get("mode") or "fit"),
                "status": fit.get("cn_status") or "diagnostic",
                "note": fit.get("warning") or "",
                "source": "fit.json",
            })
    if rows:
        return rows

    l1 = fit.get("l1") or _dig(ctx.get("structure"), "l1_first_shell")
    if not l1:
        l1 = _dig(ctx.get("meta"), "gates", "fit", "observed", "l1")
    if l1:
        n = _as_float(l1.get("N"))
        if n is not None:
            n0 = None
            cif = (_dig(ctx.get("meta"), "gates", "cif", "observed", "sites")
                   or [{}])
            if cif:
                n0 = _as_float(cif[0].get("first_shell_cn"))
            ligand = l1.get("ligand") or "X"
            geom = l1.get("geometry") or ""
            rows.append({
                "shell": f"{ligand} ({geom})" if geom else str(ligand),
                "n": n,
                "n_err": _as_float(l1.get("N_err")),
                "n0": n0,
                "n_ratio": (n / n0) if (n0 not in (None, 0)) else None,
                "r": _as_float(l1.get("R")),
                "r_err": _as_float(l1.get("R_err")),
                "sigma2": _as_float(l1.get("sigma2")),
                "rf": _as_float(l1.get("rf")),
                "mode": str(fit.get("mode") or "l1_first_shell"),
                "status": fit.get("cn_status") or "diagnostic",
                "note": l1.get("warning") or "",
                "source": "fit.json:l1",
            })
    return rows


# ---------------------------------------------------------------
# 拟合数据提取
# ---------------------------------------------------------------
def _fit_obj(ctx: dict) -> dict:
    obj = (ctx.get("fit") or {}).get("fit")
    return obj if isinstance(obj, dict) else {}


def _fit_shell_degen(ctx: dict) -> dict[str, float | None]:
    out: dict[str, float | None] = {}
    for shell in ((ctx.get("fit") or {}).get("shells") or []):
        if isinstance(shell, dict) and shell.get("label"):
            out[str(shell["label"])] = _as_float(shell.get("degen"))
    return out


def _fit_paths(ctx: dict) -> list[dict]:
    """``fit.json`` 主拟合的路径参数。"""
    obj = _fit_obj(ctx)
    degen = _fit_shell_degen(ctx)
    rows = []
    for path in (obj.get("paths") or []):
        if not isinstance(path, dict):
            continue
        label = str(path.get("label") or "path")
        n = _as_float(path.get("N"))
        r = _as_float(path.get("R"))
        reff = _as_float(path.get("reff"))
        n0 = degen.get(label)
        rows.append({
            "variant": "fit",
            "label": label,
            "n": n,
            "n_err": _as_float(path.get("N_err")),
            "r": r,
            "r_err": _as_float(path.get("R_err")),
            "reff": reff,
            "delta_r": ((r - reff) if (r is not None and reff is not None)
                        else _as_float(path.get("delta_r"))),
            "sigma2": _as_float(path.get("sigma2")),
            "sigma2_err": _as_float(path.get("sigma2_err")),
            "n0": n0,
            "n_ratio": ((n / n0) if (n is not None and n0 not in (None, 0))
                        else None),
            "rf": _as_float(obj.get("r_factor")),
            "note": "",
        })
    if rows:
        return rows
    l1 = (ctx.get("fit") or {}).get("l1") or _dig(
        ctx.get("structure"), "l1_first_shell")
    if l1:
        n = _as_float(l1.get("N"))
        r = _as_float(l1.get("R"))
        reff = _as_float(l1.get("r0"))
        rows.append({
            "variant": "l1_first_shell",
            "label": f"{l1.get('ligand') or 'X'}-{l1.get('geometry') or 'shell'}",
            "n": n,
            "n_err": _as_float(l1.get("N_err")),
            "r": r,
            "r_err": _as_float(l1.get("R_err")),
            "reff": reff,
            "delta_r": ((r - reff) if (r is not None and reff is not None)
                        else None),
            "sigma2": _as_float(l1.get("sigma2")),
            "sigma2_err": None,
            "n0": None,
            "n_ratio": None,
            "rf": _as_float(l1.get("rf")),
            "note": l1.get("warning") or "",
        })
    return rows


def _fit_variant_rows(ctx: dict) -> list[dict]:
    """固定 N / 条件 N / 自由 N 三套拟合的路径对照。"""
    rows: list[dict] = []
    fit = ctx.get("fit") or {}
    cond = fit.get("conditional_cn") or {}
    for key, variant in (("fixed_n_paths", "fixed_N"),
                         ("conditional_n_paths", "conditional_N"),
                         ("free_n_paths", "free_N")):
        for path in (cond.get(key) or []):
            if not isinstance(path, dict):
                continue
            reff = _as_float(path.get("reff"))
            dr = _as_float(path.get("delta_r"))
            r = (reff + dr) if (reff is not None and dr is not None) else None
            rows.append({
                "variant": variant,
                "label": path.get("label") or "path",
                "n": _as_float(path.get("N")),
                "n_err": _as_float(path.get("N_err")),
                "r": r,
                "r_err": None,
                "reff": reff,
                "delta_r": dr,
                "sigma2": _as_float(path.get("sigma2")),
                "sigma2_err": None,
                "n0": None,
                "n_ratio": None,
                "rf": _as_float((cond.get(key.replace("_paths", "_fit"))
                                 or {}).get("r_factor")),
                "note": "",
            })
    return rows


def _file_rows(path: Path | None) -> int | None:
    if path is None or not path.is_file():
        return None
    if path.suffix.lower() not in (".csv", ".dat", ".txt"):
        return None
    try:
        with path.open("r", encoding="utf-8-sig", errors="replace") as fh:
            n = sum(1 for line in fh if line.strip())
        return max(0, n - 1)
    except OSError:
        return None


def _file_sha(path: Path | None) -> str:
    if path is None or not path.is_file():
        return ""
    h = hashlib.sha256()
    try:
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
    except OSError:
        return ""
    return h.hexdigest()


def _panel_csv_name(panel: Panel) -> str | None:
    return f"{panel.id}.csv" if panel.rows else None


# ---------------------------------------------------------------
# 图 1 面板
# ---------------------------------------------------------------
@dataclass
class Panel:
    id: str
    label: str
    title: str
    caption: str
    source: str = ""
    columns: list[str] = field(default_factory=list)
    rows: list[dict] = field(default_factory=list)
    status: str = "ok"
    reason: str = ""


def panel_xanes(ax, ctx: dict) -> Panel:
    data = ctx.get("xanes")
    if not data:
        _draw_missing(ax, "缺少 *_xanes.csv / *_mu.csv / *_mu.dat")
        return Panel("Figure1a_xanes", "Figure 1(a)", "XANES",
                     "归一化 XANES：边位置与白线反映氧化态和配位对称性。",
                     status="skipped", reason="no XANES/μ(E) data")
    e = data["energy"]
    y = data["mu_norm"]
    e0 = _meta_e0(ctx.get("meta") or {})
    ax.plot(e, y, color="#0d47a1", lw=1.1)
    if e0 is not None:
        ax.axvline(e0, color="#c62828", ls="--", lw=0.9)
        ax.text(e0, ax.get_ylim()[1], f" E0={e0:.1f}", color="#c62828",
                fontsize=7, va="top")
        lo, hi = max(float(np.min(e)), e0 - 80.0), min(float(np.max(e)), e0 + 250.0)
        if hi > lo:
            ax.set_xlim(lo, hi)
        mask = (e >= e0 + 5.0) & (e <= e0 + 100.0)
        if np.any(mask):
            idx = int(np.nanargmax(y[mask]))
            ax.scatter([e[mask][idx]], [y[mask][idx]], s=16, color="#ef6c00",
                       zorder=4)
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel(r"Normalized $\mu(E)$")
    ax.set_title("(a) XANES", loc="left", fontweight="bold")
    ax.grid(alpha=0.18)
    rows = [{"energy_eV": float(x), "mu_norm": float(v)}
            for x, v in zip(e, y) if math.isfinite(float(v))]
    return Panel(
        "Figure1a_xanes", "Figure 1(a)", "Normalized XANES",
        "归一化 XANES：吸收边位置、白线强度和边前特征用于判断活性金属"
        "氧化态、配位对称性与电子结构变化。",
        source=_rel(ctx, data["source"]),
        columns=["energy_eV", "mu_norm"], rows=rows)


def panel_exafs(ax, ctx: dict) -> Panel:
    data = ctx.get("chik")
    if not data:
        _draw_missing(ax, "缺少 *_chik.csv / *_chi.dat")
        return Panel("Figure1b_exafs", "Figure 1(b)", "EXAFS χ(k)",
                     "k²加权 EXAFS 振荡：键长、配位数与散射体信息。",
                     status="skipped", reason="no χ(k) data")
    k, y = data["k"], data["k2chi"]
    ax.plot(k, y, color="#1b5e20", lw=1.0)
    ft_info = (ctx.get("meta") or {}).get("ft") or {}
    params = ((ctx.get("meta") or {}).get("params") or {}).get("effective") or {}
    kmin = _as_float(ft_info.get("kmin", params.get("kmin")))
    kmax = _as_float(ft_info.get("kmax", params.get("kmax")))
    if kmin is not None and kmax is not None and kmax > kmin:
        ax.axvspan(kmin, kmax, color="#c8e6c9", alpha=0.45, lw=0)
        ax.text(kmin, ax.get_ylim()[1],
                f" k={kmin:g}–{kmax:g} Å$^{{-1}}$",
                fontsize=7, va="top", color="#2e7d32")
    ax.axhline(0, color="#bdbdbd", lw=0.6)
    ax.set_xlabel(r"$k$ (Å$^{-1}$)")
    ax.set_ylabel(r"$k^2\chi(k)$")
    ax.set_title("(b) EXAFS χ(k)", loc="left", fontweight="bold")
    ax.grid(alpha=0.18)
    rows = []
    for kk, chi, k2 in zip(k, data["chi"], y):
        if math.isfinite(float(k2)):
            rows.append({"k_inv_A": float(kk), "chi": float(chi),
                         "k2chi": float(k2)})
    return Panel(
        "Figure1b_exafs", "Figure 1(b)", "k²-weighted EXAFS",
        "k²加权 χ(k)：振荡频率对应键长，振幅对应配位数与无序；高 k 振荡"
        "提示重散射体（金属–金属/金属–重元素）。",
        source=_rel(ctx, data["source"]),
        columns=["k_inv_A", "chi", "k2chi"], rows=rows)


def panel_ft(ax, ctx: dict) -> Panel:
    data = ctx.get("ft")
    if not data:
        _draw_missing(ax, "缺少 *_chir.csv / *_ft.dat")
        return Panel("Figure1c_ft", "Figure 1(c)", "Fourier transform",
                     "FT 幅值与实/虚部：壳层位置和相对强度。",
                     status="skipped", reason="no FT data")
    r, mag = data["r"], data["mag"]
    ax.plot(r, mag, color="#4a148c", lw=1.2, label=r"$|\chi(R)|$")
    ax.plot(r, data["re"], color="#90caf9", lw=0.7, ls="--", label="Re")
    ax.plot(r, data["im"], color="#ffb74d", lw=0.7, ls=":", label="Im")
    mask = (r >= 1.0) & (r <= 3.5)
    if np.any(mask):
        idx = int(np.nanargmax(mag[mask]))
        r_peak = float(r[mask][idx])
        ax.axvline(r_peak, color="#c62828", ls="--", lw=0.8)
        ax.text(r_peak, ax.get_ylim()[1], f" R={r_peak:.2f} Å",
                fontsize=7, va="top", color="#c62828")
    fit = ctx.get("fit") or {}
    rmin = _as_float(fit.get("rmin"))
    rmax = _as_float(fit.get("rmax"))
    if rmin is not None and rmax is not None and rmax > rmin:
        ax.axvspan(rmin, rmax, color="#e1bee7", alpha=0.35, lw=0)
    ax.set_xlim(0, min(float(np.max(r)), 6.0))
    ax.set_xlabel(r"$R$ (Å)")
    ax.set_ylabel(r"$|\chi(R)|$ (Å$^{-3}$)")
    ax.set_title("(c) FT magnitude", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(alpha=0.18)
    rows = []
    for rr, re, im, mm in zip(r, data["re"], data["im"], mag):
        if math.isfinite(float(mm)):
            rows.append({"R_A": float(rr), "chi_R_re": float(re),
                         "chi_R_im": float(im), "magnitude": float(mm)})
    return Panel(
        "Figure1c_ft", "Figure 1(c)", "Fourier transform |χ(R)|",
        "FT 幅值：第一壳 M–O/N 与第二壳 M–M/M–C 的位置和相对强度；"
        "实部/虚部用于与拟合模型比对。峰位未做相移校正。",
        source=_rel(ctx, data["source"]),
        columns=["R_A", "chi_R_re", "chi_R_im", "magnitude"], rows=rows)


def panel_wt(ax, ctx: dict) -> Panel:
    data = ctx.get("wt")
    if not data:
        _draw_missing(ax, "缺少 *_wavelet.csv / *_wt.dat")
        return Panel("Figure1d_wavelet", "Figure 1(d)", "Wavelet",
                     "k–R 小波：区分轻/重散射体与单原子/团簇。",
                     status="skipped", reason="no wavelet data")
    k, r, z = data["k"], data["r"], data["z"]
    mesh = ax.imshow(
        z, extent=[float(r.min()), float(r.max()),
                   float(k.min()), float(k.max())],
        origin="lower", aspect="auto", cmap="viridis",
        interpolation="nearest")
    fig = ax.figure
    cbar = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.02)
    cbar.ax.tick_params(labelsize=6)
    wire = " |W|"
    cbar.set_label(wire, fontsize=7)
    wt_info = (ctx.get("meta") or {}).get("wt") or {}
    r_peak = _as_float(wt_info.get("r_first_peak"))
    k_peak = _as_float(wt_info.get("k_at_peak"))
    if r_peak is None or k_peak is None:
        idx = np.unravel_index(int(np.nanargmax(z)), z.shape)
        k_peak = float(k[idx[0]])
        r_peak = float(r[idx[1]])
    ax.scatter([r_peak], [k_peak], marker="x", s=30, color="white", lw=1.2)
    ax.text(r_peak + 0.08, k_peak,
            f"({r_peak:.2f} Å, {k_peak:.1f} Å$^{{-1}}$)",
            color="white", fontsize=6.5, va="center")
    ax.set_xlabel(r"$R$ (Å)")
    ax.set_ylabel(r"$k$ (Å$^{-1}$)")
    ax.set_title("(d) Wavelet k–R", loc="left", fontweight="bold")
    rows = []
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            val = z[i, j]
            if math.isfinite(float(val)):
                rows.append({"k_inv_A": float(k[i]), "R_A": float(r[j]),
                             "intensity": float(val)})
    return Panel(
        "Figure1d_wavelet", "Figure 1(d)", "Continuous wavelet transform",
        "k–R 小波：轻散射体（N/O/C）在低 k 区，重散射体（Fe/Co/Ni）在"
        "高 k 区；金属–金属键在 2–3 Å 出现高 k 强信号，可用于区分单原子"
        "与金属团簇。",
        source=_rel(ctx, data["source"]),
        columns=["k_inv_A", "R_A", "intensity"], rows=rows)


def panel_fit(ax, ctx: dict) -> Panel:
    rows = _fit_rows(ctx)
    if not rows:
        _draw_missing(ax, "缺少 S7 拟合或正式 CN 结果")
        return Panel("Figure1e_fit", "Figure 1(e)", "Fit shells",
                     "拟合壳层配位数与键长。",
                     status="skipped", reason="no fit / formal CN rows")
    rows = rows[:8]
    labels = []
    for r in rows:
        label = str(r["shell"])
        bits = []
        if r.get("r") is not None:
            bits.append(f"R={r['r']:.2f} Å")
        if r.get("sigma2") is not None:
            bits.append(f"σ²={r['sigma2']:.4f}")
        labels.append(label + ("\n" + "  ".join(bits) if bits else ""))
    n = np.asarray([r["n"] for r in rows], dtype=float)
    n_err = np.asarray([r["n_err"] if r["n_err"] is not None else 0.0
                        for r in rows], dtype=float)
    ypos = np.arange(len(rows))
    colors = [_status_color(r.get("status")) for r in rows]
    ax.barh(ypos, n, color=colors, alpha=0.72, height=0.45)
    ax.errorbar(n, ypos, xerr=n_err, fmt="none", ecolor="#37474f",
                capsize=2.5, lw=0.9)
    n0 = [r["n0"] for r in rows]
    if any(v is not None for v in n0):
        ax.scatter([v if v is not None else np.nan for v in n0], ypos,
                   marker="D", s=28, facecolor="white", edgecolor="#37474f",
                   zorder=4, label=r"$N_0$ (CIF)")
        ax.legend(frameon=False, loc="lower right")
    xmax = max([float(np.nanmax(n + n_err)) if len(n) else 1.0] +
               [float(v) for v in n0 if v is not None] + [1.0])
    for i, r in enumerate(rows):
        if r.get("n_ratio") is not None:
            ax.text(n[i] + n_err[i] + 0.08, i,
                    f"N/N$_0$={r['n_ratio']:.2f}", va="center",
                    fontsize=6.5, color="#37474f")
    ax.set_xlim(0, xmax * 1.45)
    ax.set_yticks(ypos)
    ax.set_yticklabels([_display(x) for x in labels], fontsize=7.0)
    ax.invert_yaxis()
    ax.set_xlabel("Coordination number N")
    mode = str((ctx.get("fit") or {}).get("mode") or "")
    title = "(e) Fit shells" + (f" — {mode}" if mode else "")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(axis="x", alpha=0.2)
    out = []
    for r in rows:
        out.append({k: r.get(k) for k in (
            "shell", "n", "n_err", "n0", "n_ratio", "r", "r_err",
            "sigma2", "rf", "mode", "status", "note", "source")})
    return Panel(
        "Figure1e_fit", "Figure 1(e)", "EXAFS fit / formal CN shells",
        "拟合壳层：N 为 EXAFS 配位数或有效简并度；N₀ 为 CIF 晶体学值；"
        "R 为键长，σ² 为无序；颜色对应可靠性/正式 CN 状态。",
        source="fit.json / cn_table.csv",
        columns=["shell", "n", "n_err", "n0", "n_ratio", "r", "r_err",
                 "sigma2", "rf", "mode", "status", "note", "source"],
        rows=out)


def panel_structure(ax, ctx: dict) -> Panel:
    structure = ctx.get("structure") or {}
    ranked = list(structure.get("ranked") or [])
    if not ranked:
        _draw_missing(ax, "缺少结构识别结果（--identify）")
        return Panel("Figure1f_structure", "Figure 1(f)", "Structure ID",
                     "结构识别候选排名。",
                     status="skipped", reason="no structure identification")
    all_ranked = list(structure.get("ranked") or [])
    ranked = all_ranked[:6]   # 图上只画前 6 个，避免标签重叠
    labels = [str(x.get("id") or "?") for x in ranked]
    scores = [_as_float(x.get("score")) or float("nan") for x in ranked]
    ypos = np.arange(len(ranked))
    colors = ["#1565c0"] + ["#90a4ae"] * (len(ranked) - 1)
    ax.barh(ypos, scores, color=colors, alpha=0.85, height=0.55)
    limit = _as_float((structure.get("confidence") or {}).get("score_limit"))
    if limit is None:
        limit = 0.16
    ax.axvline(limit, color="#c62828", ls="--", lw=0.8)
    for i, x in enumerate(ranked):
        ax.text(float(scores[i]) + 0.01, i,
                f"{float(scores[i]):.3f}", va="center", fontsize=6.5)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Fingerprint score (lower = better)")
    status = structure.get("status") or "unknown"
    coverage = structure.get("coverage_ok")
    ax.set_title(f"(f) Structure ID — {status}, coverage={coverage}; "
                 f"score limit={limit:g}",
                 loc="left", fontweight="bold")
    ax.grid(axis="x", alpha=0.2)
    rows = []
    for i, x in enumerate(all_ranked):
        rows.append({
            "rank": i + 1,
            "id": x.get("id"),
            "kind": x.get("kind"),
            "nuclearity": x.get("nuclearity"),
            "score": _as_float(x.get("score")),
            "score_mag": _as_float(x.get("score_mag")),
            "sigma2": _as_float(x.get("sigma2")),
            "de0": _as_float(x.get("de0")),
            "status": status,
            "coverage_ok": coverage,
        })
    return Panel(
        "Figure1f_structure", "Figure 1(f)", "Structure identification",
        "结构识别候选：分数越低越接近实验谱；红色虚线为 confident 阈值。"
        "no_match 表示诊断库缺条目，不等于数据无效。",
        source=_rel(ctx, ctx.get("structure_path")) if ctx.get("structure_path")
        else "meta.json:structure_id",
        columns=["rank", "id", "kind", "nuclearity", "score", "score_mag",
                 "sigma2", "de0", "status", "coverage_ok"], rows=rows)


def panel_cn_gate(ax, ctx: dict) -> Panel:
    rows = []
    for row in (ctx.get("cn_rows") or []):
        ratio = _as_float(row.get("n_ratio"))
        if ratio is None:
            continue
        rows.append({
            "shell": row.get("shell_label") or row.get("shell_group")
            or row.get("absorber_site") or "shell",
            "n_ratio": ratio,
            "n_status": row.get("n_status"),
            "rf": _as_float(row.get("rf")),
            "gate_status": row.get("gate_status"),
            "reject_reason": row.get("reject_reason") or "",
        })
    formal = _formal_block(ctx)
    if not rows:
        observed = (_dig(ctx.get("meta"), "gates", "formal_cn", "observed")
                    or {})
        ratio = _as_float(observed.get("n_ratio"))
        if ratio is not None:
            rows.append({
                "shell": "first shell",
                "n_ratio": ratio,
                "n_status": observed.get("n_status"),
                "rf": _as_float(observed.get("conditional_rf")),
                "gate_status": formal.get("status") or observed.get("formal_status"),
                "reject_reason": formal.get("reject_reason") or "",
            })
    if not rows:
        status = formal.get("status") or "not run"
        _draw_missing(ax, f"无正式 CN 行（status={status}）")
        return Panel("Figure1g_cn", "Figure 1(g)", "Formal CN gate",
                     "正式 CN 一致性门禁。",
                     status="skipped", reason=f"no formal CN rows ({status})")
    xmax = max(1.6, max(float(r["n_ratio"]) for r in rows) + 0.15)
    ax.axvspan(0, 0.6, color=BAND_COLORS["reject"], alpha=0.6, lw=0)
    ax.axvspan(0.6, 0.7, color=BAND_COLORS["warn"], alpha=0.8, lw=0)
    ax.axvspan(0.7, 1.3, color=BAND_COLORS["ok"], alpha=0.7, lw=0)
    ax.axvspan(1.3, 1.4, color=BAND_COLORS["warn"], alpha=0.8, lw=0)
    ax.axvspan(1.4, xmax, color=BAND_COLORS["reject"], alpha=0.6, lw=0)
    for x in (0.6, 0.7, 1.3, 1.4):
        ax.axvline(x, color="#546e7a", ls=":", lw=0.7)
    ypos = np.arange(len(rows))
    xs = [r["n_ratio"] for r in rows]
    colors = [_status_color(r.get("gate_status")) for r in rows]
    ax.scatter(xs, ypos, s=70, c=colors, edgecolor="#37474f", zorder=4)
    for i, r in enumerate(rows):
        bits = [f"N/N$_0$={r['n_ratio']:.3f}"]
        if r.get("rf") is not None:
            bits.append(f"Rf={r['rf']:.4f}")
        if r.get("n_status"):
            bits.append(str(r["n_status"]))
        ax.text(r["n_ratio"] + 0.02, i, "  ".join(bits), va="center",
                fontsize=6.8)
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.6, max(1, len(rows) - 0.4))
    ax.set_yticks(ypos)
    ax.set_yticklabels([r["shell"] for r in rows], fontsize=7.5)
    ax.set_xlabel(r"$N/N_0$")
    status = formal.get("status") or "unknown"
    ax.set_title(f"(g) Formal CN gate — {status}", loc="left",
                 fontweight="bold")
    ax.grid(axis="x", alpha=0.2)
    return Panel(
        "Figure1g_cn", "Figure 1(g)", "Formal CN consistency gate",
        "正式 CN 一致性：0.7–1.3 通过，0.6–0.7/1.3–1.4 警告，越界拒绝；"
        "同时标注条件拟合 Rf。拒签不会删除诊断拟合。",
        source="cn_table.csv / meta.json:gates.formal_cn",
        columns=["shell", "n_ratio", "n_status", "rf", "gate_status",
                 "reject_reason"], rows=rows)


def _card_rows(ctx: dict) -> list[dict]:
    meta = ctx.get("meta") or {}
    fit = ctx.get("fit") or {}
    formal = _formal_block(ctx)
    params = (meta.get("params") or {}).get("effective") or {}
    source = meta.get("source") or {}
    structure = ctx.get("structure") or {}
    top = ((structure.get("ranked") or [{}])[0] if structure else {})
    model = (formal.get("model") or {}) if isinstance(formal, dict) else {}
    rows = [
        ("Sample", _sample_name(ctx)),
        ("Absorber / edge", f"{_meta_absorber(meta)} / {_meta_edge(meta, _meta_absorber(meta))}"),
        ("E0 (eV)", _fmt(_meta_e0(meta), 6)),
        ("Phase ID", model.get("model_id") or params.get("resolved_phase_id")
         or params.get("phase_id") or top.get("id") or "—"),
        ("Model level", model.get("model_level") or params.get("model_level") or "—"),
        ("Source type", model.get("source_type") or params.get("resolved_phase_source") or "—"),
        ("Absorber site", model.get("absorber_site") or params.get("resolved_absorber_site")
         or params.get("absorber_site") or "—"),
        ("S0²", f"{_fmt(fit.get('s02') or params.get('s02'), 5)} ({fit.get('s02_source') or params.get('s02_source') or '—'})"),
        ("Fit mode", fit.get("mode") or "—"),
        ("CN status", fit.get("cn_status") or "—"),
        ("Formal CN", formal.get("status") or "not run"),
        ("Data SHA256", _short(source.get("data_sha256") or meta.get("data_sha256"), 18)),
        ("CIF SHA256", _short(model.get("source_cif_sha256"), 18)),
    ]
    return [{"key": k, "value": v} for k, v in rows]


def panel_card(ax, ctx: dict) -> Panel:
    rows = _card_rows(ctx)
    _render_table(ax, ["Field", "Value"],
                  [[r["key"], _short(r["value"], 52)] for r in rows],
                  fontsize=7.2, title="(h) Structure card")
    return Panel(
        "Figure1h_card", "Figure 1(h)", "Structure card",
        "结构卡片：样品、吸收元素/边、相来源、模型级别、吸收位、S0²、"
        "拟合与正式 CN 状态、数据/CIF 哈希。",
        source="meta.json", columns=["key", "value"], rows=rows)


# ---------------------------------------------------------------
# 图 2 面板
# ---------------------------------------------------------------
def _gate_rows(ctx: dict) -> list[dict]:
    gates = (ctx.get("meta") or {}).get("gates") or {}
    rows = []
    for name in GATE_ORDER:
        payload = gates.get(name)
        if not isinstance(payload, dict):
            continue
        rows.append({
            "gate": name,
            "status": payload.get("status") or "not_available",
            "action": payload.get("action") or "—",
            "observed": _short(payload.get("observed"), 70),
            "threshold": _short(payload.get("threshold"), 70),
            "reject_reason": payload.get("reject_reason") or "",
        })
    if rows:
        return rows
    summary = _dig(ctx.get("meta"), "gates", "summary") or {}
    for status, names in (("passed", summary.get("passed")),
                          ("warning", summary.get("warned")),
                          ("failed", summary.get("failed")),
                          ("not_applicable", summary.get("not_applicable"))):
        for name in (names or []):
            rows.append({"gate": name, "status": status, "action": "—",
                         "observed": "—", "threshold": "—",
                         "reject_reason": ""})
    return rows


def panel_gate_board(ax, ctx: dict) -> Panel:
    rows = _gate_rows(ctx)
    if not rows:
        _draw_missing(ax, "缺少 gates 台账（P6.4 输出）")
        return Panel("Figure2a_gates", "Figure 2(a)", "Gate status board",
                     "门禁状态板。", status="skipped", reason="no gates ledger")
    ax.axis("off")
    cols = ["passed", "warning", "failed", "not_applicable", "not_available"]
    nrow = len(rows)
    ax.set_xlim(-2.8, len(cols) + 0.9)
    ax.set_ylim(0, nrow + 2.0)
    ax.text(-0.15, nrow + 1.35, "Gate", fontsize=7.5, fontweight="bold",
            va="center", ha="right")
    for j, col in enumerate(cols):
        ax.text(j + 0.55, nrow + 1.35, col.replace("_", "\n"),
                fontsize=6.0, fontweight="bold", va="center", ha="center")
    for i, row in enumerate(rows):
        y = nrow - i - 0.5
        ax.text(-0.15, y, row["gate"], fontsize=7, va="center", ha="right")
        for j, col in enumerate(cols):
            active = (row["status"] == col)
            color = _status_color(col) if active else "#f5f5f5"
            rect = Rectangle((j + 0.10, y - 0.34), 0.90, 0.68,
                             facecolor=color, edgecolor="#cfd8dc", lw=0.5)
            ax.add_patch(rect)
            if active:
                ax.text(j + 0.55, y, "●", ha="center", va="center",
                        fontsize=6, color="white" if col in (
                            "passed", "failed") else "#212121")
    summary = _summary_counts((ctx.get("meta") or {}).get("gates") or {})
    text = (f"passed={len(summary['passed'])}  "
            f"warning={len(summary['warned'])}  "
            f"failed={len(summary['failed'])}  "
            f"n/a={len(summary['not_applicable'])}")
    ax.text(-0.15, 0.15, text, fontsize=7, color="#37474f", ha="left")
    ax.set_title("(a) Gate status board", loc="left", fontweight="bold")
    return Panel(
        "Figure2a_gates", "Figure 2(a)", "Gate status board",
        "置信门禁状态板：通过/警告/失败/不适用一眼可见。",
        source="meta.json:gates",
        columns=["gate", "status", "action", "observed", "threshold",
                 "reject_reason"], rows=rows)


def panel_decision_flow(ax, ctx: dict) -> Panel:
    decisions = list((ctx.get("meta") or {}).get("decisions") or [])
    if not decisions:
        _draw_missing(ax, "缺少 decisions 台账（P6.4 输出）")
        return Panel("Figure2b_decisions", "Figure 2(b)", "Decision trace",
                     "决策轨迹。", status="skipped", reason="no decisions ledger")
    decisions = decisions[:12]
    ax.axis("off")
    n = len(decisions)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, n + 0.5)
    ax.plot([0.08, 0.08], [0.3, n - 0.2], color="#b0bec5", lw=1.4, zorder=1)
    for i, item in enumerate(decisions):
        y = n - i - 0.5
        status = str(item.get("status") or "info")
        color = _status_color(status)
        ax.scatter([0.08], [y], s=80, color=color, zorder=3,
                   edgecolor="white", lw=0.8)
        ax.text(0.0, y, str(item.get("step") or ""), fontsize=7,
                ha="right", va="center", fontweight="bold")
        detail = _display(_short(item.get("detail"), 54))
        ax.text(0.14, y + 0.12, _display(str(item.get("decision") or "")),
                fontsize=7, va="center", fontweight="bold")
        ax.text(0.14, y - 0.16, detail, fontsize=6.0, va="center",
                color="#546e7a")
    ax.set_title("(b) Decision trace", loc="left", fontweight="bold")
    rows = [{
        "step": item.get("step"),
        "decision": item.get("decision"),
        "status": item.get("status"),
        "detail": _short(item.get("detail"), 120),
    } for item in decisions]
    return Panel(
        "Figure2b_decisions", "Figure 2(b)", "Decision trace",
        "决策轨迹：S1→S7→正式 CN 的关键分支及其状态。",
        source="meta.json:decisions",
        columns=["step", "decision", "status", "detail"], rows=rows)


def panel_param_sources(ax, ctx: dict) -> Panel:
    params = (ctx.get("meta") or {}).get("params") or {}
    sources = params.get("sources") or {}
    stages = params.get("stages") or {}
    if not sources:
        _draw_missing(ax, "缺少 params 台账（P6.4 输出）")
        return Panel("Figure2c_sources", "Figure 2(c)", "Parameter sources",
                     "参数来源。", status="skipped", reason="no params ledger")
    counts: dict[tuple[str, str], int] = {}
    for key, source in sources.items():
        stage = stages.get(key) or "unspecified"
        counts[(stage, source)] = counts.get((stage, source), 0) + 1
    stage_totals: dict[str, int] = {}
    for (stage, _), n in counts.items():
        stage_totals[stage] = stage_totals.get(stage, 0) + n
    top_stages = [s for s, _ in sorted(stage_totals.items(),
                                       key=lambda kv: (-kv[1], kv[0]))[:10]]
    src_order = [s for s in ("cli", "config", "default", "auto",
                             "calibrated", "runtime_override")
                 if any(src == s for _, src in counts)]
    ypos = np.arange(len(top_stages))
    left = np.zeros(len(top_stages))
    for source in src_order:
        vals = np.asarray([counts.get((stage, source), 0)
                           for stage in top_stages], dtype=float)
        ax.barh(ypos, vals, left=left, color=SOURCE_COLORS.get(source, "#90a4ae"),
                alpha=0.88, height=0.55, label=source)
        left += vals
    ax.set_yticks(ypos)
    ax.set_yticklabels(top_stages, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Number of parameters")
    ax.set_title("(c) Parameter provenance", loc="left", fontweight="bold")
    ax.legend(frameon=False, ncol=3, fontsize=6.5, loc="lower right")
    ax.grid(axis="x", alpha=0.18)
    rows = [{"stage": stage, "source": source, "count": count}
            for (stage, source), count in sorted(counts.items())]
    return Panel(
        "Figure2c_sources", "Figure 2(c)", "Parameter provenance",
        "参数来源：CLI / 配置文件 / 默认值 / 自动推导 / 标定 / 运行时覆盖"
        "按阶段统计，回答“这个结果用了哪些参数、从哪来”。",
        source="meta.json:params",
        columns=["stage", "source", "count"], rows=rows)


def _gate_observed_summary(name: str, payload: dict) -> str:
    """把门禁 observed 压成表格里可读的一行；完整 JSON 留在 CSV。"""
    obs = payload.get("observed")
    if obs is None:
        return "—"
    if name == "normalization":
        return f"plateau={_fmt(obs, 4)}"
    if name == "self_absorption" and isinstance(obs, dict):
        return (f"mode={obs.get('mode')}, risk={obs.get('risk')}, "
                f"corrected={obs.get('corrected')}")
    if name == "structure_id" and isinstance(obs, dict):
        return (f"status={obs.get('status')}, best={obs.get('best_id')}, "
                f"score={_fmt(obs.get('best_score'), 4)}, "
                f"coverage={obs.get('coverage_ok')}")
    if name == "cif" and isinstance(obs, dict):
        sites = obs.get("sites") or []
        return f"ok={obs.get('ok')}, sites={len(sites)}"
    if name == "fit" and isinstance(obs, dict):
        bits = [f"mode={obs.get('mode')}", f"Rf={_fmt(obs.get('r_factor'), 4)}"]
        if obs.get("delta_e0") is not None:
            bits.append(f"dE0={_fmt(obs.get('delta_e0'), 3)}")
        bits.append(f"CN_reliable={obs.get('cn_reliable')}")
        return ", ".join(bits)
    if name == "formal_cn" and isinstance(obs, dict):
        bits = [f"status={obs.get('formal_status')}"]
        if obs.get("n_ratio") is not None:
            bits.append(f"N/N0={_fmt(obs.get('n_ratio'), 4)}")
        if obs.get("conditional_rf") is not None:
            bits.append(f"Rf={_fmt(obs.get('conditional_rf'), 4)}")
        return ", ".join(bits)
    if name == "xanes_cross_check" and isinstance(obs, dict):
        return f"warnings={obs.get('n_warnings')}"
    return _short(obs, 42)


def panel_gate_table(ax, ctx: dict) -> Panel:
    rows = _gate_rows(ctx)
    if not rows:
        _draw_missing(ax, "缺少 gates 台账")
        return Panel("Figure2d_gate_table", "Figure 2(d)", "Gate table",
                     "门禁表。", status="skipped", reason="no gates ledger")
    gates = (ctx.get("meta") or {}).get("gates") or {}
    table_rows = [[r["gate"], r["status"], r["action"],
                   _gate_observed_summary(r["gate"], gates.get(r["gate"]) or {})]
                  for r in rows]
    _render_table(ax, ["Gate", "Status", "Action", "Key observation"],
                  table_rows, status_col=1, fontsize=6.2,
                  title="(d) Gate ledger — full observed/threshold in CSV")
    return Panel(
        "Figure2d_gate_table", "Figure 2(d)", "Gate table",
        "门禁表：每个门禁的阈值、观测值、状态和动作，直接作为图片输出。",
        source="meta.json:gates",
        columns=["gate", "status", "action", "observed", "threshold"],
        rows=rows)


def panel_decision_table(ax, ctx: dict) -> Panel:
    decisions = list((ctx.get("meta") or {}).get("decisions") or [])[:12]
    if not decisions:
        _draw_missing(ax, "缺少 decisions 台账")
        return Panel("Figure2e_decision_table", "Figure 2(e)", "Decision table",
                     "决策表。", status="skipped", reason="no decisions ledger")
    rows = [{
        "step": item.get("step"),
        "decision": item.get("decision"),
        "status": item.get("status"),
        "detail": _short(item.get("detail"), 46),
    } for item in decisions]
    _render_table(ax, ["Step", "Decision", "Status", "Detail"],
                  [[r["step"], r["decision"], r["status"], r["detail"]]
                   for r in rows], status_col=2, fontsize=6.1,
                  title="(e) Decision ledger")
    return Panel(
        "Figure2e_decision_table", "Figure 2(e)", "Decision table",
        "决策表：S1→S7→正式 CN 的决策、状态与细节，直接作为图片输出。",
        source="meta.json:decisions",
        columns=["step", "decision", "status", "detail"], rows=rows)


def _key_params(ctx: dict) -> list[dict]:
    params = (ctx.get("meta") or {}).get("params") or {}
    effective = params.get("effective") or {}
    sources = params.get("sources") or {}
    units = params.get("units") or {}
    stages = params.get("stages") or {}
    formal = set(params.get("affects_formal_cn") or [])
    preferred = [
        "absorber", "edge_energy", "e0_method", "manual_e0", "post_window",
        "strict", "rbkg", "kmin", "kmax", "kweight", "s02", "s02_from_t5",
        "phase_id", "phase_cif", "absorber_site", "fit", "fit_motif",
        "fit_rmin", "fit_rmax", "formal_cn", "identify",
    ]
    keys = [k for k in preferred if k in effective]
    keys += [k for k in sorted(effective)
             if k in formal and k not in keys]
    keys = keys[:14]
    rows = []
    for key in keys:
        rows.append({
            "parameter": key,
            "value": _short(effective.get(key), 38),
            "source": sources.get(key) or "—",
            "unit": units.get(key) or "",
            "stage": stages.get(key) or "",
            "affects_formal_cn": key in formal,
        })
    return rows


def panel_param_table(ax, ctx: dict) -> Panel:
    rows = _key_params(ctx)
    if not rows:
        _draw_missing(ax, "缺少 params 台账")
        return Panel("Figure2f_param_table", "Figure 2(f)", "Key parameters",
                     "关键参数表。", status="skipped", reason="no params ledger")
    _render_table(
        ax, ["Parameter", "Value", "Source", "Unit", "Stage"],
        [[r["parameter"], r["value"], r["source"], r["unit"], r["stage"]]
         for r in rows], fontsize=6.6, title="(f) Key parameters")
    return Panel(
        "Figure2f_param_table", "Figure 2(f)", "Key parameter table",
        "关键参数表：影响正式 CN 的参数及其来源/单位/阶段，直接作为图片"
        "输出；完整参数见 CSV 与 meta.json。",
        source="meta.json:params",
        columns=["parameter", "value", "source", "unit", "stage",
                 "affects_formal_cn"], rows=rows)


# ---------------------------------------------------------------
# 图 3：拟合与误差诊断
# ---------------------------------------------------------------
def _short_param_name(name: str) -> str:
    """``Fe-Fe@2.48:ΔR`` → ``ΔR@2.48``，用于热图轴标签。"""
    text = str(name)
    if ":" in text:
        shell, param = text.rsplit(":", 1)
        if "@" in shell:
            shell = shell.split("@")[-1]
        return f"{param}@{shell}"
    return text


def panel_fit_correlation(ax, ctx: dict) -> Panel:
    obj = _fit_obj(ctx)
    corr = obj.get("correlation")
    names = [str(x) for x in (obj.get("parameter_names") or [])]
    if not corr or not names:
        _draw_missing(ax, "缺少参数相关矩阵（fit.json:fit.correlation）")
        return Panel("Figure3a_fitcorr", "Figure 3(a)", "Parameter correlation",
                     "拟合参数相关矩阵。", status="skipped",
                     reason="no fit correlation matrix")
    mat = np.asarray(corr, dtype=float)
    n = min(len(names), mat.shape[0], mat.shape[1])
    mat, names = mat[:n, :n], names[:n]
    image = ax.imshow(mat, cmap="RdBu_r", vmin=-1.0, vmax=1.0,
                      aspect="equal")
    fig = ax.figure
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.02)
    cbar.ax.tick_params(labelsize=6.5)
    cbar.set_label("Pearson ρ", fontsize=7)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    short = [_short_param_name(name) for name in names]
    ax.set_xticklabels(short, rotation=45, ha="right", fontsize=6.0)
    ax.set_yticklabels(short, fontsize=6.0)
    if n <= 8:
        for i in range(n):
            for j in range(n):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                        fontsize=5.2,
                        color="white" if abs(mat[i, j]) > 0.65 else "#212121")
    ax.set_title("(a) Parameter correlation", loc="left", fontweight="bold")
    rows = [{"param_i": names[i], "param_j": names[j],
             "rho": float(mat[i, j])}
            for i in range(n) for j in range(n)]
    return Panel(
        "Figure3a_fitcorr", "Figure 3(a)", "Parameter correlation matrix",
        "参数相关矩阵：|ρ|→1 表示两参数几乎不可独立确定；强相关时 N 与 σ²"
        "的误差应视为联合不确定度，不能单独引用。",
        source="fit.json:fit.correlation",
        columns=["param_i", "param_j", "rho"], rows=rows)


def panel_fit_paths(ax, ctx: dict) -> Panel:
    rows = _fit_paths(ctx)
    if not rows:
        _draw_missing(ax, "缺少路径参数（fit.json:fit.paths / l1）")
        return Panel("Figure3b_paths", "Figure 3(b)", "Path N–R map",
                     "拟合路径 N–R 图。", status="skipped",
                     reason="no fit paths")
    for row in rows:
        n = row.get("n")
        r = row.get("r")
        if n is None or r is None:
            continue
        sigma2 = row.get("sigma2") or 0.0
        size = 40.0 + max(0.0, min(0.02, sigma2)) * 5000.0
        ax.errorbar(r, n, xerr=row.get("r_err"), yerr=row.get("n_err"),
                    fmt="o", ms=math.sqrt(size) / 2.2, color="#1565c0",
                    ecolor="#546e7a", capsize=2.5, lw=0.9, zorder=3)
        label = row.get("label") or "path"
        dr = row.get("delta_r")
        s2 = row.get("sigma2")
        extra = []
        if dr is not None:
            extra.append(f"ΔR={float(dr):+.3f}Å")
        if s2 is not None:
            extra.append(f"σ²={float(s2) * 1e3:.2f}e-3")
        if extra:
            label += "\n" + "  ".join(extra)
        ax.annotate(label, (r, n), textcoords="offset points", xytext=(6, 5),
                    fontsize=6.0, color="#263238")
        n0 = row.get("n0")
        if n0 is not None:
            ax.scatter([r], [n0], marker="D", s=26, facecolor="white",
                       edgecolor="#c62828", zorder=4)
    ax.set_xlabel(r"$R$ (Å)")
    ax.set_ylabel("Coordination number N")
    ax.set_title("(b) Path N–R (open diamond = CIF N$_0$)", loc="left",
                 fontweight="bold")
    ax.grid(alpha=0.2)
    return Panel(
        "Figure3b_paths", "Figure 3(b)", "Fit path N–R map",
        "路径级 N–R 图：误差棒为形式误差，点大小为 σ²；空心菱形为 CIF "
        "N₀。强相关或误差过大时不能把 N 当绝对配位数。",
        source="fit.json:fit.paths",
        columns=["variant", "label", "n", "n_err", "r", "r_err", "reff",
                 "delta_r", "sigma2", "sigma2_err", "n0", "n_ratio", "rf",
                 "note"], rows=rows)


def panel_error_budget(ax, ctx: dict) -> Panel:
    systematics = (ctx.get("fit") or {}).get("cn_systematics") or {}
    components = dict(systematics.get("components") or {})
    extra = dict(systematics.get("extra_components") or {})
    items = [(str(k), _as_float(v)) for k, v in components.items()]
    items += [(f"extra:{k}", _as_float(v)) for k, v in extra.items()]
    items = [(k, v) for k, v in items if v is not None]
    total = _as_float(systematics.get("relative_1sigma"))
    if not items and total is None:
        _draw_missing(ax, "缺少 CN 系统误差预算（fit.json:cn_systematics）")
        return Panel("Figure3c_budget", "Figure 3(c)", "CN error budget",
                     "CN 系统误差预算。", status="skipped",
                     reason="no cn_systematics")
    ypos = np.arange(len(items))
    values = [100.0 * v for _, v in items]
    ax.barh(ypos, values, color="#ef6c00", alpha=0.8, height=0.55)
    ax.set_yticks(ypos)
    ax.set_yticklabels([k.replace("_", " ") for k, _ in items], fontsize=6.8)
    ax.invert_yaxis()
    if total is not None:
        ax.axvline(100.0 * total, color="#c62828", ls="--", lw=1.0)
        ax.text(100.0 * total, -0.55, f" total={100.0 * total:.1f}%",
                color="#c62828", fontsize=6.8, va="bottom")
    for i, (_, value) in enumerate(items):
        ax.text(100.0 * value + 0.5, i, f"{100.0 * value:.1f}%",
                va="center", fontsize=6.3)
    ax.set_xlabel("Relative 1σ contribution (%)")
    ax.set_title("(c) CN systematic error budget", loc="left",
                 fontweight="bold")
    ax.grid(axis="x", alpha=0.2)
    note = systematics.get("note") or ""
    rows = [{"component": k, "relative_1sigma": v,
             "percent": 100.0 * v,
             "total_relative_1sigma": total, "note": note}
            for k, v in items]
    return Panel(
        "Figure3c_budget", "Figure 3(c)", "CN systematic error budget",
        "CN 系统误差预算：各项按平方和合成总相对 1σ；即使形式误差很小，"
        "系统误差仍可能限制 N 的绝对精度。",
        source="fit.json:cn_systematics",
        columns=["component", "relative_1sigma", "percent",
                 "total_relative_1sigma", "note"], rows=rows)


def panel_rbkg(ax, ctx: dict) -> Panel:
    meta = ctx.get("meta") or {}
    table = _dig(meta, "rbkg", "table") or []
    if not table:
        _draw_missing(ax, "缺少 Rbkg 寻优表（meta.json:rbkg.table）")
        return Panel("Figure3d_rbkg", "Figure 3(d)", "Rbkg optimization",
                     "Rbkg 寻优曲线。", status="skipped",
                     reason="no rbkg table")
    rows = []
    for item in table:
        row = {"rbkg": _as_float(item.get("rbkg")),
               "r1": _as_float(item.get("r1")),
               "amplitude": _as_float(item.get("amplitude")),
               "leakage": _as_float(item.get("leakage")),
               "ratio": _as_float(item.get("ratio")),
               "score": _as_float(item.get("score")),
               "veto": item.get("veto"),
               "confidence": item.get("confidence")}
        metrics = item.get("confidence_metrics") or {}
        for key in ("valid_fraction", "roughness", "margin"):
            row[f"confidence_{key}"] = _as_float(metrics.get(key))
        rows.append(row)
    rows = [r for r in rows if r["rbkg"] is not None]
    if not rows:
        _draw_missing(ax, "Rbkg 表为空")
        return Panel("Figure3d_rbkg", "Figure 3(d)", "Rbkg optimization",
                     "Rbkg 寻优曲线。", status="skipped",
                     reason="empty rbkg table")
    x = [r["rbkg"] for r in rows]
    score = [r["score"] if r["score"] is not None else np.nan for r in rows]
    leakage = [r["leakage"] if r["leakage"] is not None else np.nan
               for r in rows]
    ax.plot(x, score, "o-", color="#1565c0", lw=1.1, ms=3.5,
            label="score (↑ better)")
    ax.set_xlabel(r"Rbkg (Å)")
    ax.set_ylabel("Rbkg score", color="#1565c0")
    ax.tick_params(axis="y", labelcolor="#1565c0")
    ax2 = ax.twinx()
    ax2.plot(x, leakage, "s--", color="#ef6c00", lw=0.9, ms=3.0,
             label="leakage (↓ better)")
    ax2.set_ylabel("Leakage", color="#ef6c00")
    ax2.tick_params(axis="y", labelcolor="#ef6c00")
    selected = _as_float(_dig(meta, "rbkg", "selected"))
    confidence = _dig(meta, "rbkg", "confidence")
    if selected is not None:
        ax.axvline(selected, color="#c62828", ls=":", lw=1.1)
        ax.text(selected, ax.get_ylim()[1],
                f" selected={selected:.2f} Å ({confidence or '—'})",
                color="#c62828", fontsize=6.5, va="top")
    ax.set_title("(d) Rbkg optimization", loc="left", fontweight="bold")
    ax.grid(alpha=0.2)
    return Panel(
        "Figure3d_rbkg", "Figure 3(d)", "Rbkg optimization curve",
        "Rbkg 寻优：score 越高、leakage 越低越好；选中值由自动判据或用户"
        "显式指定，图中红色虚线为最终采用值。",
        source="meta.json:rbkg.table",
        columns=["rbkg", "r1", "amplitude", "leakage", "ratio", "score",
                 "veto", "confidence", "confidence_valid_fraction",
                 "confidence_roughness", "confidence_margin"], rows=rows)


def _xanes_columns(ctx: dict) -> dict | None:
    data = ctx.get("xanes")
    if not data:
        return None
    path = Path(data["source"])
    try:
        header, arr = _read_numeric(path)
    except FigurePackError:
        return None
    if arr.ndim != 2 or arr.shape[0] < 3:
        return None
    idx = {str(h).strip(): i for i, h in enumerate(header or [])}
    energy = arr[:, idx.get("energy_eV", 0)]
    if header:
        mu = arr[:, idx["mu"]] if "mu" in idx else None
        pre = arr[:, idx["pre_edge_line"]] if "pre_edge_line" in idx else None
        norm = arr[:, idx["mu_norm"]] if "mu_norm" in idx else None
        flat = arr[:, idx["mu_flat"]] if "mu_flat" in idx else None
    else:
        mu = arr[:, 1] if arr.shape[1] > 1 else None
        pre = None
        norm = arr[:, 2] if arr.shape[1] > 2 else mu
        flat = norm
    return {"energy": energy, "mu": mu, "pre": pre,
            "norm": norm, "flat": flat, "source": path}


def panel_normalization(ax, ctx: dict) -> Panel:
    cols = _xanes_columns(ctx)
    if cols is None:
        _draw_missing(ax, "缺少 μ(E) 数据")
        return Panel("Figure3e_normalization", "Figure 3(e)",
                     "Normalization", "归一化分解。", status="skipped",
                     reason="no mu(E) data")
    energy = cols["energy"]
    mu = cols["mu"] if cols["mu"] is not None else cols["norm"]
    pre = cols["pre"]
    flat = cols["flat"] if cols["flat"] is not None else cols["norm"]
    if mu is None and flat is None:
        _draw_missing(ax, "μ(E) 列不可用")
        return Panel("Figure3e_normalization", "Figure 3(e)",
                     "Normalization", "归一化分解。", status="skipped",
                     reason="no usable mu columns")
    ax.plot(energy, mu, color="#37474f", lw=1.0, label=r"raw $\mu(E)$")
    if pre is not None:
        ax.plot(energy, pre, color="#1565c0", lw=0.8, ls="--",
                label="pre-edge line")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel(r"$\mu(E)$", color="#37474f")
    ax.tick_params(axis="y", labelcolor="#37474f")
    ax2 = ax.twinx()
    ax2.plot(energy, flat, color="#2e7d32", lw=1.0,
             label=r"normalized $\mu$")
    ax2.set_ylabel("Normalized μ", color="#2e7d32")
    ax2.tick_params(axis="y", labelcolor="#2e7d32")
    e0 = _meta_e0(ctx.get("meta") or {})
    if e0 is not None:
        ax.axvline(e0, color="#c62828", ls=":", lw=0.9)
        ax.text(e0, ax.get_ylim()[1], f" E0={e0:.1f}", color="#c62828",
                fontsize=6.5, va="top")
    norm_meta = (ctx.get("meta") or {}).get("normalization") or {}
    edge_step = _as_float(norm_meta.get("edge_step"))
    plateau = _as_float(norm_meta.get("plateau"))
    bits = []
    if edge_step is not None:
        bits.append(f"edge step={edge_step:.4g}")
    if plateau is not None:
        bits.append(f"plateau={plateau:.4g}")
    lines = []
    handles = ax.get_lines() + ax2.get_lines()
    for handle in handles:
        if handle.get_label() and not handle.get_label().startswith("_"):
            lines.append(handle)
    ax.legend(handles=lines, frameon=True, framealpha=0.88, fontsize=6.2,
              loc="center right")
    title = "(e) Normalization decomposition"
    if bits:
        title += " — " + ", ".join(bits)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(alpha=0.18)
    rows = []
    for i, e in enumerate(energy):
        rows.append({
            "energy_eV": float(e),
            "mu": _as_float(mu[i]) if mu is not None else None,
            "pre_edge_line": (_as_float(pre[i]) if pre is not None else None),
            "mu_norm": (_as_float(cols["norm"][i])
                        if cols["norm"] is not None else None),
            "mu_flat": (_as_float(flat[i]) if flat is not None else None),
        })
    return Panel(
        "Figure3e_normalization", "Figure 3(e)", "Normalization decomposition",
        "归一化分解：原始 μ(E)、边前线与归一化 μ；边跃高度和平台值决定"
        "定量可靠性，平台偏离 0.85–1.15 会触发归一化警告。",
        source=_rel(ctx, cols["source"]),
        columns=["energy_eV", "mu", "pre_edge_line", "mu_norm", "mu_flat"],
        rows=rows)


def panel_fit_window_table(ax, ctx: dict) -> Panel:
    fit = ctx.get("fit") or {}
    obj = _fit_obj(ctx)
    rows: list[dict] = []

    def add(key: str, value: Any, unit: str = "") -> None:
        if value is None:
            return
        rows.append({"key": key, "value": _fmt(value, 6), "unit": unit})

    add("Fit mode", fit.get("mode"))
    add("Status", fit.get("status"))
    add("kmin", obj.get("kmin", _dig((ctx.get("meta") or {}), "ft", "kmin")),
        "Å^-1")
    add("kmax", obj.get("kmax", _dig((ctx.get("meta") or {}), "ft", "kmax")),
        "Å^-1")
    kweights = obj.get("kweights")
    add("kweight", ",".join(str(x) for x in kweights) if kweights else None)
    add("Rmin", obj.get("rmin", fit.get("rmin")), "Å")
    add("Rmax", obj.get("rmax", fit.get("rmax")), "Å")
    add("S0²", fit.get("s02", obj.get("s02")), "")
    delta_e0 = _as_float(obj.get("delta_e0"))
    delta_e0_err = _as_float(obj.get("delta_e0_err"))
    if delta_e0 is not None:
        add("ΔE0 ± err",
            (f"{delta_e0:.4g} ± {delta_e0_err:.4g}"
             if delta_e0_err is not None else f"{delta_e0:.4g}"), "eV")
    add("R-factor", obj.get("r_factor"), "")
    add("N independent", obj.get("n_independent"), "")
    add("N variables", obj.get("n_variables"), "")
    chi2 = _as_float(obj.get("chi_square"))
    red_chi2 = _as_float(obj.get("reduced_chi_square"))
    if chi2 is not None:
        add("χ² / reduced",
            (f"{chi2:.4g} / {red_chi2:.4g}"
             if red_chi2 is not None else f"{chi2:.4g}"), "")
    add("N paths", len(obj.get("paths") or []), "")
    add("CN status", fit.get("cn_status"))
    add("CN reliable", fit.get("cn_reliable"))
    if not rows:
        _draw_missing(ax, "缺少拟合结果（fit.json）")
        return Panel("Figure3f_window_table", "Figure 3(f)",
                     "Fit window & settings", "拟合窗口表。",
                     status="skipped", reason="no fit result")
    _render_table(ax, ["Item", "Value", "Unit"],
                  [[r["key"], r["value"], r["unit"]] for r in rows],
                  fontsize=6.6, title="(f) Fit window & settings")
    return Panel(
        "Figure3f_window_table", "Figure 3(f)", "Fit window & settings",
        "拟合窗口与设置：k/R 范围、k 权重、S0²、ΔE0、Rf、独立点数与变量数"
        "共同决定拟合是否可信；同时给出 S7 的 CN 状态。",
        source="fit.json:fit",
        columns=["key", "value", "unit"], rows=rows)


def panel_checks_table(ax, ctx: dict) -> Panel:
    fit = ctx.get("fit") or {}
    meta = ctx.get("meta") or {}
    flags = fit.get("flags") or {}
    rows: list[dict] = []

    def add(check: str, status: str, value: Any, note: str = "") -> None:
        rows.append({"check": check, "status": status,
                     "value": _short(value, 48),
                     "note": _short(note, 28)})

    risk_flags = (
        ("overfit", "overfit"),
        ("sigma2_at_bound", "σ² at bound"),
        ("c3_at_bound", "C3 at bound"),
        ("c4_at_bound", "C4 at bound"),
        ("weight_at_bound", "weight at bound"),
        ("weight_degenerate", "weight degenerate"),
        ("sub_rayleigh_z_degeneracy", "sub-Rayleigh Z"),
        ("c3_dropped_by_nidp", "C3 dropped"),
        ("c4_dropped_by_nidp", "C4 dropped"),
    )
    for key, label in risk_flags:
        if key not in flags:
            continue
        value = bool(flags.get(key))
        if key in ("c3_dropped_by_nidp", "c4_dropped_by_nidp"):
            status = "warning" if value else "passed"
        else:
            status = "failed" if value else "passed"
        add(label, status, value)

    sa = fit.get("self_absorption") or meta.get("self_absorption") or {}
    if sa:
        cn_gate = bool(sa.get("cn_gate"))
        add("self-absorption", "warning" if cn_gate else "passed",
            f"risk={sa.get('risk')}, corrected={sa.get('corrected')}, "
            f"cn_gate={cn_gate}", sa.get("note") or "")

    resolution = fit.get("resolution_check") or {}
    if resolution:
        add("shell resolution",
            "passed" if resolution.get("resolvable") else "warning",
            f"Δk={_fmt(resolution.get('delta_k'), 4)} Å^-1, "
            f"Rayleigh={_fmt(resolution.get('rayleigh_limit_A'), 4)} Å",
            resolution.get("recommendation") or "")

    for key, label in (("mm_bond_ftest", "M–M bond F-test"),
                       ("shell_f_test", "shell F-test")):
        test = fit.get(key)
        if isinstance(test, dict):
            add(label, "info",
                f"significant={test.get('significant')}, "
                f"confidence={_fmt(test.get('confidence'), 4)}")

    systematics = fit.get("cn_systematics") or {}
    if systematics:
        add("CN systematic 1σ", "info",
            f"{100.0 * (_as_float(systematics.get('relative_1sigma')) or 0.0):.1f}%",
            systematics.get("note") or "")

    if not rows:
        _draw_missing(ax, "缺少拟合稳定性检查字段")
        return Panel("Figure3g_checks_table", "Figure 3(g)",
                     "Stability checks", "拟合稳定性检查表。",
                     status="skipped", reason="no fit stability fields")
    _render_table(ax, ["Check", "Status", "Value", "Note"],
                  [[r["check"], r["status"], r["value"], r["note"]]
                   for r in rows[:14]], status_col=1, fontsize=5.8,
                  title="(g) Fit stability checks")
    shown = rows[:14]
    return Panel(
        "Figure3g_checks_table", "Figure 3(g)", "Fit stability checks",
        "拟合稳定性检查：过拟合、σ²/C3/C4 撞边界、混合权重简并、自吸收、"
        "壳层分辨与 M–M 键检验；失败项会阻断 CN 定量。",
        source="fit.json:flags / resolution_check / cn_systematics",
        columns=["check", "status", "value", "note"], rows=shown)


def panel_scan_quality(ax, ctx: dict) -> Panel:
    meta = ctx.get("meta") or {}
    source = meta.get("source") or {}
    merge = source.get("merge_info") or {}
    norm = meta.get("normalization") or {}
    rbkg = meta.get("rbkg") or {}
    rows: list[dict] = []

    def add(metric: str, value: Any, unit: str = "", origin: str = "") -> None:
        if value is None:
            return
        rows.append({"metric": metric, "value": _short(value, 56),
                     "unit": unit, "source": origin})

    add("format", source.get("format_detected") or source.get("format"))
    add("mode", source.get("mode"))
    add("n points", source.get("n_points"))
    add("n scans", source.get("n_scans"))
    add("scan RMS", source.get("scan_rms"))
    add("merge input", merge.get("n_input"))
    add("merge kept", merge.get("n_kept"))
    add("rejected scans", merge.get("rejected_indices"))
    add("edge step", _fmt(norm.get("edge_step"), 6))
    add("plateau", _fmt(norm.get("plateau"), 6))
    add("plateau ok", norm.get("plateau_ok"))
    add("strict", norm.get("strict"))
    add("Rbkg selected", _fmt(rbkg.get("selected"), 6), "Å")
    add("Rbkg confidence", rbkg.get("confidence"))
    if not rows:
        _draw_missing(ax, "缺少数据源/归一化元数据")
        return Panel("Figure3h_scan_quality", "Figure 3(h)",
                     "Scan & data quality", "扫描与数据质量表。",
                     status="skipped", reason="no source metadata")
    shown = rows[:14]
    _render_table(ax, ["Metric", "Value", "Unit", "Source"],
                  [[r["metric"], r["value"], r["unit"], r["source"]]
                   for r in shown], fontsize=5.9,
                  title="(h) Scan & normalization quality")
    return Panel(
        "Figure3h_scan_quality", "Figure 3(h)", "Scan & data quality",
        "扫描与数据质量：格式、扫描数、合并/剔除、重复性、边跃高度、平台"
        "与 Rbkg 选择；这些是后续解释 N 的系统边界。",
        source="meta.json:source / normalization / rbkg",
        columns=["metric", "value", "unit", "source"], rows=shown)


# ---------------------------------------------------------------
# 图 4：数据特征与输入清单
# ---------------------------------------------------------------
def panel_xanes_derivative(ax, ctx: dict) -> Panel:
    cols = _xanes_columns(ctx)
    if cols is None or cols.get("norm") is None:
        _draw_missing(ax, "缺少归一化 XANES")
        return Panel("Figure4a_derivative", "Figure 4(a)",
                     "XANES derivative", "dμ/dE。", status="skipped",
                     reason="no normalized XANES")
    energy = cols["energy"]
    norm = cols["norm"]
    if len(energy) < 3:
        _draw_missing(ax, "XANES 点数不足")
        return Panel("Figure4a_derivative", "Figure 4(a)",
                     "XANES derivative", "dμ/dE。", status="skipped",
                     reason="too few XANES points")
    deriv = np.gradient(norm, energy)
    ax.plot(energy, deriv, color="#6a1b9a", lw=1.1)
    e0 = _meta_e0(ctx.get("meta") or {})
    if e0 is not None:
        ax.axvline(e0, color="#c62828", ls="--", lw=0.9)
        mask = (energy >= e0 - 20.0) & (energy <= e0 + 120.0)
        if np.any(mask):
            idx = int(np.nanargmax(deriv[mask]))
            ax.scatter([energy[mask][idx]], [deriv[mask][idx]], s=26,
                       color="#ef6c00", zorder=4)
            ax.annotate(f"max dμ/dE={energy[mask][idx]:.1f} eV",
                        (energy[mask][idx], deriv[mask][idx]),
                        textcoords="offset points", xytext=(5, 5),
                        fontsize=6.5, color="#ef6c00")
        lo = max(float(np.min(energy)), e0 - 80.0)
        hi = min(float(np.max(energy)), e0 + 200.0)
        if hi > lo:
            ax.set_xlim(lo, hi)
    ax.axhline(0, color="#bdbdbd", lw=0.6)
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("dμ/dE")
    ax.set_title("(a) XANES derivative", loc="left", fontweight="bold")
    ax.grid(alpha=0.18)
    rows = [{"energy_eV": float(e), "dmu_dE": float(d),
             "mu_norm": float(v)}
            for e, d, v in zip(energy, deriv, norm)
            if math.isfinite(float(d))]
    return Panel(
        "Figure4a_derivative", "Figure 4(a)", "XANES derivative",
        "dμ/dE：吸收边一阶导数用于定位 E0 和白线；峰形变化可提示氧化态"
        "或配位环境改变。",
        source=_rel(ctx, cols["source"]),
        columns=["energy_eV", "dmu_dE", "mu_norm"], rows=rows)


def panel_qspace(ax, ctx: dict) -> Panel:
    data = ctx.get("chiq")
    if not data:
        _draw_missing(ax, "缺少 *_chiq.csv（q 空间）")
        return Panel("Figure4b_qspace", "Figure 4(b)", "q-space",
                     "q 空间。", status="skipped", reason="no q-space data")
    q = data["q"]
    ax.plot(q, data["mag"], color="#00695c", lw=1.2, label=r"$|\chi(q)|$")
    ax.plot(q, data["re"], color="#90a4ae", lw=0.7, ls="--", label="Re")
    ax.plot(q, data["im"], color="#ffb74d", lw=0.7, ls=":", label="Im")
    ax.set_xlabel(r"$q$ (Å$^{-1}$)")
    ax.set_ylabel(r"$\chi(q)$")
    ax.set_title("(b) q-space", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.5)
    ax.grid(alpha=0.18)
    rows = [{"q_inv_A": float(qq), "chiq_re": float(re_),
             "chiq_im": float(im_), "chiq_mag": float(mag_)}
            for qq, re_, im_, mag_ in zip(q, data["re"], data["im"],
                                          data["mag"])
            if math.isfinite(float(mag_))]
    return Panel(
        "Figure4b_qspace", "Figure 4(b)", "q-space χ(q)",
        "q 空间：与 k 空间互补，可检查高频噪声、截断振铃和数据一致性；"
        "q 空间异常往往对应背景或归一化问题。",
        source=_rel(ctx, data["source"]),
        columns=["q_inv_A", "chiq_re", "chiq_im", "chiq_mag"], rows=rows)


def _ft_peak_rows(ctx: dict) -> list[dict]:
    data = ctx.get("ft")
    if not data:
        return []
    r = np.asarray(data["r"], dtype=float)
    mag = np.asarray(data["mag"], dtype=float)
    if len(r) < 5 or not np.any(np.isfinite(mag)):
        return []
    peak_idx = []
    for i in range(1, len(r) - 1):
        if not (0.8 <= r[i] <= 6.0):
            continue
        if not (math.isfinite(float(mag[i])) and
                math.isfinite(float(mag[i - 1])) and
                math.isfinite(float(mag[i + 1]))):
            continue
        if mag[i] >= mag[i - 1] and mag[i] > mag[i + 1] and mag[i] > 0:
            peak_idx.append(i)
    if not peak_idx:
        return []
    global_max = max(float(mag[i]) for i in peak_idx)
    kept: list[int] = []
    for i in sorted(peak_idx, key=lambda k: -float(mag[k])):
        if float(mag[i]) < 0.05 * global_max:
            continue
        if all(abs(float(r[i] - r[j])) > 0.15 for j in kept):
            kept.append(i)
        if len(kept) >= 8:
            break
    kept.sort(key=lambda k: float(r[k]))
    shells = [s for s in ((ctx.get("fit") or {}).get("shells") or [])
              if isinstance(s, dict) and _as_float(s.get("reff")) is not None]
    amplitude_rank = {idx: rank for rank, idx in enumerate(
        sorted(kept, key=lambda k: -float(mag[k])), 1)}
    rows = []
    for rank, i in enumerate(kept, 1):
        assignment = ""
        for shell in shells:
            if abs(float(r[i]) - float(shell["reff"])) <= 0.25:
                assignment = str(shell.get("label") or "")
                break
        rows.append({
            "rank": rank,
            "r_A": float(r[i]),
            "amplitude": float(mag[i]),
            "relative": float(mag[i]) / global_max if global_max else None,
            "amplitude_rank": amplitude_rank.get(i),
            "assignment": assignment,
        })
    return rows


def panel_ft_peaks(ax, ctx: dict) -> Panel:
    data = ctx.get("ft")
    rows = _ft_peak_rows(ctx)
    if not data or not rows:
        _draw_missing(ax, "缺少可识别 FT 峰")
        return Panel("Figure4c_ft_peaks", "Figure 4(c)", "FT peaks",
                     "FT 峰识别。", status="skipped", reason="no FT peaks")
    r = data["r"]
    mag = data["mag"]
    ax.plot(r, mag, color="#4a148c", lw=1.1)
    for item in rows:
        ax.scatter([item["r_A"]], [item["amplitude"]], s=24,
                   color="#ef6c00", zorder=4)
        ax.annotate(f"{item['rank']}", (item["r_A"], item["amplitude"]),
                    textcoords="offset points", xytext=(4, 4),
                    fontsize=6.5, color="#ef6c00")
    ax.set_xlim(0, min(float(np.max(r)), 6.0))
    ax.set_xlabel(r"$R$ (Å)")
    ax.set_ylabel(r"$|\chi(R)|$ (Å$^{-3}$)")
    ax.set_title("(c) FT peaks (numbered by R; amplitude rank in CSV)",
                 loc="left", fontweight="bold")
    ax.grid(alpha=0.18)
    return Panel(
        "Figure4c_ft_peaks", "Figure 4(c)", "FT peak detection",
        "FT 峰识别：按 R 编号的前 8 个峰，amplitude_rank 为幅值排序；"
        "assignment 为最近拟合壳层；未分配峰提示模型缺壳层或多位点结构。",
        source=_rel(ctx, data["source"]),
        columns=["rank", "r_A", "amplitude", "relative", "amplitude_rank",
                 "assignment"], rows=rows)


def _cif_shell_rows(ctx: dict) -> list[dict]:
    fit = ctx.get("fit") or {}
    model = (_dig(fit, "formal_cn", "model")
             or _dig(ctx.get("meta"), "formal_cn", "model") or {})
    shells = [s for s in (model.get("n0_shells") or [])
              if isinstance(s, dict)]
    rows = []
    source = "fit.json:formal_cn.model.n0_shells"
    if shells:
        for shell in shells:
            members = shell.get("n0_members") or []
            rows.append({
                "element": shell.get("el"),
                "r_A": _as_float(shell.get("r")),
                "n0": _as_float(shell.get("n0", shell.get("n"))),
                "n_members": len(members),
                "sub_rayleigh_grouped": shell.get("sub_rayleigh_grouped"),
                "source": source,
            })
    else:
        for shell in (fit.get("shells") or []):
            if not isinstance(shell, dict):
                continue
            rows.append({
                "element": shell.get("scatterer"),
                "r_A": _as_float(shell.get("reff")),
                "n0": _as_float(shell.get("degen")),
                "n_members": shell.get("n_members"),
                "sub_rayleigh_grouped": None,
                "source": "fit.json:shells",
            })
    return [r for r in rows
            if r["r_A"] is not None and r["n0"] is not None]


def panel_cif_shells(ax, ctx: dict) -> Panel:
    rows = _cif_shell_rows(ctx)
    if not rows:
        _draw_missing(ax, "缺少 CIF/拟合壳层 N0")
        return Panel("Figure4d_cif_shells", "Figure 4(d)", "Shell N0 vs R",
                     "CIF 配位壳层。", status="skipped",
                     reason="no CIF/fits shells")
    colors = {"O": "#c62828", "N": "#1565c0", "C": "#2e7d32",
              "Ru": "#6a1b9a", "Fe": "#ef6c00", "Ni": "#00838f",
              "Co": "#5d4037", "Cu": "#ad1457"}
    for row in rows[:18]:
        color = colors.get(str(row["element"]), "#546e7a")
        ax.hlines(row["r_A"], 0, row["n0"], color=color, lw=1.4,
                  alpha=0.85)
        ax.scatter([row["n0"]], [row["r_A"]], s=22, color=color, zorder=3)
        ax.text(row["n0"] + 0.15, row["r_A"],
                f"{row['element']}{row['n0']:g}", fontsize=6.0,
                va="center", color="#455a64",
                bbox={"facecolor": "white", "alpha": 0.65, "lw": 0})
    ax.set_xlabel("Crystallographic N$_0$")
    ax.set_ylabel(r"$R$ (Å)")
    ax.set_title("(d) CIF coordination shells", loc="left",
                 fontweight="bold")
    ax.grid(alpha=0.18)
    return Panel(
        "Figure4d_cif_shells", "Figure 4(d)", "CIF coordination shells",
        "CIF 配位壳层：每个壳层的晶体学 N₀ 与 R；用于核对拟合 N/N₀、"
        "判断模型是否覆盖实验峰，以及识别亚 Rayleigh 合并。",
        source=rows[0].get("source") or "fit.json",
        columns=["element", "r_A", "n0", "n_members",
                 "sub_rayleigh_grouped", "source"], rows=rows)


def _feature_rows(ctx: dict) -> list[dict]:
    meta = ctx.get("meta") or {}
    fit = ctx.get("fit") or {}
    obj = _fit_obj(ctx)
    cols = _xanes_columns(ctx)
    rows: list[dict] = []

    def add(feature: str, value: Any, unit: str = "", source: str = "") -> None:
        if value is None:
            return
        rows.append({"feature": feature, "value": _fmt(value, 6),
                     "unit": unit, "source": source})

    e0 = _meta_e0(meta)
    add("E0", e0, "eV", "meta:e0")
    norm = meta.get("normalization") or {}
    add("edge step", norm.get("edge_step"), "", "meta:normalization")
    add("plateau", norm.get("plateau"), "", "meta:normalization")
    if cols is not None and e0 is not None and cols.get("norm") is not None:
        energy = cols["energy"]
        values = cols["norm"]
        mask = (energy >= e0 + 5.0) & (energy <= e0 + 100.0)
        if np.any(mask):
            idx = int(np.nanargmax(values[mask]))
            add("white line E", float(energy[mask][idx]), "eV",
                "xanes")
            add("white line height", float(values[mask][idx]), "",
                "xanes")
        pre_mask = (energy >= e0 - 40.0) & (energy <= e0 - 3.0)
        if np.any(pre_mask):
            i_min = int(np.nanargmin(values[pre_mask]))
            add("pre-edge minimum E", float(energy[pre_mask][i_min]),
                "eV", "xanes")
    ft_info = meta.get("ft") or {}
    add("first FT peak R", ft_info.get("r_first_peak"), "Å", "meta:ft")
    peaks = _ft_peak_rows(ctx)
    if peaks:
        add("strongest FT peak R", peaks[0]["r_A"], "Å", "Figure4c")
        add("strongest FT peak amp", peaks[0]["amplitude"], "", "Figure4c")
    wt_info = meta.get("wt") or {}
    add("WT peak R", wt_info.get("r_peak") or wt_info.get("r_first_peak"),
        "Å", "meta:wt")
    add("WT peak k", wt_info.get("k_peak") or wt_info.get("k_at_peak"),
        "Å^-1", "meta:wt")
    add("WT sigma", wt_info.get("sigma"), "Å", "meta:wt")
    rbkg = meta.get("rbkg") or {}
    add("Rbkg", rbkg.get("selected"), "Å", "meta:rbkg")
    add("Rbkg confidence", rbkg.get("confidence"), "", "meta:rbkg")
    add("S0²", fit.get("s02", obj.get("s02")), "", "fit.json")
    add("ΔE0", obj.get("delta_e0"), "eV", "fit.json:fit")
    add("R-factor", obj.get("r_factor"), "", "fit.json:fit")
    add("CN status", fit.get("cn_status"), "", "fit.json")
    formal = _formal_block(ctx)
    add("Formal CN status", formal.get("status"), "", "formal_cn")
    structure = ctx.get("structure") or {}
    add("Structure status", structure.get("status"), "", "structure_id")
    add("N independent", obj.get("n_independent"), "", "fit.json:fit")
    add("N variables", obj.get("n_variables"), "", "fit.json:fit")
    add("strongest FT peak amp", peaks[0]["amplitude"] if peaks else None,
        "", "Figure4c")
    add("WT sigma", wt_info.get("sigma"), "Å", "meta:wt")
    add("Structure coverage", structure.get("coverage_ok"), "",
        "structure_id")
    conf = structure.get("confidence") or {}
    add("Structure best score", conf.get("best_score"), "", "structure_id")
    add("Structure margin", conf.get("margin"), "", "structure_id")
    return rows


def panel_features_table(ax, ctx: dict) -> Panel:
    rows = _feature_rows(ctx)
    if not rows:
        _draw_missing(ax, "缺少特征字段")
        return Panel("Figure4e_features_table", "Figure 4(e)",
                     "Feature summary", "特征总表。", status="skipped",
                     reason="no feature fields")
    shown = rows[:16]
    _render_table(ax, ["Feature", "Value", "Unit", "Source"],
                  [[r["feature"], r["value"], r["unit"], r["source"]]
                   for r in shown], fontsize=5.7,
                  title="(e) Feature summary")
    return Panel(
        "Figure4e_features_table", "Figure 4(e)", "Feature summary",
        "特征总表：E0、边跃、平台、白线、FT/WT 峰、Rbkg、S0²、Rf、CN 与"
        "结构识别状态；一张表回答这个结果的关键数值从哪来。",
        source="meta / fit / structure_id / XANES",
        columns=["feature", "value", "unit", "source"], rows=shown)


def panel_input_inventory(ax, ctx: dict) -> Panel:
    items = [
        ("meta", "参数/门禁台账", ctx.get("meta_path")),
        ("fit", "拟合结果", ctx.get("fit_path")),
        ("structure_id", "结构识别", ctx.get("structure_path")),
        ("cn_table", "正式 CN 表", ctx.get("cn_path")),
        ("fit_table", "拟合壳层表", ctx.get("fit_table_path")),
        ("xanes", "XANES/μ(E)", (ctx.get("xanes") or {}).get("source")),
        ("chik", "χ(k)", (ctx.get("chik") or {}).get("source")),
        ("chir", "χ(R)", (ctx.get("ft") or {}).get("source")),
        ("chiq", "χ(q)", (ctx.get("chiq") or {}).get("source")),
        ("wavelet", "小波", (ctx.get("wt") or {}).get("source")),
    ]
    rows = []
    for role, label, path in items:
        if path is None:
            continue
        path = Path(path)
        if not path.is_file():
            continue
        rows.append({
            "file": path.name,
            "role": label,
            "rows": _file_rows(path),
            "size_bytes": path.stat().st_size,
            "sha256": _file_sha(path)[:16],
            "path": str(path),
        })
    if not rows:
        _draw_missing(ax, "输出目录中没有可清点的输入文件")
        return Panel("Figure4f_inventory", "Figure 4(f)",
                     "Input inventory", "输入清单。", status="skipped",
                     reason="no input files")
    _render_table(ax, ["File", "Role", "Rows", "Size (KB)", "SHA256"],
                  [[r["file"], r["role"],
                    "" if r["rows"] is None else r["rows"],
                    f"{r['size_bytes'] / 1024:.1f}", r["sha256"]]
                   for r in rows], fontsize=5.8,
                  title="(f) Input inventory")
    return Panel(
        "Figure4f_inventory", "Figure 4(f)", "Input inventory",
        "输入清单：图包实际读取的文件、行数、大小与 SHA256 前缀；用于"
        "复现、审稿回复和判断结果是否基于同一批数据。",
        source="output directory",
        columns=["file", "role", "rows", "size_bytes", "sha256", "path"],
        rows=rows)


# ---------------------------------------------------------------
# 补充表格图
# ---------------------------------------------------------------
def _cn_table_rows(ctx: dict) -> list[dict]:
    rows = []
    for row in (ctx.get("cn_rows") or []):
        rows.append({
            "phase_id": row.get("phase_id"),
            "model_level": row.get("model_level"),
            "absorber_site": row.get("absorber_site"),
            "shell": row.get("shell_label") or row.get("shell_group"),
            "n0": row.get("n0"),
            "n_fit": row.get("n_fit"),
            "n_ratio": row.get("n_ratio"),
            "delta_r": row.get("delta_r"),
            "sigma2": row.get("sigma2"),
            "rf": row.get("rf"),
            "gate_status": row.get("gate_status"),
            "reject_reason": _short(row.get("reject_reason"), 46),
        })
    return rows


def supplementary_cn_table(ax, ctx: dict) -> Panel:
    rows = _cn_table_rows(ctx)
    if not rows:
        _draw_missing(ax, "缺少 cn_table.csv")
        return Panel("Figure2g_cn_table", "Table S1", "Formal CN table",
                     "正式 CN 表。", status="skipped", reason="no cn_table.csv")
    table_rows = [[_short(r["phase_id"], 18), _short(r["model_level"], 12),
                   _short(r["absorber_site"], 10), _short(r["shell"], 16),
                   _fmt(r["n0"], 4), _fmt(r["n_fit"], 4),
                   _fmt(r["n_ratio"], 4), _fmt(r["delta_r"], 4),
                   _fmt(r["sigma2"], 4), _fmt(r["rf"], 4),
                   str(r["gate_status"] or ""), _short(r["reject_reason"], 30)]
                  for r in rows]
    _render_table(ax, ["Phase", "Level", "Site", "Shell", "N0", "N_fit",
                       "N/N0", "dR", "sigma2", "Rf", "Gate", "Reject"],
                  table_rows, status_col=10, fontsize=5.6,
                  title="Table S1. Formal coordination-number table")
    return Panel(
        "Figure2g_cn_table", "Table S1", "Formal CN table",
        "正式 CN 表：相来源、吸收位、N₀、N_fit、N/N₀、ΔR、σ²、Rf、"
        "门禁状态与拒绝原因；表格直接出图，CSV 保留完整列。",
        source="cn_table.csv",
        columns=["phase_id", "model_level", "absorber_site", "shell", "n0",
                 "n_fit", "n_ratio", "delta_r", "sigma2", "rf",
                 "gate_status", "reject_reason"], rows=rows)


def supplementary_fit_table(ax, ctx: dict) -> Panel:
    rows = _fit_rows(ctx)
    if not rows:
        _draw_missing(ax, "缺少拟合结果")
        return Panel("Figure2h_fit_table", "Table S2", "EXAFS fit table",
                     "EXAFS 拟合表。", status="skipped", reason="no fit rows")
    table_rows = [[_short(r["shell"], 18), _fmt(r["n"], 4),
                   _fmt(r["n_err"], 3), _fmt(r["n0"], 4),
                   _fmt(r["n_ratio"], 3), _fmt(r["r"], 4),
                   _fmt(r["r_err"], 4), _fmt(r["sigma2"], 4),
                   _fmt(r["rf"], 4), _short(r["status"], 16)]
                  for r in rows[:12]]
    _render_table(ax, ["Shell", "N", "N_err", "N0", "N/N0", "R",
                       "R_err", "sigma2", "Rf", "Status"],
                  table_rows, status_col=9, fontsize=5.8,
                  title="Table S2. EXAFS fit parameters")
    return Panel(
        "Figure2h_fit_table", "Table S2", "EXAFS fit table",
        "EXAFS 拟合表：每个壳层的 N、误差、N₀、R、σ²、Rf 与可靠性状态；"
        "表格直接出图，CSV 保留完整列。",
        source="fit.json / cn_table.csv",
        columns=["shell", "n", "n_err", "n0", "n_ratio", "r", "r_err",
                 "sigma2", "rf", "mode", "status", "note", "source"],
        rows=rows)


# ---------------------------------------------------------------
# 保存 / 清单
# ---------------------------------------------------------------
def _rel(ctx: dict, path: Path | None) -> str:
    if path is None:
        return ""
    try:
        return str(Path(path).resolve().relative_to(ctx["outdir"].resolve()))
    except ValueError:
        return str(path)


def _csv_cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.10g}"
    if isinstance(value, (dict, list, tuple)):
        try:
            return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
        except Exception:  # noqa: BLE001
            return str(value)
    return value


def _write_rows_csv(path: Path, columns: list[str],
                    rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)
        for row in rows:
            writer.writerow([_csv_cell(row.get(col)) for col in columns])


def _write_panel_csv(path: Path, panel: Panel) -> None:
    _write_rows_csv(path, panel.columns, panel.rows)


def _save_standalone(ctx: dict, panel_fn: Callable, base: Path,
                     formats: tuple[str, ...], dpi: int,
                     figsize: tuple[float, float],
                     panel: Panel | None = None) -> Panel:
    fig, ax = plt.subplots(figsize=figsize)
    try:
        result = panel_fn(ax, ctx)
        fig.tight_layout()
        for fmt in formats:
            fig.savefig(base.with_suffix(f".{fmt}"), dpi=dpi)
    finally:
        plt.close(fig)
    return result if panel is None else panel


def _master_figure(ctx: dict, panels: list[Callable], *, title: str,
                   figsize: tuple[float, float], nrow: int, ncol: int,
                   panel_ids: list[Panel],
                   hspace: float = 0.32,
                   wspace: float = 0.24) -> tuple[plt.Figure, list[Panel]]:
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(
        nrow, ncol, left=0.06, right=0.985, top=0.935, bottom=0.045,
        hspace=hspace, wspace=wspace)
    out: list[Panel] = []
    for idx, fn in enumerate(panels):
        ax = fig.add_subplot(gs[idx // ncol, idx % ncol])
        out.append(fn(ax, ctx))
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    fig.text(0.5, 0.008,
             "Generated by xafs_figure_pack.py — each panel has a matching CSV; "
             "see figure_manifest.json",
             ha="center", fontsize=7, color="#607d8b")
    return fig, out


def _save_master(fig: plt.Figure, base: Path, formats: tuple[str, ...],
                 dpi: int) -> list[str]:
    paths = []
    for fmt in formats:
        path = base.with_suffix(f".{fmt}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        paths.append(path.name)
    plt.close(fig)
    return paths


def _emit_figure(ctx: dict, fns: list[Callable], *, base: str, title: str,
                 figsize: tuple[float, float], nrow: int, ncol: int,
                 formats: tuple[str, ...], dpi: int, figures_dir: Path,
                 standalone_size: tuple[float, float],
                 panels_manifest: list[dict],
                 notes: list[str],
                 hspace: float = 0.32,
                 wspace: float = 0.24,
                 standalone_panels: bool = False) -> dict:
    """落盘一张总图、逐面板 CSV（归入 panels/）及可选的独立面板图。

    总图已包含全部子面板；默认不再重复渲染散装单面板 PNG，
    避免耗时翻倍与 figures/ 根目录被淹没。需要单面板图时显式打开
    ``standalone_panels``，输出到 ``figures/panels/``。
    """
    fig, panels = _master_figure(
        ctx, fns, title=title, figsize=figsize, nrow=nrow, ncol=ncol,
        panel_ids=[], hspace=hspace, wspace=wspace)
    paths = _save_master(fig, figures_dir / base, formats, dpi)
    image_ext = formats[0]
    panels_dir = figures_dir / "panels"
    panels_dir.mkdir(parents=True, exist_ok=True)
    for fn, panel in zip(fns, panels):
        image_rel = None
        if panel.status == "ok":
            if standalone_panels:
                _save_standalone(ctx, fn, panels_dir / panel.id, formats, dpi,
                                 standalone_size, panel)
                image_rel = f"panels/{panel.id}.{image_ext}"
            else:
                # 面板图即所属总图：不重复渲染，仍可被清单定位。
                image_rel = f"{base}.{image_ext}"
        csv_name = _panel_csv_name(panel)
        csv_rel = None
        if csv_name:
            _write_panel_csv(panels_dir / csv_name, panel)
            csv_rel = f"panels/{csv_name}"
        panels_manifest.append(_manifest_panel(
            ctx, panel, figures_dir,
            image=image_rel, csv_name=csv_rel))
        if panel.status != "ok":
            notes.append(f"{panel.label} skipped: {panel.reason}")
    return {"paths": paths, "panels": panels}


# ---------------------------------------------------------------
# 数据包
# ---------------------------------------------------------------
def _write_data_bundle(ctx: dict, figures_dir: Path,
                       manifest: dict) -> list[str]:
    """跨面板汇总写入 ``figures/data/``。

    源 JSON/CSV 不在此复制：父目录已有唯一权威副本，逐面板 CSV 在
    ``figures/panels/``。本目录只保留无法从单个面板还原的汇总表。
    """
    data_dir = figures_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []

    fit = ctx.get("fit") or {}
    meta = ctx.get("meta") or {}
    obj = _fit_obj(ctx)
    cond = fit.get("conditional_cn") or {}

    # 1) 全路径表：主拟合 + fixed/conditional/free N + L1
    path_rows = _fit_paths(ctx) + _fit_variant_rows(ctx)
    path_cols = ["variant", "label", "n", "n_err", "r", "r_err", "reff",
                 "delta_r", "sigma2", "sigma2_err", "n0", "n_ratio", "rf",
                 "note"]
    _write_rows_csv(data_dir / "all_fit_paths.csv", path_cols, path_rows)
    written.append("all_fit_paths.csv")

    # 2) 拟合变体汇总：主拟合 + 固定/条件/自由 N
    variant_rows = []
    if obj:
        variant_rows.append({
            "variant": "fit",
            "r_factor": _as_float(obj.get("r_factor")),
            "n_variables": obj.get("n_variables"),
            "n_independent": _as_float(obj.get("n_independent")),
            "delta_e0": _as_float(obj.get("delta_e0")),
            "overfit": obj.get("overfit"),
            "n_ratio": _as_float(cond.get("n_ratio")),
            "n_status": cond.get("n_status"),
            "consistency_passed": cond.get("consistency_passed"),
            "hamilton_F": _as_float(_dig(cond, "hamilton_f_test", "F")),
            "hamilton_p": _as_float(_dig(cond, "hamilton_f_test", "p_value")),
            "significant": _dig(cond, "hamilton_f_test", "significant"),
        })
    for key, label in (("fixed_n_fit", "fixed_N"),
                       ("conditional_n_fit", "conditional_N"),
                       ("free_n_fit", "free_N")):
        item = cond.get(key) or {}
        if not item:
            continue
        variant_rows.append({
            "variant": label,
            "r_factor": _as_float(item.get("r_factor")),
            "n_variables": item.get("n_variables"),
            "n_independent": _as_float(item.get("n_independent")),
            "delta_e0": _as_float(item.get("delta_e0")),
            "overfit": item.get("overfit"),
            "n_ratio": _as_float(cond.get("n_ratio")),
            "n_status": cond.get("n_status"),
            "consistency_passed": cond.get("consistency_passed"),
            "hamilton_F": _as_float(_dig(cond, "hamilton_f_test", "F")),
            "hamilton_p": _as_float(_dig(cond, "hamilton_f_test", "p_value")),
            "significant": _dig(cond, "hamilton_f_test", "significant"),
        })
    if variant_rows:
        _write_rows_csv(
            data_dir / "fit_variants.csv",
            ["variant", "r_factor", "n_variables", "n_independent",
             "delta_e0", "overfit", "n_ratio", "n_status",
             "consistency_passed", "hamilton_F", "hamilton_p",
             "significant"], variant_rows)
        written.append("fit_variants.csv")

    # 3) 参数全表
    params = meta.get("params") or {}
    effective = params.get("effective") or {}
    sources = params.get("sources") or {}
    units = params.get("units") or {}
    stages = params.get("stages") or {}
    formal = set(params.get("affects_formal_cn") or [])
    param_rows = [{
        "parameter": key,
        "value": _csv_cell(effective.get(key)),
        "source": sources.get(key) or "",
        "unit": units.get(key) or "",
        "stage": stages.get(key) or "",
        "affects_formal_cn": key in formal,
    } for key in sorted(effective)]
    if param_rows:
        _write_rows_csv(
            data_dir / "all_parameters.csv",
            ["parameter", "value", "source", "unit", "stage",
             "affects_formal_cn"], param_rows)
        written.append("all_parameters.csv")

    # 4) 门禁全表（保留完整 observed/threshold JSON；面板表会截断）
    gates = meta.get("gates") or {}
    gate_rows = []
    for name, payload in gates.items():
        if name == "summary" or not isinstance(payload, dict):
            continue
        gate_rows.append({
            "gate": name,
            "status": payload.get("status"),
            "action": payload.get("action"),
            "threshold": _csv_cell(payload.get("threshold")),
            "observed": _csv_cell(payload.get("observed")),
            "reject_reason": payload.get("reject_reason") or "",
            "location": payload.get("location") or "",
        })
    if gate_rows:
        _write_rows_csv(
            data_dir / "gate_checks.csv",
            ["gate", "status", "action", "threshold", "observed",
             "reject_reason", "location"], gate_rows)
        written.append("gate_checks.csv")

    # 5) 决策全表（不截断 detail；面板表只画前 12 条）
    decision_rows = []
    for item in (meta.get("decisions") or []):
        if not isinstance(item, dict):
            continue
        decision_rows.append({
            "step": item.get("step"),
            "decision": item.get("decision"),
            "status": item.get("status"),
            "detail": item.get("detail"),
        })
    if decision_rows:
        _write_rows_csv(
            data_dir / "decisions.csv",
            ["step", "decision", "status", "detail"], decision_rows)
        written.append("decisions.csv")

    # 6) 面板索引
    index_rows = [{
        "label": p.get("label"),
        "id": p.get("id"),
        "title": p.get("title"),
        "status": p.get("status"),
        "image": p.get("image") or "",
        "csv": p.get("csv") or "",
        "n_rows": p.get("n_rows"),
        "source": p.get("source") or "",
        "caption": p.get("caption") or "",
    } for p in manifest.get("panels", [])]
    _write_rows_csv(
        data_dir / "panel_index.csv",
        ["label", "id", "title", "status", "image", "csv", "n_rows",
         "source", "caption"], index_rows)
    written.append("panel_index.csv")

    # 7) 汇总 JSON
    gates = meta.get("gates") or {}
    summary_counts = _summary_counts(gates)
    structure = ctx.get("structure") or {}
    formal_block = _formal_block(ctx)
    cif_rows = _cif_shell_rows(ctx)
    summary = {
        "schema": "xafs-figure-pack-data/v1",
        "generated_at": manifest.get("generated_at"),
        "sample": manifest.get("sample"),
        "stem": ctx["stem"],
        "absorber": manifest.get("absorber"),
        "edge": manifest.get("edge"),
        "e0_eV": _meta_e0(meta),
        "fit": {
            "mode": fit.get("mode"),
            "status": fit.get("status"),
            "r_factor": _as_float(obj.get("r_factor")),
            "delta_e0": _as_float(obj.get("delta_e0")),
            "s02": _as_float(fit.get("s02")),
            "s02_source": fit.get("s02_source"),
            "cn_status": fit.get("cn_status"),
            "cn_reliable": fit.get("cn_reliable"),
            "n_paths": len(obj.get("paths") or []),
            "n_independent": _as_float(obj.get("n_independent")),
            "n_variables": obj.get("n_variables"),
        },
        "formal_cn": {
            "status": formal_block.get("status"),
            "reject_reason": formal_block.get("reject_reason"),
        },
        "structure": {
            "status": structure.get("status"),
            "best_id": (((structure.get("ranked") or [{}])[0] or {}).get("id")
                        if structure.get("ranked") else None),
            "best_score": _dig(structure, "confidence", "best_score"),
            "margin": _dig(structure, "confidence", "margin"),
            "coverage_ok": structure.get("coverage_ok"),
        },
        "counts": {
            "gates": {k: len(v) for k, v in summary_counts.items()},
            "decisions": len(meta.get("decisions") or []),
            "parameters": len(effective),
            "fit_paths": len(path_rows),
            "cif_shells": len(cif_rows),
            "panels_ok": sum(1 for p in manifest.get("panels", [])
                             if p.get("status") == "ok"),
        },
        "inputs": manifest.get("inputs"),
    }
    (data_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    written.append("summary.json")

    readme = """# 图包数据目录

本目录由 `xafs_figure_pack.py` 自动生成，对应 `figures/` 中的 Figure 1–4
与 Table S1/S2。所有文件均为只读快照，不重新拟合、不改变正式输出。

目录分工（不重复拷贝）：

- `figures/` 根目录：4 张总图 + Table S1/S2 表图 + 清单/索引。
- `figures/panels/`：每个子面板的唯一 CSV（需要时还含独立面板 PNG）。
- `figures/data/`：本目录，只放无法从单个面板还原的跨面板汇总。
- 父目录：`<stem>_meta.json` / `<stem>_fit.json` /
  `<stem>_structure_id.json` / `<stem>_cn_table.csv` 等权威原件，
  不在本目录重复复制。

## 跨面板汇总

| 文件 | 内容 |
|---|---|
| `all_fit_paths.csv` | 主拟合 + fixed/conditional/free N + L1 的全部路径参数 |
| `fit_variants.csv` | 主拟合与 fixed/conditional/free N 三套变体的 Rf、变量数、Hamilton F |
| `all_parameters.csv` | `meta.params.effective` 全参数：值、来源、单位、阶段、是否影响正式 CN |
| `gate_checks.csv` | 全部门禁：状态、动作、完整阈值/观测值、拒绝原因、代码位置 |
| `decisions.csv` | S1 → S7 → 正式 CN 的完整决策轨迹（detail 不截断） |
| `panel_index.csv` | 全部面板 ↔ 图 ↔ 面板 CSV ↔ 源文件 ↔ 图注 ↔ 状态 |
| `summary.json` | 机器可读汇总：样品、拟合、正式 CN、结构识别与计数 |

参数相关矩阵、CIF 壳层表、结构候选等逐面板数据只在
`figures/panels/` 保留一份，避免同一张表多处落盘。

## 使用建议

- 重绘某个面板：读 `panel_index.csv` 找到 `panels/` 下对应 CSV。
- 回复审稿：优先给 `summary.json` + `gate_checks.csv` +
  `panels/Figure3a_fitcorr.csv`。
- 复核拟合：`all_fit_paths.csv` 与 `fit_variants.csv` 提供路径级和变体级数据。
- 缺失数据不会伪造：对应面板在 `figure_manifest.json` 中标记
  `status=skipped` 并写明原因。
"""
    (data_dir / "README.md").write_text(readme, encoding="utf-8")
    written.append("README.md")
    return written


def _manifest_panel(ctx: dict, panel: Panel, figures_dir: Path,
                    image: str | None = None,
                    csv_name: str | None = None) -> dict:
    return {
        "label": panel.label,
        "id": panel.id,
        "title": panel.title,
        "status": panel.status,
        "reason": panel.reason,
        "image": image,
        "csv": csv_name,
        "source": panel.source,
        "columns": list(panel.columns),
        "n_rows": len(panel.rows),
        "caption": panel.caption,
    }


def _write_manifest(ctx: dict, figures_dir: Path, manifest: dict) -> Path:
    path = figures_dir / "figure_manifest.json"
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2),
                    encoding="utf-8")
    return path


def _write_markdown(ctx: dict, figures_dir: Path, manifest: dict) -> Path:
    lines = [
        f"# Figure pack — {manifest.get('sample', ctx['stem'])}",
        "",
        f"- Generated: {manifest['generated_at']}",
        f"- Source directory: `{ctx['outdir']}`",
        f"- Stem: `{ctx['stem']}`",
        f"- Absorber / edge: {manifest.get('absorber') or '—'} / "
        f"{manifest.get('edge') or '—'}",
        "",
        "## Master figures",
        "",
    ]
    for key, label in (("Figure1", "Figure 1 — structure overview"),
                       ("Figure2", "Figure 2 — confidence & process"),
                       ("Figure3", "Figure 3 — fit & data quality"),
                       ("Figure4", "Figure 4 — data features & inputs")):
        info = manifest.get("masters", {}).get(key)
        if info:
            lines.append(f"- **{label}**: `{info['image']}`")
    lines += ["", "## Panels and CSV files", "",
              "| Panel | Image | CSV | What it explains |",
              "|---|---|---|---|"]
    for item in manifest.get("panels", []):
        image = f"`{item['image']}`" if item.get("image") else "—"
        csv_name = f"`{item['csv']}`" if item.get("csv") else "—"
        caption = item.get("caption") or ""
        lines.append(f"| {item['label']} | {image} | {csv_name} | {caption} |")
    data_files = manifest.get("data_files") or []
    if data_files:
        lines += ["", "## Cross-panel data bundle", "",
                  "`data/` holds cross-panel summary tables; the authoritative",
                  "source JSON/CSV stay in the parent output directory:",
                  ""]
        for name in data_files:
            lines.append(f"- `{name}`")
    lines += ["", "## Notes", ""]
    for note in manifest.get("notes", []):
        lines.append(f"- {note}")
    path = figures_dir / "figure_pack.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def build_pack(outdir: Path | str, *, stem: str | None = None,
               formats: tuple[str, ...] = ("png",),
               dpi: int = 300, standalone_panels: bool = False,
               verbose: bool = True) -> dict:
    """生成完整图包；返回 manifest。

    ``standalone_panels=False``（默认）只画 4 张总图 + Table S1/S2，
    子面板 CSV 收口到 ``figures/panels/``；打开时额外渲染单面板 PNG。
    """
    _style()
    ctx = load_context(Path(outdir), stem)
    figures_dir = ctx["outdir"] / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    formats = tuple(dict.fromkeys(str(f).lower().lstrip(".") for f in formats))
    if not formats:
        formats = ("png",)
    meta = ctx.get("meta") or {}
    absorber = _meta_absorber(meta)
    edge = _meta_edge(meta, absorber)
    sample = _sample_name(ctx)
    notes: list[str] = []
    panels_manifest: list[dict] = []

    # ---- Figure 1 ----
    fig1 = _emit_figure(
        ctx,
        [panel_xanes, panel_exafs, panel_ft, panel_wt,
         panel_fit, panel_structure, panel_cn_gate, panel_card],
        base="Figure1_master",
        title=f"Figure 1. {sample} — {absorber} {edge}",
        figsize=(13.5, 16.5), nrow=4, ncol=2,
        formats=formats, dpi=dpi, figures_dir=figures_dir,
        standalone_size=(6.4, 4.2),
        panels_manifest=panels_manifest, notes=notes,
        standalone_panels=standalone_panels)

    # ---- Figure 2 ----
    fig2 = _emit_figure(
        ctx,
        [panel_gate_board, panel_decision_flow, panel_param_sources,
         panel_gate_table, panel_decision_table, panel_param_table],
        base="Figure2_confidence",
        title=f"Figure 2. Confidence & process — {sample}",
        figsize=(14.0, 15.0), nrow=3, ncol=2,
        formats=formats, dpi=dpi, figures_dir=figures_dir,
        standalone_size=(7.2, 4.8),
        panels_manifest=panels_manifest, notes=notes,
        standalone_panels=standalone_panels)

    # ---- Figure 3 ----
    fig3 = _emit_figure(
        ctx,
        [panel_fit_correlation, panel_fit_paths, panel_error_budget,
         panel_rbkg, panel_normalization, panel_fit_window_table,
         panel_checks_table, panel_scan_quality],
        base="Figure3_diagnostics",
        title=f"Figure 3. Fit & data quality — {sample}",
        figsize=(13.5, 16.5), nrow=4, ncol=2,
        formats=formats, dpi=dpi, figures_dir=figures_dir,
        standalone_size=(7.2, 4.8),
        panels_manifest=panels_manifest, notes=notes,
        hspace=0.55, standalone_panels=standalone_panels)

    # ---- Figure 4 ----
    fig4 = _emit_figure(
        ctx,
        [panel_xanes_derivative, panel_qspace, panel_ft_peaks,
         panel_cif_shells, panel_features_table, panel_input_inventory],
        base="Figure4_features",
        title=f"Figure 4. Data features & inputs — {sample}",
        figsize=(14.0, 13.0), nrow=3, ncol=2,
        formats=formats, dpi=dpi, figures_dir=figures_dir,
        standalone_size=(7.2, 4.8),
        panels_manifest=panels_manifest, notes=notes,
        hspace=0.40, standalone_panels=standalone_panels)

    # ---- 补充表格图（独立成图；逐面板 CSV 归入 panels/）----
    image_ext = formats[0]
    panels_dir = figures_dir / "panels"
    panels_dir.mkdir(parents=True, exist_ok=True)
    for fn, figsize in ((supplementary_cn_table, (12.0, 4.8)),
                        (supplementary_fit_table, (10.0, 4.2))):
        fig, ax = plt.subplots(figsize=figsize)
        try:
            panel = fn(ax, ctx)
            if panel.status == "ok":
                fig.tight_layout()
                for fmt in formats:
                    fig.savefig(figures_dir / f"{panel.id}.{fmt}", dpi=dpi)
        finally:
            plt.close(fig)
        if panel.status != "ok":
            panels_manifest.append(_manifest_panel(
                ctx, panel, figures_dir, image=None, csv_name=None))
            notes.append(f"{panel.label} skipped: {panel.reason}")
            continue
        csv_name = _panel_csv_name(panel)
        csv_rel = None
        if csv_name:
            _write_panel_csv(panels_dir / csv_name, panel)
            csv_rel = f"panels/{csv_name}"
        panels_manifest.append(_manifest_panel(
            ctx, panel, figures_dir,
            image=f"{panel.id}.{image_ext}", csv_name=csv_rel))

    manifest = {
        "schema": "xafs-figure-pack/v1",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "source_dir": str(ctx["outdir"]),
        "stem": ctx["stem"],
        "sample": sample,
        "absorber": absorber,
        "edge": edge,
        "masters": {
            "Figure1": {"image": fig1["paths"][0],
                        "images": fig1["paths"],
                        "label": "Figure 1 — structure overview"},
            "Figure2": {"image": fig2["paths"][0],
                        "images": fig2["paths"],
                        "label": "Figure 2 — confidence & process"},
            "Figure3": {"image": fig3["paths"][0],
                        "images": fig3["paths"],
                        "label": "Figure 3 — fit & data quality"},
            "Figure4": {"image": fig4["paths"][0],
                        "images": fig4["paths"],
                        "label": "Figure 4 — data features & inputs"},
        },
        "panels": panels_manifest,
        "csv_files": [p["csv"] for p in panels_manifest if p.get("csv")],
        "inputs": {
            "meta": _rel(ctx, ctx.get("meta_path")),
            "fit": _rel(ctx, ctx.get("fit_path")),
            "structure_id": _rel(ctx, ctx.get("structure_path")),
            "cn_table": _rel(ctx, ctx.get("cn_path")),
            "fit_table": _rel(ctx, ctx.get("fit_table_path")),
            "xanes": _rel(ctx, (ctx.get("xanes") or {}).get("source")),
            "chik": _rel(ctx, (ctx.get("chik") or {}).get("source")),
            "ft": _rel(ctx, (ctx.get("ft") or {}).get("source")),
            "chiq": _rel(ctx, (ctx.get("chiq") or {}).get("source")),
            "wavelet": _rel(ctx, (ctx.get("wt") or {}).get("source")),
        },
        "notes": notes,
    }
    data_files = _write_data_bundle(ctx, figures_dir, manifest)
    manifest["data_files"] = [f"data/{name}" for name in data_files]
    manifest_path = _write_manifest(ctx, figures_dir, manifest)
    md_path = _write_markdown(ctx, figures_dir, manifest)
    if verbose:
        print(f"图包 → {figures_dir}")
        print(f"  总图 1: {figures_dir / fig1['paths'][0]}")
        print(f"  总图 2: {figures_dir / fig2['paths'][0]}")
        print(f"  总图 3: {figures_dir / fig3['paths'][0]}")
        print(f"  总图 4: {figures_dir / fig4['paths'][0]}")
        print(f"  面板 CSV: {figures_dir / 'panels'}（{len(manifest['csv_files'])} 个）")
        print(f"  数据包: {figures_dir / 'data'}（{len(data_files)} 个文件）")
        print(f"  清单: {manifest_path}")
        print(f"  索引: {md_path}")
        if notes:
            print(f"  跳过/提示: {len(notes)} 条（见 manifest.notes）")
    return manifest


# ---------------------------------------------------------------
# CLI
# ---------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="催化剂 XAFS 结果图包：总图 + 逐面板 CSV + "
                    "置信/过程仪表盘 + 表格图")
    parser.add_argument("outdir", help="run_xafs.py / xafs_report.py 的输出目录")
    parser.add_argument("--stem", default=None,
                        help="样品 stem；默认自动识别 *_meta.json")
    parser.add_argument("--formats", default="png",
                        help="输出格式，逗号分隔（默认 png；可 png,pdf）")
    parser.add_argument("--dpi", type=int, default=300, help="图片 DPI（默认 300）")
    parser.add_argument("--standalone-panels", action="store_true",
                        help="额外把子面板单独渲染为 panels/<id>.png"
                             "（默认关闭：只输出 4 张总图 + Table S1/S2）")
    parser.add_argument("--quiet", action="store_true", help="只打印错误")
    args = parser.parse_args(argv)
    formats = tuple(x.strip() for x in args.formats.split(",") if x.strip())
    try:
        build_pack(args.outdir, stem=args.stem, formats=formats,
                   dpi=args.dpi, standalone_panels=args.standalone_panels,
                   verbose=not args.quiet)
    except FigurePackError as exc:
        print(f"图包生成失败: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
