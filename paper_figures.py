# -*- coding: utf-8 -*-
"""paper_figures.py — 通用 XAFS 论文出图模块
================================================

目标：把一次 XAFS 分析的“论文主图 + CN 拟合表”标准化，同时**不绑定
吸收元素、边、壳层或某个具体拟合模型**。所有元素名、壳层名、标线位置
都由调用者传入；模块只负责处理、重建模型和排版。

产出的主图固定为六个面板（有拟合时）：

    (a) 归一化 XANES
    (b) k^n chi(k)
    (c) |FT(k^n chi(k))|
    (d) R 空间数据 vs 拟合
    (e) k 空间数据 vs 拟合
    (f) CN 拟合表

没有拟合时退化为一行三面板 (a)-(c)。

设计边界
--------
- **不重新拟合**：S7 的结果必须已有（``FitResult`` 或 fit JSON）。
  (d)(e) 面板只是按已有参数重建模型，用来核对，不产生新的物理判定。
- **不改 S0^2 / E0 / 拟合窗口**：这些值原样从输入读取；本模块不做
  “为了让图好看”的参数调整。
- **不猜壳层归属**：|FT| 的标线由 ``ft_guides`` 传入；R 空间拟合面板的
  标线来自 ``fit.paths``，并明确标注“未做相位修正”。
- **CN 表分级**：formal / diagnostic / upper limit / sigma^2-at-bound
  分开显示，避免把诊断级 CN 当正式配位数写进论文。

用法（库）
----------
>>> import paper_figures as pf
>>> out = pf.make_paper_figure(
...     sample="sample_mu.dat",
...     absorber="M",
...     refs=[("Reference", "ref_mu.dat")],
...     ft_guides=[pf.FTGuide("M-O", 1.66), pf.FTGuide("M-M", 2.58)],
...     fit="sample_fit.json",
... )
>>> out["figure"]
PosixPath('.../sample_paper.png')

用法（命令行）
--------------
python paper_figures.py --sample sample_mu.dat --absorber M --edge K \
    --ref "Reference=ref_mu.dat" \
    --kmin 3 --kmax 12 --kweight 3 --rbkg 1.3 \
    --guide "M-O=1.66" --guide "M-M=2.58" --guide "M-Ce=2.70..3.30" \
    --fit-json sample_fit.json \
    --out paper_out

fit JSON 的兼容性
-----------------
优先读取 ``run_xafs.py`` / ``xafs_report.py`` 的 S7 输出结构（``fit.paths``
+ ``path_files`` + ``cn_status`` + ``cn_systematics``）。也接受本模块的
``fit_result_to_json(...)`` 产物。缺失字段只影响对应面板，不会伪造数值。

结构 / CIF provenance
---------------------
fit JSON 可带顶层 ``"structure"``；没有时可用 ``structure_provenance(...)``
构造并通过 ``make_paper_figure(..., structure=...)`` 或 CLI
``--structure-json`` 补录。provenance 只影响图注、CN 表、CSV 与 fit JSON
的披露文字，不参与任何拟合或数值判定。没有使用 CIF 时用 ``cif=None``
显式记录，输出写 ``CIF: none (no CIF used)``，不会拿无关 CIF 顶替。
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import artemis          # noqa: E402
import athena           # noqa: E402
import xafs_pipeline    # noqa: E402
from absorber import get_absorber  # noqa: E402
from xafs_core import FitResult, PipelineResult  # noqa: E402
from cn_disclosure import check_fit_payload, cif_status  # noqa: E402

# 控制台默认 GBK，打印 Å / σ² / Δ 会直接崩；入口统一改成 UTF-8。
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8", errors="replace")

SERIES_COLORS = ("#c0392b", "#2c7fb8", "#4daf4a", "#984ea3",
                 "#ff7f00", "#8d6e3a", "#17becf", "#bcbd22")
FIT_COLOR = "#d62728"
DATA_COLOR = "#111111"


# ================================================================
# 1. 公共小工具
# ================================================================

@dataclass(frozen=True)
class FTGuide:
    """|FT| 面板的一根标线或一个区间。

    Parameters
    ----------
    label : str
        标线文字，例如 ``"M-O"``。
    x : float
        标线位置（Å），或区间起点。
    x_end : float | None
        给定时画双向箭头区间 ``[x, x_end]``；否则画单根竖线。
    text_y : float | None
        文字的 y 位置（与 |FT| 同单位）。None 时自动放在曲线峰值上方。
    """
    label: str
    x: float
    x_end: float | None = None
    text_y: float | None = None


def parse_guide(text: str) -> FTGuide:
    """解析 CLI 的 ``LABEL=X`` 或 ``LABEL=X..Y``。"""
    label, sep, value = text.partition("=")
    label = label.strip()
    value = value.strip()
    if not sep or not label or not value:
        raise ValueError(f"标线格式应为 LABEL=X 或 LABEL=X..Y，收到: {text!r}")
    if ".." in value:
        lo_text, hi_text = value.split("..", 1)
        lo, hi = float(lo_text), float(hi_text)
        if lo > hi:
            lo, hi = hi, lo
        return FTGuide(label, lo, hi)
    return FTGuide(label, float(value))


def _coerce_guides(guides: Iterable[Any] | None) -> list[FTGuide]:
    out: list[FTGuide] = []
    for item in guides or ():
        if isinstance(item, FTGuide):
            out.append(item)
        elif isinstance(item, str):
            out.append(parse_guide(item))
        elif isinstance(item, dict):
            out.append(FTGuide(
                label=str(item.get("label", "")),
                x=float(item["x"]),
                x_end=(None if item.get("x_end") is None
                       else float(item["x_end"])),
                text_y=(None if item.get("text_y") is None
                        else float(item["text_y"])),
            ))
        elif isinstance(item, (tuple, list)) and len(item) >= 2:
            out.append(FTGuide(str(item[0]), float(item[1]),
                               None if len(item) < 3 else float(item[2])))
        else:
            raise TypeError(f"无法识别的标线: {item!r}")
    return [g for g in out if g.label]


def _fmt_pm(value: Any, err: Any, nd: int = 2,
            fixed_text: str | None = "fixed") -> str:
    """``value ± err``；err 为 None 时按需要标注 fixed。"""
    if value is None:
        return "—"
    try:
        v = float(value)
    except (TypeError, ValueError):
        return str(value)
    if err is None:
        return f"{v:.{nd}f}" + (f" ({fixed_text})" if fixed_text else "")
    try:
        e = float(err)
    except (TypeError, ValueError):
        return f"{v:.{nd}f}"
    return f"{v:.{nd}f} ± {e:.{nd}f}"


def _fmt_float(value: Any, nd: int = 3) -> str:
    if value is None:
        return "—"
    try:
        return f"{float(value):.{nd}f}"
    except (TypeError, ValueError):
        return str(value)


def _style() -> None:
    """出版风格：无衬线、内向刻度、四边刻度、粗轴线。"""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "axes.linewidth": 1.4,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9.5,
        "font.size": 11,
        "figure.dpi": 120,
        "savefig.dpi": 300,
    })


def _panel_title(ax, letter: str, title: str) -> None:
    ax.set_title(f"({letter}) {title}", loc="left", fontsize=12,
                 fontweight="bold", pad=8)


def _tidy_axis(ax) -> None:
    ax.tick_params(which="major", length=5, width=1.2)
    ax.tick_params(which="minor", length=2.5, width=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(1.4)


def _clean_stem(text: str) -> str:
    return (str(text).replace("$", "").replace("/", "-")
            .replace("\\", "-").replace("@", "_at_").replace(" ", "_"))


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (str, Path)):
        return [str(value)]
    return [str(v) for v in value]


def _short_path(value: str | Path) -> str:
    p = Path(str(value))
    parent = p.parent.name
    return f"{parent}/{p.name}" if parent and parent not in (".", "/") else p.name


def _normalize_structure(raw: Any, base: Path) -> dict:
    """把结构/CIF provenance 统一成可显示、可落盘的字段。

    调用者可以直接给 ``fit_result_to_json(..., structure=...)`` 或
    ``structure_provenance(...)``。本函数不猜结构，只做字段别名、
    相对路径解析和 SHA256 补齐；缺失的字段保持缺失，输出时明确写
    "not stated"。
    """
    if raw is None:
        return {}
    if isinstance(raw, str):
        return {"description": raw}
    if not isinstance(raw, dict):
        return {"description": str(raw)}

    s = dict(raw)
    src = s.get("source")
    if isinstance(src, dict):
        s.setdefault("source_type", src.get("type"))
        s.setdefault("source_ref", src.get("ref"))
        s.setdefault("provenance_status", src.get("status"))

    aliases = {
        "type": "source_type",
        "ref": "source_ref",
        "structure_type": "source_type",
        "structure_ref": "source_ref",
        "cif_path": "cif",
        "cif_file": "cif",
        "source_cif": "cif",
        "source_cif_sha256": "cif_sha256",
        "motif": "phase_id",
        "phase": "phase_id",
    }
    for old, new in aliases.items():
        if new not in s and s.get(old) not in (None, ""):
            s[new] = s[old]

    cif_declared = bool(s.get("cif_declared") or s.get("_cif_declared")) or any(
        k in s for k in (
            "cif", "cif_path", "cif_file", "source_cif",
            "cif_sha256", "source_cif_sha256"))
    s.pop("_cif_declared", None)
    s["cif_declared"] = cif_declared

    cif = s.get("cif")
    if cif not in (None, "") and str(cif).strip().lower() not in (
            "none", "null", "n/a", "na", "no"):
        p = Path(str(cif))
        if not p.is_absolute():
            q = base / p
            if q.is_file():
                p = q
        s["cif"] = str(p)
        if not s.get("cif_sha256") and p.is_file():
            s["cif_sha256"] = _sha256_file(p)
    else:
        s.pop("cif", None)

    s["feff_inputs"] = _as_list(
        s.get("feff_inputs") or s.get("feff_input"))
    if s.get("feff_inputs_sha256") and not s.get("feff_input_sha256"):
        s["feff_input_sha256"] = s["feff_inputs_sha256"]
    feff_sha = _as_list(s.get("feff_input_sha256"))
    resolved: list[str] = []
    for i, item in enumerate(s["feff_inputs"]):
        p = Path(item)
        if not p.is_absolute():
            q = base / p
            if q.is_file():
                p = q
        resolved.append(str(p))
        if i >= len(feff_sha):
            feff_sha.append(_sha256_file(p) if p.is_file() else "")
        elif not feff_sha[i] and p.is_file():
            feff_sha[i] = _sha256_file(p)
    s["feff_inputs"] = resolved
    s["feff_input_sha256"] = feff_sha

    if not s.get("model_level"):
        if s.get("cif"):
            s["model_level"] = "crystal_cif"
        elif str(s.get("source_type") or "").lower() in (
                "constructed", "constructed_cluster", "cluster"):
            s["model_level"] = "cluster"
    return s


def structure_provenance(
    *,
    source_type: str | None = None,
    source_ref: str | None = None,
    short_label: str | None = None,
    phase_id: str | None = None,
    model_level: str | None = None,
    absorber_site: str | None = None,
    provenance_status: str | None = None,
    description: str | None = None,
    cif: str | Path | None = None,
    cif_sha256: str | None = None,
    feff_inputs: Sequence[str | Path] = (),
    base: str | Path | None = None,
) -> dict:
    """构造结构/CIF provenance，并自动补齐 CIF/FEFF 输入的 SHA256。

    ``cif=None`` 会显式记录“未使用 CIF”，而不是让读者误以为漏写。
    所有参数都只是描述性元数据，不参与任何拟合或数值判定。
    """
    raw = {
        "source_type": source_type,
        "source_ref": source_ref,
        "short_label": short_label,
        "phase_id": phase_id,
        "model_level": model_level,
        "absorber_site": absorber_site,
        "provenance_status": provenance_status,
        "description": description,
        "cif": cif,
        "cif_sha256": cif_sha256,
        "feff_inputs": list(feff_inputs or ()),
    }
    raw = {k: v for k, v in raw.items()
           if v not in (None, "", [], ()) or k == "cif"}
    return _normalize_structure(raw, Path(base) if base else Path.cwd())


# ================================================================
# 2. 拟合结果的统一读取
# ================================================================

@dataclass
class FitData:
    """与具体拟合实现解耦的拟合摘要。"""
    paths: list[dict]
    path_files: list[str]
    s02: float
    delta_e0: float
    delta_e0_err: float | None
    r_factor: float
    kmin: float
    kmax: float
    kweight: int
    rmin: float
    rmax: float
    n_independent: float | None = None
    n_variables: int | None = None
    cn_status: str | None = None
    cn_systematics: dict | None = None
    diagnostic_only: bool = False
    flags: dict = field(default_factory=dict)
    structure: dict = field(default_factory=dict)


def _path_to_dict(p) -> dict:
    """FitResult.paths 的 PathFit -> JSON 兼容 dict。"""
    rec = {
        "label": p.label,
        "reff": p.reff,
        "N": p.degeneracy,
        "N_err": p.n_err,
        "R": p.r,
        "R_err": p.delta_r_err,
        "sigma2": p.sigma2,
        "sigma2_err": p.sigma2_err,
        "sigma2_at_bound": p.sigma2_at_bound,
    }
    if getattr(p, "c3", None) is not None:
        rec["c3"] = p.c3
    if getattr(p, "c4", None) is not None:
        rec["c4"] = p.c4
    return rec


def fit_result_to_json(
    fit: FitResult,
    path_files: Sequence[str | Path],
    *,
    sample: str = "",
    cn_status: str | None = None,
    cn_systematics: dict | None = None,
    diagnostic_only: bool | None = None,
    flags: dict | None = None,
    notes: dict[str, str] | None = None,
    structure: dict | None = None,
) -> dict:
    """把 ``FitResult`` 转成 ``run_xafs`` 风格的 fit JSON。

    这个函数不改变任何拟合数值；只是把路径文件、N/R/sigma^2 和
    S0^2/Delta E0/R-factor 写成可复现的 JSON，供本模块和下游报告读取。
    ``structure`` 用于披露本次 CN 拟合采用的结构/CIF 来源；没有就留空，
    输出会明确写 "not stated"，不会伪造 CIF。
    """
    if len(path_files) < len(fit.paths):
        raise ValueError(
            f"path_files 数量不足: {len(path_files)} < paths {len(fit.paths)}")
    notes = notes or {}
    paths = []
    for p in fit.paths:
        rec = _path_to_dict(p)
        rec["path_file"] = str(path_files[len(paths)])
        if p.label in notes:
            rec["note"] = str(notes[p.label])
        if p.n_err is None:
            rec["n_status"] = "fixed"
        paths.append(rec)
    kweights = list(getattr(fit, "kweights", None) or (2,))
    return {
        "sample": sample,
        "mode": "paper",
        "status": "diagnostic" if diagnostic_only else "fit",
        "fit": {
            "delta_e0": fit.delta_e0,
            "delta_e0_err": fit.delta_e0_err,
            "r_factor": fit.r_factor,
            "n_independent": fit.n_independent,
            "n_variables": fit.n_variables,
            "rmin": fit.rmin,
            "rmax": fit.rmax,
            "kmin": fit.kmin,
            "kmax": fit.kmax,
            "s02": fit.s02,
            "kweights": kweights,
            "paths": paths,
        },
        "path_files": [str(p) for p in path_files[:len(fit.paths)]],
        "cn_status": cn_status or ("diagnostic_only" if diagnostic_only
                                   else "not_stated"),
        "diagnostic_only": bool(diagnostic_only),
        "cn_systematics": cn_systematics or {},
        "flags": flags or {},
        "structure": (_normalize_structure(structure, Path.cwd())
                      if structure else {}),
    }


def _resolve_path_file(raw: str | Path, base: Path) -> Path:
    p = Path(raw)
    if p.is_file():
        return p
    if not p.is_absolute():
        q = (base / p)
        if q.is_file():
            return q
    q = base / p.name
    if q.is_file():
        return q
    hits = list(base.rglob(p.name))
    if hits:
        return hits[0]
    raise FileNotFoundError(
        f"找不到 FEFF 路径文件 {raw!r}（基准目录 {base}）；"
        "fit JSON 必须带 path_files，或把 path_files 改成相对路径")


def _infer_structure(payload: dict, f: dict) -> dict:
    """从常见 S7 字段推断结构来源；只做映射，不猜 CIF。"""
    out: dict = {}
    motif = payload.get("motif") or f.get("motif")
    if isinstance(motif, dict):
        out["phase_id"] = motif.get("id")
        src = motif.get("source") or {}
        if isinstance(src, dict):
            out["source_type"] = src.get("type")
            out["source_ref"] = src.get("ref")
        out["model_level"] = motif.get("model_level")
        out["absorber_site"] = motif.get("absorber_site")
    elif isinstance(motif, str):
        out["phase_id"] = motif
        out["source_type"] = "motif"
        out["source_ref"] = f"motif {motif}; provenance not stated"

    for key in ("phase_id", "phase_source", "source_type", "source_ref",
                "source_cif_sha256", "cif", "cif_sha256", "model_level",
                "absorber_site", "provenance_status"):
        if key in out:
            continue
        val = payload.get(key)
        if val is None:
            val = f.get(key)
        if val not in (None, ""):
            out[key] = val
    return out


def load_fit(
    source: FitResult | dict | str | Path | None,
    *,
    path_files: Sequence[str | Path] | None = None,
    fit_kweight: int | None = None,
    structure: dict | None = None,
) -> FitData | None:
    """读取 fit JSON / dict / FitResult，返回统一的 FitData。"""
    if source is None:
        return None

    base = Path.cwd()
    payload: dict
    if isinstance(source, FitResult):
        if not path_files:
            raise ValueError("传入 FitResult 时必须同时给 path_files（FEFF 路径文件）")
        payload = fit_result_to_json(source, path_files, structure=structure)
        payload["path_files"] = [str(p) for p in path_files]
    elif isinstance(source, (str, Path)):
        p = Path(source)
        base = p.resolve().parent
        payload = json.loads(p.read_text(encoding="utf-8"))
    elif isinstance(source, dict):
        payload = dict(source)
        if path_files:
            payload["path_files"] = [str(p) for p in path_files]
    else:
        raise TypeError(f"不支持的 fit 来源: {type(source)!r}")

    if structure:
        payload["structure"] = {
            **(payload.get("structure") or {}), **structure}

    f = payload.get("fit") or payload
    paths = f.get("paths") or payload.get("paths") or []
    if not paths:
        raise ValueError("fit JSON 里没有 fit.paths；无法生成拟合面板或 CN 表")

    struct_raw = payload.get("structure") or f.get("structure")
    if not struct_raw:
        struct_raw = _infer_structure(payload, f)
    structure_data = _normalize_structure(struct_raw, base)

    s02 = f.get("s02", payload.get("s02"))
    if s02 is None:
        raise ValueError("fit JSON 缺少 s02；本模块不会猜 S0^2")
    kws = f.get("kweights")
    if not kws:
        kws = [f.get("kweight", payload.get("kweight", 2))]
    kws = [int(k) for k in kws]
    kweight = int(fit_kweight) if fit_kweight is not None else kws[0]

    file_list = payload.get("path_files") or f.get("path_files") or []
    if not file_list:
        # 兼容“每条路径自带 path_file”的写法
        file_list = [str(r.get("path_file", "")) for r in paths]
    file_list = [str(x) for x in file_list if str(x).strip()]

    return FitData(
        paths=[dict(r) for r in paths],
        path_files=[str(_resolve_path_file(x, base)) for x in file_list],
        s02=float(s02),
        delta_e0=float(f.get("delta_e0", payload.get("delta_e0", 0.0))),
        delta_e0_err=(None if f.get("delta_e0_err") is None
                      else float(f["delta_e0_err"])),
        r_factor=float(f.get("r_factor", payload.get("r_factor", float("nan")))),
        kmin=float(f.get("kmin", payload.get("kmin", 3.0))),
        kmax=float(f.get("kmax", payload.get("kmax", 12.0))),
        kweight=kweight,
        rmin=float(f.get("rmin", payload.get("rmin", 1.0))),
        rmax=float(f.get("rmax", payload.get("rmax", 3.0))),
        n_independent=(None if f.get("n_independent") is None
                       else float(f["n_independent"])),
        n_variables=(None if f.get("n_variables") is None
                     else int(f["n_variables"])),
        cn_status=(payload.get("cn_status") or f.get("cn_status")),
        cn_systematics=(payload.get("cn_systematics")
                        or f.get("cn_systematics") or {}),
        diagnostic_only=bool(payload.get("diagnostic_only", False)),
        flags=(payload.get("flags") or f.get("flags") or {}),
        structure=structure_data,
    )


def rebuild_model_k(fit: FitData, k: np.ndarray) -> np.ndarray:
    """按 fit JSON 的 N/R/sigma^2/Delta E0 重建模型 chi(k)。

    不重新优化任何参数。路径文件必须是生成该拟合时使用的 FEFF 文件。
    """
    if len(fit.path_files) < len(fit.paths):
        raise ValueError(
            f"FEFF 路径文件不足: {len(fit.path_files)} < {len(fit.paths)}")
    k = np.asarray(k, dtype=float)
    paths = []
    for i, rec in enumerate(fit.paths):
        p = artemis.load_feff_path(fit.path_files[i])
        n = rec.get("N", rec.get("degeneracy"))
        if n is not None:
            p.degen = float(n)
        reff = rec.get("reff")
        if reff is None:
            reff = p.reff
        if rec.get("R") is not None:
            delta_r = float(rec["R"]) - float(reff)
        else:
            delta_r = float(rec.get("delta_r", 0.0))
        if rec.get("sigma2") is None:
            raise ValueError(f"路径 {rec.get('label', i)} 缺少 sigma2")
        artemis.set_path_parameters(
            p,
            s02=fit.s02,
            delta_r=delta_r,
            sigma2=float(rec["sigma2"]),
            delta_e0=fit.delta_e0,
            third=rec.get("c3"),
            fourth=rec.get("c4"),
        )
        paths.append(p)
    return artemis.exafs_model_k(paths, k)


# ================================================================
# 3. CN 表
# ================================================================

_STATUS_LABELS = {
    "ok": "formal",
    "formal": "formal",
    "calibrated": "formal",
    "pending_consistency": "pending",
    "diagnostic_only": "diagnostic",
    "l1_diagnostic_only": "diagnostic",
    "hypothesis_model_only": "hypothesis",
    "no_cif": "no CIF",
    "consistency_failed": "rejected",
    "not_stated": "not stated",
}


def _status_label(fit: FitData) -> str:
    raw = str(fit.cn_status or "").strip().lower()
    if raw in _STATUS_LABELS:
        return _STATUS_LABELS[raw]
    if raw:
        return raw
    return "diagnostic" if fit.diagnostic_only else "not stated"


def _path_flags(fit: FitData, rec: dict) -> list[str]:
    flags: list[str] = []
    label = str(rec.get("label", ""))
    bound = set(fit.flags.get("sigma2_bound_shells") or [])
    if rec.get("sigma2_at_bound") or label in bound:
        flags.append("σ² at bound")
    n_err = rec.get("N_err")
    n_val = rec.get("N", rec.get("degeneracy"))
    if n_err is None:
        flags.append("N fixed")
    else:
        try:
            if n_val is not None and abs(float(n_err)) >= abs(float(n_val)):
                flags.append("N unconstrained")
        except (TypeError, ValueError):
            pass
    n_status = str(rec.get("n_status", "")).lower()
    if rec.get("upper_limit") or "upper" in n_status:
        flags.append("upper limit")
    # 自由文本 note 只进 CSV/Markdown，不塞进紧凑表的 Status 列。
    seen = set()
    return [x for x in flags if not (x in seen or seen.add(x))]


def cn_rows(fit: FitData) -> list[dict]:
    """生成 CN 表的行（数值、误差、状态、备注）。"""
    global_status = _status_label(fit)
    rows = []
    for rec in fit.paths:
        flags = _path_flags(fit, rec)
        if "upper limit" in flags:
            # “上限”本身就是一个独立结论，不再冒充 formal。
            status = "upper limit"
            rest = [f for f in flags if f != "upper limit"]
            if rest:
                status += " · " + " · ".join(rest)
        elif global_status == "formal" and not flags:
            status = "formal"
        elif global_status == "formal" and any(
                f in ("σ² at bound", "N unconstrained") for f in flags):
            # 形式上是 formal，但该壳层参数不可用于定量 -> 降级。
            status = " · ".join(["diagnostic"] + flags)
        elif flags:
            status = " · ".join([global_status] + flags)
        else:
            status = global_status
        rows.append({
            "shell": str(rec.get("label", "")),
            "n": rec.get("N", rec.get("degeneracy")),
            "n_err": rec.get("N_err"),
            "r": rec.get("R"),
            "r_err": rec.get("R_err"),
            "sigma2": rec.get("sigma2"),
            "sigma2_err": rec.get("sigma2_err"),
            "n_text": _fmt_pm(rec.get("N", rec.get("degeneracy")),
                              rec.get("N_err"), 2, fixed_text="fixed"),
            "r_text": _fmt_pm(rec.get("R"), rec.get("R_err"), 3,
                              fixed_text=None),
            "sigma2_text": _fmt_pm(rec.get("sigma2"),
                                   rec.get("sigma2_err"), 4,
                                   fixed_text=None),
            "cn_status": global_status,
            "flags": flags,
            "status": status,
            "note": str(rec.get("note", "") or ""),
        })
    return rows


def validate_cn_disclosure(fit: FitData) -> list[str]:
    """检查 CN 输出是否达到项目披露标准；只报告缺口，不改数值。

    规则见 ``docs/CN_REPORTING_STANDARD.md``。空列表表示最低披露要求
    已满足；任何问题都只影响报告文字/门禁，不影响拟合结果。
    """
    payload = {
        "structure": fit.structure,
        "path_files": list(fit.path_files or []),
        "cn_status": fit.cn_status,
        "cn_systematics": fit.cn_systematics or {},
        "n_independent": fit.n_independent,
        "n_variables": fit.n_variables,
        "diagnostic_only": fit.diagnostic_only,
        "fit": {
            "paths": fit.paths,
            "s02": fit.s02,
            "delta_e0": fit.delta_e0,
            "kmin": fit.kmin,
            "kmax": fit.kmax,
            "kweight": fit.kweight,
            "rmin": fit.rmin,
            "rmax": fit.rmax,
            "r_factor": fit.r_factor,
        },
    }
    return check_fit_payload(payload)


def _structure_fields(fit: FitData) -> list[tuple[str, str]]:
    """结构/CIF 来源的完整字段；缺失时明确写 not stated。"""
    s = fit.structure or {}
    if not s:
        return [("Structure / CIF", "not stated in fit JSON")]
    out: list[tuple[str, str]] = []
    for key, label in (
            ("source_type", "Source type"),
            ("source_ref", "Structure"),
            ("description", "Description"),
            ("phase_id", "Phase ID"),
            ("model_level", "Model level"),
            ("absorber_site", "Absorber site"),
            ("provenance_status", "Provenance")):
        val = s.get(key)
        if val not in (None, ""):
            out.append((label, str(val)))
    if s.get("cif"):
        out.append(("CIF", str(s["cif"])))
    elif s.get("cif_sha256"):
        out.append(("CIF", "hash only (path not stated)"))
    elif s.get("cif_declared") or s.get("_cif_declared"):
        out.append(("CIF", "none (no CIF used)"))
    else:
        out.append(("CIF", "not stated in fit JSON"))
    if s.get("cif_sha256"):
        out.append(("CIF SHA256", str(s["cif_sha256"])))
    feff = s.get("feff_inputs") or []
    if feff:
        out.append(("FEFF inputs", "; ".join(str(x) for x in feff)))
    feff_sha = [x for x in (s.get("feff_input_sha256") or []) if x]
    if feff_sha:
        out.append(("FEFF SHA256", "; ".join(str(x) for x in feff_sha)))
    return out


def _structure_lines(fit: FitData, *, compact: bool = False) -> list[str]:
    s = fit.structure or {}
    if compact:
        if not s:
            return ["Structure/CIF: not stated in fit JSON"]
        short = (s.get("short_label") or s.get("phase_id")
                 or s.get("source_type") or "not stated")
        line = f"Structure: {short}"
        if s.get("cif"):
            line += f"; CIF: {Path(str(s['cif'])).name}"
            if s.get("cif_sha256"):
                line += f" (sha256 {str(s['cif_sha256'])[:12]}…)"
        elif s.get("cif_sha256"):
            line += f"; CIF sha256: {str(s['cif_sha256'])[:12]}…"
        elif s.get("cif_declared") or s.get("_cif_declared"):
            line += "; CIF: none"
        else:
            line += "; CIF: not stated"
        lines = [line]
        feff = s.get("feff_inputs") or []
        if feff:
            lines.append("FEFF: " + ", ".join(_short_path(x) for x in feff))
        return lines
    return [f"{label}: {value}" for label, value in _structure_fields(fit)]


def _cn_footnote(fit: FitData, *, compact: bool = False) -> str:
    s02 = f"$S_0^2$ = {fit.s02:g} (fixed)"
    de0 = f"$\\Delta E_0$ = {fit.delta_e0:+.2f} eV"
    rf = f"R-factor = {fit.r_factor:.4f}"
    idp = ""
    if fit.n_independent is not None:
        idp = f"$N_{{idp}}$ = {fit.n_independent:.1f}"
        if fit.n_variables is not None:
            idp += f", variables = {fit.n_variables}"
    line1 = f"{s02}; {de0}; {rf}" + (f"; {idp}" if idp else "")
    line2 = (f"Fit: $k$ = {fit.kmin:g}–{fit.kmax:g} Å$^{{-1}}$, "
             f"$k_{{\\rm weight}}$ = {fit.kweight}; "
             f"$R$ = {fit.rmin:g}–{fit.rmax:g} Å")
    status = f"CN status: {fit.cn_status or 'not stated'}"
    budget = (fit.cn_systematics or {}).get("relative_1sigma")
    if budget is not None:
        status += f"; systematic budget ≈ ±{100 * float(budget):.0f}% (1σ)"
    if fit.diagnostic_only and "diagnostic" not in status.lower():
        status += "; diagnostic only, not for formal CN"
    lines = [line1, line2, status]
    if not compact:
        lines.append("Errors are formal only. See the fit JSON for the "
                     "full parameter set and correlation matrix.")
        lines.append("Structure / CIF provenance:")
    lines.extend(_structure_lines(fit, compact=compact))
    return "\n".join(lines)


def _table_cell_style(tbl, rows: list[dict]) -> None:
    for (r, _c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.6)
        cell.set_edgecolor("#666666")
        if r == 0:
            cell.set_facecolor("#e8e8e8")
            cell.set_text_props(fontweight="bold")
        else:
            flags = " ".join(rows[r - 1]["flags"]).lower()
            if "bound" in flags:
                cell.set_facecolor("#fff3e0")
            elif "diagnostic" in rows[r - 1]["status"].lower():
                cell.set_facecolor("#f5f8fc")


def plot_cn_table(
    fit: FitData,
    out: str | Path,
    stem: str = "cn_table",
    *,
    title: str | None = None,
    dpi: int = 300,
) -> Path:
    """单独输出一张 CN 表 PNG（论文 Supporting Information 用）。"""
    _style()
    rows = cn_rows(fit)
    headers = ["Shell", "N", "R (Å)", r"$\sigma^2$ (Å$^2$)",
               "CN status"]
    cell_text = [[r["shell"], r["n_text"], r["r_text"],
                  r["sigma2_text"], r["status"]] for r in rows]
    notes = [f"{r['shell']}: {r['note']}" for r in rows if r["note"]]
    foot = _cn_footnote(fit, compact=False)
    if notes:
        foot += "\nNotes: " + "\n       ".join(notes)
    foot_lines = foot.count("\n") + 1
    height = (1.25 + 0.42 * (len(rows) + 2)
              + 0.15 * foot_lines + (0.25 if notes else 0.0))
    fig, ax = plt.subplots(figsize=(11.0, height))
    ax.axis("off")
    ax.set_title(title or "Coordination-number fit table",
                 loc="left", fontsize=13, fontweight="bold", pad=12)
    n_rows = len(rows) + 1
    table_h = min(0.62, 0.12 * n_rows + 0.10)
    table_y = 0.90 - table_h
    tbl = ax.table(cellText=cell_text, colLabels=headers,
                   cellLoc="center", bbox=[0.0, table_y, 1.0, table_h],
                   colWidths=[0.22, 0.18, 0.20, 0.22, 0.18])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1, 1.75)
    _table_cell_style(tbl, rows)
    ax.text(0.0, max(table_y - 0.04, 0.01), foot,
            transform=ax.transAxes, va="top", fontsize=8.5, color="0.25")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{stem}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def write_cn_csv(fit: FitData, out: str | Path, stem: str = "cn_table") -> Path:
    """CN 表 CSV；每行重复全局拟合信息，便于直接进 Excel/Origin。"""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{stem}.csv"
    rows = cn_rows(fit)
    header = [
        "shell", "N", "N_err", "R_angstrom", "R_err",
        "sigma2_angstrom2", "sigma2_err", "cn_status", "flags", "note",
        "S02", "delta_e0_eV", "delta_e0_err_eV", "R_factor",
        "kmin", "kmax", "kweight", "rmin", "rmax",
        "n_independent", "n_variables", "cn_systematic_1sigma",
        "structure_source_type", "structure_source_ref", "phase_id",
        "model_level", "absorber_site", "provenance_status",
        "cif_status", "cif_file", "cif_sha256",
        "feff_inputs", "feff_input_sha256",
        "disclosure_status", "disclosure_issues",
    ]
    budget = (fit.cn_systematics or {}).get("relative_1sigma")
    s = fit.structure or {}
    disclosure_issues = validate_cn_disclosure(fit)
    struct_cols = [
        s.get("source_type") or "", s.get("source_ref") or "",
        s.get("phase_id") or "", s.get("model_level") or "",
        s.get("absorber_site") or "", s.get("provenance_status") or "",
        cif_status(s), s.get("cif") or "", s.get("cif_sha256") or "",
        "; ".join(str(x) for x in (s.get("feff_inputs") or [])),
        "; ".join(str(x) for x in (s.get("feff_input_sha256") or []) if x),
        ("warn" if disclosure_issues else "ok"),
        "; ".join(disclosure_issues),
    ]
    with open(path, "w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for r in rows:
            writer.writerow([
                r["shell"], r["n"], r["n_err"], r["r"], r["r_err"],
                r["sigma2"], r["sigma2_err"], r["cn_status"],
                "; ".join(r["flags"]), r["note"],
                fit.s02, fit.delta_e0, fit.delta_e0_err, fit.r_factor,
                fit.kmin, fit.kmax, fit.kweight, fit.rmin, fit.rmax,
                fit.n_independent, fit.n_variables, budget,
            ] + struct_cols)
    return path


def write_cn_markdown(fit: FitData, out: str | Path,
                      stem: str = "cn_table") -> Path:
    """CN 表 Markdown；便于直接贴进论文草稿或 Supplementary。"""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{stem}.md"
    rows = cn_rows(fit)
    lines = [
        "# Coordination-number fit table",
        "",
        "| Shell | N | R (Å) | σ² (Å²) | Status |",
        "|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(f"| {r['shell']} | {r['n_text']} | {r['r_text']} | "
                     f"{r['sigma2_text']} | {r['status']} |")
    lines += ["", "```", _cn_footnote(fit, compact=False), "```"]
    notes = [f"- **{r['shell']}**: {r['note']}" for r in rows if r["note"]]
    if notes:
        lines += ["", "**Notes**", ""] + notes
    issues = validate_cn_disclosure(fit)
    if issues:
        lines += ["", "**Disclosure warnings**", ""]
        lines += [f"- {issue}" for issue in issues]
    lines += ["", "**Structure / CIF provenance**", "",
              "| Field | Value |", "|---|---|"]
    for label, value in _structure_fields(fit):
        lines.append(f"| {label} | {value} |")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


# ================================================================
# 4. 数据处理
# ================================================================

def _as_result(
    source: PipelineResult | str | Path,
    *,
    absorber: str,
    kmin: float,
    kmax: float,
    kweight: int,
    rbkg: float | None,
) -> PipelineResult:
    if isinstance(source, PipelineResult):
        return source
    return xafs_pipeline.process_sample(
        [str(source)], absorber=absorber, strict=False, verbose=False,
        kmin=kmin, kmax=kmax, kweight=kweight, rbkg=rbkg)


def _spectrum(res: PipelineResult, kmin: float, kmax: float,
              kweight: int, dk: float):
    k = np.asarray(res.chi_k.k, dtype=float)
    chi = np.asarray(res.chi_k.chi, dtype=float)
    ft = athena.ft_to_r_space(k, chi, None, kmin, kmax, kweight, dk)
    return k, chi, np.asarray(ft.r, dtype=float), np.abs(ft.chi_r)


def _xanes_arrays(res: PipelineResult):
    norm = res.norm
    energy = np.asarray(norm.energy, dtype=float)
    # NormalizedScan.mu 是原始 μ；mu_norm 才是归一化 μ。
    mu = np.asarray(norm.mu_norm, dtype=float)
    return energy, mu


def _sample_e0(res: PipelineResult) -> float:
    e0 = res.e0
    for key in ("e0_sample", "e0", "E0"):
        if key in e0 and e0[key] is not None:
            return float(e0[key])
    return float(res.norm.e0)


# ================================================================
# 5. 绘图
# ================================================================

def _plot_xanes(ax, series, e0: float, window: tuple[float, float]) -> None:
    for s in series:
        ax.plot(s["energy"], s["mu"], color=s["color"], lw=1.7,
                label=s["label"], solid_capstyle="round")
    lo, hi = float(window[0]), float(window[1])
    if lo > hi:
        lo, hi = hi, lo
    ax.axvline(e0, color="0.65", lw=0.9, ls=(0, (2, 3)), zorder=0)
    ax.text(e0, 0.02, " $E_0$", transform=ax.get_xaxis_transform(),
            ha="left", va="bottom", fontsize=8.5, color="0.45")
    ax.set_xlim(e0 + lo, e0 + hi)
    vals = []
    for s in series:
        m = (s["energy"] >= e0 + lo) & (s["energy"] <= e0 + hi)
        if np.any(m):
            vals.append(s["mu"][m])
    if vals:
        allv = np.concatenate(vals)
        ymin = min(-0.05, float(np.nanmin(allv)) - 0.03)
        ymax = float(np.nanmax(allv)) + 0.08
        ax.set_ylim(ymin, ymax)
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel(r"Normalized $\mu(E)$")
    ax.legend(loc="best", frameon=False, fontsize=9)
    _tidy_axis(ax)


def _plot_k(ax, series, kmin: float, kmax: float,
            kweight: int) -> None:
    n = int(kweight)
    ymax = 0.0
    xlo = max(0.0, kmin - 1.0)
    xhi = kmax + 1.0
    for s in series:
        m = (s["k"] >= xlo) & (s["k"] <= xhi)
        y = s["k"][m] ** n * s["chi"][m]
        ax.plot(s["k"][m], y, color=s["color"], lw=1.25,
                label=s["label"], solid_capstyle="round")
        if y.size:
            ymax = max(ymax, float(np.nanmax(np.abs(y))))
    ax.axhline(0.0, color="0.65", lw=0.8, zorder=0)
    ax.set_xlim(xlo, xhi)
    if ymax > 0:
        ax.set_ylim(-1.12 * ymax, 1.12 * ymax)
    ax.set_xlabel(r"$k$ (Å$^{-1}$)")
    ax.set_ylabel(rf"$k^{{{n}}}\chi(k)$ (Å$^{{-{n}}}$)")
    _tidy_axis(ax)


def _draw_ft_guides(ax, guides: list[FTGuide], series,
                    ymax: float) -> None:
    """画 |FT| 标线；标签自动避让，允许调用者给 text_y。"""
    used: list[tuple[float, float]] = []
    for g in guides:
        if g.x_end is None:
            ax.axvline(g.x, color="0.45", lw=1.0, ls=(0, (5, 4)), zorder=0)
            vals = []
            for s in series:
                r, mag = s["r"], s["mag"]
                if r.size and r[0] <= g.x <= r[-1]:
                    vals.append(float(np.interp(g.x, r, mag)))
            y = max(vals) if vals else 0.0
            ty = g.text_y if g.text_y is not None else y + 0.07 * ymax
            ty = min(ty, 0.96 * ymax)
            for ux, uy in used:
                while abs(g.x - ux) < 0.38 and abs(ty - uy) < 0.085 * ymax:
                    ty -= 0.085 * ymax
            ax.text(g.x, ty, g.label, ha="center", va="bottom",
                    fontsize=11, fontweight="bold", color="0.15")
            used.append((g.x, ty))
        else:
            lo, hi = sorted((g.x, g.x_end))
            y = g.text_y if g.text_y is not None else 0.10 * ymax
            ax.annotate("", xy=(hi, y), xytext=(lo, y), zorder=6,
                        arrowprops=dict(arrowstyle="<|-|>", color="0.35",
                                        lw=1.2, shrinkA=0, shrinkB=0))
            ax.text((lo + hi) / 2, y + 0.022 * ymax, g.label,
                    ha="center", va="bottom", fontsize=11,
                    fontweight="bold", color="0.15", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white",
                              ec="none", alpha=0.8))


def _plot_ft(ax, series, guides: list[FTGuide], rmax_plot: float,
             rbkg: float | None, shade_rbkg: bool) -> None:
    ymax = 0.0
    for s in series:
        m = (s["r"] >= 0.0) & (s["r"] <= rmax_plot)
        ax.plot(s["r"][m], s["mag"][m], color=s["color"], lw=1.7,
                label=s["label"], solid_capstyle="round")
        if np.any(m):
            ymax = max(ymax, float(np.nanmax(s["mag"][m])))
    ymax = ymax if ymax > 0 else 1.0
    if shade_rbkg and rbkg is not None and 0 < rbkg < rmax_plot:
        ax.axvspan(0.0, rbkg, color="0.93", zorder=0)
        ax.text(rbkg * 0.5, 0.018 * ymax, r"$R < R_{\rm bkg}$",
                ha="center", va="bottom", fontsize=8, color="0.45",
                zorder=4,
                bbox=dict(boxstyle="round,pad=0.1", fc="white",
                          ec="none", alpha=0.7))
    _draw_ft_guides(ax, guides, series, ymax)
    ax.set_xlim(0.0, rmax_plot)
    ax.set_ylim(0.0, 1.28 * ymax)
    ax.set_xlabel("Radial distance (Å)")
    ax.set_ylabel(rf"$|\mathrm{{FT}}\,(k^{{{int(series[0]['kweight'])}}}"
                  rf"\chi(k))|$ (Å$^{{-{int(series[0]['kweight']) + 1}}}$)")
    _tidy_axis(ax)


def _fit_guide_positions(fit: FitData, r: np.ndarray,
                         mag: np.ndarray, rmax_plot: float) -> list[tuple[str, float]]:
    """为每条拟合路径找模型 |FT| 的局部峰，作为 R 空间标线。

    未相位修正的 FT 峰位比真实键长短 ~0.2–0.5 Å，因此搜索窗口取
    ``[R-0.75, R-0.15]``。找不到峰时退回 ``R-0.35``。
    """
    out = []
    for rec in fit.paths:
        try:
            center = float(rec.get("R", rec.get("reff", 0.0)))
        except (TypeError, ValueError):
            continue
        if center <= 0 or center > rmax_plot + 0.5:
            continue
        lo = max(0.0, center - 0.75)
        hi = max(lo + 0.05, center - 0.15)
        m = (r >= lo) & (r <= hi)
        if np.any(m):
            local = np.where(m)[0]
            idx = local[int(np.argmax(mag[local]))]
            pos = float(r[idx])
        else:
            pos = max(0.0, center - 0.35)
        out.append((str(rec.get("label", "")), pos))
    return out


def _plot_fit_r(ax, r: np.ndarray, data_mag: np.ndarray,
                model_mag: np.ndarray, fit: FitData,
                rmax_plot: float) -> None:
    ax.axvspan(fit.rmin, fit.rmax, color="#e8f0fe", zorder=0)
    ax.plot(r, data_mag, color=DATA_COLOR, lw=1.7, label="data",
            solid_capstyle="round")
    ax.plot(r, model_mag, color=FIT_COLOR, lw=1.5, ls="--", label="fit")
    ymax = 1.0
    m = (r >= 0) & (r <= rmax_plot)
    if np.any(m):
        ymax = max(float(np.nanmax(data_mag[m])),
                   float(np.nanmax(model_mag[m])), 1e-12)
    guides = _fit_guide_positions(fit, r, model_mag, rmax_plot)
    for i, (label, pos) in enumerate(guides):
        ax.axvline(pos, color="0.55", lw=0.8, ls=(0, (2, 3)),
                   zorder=0)
        ax.text(pos, ymax * (0.93 - 0.10 * i), label, ha="center",
                va="bottom", fontsize=9.5, fontweight="bold",
                color="0.2",
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec="none", alpha=0.75))
    ax.set_xlim(0.0, rmax_plot)
    ax.set_ylim(0.0, 1.28 * ymax)
    ax.set_xlabel("Radial distance (Å)")
    ax.set_ylabel(rf"$|\mathrm{{FT}}|$ (Å$^{{-{fit.kweight + 1}}}$)")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    # 放在左下空白区，避开右侧高壳层和 “not phase-corrected” 标注。
    ax.text(0.03, 0.62,
            f"R-factor = {fit.r_factor:.4f}\n"
            f"$\\Delta E_0$ = {fit.delta_e0:+.1f} eV",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.5,
            color="0.25", zorder=5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white",
                      ec="0.8", alpha=0.85))
    ax.text(0.98, 0.02, "not phase-corrected", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7.5, color="0.45",
            bbox=dict(boxstyle="round,pad=0.15", fc="white",
                      ec="none", alpha=0.75))
    _tidy_axis(ax)


def _plot_fit_k(ax, k: np.ndarray, data_chi: np.ndarray,
                model_chi: np.ndarray, fit: FitData) -> None:
    n = int(fit.kweight)
    yd = k ** n * data_chi
    ym = k ** n * model_chi
    ymax = max(float(np.nanmax(np.abs(yd))),
               float(np.nanmax(np.abs(ym))), 1e-12)
    ax.plot(k, yd, color=DATA_COLOR, lw=1.5, label="data")
    ax.plot(k, ym, color=FIT_COLOR, lw=1.3, ls="--", label="fit")
    ax.axhline(0.0, color="0.65", lw=0.8, zorder=0)
    ax.set_xlim(float(k.min()), float(k.max()))
    ax.set_ylim(-1.15 * ymax, 1.15 * ymax)
    ax.set_xlabel(r"$k$ (Å$^{-1}$)")
    ax.set_ylabel(rf"$k^{{{n}}}\chi(k)$ (Å$^{{-{n}}}$)")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    _tidy_axis(ax)


def _two_line_pm(text: str) -> str:
    """紧凑表里把 ``a ± b`` 折成两行，避免 R / sigma^2 列互相挤压。"""
    if " ± " in text:
        value, err = text.split(" ± ", 1)
        return f"{value}\n± {err}"
    return text


def _draw_cn_panel(ax, fit: FitData) -> None:
    ax.axis("off")
    rows = cn_rows(fit)
    headers = ["Shell", "N", "R (Å)", r"$\sigma^2$ (Å$^2$)",
               "CN status"]
    cell_text = [[r["shell"], _two_line_pm(r["n_text"]),
                  _two_line_pm(r["r_text"]),
                  _two_line_pm(r["sigma2_text"]), r["status"]]
                 for r in rows]
    n_rows = len(rows) + 1
    table_h = min(0.70, 0.17 * n_rows + 0.12)
    table_y = 0.94 - table_h
    tbl = ax.table(cellText=cell_text, colLabels=headers,
                   cellLoc="center", bbox=[0.0, table_y, 1.0, table_h],
                   colWidths=[0.22, 0.18, 0.20, 0.22, 0.18])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.0)
    tbl.scale(1, 1.35)
    _table_cell_style(tbl, rows)
    ax.text(0.0, max(table_y - 0.035, 0.01),
            _cn_footnote(fit, compact=True),
            transform=ax.transAxes, va="top", fontsize=7.2,
            color="0.25")


def _series_payload(label: str, res: PipelineResult, color: str,
                    kmin: float, kmax: float, kweight: int,
                    dk: float) -> dict:
    energy, mu = _xanes_arrays(res)
    k, chi, r, mag = _spectrum(res, kmin, kmax, kweight, dk)
    return {
        "label": label,
        "color": color,
        "energy": energy,
        "mu": mu,
        "k": k,
        "chi": chi,
        "r": r,
        "mag": mag,
        "kweight": int(kweight),
    }


def _write_ft_csv(series: list[dict], out: Path, stem: str) -> Path:
    """把 |FT| 面板实际绘制的数据存成 CSV（便于复现/重绘）。"""
    r_ref = series[0]["r"]
    cols = {"R": r_ref}
    for s in series:
        cols[_clean_stem(s["label"])] = np.interp(r_ref, s["r"], s["mag"])
    path = out / f"{stem}_ft.csv"
    keys = list(cols)
    with open(path, "w", encoding="utf-8-sig", newline="") as fh:
        fh.write(",".join(keys) + "\n")
        for i in range(len(r_ref)):
            fh.write(",".join(f"{cols[k][i]:.6g}" for k in keys) + "\n")
    return path


# ================================================================
# 6. 主入口
# ================================================================

def make_paper_figure(
    sample: PipelineResult | str | Path,
    *,
    absorber: str,
    refs: Sequence[tuple[str, PipelineResult | str | Path]] = (),
    sample_label: str | None = None,
    edge: str | None = None,
    kmin: float = 3.0,
    kmax: float = 12.0,
    kweight: int = 3,
    dk: float = 1.0,
    rbkg: float | None = None,
    xanes_window: tuple[float, float] = (-30.0, 60.0),
    ft_guides: Iterable[Any] = (),
    fit: FitResult | dict | str | Path | None = None,
    fit_path_files: Sequence[str | Path] | None = None,
    fit_kweight: int | None = None,
    structure: dict | None = None,
    strict_disclosure: bool = False,
    fit_r_plot: tuple[float, float] | None = None,
    rmax_plot: float = 6.0,
    shade_rbkg: bool = True,
    title: str | None = None,
    out: str | Path = ".",
    stem: str = "paper",
    dpi: int = 300,
    verbose: bool = True,
) -> dict[str, Path]:
    """生成论文主图 + CN 表（CSV/Markdown/PNG）。

    Parameters
    ----------
    sample :
        ``PipelineResult``，或原始谱文件路径（走 ``process_sample``）。
    absorber :
        吸收元素符号；只用于 S2 定边、S0^2 查表和默认边标签。
    refs :
        ``[(label, PipelineResult | path), ...]``；参考谱与样品使用同一
        kmin/kmax/kweight 和 Rbkg（未显式给 Rbkg 时继承样品值）。
    ft_guides :
        |FT| 面板的标线，接受 ``FTGuide`` / dict / ``"LABEL=X"`` /
        ``"LABEL=X..Y"``。
    fit :
        ``FitResult`` / fit JSON 路径 / dict。None 时不画 (d)(e) 和表。
    fit_path_files :
        传 ``FitResult`` 时必填；传 JSON 时忽略（JSON 自带 path_files）。
    fit_kweight :
        覆盖 fit JSON 里的 k 权重（只影响 (d)(e) 的显示和 FT，不改拟合）。
    structure :
        结构/CIF provenance dict；与 fit JSON 自带的 ``structure`` 合并
        （本参数优先）。可用 ``structure_provenance(...)`` 构造。它只影响
        图注和 CN 表的披露文字，不参与任何拟合或数值判定。
    strict_disclosure :
        True 时若 CN 披露不完整（见 ``docs/CN_REPORTING_STANDARD.md``）
        直接抛 ``ValueError``；默认只打印警告。
    fit_r_plot :
        R 空间拟合面板的横轴范围；默认 ``0 .. max(3.5, rmax+0.8)``。
    """
    _style()
    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = get_absorber(absorber)
    edge_label = (edge or cfg.edge or "").strip()
    sample_stem = _clean_stem(Path(sample).stem if isinstance(sample, (str, Path))
                              else (getattr(sample, "stem", "") or "sample"))
    sample_label = sample_label or sample_stem
    stem = stem or f"{sample_stem}_paper"
    guides = _coerce_guides(ft_guides)

    # 样品：显式 Rbkg 优先；否则让管线自动选。
    res_sample = _as_result(sample, absorber=absorber, kmin=kmin, kmax=kmax,
                            kweight=kweight, rbkg=rbkg)
    actual_rbkg = (float(rbkg) if rbkg is not None
                   else float(res_sample.rbkg["rbkg"]))
    series = [_series_payload(sample_label, res_sample, SERIES_COLORS[0],
                              kmin, kmax, kweight, dk)]
    for i, (label, source) in enumerate(refs, start=1):
        res_ref = _as_result(source, absorber=absorber, kmin=kmin, kmax=kmax,
                             kweight=kweight, rbkg=actual_rbkg)
        series.append(_series_payload(
            label, res_ref, SERIES_COLORS[i % len(SERIES_COLORS)],
            kmin, kmax, kweight, dk))

    fit_data = load_fit(fit, path_files=fit_path_files,
                        fit_kweight=fit_kweight, structure=structure)

    has_fit = fit_data is not None
    if has_fit and len(fit_data.path_files) < len(fit_data.paths):
        raise ValueError(
            "fit JSON 缺少足够的 path_files；无法重建 (d)(e) 面板。"
            "请确认 fit JSON 含 path_files，或传 fit_path_files。")

    disclosure_issues = validate_cn_disclosure(fit_data) if has_fit else []
    if disclosure_issues:
        if strict_disclosure:
            raise ValueError(
                "CN disclosure standard not met:\n  - "
                + "\n  - ".join(disclosure_issues))
        if verbose:
            print("[paper] CN disclosure warnings:")
            for issue in disclosure_issues:
                print(f"        - {issue}")

    # ---- 画布 ----
    if has_fit:
        fig = plt.figure(figsize=(15.5, 8.8))
        gs = GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.32)
        ax_a = fig.add_subplot(gs[0, 0])
        ax_b = fig.add_subplot(gs[0, 1])
        ax_c = fig.add_subplot(gs[0, 2])
        ax_d = fig.add_subplot(gs[1, 0])
        ax_e = fig.add_subplot(gs[1, 1])
        ax_f = fig.add_subplot(gs[1, 2])
    else:
        fig = plt.figure(figsize=(15.0, 4.8))
        gs = GridSpec(1, 3, figure=fig, wspace=0.32)
        ax_a = fig.add_subplot(gs[0, 0])
        ax_b = fig.add_subplot(gs[0, 1])
        ax_c = fig.add_subplot(gs[0, 2])
        ax_d = ax_e = ax_f = None

    e0 = _sample_e0(res_sample)
    _plot_xanes(ax_a, series, e0, xanes_window)
    _panel_title(ax_a, "a", "XANES")
    if edge_label:
        ax_a.text(0.02, 0.95, f"{cfg.symbol} {edge_label}-edge",
                  transform=ax_a.transAxes, ha="left", va="top",
                  fontsize=9, color="0.35")

    _plot_k(ax_b, series, kmin, kmax, kweight)
    _panel_title(ax_b, "b", rf"$k^{{{kweight}}}\chi(k)$")

    _plot_ft(ax_c, series, guides, rmax_plot, actual_rbkg, shade_rbkg)
    _panel_title(ax_c, "c", rf"$|\mathrm{{FT}}\,(k^{{{kweight}}}\chi(k))|$")

    if has_fit:
        k_fit = np.asarray(res_sample.chi_k.k, dtype=float)
        chi_fit = np.asarray(res_sample.chi_k.chi, dtype=float)
        sel = (k_fit >= fit_data.kmin) & (k_fit <= fit_data.kmax)
        k_fit = k_fit[sel]
        chi_fit = chi_fit[sel]
        model_k = rebuild_model_k(fit_data, k_fit)
        ft_data = athena.ft_to_r_space(
            k_fit, chi_fit, None, fit_data.kmin, fit_data.kmax,
            fit_data.kweight, dk)
        ft_model = athena.ft_to_r_space(
            k_fit, model_k, None, fit_data.kmin, fit_data.kmax,
            fit_data.kweight, dk)
        r_fit = np.asarray(ft_data.r, dtype=float)
        data_mag = np.abs(ft_data.chi_r)
        model_mag = np.abs(ft_model.chi_r)

        r_lo, r_hi = (fit_r_plot if fit_r_plot is not None
                      else (0.0, max(3.5, fit_data.rmax + 0.8)))
        if r_hi <= r_lo:
            r_hi = r_lo + 4.0
        _plot_fit_r(ax_d, r_fit, data_mag, model_mag, fit_data, r_hi)
        _panel_title(ax_d, "d", "R-space fit")

        _plot_fit_k(ax_e, k_fit, chi_fit, model_k, fit_data)
        _panel_title(ax_e, "e", rf"$k^{{{fit_data.kweight}}}\chi(k)$ fit")

        _draw_cn_panel(ax_f, fit_data)
        _panel_title(ax_f, "f", "Coordination numbers")

    if title:
        fig.suptitle(title, fontsize=14, fontweight="bold", y=0.985)
        rect_top = 0.93
    else:
        rect_top = 0.97
    fig.tight_layout(rect=(0.0, 0.02, 1.0, rect_top))

    fig_path = out_dir / f"{stem}.png"
    fig.savefig(fig_path, dpi=dpi)
    plt.close(fig)

    outputs: dict[str, Path] = {
        "figure": fig_path,
        "ft_csv": _write_ft_csv(series, out_dir, stem),
    }
    if has_fit:
        outputs["cn_table_png"] = plot_cn_table(
            fit_data, out_dir, f"{stem}_cn_table", dpi=dpi)
        outputs["cn_table_csv"] = write_cn_csv(
            fit_data, out_dir, f"{stem}_cn_table")
        outputs["cn_table_md"] = write_cn_markdown(
            fit_data, out_dir, f"{stem}_cn_table")

    if verbose:
        print(f"[paper] sample = {sample_label}  E0 = {e0:.2f} eV  "
              f"Rbkg = {actual_rbkg:.2f} Å  "
              f"k = {kmin:g}–{kmax:g} Å^-1  kweight = {kweight}")
        for s in series:
            m = (s["r"] >= 0) & (s["r"] <= rmax_plot)
            peak = (s["r"][m][int(np.argmax(s["mag"][m]))]
                    if np.any(m) else float("nan"))
            print(f"        {s['label']:<24} |FT| max = {np.nanmax(s['mag'][m]):7.2f}  "
                  f"main peak = {peak:.3f} Å")
        if has_fit:
            print(f"[paper] fit = {len(fit_data.paths)} path(s), "
                  f"Rf = {fit_data.r_factor:.4f}, "
                  f"CN status = {fit_data.cn_status or 'not stated'}")
        for key, path in outputs.items():
            print(f"        wrote {path}")
    return outputs


# ================================================================
# 7. CLI
# ================================================================

def _parse_ref(text: str) -> tuple[str, str]:
    label, sep, path = text.partition("=")
    if not sep or not label.strip() or not path.strip():
        raise argparse.ArgumentTypeError(
            f"--ref 格式应为 LABEL=PATH，收到: {text!r}")
    return label.strip(), path.strip()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="通用 XAFS 论文出图：XANES / k^n chi(k) / |FT| / "
                    "R 空间拟合 / k 空间拟合 / CN 表",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--sample", required=True,
                   help="样品谱文件（XDI / prj / 两列 μ(E) 等）")
    p.add_argument("--sample-label", default=None, help="样品图例名")
    p.add_argument("--ref", action="append", default=[], type=_parse_ref,
                   metavar="LABEL=PATH",
                   help="参考谱，可重复；例如 --ref 'M foil=foil.dat'")
    p.add_argument("--absorber", required=True, help="吸收元素符号")
    p.add_argument("--edge", default=None, help="边标签（默认读 absorber 表）")
    p.add_argument("--kmin", type=float, default=3.0)
    p.add_argument("--kmax", type=float, default=12.0)
    p.add_argument("--kweight", type=int, default=3,
                   help="(a)-(c) 面板的 k 权重")
    p.add_argument("--dk", type=float, default=1.0, help="FT 窗过渡宽度")
    p.add_argument("--rbkg", type=float, default=None,
                   help="背景扣除 Rbkg；默认样品自动、参考谱继承样品值")
    p.add_argument("--xanes-window", type=float, nargs=2,
                   default=(-30.0, 60.0), metavar=("LO", "HI"),
                   help="XANES 横轴相对 E0 的范围 (eV)")
    p.add_argument("--guide", action="append", default=[],
                   metavar="LABEL=X[..Y]",
                   help="|FT| 标线；区间写作 LABEL=X..Y，可重复")
    p.add_argument("--fit-json", default=None,
                   help="S7 fit JSON；提供后增加 (d)(e) 和 CN 表")
    p.add_argument("--fit-kweight", type=int, default=None,
                   help="覆盖 fit JSON 的 k 权重（只影响显示）")
    p.add_argument("--structure-json", default=None,
                   help="结构/CIF provenance JSON；与 fit JSON 自带的 "
                        "structure 合并后写入图注、CN 表和 CSV")
    p.add_argument("--strict-disclosure", action="store_true",
                   help="CN 披露不完整时直接报错；默认只打印警告")
    p.add_argument("--fit-r-plot", type=float, nargs=2, default=None,
                   metavar=("LO", "HI"), help="R 空间拟合面板横轴范围")
    p.add_argument("--rmax-plot", type=float, default=6.0,
                   help="(c) |FT| 面板横轴上限")
    p.add_argument("--no-rbkg-shade", action="store_true",
                   help="不标注 R < Rbkg 的非解释区")
    p.add_argument("--title", default=None, help="整图标题")
    p.add_argument("--out", default=".", help="输出目录")
    p.add_argument("--stem", default=None, help="输出文件名前缀")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--quiet", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    guides = [parse_guide(g) for g in args.guide]
    structure = None
    if args.structure_json:
        sp = Path(args.structure_json)
        structure = _normalize_structure(
            json.loads(sp.read_text(encoding="utf-8")), sp.resolve().parent)
    sample_stem = _clean_stem(Path(args.sample).stem)
    stem = args.stem or f"{sample_stem}_paper"
    make_paper_figure(
        sample=args.sample,
        absorber=args.absorber,
        refs=[(label, path) for label, path in args.ref],
        sample_label=args.sample_label,
        edge=args.edge,
        kmin=args.kmin,
        kmax=args.kmax,
        kweight=args.kweight,
        dk=args.dk,
        rbkg=args.rbkg,
        xanes_window=tuple(args.xanes_window),
        ft_guides=guides,
        fit=args.fit_json,
        fit_kweight=args.fit_kweight,
        structure=structure,
        strict_disclosure=args.strict_disclosure,
        fit_r_plot=(tuple(args.fit_r_plot) if args.fit_r_plot else None),
        rmax_plot=args.rmax_plot,
        shade_rbkg=not args.no_rbkg_shade,
        title=args.title,
        out=args.out,
        stem=stem,
        dpi=args.dpi,
        verbose=not args.quiet,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
