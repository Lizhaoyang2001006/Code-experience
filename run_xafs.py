# -*- coding: utf-8 -*-
"""
run_xafs.py — 命令行入口：解析参数、调 xafs_pipeline、落盘
================================================================
本文件只做三件事：解析参数、调用 xafs_pipeline.process_sample、写文件。
**不含任何物理计算**——数学在 athena / artemis / larch_wt，
编排在 xafs_pipeline，诊断图在 xafs_diagnostics。

用法:
    python run_xafs.py <数据文件...> [--absorber Fe] [-o 输出目录] [--align]
                        [--identify] [--fit] [--fit-motif ID]
                        [--fit-rmin A] [--fit-rmax A] [--allow-ambiguous]
                        [--fit-weights]
                        [--l1-ligand O] [--l1-geometry oct]
                        [--strict/--no-strict] [--diagnostics]
    python run_xafs.py                     # 缺参数时跑内置示例（Fe 标样）

    --fit 在识别 confident 时自动做多壳层 CN；partial/ambiguous/no_match
    时只做 L1 第一壳。S₀² 未标定的元素一律拒绝 CN 定量。

处理链（xafs_pipeline.process_sample 内部）:
    S1   xafs_io.read_scan          → RawScan
    S1.5 athena.merge_scans         → 合并（可选 align_by_reference）
    S2   athena.determine_e0        → E0
    S3   athena.normalize_absorption→ NormalizedScan
    S4   athena.optimize_rbkg       → Rbkg
    S5   athena.extract_chi_k / ft_to_r_space → ChiK / FourierTransform
    S6   larch_wt.continuous_wavelet_transform → WaveletMap
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

import absorber as absorber_mod
import cn_table
import cif_gate
import phase_resolver
import xafs_config
import xafs_diagnostics
import xafs_fit
import xafs_io
import xafs_lcf
import xafs_pipeline
from xafs_core import XAFSConfigError

PROJECT_ROOT = Path(__file__).resolve().parent
T5_CALIBRATION_DIR = PROJECT_ROOT / "data" / "s02_calibration"

# Windows 控制台默认 GBK；帮助文本含 Å/σ²/S₀²/α 等字符，
# 不重配会在 --help 时抛 UnicodeEncodeError。改 UTF-8 + replace
# 只影响控制台，不改任何数据文件编码。
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")


def _fit_outputs(fit_out, outdir: Path, stem: str) -> dict:
    """把 S7 结果写成 CSV + JSON，返回可并入 meta 的摘要。"""
    mode = fit_out["mode"]
    if mode == "l1_first_shell":
        l1 = fit_out["l1"]
        rows = [[l1.get("ligand", ""), l1.get("N"), l1.get("N_err"),
                 l1.get("R"), l1.get("R_err"), l1.get("sigma2"), "",
                 l1.get("r0"), l1.get("delta_e0"), l1.get("rf"),
                 l1.get("n_independent"), l1.get("n_variables")]]
        meta = {k: v for k, v in fit_out.items() if k != "structure"}
    else:
        fit = fit_out["fit"]
        rows = []
        for p in fit.paths:
            rows.append([p.label, p.degeneracy, p.n_err, p.r, p.delta_r_err,
                         p.sigma2, p.sigma2_err, p.reff, fit.delta_e0,
                         fit.r_factor, fit.n_independent, fit.n_variables])
        meta = {k: v for k, v in fit_out.items() if k != "structure"}
        meta["fit"] = {
            "delta_e0": fit.delta_e0, "delta_e0_err": fit.delta_e0_err,
            "r_factor": fit.r_factor, "n_independent": fit.n_independent,
            "n_variables": fit.n_variables, "overfit": fit.overfit,
            "rmin": fit.rmin, "rmax": fit.rmax,
            "kmin": fit.kmin, "kmax": fit.kmax, "s02": fit.s02,
            "kweights": list(fit.kweights or ()),
            "r_factors_by_kweight": (fit.r_factors_by_kweight or {}),
            "epsilon_r": fit.epsilon_r,
            "chi_square": fit.chi_square,
            "reduced_chi_square": fit.reduced_chi_square,
            "mixture_weights": fit.mixture_weights,
            "parameter_names": fit.parameter_names,
            "correlation": ([[round(float(v), 5) for v in row]
                             for row in fit.correlation]
                            if fit.correlation is not None else None),
            "paths": [{"label": p.label, "reff": p.reff,
                       "N": p.degeneracy, "N_err": p.n_err,
                       "R": p.r, "R_err": p.delta_r_err,
                       "sigma2": p.sigma2, "sigma2_err": p.sigma2_err,
                       "sigma2_at_bound": p.sigma2_at_bound}
                      for p in fit.paths]}
    header = ["shell", "N", "N_err", "R_A", "R_err", "sigma2",
              "sigma2_err", "Reff_A", "delta_e0", "R_factor",
              "N_independent", "N_variables"]
    with open(outdir / f"{stem}_fit.csv", "w", newline="",
              encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for row in rows:
            w.writerow(["" if v is None else
                        (f"{v:.6g}" if isinstance(v, float) else v)
                        for v in row])
    xafs_io.write_json(outdir / f"{stem}_fit.json", meta)
    return meta


_ELEMENT_RE = re.compile(r"^[A-Za-z]{1,2}$")


def _normalize_element(value) -> str | None:
    """元素符号归一化并验证是否在 absorber 表内；不猜未知符号。"""
    if not value:
        return None
    s = str(value).strip()
    if not _ELEMENT_RE.match(s):
        return None
    el = s[0].upper() + s[1:].lower()
    try:
        absorber_mod.get_absorber(el)
    except Exception:  # noqa: BLE001 - 未知元素按“未检测到”处理
        return None
    return el


def _sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _project_rel(path) -> str:
    try:
        return str(Path(path).resolve().relative_to(PROJECT_ROOT))
    except ValueError:
        return str(Path(path).resolve())


def load_t5_s02(element, root=None) -> tuple[float, str]:
    """读取 T5 独立标定产物中的 S₀²，并返回可写入正式表的溯源串。

    只接受 ``data/s02_calibration/<元素>/calibration.json`` 中
    ``status="calibrated"`` 的结果。未通过元素**明确拒绝**，不静默
    回退到静态表或猜测值。该函数不修改 ``absorber.py``（#60）。
    """
    el = _normalize_element(element)
    if el is None:
        raise XAFSConfigError(
            "T5", f"未知吸收元素 {element!r}；无法读取 T5 标定产物")
    root = Path(root) if root is not None else T5_CALIBRATION_DIR
    path = root / el / "calibration.json"
    if not path.is_file():
        raise XAFSConfigError(
            "T5", f"{el} 缺少 T5 标定产物 {_project_rel(path)}；"
                   f"先运行 calibrate_t5.py，或显式给 --s02")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise XAFSConfigError(
            "T5", f"{_project_rel(path)} 读取/解析失败: {exc}") from exc
    if str(payload.get("absorber", "")).capitalize() != el:
        raise XAFSConfigError(
            "T5", f"{_project_rel(path)} 的 absorber="
                   f"{payload.get('absorber')!r} 与请求元素 {el!r} 不一致")
    status = payload.get("status")
    if status != "calibrated":
        reason = payload.get("reason") or "未通过标定门禁"
        raise XAFSConfigError(
            "T5", f"{el} 的 T5 标定状态={status!r}，拒绝自动注入 S₀²："
                   f"{reason}；请显式 --s02 或先补标定")
    if payload.get("static_absorber_py_unchanged") is not True:
        raise XAFSConfigError(
            "T5", f"{_project_rel(path)} 未声明 "
                   f"static_absorber_py_unchanged=True；拒绝作为 T5 产物使用")
    try:
        s02 = float(payload.get("s02"))
    except (TypeError, ValueError) as exc:
        raise XAFSConfigError(
            "T5", f"{_project_rel(path)} 的 s02={payload.get('s02')!r} "
                   f"不是数值") from exc
    if not (0.4 <= s02 <= 1.3):
        raise XAFSConfigError(
            "T5", f"{el} 的 T5 s02={s02} 超出 [0.4, 1.3]")
    calibration_sha = _sha256_file(path)
    source = (f"s02_calibration:{_project_rel(path)}"
              f"#sha256:{calibration_sha[:12]}")
    standard_sha = payload.get("standard_sha256")
    if standard_sha:
        source += f";standard_sha256:{str(standard_sha)[:12]}"
    return s02, source


def _resolve_s02_for_paths(args, paths, absorber_hint=None):
    """按 CLI 优先级解析 S₀²：显式 --s02 > --s02-from-t5 > 静态表。

    返回 ``(s02, source)``；未启用任何覆盖时返回 ``(None, None)``，
    由 ``xafs_fit`` 沿用 absorber.py 静态表。
    """
    if getattr(args, "s02", None) is not None:
        return float(args.s02), "--s02 runtime override"
    if not getattr(args, "s02_from_t5", False):
        return None, None
    hint = absorber_hint if absorber_hint is not None \
        else getattr(args, "absorber", None)
    if hint:
        # 显式 --absorber 优先：避免文件头与用户指定元素不一致时，
        # S₀² 与真正用于拟合的元素错配。
        element = _normalize_element(hint)
    else:
        element, _source = detect_absorber(paths[0])
    if element is None:
        raise XAFSConfigError(
            "T5", f"{Path(paths[0]).name}: 启用 --s02-from-t5 但无法判定"
                   f"吸收元素；请显式给 --absorber")
    return load_t5_s02(element)


def detect_absorber(path, fallback: str | None = None) -> tuple[str | None, str]:
    """从 XDI/纯文本头自动判定吸收元素。

    优先级：XDI ``Element.symbol`` / 纯文本 ``Absorber:`` 头 >
    ``RawScan.meta['element']`` > CLI ``--absorber`` 兜底。
    返回 ``(element, source)``；``source`` 用于报告溯源。
    """
    p = Path(path)
    try:
        head = p.read_text(encoding="utf-8", errors="replace")[:8192]
    except OSError:
        head = ""
    patterns = [
        ("xdi_element", r"(?im)^\s*#\s*Element\.symbol\s*:\s*([A-Za-z]{1,2})\b"),
        ("xdi_element", r"(?im)^\s*#\s*Element\s*:\s*([A-Za-z]{1,2})\b"),
        ("header_absorber", r"(?im)^\s*#\s*Absorber\s*[:=]\s*([A-Za-z]{1,2})\b"),
    ]
    for source, pat in patterns:
        m = re.search(pat, head)
        if m:
            el = _normalize_element(m.group(1))
            if el:
                return el, source
    try:
        scan = xafs_io.read_scan(p)
        el = _normalize_element((scan.meta or {}).get("element"))
        if el:
            return el, "scan_meta"
    except Exception:  # noqa: BLE001 - 头信息已失败，读整谱失败不致命
        pass
    el = _normalize_element(fallback)
    if el:
        return el, "cli_fallback"
    return None, "not_found"


def _multi_row(element, path, res, outdir) -> dict:
    """把一个元素的 PipelineResult 压成综合报告行（不排名、不推荐）。"""
    structure = getattr(res, "structure", None) or {}
    fit = getattr(res, "fit", None) or {}
    ranked = structure.get("ranked") or []
    candidates = [{"id": r.get("id"), "score": r.get("score")}
                  for r in ranked[:3]]
    row = {
        "absorber": element,
        "file": str(path),
        "output_dir": str(outdir),
        "identify_status": structure.get("status"),
        "candidates": candidates,
        "fit_mode": fit.get("mode"),
        "fit_motif": fit.get("motif"),
        "cn_reliable": fit.get("cn_reliable"),
        "s02_source": fit.get("s02_source"),
        "s02": fit.get("s02"),
        "r_factor": None,
        "first_shell": None,
        "warnings": [],
        "xanes_cross_check": getattr(res, "xanes_cross_check", None),
        "status": "ok",
    }
    first = xafs_fit._first_shell_cn(fit) if fit else None
    if first:
        row["first_shell"] = {
            "label": first.get("label"), "n_fit": first.get("n_fit"),
            "n0": first.get("n0"), "reff": first.get("reff")}
    fit_obj = fit.get("fit") if isinstance(fit, dict) else None
    if fit_obj is not None:
        row["r_factor"] = getattr(fit_obj, "r_factor", None)
    for section in (structure, fit):
        if isinstance(section, dict) and section.get("warning"):
            row["warnings"].append(section["warning"])
    if row["xanes_cross_check"]:
        row["warnings"].extend(
            w["message"] for w in
            (row["xanes_cross_check"].get("warnings") or []))
    return row


def run_multi_element(files, outdir, *, fallback_absorber=None,
                      s02_map=None, run_one=None, verbose=True) -> dict:
    """T7：逐文件自动判定吸收元素，分元素独立处理并输出综合报告。

    ``run_one(file, element, sub_outdir) -> PipelineResult`` 是注入点：
    生产入口传入绑定好 CLI 参数的 ``process_sample``；回归测试可传入
    轻量假实现，避免为多元素路由重复跑 FEFF。

    报告**并列**各元素结果：不含 ``recommended``/``best_element``，
    也不按分数跨元素排序——多元素样品的“哪个元素最重要”必须由
    用户依据科学问题决定，工具不自动选择。
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    s02_map = dict(s02_map or {})
    rows, errors = [], []
    seen = set()
    for f in files:
        element, source = detect_absorber(f, fallback_absorber)
        if element is None:
            errors.append({"file": str(f),
                           "error": "无法判定吸收元素：头信息缺失且未给 --absorber"})
            if verbose:
                print(f"  ❌ {Path(f).name}: 无法判定吸收元素")
            continue
        stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", Path(f).stem)[:60]
        sub = outdir / f"{element}_{stem}"
        # 同名同元素文件不应互相覆盖；出现时加序号。
        base, i = sub, 1
        while sub in seen:
            i += 1
            sub = Path(f"{base}_{i}")
        seen.add(sub)
        if verbose:
            print(f"[multi {len(rows) + len(errors) + 1}/{len(files)}] "
                  f"{element} <- {Path(f).name} ({source})")
        try:
            res = (run_one(f, element, sub) if run_one is not None
                   else process_sample([f], sub, absorber=element,
                                       s02=s02_map.get(element)))
            row = _multi_row(element, f, res, sub)
            row["absorber_source"] = source
            row["s02_override"] = s02_map.get(element)
            rows.append(row)
        except Exception as exc:  # noqa: BLE001 - 单元素失败不中断整批
            errors.append({"file": str(f), "absorber": element,
                           "error": f"{type(exc).__name__}: {exc}"})
            if verbose:
                print(f"  ❌ 失败: {exc}")
    elements = sorted({r["absorber"] for r in rows})
    summary = {
        "mode": "multi_element",
        "n_files": len(files),
        "n_ok": len(rows),
        "n_failed": len(errors),
        "elements": elements,
        "recommended": None,
        "best_element": None,
        "recommendation_policy": (
            "none: 并列各元素结果，不自动选择/推荐 Top-1 元素或相"),
        "s02_overrides": s02_map,
        "rows": rows,
        "errors": errors,
    }
    xafs_io.write_json(outdir / "multi_element_summary.json", summary)
    md = _multi_markdown(summary)
    (outdir / "multi_element_summary.md").write_text(md, encoding="utf-8")
    if verbose:
        print("-" * 66)
        print(f"多元素完成：成功 {len(rows)}/{len(files)}，"
              f"元素 {', '.join(elements) if elements else '—'}")
        print(f"  → {outdir / 'multi_element_summary.json'}")
    return summary


def _multi_markdown(summary: dict) -> str:
    lines = [
        "# 多元素联合报告（T7）",
        "",
        f"- 输入文件：{summary['n_files']}；成功：{summary['n_ok']}；"
        f"失败：{summary['n_failed']}",
        f"- 元素：{', '.join(summary['elements']) or '—'}",
        "- **不自动选择/推荐 Top-1 元素或相**；各元素结果并列，"
        "由用户依据科学问题判断。",
        "",
        "| 元素 | 文件 | 识别状态 | 拟合 motif | 第一壳 | N_fit | N0 | Rf | CN 可靠 | S₀²来源 | 警告 |",
        "|------|------|----------|------------|--------|-------|----|----|---------|----------|------|",
    ]
    for r in summary["rows"]:
        fs = r.get("first_shell") or {}
        def fmt(v, nd=3):
            return "—" if v is None else f"{float(v):.{nd}f}"
        warns = "；".join(r.get("warnings") or []) or "—"
        lines.append(
            f"| {r['absorber']} | {Path(r['file']).name} | "
            f"{r.get('identify_status') or '—'} | {r.get('fit_motif') or '—'} | "
            f"{fs.get('label') or '—'} | {fmt(fs.get('n_fit'))} | "
            f"{fmt(fs.get('n0'))} | {fmt(r.get('r_factor'), 4)} | "
            f"{r.get('cn_reliable')} | {r.get('s02_source') or '—'} | "
            f"{warns} |")
    if summary["errors"]:
        lines += ["", "## 失败文件", ""]
        for e in summary["errors"]:
            lines.append(f"- `{e['file']}`：{e['error']}")
    lines += ["", "> 顺序为输入文件顺序，不代表优先级或推荐。" ]
    return "\n".join(lines) + "\n"


def process_sample(paths, outdir, align: bool = False,
                   absorber: str | None = None,
                   edge_energy: float | None = None,
                   e0_method: str = "derivative", manual_e0: float | None = None,
                   rbkg: float | None = None, kmin: float = 3.0, kmax: float = 12.0,
                   kweight: int = 2, post_window=(150, 600),
                   strict: bool = True,
                   diagnostics: bool = False, identify: bool = False,
                   fit: bool = False, fit_motif: str | None = None,
                   fit_rmin: float | None = None, fit_rmax: float | None = None,
                   allow_ambiguous: bool = False,
                   l1_ligand: str = "O", l1_geometry: str = "oct",
                   fit_weights: bool = False,
                   fit_kweights: tuple | None = None,
                   include_ms: bool = False,
                   despike: bool = False,
                   deglitch: list[tuple[float, float]] | None = None,
                   deadtime: bool = False,
                   reject_outlier_scans: bool = False,
                   lcf_standards: list[str] | None = None,
                   vary_c3: bool | list | None = None,
                   vary_c4: bool | list | None = None,
                   isotropic_deltar: bool = False,
                   shared_sigma2: str | list | None = None,
                   mback: bool = False,
                   chi_std=None, chi_std_k=None,
                   clamp_lo: float | None = None,
                   clamp_hi: float | None = None,
                   fit_space: str = "r",
                   bvs_restraint: dict | None = None,
                   catalytic_model=None,
                   hetero_cocupancy: bool = False,
                   constraints=None,
                   s02: float | None = None,
                   s02_source: str | None = None,
                   phase_correct: bool = False,
                   chi2_contour: bool = False,
                   self_absorption_formula: str | None = None,
                   self_absorption_thickness: float | None = None,
                   self_absorption_density: float | None = None,
                   self_absorption_angles: tuple = (45.0, 45.0),
                   self_absorption_line: str | None = None,
                   verbose: bool = True,
                   formal_cn: bool = False,
                   legacy_auto_fit: bool = False,
                   xanes_cross_check: bool = False,
          # ---- P1b 相解析（正式 CN 结构来源）----
          phase_cif: list[str] | None = None,
          phase_id: str | None = None,
          phase_spec: str | None = None,
          phase_index_path: str | Path | None = None,
          phase_pack_path: str | Path | None = None,
          absorber_site: str | None = None,
          offline: bool = False):
    """
    跑完 S1–S6 并落盘。返回 PipelineResult。

    Rbkg 默认由 optimize_rbkg 自动判定；显式给 rbkg 时跳过寻优。
    diagnostics=True 时在 <outdir>/diagnostics/ 下画各阶段诊断图。
    identify=True 时经 xafs_fit 追加 S6.5 结构识别（拟合链唯一触发点）。
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    result = xafs_pipeline.process_sample(
        paths, align=align, absorber=absorber, edge_energy=edge_energy,
        e0_method=e0_method, manual_e0=manual_e0, rbkg=rbkg,
        kmin=kmin, kmax=kmax, kweight=kweight,
        post_window=post_window, strict=strict,
        despike=despike, deglitch=deglitch, deadtime=deadtime,
        reject_outlier_scans=reject_outlier_scans,
        mback=mback, chi_std=chi_std, chi_std_k=chi_std_k,
        clamp_lo=clamp_lo, clamp_hi=clamp_hi,
        self_absorption_formula=self_absorption_formula,
        self_absorption_thickness=self_absorption_thickness,
        self_absorption_density=self_absorption_density,
        self_absorption_angles=self_absorption_angles,
        self_absorption_line=self_absorption_line,
        verbose=verbose)
    result.outdir = outdir
    result.stem = Path(paths[0]).stem
    # 正式 CN 溯源：输入文件哈希（不含绝对路径与 mtime）。
    result.data_sha256 = cn_table.data_sha256(paths)

    # ---- B2-2 XANES LCF（可选；端元文件先各自走 S1–S3）----
    if lcf_standards:
        std_norms = {}
        for std_path in lcf_standards:
            std_res = xafs_pipeline.process_sample(
                [std_path], absorber=absorber, strict=False, verbose=False)
            std_norms[Path(std_path).stem] = std_res.norm
        lcf_res = xafs_lcf.lcf_from_scans(
            result.norm, std_norms, mode="xanes")
        xafs_io.write_json(outdir / f"{result.stem}_lcf.json",
                           xafs_lcf.lcf_summary(lcf_res))
        if verbose:
            w = ", ".join(f"{k}={v:.3f}" for k, v in
                          lcf_res["weights"].items())
            print(f"  [B2-2] XANES LCF: {w}   Rf={lcf_res['r_factor']:.2e}")

    # P1b phase resolution: parse --phase-* into StructureModel (formal CN strict track)
    structure_model = None
    structure_models: list = []
    gate_result = None
    gate_results: list = []
    if formal_cn and (phase_cif or phase_id or phase_spec or phase_index_path):
        if not absorber:
            raise XAFSConfigError(
                "FIT", "--formal-cn + --phase-* requires --absorber to specify absorber element")
        resolution = phase_resolver.resolve_phase(
            absorber,
            phase_cif=phase_cif,
            phase_id=phase_id,
            phase_spec=phase_spec,
            phase_index_path=phase_index_path,
            phase_pack_path=phase_pack_path,
            absorber_site=absorber_site,
            offline=offline)
        if resolution.status == "ambiguous":
            raise XAFSConfigError(
                "FIT", f"Phase ambiguous: {resolution.reason}; "
                f"please use --phase-id to specify precisely or --absorber-site to disambiguate")
        if resolution.status == "offline_miss":
            raise XAFSConfigError(
                "FIT", f"Phase offline miss: {resolution.reason}; does not fall back to motif")
        structure_model = resolution.structure_model
        if structure_model is None:
            raise XAFSConfigError(
                "FIT", f"Phase resolution did not produce a usable structure model: {resolution.reason}")
        site_models = list(getattr(resolution, "site_models", []) or [])
        multi_site = len(site_models) > 1
        if multi_site:
            structure_models = site_models
            structure_model = site_models[0]
        if verbose:
            sel = resolution.selected
            print(f"  [P1b] phase resolved: {sel.source}:{sel.id} "
                  f"sha256={sel.cif_sha256[:12]}... site={resolution.absorber_site}")

        # P2a CIF gate: validate the resolved CIF and compute N0 shell grouping
        # (Rayleigh). The gate result carries n0_shells needed for P3 conditional CN.
        gate_kwargs = dict(
            kmin=kmin, kmax=kmax, rmax=6.0,
            allow_computed=False,
            allow_unusual_coordination=False,
            is_computed=not resolution.selected.experimental)
        if multi_site:
            # P5-D：每个不等价吸收位独立过门禁；任一位点失败即拒绝，
            # 不做“部分位点成功即放行”。
            gate_results = cif_gate.gate_cif_all_sites(
                resolution.selected.cif_path, absorber, **gate_kwargs)
            failed = [g for g in gate_results if not g.ok]
            if failed:
                labels = ", ".join(
                    f"{g.absorber_site or '?'}({g.code})" for g in failed)
                raise XAFSConfigError(
                    "FIT", f"CIF gate rejected multi-site: {labels} — "
                           f"{failed[0].reason}")
            gate_result = gate_results[0]
            for m, g in zip(structure_models, gate_results):
                m.occupancy_policy = g.occupancy_source
                m.n0_shells = list(g.n0_shells)
                m.warnings = list(g.warnings)
        else:
            gate_result = cif_gate.gate_cif(
                resolution.selected.cif_path, absorber,
                absorber_site=resolution.absorber_site,
                **gate_kwargs)
            if not gate_result.ok:
                raise XAFSConfigError(
                    "FIT", f"CIF gate rejected: {gate_result.code} — {gate_result.reason}")
            if structure_model is not None:
                structure_model.occupancy_policy = gate_result.occupancy_source
                structure_model.n0_shells = list(gate_result.n0_shells)
                structure_model.warnings = list(gate_result.warnings)
        if verbose:
            if multi_site:
                weights = ", ".join(
                    f"{g.absorber_site}={g.site_weight:.3f}"
                    for g in gate_results)
                print(f"  [P2a] CIF gate passed for {len(gate_results)} sites: "
                      f"{weights}")
            else:
                print(f"  [P2a] CIF gate passed: occ_src={gate_result.occupancy_source}, "
                      f"first_shell_cn={gate_result.first_shell_cn}, "
                      f"n0_shells={len(gate_result.n0_shells)}")

    # S6.5/S7 only triggered via the fitting chain (entry -> xafs_fit -> structure_id, linear no cycle)
    multi_site_formal = bool(formal_cn and len(structure_models) > 1)
    if identify or fit or xanes_cross_check:
        xafs_fit.run(
            result, identify=identify or fit, fit=fit,
            absorber_symbol=absorber, fit_motif=fit_motif,
            fit_rmin=fit_rmin, fit_rmax=fit_rmax,
            allow_ambiguous=allow_ambiguous,
            l1_ligand=l1_ligand, l1_geometry=l1_geometry,
            fit_weights=fit_weights,
            shared_sigma2=shared_sigma2, kweights=fit_kweights,
            include_ms=include_ms, vary_c3=vary_c3, vary_c4=vary_c4,
            isotropic_deltar=isotropic_deltar,
            fit_space=fit_space, bvs_restraint=bvs_restraint,
            catalytic_model=catalytic_model,
            hetero_cocupancy=hetero_cocupancy,
            constraints=constraints,
            s02=s02,
            s02_source=s02_source,
            verbose=verbose,
            # 多位点时单一位点的 formal 块会被 P5-D 结果覆盖；先跳过，
            # 避免重复跑一套 FEFF/条件 CN。
            formal_cn=formal_cn and not multi_site_formal,
            structure_model=structure_model,
            gate_result=gate_result,
            xanes_cross_check=xanes_cross_check)

        # P5-D：--absorber-site all 时，用每位点独立的正式 CN 结果
        # 覆盖单一位点输出；只有全部位点 ok 才追加 ALL_WEIGHTED 汇总行。
        if multi_site_formal:
            multi = xafs_fit.run_multi_site_formal_cn(
                result,
                site_models=structure_models,
                gate_results=gate_results,
                s02=s02,
                s02_source=s02_source,
                data_sha256=result.data_sha256,
                verbose=verbose)
            result.fit["formal_cn"] = multi
            result.fit["cn_status"] = multi.get("status", "multi_site_failed")

    # ---- 落盘 ----
    stem, norm, chi_k, ft, wave_map = (result.stem, result.norm, result.chi_k,
                                       result.ft, result.wt)
    best, r1, r_wt = result.rbkg, result.r1, result.r_wt
    xafs_io.write_columns(
        outdir / f"{stem}_mu.dat",
        f"# {stem}   E0={norm.e0:.2f} eV   edge_step={norm.edge_step:.4f}",
        [norm.energy, norm.mu, norm.mu_norm]
        + ([result.merged.sigma_mu] if result.merged.sigma_mu is not None else []),
        ["E(eV)", "mu", "mu_norm"]
        + (["sigma_mu"] if result.merged.sigma_mu is not None else []))
    xafs_io.write_columns(
        outdir / f"{stem}_chi.dat",
        f"# Rbkg={best['rbkg']:.2f} A   kmin/kmax={kmin}/{kmax}   kweight={kweight}",
        [chi_k.k, chi_k.chi, chi_k.k ** kweight * chi_k.chi],
        ["k(A^-1)", "chi", f"k{kweight}chi"])
    xafs_io.write_columns(
        outdir / f"{stem}_ft.dat",
        f"# R(A)  FT kmin={kmin} kmax={kmax}",
        [ft.r, ft.real, ft.imag, ft.magnitude],
        ["R(A)", "Re", "Im", "|chi(R)|"])
    with open(outdir / f"{stem}_wt.dat", "w", encoding="utf-8") as f:
        f.write(f"# CWT sigma={wave_map.sigma}  row=k, col=R\n")
        f.write("# k \\ R  " + "  ".join(f"{v:.2f}" for v in wave_map.r) + "\n")
        for i, kc in enumerate(wave_map.k_centers):
            f.write(f"{kc:.2f}  "
                    + "  ".join(f"{v:.4f}" for v in wave_map.magnitude[i]) + "\n")
    meta = {
        "source": result.merged.meta,
        "absorber": result.absorber_symbol,
        "e0": result.e0,
        "normalization": {
            "edge_step": norm.edge_step, "plateau": norm.plateau,
            "plateau_ok": bool(0.85 <= norm.plateau <= 1.15),
            "strict": bool(strict)},
        "rbkg": {"selected": best["rbkg"], "confidence": best.get("confidence"),
                 "table": result.rbkg_table},
        "ft": {"kmin": kmin, "kmax": kmax, "kweight": kweight,
               "r_first_peak": r1},
        "wt": {"sigma": wave_map.sigma, "r_first_peak": r_wt,
               "k_at_peak": result.k_wt},
        "self_absorption": result.self_absorption,
    }
    if result.structure:
        meta["structure_id"] = {
            "status": result.structure.get("status"),
            "confidence": result.structure.get("confidence"),
            "top3": result.structure["ranked"][:3],
            "indistinguishable": result.structure["indistinguishable"],
            "uncovered_peaks": result.structure["uncovered_peaks"],
            "coverage_ok": result.structure.get("coverage_ok"),
        }
        if result.structure.get("warning"):
            meta["structure_id"]["warning"] = result.structure["warning"]
        if result.structure.get("l1_first_shell"):
            meta["structure_id"]["l1_first_shell"] = \
                result.structure["l1_first_shell"]
        if result.structure.get("l1_error"):
            meta["structure_id"]["l1_error"] = result.structure["l1_error"]
    if result.fit:
        meta["fit"] = _fit_outputs(result.fit, outdir, stem)
        if formal_cn:
            formal = result.fit.get("formal_cn") or {}
            if legacy_auto_fit:
                formal["status"] = "legacy_unverified"
                formal["reject_reason"] = formal.get("reject_reason") or "legacy_auto_fit: 旧版自动拟合，不进正式表"
                formal["legacy_auto_fit"] = True
            rows = list(formal.get("rows") or [])
            model = formal.get("model") or {}
            if not rows:
                reason = formal.get("reject_reason") or "未生成正式 CN"
                rows = [cn_table.refusal_row(
                    formal.get("status", "no_cif"), reason,
                    phase_id=model.get("model_id", ""),
                    model_level=model.get("model_level", ""),
                    source_type=model.get("source_type", ""),
                    source_ref=model.get("source_ref", ""),
                    source_cif_sha256=model.get("source_cif_sha256", ""),
                    absorber_site=model.get("absorber_site", ""),
                    s02=result.fit.get("s02", ""),
                    kmin=kmin, kmax=kmax, kweight=kweight,
                    rmin=result.fit.get("rmin", ""),
                    rmax=result.fit.get("rmax", ""),
                )]
            cn_table.write_cn_table(outdir / "cn_table.csv", rows)
            meta["formal_cn"] = formal
        # B9-1：只有真正放开 C3/C4 才重构非高斯 g(R)，不凭空造分布。
        fit_obj = result.fit.get("fit")
        if fit_obj is not None:
            rdf_tab = xafs_fit.rdf_table(fit_obj)
            if rdf_tab["rows"]:
                xafs_fit.write_rdf_csv(outdir / f"{stem}_rdf.csv", rdf_tab)
                meta["rdf"] = rdf_tab["summary"]
    if result.xanes_cross_check:
        meta["xanes_cross_check"] = result.xanes_cross_check
    xafs_io.write_json(outdir / f"{stem}_meta.json", meta)

    if diagnostics:
        made = xafs_diagnostics.write_diagnostics(result,
                                                  outdir / "diagnostics")
        if verbose:
            print(f"  [诊断] {len(made)} 张阶段图 → {outdir / 'diagnostics'}")

    # ---- B14-2 审稿级诊断：相移校正真键长 / 二维 χ² 等高线 ----
    if phase_correct or chi2_contour:
        diag_dir = outdir / "diagnostics"
        diag_dir.mkdir(parents=True, exist_ok=True)
        payload = result.fit
        if payload is None and result.structure:
            payload = {"l1": result.structure.get("l1_first_shell"),
                       "absorber": result.absorber_symbol}
        if phase_correct:
            pc = xafs_diagnostics.write_phase_corrected_diagnostics(
                result, payload, diag_dir, stem=result.stem)
            if verbose:
                print(f"  [B14-2] 相移校正：R_raw={pc['r_peak_raw']:.3f} → "
                      f"R_true={pc['r_peak_corrected']:.3f} Å "
                      f"（拟合 R={pc['r_fit']:.3f}）")
        if chi2_contour:
            cc = xafs_diagnostics.write_chi2_contour_diagnostics(
                result, payload, diag_dir, stem=result.stem)
            if verbose and cc is not None:
                print(f"  [B14-2] χ² 等高线：min_chi2={cc['min_chi2']:.3f} "
                      f"→ {cc['png']}")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="XAFS 自动处理：μ(E) → 归一化 → χ(k) → |χ(R)| → 小波")
    parser.add_argument("files", nargs="*", help="原始数据文件（可多个扫描）")
    parser.add_argument("-o", "--outdir", default="./xafs_out", help="输出目录")
    parser.add_argument("--align", action="store_true", help="用参考通道做能量对齐")
    parser.add_argument("--strict", action=argparse.BooleanOptionalAction,
                        default=True,
                        help="边后平台不在 0.85–1.15 时是否中断（默认中断；"
                             "--no-strict 降级为警告，供已知厚样/标样使用）")
    parser.add_argument("--diagnostics", action="store_true",
                        help="输出 S1–S6 各阶段诊断图（<outdir>/diagnostics/）")
    parser.add_argument("--phase-correct", action="store_true",
                        help="B14-2 相移校正真键长 FT/WT 诊断"
                             "（默认关闭；需 --identify 或 --fit）")
    parser.add_argument("--chi2-contour", action="store_true",
                        help="B14-2 第一壳 N–σ² 二维 χ² 置信椭圆"
                             "（默认关闭；需 --fit）")
    parser.add_argument("--operando", action="store_true",
                        help="B15-1 原位序列模式：每个文件视为一个时间点"
                             "（至少 3 个），输出 PCA/SPOIL/LCF/MCR-ALS "
                             "演变总表与四宫格图")
    parser.add_argument("--identify", action="store_true",
                        help="追加 S6.5 结构识别（motif 库指纹匹配）")
    parser.add_argument("--fit", action="store_true",
                        help="追加 S7 拟合：confident 时自动多壳层 CN，"
                             "partial/ambiguous/no_match 时 L1 第一壳；"
                             "需要 --absorber 且 S₀² 已标定")
    parser.add_argument("--formal-cn", action="store_true",
                        help="P0a 正式 CN 链：无合格结构来源时写拒绝行到 "
                             "cn_table.csv；诊断拟合照常输出（需 --fit）")
    parser.add_argument("--xanes-cross-check", action="store_true",
                        help="T8 XANES–EXAFS 联合软校验：白线/边步/归一化"
                             "异常时警告；只警告，不改变 cn_reliable/正式 CN"
                             "（可与 --identify/--fit 同用）")
    parser.add_argument("--legacy-auto-fit", action="store_true",
                        help="P4 旧版自动拟合：标记 legacy_unverified，不进正式 CN 表")
    parser.add_argument("--phase-cif", action="append", default=None,
                        help="P1b 用户真实 CIF 路径（可重复，最高优先级）；"
                             "必须可被 gemmi 解析，且不得是现场生成 CIF")
    parser.add_argument("--phase-id", default=None,
                        help="P1b 精确数据库编号 SOURCE:ID，如 COD:1100108；"
                             "本地缓存/phase pack 优先；未命中且非离线才允许 fetcher")
    parser.add_argument("--phase", default=None,
                        help="P1b 相 SPEC 字符串，如 name=fcc Fe;formula=Fe;"
                             "spacegroup=Im-3m;cell=a=2.868,b=2.868,c=2.868,"
                             "alpha=90,beta=90,gamma=90；仅 split(';')/split('=',1)，禁止 eval/exec")
    parser.add_argument("--phase-index", default=None,
                        help="P1b 本地 phase_index.json 路径（默认 _structure_lib/phase_index.json）")
    parser.add_argument("--phase-pack", default=None,
                        help="P1b phase pack 根目录（默认 _structure_lib/_phase_pack）")
    parser.add_argument("--absorber-site", default=None,
                        help="P1b 吸收位标签，如 Fe1 或 frac:(0,0,0)；"
                             "多位点相未显式指定时仍拒绝（absorber_site_ambiguous），"
                             "传 'all' 表示对每个不等价吸收位分别做正式 CN "
                             "并输出化学计量加权汇总（P5-D）")
    parser.add_argument("--offline", action="store_true",
                        help="P1b 离线模式：不允许 fetcher 联网；未命中时报 offline_miss，不退回 motif")
    parser.add_argument("--fit-motif", default=None,
                        help="强制用指定 motif 做 S7（跳过识别歧义）")
    parser.add_argument("--catalytic-model", default=None,
                        help="B13-1 显式催化位点模型：'sac:N4@2.00,axial=O@2.10,"
                             "C@2.90' 或 'dac:Co@2.50,N4@2.00'，或 JSON 字典；"
                             "绕过结构库搜索直接做 S7（需 --fit）")
    parser.add_argument("--s02", type=float, default=None,
                        help="B18-2 运行时覆盖 S₀²（0.4–1.3）：解锁 absorber.py "
                             "未标定元素（Au/Mo/Cr/Zn 等）的 S6.5/S7，"
                             "不修改静态标定表")
    parser.add_argument("--s02-from-t5", action="store_true",
                        help="从 data/s02_calibration/<元素>/calibration.json "
                             "读取 T5 独立标定值；仅接受 status=calibrated "
                             "的元素，未通过元素明确拒绝；显式 --s02/"
                             "--multi-element-s02 优先。不修改静态表")
    parser.add_argument("--hetero-cocupancy", action="store_true",
                        help="B13-2 异核共占位双分支精修：对 motif.cif_file 中"
                             "occ<1 的异核位点拆成逐元素路径并共享 ΔR/σ²"
                             "（需 --fit 与含占据度的 CIF motif）")
    parser.add_argument("--constraints", default=None,
                        help="B14-1 化学计量约束：'N(Fe-O)+N(Fe-Cl)=6.0' 或 "
                             "JSON 字典列表；被约束的 N 不再作为自由变量")
    parser.add_argument("--fit-rmin", type=float, default=None,
                        help="S7 拟合 R 下限（默认按路径自动）")
    parser.add_argument("--fit-rmax", type=float, default=None,
                        help="S7 拟合 R 上限（默认按路径自动）")
    parser.add_argument("--allow-ambiguous", action="store_true",
                        help="识别 ambiguous/partial 时仍用 Top-1 拟合"
                             "（默认只用 L1）")
    parser.add_argument("--fit-weights", action="store_true",
                        help="mixture 多位点模型：放开各组分权重做相分数拟合"
                             "（N 列变为有效简并度 N_crystal·w，不是配位数）")
    parser.add_argument("--l1-ligand", default="O",
                        help="L1 降级的第一壳配体元素（默认 O）")
    parser.add_argument("--l1-geometry", default="oct", choices=["oct", "tet"],
                        help="L1 降级的第一壳几何（默认 oct）")
    parser.add_argument("--absorber", default=None,
                        help="吸收元素符号（如 Fe/Co/Pt）——S6.5 识别与 S₀² 查表必需")
    parser.add_argument("--edge-energy", type=float, default=None,
                        help="参考边能量（eV）；缺省按 --absorber 查表（对齐时需要）")
    parser.add_argument("--e0-method", default="derivative",
                        choices=["derivative", "half_height", "manual"])
    parser.add_argument("--manual-e0", type=float, default=None)
    parser.add_argument("--rbkg", type=float, default=None, help="指定 Rbkg，跳过寻优")
    parser.add_argument("--self-absorption-formula", default=None,
                        help="荧光自吸收校正的样品总化学式（如 Fe3O4）；"
                             "仅荧光模式，且必须同时给厚度/密度")
    parser.add_argument("--self-absorption-thickness", type=float, default=None,
                        help="样品厚度（μm）；Booth 校正必填")
    parser.add_argument("--self-absorption-density", type=float, default=None,
                        help="有效密度（g/cm³）；粉末压片用几何密度，Booth 校正必填")
    parser.add_argument("--self-absorption-angles", default="45,45",
                        help="入射角,出射角（相对样品表面，度；默认 45,45）")
    parser.add_argument("--self-absorption-line", default=None,
                        help="荧光线（默认 K→Ka、L3→La；可显式 Kb/Lb1 等）")
    parser.add_argument("--kmin", type=float, default=3.0)
    parser.add_argument("--kmax", type=float, default=12.0)
    parser.add_argument("--kweight", type=int, default=2)
    parser.add_argument("--shared-sigma2", default=None,
                        help="σ² 共享组（对标 Artemis 路径约束）：'all' 让所有路径共享一个 σ²；"
                             "或 JSON 列表如 '[[0,1],[2]]' 显式分组；默认每条路径独立 σ²")
    parser.add_argument("--fit-kweights", default=None,
                        help="多 k 权重联合拟合（对标 Artemis），如 '1,2,3'；"
                             "默认沿用单权 --kweight（2）")
    parser.add_argument("--include-ms", action="store_true",
                        help="S7 纳入近共线多散射路径（nleg=3/4），参数由母壳层"
                             "几何比例约束，不新增自由变量（默认关闭）")
    parser.add_argument("--vary-c3", default=None,
                        choices=["first", "all"],
                        help="放开三阶累积量 C3（非谐）：first=仅第一壳，"
                             "all=全部壳层；N_idp 不足时自动降级/关闭")
    parser.add_argument("--vary-c4", default=None,
                        choices=["first", "all"],
                        help="放开四阶累积量 C4（非高斯展宽）：first=仅第一壳，"
                             "all=全部壳层；N_idp 不足时自动降级/关闭")
    parser.add_argument("--mback", action="store_true",
                        help="用理论截面（Cromer–Liberman）做 MBACK 归一化，"
                             "适用于边后窗口仅 60–150 eV 的短谱")
    parser.add_argument("--chi-std", default=None,
                        help="AUTOBK 先验标准谱文件（两列 k,chi(k)），"
                             "用于保护低 R 真实信号")
    parser.add_argument("--clamp-lo", type=float, default=None,
                        help="AUTOBK 低 k 端钳制权重（0–1）")
    parser.add_argument("--clamp-hi", type=float, default=None,
                        help="AUTOBK 高 k 端钳制权重（0–1）")
    parser.add_argument("--fit-space", default=None,
                        choices=["r", "q", "wt"],
                        help="S7 拟合空间：r=R 空间（默认，历史口径）/ "
                             "q=R 窗反 FT 的波数空间 / wt=二维小波空间")
    parser.add_argument("--bvs-target", type=float, default=None,
                        help="键价和（BVS）软约束目标，如 Fe³⁺=3.0；"
                             "只增加一行约束、不增加自由变量")
    parser.add_argument("--bvs-sigma", type=float, default=0.1,
                        help="BVS 软约束标准误差（默认 0.1 价单位）")
    parser.add_argument("--bvs-b", type=float, default=0.37,
                        help="BVS 键价参数 b（Å，默认 0.37，Brown & Altermatt）")
    parser.add_argument("--isotropic-deltar", action="store_true",
                        help="各向同性晶格膨胀约束：多壳层 ΔR 合并为唯一全局 "
                             "α（ΔR_i=α·Reff_i），省 M−1 个自由变量")
    parser.add_argument("--config", default=None,
                        help="配置文件（.yaml/.json），CLI 参数优先级高于配置文件")
    parser.add_argument("--despike", action="store_true",
                        help="启用去毛刺（S1.5–S3 之间，中值滑窗+标准差门禁，默认关闭）")
    parser.add_argument("--deglitch", default=None,
                        help="显式能量区间去 glitch，如 '7250:7262,7410:7415'；"
                             "直接删点，不插值（默认关闭）")
    parser.add_argument("--deadtime", action="store_true",
                        help="按 dtime 列做荧光探测器死时间校正（默认关闭）")
    parser.add_argument("--reject-outlier-scans", action="store_true",
                        help="多扫描合并前按中位谱 RMS+MAD 剔除离群扫描，"
                             "至少保留 2 条（默认关闭；默认均值逐位兼容）")
    parser.add_argument("--lcf-standards", default=None,
                        help="XANES LCF 端元文件（逗号分隔），输出 <stem>_lcf.json")
    parser.add_argument("--batch-dir", default=None,
                        help="批处理目录：对该目录下所有匹配的数据文件循环处理，"
                             "失败隔离到 batch_errors.json，与 files 互斥；"
                             "输出目录用 -o 指定（否则取配置 batch.outdir）")
    parser.add_argument("--multi-element", action="store_true",
                        help="T7 多元素模式：逐文件自动判定吸收元素，分元素"
                             "独立跑批并输出 multi_element_summary.json/.md；"
                             "并列结果，不自动选择/推荐 Top-1")
    parser.add_argument("--multi-element-s02", default=None,
                        help="T7 逐元素运行时 S₀² 覆盖，如 'Fe=0.803,"
                             "Co=0.7308'（未列出的元素沿用 --s02/静态表）")
    args = parser.parse_args()

    # ---- 解析 --shared-sigma2（'all' 或 JSON 列表如 '[[0,1],[2]]'）----
    shared_sigma2 = None
    if args.shared_sigma2 is not None:
        s = args.shared_sigma2.strip()
        if s == "all":
            shared_sigma2 = "all"
        else:
            try:
                shared_sigma2 = json.loads(s)
                if not isinstance(shared_sigma2, list):
                    parser.error("--shared-sigma2 必须是 'all' 或 JSON 路径下标列表")
            except json.JSONDecodeError:
                parser.error("--shared-sigma2 格式无效：用 'all' 或 JSON 列表如 '[[0,1],[2]]'")

    try:
        fit_kweights = xafs_config.parse_kweights(args.fit_kweights)
    except ValueError as exc:
        parser.error(str(exc))
    try:
        deglitch = xafs_config.parse_deglitch(args.deglitch)
    except ValueError as exc:
        parser.error(str(exc))
    lcf_standards = None
    if args.lcf_standards:
        lcf_standards = [s.strip() for s in args.lcf_standards.split(",")
                         if s.strip()]
        if not lcf_standards:
            parser.error("--lcf-standards 至少需要一个端元文件")
    vary_c3 = (None if args.vary_c3 is None
               else ([0] if args.vary_c3 == "first" else True))
    vary_c4 = (None if args.vary_c4 is None
               else ([0] if args.vary_c4 == "first" else True))

    # ---- 配置文件加载与 CLI 合并（CLI 显式给出的优先；仅代码默认值时被覆盖）----
    config = xafs_config.load_config(args.config)
    _cli_defaults = {
        "kmin": 3.0, "kmax": 12.0, "kweight": 2,
        "strict": True, "rbkg": None, "wt_sigma": 1.5,
        "self_absorption_formula": None,
        "self_absorption_thickness": None,
        "self_absorption_density": None,
        "self_absorption_angles": "45,45",
        "self_absorption_line": None,
    }
    for _key, _default in _cli_defaults.items():
        if getattr(args, _key, None) == _default:
            if _key.startswith("self_absorption_"):
                _sub = _key[len("self_absorption_"):]
                _cfg_val = config.get("self_absorption", {}).get(_sub)
                if _key == "self_absorption_angles" and _cfg_val:
                    _cfg_val = ",".join(str(a) for a in _cfg_val)
            else:
                _cfg_val = config.get(_key)
            if _cfg_val is not None:
                setattr(args, _key, _cfg_val)
    if not args.despike and config.get("despike"):
        args.despike = True
    if not args.formal_cn and config.get("formal_cn"):
        args.formal_cn = True
    if args.catalytic_model is None and config.get("catalytic_model") is not None:
        args.catalytic_model = config["catalytic_model"]
    if args.s02 is None and config.get("s02") is not None:
        args.s02 = float(config["s02"])
    if args.constraints is None and config.get("constraints") is not None:
        args.constraints = config["constraints"]

    # ---- Phase I–III 开关：CLI 显式 > 配置文件 > 默认 ----
    for _key in ("deadtime", "reject_outlier_scans", "include_ms",
                 "isotropic_deltar", "mback", "hetero_cocupancy",
                 "phase_correct", "chi2_contour", "operando"):
        if not getattr(args, _key, False) and config.get(_key):
            setattr(args, _key, True)
    if args.deglitch is None and config.get("deglitch") is not None:
        args.deglitch = config["deglitch"]
    if args.fit_kweights is None and config.get("fit_kweights") is not None:
        args.fit_kweights = config["fit_kweights"]
    if args.shared_sigma2 is None and config.get("shared_sigma2") is not None:
        args.shared_sigma2 = config["shared_sigma2"]
    if args.lcf_standards is None and config.get("lcf_standards") is not None:
        args.lcf_standards = config["lcf_standards"]
    if args.vary_c3 is None and config.get("vary_c3") is not None:
        args.vary_c3 = config["vary_c3"]
    if args.vary_c4 is None and config.get("vary_c4") is not None:
        args.vary_c4 = config["vary_c4"]
    if args.chi_std is None and config.get("chi_std") is not None:
        args.chi_std = config["chi_std"]
    if args.clamp_lo is None and config.get("clamp_lo") is not None:
        args.clamp_lo = config["clamp_lo"]
    if args.clamp_hi is None and config.get("clamp_hi") is not None:
        args.clamp_hi = config["clamp_hi"]
    if args.fit_space is None:
        args.fit_space = str(config.get("fit_space") or "r").strip().lower()
    if args.fit_space not in ("r", "q", "wt"):
        parser.error(f"--fit-space 只支持 r/q/wt，收到 {args.fit_space!r}")
    try:
        fit_kweights = xafs_config.parse_kweights(args.fit_kweights) \
            if not isinstance(args.fit_kweights, (list, tuple)) \
            else tuple(args.fit_kweights)
    except ValueError as exc:
        parser.error(str(exc))
    try:
        deglitch = xafs_config.parse_deglitch(args.deglitch) \
            if not isinstance(args.deglitch, (list, tuple)) \
            else xafs_config._coerce_config_value("deglitch", args.deglitch)
    except ValueError as exc:
        parser.error(str(exc))
    if isinstance(args.lcf_standards, str):
        lcf_standards = [s.strip() for s in args.lcf_standards.split(",")
                         if s.strip()] or None
    else:
        lcf_standards = args.lcf_standards
    if isinstance(args.shared_sigma2, str) and args.shared_sigma2.strip() == "all":
        shared_sigma2 = "all"
    elif isinstance(args.shared_sigma2, str) and args.shared_sigma2:
        try:
            shared_sigma2 = json.loads(args.shared_sigma2)
        except json.JSONDecodeError:
            parser.error("--shared-sigma2 格式无效：用 'all' 或 JSON 列表如 '[[0,1],[2]]'")
    else:
        shared_sigma2 = args.shared_sigma2
    vary_c3 = xafs_config.normalize_cumulant_switch(args.vary_c3)
    vary_c4 = xafs_config.normalize_cumulant_switch(args.vary_c4)

    # ---- B10-1 chi_std 先验谱：两列 k, chi(k) ----
    chi_std = chi_std_k = None
    if args.chi_std:
        std_path = Path(args.chi_std)
        if not std_path.is_file():
            parser.error(f"--chi-std 文件不存在: {std_path}")
        try:
            std_arr = np.loadtxt(std_path, comments=("#", ";"), ndmin=2)
        except Exception as exc:  # noqa: BLE001 - 输入格式必须可解释
            parser.error(f"--chi-std 读取失败: {exc}")
        if std_arr.shape[1] < 2:
            parser.error("--chi-std 需要至少两列：k(Å⁻¹), chi(k)")
        chi_std_k = np.asarray(std_arr[:, 0], dtype=float)
        chi_std = np.asarray(std_arr[:, 1], dtype=float)
        order = np.argsort(chi_std_k)
        chi_std_k, chi_std = chi_std_k[order], chi_std[order]
    for _name in ("clamp_lo", "clamp_hi"):
        _val = getattr(args, _name)
        if _val is not None and not (0.0 <= float(_val) <= 1.0):
            parser.error(f"--{_name.replace('_', '-')} 必须在 0–1 之间")

    # ---- B11-2 BVS 软约束：CLI 显式 target 覆盖配置文件 ----
    bvs_restraint = config.get("bvs_restraint")
    if bvs_restraint is not None and not isinstance(bvs_restraint, dict):
        parser.error("配置里的 bvs_restraint 必须是字典 "
                     "{target, sigma, r0_table?, b?}")
    if args.bvs_target is not None:
        if not np.isfinite(args.bvs_target):
            parser.error("--bvs-target 必须是有限数")
        if not np.isfinite(args.bvs_sigma) or args.bvs_sigma <= 0:
            parser.error("--bvs-sigma 必须 > 0")
        if not np.isfinite(args.bvs_b) or args.bvs_b <= 0:
            parser.error("--bvs-b 必须 > 0")
        bvs_restraint = {"target": float(args.bvs_target),
                         "sigma": float(args.bvs_sigma),
                         "b": float(args.bvs_b),
                         "absorber": (args.absorber or
                                      (config.get("absorber") or "Fe"))}
    if bvs_restraint is not None and not args.fit:
        parser.error("BVS 软约束只在 S7 拟合中有意义，请同时给 --fit")

    if args.fit_weights and not args.fit:
        parser.error("--fit-weights 需要同时给 --fit（否则不会运行 S7）")
    if args.formal_cn and not args.fit:
        parser.error("--formal-cn 需要同时给 --fit（正式 CN 依赖 S7 拟合结果）")
    if args.s02 is not None and not (0.4 <= args.s02 <= 1.3):
        parser.error("--s02 必须在 [0.4, 1.3] 内")
    if args.s02_from_t5 and not (args.identify or args.fit):
        parser.error("--s02-from-t5 需要 --identify 或 --fit（S₀² 只用于拟合链）")
    for _flag, _label in (("catalytic_model", "--catalytic-model"),
                          ("hetero_cocupancy", "--hetero-cocupancy"),
                          ("constraints", "--constraints")):
        if getattr(args, _flag, None) and not args.fit:
            parser.error(f"{_label} 需要同时给 --fit（否则不会运行 S7）")
    constraints = args.constraints
    if isinstance(constraints, str):
        _c = constraints.strip()
        if _c.startswith("[") or _c.startswith("{"):
            try:
                constraints = json.loads(_c)
            except json.JSONDecodeError as exc:
                parser.error(f"--constraints JSON 解析失败: {exc}")
        else:
            constraints = _c
    if args.phase_correct and not (args.identify or args.fit):
        parser.error("--phase-correct 需要 --identify 或 --fit（要第一壳路径）")
    if args.chi2_contour and not args.fit:
        parser.error("--chi2-contour 需要 --fit（要 S7 最优参数）")
    try:
        sa_angles = tuple(float(x) for x in
                          args.self_absorption_angles.replace(" ", "").split(","))
    except ValueError:
        parser.error("--self-absorption-angles 格式为 入射,出射（度），如 45,45")
    if len(sa_angles) != 2:
        parser.error("--self-absorption-angles 需要两个角度：入射,出射")
    if args.self_absorption_formula:
        if (args.self_absorption_thickness is None
                or args.self_absorption_density is None):
            parser.error("Booth 自吸收校正需要 --self-absorption-thickness "
                         "与 --self-absorption-density")
    elif (args.self_absorption_thickness is not None
          or args.self_absorption_density is not None
          or args.self_absorption_line is not None):
        parser.error("厚度/密度/荧光线只在给 --self-absorption-formula 时使用")
    # 文件来源：--batch-dir 与 files 互斥
    if args.batch_dir and args.files:
        parser.error("--batch-dir 与 files 互斥，请二选一")
    if args.batch_dir:
        batch_dir = Path(args.batch_dir)
        if not batch_dir.is_dir():
            parser.error(f"--batch-dir 不是目录: {batch_dir}")
        pattern = config.get("batch", {}).get("pattern", "*.xdi")
        files = sorted(str(p) for p in batch_dir.glob(pattern))
        if not files:
            parser.error(f"--batch-dir {batch_dir} 下未匹配到 {pattern}")
    elif args.files:
        files = args.files
    elif args.multi_element:
        parser.error("--multi-element 需要 files 或 --batch-dir")
    else:   # 内置演示 = Fe 标样，元素已知
        files = [str(Path(__file__).resolve().parent / "data" / "练手数据" /
                     "02_Fe_standards" / "Fe_lepidocrocite.000")]
        args.absorber = args.absorber or "Fe"

    # ---- T7 多元素布局：逐元素 s02 覆盖表 ----
    multi_s02 = {}
    if args.multi_element_s02:
        for part in args.multi_element_s02.split(","):
            part = part.strip()
            if not part:
                continue
            if "=" not in part:
                parser.error("--multi-element-s02 格式为 'Fe=0.803,Co=0.73'")
            key, val = part.split("=", 1)
            el = _normalize_element(key)
            if el is None:
                parser.error(f"--multi-element-s02 未知元素: {key!r}")
            try:
                v = float(val)
            except ValueError:
                parser.error(f"--multi-element-s02 {el} 的值不是数字: {val!r}")
            if not (0.4 <= v <= 1.3):
                parser.error(f"--multi-element-s02 {el}={v} 必须在 [0.4, 1.3]")
            multi_s02[el] = v

    xafs_pipeline.configure_console()

    # ---- T7 多元素：逐文件自动判定吸收元素，各自独立跑批 ----
    if args.multi_element:
        if args.operando:
            parser.error("--multi-element 与 --operando 不能同时使用")
        if not files:
            parser.error("--multi-element 没有可处理文件")

        def _run_one_element(fpath, element, sub_out):
            # 优先级：--multi-element-s02 > 全局 --s02 > T5 独立标定产物。
            s02_val = multi_s02.get(element)
            s02_src = (f"--multi-element-s02 ({element})"
                       if s02_val is not None else None)
            if s02_val is None and args.s02 is not None:
                s02_val = args.s02
                s02_src = "--s02 runtime override"
            elif s02_val is None and args.s02_from_t5:
                s02_val, s02_src = load_t5_s02(element)
            return process_sample(
                [fpath], sub_out, align=args.align, absorber=element,
                edge_energy=args.edge_energy, e0_method=args.e0_method,
                manual_e0=args.manual_e0, rbkg=args.rbkg, kmin=args.kmin,
                kmax=args.kmax, kweight=args.kweight, strict=args.strict,
                diagnostics=args.diagnostics, identify=args.identify,
                fit=args.fit, fit_motif=args.fit_motif,
                fit_rmin=args.fit_rmin, fit_rmax=args.fit_rmax,
                allow_ambiguous=args.allow_ambiguous,
                l1_ligand=args.l1_ligand, l1_geometry=args.l1_geometry,
                fit_weights=args.fit_weights,
                fit_kweights=fit_kweights, include_ms=args.include_ms,
                despike=args.despike, deglitch=deglitch,
                deadtime=args.deadtime,
                reject_outlier_scans=args.reject_outlier_scans,
                lcf_standards=lcf_standards,
                vary_c3=vary_c3, vary_c4=vary_c4,
                isotropic_deltar=args.isotropic_deltar,
                shared_sigma2=shared_sigma2,
                mback=args.mback, chi_std=chi_std, chi_std_k=chi_std_k,
                clamp_lo=args.clamp_lo, clamp_hi=args.clamp_hi,
                fit_space=args.fit_space, bvs_restraint=bvs_restraint,
                catalytic_model=args.catalytic_model,
                hetero_cocupancy=args.hetero_cocupancy,
                constraints=constraints,
                s02=s02_val,
                s02_source=s02_src,
                phase_correct=args.phase_correct,
                chi2_contour=args.chi2_contour,
                self_absorption_formula=args.self_absorption_formula,
                self_absorption_thickness=args.self_absorption_thickness,
                self_absorption_density=args.self_absorption_density,
                self_absorption_angles=sa_angles,
                self_absorption_line=args.self_absorption_line,
                verbose=False,
                formal_cn=args.formal_cn,
                legacy_auto_fit=args.legacy_auto_fit,
                xanes_cross_check=args.xanes_cross_check,
                phase_cif=args.phase_cif,
                phase_id=args.phase_id,
                phase_spec=args.phase,
                phase_index_path=args.phase_index,
                phase_pack_path=args.phase_pack,
                absorber_site=args.absorber_site,
                offline=args.offline)

        run_multi_element(
            files, args.outdir, fallback_absorber=args.absorber,
            s02_map=multi_s02, run_one=_run_one_element)
        return

    # ---- B15-1 原位序列：每个文件=一个时间点，不合并 ----
    if args.operando:
        if args.batch_dir:
            parser.error("--operando 与 --batch-dir 不能同时使用")
        if len(files) < 3:
            parser.error("--operando 至少需要 3 个样品文件")
        op_s02, op_s02_src = _resolve_s02_for_paths(args, files, args.absorber)
        results = []
        for i, fpath in enumerate(files, 1):
            print(f"[operando {i}/{len(files)}] {Path(fpath).name}")
            results.append(process_sample(
                [fpath], Path(args.outdir) / f"scan_{i - 1:02d}",
                align=args.align, absorber=args.absorber,
                edge_energy=args.edge_energy, e0_method=args.e0_method,
                manual_e0=args.manual_e0, rbkg=args.rbkg, kmin=args.kmin,
                kmax=args.kmax, kweight=args.kweight, strict=args.strict,
                despike=args.despike,
                deadtime=args.deadtime,
                reject_outlier_scans=args.reject_outlier_scans,
                mback=args.mback, chi_std=chi_std, chi_std_k=chi_std_k,
                clamp_lo=args.clamp_lo, clamp_hi=args.clamp_hi,
                self_absorption_formula=args.self_absorption_formula,
                self_absorption_thickness=args.self_absorption_thickness,
                self_absorption_density=args.self_absorption_density,
                self_absorption_angles=sa_angles,
                self_absorption_line=args.self_absorption_line,
                s02=op_s02,
                s02_source=op_s02_src,
                verbose=False))
        std_scans = None
        if lcf_standards:
            std_scans = {}
            for p in lcf_standards:
                std = xafs_pipeline.process_sample(
                    [p], absorber=args.absorber, strict=False, verbose=False)
                std_scans[Path(p).stem] = std.norm
        op = xafs_lcf.run_operando_series(results, standards=std_scans)
        op_dir = Path(args.outdir) / "operando"
        made = xafs_diagnostics.write_operando_diagnostics(
            op, op_dir, stem="operando")
        print(f"原位序列：{len(results)} 条谱，n_components="
              f"{op['n_components']}，accepted="
              f"{op['accepted_standards']}，rejected="
              f"{op['rejected_standards']}")
        for key, path in made.items():
            print(f"  → {path}")
        return

    is_batch = args.batch_dir is not None
    if is_batch:
        # 输出目录优先级：显式 -o > 配置文件 batch.outdir > 默认 xafs_out_batch
        if args.outdir and args.outdir != parser.get_default("outdir"):
            batch_out = Path(args.outdir)
        else:
            batch_out = Path(config.get("batch", {}).get("outdir",
                                                         "xafs_out_batch"))
        batch_out.mkdir(parents=True, exist_ok=True)
        errors = []
        print("=" * 66)
        print(f"批处理 {len(files)} 个文件 → {batch_out}")
        print("=" * 66)
        for i, fpath in enumerate(files, 1):
            stem = Path(fpath).stem
            sub_out = batch_out / stem
            print(f"\n[{i}/{len(files)}] {Path(fpath).name}")
            try:
                file_s02, file_s02_src = _resolve_s02_for_paths(
                    args, [fpath], args.absorber)
                process_sample(
                    [fpath], sub_out, align=args.align, absorber=args.absorber,
                    edge_energy=args.edge_energy, e0_method=args.e0_method,
                    manual_e0=args.manual_e0, rbkg=args.rbkg, kmin=args.kmin,
                    kmax=args.kmax, kweight=args.kweight, strict=args.strict,
                    diagnostics=args.diagnostics, identify=args.identify,
                    fit=args.fit, fit_motif=args.fit_motif,
                    fit_rmin=args.fit_rmin, fit_rmax=args.fit_rmax,
                    allow_ambiguous=args.allow_ambiguous,
                    l1_ligand=args.l1_ligand, l1_geometry=args.l1_geometry,
                    fit_weights=args.fit_weights,
                    fit_kweights=fit_kweights, include_ms=args.include_ms,
                    despike=args.despike, deglitch=deglitch,
                    deadtime=args.deadtime,
                    reject_outlier_scans=args.reject_outlier_scans,
                    lcf_standards=lcf_standards,
                    vary_c3=vary_c3,
                    vary_c4=vary_c4,
                    isotropic_deltar=args.isotropic_deltar,
                    shared_sigma2=shared_sigma2,
                    mback=args.mback, chi_std=chi_std, chi_std_k=chi_std_k,
                    clamp_lo=args.clamp_lo, clamp_hi=args.clamp_hi,
                    catalytic_model=args.catalytic_model,
                    hetero_cocupancy=args.hetero_cocupancy,
                    constraints=constraints,
                    s02=file_s02,
                    s02_source=file_s02_src,
                    phase_correct=args.phase_correct,
                    chi2_contour=args.chi2_contour,
                    self_absorption_formula=args.self_absorption_formula,
                    self_absorption_thickness=args.self_absorption_thickness,
                    self_absorption_density=args.self_absorption_density,
                    self_absorption_angles=sa_angles,
                    self_absorption_line=args.self_absorption_line,
                    verbose=False,
                    formal_cn=args.formal_cn,
                    legacy_auto_fit=args.legacy_auto_fit,
                    xanes_cross_check=args.xanes_cross_check,
                    phase_cif=args.phase_cif,
                    phase_id=args.phase_id,
                    phase_spec=args.phase,
                    phase_index_path=args.phase_index,
                    phase_pack_path=args.phase_pack,
                    absorber_site=args.absorber_site,
                    offline=args.offline)
            except Exception as exc:  # 批处理：单个失败不中断整批
                errors.append({"file": fpath, "error": str(exc)})
                print(f"  ❌ 失败: {exc}")
        # 汇总错误
        err_file = batch_out / "batch_errors.json"
        xafs_io.write_json(err_file, errors)
        print("-" * 66)
        print(f"完成：成功 {len(files) - len(errors)}/{len(files)}, "
              f"失败 {len(errors)} → {err_file}")
        return
    print("=" * 66)
    print(f"处理 {len(files)} 个文件 → {args.outdir}")
    print("=" * 66)
    run_s02, run_s02_src = _resolve_s02_for_paths(args, files, args.absorber)
    result = process_sample(
        files, args.outdir, align=args.align, absorber=args.absorber,
        edge_energy=args.edge_energy, e0_method=args.e0_method,
        manual_e0=args.manual_e0, rbkg=args.rbkg, kmin=args.kmin,
        kmax=args.kmax, kweight=args.kweight, strict=args.strict,
        diagnostics=args.diagnostics, identify=args.identify, fit=args.fit,
        fit_motif=args.fit_motif, fit_rmin=args.fit_rmin,
        fit_rmax=args.fit_rmax, allow_ambiguous=args.allow_ambiguous,
        l1_ligand=args.l1_ligand, l1_geometry=args.l1_geometry,
        fit_weights=args.fit_weights,
        fit_kweights=fit_kweights, include_ms=args.include_ms,
        despike=args.despike, deglitch=deglitch, deadtime=args.deadtime,
        reject_outlier_scans=args.reject_outlier_scans,
        lcf_standards=lcf_standards, vary_c3=vary_c3, vary_c4=vary_c4,
        isotropic_deltar=args.isotropic_deltar,
        shared_sigma2=shared_sigma2,
        mback=args.mback, chi_std=chi_std, chi_std_k=chi_std_k,
        clamp_lo=args.clamp_lo, clamp_hi=args.clamp_hi,
        fit_space=args.fit_space, bvs_restraint=bvs_restraint,
        catalytic_model=args.catalytic_model,
        hetero_cocupancy=args.hetero_cocupancy,
        constraints=constraints,
        s02=run_s02,
        s02_source=run_s02_src,
        phase_correct=args.phase_correct,
        chi2_contour=args.chi2_contour,
        self_absorption_formula=args.self_absorption_formula,
        self_absorption_thickness=args.self_absorption_thickness,
        self_absorption_density=args.self_absorption_density,
        self_absorption_angles=sa_angles,
        self_absorption_line=args.self_absorption_line,
        formal_cn=args.formal_cn,
                    legacy_auto_fit=args.legacy_auto_fit,
                    xanes_cross_check=args.xanes_cross_check,
                    phase_cif=args.phase_cif,
                    phase_id=args.phase_id,
                    phase_spec=args.phase,
                    phase_index_path=args.phase_index,
                    phase_pack_path=args.phase_pack,
                    absorber_site=args.absorber_site,
                    offline=args.offline)

    print("-" * 66)
    for f in sorted(result.outdir.glob(f"{result.stem}*")):
        print(f"  {f.name:<34} {f.stat().st_size / 1024:>8.1f} KB")
    return result


if __name__ == "__main__":
    main()
