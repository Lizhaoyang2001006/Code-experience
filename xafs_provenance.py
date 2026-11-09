"""P6.4：输出台账 —— 参数、门禁、决策轨迹与摘要。

本模块**只记录**，不参与任何物理判定，也不改变 ``cn_table.csv``、
``fit.json`` 或任何门禁行为。核心阈值仍由各物理模块决定；这里只把
“本次运行实际用了什么、命中了哪条门禁、为什么得出这个结论”写成
可复核的 JSON 块。

设计原则：
- 只增不改：现有 ``_meta.json`` 字段保持原样，新块追加在顶层；
- 来源可追溯：每个参数标记 ``cli / config / default / auto``；
- 门禁可复核：统一记录 ``threshold / observed / status / action``；
- 决策可解释：按 S1 → S7 → 正式 CN 的顺序记录关键分支；
- 结果可汇总：``summary`` 给出通过/警告/失败的门禁清单。

与 ``regression_test.py`` 的契约 #132 对应。
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import cif_gate
import structure_id
import xafs_config
import xafs_fit

# ---------------------------------------------------------------
# 阈值镜像
# ---------------------------------------------------------------
# 这些阈值在核心模块里是硬编码的；台账只做镜像，不改变判定。
# 契约 #132 会锁定本文件与核心行为一致。
PLATEAU_LO, PLATEAU_HI = 0.85, 1.15
RF_LIMIT = 0.08
DELTA_E0_LIMIT = 15.0
N_ERR_LIMIT = 1.0
N_ERR_FRACTION_LIMIT = 0.5
SELF_ABS_HIGH = xafs_fit.SELF_ABSORPTION_HIGH_EDGE_STEP
SELF_ABS_WARN = xafs_fit.SELF_ABSORPTION_WARN_EDGE_STEP
CIF_CN_NORMAL = (cif_gate.CN_NORMAL_LO, cif_gate.CN_NORMAL_HI)
CIF_CN_WARN_HI = cif_gate.CN_WARN_HI
CIF_BOND_SLACK = cif_gate.PAIR_BOND_SLACK
CIF_OCC_TOL = cif_gate.OCCUPANCY_TOL
SID_SCORE_LIMIT = structure_id.CONFIDENT_SCORE
SID_INDIST = structure_id.INDIST_TOL
SID_MARGIN = structure_id.CONFIDENT_MARGIN
CN_RATIO_OK = (xafs_fit.CN_RATIO_OK_LO, xafs_fit.CN_RATIO_OK_HI)
CN_RATIO_WARN = (xafs_fit.CN_RATIO_WARN_LO, xafs_fit.CN_RATIO_WARN_HI)

# ---------------------------------------------------------------
# 参数元数据：单位、所在阶段、是否影响正式 CN
# ---------------------------------------------------------------
# 未列出的参数在台账里 unit/stage 为空，不影响“全部参数已记录”。
PARAM_META: dict[str, dict[str, Any]] = {
    # S1 读取与预处理
    "align": {"unit": "", "stage": "S1", "formal": False},
    "batch_dir": {"unit": "", "stage": "S1", "formal": False},
    "deadtime": {"unit": "", "stage": "S1", "formal": False},
    "deglitch": {"unit": "eV", "stage": "S1", "formal": False},
    "despike": {"unit": "", "stage": "S1", "formal": False},
    "reject_outlier_scans": {"unit": "", "stage": "S1", "formal": False},
    "multi_element": {"unit": "", "stage": "S1", "formal": False},
    "multi_element_s02": {"unit": "", "stage": "S1", "formal": True},
    "operando": {"unit": "", "stage": "S1–S6", "formal": False},
    # S2 元素与 E0
    "absorber": {"unit": "", "stage": "S2", "formal": True},
    "edge_energy": {"unit": "eV", "stage": "S2", "formal": True},
    "e0_method": {"unit": "", "stage": "S2", "formal": True},
    "manual_e0": {"unit": "eV", "stage": "S2", "formal": True},
    # S3 归一化与自吸收
    "post_window": {"unit": "eV", "stage": "S3", "formal": True},
    "strict": {"unit": "", "stage": "S3", "formal": True},
    "self_absorption_formula": {"unit": "", "stage": "S3", "formal": True},
    "self_absorption_thickness": {"unit": "μm", "stage": "S3", "formal": True},
    "self_absorption_density": {"unit": "g/cm³", "stage": "S3", "formal": True},
    "self_absorption_angles": {"unit": "deg", "stage": "S3", "formal": True},
    "self_absorption_line": {"unit": "", "stage": "S3", "formal": True},
    # S4/S5 背景与 FT
    "rbkg": {"unit": "Å", "stage": "S4", "formal": True},
    "kmin": {"unit": "Å⁻¹", "stage": "S5", "formal": True},
    "kmax": {"unit": "Å⁻¹", "stage": "S5", "formal": True},
    "kweight": {"unit": "", "stage": "S5", "formal": True},
    "mback": {"unit": "", "stage": "S5", "formal": False},
    "chi_std": {"unit": "", "stage": "S5", "formal": False},
    "clamp_lo": {"unit": "", "stage": "S5", "formal": False},
    "clamp_hi": {"unit": "", "stage": "S5", "formal": False},
    # S6 小波
    "wt_sigma": {"unit": "Å", "stage": "S6", "formal": False},
    # S6.5 结构识别/相解析
    "identify": {"unit": "", "stage": "S6.5", "formal": False},
    "allow_ambiguous": {"unit": "", "stage": "S6.5", "formal": False},
    "fit_motif": {"unit": "", "stage": "S6.5", "formal": True},
    "l1_ligand": {"unit": "", "stage": "S6.5", "formal": False},
    "l1_geometry": {"unit": "", "stage": "S6.5", "formal": False},
    "catalytic_model": {"unit": "", "stage": "S6.5", "formal": False},
    "phase_cif": {"unit": "", "stage": "S6.5", "formal": True},
    "phase_id": {"unit": "", "stage": "S6.5", "formal": True},
    "phase": {"unit": "", "stage": "S6.5", "formal": True},
    "phase_index": {"unit": "", "stage": "S6.5", "formal": True},
    "phase_pack": {"unit": "", "stage": "S6.5", "formal": True},
    "absorber_site": {"unit": "", "stage": "S6.5", "formal": True},
    "offline": {"unit": "", "stage": "S6.5", "formal": True},
    "s02_from_t5": {"unit": "", "stage": "S6.5", "formal": True},
    # S7 拟合
    "fit": {"unit": "", "stage": "S7", "formal": True},
    "fit_rmin": {"unit": "Å", "stage": "S7", "formal": True},
    "fit_rmax": {"unit": "Å", "stage": "S7", "formal": True},
    "fit_weights": {"unit": "", "stage": "S7", "formal": True},
    "fit_kweights": {"unit": "", "stage": "S7", "formal": True},
    "fit_space": {"unit": "", "stage": "S7", "formal": True},
    "include_ms": {"unit": "", "stage": "S7", "formal": True},
    "vary_c3": {"unit": "", "stage": "S7", "formal": True},
    "vary_c4": {"unit": "", "stage": "S7", "formal": True},
    "isotropic_deltar": {"unit": "", "stage": "S7", "formal": True},
    "shared_sigma2": {"unit": "", "stage": "S7", "formal": True},
    "constraints": {"unit": "", "stage": "S7", "formal": True},
    "bvs_target": {"unit": "", "stage": "S7", "formal": True},
    "bvs_sigma": {"unit": "", "stage": "S7", "formal": True},
    "bvs_b": {"unit": "Å", "stage": "S7", "formal": True},
    "hetero_cocupancy": {"unit": "", "stage": "S7", "formal": True},
    "s02": {"unit": "", "stage": "S7", "formal": True},
    # T8 与输出
    "xanes_cross_check": {"unit": "", "stage": "T8", "formal": False},
    "formal_cn": {"unit": "", "stage": "P3", "formal": True},
    "legacy_auto_fit": {"unit": "", "stage": "P3", "formal": True},
    "allow_unusual_coordination": {"unit": "", "stage": "P2a", "formal": True},
    "diagnostics": {"unit": "", "stage": "output", "formal": False},
    "phase_correct": {"unit": "", "stage": "output", "formal": False},
    "chi2_contour": {"unit": "", "stage": "output", "formal": False},
    "lcf_standards": {"unit": "", "stage": "XANES", "formal": False},
    # report 入口
    "demo": {"unit": "", "stage": "output", "formal": False},
    "edge": {"unit": "", "stage": "S2", "formal": False},
    "no_fit": {"unit": "", "stage": "S7", "formal": False},
    "no_html": {"unit": "", "stage": "output", "formal": False},
    "no_window_tuner": {"unit": "", "stage": "output", "formal": False},
    "sample": {"unit": "", "stage": "output", "formal": False},
    "shells": {"unit": "", "stage": "output", "formal": False},
    "stem": {"unit": "", "stage": "output", "formal": False},
    "wt_sweep": {"unit": "", "stage": "output", "formal": False},
    "print_gates": {"unit": "", "stage": "output", "formal": False},
    "print_params": {"unit": "", "stage": "output", "formal": False},
}

_EXCLUDE_FROM_PARAMS = {"files"}


def _jsonable(value: Any) -> Any:
    """把 argparse / numpy / Path 值转成 JSON 可写对象。"""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(v) for v in value]
    # numpy 标量
    if hasattr(value, "item"):
        try:
            return _jsonable(value.item())
        except Exception:  # noqa: BLE001
            pass
    return str(value)


def _get(obj: Any, name: str, default: Any = None) -> Any:
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(name, default)
    return getattr(obj, name, default)


def _status_from_bool(ok: bool, *, strict: bool = False) -> str:
    if ok:
        return "passed"
    return "failed" if strict else "warning"


def _action_from_status(status: str, *, strict: bool = False) -> str:
    if status == "passed":
        return "none"
    if status == "warning":
        return "downgraded_to_warning" if strict else "warning_only"
    if status == "failed":
        return "abort" if strict else "block_quantification"
    return "none"


def _cli_dests(parser, argv: list[str]) -> set[str]:
    """从 argv 反推哪些 dest 是 CLI 显式给出的（支持 --x=y 与 -x）。"""
    if parser is None:
        return set()
    mapping: dict[str, str] = {}
    for action in getattr(parser, "_actions", []):
        for opt in getattr(action, "option_strings", []):
            mapping[opt] = action.dest
    dests: set[str] = set()
    for token in argv:
        if not token.startswith("-") or token == "-":
            continue
        opt = token.split("=", 1)[0]
        if opt in mapping:
            dests.add(mapping[opt])
    return dests


@dataclass
class Provenance:
    """一次运行的有效参数与来源快照。"""

    effective: dict[str, Any] = field(default_factory=dict)
    sources: dict[str, str] = field(default_factory=dict)
    argv: list[str] = field(default_factory=list)
    config_path: str | None = None
    config_keys: list[str] = field(default_factory=list)
    inputs: list[str] = field(default_factory=list)

    @classmethod
    def from_args(cls, args, config: dict | None, argv: list[str],
                  parser=None, inputs=None,
                  overrides: dict | None = None,
                  source_overrides: dict | None = None) -> "Provenance":
        effective: dict[str, Any] = {}
        for key, value in vars(args).items():
            if key in _EXCLUDE_FROM_PARAMS:
                continue
            effective[key] = _jsonable(value)

        cli_dests = _cli_dests(parser, argv)
        cfg_path = getattr(args, "config", None)
        # 只把**用户配置文件里真实出现的键**标成 config；load_config(None)
        # 返回的是 DEFAULTS 副本，不能把默认值误标成配置来源。
        config_keys: list[str] = []
        if cfg_path:
            try:
                raw_cfg = xafs_config._load_text(Path(str(cfg_path)))
                if isinstance(raw_cfg, dict):
                    config_keys = sorted(set(raw_cfg) & set(xafs_config.DEFAULTS))
            except Exception:  # noqa: BLE001 - 台账不应因配置读取失败而中断
                config_keys = []
        sources: dict[str, str] = {}
        for key in effective:
            if key in cli_dests:
                sources[key] = "cli"
            elif key in config_keys:
                sources[key] = "config"
            else:
                sources[key] = "default"

        for key, value in (overrides or {}).items():
            effective[key] = _jsonable(value)
            sources.setdefault(key, "auto")
        for key, value in (source_overrides or {}).items():
            sources[key] = str(value)

        if inputs is None:
            inputs = getattr(args, "files", None) or []
        return cls(
            effective=effective,
            sources=sources,
            argv=[str(x) for x in argv],
            config_path=(str(cfg_path) if cfg_path else None),
            config_keys=config_keys,
            inputs=[str(x) for x in inputs],
        )

    def with_overrides(self, values: dict | None = None,
                       sources: dict | None = None) -> "Provenance":
        eff = dict(self.effective)
        src = dict(self.sources)
        for key, value in (values or {}).items():
            eff[key] = _jsonable(value)
            src.setdefault(key, "auto")
        for key, value in (sources or {}).items():
            src[key] = str(value)
        return replace(self, effective=eff, sources=src)

    def param_block(self, values: dict | None = None,
                    sources: dict | None = None) -> dict:
        prov = self.with_overrides(values, sources)
        return {
            "effective": prov.effective,
            "sources": prov.sources,
            "units": {k: PARAM_META.get(k, {}).get("unit", "")
                      for k in prov.effective},
            "stages": {k: PARAM_META.get(k, {}).get("stage", "")
                       for k in prov.effective},
            "affects_formal_cn": sorted(
                k for k in prov.effective
                if PARAM_META.get(k, {}).get("formal")),
            "cli_argv": prov.argv,
            "inputs": prov.inputs,
            "config_file": prov.config_path,
            "config_keys_applied": prov.config_keys,
        }


# ---------------------------------------------------------------
# 门禁台账
# ---------------------------------------------------------------
def _gate(status: str, *, threshold: Any = None, observed: Any = None,
          action: str = "none", **extra: Any) -> dict:
    payload = {
        "status": status,
        "threshold": threshold,
        "observed": observed,
        "action": action,
    }
    payload.update(extra)
    return payload


def _normalization_gate(norm, strict: bool) -> dict:
    plateau = _get(norm, "plateau")
    if plateau is None or not math.isfinite(float(plateau)):
        return _gate("not_available", threshold=[PLATEAU_LO, PLATEAU_HI],
                     observed=None, action="not_evaluated")
    plateau = float(plateau)
    ok = PLATEAU_LO <= plateau <= PLATEAU_HI
    status = _status_from_bool(ok, strict=strict)
    return _gate(
        status,
        threshold=[PLATEAU_LO, PLATEAU_HI],
        observed=round(plateau, 6),
        action=_action_from_status(status, strict=strict),
        strict=bool(strict),
        location="xafs_pipeline.py:196",
    )


def _self_absorption_gate(self_absorption) -> dict:
    sa = self_absorption or {}
    if not sa:
        return _gate("not_available", observed=None, action="not_evaluated")
    risk = str(sa.get("risk") or "unknown")
    cn_gate = bool(sa.get("cn_gate"))
    if risk in ("none", "low", "corrected"):
        status = "passed"
    elif risk == "possible":
        status = "warning"
    elif risk in ("high", "corrected_unreliable", "unknown"):
        status = "failed" if cn_gate else "warning"
    else:
        status = "warning"
    action = "cn_quantification_blocked" if cn_gate else (
        "warning_only" if status == "warning" else "none")
    return _gate(
        status,
        threshold={"high": SELF_ABS_HIGH, "warn": SELF_ABS_WARN},
        observed={
            "mode": sa.get("mode"),
            "edge_step": sa.get("edge_step"),
            "corrected": bool(sa.get("corrected")),
            "reliable": bool(sa.get("reliable")),
            "risk": risk,
            "correction_factor_median": sa.get("correction_factor_median"),
        },
        action=action,
        cn_gate=cn_gate,
        note=sa.get("note"),
        location="xafs_fit.py:60-116",
    )


def _structure_id_gate(structure) -> dict:
    if not structure:
        return _gate("not_applicable", observed=None, action="not_run")
    conf = structure.get("confidence") or {}
    sid_status = str(structure.get("status") or "unknown")
    if sid_status == "confident":
        status = "passed"
        action = "none"
    elif sid_status == "ambiguous":
        status = "warning"
        action = "require_user_disambiguation"
    elif sid_status == "partial":
        status = "warning"
        action = "model_incomplete"
    elif sid_status == "no_match":
        # 库未命中不是数据失败：S7 会退回 L1 第一壳，正式 CN 仍可走
        # 用户 CIF。这里只警告，避免把“诊断库缺条目”误读成“数据无效”。
        status = "warning"
        action = "library_miss_l1_fallback"
    else:
        status = "warning"
        action = "review_required"
    ranked = structure.get("ranked") or []
    return _gate(
        status,
        threshold={
            "score_limit": SID_SCORE_LIMIT,
            "indist_tol": SID_INDIST,
            "margin_limit": SID_MARGIN,
            "peak_width_min": structure_id.PEAK_WIDTH_MIN,
            "peak_prominence_min": structure_id.PEAK_PROM,
        },
        observed={
            "status": sid_status,
            "best_id": ranked[0].get("id") if ranked else None,
            "best_score": conf.get("best_score"),
            "second_score": conf.get("second_score"),
            "margin": conf.get("margin"),
            "coverage_ok": structure.get("coverage_ok"),
            "uncovered_peaks": structure.get("uncovered_peaks"),
            "indistinguishable": structure.get("indistinguishable"),
        },
        action=action,
        warning=structure.get("warning"),
        location="structure_id.py:338-356",
    )


def _cif_gate(gate_result, gate_results) -> dict:
    results = list(gate_results or [])
    if not results and gate_result is not None:
        results = [gate_result]
    if not results:
        return _gate("not_applicable", observed=None, action="not_run")
    ok = all(bool(_get(g, "ok", False)) for g in results)
    status = "passed" if ok else "failed"
    sites = []
    for g in results:
        sites.append({
            "absorber_site": _get(g, "absorber_site"),
            "ok": bool(_get(g, "ok", False)),
            "code": _get(g, "code"),
            "first_shell_cn": _get(g, "first_shell_cn"),
            "occupancy_source": _get(g, "occupancy_source"),
            "n0_shells": len(_get(g, "n0_shells", []) or []),
            "site_weight": _get(g, "site_weight"),
            "warnings": list(_get(g, "warnings", []) or []),
        })
    return _gate(
        status,
        threshold={
            "cn_normal": list(CIF_CN_NORMAL),
            "cn_warn_hi": CIF_CN_WARN_HI,
            "pair_bond_slack_A": CIF_BOND_SLACK,
            "occupancy_tol_A": CIF_OCC_TOL,
        },
        observed={"ok": ok, "sites": sites},
        action="none" if ok else "reject_formal_cn",
        location="cif_gate.py:616-660",
    )


def _fit_gate(fit_out) -> dict:
    if not fit_out or (fit_out.get("fit") is None and not fit_out.get("l1")):
        return _gate("not_applicable", observed=None,
                     action="not_run", location="xafs_fit.py:1500-1564")
    flags = fit_out.get("flags") or {}
    fit = fit_out.get("fit")
    l1 = fit_out.get("l1") or {}
    if fit is not None:
        paths = list(getattr(fit, "paths", []) or [])
        n_err_bad = any(
            p.n_err is not None and (p.n_err > N_ERR_LIMIT
                                     or p.n_err > N_ERR_FRACTION_LIMIT * p.degeneracy)
            for p in paths)
        sigma2_bad = any(bool(getattr(p, "sigma2_at_bound", False))
                         for p in paths)
        c3_bad = any(bool(getattr(p, "c3_at_bound", False)) for p in paths)
        c4_bad = any(bool(getattr(p, "c4_at_bound", False)) for p in paths)
        r_factor = getattr(fit, "r_factor", None)
        delta_e0 = getattr(fit, "delta_e0", None)
        overfit = bool(getattr(fit, "overfit", False))
        n_independent = getattr(fit, "n_independent", None)
        n_variables = getattr(fit, "n_variables", None)
    else:
        # L1 第一壳降级：没有 motif fit 对象，但仍有 rf / N_err 等
        # 可复核指标；不能因为“没有 fit 对象”就整段标成不适用。
        paths = []
        n_val = l1.get("N")
        n_err = l1.get("N_err")
        n_err_bad = bool(
            n_err is not None and (float(n_err) > N_ERR_LIMIT
                                   or (n_val is not None
                                       and float(n_err) > N_ERR_FRACTION_LIMIT
                                       * float(n_val))))
        sigma2_bad = False
        c3_bad = bool(l1.get("c3_at_bound", False))
        c4_bad = False
        r_factor = l1.get("rf")
        delta_e0 = l1.get("delta_e0")
        overfit = False
        n_independent = l1.get("n_independent")
        n_variables = l1.get("n_variables")
    coverage_ok = bool(fit_out.get("coverage_ok", True))
    sa_gate = bool(flags.get("self_absorption_gate"))

    failed: list[str] = []
    if r_factor is not None and float(r_factor) > RF_LIMIT:
        failed.append("r_factor")
    if delta_e0 is not None and abs(float(delta_e0)) >= DELTA_E0_LIMIT:
        failed.append("delta_e0")
    if n_err_bad:
        failed.append("n_err")
    if sigma2_bad:
        failed.append("sigma2_at_bound")
    if c3_bad:
        failed.append("c3_at_bound")
    if c4_bad:
        failed.append("c4_at_bound")
    if overfit:
        failed.append("overfit")
    if not coverage_ok:
        failed.append("coverage")
    if sa_gate:
        failed.append("self_absorption")

    cn_reliable = fit_out.get("cn_reliable")
    if failed:
        status = "failed"
    elif cn_reliable is None:
        status = "warning"
    elif bool(cn_reliable):
        status = "passed"
    else:
        status = "warning"
    action = ("cn_quantification_blocked"
              if (not bool(cn_reliable) or failed) else "none")
    return _gate(
        status,
        threshold={
            "r_factor": RF_LIMIT,
            "abs_delta_e0_eV": DELTA_E0_LIMIT,
            "n_err": N_ERR_LIMIT,
            "n_err_fraction_of_degeneracy": N_ERR_FRACTION_LIMIT,
        },
        observed={
            "mode": fit_out.get("mode"),
            "r_factor": r_factor,
            "delta_e0": delta_e0,
            "n_err_bad": n_err_bad,
            "sigma2_at_bound": sigma2_bad,
            "c3_at_bound": c3_bad,
            "c4_at_bound": c4_bad,
            "overfit": overfit,
            "coverage_ok": coverage_ok,
            "self_absorption_gate": sa_gate,
            "cn_reliable": cn_reliable,
            "cn_status": fit_out.get("cn_status"),
            "n_independent": n_independent,
            "n_variables": n_variables,
            "l1": (l1 if l1 else None),
        },
        action=action,
        failed_checks=failed,
        flags=flags,
        location="xafs_fit.py:1500-1564",
    )


def _formal_cn_gate(fit_out, formal_block) -> dict:
    if formal_block is None:
        return _gate("not_applicable", observed=None,
                     action="not_run", location="xafs_fit.py:2675-2714")
    status_in = str(formal_block.get("status") or "")
    cond = (fit_out or {}).get("conditional_cn") or {}
    cond_fit = cond.get("conditional_n_fit") or {}
    ftest = cond.get("hamilton_f_test") or {}
    n_status = cond.get("n_status")
    if status_in == "ok":
        if n_status == "ok":
            status = "passed"
            action = "none"
        elif n_status == "warning":
            status = "warning"
            action = "warning_emitted"
        else:
            status = "failed"
            action = "formal_cn_refused"
    elif status_in == "hypothesis_model_only":
        status = "not_applicable"
        action = "diagnostic_only"
    elif status_in == "no_cif":
        status = "failed"
        action = "formal_cn_refused"
    elif status_in in ("legacy_unverified", "multi_site_failed"):
        status = "failed"
        action = "formal_cn_refused"
    else:
        status = "failed"
        action = "formal_cn_refused"
    return _gate(
        status,
        threshold={
            "n_ratio_ok": list(CN_RATIO_OK),
            "n_ratio_warn": list(CN_RATIO_WARN),
            "conditional_rf": RF_LIMIT,
            "hamilton_significant": False,
        },
        observed={
            "formal_status": status_in,
            "n_ratio": cond.get("n_ratio"),
            "n_status": n_status,
            "conditional_rf": cond_fit.get("r_factor"),
            "hamilton_significant": ftest.get("significant"),
            "hamilton_p_value": ftest.get("p_value"),
            "consistency_passed": cond.get("consistency_passed"),
        },
        action=action,
        reject_reason=formal_block.get("reject_reason"),
        location="xafs_fit.py:2694-2714",
    )


def _xanes_cross_check_gate(cross_check) -> dict:
    if not cross_check:
        return _gate("not_applicable", observed=None, action="not_run")
    warnings = list(cross_check.get("warnings") or [])
    return _gate(
        "warning" if warnings else "passed",
        threshold={"warnings": 0},
        observed={"n_warnings": len(warnings),
                  "codes": [w.get("code") for w in warnings]},
        action="warning_only" if warnings else "none",
        warnings=warnings,
        location="xafs_fit.py:1912-1920",
    )


def build_gate_ledger(*, norm=None, strict: bool = True,
                      self_absorption=None, structure=None,
                      gate_result=None, gate_results=None,
                      fit_out=None, formal_block=None,
                      xanes_cross_check=None) -> dict:
    """汇总本次运行命中的全部门禁。"""
    gates = {
        "normalization": _normalization_gate(norm, strict),
        "self_absorption": _self_absorption_gate(self_absorption),
        "structure_id": _structure_id_gate(structure),
        "cif": _cif_gate(gate_result, gate_results),
        "fit": _fit_gate(fit_out),
        "formal_cn": _formal_cn_gate(fit_out, formal_block),
        "xanes_cross_check": _xanes_cross_check_gate(xanes_cross_check),
    }
    gates["summary"] = _summarize_gates(gates)
    return gates


def _summarize_gates(gates: dict) -> dict:
    passed, warned, failed, not_applicable = [], [], [], []
    for name, payload in gates.items():
        if name == "summary" or not isinstance(payload, dict):
            continue
        status = payload.get("status")
        if status == "passed":
            passed.append(name)
        elif status == "warning":
            warned.append(name)
        elif status == "failed":
            failed.append(name)
        else:
            not_applicable.append(name)
    return {
        "passed": passed,
        "warned": warned,
        "failed": failed,
        "not_applicable": not_applicable,
    }


# ---------------------------------------------------------------
# 决策轨迹
# ---------------------------------------------------------------
def build_decision_trace(context: dict | None = None) -> list[dict]:
    """按 S1 → S7 → 正式 CN 记录关键分支。

    ``context`` 只接受原始值（数字/字符串/布尔），避免把大对象写进 JSON。
    """
    c = dict(context or {})
    rows: list[dict] = []

    def add(step: str, decision: str, detail: Any = None,
            status: str = "info") -> None:
        if detail is None:
            return
        rows.append({"step": step, "decision": decision,
                     "detail": _jsonable(detail), "status": status})

    n_scans = c.get("n_scans")
    if n_scans is not None:
        add("S1", "merged_scans",
            {"n_scans": n_scans,
             "n_rejected": c.get("n_rejected_scans", 0),
             "align": c.get("align")})
    if c.get("absorber") is not None or c.get("e0") is not None:
        add("S2", "absorber_and_e0",
            {"absorber": c.get("absorber"),
             "e0_method": c.get("e0_method"),
             "e0_eV": c.get("e0"),
             "edge_energy_eV": c.get("edge_energy")})
    if c.get("plateau") is not None:
        add("S3", "normalization",
            {"plateau": c.get("plateau"),
             "strict": c.get("strict"),
             "post_window_eV": c.get("post_window")},
            status=("passed" if PLATEAU_LO <= float(c["plateau"]) <= PLATEAU_HI
                    else ("failed" if c.get("strict") else "warning")))
    if c.get("rbkg") is not None:
        add("S4", "rbkg",
            {"rbkg_A": c.get("rbkg"),
             "confidence": c.get("rbkg_confidence"),
             "manual": c.get("rbkg_manual")})
    if c.get("kmin") is not None:
        add("S5", "ft_window",
            {"kmin": c.get("kmin"), "kmax": c.get("kmax"),
             "kweight": c.get("kweight"),
             "r_first_peak_A": c.get("r_first_peak")})
    if c.get("wt_sigma") is not None:
        add("S6", "wavelet",
            {"sigma": c.get("wt_sigma"),
             "r_peak_A": c.get("r_wt_peak"),
             "k_peak": c.get("k_wt_peak")})
    if c.get("structure_status") is not None:
        add("S6.5", "structure_identification",
            {"status": c.get("structure_status"),
             "top1": c.get("structure_top1"),
             "margin": c.get("structure_margin"),
             "coverage_ok": c.get("structure_coverage_ok")})
    if c.get("phase_id") is not None or c.get("phase_source") is not None:
        add("S6.5", "phase_resolution",
            {"phase_id": c.get("phase_id"),
             "source": c.get("phase_source"),
             "model_level": c.get("model_level"),
             "absorber_site": c.get("absorber_site")})
    if c.get("catalytic_model") is not None:
        add("S6.5", "catalytic_model",
            {"model": c.get("catalytic_model")})
    if c.get("fit_mode") is not None:
        add("S7", "fit",
            {"mode": c.get("fit_mode"),
             "motif": c.get("fit_motif"),
             "rmin_A": c.get("fit_rmin"),
             "rmax_A": c.get("fit_rmax"),
             "s02": c.get("s02"),
             "s02_source": c.get("s02_source")})
    if c.get("cn_reliable") is not None:
        add("S7", "cn_reliability",
            {"cn_reliable": c.get("cn_reliable"),
             "cn_status": c.get("cn_status"),
             "flags": c.get("fit_flags")},
            status=("passed" if c.get("cn_reliable") else "warning"))
    if c.get("formal_status") is not None:
        add("P3", "formal_cn",
            {"status": c.get("formal_status"),
             "reject_reason": c.get("formal_reject_reason"),
             "n_ratio": c.get("n_ratio")},
            status=("passed" if c.get("formal_status") == "ok"
                    else "warning" if c.get("formal_status") == "hypothesis_model_only"
                    else "failed"))
    if c.get("xanes_cross_check_warnings") is not None:
        add("T8", "xanes_cross_check",
            {"n_warnings": c.get("xanes_cross_check_warnings")},
            status=("warning" if c.get("xanes_cross_check_warnings") else "passed"))
    return rows


def build_operando_meta(prov: Provenance, op: dict, *,
                        overrides: dict | None = None) -> dict:
    """原位序列的参数/门禁台账（PCA / SPOIL / MCR-ALS）。"""
    pca = op.get("pca") or {}
    mcr = op.get("mcr_als") or {}
    cum = list(pca.get("cumulative_variance") or [])
    rejected = list(op.get("rejected_standards") or [])
    pca_gate = {
        "status": "passed" if int(op.get("n_components") or 0) >= 2 else "warning",
        "threshold": {"cumulative_variance": 0.99, "min_components": 2},
        "observed": {
            "n_components": op.get("n_components"),
            "cumulative_variance": (cum[-1] if cum else None),
        },
        "action": ("none" if int(op.get("n_components") or 0) >= 2
                   else "review_component_count"),
        "location": "xafs_lcf.py:244-274",
    }
    spoil_gate = {
        "status": "warning" if rejected else "passed",
        "threshold": {"spoil_good": 3.0, "spoil_poor": 6.0,
                      "reject_above": "2 × spoil_poor"},
        "observed": {"accepted": list(op.get("accepted_standards") or []),
                     "rejected": rejected},
        "action": "warning_only" if rejected else "none",
        "location": "xafs_lcf.py:295-344",
    }
    mcr_rf = mcr.get("r_factor")
    mcr_gate = {
        "status": ("passed" if (mcr_rf is not None and float(mcr_rf) <= RF_LIMIT)
                   else "warning" if mcr_rf is not None else "not_available"),
        "threshold": {"r_factor": RF_LIMIT},
        "observed": {"r_factor": mcr_rf,
                     "n_components": mcr.get("n_components")},
        "action": ("none" if (mcr_rf is not None and float(mcr_rf) <= RF_LIMIT)
                   else "review_mcr_fit"),
        "location": "xafs_lcf.py:348-458",
    }
    gates = {
        "pca_variance": pca_gate,
        "standards_spoil": spoil_gate,
        "mcr_als": mcr_gate,
    }
    gates["summary"] = _summarize_gates(gates)
    return {
        "params": prov.param_block(overrides),
        "gates": gates,
        "decisions": [
            {"step": "S1–S6", "decision": "operando_series",
             "detail": {"n_spectra": len(op.get("names") or []),
                        "mode": op.get("mode"),
                        "x_range": op.get("x_range")},
             "status": "info"},
            {"step": "PCA/MCR", "decision": "blind_decomposition",
             "detail": {"n_components": op.get("n_components"),
                        "accepted_standards": op.get("accepted_standards"),
                        "rejected_standards": rejected},
             "status": "info"},
        ],
        "summary": build_summary(gates, fit_out=None, formal_block=None),
    }


def build_summary(gates: dict, fit_out=None, formal_block=None) -> dict:
    summary = dict(gates.get("summary") or _summarize_gates(gates))
    summary.update({
        "cn_reliable": (fit_out or {}).get("cn_reliable"),
        "cn_status": (fit_out or {}).get("cn_status"),
        "formal_cn_status": (formal_block or {}).get("status"),
        "formal_cn_reject_reason": (formal_block or {}).get("reject_reason"),
    })
    return summary


def build_output_meta(prov: Provenance, *, norm=None, strict: bool = True,
                      self_absorption=None, structure=None,
                      gate_result=None, gate_results=None,
                      fit_out=None, formal_block=None,
                      xanes_cross_check=None,
                      overrides: dict | None = None,
                      source_overrides: dict | None = None,
                      decision_context: dict | None = None) -> dict:
    """生成可直接合并进 ``_meta.json`` 的四个顶层块。"""
    params = prov.param_block(overrides, source_overrides)
    gates = build_gate_ledger(
        norm=norm, strict=strict, self_absorption=self_absorption,
        structure=structure, gate_result=gate_result,
        gate_results=gate_results, fit_out=fit_out,
        formal_block=formal_block, xanes_cross_check=xanes_cross_check)
    return {
        "params": params,
        "gates": gates,
        "decisions": build_decision_trace(decision_context),
        "summary": build_summary(gates, fit_out, formal_block),
    }


# ---------------------------------------------------------------
# 终端格式化
# ---------------------------------------------------------------
def _short(value: Any, limit: int = 52) -> str:
    try:
        text = json.dumps(_jsonable(value), ensure_ascii=False, default=str)
    except Exception:  # noqa: BLE001
        text = str(value)
    return text if len(text) <= limit else text[:limit - 3] + "..."


def format_params(params: dict) -> str:
    eff = params.get("effective") or {}
    src = params.get("sources") or {}
    units = params.get("units") or {}
    stages = params.get("stages") or {}
    lines = [
        "生效参数（全部 CLI / 配置 / 默认 / 自动）",
        f"{'参数':<30} {'生效值':<30} {'来源':<8} {'单位':<8} 阶段",
        "-" * 96,
    ]
    for key in sorted(eff):
        lines.append(
            f"{key:<30} {_short(eff[key], 28):<30} "
            f"{src.get(key, ''):<8} {units.get(key, ''):<8} "
            f"{stages.get(key, '')}")
    lines.append("-" * 96)
    lines.append(f"输入文件: {params.get('inputs') or []}")
    lines.append(f"配置文件: {params.get('config_file') or '（无）'}")
    lines.append(f"配置键: {params.get('config_keys_applied') or []}")
    lines.append(f"命令行: {' '.join(params.get('cli_argv') or [])}")
    return "\n".join(lines)


def format_gates(gates: dict) -> str:
    lines = [
        "置信门禁台账（阈值 / 观测 / 状态 / 动作）",
        f"{'门禁':<22} {'状态':<14} {'动作':<28} 观测 / 阈值",
        "-" * 120,
    ]
    for name in ("normalization", "self_absorption", "structure_id", "cif",
                 "fit", "formal_cn", "xanes_cross_check"):
        payload = gates.get(name) or {}
        lines.append(
            f"{name:<22} {str(payload.get('status', '')):<14} "
            f"{str(payload.get('action', '')):<28} "
            f"{_short(payload.get('observed'), 40)}")
        if payload.get("threshold") is not None:
            lines.append(f"{'':<22} {'':<14} {'':<28} "
                         f"阈值 {_short(payload.get('threshold'), 52)}")
        if payload.get("reject_reason"):
            lines.append(f"{'':<22} {'':<14} {'':<28} "
                         f"拒绝原因 {_short(payload.get('reject_reason'), 52)}")
    summary = gates.get("summary") or {}
    lines.append("-" * 120)
    lines.append(f"通过: {summary.get('passed') or []}")
    lines.append(f"警告: {summary.get('warned') or []}")
    lines.append(f"失败: {summary.get('failed') or []}")
    lines.append(f"不适用: {summary.get('not_applicable') or []}")
    return "\n".join(lines)
