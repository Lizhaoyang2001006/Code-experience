# -*- coding: utf-8 -*-
"""build_source_pool.py — S₀² 标准谱源隔离池（P5 遗留项准备）

把 T5 未通过的 5 个元素（Rh/Pd/Ir/Ce/U）的谱源状态固化为
``data/s02_calibration/source_pool/index.json``：

- 只读取本地标定产物与谱源元数据，不运行拟合；
- 不修改 ``absorber.py`` 静态表、不修改任何门禁；
- 隔离是**逻辑隔离**：不移动/删除原始谱，避免破坏既有标定溯源；
- 隔离池中的谱源永远不产生 S₀²，运行期不得用于 ``--s02-from-t5``。

用法::

    python tools/build_source_pool.py
    python tools/build_source_pool.py --check
    python tools/build_source_pool.py --json
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CAL_ROOT = ROOT / "data" / "s02_calibration"
POOL_ROOT = CAL_ROOT / "source_pool"
POOL_INDEX = POOL_ROOT / "index.json"
GATE_SOURCE = ROOT / "calibrate_s02.py"

# 池内容变更时手动更新；不取系统时间，保证产物可复现。
POOL_DATE = "2026-11-09"

ELEMENTS = ("Rh", "Pd", "Ir", "Ce", "U")

# k = sqrt((E - E0) / 3.81)，与项目既有 EXAFS k 换算一致。
K_CONVERSION = 3.81

# 历史隔离证据（2026-11-09 首批观测；旧谱已被替换，文件不再引用）：
# 固定保存，保证 index.json 在没有未入库谱文件的机器上也可复现。
LEGACY_SOURCE_EVIDENCE = {
    "Rh": {"source": "_standards_lib/Rh/Rh_metal_foil.dat",
           "source_sha256": "074792a1b5de859f72921a69f0ef90b0f5fa3a023f8dafdfda1c8204fa7a9d85",
           "points": 137, "e_min": 23169.73, "e_max": 23342.59,
           "edge_e0": 23220.0, "k_max": 5.67, "file_format": "IIT-Lytle .dat",
           "attempt_status": "gated_out",
           "gate_failures": ["wt_no_valid_grid"],
           "failure_reason": "标准谱过短：k_max 不足（IIT-Lytle 仅约 170–195 eV above E0），S5/S6 WT 无有效网格",
           "failure_error": "ValueError: zero-size array to reduction operation maximum which has no identity"},
    "Pd": {"source": "_standards_lib/Pd/Pd_metal_foil.dat",
           "source_sha256": "9f07525aaa52e0ac6ebbc4344d9c715e3c50cea6fd2654729b3909e13d77a21e",
           "points": 195, "e_min": 24275.21, "e_max": 24472.85,
           "edge_e0": 24350.0, "k_max": 5.68, "file_format": "IIT-Lytle .dat",
           "attempt_status": "gated_out",
           "gate_failures": ["wt_no_valid_grid"],
           "failure_reason": "标准谱过短：k_max 不足（IIT-Lytle 仅约 170–195 eV above E0），S5/S6 WT 无有效网格",
           "failure_error": "ValueError: zero-size array to reduction operation maximum which has no identity"},
    "Ir": {"source": "_standards_lib/Ir/Ir_metal_foil.dat",
           "source_sha256": "959da06b723b3342c26b5e349dce2ce1aa6ac3dda9fa62479599e5a0387d341b",
           "points": 272, "e_min": 11142.70, "e_max": 11336.90,
           "edge_e0": 11215.0, "k_max": 5.66, "file_format": "IIT-Lytle .dat",
           "attempt_status": "gated_out",
           "gate_failures": ["wt_no_valid_grid"],
           "failure_reason": "标准谱过短：k_max 不足（IIT-Lytle 仅约 170–195 eV above E0），S5/S6 WT 无有效网格",
           "failure_error": "ValueError: zero-size array to reduction operation maximum which has no identity"},
    "Ce": {"source": "data/bmm_standards/Data/Ce/Ce-L3-CeO2.xdi",
           "source_sha256": "1f679ebd0fb31ed37bca9e8ba125ad06fc0b08554052f8329a2ce7386b8c1e7e",
           "points": 259, "e_min": 5523.405, "e_max": 6152.412,
           "edge_e0": 5726.0, "k_max": 10.58, "file_format": "NIST BMM XDI",
           "attempt_status": "gated_out",
           "gate_failures": ["window_runs_interior", "s02_in_range",
                             "s02_err_ok", "window_spread_ok",
                             "sigma2_not_at_bound"],
           "failure_reason": "门禁未通过: window_runs_interior, s02_in_range, s02_err_ok, window_spread_ok, sigma2_not_at_bound",
           "failure_error": None},
    "U": {"source": "_standards_lib/U/U_metal_foil.dat",
          "source_sha256": "5b31461e25fe6379d30197199d31a59bb44af6e8f4855e3cb2c0b995c9510dfe",
          "points": 455, "e_min": 16966.00, "e_max": 18759.65,
          "edge_e0": 17163.1, "k_max": 20.47, "file_format": "IIT-Lytle .dat",
          "attempt_status": "gated_out",
          "gate_failures": ["window_runs_interior", "s02_in_range",
                            "s02_err_ok", "window_spread_ok",
                            "sigma2_not_at_bound"],
          "failure_reason": "门禁未通过: window_runs_interior, s02_in_range, s02_err_ok, window_spread_ok, sigma2_not_at_bound",
          "failure_error": None},
}

# 替换后通过冻结门禁的谱源观测值（来自 _standards_lib 本地文件头，
# 固定保存以保证 build_pool() 与文件是否在本机无关）。
RESOLVED_EVIDENCE = {
    "Rh": {"points": 514, "e_min": 22694.49, "e_max": 24807.16},
    "Pd": {"points": 520, "e_min": 24156.34, "e_max": 26352.59},
    "Ir": {"points": 645, "e_min": 11015.00, "e_max": 12820.32},
    "Ce": {"points": 327, "e_min": 5523.00, "e_max": 6158.49},
    "U": {"points": 381, "e_min": 16966.34, "e_max": 18170.35},
}

# 人工判定的分类与下一步。user_class 保留 Gemini 建议的原始分类，
# assessment 是本项目依据实际数据给出的评估（可以不同意建议）。
POOL_META = {
    "Rh": {
        "edge": "K",
        "source_license": "IIT-Lytle；许可未明确，不提交",
        "phase_baseline": "Rh fcc（COD 9008482）",
        "user_class": "xanes_truncation_replace_spectrum",
        "classification": "short_spectrum_replace",
        "assessment": "同意：谱长不足，无法进入拟合。",
        "resolution": ("已用 IIT-Lytle RAW/Rh/rhfoil.003（514 点，22694–24807 eV）"
                       "替换；half_height E0 + strict 门禁 5/5 通过，"
                       "S₀²=0.7022±0.0191。"),
        "pool_state": "quarantined",
        "action_class": "replace_source",
        "phase_switch_needed": False,
        "next_actions": [
            "已完成：替换为 RAW 长谱并通过全部冻结门禁",
            "维护：运行期仍只通过 --s02/--s02-from-t5 使用标定值，不改静态表",
        ],
        "replacement_candidates": [
            "Rh K 边 Rh 箔透射/荧光谱（k_max ≥ 10 Å⁻¹）",
            "Rh2O3 corundum（P5-C COD 1010584）配套谱",
        ],
    },
    "Pd": {
        "edge": "K",
        "source_license": "IIT-Lytle；许可未明确，不提交",
        "phase_baseline": "Pd fcc（COD 9008478）",
        "user_class": "xanes_truncation_replace_spectrum",
        "classification": "short_spectrum_replace",
        "assessment": "同意：谱长不足，无法进入拟合。",
        "resolution": ("已用 IIT-Lytle RAW/Pd/pdfl.404（520 点，24156–26353 eV）"
                       "替换；strict=False 记录 plateau=0.8138，"
                       "五项拟合门禁 5/5 通过，S₀²=0.6559±0.0141。"),
        "pool_state": "quarantined",
        "action_class": "replace_source",
        "phase_switch_needed": False,
        "next_actions": [
            "已完成：替换为 RAW 长谱并通过全部冻结门禁",
            "维护：运行期仍只通过 --s02/--s02-from-t5 使用标定值，不改静态表",
        ],
        "replacement_candidates": [
            "Pd K 边 Pd 箔透射/荧光谱（k_max ≥ 10 Å⁻¹）",
            "PdO（P5-C COD 1009031）配套谱",
        ],
    },
    "Ir": {
        "edge": "L3",
        "source_license": "IIT-Lytle；许可未明确，不提交",
        "phase_baseline": "Ir fcc（COD 9008470）",
        "user_class": "parameter_window_tuning_raise_kmin",
        "classification": "short_spectrum_replace",
        "assessment": (
            "不同意‘仅调 k_min 可解决’。当前谱 k_max≈5.66 Å⁻¹，"
            "即使 k_min=2.5 Å⁻¹，Δk≈3.16 Å⁻¹，3 窗口的 N_idp 仍小于 "
            "4 个自由参数；且失败发生在 WT 网格为空，尚未进入参数窗口。"
            "k_min 只能在换到长谱后作为诊断项，不能替代换谱，也不能"
            "用来绕过门禁。"
        ),
        "pool_state": "quarantined",
        "action_class": "replace_source_then_kmin_diagnostic",
        "phase_switch_needed": False,
        "resolution": ("已用 IIT-Lytle RAW/Ir/irfoil.146（645 点，11015–12820 eV）"
                       "替换；kmin=3.0 仅作 L3 白线规避诊断，strict 门禁 5/5 "
                       "通过，S₀²=0.4752±0.0149。"),
        "kmin_tuning_assessment": {
            "suggested": True,
            "sufficient_with_current_source": False,
            "reason": (
                "k_max≈5.66 Å⁻¹；k_min≥2.5 Å⁻¹ 时 Δk≈3.16 Å⁻¹，"
                "N_idp 不足以支撑 s02+dE0+dR+sigma2 四个参数；"
                "当前失败点是 WT 无有效网格，不是白线污染。"
            ),
            "allowed_after_replacement": True,
            "may_relax_gates": False,
        },
        "next_actions": [
            "已完成：替换为 RAW 长谱并通过全部冻结门禁",
            "维护：kmin 敏感性已记录在 calibration.json（kw_sensitivity/runs），"
            "不得用于放宽门禁",
        ],
        "replacement_candidates": [
            "Ir L3 边 Ir 箔透射/荧光谱（k_max ≥ 10 Å⁻¹）",
            "IrO2 rutile（P5-C COD 2101854）配套谱",
        ],
    },
    "Ce": {
        "edge": "L3",
        "source_license": "NIST BMM 公共数据（DOI 10.18434/mds2-4032）",
        "phase_baseline": "CeO2 fluorite（COD 1562989）",
        "user_class": "phase_recalibration_oxide_baseline",
        "classification": "oxide_baseline_reprocess",
        "assessment": (
            "物相基准已经是 CeO2：T5 主谱为 NIST BMM "
            "Ce-L3-CeO2.xdi，motif 为 ceo2_fluorite_cod_1562989，"
            "不是 Ce 金属箔。问题不在‘换 CeO2’，而在归一化/ΔE₀/σ²："
            "三个窗口 ΔE₀=9.78–10.20 eV 贴上限，σ² 压到 1e-4 下界，"
            "无 interior 窗口，S₀² 无唯一解。"
        ),
        "pool_state": "quarantined",
        "action_class": "reprocess_normalization_e0",
        "phase_switch_needed": False,
        "resolution": ("已改用 IIT-Lytle RAW/Ce/ceo2.004（327 点，5523–6158 eV，"
                       "Ce L2 前截断）+ rmin=1.5，strict 门禁 5/5 通过，"
                       "S₀²=0.6279±0.0268；BMM 两条失败谱保留在 attempts 审计链。"),
        "next_actions": [
            "已完成：替换为 RAW 长谱并通过全部冻结门禁",
            "维护：BMM 负例（σ² 贴界/平台 2.271）保留在 calibration.json 供复核",
        ],
        "replacement_candidates": [
            "其它来源的 CeO2 L3 透射谱（平台 0.85–1.15）",
            "NIST BMM Ce-L3-CeO2.xdi 重新归一化/换 E0 方法",
        ],
        "rejected_sources": [
            {
                "source": "data/bmm_standards/Data/Ce/Ce-L3-CeriumDioxide.xdi",
                "sha256": "e4d10395c858d41f56a1cddb173745fe71ca45405d532bd8e9cf603894db6603",
                "reason": "边后平台 2.271，归一化伪迹；严格门禁拒绝",
            }
        ],
    },
    "U": {
        "edge": "L3",
        "source_license": "IIT-Lytle；许可未明确，不提交",
        "phase_baseline": "UO2 fluorite（COD 9009049，标定专用 uo2_fluorite motif）",
        "user_class": "phase_recalibration_oxide_baseline",
        "classification": "oxide_phase_verification",
        "assessment": (
            "同意核实 UO2 基准，但当前本地没有 UO2 标准谱。现有 "
            "U_metal_foil.dat 是 IIT-Lytle RAW 荧光谱，α-U 第一壳模型"
            "（4×U–U 2.7425 Å）拟合 s02 钉在 0.3、ΔE₀≈9.4–10.9 eV、"
            "Rf≈0.55–0.60、N_idp=2.48<4 参数。UO2 结构 CIF 已在 P5-C "
            "入库，但还缺配套 UO2 L3 谱；先核实谱源，再决定是否切换。"
        ),
        "pool_state": "quarantined",
        "action_class": "verify_uo2_phase_switch",
        "phase_switch_needed": False,
        "resolution": ("已完成 UO2 基准切换：IIT-Lytle RAW/U/uo2.379（381 点，"
                       "16966–18170 eV）+ COD 9009049 标定专用 uo2_fluorite "
                       "motif（8×U–O 2.3678 Å），strict 门禁 5/5 通过，"
                       "S₀²=0.7800±0.0229。"),
        "next_actions": [
            "已完成：UO2 基准切换并通过全部冻结门禁",
            "维护：uo2_fluorite 仅存在于 data/s02_calibration/motifs，不进生产 _structure_lib",
        ],
        "replacement_candidates": [
            "UO2 L3 透射/荧光谱（许可明确、k_max ≥ 10 Å⁻¹）",
            "更干净的 α-U 金属箔谱 + 多壳层 α-U 模型（备选，不优先）",
        ],
    },
}

GATE_NAMES = (
    "S02_PASS", "S02_ERR_MAX", "RF_MAX", "E0_MAX", "WINDOW_FRACS",
    "WINDOW_SPREAD_MAX", "DELTAR_LIM_CAL",
)


def _literal_assignments(path: Path) -> dict:
    """用 AST 读取模块级字面量常量，不执行被导入模块。"""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    found = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name) or target.id not in GATE_NAMES:
            continue
        try:
            found[target.id] = ast.literal_eval(node.value)
        except (ValueError, SyntaxError):
            continue
    return found


def gate_policy() -> dict:
    raw = _literal_assignments(GATE_SOURCE)
    missing = [name for name in GATE_NAMES if name not in raw]
    if missing:
        raise SystemExit(f"calibrate_s02.py 缺少门禁常量: {missing}")
    return {
        "s02_pass_range": list(raw["S02_PASS"]),
        "s02_err_max": float(raw["S02_ERR_MAX"]),
        "rf_max": float(raw["RF_MAX"]),
        "abs_delta_e0_max_ev": float(raw["E0_MAX"]),
        "window_fracs": list(raw["WINDOW_FRACS"]),
        "window_spread_max": float(raw["WINDOW_SPREAD_MAX"]),
        "delta_r_limits": list(raw["DELTAR_LIM_CAL"]),
        "sigma2_optimizer_floor": 0.0001,
        "sigma2_at_bound_threshold": 0.0005,
        "require_vars_lt_n_idp": True,
        "require_three_interior_windows": True,
        "gate_relaxation": "forbidden",
        "manual_window_trimming": "forbidden",
        "auto_top1_selection": "forbidden",
        "static_absorber_py": "unchanged",
        "runtime_effect": "none",
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _normalize_path(raw) -> str:
    text = str(raw or "").replace("\\", "/")
    while text.startswith("./"):
        text = text[2:]
    return text


def _first_attempt(calibration: dict) -> dict:
    attempts = calibration.get("attempts") or []
    return attempts[0] if attempts else {}


def _gate_failures(attempt: dict) -> list[str]:
    result = attempt.get("result") or {}
    checks = result.get("checks") or {}
    failures = [key for key, value in checks.items() if value is False]
    if failures:
        return failures
    if attempt.get("status") == "error":
        return ["source_or_normalization_error"]
    if attempt.get("error"):
        return ["wt_no_valid_grid"]
    return []


def build_element(element: str) -> dict:
    cal_path = CAL_ROOT / element / "calibration.json"
    calibration = _load_json(cal_path)
    attempts = calibration.get("attempts") or []
    winner = next((a for a in attempts if a.get("status") == "passed"), None)
    meta = POOL_META[element]

    legacy = dict(LEGACY_SOURCE_EVIDENCE[element])
    legacy["k_max_estimate"] = legacy.pop("k_max")

    row = {
        "element": element,
        "edge": meta["edge"],
        "calibration_artifact":
            f"data/s02_calibration/{element}/calibration.json",
        "last_status": calibration.get("status"),
        "classification": meta["classification"],
        "user_suggested_class": meta["user_class"],
        "assessment": meta["assessment"],
        "resolution": meta.get("resolution"),
        "action_class": meta["action_class"],
        "phase_baseline": meta["phase_baseline"],
        "phase_switch_needed": meta["phase_switch_needed"],
        "next_actions": list(meta["next_actions"]),
        "replacement_candidates": list(meta["replacement_candidates"]),
        "quarantined_legacy_source": legacy,
        "acceptance": (
            "替换谱必须通过冻结门禁：三窗口 interior、S₀²∈[0.4,1.15]、"
            "s02_err<0.06、Rf≤0.12、|ΔE₀|<10 eV、窗口极差≤0.05、"
            "σ² 不贴界、vars<N_idp；任一不满足则保持 not_calibrated"
        ),
    }
    if "kmin_tuning_assessment" in meta:
        row["kmin_tuning_assessment"] = meta["kmin_tuning_assessment"]
    if "rejected_sources" in meta:
        row["rejected_sources"] = meta["rejected_sources"]

    if winner is not None:
        # 已解决：所有数值都来自已提交的 calibration.json + 固定观测常量，
        # 不依赖 _standards_lib 原始谱是否在本机。
        result = winner.get("result") or {}
        provenance = winner.get("provenance") or {}
        evidence = RESOLVED_EVIDENCE[element]
        e0 = result.get("e0")
        if e0 is None:
            e0 = winner.get("manual_e0") or legacy["edge_e0"]
        post_edge_span = float(evidence["e_max"]) - float(e0)
        k_max = (math.sqrt(post_edge_span / K_CONVERSION)
                 if post_edge_span > 0 else 0.0)
        source = _normalize_path(
            provenance.get("source_path") or winner.get("standard"))
        source_sha = str(provenance.get("used_sha256")
                         or provenance.get("source_sha256") or "").lower()
        resolved_by = {
            "standard": source,
            "source_sha256": source_sha,
            "note": winner.get("note", ""),
            "motif": winner.get("motif"),
            "motif_file": winner.get("motif_file"),
            "kmin": winner.get("kmin"),
            "kmax": winner.get("kmax"),
            "rmin": winner.get("rmin"),
            "e0_method": winner.get("e0_method"),
            "s02": calibration.get("s02"),
            "s02_err": calibration.get("s02_err"),
            "rf": calibration.get("rf"),
            "spread": calibration.get("spread"),
            "gate_failures": [],
        }
        row.update({
            "e0_ev": round(float(e0), 3),
            "current_source": source,
            "source_sha256": source_sha,
            "source_license": meta["source_license"],
            "source_committed": bool(provenance.get("committed")),
            "data_points": int(evidence["points"]),
            "energy_range_ev": [round(float(evidence["e_min"]), 3),
                                round(float(evidence["e_max"]), 3)],
            "post_edge_span_ev": round(post_edge_span, 3),
            "k_max_estimate": round(k_max, 2),
            "calibration_motif": winner.get("motif"),
            "attempt_status": winner.get("status"),
            "failure_reason": None,
            "failure_error": None,
            "gate_failures": [],
            "s02": calibration.get("s02"),
            "s02_err": calibration.get("s02_err"),
            "rf": calibration.get("rf"),
            "spread": calibration.get("spread"),
            "delta_e0": calibration.get("delta_e0"),
            "sigma2": calibration.get("sigma2"),
            "plateau": calibration.get("plateau"),
            "k_window": calibration.get("k_window"),
            "e0_method": calibration.get("e0_method"),
            "n_variables": calibration.get("n_variables"),
            "n_independent": calibration.get("n_independent"),
            "usable_for_s02": True,
            "pool_state": "resolved",
            "resolved_by": resolved_by,
        })
        return row

    # ---- 未解决分支（保留给未来回归；当前 5 元素均已 resolved） ----
    attempt = attempts[0] if attempts else {}
    result = attempt.get("result") or {}
    provenance = attempt.get("provenance") or {}
    evidence = LEGACY_SOURCE_EVIDENCE[element]

    source = _normalize_path(
        provenance.get("source_path") or attempt.get("standard"))
    source_sha = str(provenance.get("source_sha256") or "").lower()
    source_path = ROOT / source
    if source_path.is_file() and source_sha:
        actual_sha = _sha256(source_path)
        if actual_sha != source_sha:
            raise SystemExit(
                f"{element} 源文件 sha256 不匹配: {source} 实际 {actual_sha}")

    e0 = result.get("e0")
    if e0 is None:
        e0 = attempt.get("manual_e0") or evidence["edge_e0"]
    post_edge_span = float(evidence["e_max"]) - float(e0)
    k_max = (math.sqrt(post_edge_span / K_CONVERSION)
             if post_edge_span > 0 else 0.0)

    row.update({
        "e0_ev": round(float(e0), 3),
        "current_source": source,
        "source_sha256": source_sha,
        "source_license": meta["source_license"],
        "source_committed": bool(provenance.get("committed")),
        "data_points": int(evidence["points"]),
        "energy_range_ev": [round(float(evidence["e_min"]), 3),
                            round(float(evidence["e_max"]), 3)],
        "post_edge_span_ev": round(post_edge_span, 3),
        "k_max_estimate": round(k_max, 2),
        "calibration_motif": attempt.get("motif"),
        "attempt_status": attempt.get("status"),
        "failure_reason": calibration.get("reason"),
        "failure_error": attempt.get("error"),
        "gate_failures": _gate_failures(attempt),
        "s02": calibration.get("s02"),
        "usable_for_s02": False,
        "pool_state": "quarantined",
    })
    return row


def build_pool() -> dict:
    elements = [build_element(element) for element in ELEMENTS]
    counts: dict[str, int] = {}
    for row in elements:
        counts[row["classification"]] = counts.get(row["classification"], 0) + 1
    state_counts: dict[str, int] = {}
    for row in elements:
        state_counts[row["pool_state"]] = (
            state_counts.get(row["pool_state"], 0) + 1)
    return {
        "schema_version": 2,
        "pool_id": "s02-standard-source-isolation",
        "generated_at": POOL_DATE,
        "generated_by": "tools/build_source_pool.py",
        "scope": (
            "T5 5 元素谱源隔离池（历史审计 + 替换准入）；"
            "替换谱通过冻结门禁后标记 resolved 并记录 resolved_by，"
            "旧谱保留在 quarantined_legacy_source；不产生新 S₀²、"
            "不改静态表、不改门禁"
        ),
        "policy": gate_policy(),
        "user_classification_map": {
            "xanes_truncation_replace_spectrum": ["Rh", "Pd"],
            "phase_recalibration_oxide_baseline": ["Ce", "U"],
            "parameter_window_tuning": ["Ir"],
        },
        "elements_total": len(elements),
        "classification_counts": dict(sorted(counts.items())),
        "pool_state_counts": dict(sorted(state_counts.items())),
        "elements": elements,
    }


def _write_json_if_changed(path: Path, data: dict) -> bool:
    payload = (json.dumps(data, ensure_ascii=False, indent=2) + "\n").encode(
        "utf-8")
    if path.is_file() and path.read_bytes() == payload:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        fh.write(payload)
    return True


def _print_summary(data: dict) -> None:
    print(f"pool: {data['pool_id']}  elements={data['elements_total']}")
    for row in data["elements"]:
        print(f"  {row['element']:3s} {row['classification']:30s} "
              f"k_max={row['k_max_estimate']:5.2f}  "
              f"status={row['last_status']}  "
              f"usable={row['usable_for_s02']}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="生成/校验 S₀² 标准谱源隔离池")
    parser.add_argument("--check", action="store_true",
                        help="只校验现有 index.json 是否最新")
    parser.add_argument("--json", action="store_true",
                        help="把池索引输出到 stdout")
    args = parser.parse_args(argv)

    data = build_pool()
    if args.json:
        print(json.dumps(data, ensure_ascii=False, indent=2))
        return 0
    if args.check:
        if not POOL_INDEX.is_file():
            print(f"缺少 {POOL_INDEX.relative_to(ROOT)}", file=sys.stderr)
            return 1
        stored = _load_json(POOL_INDEX)
        if stored != data:
            print("source_pool/index.json 已过期；请运行 "
                  "python tools/build_source_pool.py", file=sys.stderr)
            return 1
        _print_summary(data)
        print("OK: source pool index 与标定产物/门禁常量一致")
        return 0

    changed = _write_json_if_changed(POOL_INDEX, data)
    _print_summary(data)
    print(f"{'更新' if changed else '未变化'}: "
          f"{POOL_INDEX.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    raise SystemExit(main())
