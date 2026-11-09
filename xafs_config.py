# -*- coding: utf-8 -*-
"""
xafs_config.py — 管线配置文件支持（YAML / JSON / INI）
==============================================================

职责：
    把 CLI 频繁重复的处理参数固化到配置文件，CLI 参数优先级高于配置文件。
    支持三种格式，无需额外依赖（YAML 用 Python 标准库解析子集，
    遇到复杂 YAML 自动退回 JSON）。

读取优先级：
    CLI 参数 > 配置文件 > 代码默认值

配置文件示例（xafs.yaml）：
    kmin: 3.0
    kmax: 12.0
    kweight: 2
    post_window: [150, 600]
    rbkg: null
    strict: true
    wt_sigma: 1.5
    self_absorption:
        formula: null
        thickness: null
        density: null
        angles: [45, 45]
        line: null
    despike: false
    batch:
        pattern: "*.xdi"
        outdir: "xafs_out_batch"
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

try:
    import yaml  # type: ignore

    def _load_text(path: Path) -> dict[str, Any]:
        with path.open("r", encoding="utf-8") as fh:
            return yaml.safe_load(fh) or {}

except ImportError:  # pragma: no cover - 无 PyYAML 时退回 JSON
    yaml = None

    def _load_text(path: Path) -> dict[str, Any]:
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)


# ------------------------------------------------------------ 默认值
# Phase III（B9-1）：这里必须与 CLI 开关 100% 镜像。任何新增开关
# 如果只加在 argparse 而不加在这里，配置文件就会被静默忽略。
DEFAULTS: dict[str, Any] = {
    "kmin": 3.0,
    "kmax": 12.0,
    "kweight": 2,
    "post_window": [150, 600],
    "rbkg": None,
    "strict": True,
    "wt_sigma": 1.5,
    "despike": False,
    # ---- Phase I/II ----
    "deglitch": None,
    "deadtime": False,
    "reject_outlier_scans": False,
    "fit_kweights": None,
    "include_ms": False,
    "vary_c3": None,
    "vary_c4": None,
    "isotropic_deltar": False,
    "lcf_standards": None,
    "shared_sigma2": None,
    # ---- Phase III ----
    "mback": False,
    "chi_std": None,
    "clamp_lo": None,
    "clamp_hi": None,
    "fit_space": "r",
    "bvs_restraint": None,
    "target_transform": None,
    "mcr_als": None,
    # ---- Phase IV（B13）----
    "catalytic_model": None,
    "hetero_cocupancy": False,
    "constraints": None,
    "phase_correct": False,
    "chi2_contour": False,
    "operando": False,
    # ---- Phase V（B18-2）----
    "s02": None,
    # P0a：正式 CN 双轨开关；默认关闭，保持诊断输出零退化。
    "formal_cn": False,
    # ---- P1b / P5：独立 S₀²、联合校验与相来源 ----
    "s02_from_t5": False,
    "xanes_cross_check": False,
    "phase_cif": None,
    "phase_id": None,
    "phase": None,
    "phase_index": None,
    "phase_pack": None,
    "absorber_site": None,
    "offline": False,
    "self_absorption": {
        "formula": None,
        "thickness": None,
        "density": None,
        "angles": [45.0, 45.0],
        "line": None,
    },
    "batch": {
        "pattern": "*.xdi",
        "outdir": "xafs_out_batch",
    },
}


def normalize_cumulant_switch(value):
    """把 CLI/YAML/API 的 C3/C4 开关统一成 None / True / [idx...]。"""
    if value is None or value is False:
        return None
    if value is True:
        return True
    if isinstance(value, str):
        low = value.strip().lower()
        if low in ("first", "1", "shell1"):
            return [0]
        if low in ("all", "true", "yes", "on"):
            return True
        if low in ("none", "false", "no", "off", ""):
            return None
    if isinstance(value, (list, tuple, set)):
        return [int(v) for v in value]
    return value


def _coerce_config_value(key: str, value):
    """把 YAML/JSON 里的自然写法转成 API 需要的类型（不猜物理量）。"""
    if value is None:
        return None
    if key == "deglitch":
        if isinstance(value, str):
            return parse_deglitch(value)
        if isinstance(value, (list, tuple)):
            out = []
            for item in value:
                if isinstance(item, str):
                    out.extend(parse_deglitch(item) or [])
                else:
                    lo, hi = float(item[0]), float(item[1])
                    out.append((lo, hi) if hi >= lo else (hi, lo))
            return out or None
    if key == "fit_kweights":
        if isinstance(value, str):
            return parse_kweights(value)
        if isinstance(value, (list, tuple)):
            return tuple(int(v) for v in value)
    if key in ("vary_c3", "vary_c4"):
        if isinstance(value, str):
            low = value.strip().lower()
            if low in ("first", "1", "shell1"):
                return [0]
            if low in ("all", "true", "yes", "on"):
                return True
            if low in ("none", "false", "no", "off", ""):
                return None
        if isinstance(value, bool):
            return value
        if isinstance(value, (list, tuple)):
            return [int(v) for v in value]
    if key == "lcf_standards":
        if isinstance(value, str):
            return [s.strip() for s in value.split(",") if s.strip()]
        if isinstance(value, (list, tuple)):
            return [str(v) for v in value]
    if key == "shared_sigma2":
        if isinstance(value, str) and value.strip().lower() == "all":
            return "all"
        if isinstance(value, str):
            return json.loads(value)
        return value
    if key == "post_window":
        if isinstance(value, str):
            parts = [p for p in re.split(r"[,:;\s]+", value.strip()) if p]
            return tuple(float(p) for p in parts)
        return tuple(float(v) for v in value)
    if key == "bvs_restraint":
        if isinstance(value, str):
            return json.loads(value)
        return value
    if key == "fit_space":
        low = str(value).strip().lower()
        if low not in ("r", "q", "wt"):
            raise ValueError(f"fit_space 只支持 r/q/wt，收到 {value!r}")
        return low
    if key == "phase_cif":
        if isinstance(value, str):
            return [s.strip() for s in value.split(",") if s.strip()]
        if isinstance(value, (list, tuple)):
            return [str(v) for v in value]
    return value


def load_config(path: str | Path | None) -> dict[str, Any]:
    """从文件加载配置，与 DEFAULTS 深度合并；path 为 None 时返回默认值。"""
    if path is None:
        return json.loads(json.dumps(DEFAULTS))
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"配置文件不存在: {p}")
    if p.suffix.lower() in (".yaml", ".yml"):
        if yaml is None:
            raise ImportError(
                f"配置文件 {p} 是 YAML 格式，但未安装 PyYAML；"
                "请改用 .json 或 pip install pyyaml"
            )
        user = _load_text(p)
    elif p.suffix.lower() == ".json":
        user = _load_text(p)
    else:
        raise ValueError(
            f"不支持的配置文件格式: {p.suffix}；仅支持 .yaml / .yml / .json"
        )
    return _deep_merge(DEFAULTS, user)


def parse_kweights(text: str | None) -> tuple | None:
    """解析 ``--fit-kweights 1,2,3`` → (1, 2, 3)；None/空串返回 None。

    去重保序、拒绝负数。单值（"2"）等价于单权模式。
    """
    if text is None:
        return None
    s = str(text).replace(" ", "").strip()
    if not s:
        return None
    try:
        vals = [int(x) for x in s.split(",") if x != ""]
    except ValueError as exc:
        raise ValueError(f"kweights 必须是逗号分隔整数，收到 {text!r}") from exc
    if not vals:
        return None
    if any(v < 0 for v in vals):
        raise ValueError(f"kweights 不能为负: {vals}")
    out = []
    for v in vals:
        if v not in out:
            out.append(v)
    return tuple(out)


def parse_deglitch(text: str | None) -> list[tuple[float, float]] | None:
    """解析 ``--deglitch "7250:7262,7410:7415"`` → [(7250,7262),(7410,7415)]。

    也接受分号分隔。区间自动按小→大排序；格式错误抛 ValueError。
    """
    if text is None:
        return None
    s = str(text).strip()
    if not s:
        return None
    ranges = []
    for part in s.replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise ValueError(
                f"deglitch 区间必须写成 e_lo:e_hi，收到 {part!r}")
        lo_s, hi_s = part.split(":", 1)
        try:
            lo, hi = float(lo_s), float(hi_s)
        except ValueError as exc:
            raise ValueError(
                f"deglitch 区间必须是数字: {part!r}") from exc
        if hi < lo:
            lo, hi = hi, lo
        ranges.append((lo, hi))
    if not ranges:
        return None
    return ranges


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """递归合并两层字典，override 覆盖 base。"""
    out = dict(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


# P1b/P5 入口配置键：两个 CLI（run_xafs / xafs_report）必须完全一致，
# 否则配置文件里的相来源、T5 与离线开关会被静默忽略。
P1B_BOOL_KEYS = ("s02_from_t5", "xanes_cross_check", "offline")
P1B_VALUE_KEYS = ("phase_id", "phase", "phase_index", "phase_pack",
                  "absorber_site")


def apply_p1b_config(args, config: dict[str, Any]):
    """把 P1b/P5 配置键合并进 argparse.Namespace（CLI 显式值优先）。

    ``args`` 里已有的非默认值（布尔 True / 非 None）不会被配置覆盖；
    因此两个入口都可以在解析完 CLI 后无条件调用。
    """
    for key in P1B_BOOL_KEYS:
        if not getattr(args, key, False) and config.get(key):
            setattr(args, key, True)
    for key in P1B_VALUE_KEYS:
        if getattr(args, key, None) is None and config.get(key) is not None:
            setattr(args, key, config[key])
    if (getattr(args, "phase_cif", None) is None
            and config.get("phase_cif") is not None):
        args.phase_cif = _coerce_config_value("phase_cif", config["phase_cif"])
    if (getattr(args, "post_window", None) is None
            and config.get("post_window") is not None):
        args.post_window = _coerce_config_value(
            "post_window", config["post_window"])
    return args


def operando_kwargs_from_config(config: dict[str, Any]) -> dict[str, Any]:
    """把 B10 的 target_transform / mcr_als 配置映射到 run_operando_series。

    DEFAULTS 里的这两个键此前只是"声明"，从未被消费；这里接线后
    ``target_transform.spoil_good/spoil_poor`` 与
    ``mcr_als.n_components/max_iter/closure`` 才会真正生效。
    """
    kw: dict[str, Any] = {}
    tt = config.get("target_transform") or {}
    mcr = config.get("mcr_als") or {}
    if tt.get("spoil_good") is not None:
        kw["spoil_good"] = float(tt["spoil_good"])
    if tt.get("spoil_poor") is not None:
        kw["spoil_poor"] = float(tt["spoil_poor"])
    if mcr.get("n_components") is not None:
        kw["n_components"] = int(mcr["n_components"])
    if mcr.get("max_iter") is not None:
        kw["max_iter"] = int(mcr["max_iter"])
    if mcr.get("closure") is not None:
        kw["closure"] = bool(mcr["closure"])
    return kw


def apply_config_to_kwargs(config: dict[str, Any], kwargs: dict[str, Any],
                           explicit: set[str] | None = None) -> dict[str, Any]:
    """把配置文件值填入 kwargs，但 CLI 显式给出的键永远优先。

    优先级：``explicit``（CLI 显式 dest） > 配置文件 > 代码默认值。

    ``explicit`` 用来区分"用户真的传了 False/None"和"argparse 的默认
    False/None"；缺省时用值等价判断（仍等于代码默认值才允许配置覆盖），
    因此直接调用本函数的脚本也能得到正确行为。
    """
    explicit = set(explicit or ())
    merged = dict(kwargs)
    for key, default in DEFAULTS.items():
        if key in ("self_absorption", "batch") or key in explicit:
            continue
        cfg = config.get(key)
        if cfg is None:
            continue
        if key not in merged or merged[key] == default:
            merged[key] = _coerce_config_value(key, cfg)
    # 自吸收子表（CLI 显式 > 配置 > 默认）
    sa = config.get("self_absorption", {}) or {}
    for key in ("formula", "thickness", "density", "line"):
        dest = f"self_absorption_{key}"
        if dest in explicit:
            continue
        if merged.get(dest) is None and sa.get(key) is not None:
            merged[dest] = sa[key]
    if "self_absorption_angles" not in explicit:
        if (merged.get("self_absorption_angles") in (None, (45.0, 45.0), "45,45")
                and sa.get("angles")):
            merged["self_absorption_angles"] = tuple(float(a) for a in sa["angles"])
    return merged
