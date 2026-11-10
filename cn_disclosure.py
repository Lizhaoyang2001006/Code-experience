# -*- coding: utf-8 -*-
"""cn_disclosure.py — CN 输出披露标准检查（项目级）

本模块把 XAFS CN 报告中总结出的通用披露要求固化成可执行检查：

- 结构与 CIF 来源必须写清；没有 CIF 时显式声明 ``cif=None``；
- 构造团簇必须给出构造方法，不能用无关 CIF 顶替；
- 必须披露 S0^2、Delta E0、k/R 窗口、kweight、N_idp、变量数、Rf；
- 必须给出 CN 状态与系统误差预算；
- 诊断级 CN 不能伪装成正式 CN；
- ``|Delta E0| > 5 eV`` 只给非阻断提示：检查 FEFF 路径零点、实验 E0
  约定和模型完整性，不要直接当成吸收边位移；不改变披露检查的通过/失败。

检查器只读 fit JSON，不重新拟合、不修改任何数值，也不猜测结构。
它可以独立运行，也可以被 ``paper_figures.py`` 调用。

用法::

    python cn_disclosure.py fit.json
    python cn_disclosure.py fit1.json fit2.json --warn-only

返回码：发现问题返回 1；``--warn-only`` 时始终返回 0。
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

# 这些字段可以单独缺失，但至少要有一种结构描述。
_STRUCTURE_TEXT_FIELDS = (
    "source_type", "source_ref", "description", "short_label", "phase_id",
)

# fit JSON 里必须能读到的全局拟合字段。
_FIT_FIELDS = (
    "s02", "delta_e0", "kmin", "kmax", "kweight", "rmin", "rmax",
    "r_factor",
)

# 每条路径必须能读到的数值字段。
_PATH_FIELDS = ("N", "R", "sigma2")

# |Delta E0| 超过这个值只给非阻断提示，不改变披露门禁。
# 5 eV 是 EXAFS 拟合的经验舒适区；项目硬门禁仍是 |Delta E0| < 15 eV。
DELTA_E0_COMFORT_EV = 5.0


def _cif_path(structure: dict) -> str | None:
    for key in ("cif", "cif_path", "cif_file", "source_cif"):
        value = structure.get(key)
        if value not in (None, ""):
            return str(value)
    return None


def cif_status(structure: dict | None) -> str:
    """返回 ``provided`` / ``hash_only`` / ``none`` / ``not_stated``。

    ``cif=None`` 是显式声明“没有使用 CIF”，不是缺字段；这与
    ``paper_figures.structure_provenance(cif=None)`` 的约定一致。
    """
    if not isinstance(structure, dict):
        return "not_stated"
    s = structure
    if _cif_path(s):
        return "provided"
    if s.get("cif_sha256") or s.get("source_cif_sha256"):
        return "hash_only"
    if s.get("cif_declared") or s.get("_cif_declared"):
        return "none"
    if "cif" in s and s.get("cif") in (None, ""):
        return "none"
    return "not_stated"


def structure_issues(structure: dict | None, *, motif: Any = None) -> list[str]:
    """检查结构/CIF provenance；缺失时明确说“未声明”，不猜 CIF。"""
    s = structure or {}
    if not isinstance(s, dict):
        return [f"structure provenance is not an object: {s!r}"]
    if not s:
        if motif:
            motif_label = motif.get("id") if isinstance(motif, dict) else motif
            return [f"structure/CIF provenance not stated (motif={motif_label})"]
        return ["structure/CIF provenance not stated"]

    issues: list[str] = []
    if not any(s.get(key) not in (None, "") for key in _STRUCTURE_TEXT_FIELDS):
        issues.append(
            "structure provenance lacks source_type/source_ref/description")

    status = cif_status(s)
    if status == "not_stated":
        issues.append(
            "CIF status not stated; add a CIF path + sha256 or explicit cif=None")
    elif status == "provided" and not (
            s.get("cif_sha256") or s.get("source_cif_sha256")):
        issues.append("CIF path present but CIF sha256 missing")

    if not (s.get("feff_inputs") or s.get("feff_input")):
        issues.append("FEFF input files not stated")
    return issues


def _get(payload: dict, f: dict, key: str, default=None):
    value = payload.get(key)
    if value is None:
        value = f.get(key)
    return default if value is None else value


def check_fit_payload(payload: dict) -> list[str]:
    """检查一份 fit JSON（``run_xafs`` / ``paper_figures`` 风格）。

    返回人类可读的问题列表；空列表表示达到最低披露要求。
    本函数不读取外部文件，只检查 JSON 里是否写清了来源与参数。
    """
    if not isinstance(payload, dict):
        return ["fit payload is not a JSON object"]

    f = payload.get("fit")
    if not isinstance(f, dict):
        f = payload if "paths" in payload else {}

    structure = payload.get("structure") or f.get("structure") or {}
    motif = payload.get("motif") or f.get("motif")
    issues = structure_issues(structure, motif=motif)

    paths = f.get("paths") or payload.get("paths") or []
    if not isinstance(paths, (list, tuple)):
        issues.append("fit paths is not a list")
        paths = []
    if not paths:
        issues.append("fit paths missing")
    path_files = payload.get("path_files") or f.get("path_files") or []
    if isinstance(path_files, (str, Path)):
        path_files = [str(path_files)]
    if len(path_files) < len(paths):
        issues.append(
            f"path_files incomplete ({len(path_files)} for {len(paths)} paths)")

    s02 = _get(payload, f, "s02")
    if s02 is None:
        issues.append("S0^2 missing")
    else:
        try:
            if not 0.0 < float(s02) <= 1.5:
                issues.append(f"S0^2 outside plausible range: {s02!r}")
        except (TypeError, ValueError):
            issues.append(f"S0^2 is not numeric: {s02!r}")

    for key in _FIT_FIELDS:
        if key in ("s02", "kweight"):
            continue
        if _get(payload, f, key) is None:
            issues.append(f"{key} missing")

    kweight = _get(payload, f, "kweight")
    if kweight is None:
        kws = f.get("kweights") or payload.get("kweights")
        if isinstance(kws, (list, tuple)) and kws:
            kweight = kws[0]
    if kweight is None:
        issues.append("kweight missing (or kweights list)")

    kmin = _get(payload, f, "kmin")
    kmax = _get(payload, f, "kmax")
    if kmin is not None and kmax is not None:
        try:
            if float(kmin) >= float(kmax):
                issues.append("kmin must be smaller than kmax")
        except (TypeError, ValueError):
            issues.append("kmin/kmax are not numeric")

    rmin = _get(payload, f, "rmin")
    rmax = _get(payload, f, "rmax")
    if rmin is not None and rmax is not None:
        try:
            if float(rmin) >= float(rmax):
                issues.append("rmin must be smaller than rmax")
        except (TypeError, ValueError):
            issues.append("rmin/rmax are not numeric")

    cn_status = _get(payload, f, "cn_status")
    if not cn_status:
        issues.append("CN status not stated (formal/diagnostic/upper limit)")

    diagnostic = _get(payload, f, "diagnostic_only", False)
    if diagnostic and cn_status and "diagnostic" not in str(cn_status).lower():
        issues.append(
            "diagnostic_only=True but cn_status does not say diagnostic")

    systematics = _get(payload, f, "cn_systematics", {}) or {}
    if not isinstance(systematics, dict):
        issues.append("cn_systematics is not an object")
    else:
        budget = systematics.get("relative_1sigma")
        if budget is None:
            issues.append("CN systematic budget missing (relative_1sigma)")
        else:
            try:
                value = float(budget)
                if not 0.0 <= value < 1.0:
                    issues.append(
                        f"CN systematic budget outside 0-1: {budget!r}")
            except (TypeError, ValueError):
                issues.append(
                    f"CN systematic budget is not numeric: {budget!r}")

    n_idp = _get(payload, f, "n_independent")
    if n_idp is None:
        issues.append("N_idp missing")
    n_variables = _get(payload, f, "n_variables")
    if n_variables is None:
        issues.append("number of fit variables missing")

    for index, path in enumerate(paths, start=1):
        if not isinstance(path, dict):
            issues.append(f"path {index} is not an object")
            continue
        label = path.get("label") or path.get("shell") or f"path {index}"
        for key in _PATH_FIELDS:
            if path.get(key) is None:
                issues.append(f"{label}: {key} missing")

    return issues


def delta_e0_warnings(payload: dict) -> list[str]:
    """返回非阻断的 Delta E0 提示，不参与 `check_fit_payload` 的失败判定。

    ``|Delta E0| > 5 eV`` 时提示：检查 FEFF 路径能量零点、实验 E0
    约定和模型完整性，不要直接当成吸收边位移。项目硬门禁仍是
    ``|Delta E0| < 15 eV``。
    """
    if not isinstance(payload, dict):
        return []
    f = payload.get("fit")
    if not isinstance(f, dict):
        f = payload
    value = f.get("delta_e0")
    if value is None:
        value = payload.get("delta_e0")
    if value is None:
        return []
    try:
        value = float(value)
    except (TypeError, ValueError):
        return []
    if abs(value) <= DELTA_E0_COMFORT_EV:
        return []
    return [
        f"large |Delta E0| = {abs(value):.2f} eV (> "
        f"{DELTA_E0_COMFORT_EV:.0f} eV): check the FEFF path energy zero, "
        "experimental E0 convention, and model completeness; do not report "
        "it as an absorption-edge shift without independent XANES evidence."
    ]


def _print_report(path: Path, issues: list[str],
                  advisories: list[str] | None = None) -> None:
    if issues:
        print(f"[WARN] {path}")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print(f"[OK]   {path}")
    for advisory in advisories or []:
        print(f"  [NOTE] {advisory}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="检查 XAFS fit JSON 的 CN 披露是否达到项目标准")
    parser.add_argument("fit_json", nargs="+", help="一个或多个 fit JSON")
    parser.add_argument(
        "--warn-only", action="store_true",
        help="只打印问题，不因问题返回非零退出码")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    failed = False
    for item in args.fit_json:
        path = Path(item)
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            print(f"[ERROR] {path}: {exc}")
            failed = True
            continue
        issues = check_fit_payload(payload)
        advisories = delta_e0_warnings(payload)
        _print_report(path, issues, advisories)
        failed = failed or bool(issues)
    if failed and not args.warn_only:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
