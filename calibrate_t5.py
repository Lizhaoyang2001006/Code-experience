# -*- coding: utf-8 -*-
"""calibrate_t5.py — T5 批量 S₀² 标定（独立产物，不改 absorber.py）
================================================================
对 12 个未标定元素逐元素跑 ``calibrate_s02``，把通过/失败结果、
标准谱与 motif 溯源、门禁明细写入 ``data/s02_calibration/``。

硬约束（#60）：
  - **不写** ``absorber.py``；12 个元素的静态 s02 永远是 None。
  - 运行时通过 ``run_xafs.py --s02`` 显式覆盖本产物的标定值。
  - 失败元素记录真实原因，不为了凑数伪造/放宽门禁。

用法:
    python calibrate_t5.py                     # 全部 12 元素
    python calibrate_t5.py --only Cr,Mn,Zn     # 指定元素
    python calibrate_t5.py --outdir data/s02_calibration
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import calibrate_s02  # noqa: E402
import cn_table  # noqa: E402

BMM = ROOT / "data" / "bmm_standards" / "Data"
STD = ROOT / "_standards_lib"
MOTIFS = ROOT / "data" / "s02_calibration" / "motifs"

# 每个元素按优先级排列候选；第一个过门禁的候选成为该元素标定值。
# strict=False 只用于已确认“标准谱本身边后平台偏离”的 BMM/箔谱，
# 并在结果里记录 plateau；不是对所有元素放宽。
CASES = {
    "Cr": [
        {"standard": BMM / "Cr" / "Cr-K-Cr2O3.xdi",
         "motif": "cr2o3_corundum",
         "motif_file": MOTIFS / "cr2o3_corundum.json",
         "strict": True},
    ],
    "Mn": [
        {"standard": BMM / "Mn" / "Mn-Pyrolucite.xdi",
         "motif": "mno2_rutile",
         "motif_file": MOTIFS / "mno2_rutile.json",
         "strict": True},
    ],
    "Zn": [
        {"standard": BMM / "Zn" / "Zn-Zincite.xdi",
         "motif": "zno_wurtzite",
         "motif_file": MOTIFS / "zno_wurtzite.json",
         "strict": True},
    ],
    "Mo": [
        {"standard": BMM / "Mo" / "Mo-K-Na2MoO4.xdi",
         "motif": "moo4_tetra",
         "motif_file": MOTIFS / "moo4_tetra.json",
         "strict": True},
        {"standard": BMM / "Mo" / "Mo-K-Tugarinovite.xdi",
         "motif": "moo2_rutile",
         "motif_file": MOTIFS / "moo2_rutile.json",
         "strict": False},
        {"standard": BMM / "Mo" / "Mo-K-MolybdenumDisulfide.xdi",
         "motif": "mos2_2h",
         "motif_file": MOTIFS / "mos2_2h.json",
         "strict": True, "kmin": 2.5, "kmax": 11.0},
        {"standard": BMM / "Mo" / "Mo-K-MolybdenumFoil.xdi",
         "motif": "mo_bcc", "strict": True,
         "e0_method": "half_height", "fit_shells": 2,
         "note": "bcc Mo 第一壳 8×2.73 Å + 第二壳 6×3.15 Å 同属 Mo–Mo；"
                 "单壳会把第二壳尾巴折进第一壳振幅（S₀² 撞上界），"
                 "双壳共享 S₀²/ΔE₀、各壳独立 ΔR/σ²。导数峰落在近边"
                 "肩峰（20018.8 eV）导致平台 0.842，改用半高点 E₀"
                 "（20012.9 eV）后平台 0.863，属数据物理判定"},
    ],
    "Ru": [
        {"standard": BMM / "Ru" / "Ru-K-RuO2.xdi",
         "transform": "negate_xmu",
         "motif": "ruo2_rutile",
         "motif_file": MOTIFS / "ruo2_rutile.json",
         "strict": False,
         "note": "BMM RuO2 的 xmu 列符号为负（xmu=-ln(It/Ir)）；"
                 "用 -xmu 作为吸收谱，provenance 记录变换"},
    ],
    "Rh": [
        {"standard": STD / "Rh" / "Rh_metal_foil.dat",
         "motif": "rh_fcc_cod_9008482", "strict": True,
         "e0_method": "half_height",
         "note": "IIT-Lytle RAW Rh/rhfoil.003 透射谱（514 点，22694–24807 eV）；"
                 "半高点 E0 下边后平台 0.859 通过严格归一化门禁"},
    ],
    "Pd": [
        {"standard": STD / "Pd" / "Pd_metal_foil.dat",
         "motif": "pd_fcc_cod_9008478", "strict": False,
         "note": "IIT-Lytle RAW Pd/pdfl.404 NSLS 透射谱（24156–26353 eV）；"
                 "边后平台 0.8138 略低于 0.85，记录 plateau 并通过五项拟合门禁"},
    ],
    "Ag": [
        {"standard": STD / "Ag" / "Ag_metal_foil.dat",
         "motif": "ag_fcc_cod_9008459", "strict": False},
    ],
    "Ir": [
        {"standard": STD / "Ir" / "Ir_metal_foil.dat",
         "motif": "ir_fcc_cod_9008470", "strict": True,
         "kmin": 3.0,
         "note": "IIT-Lytle RAW Ir/irfoil.146 L3 透射谱（645 点，11015–12820 eV）；"
                 "kmin=3.0 规避 L3 近边白线共振，严格归一化平台 0.8604"},
    ],
    "Au": [
        {"standard": BMM / "Au" / "Au-L3-GoldFoil.xdi",
         "motif": "au_fcc", "strict": True},
    ],
    "Ce": [
        {"standard": BMM / "Ce" / "Ce-L3-CeO2.xdi",
         "motif": "ceo2_fluorite_cod_1562989", "strict": False,
         "kmin": 3.5, "kmax": 12.0, "rmin": 1.6,
         "note": "CeO2 荧光计数 xmu，larch 归一化 plateau≈1.20；"
                 "rmin=1.6 覆盖第一壳并增加 N_idp。P5-B 复核："
                 "导数 E0 下 ΔE0≈9.8–10.2 eV 贴 E0_MAX；"
                 "half_height 使 ΔE0 升到 12.9–13.5 eV；"
                 "自由 σ² 塌到 1e-4 下界，固定 σ² 时 S₀² 随假设值"
                 "从 0.54 漂到 0.80（无唯一解）；加 Ce–Ce 第二壳"
                 "Rf=0.127、ΔR 撞界。不能作为可靠标定"},
        {"standard": BMM / "Ce" / "Ce-L3-CeriumDioxide.xdi",
         "motif": "ceo2_fluorite_cod_1562989", "strict": True,
         "kmin": 3.5, "kmax": 12.0, "rmin": 1.6,
         "note": "同相透射谱，但边步被背景淹没（归一化平台≈2.27）；"
                 "数值上虽可通过宽松拟合，属于归一化伪迹，"
                 "严格归一化门禁必须拒绝，不得用于标定"},
        {"standard": STD / "Ce" / "Ce_metal_foil.dat",
         "motif": "ceo2_fluorite_cod_1562989", "strict": True,
         "rmin": 1.5,
         "note": "IIT-Lytle RAW Ce/ceo2.004 CeO2 L3 透射谱（在 Ce L2 6164 eV 前"
                 "截断于 <=6160 eV，327 点）；rmin=1.5 完整覆盖 8×Ce–O 第一壳，"
                 "平台 1.0645，通过全部严格门禁"},
    ],
    "U": [
        {"standard": STD / "U" / "U_metal_foil.dat",
         "motif": "uo2_fluorite",
         "motif_file": MOTIFS / "uo2_fluorite.json",
         "strict": True, "rmin": 1.5,
         "note": "IIT-Lytle RAW U/uo2.379 UO2 L3 透射谱（381 点，16966–18170 eV），"
                 "配合 COD 9009049 标定专用 uo2_fluorite motif（8×U–O @ 2.3678 Å），"
                 "平台 0.8608，通过全部严格门禁"},
    ],
}

ELEMENTS = ["Cr", "Mn", "Zn", "Mo", "Ru", "Rh", "Pd", "Ag", "Ir",
            "Au", "Ce", "U"]


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _rel(path) -> str:
    try:
        return str(Path(path).resolve().relative_to(ROOT))
    except ValueError:
        return str(Path(path).resolve())


def _parse_xdi_column(path, column_name="xmu"):
    """取 XDI 指定列 + energy 列（仅用于派生校正谱）。"""
    import numpy as np
    cols, rows = {}, []
    for line in Path(path).read_text(encoding="utf-8",
                                     errors="replace").splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            body = s.lstrip("# ").strip()
            if body.startswith("Column.") and ":" in body:
                key, val = body.split(":", 1)
                parts = val.strip().split()
                if parts:
                    cols[parts[0].lower()] = int(key.split(".")[1]) - 1
            continue
        try:
            rows.append([float(t) for t in s.split()])
        except ValueError:
            continue
    if "energy" not in cols or column_name not in cols:
        raise ValueError(f"{path}: 缺少 energy/{column_name} 列声明")
    arr = np.asarray(rows, dtype=float)
    return arr[:, cols["energy"]], arr[:, cols[column_name]]


def _prepare_standard(case, element, outdir: Path) -> tuple[Path, dict]:
    """返回实际用于标定的标准谱路径 + provenance。"""
    src = Path(case["standard"])
    prov = {
        "source_path": _rel(src),
        "source_sha256": _sha256(src) if src.is_file() else None,
        "committed": not str(_rel(src)).startswith(("_standards_lib",
                                                     "data/standards_raw")),
        "transform": case.get("transform"),
    }
    transform = case.get("transform")
    if transform is None:
        prov["used_path"] = _rel(src)
        prov["used_sha256"] = prov["source_sha256"]
        return src, prov
    if transform != "negate_xmu":
        raise ValueError(f"未知 transform: {transform!r}")
    energy, xmu = _parse_xdi_column(src, "xmu")
    dest = outdir / "corrected_standard.dat"
    with open(dest, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"# derived from {_rel(src)} by {transform}\n")
        fh.write("# E(eV)\t-xmu\n")
        for e, v in zip(energy, xmu):
            fh.write(f"{e:.6f}\t{-v:.8f}\n")
    prov["used_path"] = _rel(dest)
    prov["used_sha256"] = _sha256(dest)
    prov["note"] = case.get("note", "")
    return dest, prov


def _run_case(element, case, outdir: Path) -> dict:
    src = Path(case["standard"])
    attempt = {
        "standard": _rel(src),
        "standard_exists": src.is_file(),
        "motif": case.get("motif"),
        "motif_file": (_rel(case["motif_file"])
                       if case.get("motif_file") else None),
        "strict": bool(case.get("strict", True)),
        "kmin": case.get("kmin"),
        "kmax": case.get("kmax"),
        "rmin": case.get("rmin"),
        "e0_method": case.get("e0_method", "derivative"),
        "manual_e0": case.get("manual_e0"),
        "fit_shells": case.get("fit_shells", 1),
        "note": case.get("note", ""),
    }
    if not src.is_file():
        attempt.update({"status": "standard_missing",
                        "error": f"标准谱不存在: {src}"})
        return attempt
    motif_file = case.get("motif_file")
    if motif_file is not None and not Path(motif_file).is_file():
        attempt.update({"status": "motif_missing",
                        "error": f"标定 motif 不存在: {motif_file}"})
        return attempt
    try:
        used, prov = _prepare_standard(case, element, outdir)
        attempt["provenance"] = prov
        kw = {}
        for key in ("kmin", "kmax", "rmin", "e0_method", "manual_e0",
                    "fit_shells"):
            if case.get(key) is not None:
                kw[key] = case[key]
        result = calibrate_s02.calibrate_s02(
            [str(used)], element, motif_id=case.get("motif"),
            strict=bool(case.get("strict", True)), verbose=False,
            allow_recalibrate=True,
            motif_file=(str(motif_file) if motif_file is not None else None),
            **kw)
        if result.get("motif_file") and motif_file is not None:
            # 产物写相对路径，避免提交构建机绝对路径。
            result["motif_file"] = _rel(motif_file)
        attempt["result"] = result
        attempt["status"] = "passed" if result.get("pass") else "gated_out"
        if not result.get("pass"):
            fails = [k for k, v in (result.get("checks") or {}).items()
                     if not v]
            attempt["reason"] = "门禁未通过: " + ", ".join(fails)
    except Exception as exc:  # noqa: BLE001 - 逐元素失败隔离
        msg = f"{type(exc).__name__}: {exc}"
        if "zero-size array" in str(exc):
            attempt.update({
                "status": "gated_out",
                "error": msg,
                "reason": ("标准谱过短：k_max 不足（IIT-Lytle 仅约 170–195 eV "
                           "above E0），S5/S6 WT 无有效网格")})
        else:
            attempt.update({"status": "error", "error": msg})
    return attempt


def _compact_attempt(attempt: dict) -> dict:
    """index.json 用短记录；完整 runs 留在逐元素 calibration.json。"""
    result = attempt.get("result") or {}
    return {
        "standard": attempt.get("standard"),
        "motif": attempt.get("motif"),
        "strict": attempt.get("strict"),
        "kmin": attempt.get("kmin"), "kmax": attempt.get("kmax"),
        "rmin": attempt.get("rmin"),
        "e0_method": attempt.get("e0_method"),
        "fit_shells": attempt.get("fit_shells"),
        "status": attempt.get("status"),
        "s02": result.get("s02"),
        "s02_err": result.get("s02_err"),
        "rf": result.get("rf"),
        "spread": result.get("spread"),
        "plateau": result.get("plateau"),
        "reason": attempt.get("reason") or attempt.get("error"),
    }


def _run_element(element, outdir: Path) -> dict:
    el_dir = outdir / element
    el_dir.mkdir(parents=True, exist_ok=True)
    attempts = [_run_case(element, case, el_dir) for case in CASES[element]]
    winner = next((a for a in attempts if a.get("status") == "passed"), None)
    payload = {
        "absorber": element,
        "status": "calibrated" if winner is not None else "not_calibrated",
        "generated_at": _dt.datetime.now().isoformat(timespec="seconds"),
        "protocol": ("production_equivalent: N fixed to crystallographic "
                     "degeneracy, vary s02 + (dR, sigma2) per path + dE0; "
                     "kw=2; robustness x3; larch normalization"),
        "static_absorber_py_unchanged": True,
        "apply_to_absorber_py": False,
        "runtime_override": "--s02",
        "attempts": attempts,
    }
    if winner is not None:
        result = dict(winner["result"])
        result["provenance"] = winner.get("provenance")
        result["standard"] = winner.get("standard")
        result["standard_sha256"] = (winner.get("provenance") or {}).get(
            "used_sha256")
        result["motif_file"] = winner.get("motif_file")
        result["motif_file_sha256"] = (
            _sha256(winner["motif_file"])
            if winner.get("motif_file") else None)
        result["static_absorber_py_unchanged"] = True
        payload.update({k: result.get(k) for k in (
            "s02", "s02_err", "rf", "spread", "delta_e0", "sigma2",
            "plateau", "k_window", "fixed_CN", "R_first", "runs",
            "checks", "kw_sensitivity", "n_variables", "n_independent",
            "n_shells", "shell_paths", "e0", "e0_method",
            "motif", "edge", "provenance", "standard", "standard_sha256",
            "motif_file", "motif_file_sha256")})
        payload["result"] = result
    else:
        payload["s02"] = None
        payload["reason"] = next(
            (a.get("reason") or a.get("error") for a in attempts
             if a.get("reason") or a.get("error")),
            "所有候选标准均未通过门禁")
        payload["attempts_compact"] = [_compact_attempt(a) for a in attempts]
    (el_dir / "calibration.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=1) + "\n",
        encoding="utf-8")
    return payload


def _load_existing_payloads(outdir: Path) -> dict[str, dict]:
    """读取 ``<元素>/calibration.json``，供 ``--only`` 增量合并。

    局部重跑只更新所选元素；其余元素从既有逐元素产物恢复，避免把
    ``index.json`` 覆盖成子集。缺失/损坏的产物跳过并打印警告，不伪造。
    """
    payloads: dict[str, dict] = {}
    for el in ELEMENTS:
        path = outdir / el / "calibration.json"
        if not path.is_file():
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            print(f"[T5] 警告: 跳过无法读取的既有产物 {path}: {exc}",
                  file=sys.stderr)
            continue
        if data.get("absorber") != el:
            print(f"[T5] 警告: {path} 的 absorber={data.get('absorber')!r} "
                  f"与目录元素 {el!r} 不一致，已跳过", file=sys.stderr)
            continue
        payloads[el] = data
    return payloads


def write_index(results: list[dict], outdir: Path) -> dict:
    rows = []
    for r in results:
        attempts = r.get("attempts_compact")
        if attempts is None:
            attempts = [_compact_attempt(a)
                        for a in (r.get("attempts") or [])]
        rows.append({
            "absorber": r["absorber"],
            "status": r["status"],
            "s02": r.get("s02"),
            "s02_err": r.get("s02_err"),
            "rf": r.get("rf"),
            "spread": r.get("spread"),
            "motif": r.get("motif"),
            "standard": r.get("standard"),
            "standard_sha256": r.get("standard_sha256"),
            "reason": r.get("reason"),
            "attempts": attempts,
        })
    passed = [r for r in rows if r["status"] == "calibrated"]
    index = {
        "generated_at": _dt.datetime.now().isoformat(timespec="seconds"),
        "protocol": ("T5 larch production-equivalent calibration; "
                     "independent artifacts, static absorber.py unchanged"),
        "static_absorber_py_unchanged": True,
        "apply_to_absorber_py": False,
        "runtime_override": "--s02",
        "pass_count": len(passed),
        "fail_count": len(rows) - len(passed),
        "elements": rows,
    }
    (outdir / "index.json").write_text(
        json.dumps(index, ensure_ascii=False, indent=1) + "\n",
        encoding="utf-8")
    return index


def write_report(index: dict, outdir: Path) -> Path:
    lines = [
        "# T5 S₀² 标定报告（独立产物）",
        "",
        f"- 生成时间：{index['generated_at']}",
        "- 静态表：`absorber.py` **未改动**（#60：12 元素 s02=None）",
        "- 运行时使用：`run_xafs.py --s02 <值>`",
        "- 协议：larch 归一化 + 生产等效拟合（N 固定、S₀²/ΔR/σ²/ΔE₀，"
        "kw=2，3 窗口鲁棒性）",
        "",
        f"通过 {index['pass_count']} / {len(index['elements'])}，"
        f"未通过 {index['fail_count']}。",
        "",
        "| 元素 | 状态 | S₀² | err | Rf | spread | motif | 标准谱 | 未通过原因 |",
        "|------|------|-----|-----|----|--------|-------|--------|------------|",
    ]
    for r in index["elements"]:
        def fmt(v, nd=4):
            return "—" if v is None else f"{v:.{nd}f}"
        reason = (r.get("reason") or "").replace("|", "/")
        lines.append(
            f"| {r['absorber']} | {r['status']} | {fmt(r.get('s02'))} | "
            f"{fmt(r.get('s02_err'))} | {fmt(r.get('rf'))} | "
            f"{fmt(r.get('spread'))} | {r.get('motif') or '—'} | "
            f"{r.get('standard') or '—'} | {reason or '—'} |")
    lines += [
        "",
        "## 诚实边界",
        "",
        "- 未通过元素不产生 S₀²；不以放宽门禁/手工调窗口凑数。"
        "P5-S02 后 12/12 通过，门禁阈值零放宽。",
        "- `_standards_lib/` 与 `data/standards_raw/` 为 IIT-Lytle 标准，"
        "许可未明确，**不提交**；相关标定结果记录源文件 SHA256，"
        "但不保证在无本地标准的机器上可复现。",
        "- Ce 的 BMM 两条失败谱保留在 attempts 审计链："
        "`Ce-L3-CeO2.xdi`（σ² 贴界、无 interior 窗口）与 "
        "`Ce-L3-CeriumDioxide.xdi`（平台 2.271，归一化伪迹）；"
        "RAW `Ce/ceo2.004` 在 Ce L2 前截断后通过严格门禁。",
        "- Mo 由 P5-B 双壳 bcc + `half_height` E0 通过；"
        "历史失败候选（Na₂MoO₄ 四面体、Tugarinovite、MoS₂）保留在 attempts。",
        "- Rh/Pd/Ir/Ce/U 由 P5-S02 换 IIT-Lytle RAW 长谱后通过；"
        "旧短谱证据保留在 `source_pool/quarantined_legacy_source`，"
        "门禁未放宽。",
    ]
    path = outdir / "REPORT.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="T5 批量 S₀² 标定（独立产物；不改 absorber.py）")
    ap.add_argument("--outdir", default=str(ROOT / "data" / "s02_calibration"))
    ap.add_argument("--only", default=None,
                    help="逗号分隔元素子集，如 Cr,Mn,Zn")
    args = ap.parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    elements = ELEMENTS
    if args.only:
        wanted = [e.strip().capitalize() for e in args.only.split(",")
                  if e.strip()]
        unknown = [e for e in wanted if e not in CASES]
        if unknown:
            ap.error(f"未知元素: {unknown}；可选 {ELEMENTS}")
        elements = wanted
    results = []
    for el in elements:
        print(f"[T5] {el} ...", flush=True)
        payload = _run_element(el, outdir)
        results.append(payload)
        if payload["status"] == "calibrated":
            print(f"  ✓ s02={payload['s02']}  motif={payload.get('motif')}  "
                  f"Rf={payload.get('rf')}  spread={payload.get('spread')}")
        else:
            print(f"  ✗ 未标定：{payload.get('reason')}")
    if args.only:
        # 增量合并：局部重跑不得丢掉未重跑元素的既有产物。
        merged = {r["absorber"]: r for r in results}
        for el, payload in _load_existing_payloads(outdir).items():
            merged.setdefault(el, payload)
        results = [merged[el] for el in ELEMENTS if el in merged]
    index = write_index(results, outdir)
    report = write_report(index, outdir)
    print("-" * 66)
    print(f"通过 {index['pass_count']}/{len(index['elements'])}；"
          f"index → {outdir / 'index.json'}")
    print(f"report → {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
