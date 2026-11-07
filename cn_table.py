# -*- coding: utf-8 -*-
"""cn_table.py — 正式 CN 表的读写与拒绝行（P0a）

``cn_table.csv`` 是报告唯一的正式 CN 来源；诊断拟合表继续由
``xafs_report.export_fit_table`` 输出，两者不得混用。

P0a 只定义表结构、拒绝行和读写往返；正式 CN 行从 P3 才开始产生。
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

# 与 REVIEW_AND_REVISED_PLAN.md §4.4 保持一致；顺序即 CSV 列顺序。
CN_TABLE_COLUMNS = [
    "phase_id", "model_level", "source_type", "source_ref",
    "source_cif_sha256", "absorber_site", "shell_group", "shell_label",
    "scatterer", "r_cif", "n0", "n0_members", "n_fit", "n_status",
    "n_ratio", "n_err_formal", "n_err_systematic", "delta_r", "sigma2",
    "delta_e0", "rf", "s02", "s02_source", "s02_calibrated",
    "kmin", "kmax", "kweight", "rmin", "rmax", "n_idp", "n_variables",
    "cn_label", "gate_status", "reject_reason", "flags", "assumptions",
    "data_sha256", "software_version",
]


def _cell(value) -> str:
    """把单元格值稳定地写成字符串；None/NaN 一律空串。"""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return repr(value)
    return str(value)


def refusal_row(gate_status: str, reject_reason: str, **meta) -> dict:
    """构造一行“拒绝出 CN”的记录；CN 数值列全部为空。

    拒绝行不是 CN 行：``n_fit`` 必须为空，``gate_status`` 记录拒绝阶段，
    ``reject_reason`` 记录可复核原因。
    """
    row = {column: "" for column in CN_TABLE_COLUMNS}
    row["gate_status"] = str(gate_status or "no_cif")
    row["reject_reason"] = str(reject_reason or "未生成正式 CN")
    row["n_status"] = "diagnostic"
    row["cn_label"] = "none"
    row["flags"] = "refused"
    for key, value in meta.items():
        if key in row and value is not None:
            row[key] = value
    return row


def write_cn_table(path: str | Path, rows: list[dict]) -> None:
    """写正式 CN 表；UTF-8 BOM 便于 Excel 直接打开。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=CN_TABLE_COLUMNS,
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: _cell(row.get(column, ""))
                             for column in CN_TABLE_COLUMNS})


def read_cn_table(path: str | Path) -> list[dict]:
    """读回正式 CN 表；所有值保持字符串，便于逐字节复核。"""
    with Path(path).open("r", newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))
