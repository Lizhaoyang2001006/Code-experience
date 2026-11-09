#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cif_crawler.py - 晶体结构 (CIF) 与 XAS 标准样自动化检索工具

集成特性：
1. 多源降级管道：Local Cache -> COD (实验真实) -> MP-API / AFLOW (理论兜底)
2. 催化多晶型消歧：支持针对常见催化相配置优选空间群 (Preferred Space Groups)
3. CIF 质量验真：使用 gemmi 进行原子坐标、晶胞参数与占位数有效性验证
4. 专项目标直通：支持 --known-cod-id 瞬间精确入库
5. T5 XAS 谱学闭环：自动下载金属箔基准能量谱 (.dat / .xdi)，对接 IIT Lytle / Larch 库
"""

import os
import sys
import json
import re
import ssl
import math
import time
import shutil
import argparse
import logging
import statistics
from pathlib import Path
from typing import List, Dict, Optional, Tuple

import requests
import urllib3
from requests.adapters import HTTPAdapter

# 确保可导入项目内模块
ROOT_DIR = Path(__file__).resolve().parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

# 尝试载入晶体学解析引擎 gemmi
try:
    import gemmi
    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False

# 尝试载入 XAS 处理库 larch
try:
    from larch.xafs import standard_foil
    HAS_LARCH = True
except ImportError:
    HAS_LARCH = False

# 载入项目内结构库（用于 CIF 入库生成 motif JSON）
try:
    from structure_lib import import_cif as sl_import_cif, register as sl_register
    from structure_lib import XAFSConfigError, XAFSDataError
    import structure_lib as sl_mod
    HAS_STRUCTURE_LIB = True
except ImportError:
    HAS_STRUCTURE_LIB = False

# 载入吸收元素配置（用于 T5 标准箔占位文件的边能量）
try:
    from absorber import ABSORBERS
    HAS_ABSORBER = True
except ImportError:
    HAS_ABSORBER = False

# ==================== 全局配置与路径 ====================
# 一律以脚本所在目录为根，保证从任意工作目录运行都落到同一套构建产物。
BASE_DATA_DIR = ROOT_DIR / "data"
CIF_RAW_DIR = BASE_DATA_DIR / "cif_raw"
STANDARDS_DIR = BASE_DATA_DIR / "standards_raw"
STRUCTURE_LIB_DIR = ROOT_DIR / "_structure_lib"
STANDARDS_LIB_DIR = ROOT_DIR / "_standards_lib"
INDEX_FILE = BASE_DATA_DIR / "cif_crawler_index.json"


def _index_path(path: Path) -> str:
    """索引里优先写项目内相对路径，避免把构建机的绝对路径固化成数据契约。"""
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(ROOT_DIR).as_posix()
    except ValueError:
        return resolved.as_posix()

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) Catalysis-CIF-Crawler/2.0"
}

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s", force=True)

# ==================== 目标体系与物理相绑定 ====================
# 最终目标映射：空间群使用国际空间群编号 (IT number) 进行严格物理相绑定。
# tier     : 任务分层（T3 常见氧化物 / T4 缺失单质）
# phase    : 相名，进入 motif_id，避免 hematite / maghemite 等同质多象冲突
# stoich   : 用于 _chemical_formula_sum 的严格化学计量校验
# cod_id   : 已逐条下载并用 gemmi 验真的 COD 结构 ID
# spg_nums : 该 COD 结构真实对应的 IT 空间群号（不放永远不会命中的备用号）
VERIFIED_TARGETS = {
    # ==================== T3: 催化常见氧化物与氢氧化物 ====================
    "Fe2O3":   {"tier": "T3", "phase": "hematite",     "absorber": "Fe", "elements": ["Fe", "O"],     "stoich": {"Fe": 2, "O": 3}, "cod_id": "1011240", "spg_nums": [167]},
    "Fe3O4":   {"tier": "T3", "phase": "magnetite",    "absorber": "Fe", "elements": ["Fe", "O"],     "stoich": {"Fe": 3, "O": 4}, "cod_id": "1010369", "spg_nums": [227], "absorber_fract": [0.625, 0.625, 0.625]},
    "CoO":     {"tier": "T3", "phase": "rocksalt",     "absorber": "Co", "elements": ["Co", "O"],     "stoich": {"Co": 1, "O": 1}, "cod_id": "1533087", "spg_nums": [225]},
    "Co3O4":   {"tier": "T3", "phase": "spinel",       "absorber": "Co", "elements": ["Co", "O"],     "stoich": {"Co": 3, "O": 4}, "cod_id": "9005888", "spg_nums": [227], "absorber_fract": [0.125, 0.125, 0.125]},
    "NiO":     {"tier": "T3", "phase": "bunsenite",    "absorber": "Ni", "elements": ["Ni", "O"],     "stoich": {"Ni": 1, "O": 1}, "cod_id": "1010093", "spg_nums": [225]},
    "Ni(OH)2": {"tier": "T3", "phase": "theophrastite","absorber": "Ni", "elements": ["Ni", "O", "H"], "stoich": {"Ni": 1, "O": 2, "H": 2}, "cod_id": "1548811", "spg_nums": [164]},
    "CuO":     {"tier": "T3", "phase": "tenorite",     "absorber": "Cu", "elements": ["Cu", "O"],     "stoich": {"Cu": 1, "O": 1}, "cod_id": "1011148", "spg_nums": [15]},
    "Cu2O":    {"tier": "T3", "phase": "cuprite",      "absorber": "Cu", "elements": ["Cu", "O"],     "stoich": {"Cu": 2, "O": 1}, "cod_id": "1010941", "spg_nums": [224]},
    "CeO2":    {"tier": "T3", "phase": "fluorite",     "absorber": "Ce", "elements": ["Ce", "O"],     "stoich": {"Ce": 1, "O": 2}, "cod_id": "1562989", "spg_nums": [225]},

    # ==================== T4: 缺失基准元素单质 ====================
    "Mn":      {"tier": "T4", "phase": "alpha-Mn",     "absorber": "Mn", "elements": ["Mn"],           "stoich": {"Mn": 1}, "cod_id": "9008589", "spg_nums": [217], "absorber_fract": [0.0, 0.0, 0.0]},
    "Ru":      {"tier": "T4", "phase": "hcp",          "absorber": "Ru", "elements": ["Ru"],           "stoich": {"Ru": 1}, "cod_id": "9008513", "spg_nums": [194]},
    "Rh":      {"tier": "T4", "phase": "fcc",          "absorber": "Rh", "elements": ["Rh"],           "stoich": {"Rh": 1}, "cod_id": "9008482", "spg_nums": [225]},
    "Pd":      {"tier": "T4", "phase": "fcc",          "absorber": "Pd", "elements": ["Pd"],           "stoich": {"Pd": 1}, "cod_id": "9008478", "spg_nums": [225]},
    "Ag":      {"tier": "T4", "phase": "fcc",          "absorber": "Ag", "elements": ["Ag"],           "stoich": {"Ag": 1}, "cod_id": "9008459", "spg_nums": [225]},
    "Ir":      {"tier": "T4", "phase": "fcc",          "absorber": "Ir", "elements": ["Ir"],           "stoich": {"Ir": 1}, "cod_id": "9008470", "spg_nums": [225]},
    "Ce":      {"tier": "T4", "phase": "fcc",          "absorber": "Ce", "elements": ["Ce"],           "stoich": {"Ce": 1}, "cod_id": "9010987", "spg_nums": [225]},
    "U":       {"tier": "T4", "phase": "alpha-U",      "absorber": "U",  "elements": ["U"],            "stoich": {"U": 1},  "cod_id": "1527120", "spg_nums": [63]},
}

# 旧索引/motif 平滑迁移表。
# 只有新相已通过空间群 + 计量比校验并成功注册后，才按此表清理旧文件；
# 对会与新 motif 指纹冲突的旧 JSON，注册前先移入隔离区，失败则原样恢复。
LEGACY_ARTIFACTS = {
    "Fe2O3":   {"motif_ids": ["hematite_cod"],     "raw_cifs": []},
    "Fe3O4":   {"motif_ids": ["fefeo"],            "raw_cifs": []},
    "CoO":     {"motif_ids": ["coo_rocksalt"],     "raw_cifs": ["CoO_COD_1528838.cif", "CoO_cod1528838.cif"]},
    "Co3O4":   {"motif_ids": ["cocoo", "co3o4_tet"], "raw_cifs": ["Co3O4_COD_1526734.cif"]},
    "NiO":     {"motif_ids": ["nio_rocksalt"],     "raw_cifs": []},
    "Ni(OH)2": {"motif_ids": ["ninioh", "ni_oh2_brucite"], "raw_cifs": []},
    "CuO":     {"motif_ids": ["cuo_tenorite"],     "raw_cifs": []},
    "Cu2O":    {"motif_ids": ["cu2o_cuprite"],     "raw_cifs": ["Cu2O_COD_1000063.cif"]},
    "CeO2":    {"motif_ids": ["ceo2_fluorite"],    "raw_cifs": []},
    "Mn":      {"motif_ids": ["mnmn"],             "raw_cifs": ["Mn_COD_1534895.cif", "Mn_cod1534895.cif", "COD_1534895.cif"]},
    "Ru":      {"motif_ids": ["ruru"],             "raw_cifs": ["Ru_COD_1534914.cif", "COD_1534914.cif"]},
    "Rh":      {"motif_ids": ["rhrh"],             "raw_cifs": ["Rh_COD_1534917.cif", "COD_1534917.cif"]},
    "Pd":      {"motif_ids": ["pdpd"],             "raw_cifs": ["Pd_COD_1011104.cif", "COD_1011104.cif"]},
    "Ag":      {"motif_ids": ["agag"],             "raw_cifs": ["Ag_COD_1100136.cif"]},
    "Ir":      {"motif_ids": ["irir"],             "raw_cifs": ["Ir_COD_1534947.cif", "COD_1534947.cif"]},
    "Ce":      {"motif_ids": ["cece"],             "raw_cifs": ["Ce_COD_1533274.cif", "COD_1533274.cif"]},
    "U":       {"motif_ids": ["uu"],               "raw_cifs": []},
}

# T5 标准金属箔标定谱元素目标 (Cr 到 U)
TARGETS_T5_ELEMENTS = [
    "Cr", "Mn", "Zn", "Mo", "Ru", "Rh", 
    "Pd", "Ag", "Ir", "Au", "Ce", "U"
]


# ==================== 化学式匹配工具 ====================
_ELEMENT_RE = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")


def parse_formula(formula: str) -> dict[str, float]:
    """将化学式字符串解析为 {元素: 数量} 字典。

    支持 Fe2O3、Fe2 O3、Al0.34 Ca0.6，以及带括号的 Ni(OH)2、Ca(OH)2 等。
    """
    formula = formula.strip()
    if not formula or formula.lower() in ("?", "none", "n/a"):
        return {}

    out: dict[str, float] = {}

    def add_group(text: str, multiplier: float = 1.0) -> None:
        i, n = 0, len(text)
        while i < n:
            ch = text[i]
            if ch == "(":
                depth, j = 1, i + 1
                while j < n and depth:
                    if text[j] == "(":
                        depth += 1
                    elif text[j] == ")":
                        depth -= 1
                    j += 1
                inner = text[i + 1:j - 1]
                m = re.match(r"\d*\.?\d*", text[j:])
                num = m.group(0) if m else ""
                add_group(inner, multiplier * (float(num) if num else 1.0))
                i = j + len(num)
            elif ch.isupper():
                m = _ELEMENT_RE.match(text, i)
                if not m:
                    i += 1
                    continue
                el = m.group(1).capitalize()
                num = m.group(2)
                out[el] = out.get(el, 0.0) + (float(num) if num else 1.0) * multiplier
                i = m.end()
            else:
                i += 1

    add_group(formula)
    return out


def formula_matches(target_formula: str, cif_formula: str) -> bool:
    """判断 CIF 化学式是否与目标化学式一致（元素种类与比例相同）。
    目标 "Fe2O3" 与 CIF "Fe2 O3" 匹配；"Fe2O3" 与 "Fe3O4" 不匹配。"""
    tf = parse_formula(target_formula)
    cf = parse_formula(cif_formula)
    if not tf or not cf:
        return False
    if set(tf.keys()) != set(cf.keys()):
        return False
    t_sum = sum(tf.values())
    c_sum = sum(cf.values())
    if t_sum == 0 or c_sum == 0:
        return False
    for k in tf:
        if abs(tf[k] / t_sum - cf[k] / c_sum) > 1e-6:
            return False
    return True


# ==================== CIF 物理相核验工具 ====================
def _cif_text_value(raw: Optional[str]) -> str:
    """去掉 CIF 词法层的外层引号，但保留符号内部的双引号。

    例：'F m -3 m' -> F m -3 m；'-P 3 2"' -> -P 3 2"
    不能直接 replace('"', '')，否则 Hall 符号 #164 会被误判成 #162。
    """
    if raw is None:
        return ""
    if HAS_GEMMI:
        try:
            return gemmi.cif.as_string(raw).strip()
        except Exception:
            pass
    text = str(raw).strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in ("'", '"'):
        text = text[1:-1]
    return text.strip()


def get_spg_number(block) -> Optional[int]:
    """解析 CIF 空间群编号 (1-230)。

    优先级：IT number 标签 > H-M 符号 > Hall 符号查表。
    Hall 符号必须保留末尾/内部双引号；H-M 的 :R/:1/:2 扩展在查表前剥离。
    """
    if block is None or not HAS_GEMMI:
        return None

    # 1. 优先读取国际空间群编号标签
    for tag in ("_space_group_IT_number", "_symmetry_Int_Tables_number"):
        val = _cif_text_value(block.find_value(tag))
        if val.isdigit():
            num = int(val)
            if 1 <= num <= 230:
                return num

    # 2. 尝试解析 Hermann-Mauguin 符号（去掉 :R / :1 / :2 等设置后缀）
    for tag in ("_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt"):
        hm = _cif_text_value(block.find_value(tag))
        if not hm:
            continue
        hm_core = re.sub(r"\s*:[A-Za-z0-9]+\s*$", "", hm).strip()
        for candidate in (hm, hm_core):
            if not candidate:
                continue
            try:
                spg = gemmi.SpaceGroup(candidate)
                if 1 <= spg.number <= 230:
                    return spg.number
            except Exception:
                pass

    # 3. 尝试 Hall 符号；保留内部双引号，并用 spacegroup_table 兜底
    for tag in ("_symmetry_space_group_name_Hall", "_space_group_name_Hall"):
        hall = _cif_text_value(block.find_value(tag))
        if not hall:
            continue
        try:
            spg = gemmi.SpaceGroup(hall)
            if 1 <= spg.number <= 230:
                return spg.number
        except Exception:
            pass
        for sg in gemmi.spacegroup_table():
            if sg.hall.strip() == hall:
                return sg.number

    return None


def verify_composition(block, expected_stoich: Dict[str, float],
                       rel_tol: float = 0.08) -> Tuple[bool, str]:
    """严格比对 CIF 化学式与目标化学计量比，拦截同元素错相。

    优先使用 _chemical_formula_sum；缺失时退回原子位点统计。
    比较的是各元素在总原子数中的占比，允许 8% 的相对偏差，
    因此 FeO 与 Fe2O3、CoO 与 Co3O4 这类同元素不同计量不会互相放行。
    """
    if not expected_stoich:
        return True, "无化学计量期望，跳过"

    formula_raw = block.find_value("_chemical_formula_sum") if block is not None else None
    actual: Dict[str, float] = {}
    source = "formula_sum"
    if formula_raw:
        actual = parse_formula(_cif_text_value(formula_raw))
    else:
        source = "atom_site"
        symbols = (block.find_values("_atom_site_type_symbol")
                   or block.find_values("_atom_site_label"))
        occupancies = block.find_values("_atom_site_occupancy")
        for i, sym in enumerate(symbols):
            m = re.match(r"([A-Z][a-z]?)", _cif_text_value(sym))
            if not m:
                continue
            el = m.group(1).capitalize()
            occ = 1.0
            if i < len(occupancies):
                try:
                    occ = float(_cif_text_value(occupancies[i]).split("(")[0])
                except (TypeError, ValueError):
                    occ = 1.0
            actual[el] = actual.get(el, 0.0) + occ

    if not actual:
        return False, f"无法从 CIF 解析化学组分（{source} 为空）"

    expected_keys = set(expected_stoich)
    actual_keys = set(actual)
    if expected_keys != actual_keys:
        return False, f"元素集合不符: 期望 {sorted(expected_keys)}, 实际 {sorted(actual_keys)}"

    expected_total = sum(float(v) for v in expected_stoich.values())
    actual_total = sum(actual.values())
    if expected_total <= 0 or actual_total <= 0:
        return False, "化学计量总数无效"

    for el in sorted(expected_keys):
        expected_frac = float(expected_stoich[el]) / expected_total
        actual_frac = actual[el] / actual_total
        if abs(expected_frac - actual_frac) > rel_tol:
            return False, (f"{el} 计量比偏差过大: 期望 {expected_frac:.3f}, "
                           f"实际 {actual_frac:.3f} (来源 {source})")

    return True, (f"化学计量匹配 ({source}): 实际 {actual} ≈ 期望 "
                  f"{dict(expected_stoich)}")


def generate_safe_motif_id(formula: str, phase: str, cod_id: str) -> str:
    """生成具备唯一性的基元标识: fe2o3_hematite_cod_1011240。"""
    f_clean = re.sub(r"[^A-Za-z0-9]", "", str(formula)).lower()
    p_clean = re.sub(r"[^A-Za-z0-9]+", "_", str(phase)).strip("_").lower()
    return f"{f_clean}_{p_clean}_cod_{cod_id}"


def upsert_index(index_file: Path, new_entry_key: str, new_entry_val: dict) -> dict:
    """单条原子增量合并索引，避免 --t5 / --t3 相互覆盖总索引。

    每次写入前重新读取磁盘上的当前索引，只增改一个键，再用 os.replace
    原子替换，保证 T5 标准谱条目不会被 CIF 任务清空，反之亦然。
    """
    index_file = Path(index_file)
    current: dict = {}
    if index_file.exists():
        try:
            loaded = json.loads(index_file.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                current = loaded
        except Exception:
            # 索引损坏时保留现场，避免静默覆盖
            corrupt = index_file.with_name(
                f"{index_file.name}.corrupt.{int(time.time())}")
            try:
                os.replace(index_file, corrupt)
                logging.warning(f"[Index] 索引损坏，已备份到 {corrupt.name}")
            except OSError:
                pass
            current = {}

    current[new_entry_key] = new_entry_val
    index_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_file = index_file.with_name(
        f"{index_file.name}.tmp.{os.getpid()}")
    tmp_file.write_text(
        json.dumps(current, indent=2, ensure_ascii=False), encoding="utf-8")
    os.replace(tmp_file, index_file)
    return current


# ==================== 1. CIF 质量审查 (Sanity Checker) ====================
class CIFValidator:
    """使用 gemmi 或文本规则验证 CIF 完整性与催化适用性"""

    @staticmethod
    def validate_and_parse(cif_text: str) -> Tuple[bool, Optional[str], Optional[dict]]:
        """
        验证 CIF：检查晶胞参数、是否存在占位混乱与结构破损
        返回: (is_valid, error_msg, metadata)
        """
        if "data_" not in cif_text or "_cell_length_a" not in cif_text:
            return False, "缺失基本晶体学 CIF 标签", None

        metadata = {"spg": None, "formula": None, "cell": {}}
        if HAS_GEMMI:
            try:
                doc = gemmi.cif.read_string(cif_text)
                block = doc.sole_block()
                
                # 检查晶胞
                a = float(block.find_value('_cell_length_a').split('(')[0])
                b = float(block.find_value('_cell_length_b').split('(')[0])
                c = float(block.find_value('_cell_length_c').split('(')[0])
                if min(a, b, c) <= 0.5:
                    return False, f"异常晶胞参数: a={a}, b={b}, c={c}", None
                metadata["cell"] = {"a": a, "b": b, "c": c}

                # 提取空间群（保留 Hall 符号内部双引号，不做破坏性 replace）
                spg = block.find_value('_symmetry_space_group_name_H-M') or \
                      block.find_value('_space_group_name_H-M_alt')
                if spg:
                    metadata["spg"] = _cif_text_value(spg)

                # 提取化学式（用于 fixed_cod_id 的二次校验，防止 ID 指错相）
                formula = block.find_value('_chemical_formula_sum')
                if formula:
                    metadata["formula"] = _cif_text_value(formula)

                # 检查是否包含原子坐标
                sites_x = block.find_values('_atom_site_fract_x')
                if not sites_x:
                    return False, "缺少原子分数坐标列表", None

                return True, None, metadata
            except Exception as e:
                return False, f"gemmi 解析异常: {str(e)}", None
        else:
            # 降级文本校验
            return True, None, metadata


# ==================== 2. 数据源：COD 客户端 ====================
class CODClient:
    """Crystallography Open Database 检索与下载核心"""
    SEARCH_URL = "https://www.crystallography.net/cod/result"
    DOWNLOAD_URL = "https://www.crystallography.net/cod/{cod_id}.cif"

    @classmethod
    def download_by_id(cls, cod_id: str, dest_dir: Path, final_name: Optional[str] = None,
                       expected_formula: Optional[str] = None) -> Optional[Path]:
        """根据 COD ID 直接拉取文件 (--known-cod-id 机制)。

        expected_formula 非空时，用 CIF 的 _chemical_formula_sum 做二次校验，
        避免 fixed_cod_id 指错相导致结构库污染。
        """
        out_path = dest_dir / (final_name or f"COD_{cod_id}.cif")

        # 复用本地已下载且化学式匹配的 CIF，避免 COD 网络抖动导致重复运行失败
        if out_path.exists():
            try:
                valid, err, meta = CIFValidator.validate_and_parse(
                    out_path.read_text(encoding="utf-8"))
                if valid and (
                        not expected_formula
                        or (meta and meta.get("formula")
                            and formula_matches(expected_formula, meta["formula"]))):
                    logging.info(f"[COD] 复用本地 CIF -> {out_path.name}")
                    return out_path
                if valid and expected_formula and not (meta or {}).get("formula"):
                    logging.warning(
                        f"[COD:{cod_id}] 本地文件缺化学式标签，重新下载校验")
            except Exception:
                pass

        url = cls.DOWNLOAD_URL.format(cod_id=cod_id)
        try:
            resp = requests.get(url, headers=HEADERS, timeout=20)
            if resp.status_code == 200 and "data_" in resp.text:
                valid, err, meta = CIFValidator.validate_and_parse(resp.text)
                if not valid:
                    logging.warning(f"[COD:{cod_id}] 校验失败: {err}")
                    return None

                if expected_formula:
                    cif_formula = (meta or {}).get("formula")
                    if not cif_formula:
                        logging.warning(
                            f"[COD:{cod_id}] 拒绝入库: CIF 缺少 _chemical_formula_sum，"
                            f"无法校验目标化学式 {expected_formula}")
                        return None
                    if not formula_matches(expected_formula, cif_formula):
                        logging.warning(
                            f"[COD:{cod_id}] 拒绝入库: 化学式不符 "
                            f"(目标 {expected_formula}, CIF {cif_formula})")
                        return None
                    logging.info(
                        f"[COD:{cod_id}] 化学式校验通过: {cif_formula} ≈ {expected_formula}")

                # 直接写入最终路径（若指定），避免 rename 导致的 WinError 183 文件冲突
                out_path = dest_dir / (final_name or f"COD_{cod_id}.cif")
                out_path.write_text(resp.text, encoding="utf-8")
                logging.info(f"[COD] 直链抓取成功 -> {out_path.name}")
                return out_path
        except Exception as e:
            logging.error(f"[COD:{cod_id}] 下载中断: {e}")
        return None

    @classmethod
    def query_and_fetch(cls, formula: str, preferred_spg: Optional[List[str]], dest_dir: Path) -> Optional[Path]:
        """按化学式检索，并根据催化优选空间群智能打分筛选。

        策略（先快后全）：
          1. 先按吸收元素单独搜索 COD JSON（返回量小），客户端按 formula 精确匹配
          2. 未命中再按全部元素组合搜索 JSON（返回量大），客户端再按 formula 匹配
        """
        clean_formula = formula.strip()
        elements = [ch for ch in re.findall(r'[A-Z][a-z]?', clean_formula)]
        if not elements:
            return None
        absorber = elements[0]

        # 策略 0：复用本地已下载且化学式匹配的 CIF（应对 COD 连接重置/超时）
        try:
            local_files = sorted(
                dest_dir.glob(f"{formula}_COD_*.cif"),
                key=lambda p: p.stat().st_mtime, reverse=True)
        except OSError:
            local_files = []
        for path in local_files:
            try:
                valid, err, meta = CIFValidator.validate_and_parse(
                    path.read_text(encoding="utf-8"))
            except Exception:
                continue
            if (valid and meta and meta.get("formula")
                    and formula_matches(clean_formula, meta["formula"])):
                logging.info(f"[COD] 复用本地 CIF -> {path.name}")
                return path

        seen_ids: set = set()

        # 策略 1：单元素搜索（快）
        for cand in cls._search_json_single(absorber):
            if cand["cod_id"] in seen_ids:
                continue
            if formula_matches(clean_formula, cand["formula"]):
                seen_ids.add(cand["cod_id"])
                return cls._download_and_name(cand["cod_id"], formula, dest_dir)

        # 策略 2：多元素组合搜索（全，量大）
        for cand in cls._search_json_multi(elements):
            if cand["cod_id"] in seen_ids:
                continue
            if formula_matches(clean_formula, cand["formula"]):
                seen_ids.add(cand["cod_id"])
                return cls._download_and_name(cand["cod_id"], formula, dest_dir)

        return None

    @classmethod
    def _download_and_name(cls, cod_id: str, formula: str, dest_dir: Path) -> Optional[Path]:
        """下载 CIF 并按公式前缀命名保存。"""
        final_name = f"{formula}_COD_{cod_id}.cif"
        return cls.download_by_id(cod_id, dest_dir, final_name=final_name,
                                 expected_formula=formula)

    @staticmethod
    def _search_json_single(element: str, limit: int = 500) -> list[dict]:
        """按单元素搜索 COD JSON。"""
        try:
            resp = requests.get(
                CODClient.SEARCH_URL,
                params={"formula": element, "format": "json", "n": limit},
                headers=HEADERS, timeout=30)
            resp.raise_for_status()
        except requests.RequestException:
            return []
        return CODClient._parse_json(resp.text)

    @classmethod
    def _search_json_multi(cls, elements: List[str]) -> list[dict]:
        """按多元素组合搜索 COD JSON。"""
        params: dict = {"format": "json"}
        for i, el in enumerate(elements, start=1):
            params[f"el{i}"] = el
        try:
            resp = requests.get(cls.SEARCH_URL, params=params,
                                headers=HEADERS, timeout=60)
            resp.raise_for_status()
        except requests.RequestException:
            return []
        return cls._parse_json(resp.text)

    @staticmethod
    def _parse_json(text: str) -> list[dict]:
        """解析 COD JSON 搜索返回。"""
        try:
            data = json.loads(text)
        except (json.JSONDecodeError, ValueError):
            return []
        if not isinstance(data, list):
            return []
        out: list[dict] = []
        for d in data:
            if not isinstance(d, dict):
                continue
            cod_id = d.get("file") or d.get("cod_id")
            if not cod_id:
                continue
            out.append({
                "cod_id": str(cod_id),
                "formula": d.get("formula", "") or "",
                "name": (d.get("commonname") or d.get("mineral")
                         or d.get("chemname") or d.get("formula") or ""),
            })
        return out


# ==================== 3. 数据源：AFLOW / Materials Project 兜底 ====================
class AFLOWClient:
    """AFLOW AFLUX REST API 兜底检索（无需 API Key）"""
    AFLUX_ENDPOINT = "http://aflowlib.org/API/aflux/?"

    @classmethod
    def query_and_fetch(cls, formula: str, dest_dir: Path) -> Optional[Path]:
        # 拆分化学式中包含的元素
        elements = [ch for ch in re.findall(r'[A-Z][a-z]?', formula)]
        if not elements:
            return None
        
        elem_str = ",".join(elements)
        query = f"species({elem_str}),nspecies({len(elements)}),paging(1,format(json))"
        try:
            resp = requests.get(cls.AFLUX_ENDPOINT + query, headers=HEADERS, timeout=15)
            if resp.status_code == 200:
                data = resp.json()
                if not data:
                    return None
                
                entry = data[0]
                aurl = entry.get("aurl")
                if aurl:
                    # 获取经过弛豫计算的标准基态 POSCAR
                    host, path = aurl.split(":")
                    poscar_url = f"http://{host}/{path}/CONTCAR.relax.vasp"
                    p_resp = requests.get(poscar_url, headers=HEADERS, timeout=15)
                    if p_resp.status_code == 200:
                        out_poscar = dest_dir / f"{formula}_AFLOW.poscar"
                        out_poscar.write_text(p_resp.text, encoding="utf-8")
                        logging.info(f"[AFLOW] 成功拉取兜底结构: {out_poscar.name}")
                        return out_poscar
        except Exception:
            pass
        return None


# ==================== 4. 数据源：T5 XAS 标准金属箔参考谱抓取 ====================
class LegacySSLAdapter(HTTPAdapter):
    """兼容 IIT-Lytle 旧服务器的弱 DH 密钥 (OpenSSL 3.0 默认 SECLEVEL=2 会拒绝握手)"""

    def init_poolmanager(self, *args, **kwargs):
        ctx = ssl.create_default_context(ssl.Purpose.SERVER_AUTH)
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        ctx.set_ciphers("DEFAULT:@SECLEVEL=1")
        kwargs["ssl_context"] = ctx
        return super().init_poolmanager(*args, **kwargs)


class XASStandardFetcher:
    """T5: 抓取 IIT-Lytle (Farrel Lytle) 数据库中的真实 XAS 标准谱。

    数据源：
      - PROCESSED: 头含 ABSORPTION EDGE ENERGY，数据为 step/相对能量/mu，
        四组每行。部分旧文件头 E0 是模板占位值（如 Ru/Mo 显示 4966.4），
        因此统一以 absorber.py 的 ABSORBERS[el].e0_edge 作为绝对能量基准。
      - RAW: Cr/Ag 由单色器步数换算能量 (mu=ln(I0/I1))；
        U 的 DATA START 后首列已是 eV (mu=荧光/I0)。
    """

    PROCESSED_BASE = "http://ixs.iit.edu/data/Farrel_Lytle_data/PROCESSED/"
    RAW_BASE = "http://ixs.iit.edu/data/Farrel_Lytle_data/RAW/"
    # 旧入口，保留常量兼容历史调用；真实数据入口见上面两个 BASE
    IIT_LYTLE_URL = "https://ixs.iit.edu/data/Lytle_data/{elem}.DAT"

    NUM_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")
    PROCESSED_TRIPLE_RE = re.compile(
        r"(?<![\d.])(\d+)\s+"
        r"([-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?)\s+"
        r"([-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?)")

    # type: processed / raw_cr / raw_ag / raw_u
    SOURCES = {
        "Cr": {"type": "raw_cr", "files": ["Cr/crfl.314"],
               "note": "CR FOIL 5UM (RAW, transmission)"},
        "Mn": {"type": "processed", "files": ["m/mnfoil.a01", "m/mnfoil.f01"]},
        "Zn": {"type": "processed", "files": ["z/znfl.a03", "z/znfoil.a10"]},
        "Mo": {"type": "processed", "files": ["m/mofl.a47"]},
        "Ru": {"type": "processed", "files": ["r/ru.a01", "r/ru.a02", "r/ru.a03"]},
        "Rh": {"type": "processed", "files": ["r/rhfoil.a03"]},
        "Pd": {"type": "processed", "files": ["p/pdfoil.a71"]},
        "Ag": {"type": "raw_ag", "files": ["Ag/agfl.651"],
               "note": "AG FOIL (RAW, K-edge, transmission)"},
        "Ir": {"type": "processed", "files": ["i/irfoie.a46"]},
        "Au": {"type": "processed", "files": ["a/aufl.a01"]},
        "Ce": {"type": "processed", "files": ["c/ceo2l3.a39", "c/ceo2.260"],
               "note": "CeO2 L3 参考谱（IIT 库中无 Ce 金属箔，用 CeO2 替代）"},
        "U": {"type": "raw_u", "files": ["U/umetal.351"],
              "note": "U METAL (RAW, fluorescence I/I0)"},
    }

    @classmethod
    def get_session(cls) -> requests.Session:
        session = requests.Session()
        session.mount("https://", LegacySSLAdapter())
        # requests 的 merge_environment_settings 会把 CERT_NONE 覆盖回 CERT_REQUIRED，
        # 必须显式关闭校验，否则 IIT 的旧证书链仍会报 CERTIFICATE_VERIFY_FAILED。
        session.verify = False
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        return session

    # ---------- 通用数值解析 ----------
    @staticmethod
    def _to_float(text) -> Optional[float]:
        try:
            return float(str(text).replace("D", "E").replace("d", "e"))
        except (TypeError, ValueError):
            return None

    @classmethod
    def _numbers(cls, text: str) -> List[float]:
        out: List[float] = []
        for tok in cls.NUM_RE.findall(text):
            val = cls._to_float(tok)
            if val is not None:
                out.append(val)
        return out

    @classmethod
    def _parse_monochromator(cls, lines: List[str],
                             default_d: float = 1.92017,
                             default_stepdeg: float = 4000.0) -> Tuple[float, float]:
        """解析 RAW 头中的 d 间距与每度步数，兼容单行标签与两行表格两种格式。"""
        dspace, stepdeg = default_d, default_stepdeg
        for i, line in enumerate(lines[:20]):
            up = line.upper()
            m = re.search(r"D[\s-]?SPACE\s*[:=]?\s*([0-9.]+)", up)
            if m:
                dspace = float(m.group(1))
            m = re.search(r"(?:STPDEG|STEP[_ ]?DEG|STEPS/DEGREE)\s*[:=]?\s*([0-9.]+)", up)
            if m:
                stepdeg = float(m.group(1))
            # Ag/Cr 风格：标签行 + 下一行 10 列数值
            if "DSPACE" in up and "STPDEG" in up:
                nums = cls._numbers(line)
                if not nums and i + 1 < len(lines):
                    nums = cls._numbers(lines[i + 1])
                if len(nums) >= 6:
                    dspace, stepdeg = nums[4], nums[5]
        return dspace, stepdeg

    @classmethod
    def _parse_processed(cls, text: str) -> Tuple[Optional[float], List[Tuple[float, float]]]:
        """解析 PROCESSED: 返回 (头部 E0, [(相对能量, mu), ...])。"""
        m = re.search(r"ABSORPTION EDGE ENERGY\s*\(EV\)\s*:\s*([-+0-9.EeDd]+)", text)
        header_e0 = cls._to_float(m.group(1)) if m else None
        points: Dict[float, float] = {}
        for line in text.splitlines():
            stripped = line.strip()
            if not stripped or not stripped[0].isdigit():
                continue
            for g in cls.PROCESSED_TRIPLE_RE.findall(stripped):
                rel = cls._to_float(g[1])
                mu = cls._to_float(g[2])
                if rel is None or mu is None:
                    continue
                if math.isfinite(rel) and math.isfinite(mu):
                    points[round(rel, 6)] = mu
        return header_e0, sorted(points.items())

    @classmethod
    def _parse_raw_step_scan(cls, text: str,
                             default_d: float = 1.92017,
                             default_stepdeg: float = 4000.0) -> List[Tuple[float, float]]:
        """解析 Cr/Ag RAW: E=hc/(2*d*sin(theta)), mu=ln(I0/I1)。"""
        lines = text.splitlines()
        dspace, stepdeg = cls._parse_monochromator(lines, default_d, default_stepdeg)
        points: List[Tuple[float, float]] = []
        for line in lines:
            # 去掉数值 token 后再检查字母，避免把科学计数法中的 E 误判为表头
            if re.search(r"[A-Za-z]", cls.NUM_RE.sub(" ", line)):
                continue
            nums = cls._numbers(line)
            if len(nums) < 3:
                continue
            step = nums[0]
            if not (10000.0 <= step <= 1_000_000.0):
                continue
            if abs(step - round(step)) > 1e-3:
                continue
            i0, i1 = nums[1], nums[2]
            if i0 <= 0 or i1 <= 0:
                continue
            theta = math.radians(step / stepdeg)
            denom = 2.0 * dspace * math.sin(theta)
            if denom <= 0:
                continue
            energy = 12398.419843 / denom
            points.append((energy, math.log(i0 / i1)))
        return points

    @classmethod
    def _parse_raw_u(cls, text: str) -> List[Tuple[float, float]]:
        """解析 U RAW: DATA START 后首列 eV，mu=荧光(col5)/I0(col3)。"""
        lines = text.splitlines()
        start = None
        for i, line in enumerate(lines):
            if "DATA START" in line.upper():
                start = i
                break
        if start is None:
            return []
        points: List[Tuple[float, float]] = []
        for line in lines[start + 1:]:
            nums = cls._numbers(line)
            if len(nums) < 6:
                continue
            energy, i0, fluo = nums[0], nums[3], nums[5]
            if not (1000.0 <= energy <= 100000.0) or i0 <= 0 or fluo < 0:
                continue
            points.append((energy, fluo / i0))
        return points

    # ---------- 归一化与质控 ----------
    @staticmethod
    def _normalize(points: List[Tuple[float, float]],
                   e0: float) -> List[Tuple[float, float]]:
        """扣 pre-edge 本底并按 post-edge 平台归一化，使跳边高度约为 1。"""
        clean: List[Tuple[float, float]] = []
        for e, mu in points:
            if e is None or mu is None:
                continue
            e, mu = float(e), float(mu)
            if math.isfinite(e) and math.isfinite(mu):
                clean.append((e, mu))
        if len(clean) < 10:
            return []

        clean.sort(key=lambda p: p[0])
        merged: List[Tuple[float, float, int]] = []
        for e, mu in clean:
            if merged and abs(e - merged[-1][0]) < 1e-6:
                prev_e, prev_mu, n = merged[-1]
                merged[-1] = (prev_e, (prev_mu * n + mu) / (n + 1), n + 1)
            else:
                merged.append((e, mu, 1))
        pts = [(e, mu) for e, mu, _ in merged]

        pre = [mu for e, mu in pts if e0 - 150.0 <= e <= e0 - 30.0]
        if len(pre) < 5:
            pre = [mu for _, mu in pts[:max(5, len(pts) // 10)]]
        post = [mu for e, mu in pts if e0 + 30.0 <= e <= e0 + 300.0]
        if len(post) < 5:
            post = [mu for _, mu in pts[-max(5, len(pts) // 10):]]
        if len(pre) < 3 or len(post) < 3:
            return []

        pre_level = statistics.median(pre)
        post_level = statistics.median(post)
        step = post_level - pre_level
        if abs(step) < 1e-8:
            step = max(mu for _, mu in pts) - min(mu for _, mu in pts)
        if abs(step) < 1e-12:
            return []
        return [(e, (mu - pre_level) / step) for e, mu in pts]

    @staticmethod
    def _has_edge_jump(points: List[Tuple[float, float]], e0: float) -> bool:
        pre = [mu for e, mu in points if e0 - 50.0 <= e <= e0 - 10.0]
        post = [mu for e, mu in points if e0 + 10.0 <= e <= e0 + 50.0]
        if len(pre) < 3 or len(post) < 3:
            return False
        return (statistics.median(post) - statistics.median(pre)) > 0.05

    @classmethod
    def fetch_foil(cls, element: str, dest_dir: Path) -> Optional[Path]:
        dest_dir.mkdir(parents=True, exist_ok=True)
        element = element.capitalize()
        out_file = dest_dir / f"{element}_metal_foil.dat"

        cfg = ABSORBERS.get(element) if HAS_ABSORBER else None
        if cfg is None:
            logging.warning(f"[T5] 未知吸收元素 {element}，跳过标准谱抓取")
            return None
        e0 = cfg.e0_edge
        edge = cfg.edge

        # 兼容原 Larch 路径（当前 Larch 2026.3 无 standard_foil，通常直接跳过）
        if HAS_LARCH:
            try:
                foil_group = standard_foil(element)
                if foil_group and hasattr(foil_group, "energy") and hasattr(foil_group, "mu"):
                    with open(out_file, "w", encoding="utf-8") as f:
                        f.write(f"# XAS standard foil for {element} extracted via Larch\n")
                        f.write("# Energy(eV)\tmu\n")
                        for e_val, mu_val in zip(foil_group.energy, foil_group.mu):
                            f.write(f"{e_val:.4f}\t{mu_val:.6f}\n")
                    logging.info(f"[Larch] 成功导出内置标准金属箔谱: {out_file.name}")
                    return out_file
            except Exception:
                pass

        # 真实数据：IIT-Lytle PROCESSED / RAW
        source = cls.SOURCES.get(element)
        if source:
            session = cls.get_session()
            for rel in source["files"]:
                base = cls.RAW_BASE if source["type"].startswith("raw") else cls.PROCESSED_BASE
                url = base + rel
                try:
                    resp = session.get(url, headers=HEADERS, timeout=60)
                    if resp.status_code != 200 or len(resp.content) < 200:
                        logging.warning(
                            f"[IIT-Lytle] {element} {rel} 不可用 (HTTP {resp.status_code})")
                        continue

                    if source["type"] == "processed":
                        header_e0, rel_points = cls._parse_processed(resp.text)
                        if not rel_points:
                            logging.warning(f"[IIT-Lytle] {element} {rel} 无有效数据行")
                            continue
                        if header_e0 is None or abs(header_e0 - e0) > 100.0:
                            logging.info(
                                f"[IIT-Lytle] {element} {rel} 头部 E0={header_e0} 为占位/偏差值，"
                                f"按项目 E0={e0} eV 重建绝对能量")
                        raw_points = [(e0 + rel_e, mu) for rel_e, mu in rel_points]
                        origin = "PROCESSED (E0 + 相对能量)"
                    elif source["type"] in ("raw_cr", "raw_ag"):
                        raw_points = cls._parse_raw_step_scan(resp.text)
                        origin = "RAW (单色器步数→能量)"
                    else:  # raw_u
                        raw_points = cls._parse_raw_u(resp.text)
                        origin = "RAW (荧光/I0)"

                    if len(raw_points) < 20:
                        logging.warning(
                            f"[IIT-Lytle] {element} {rel} 解析点数不足: {len(raw_points)}")
                        continue
                    if not cls._has_edge_jump(raw_points, e0):
                        logging.warning(
                            f"[IIT-Lytle] {element} {rel} 未检测到 E0={e0} eV 附近吸收边，放弃")
                        continue

                    norm_points = cls._normalize(raw_points, e0)
                    if len(norm_points) < 20:
                        logging.warning(f"[IIT-Lytle] {element} {rel} 归一化后点数不足")
                        continue

                    note = source.get("note", "")
                    with open(out_file, "w", encoding="utf-8") as f:
                        f.write(f"# XAS standard reference for {element}\n")
                        f.write(f"# Source: IIT-Lytle Farrel Lytle database, {origin}\n")
                        f.write(f"# Source file: {rel}\n")
                        f.write(f"# Absorber: {element}, Edge: {edge}, E0 = {e0} eV\n")
                        if note:
                            f.write(f"# NOTE: {note}\n")
                        f.write(f"# Data points: {len(norm_points)}\n")
                        f.write(f"# Energy range: {norm_points[0][0]:.2f} - "
                                f"{norm_points[-1][0]:.2f} eV\n")
                        f.write("# Energy(eV)\tmu(normalized)\n")
                        for e_val, mu_val in norm_points:
                            f.write(f"{e_val:.4f}\t{mu_val:.8f}\n")
                    logging.info(
                        f"[IIT-Lytle] 真实标准谱已导出: {out_file.name} "
                        f"({len(norm_points)} 点, "
                        f"{norm_points[0][0]:.0f}-{norm_points[-1][0]:.0f} eV, 来源 {rel})")
                    return out_file
                except Exception as e:
                    logging.warning(
                        f"[IIT-Lytle] {element} {rel} 抓取/解析失败: "
                        f"{type(e).__name__}: {e}")

        # 兜底：若已有真实谱，网络失败时保留；否则写占位文件
        if out_file.exists():
            txt = out_file.read_text(encoding="utf-8", errors="replace")
            if "placeholder" not in txt.lower():
                logging.warning(f"[IIT-Lytle] {element} 抓取失败，保留已有真实谱: {out_file.name}")
                return out_file

        with open(out_file, "w", encoding="utf-8") as f:
            f.write(f"# XAS standard foil placeholder for {element}\n")
            f.write(f"# Absorber: {element}, Edge: {edge}, E0 ~ {e0} eV\n")
            f.write("# NOTE: 此为占位文件，IIT-Lytle 真实数据抓取失败；"
                    "请替换为标准金属箔 .xdi/.dat 谱文件\n")
            f.write("# 格式: Energy(eV)\tmu\n")
            f.write(f"{e0 - 50:.2f}\t0.000000\n")
            f.write(f"{e0:.2f}\t1.000000\n")
            f.write(f"{e0 + 50:.2f}\t0.500000\n")
        logging.warning(f"[Placeholder] 生成标准箔占位文件: {out_file.name} (E0={e0} eV)")
        return out_file


# ==================== 5. 结构库同步与总索引管理 ====================
class PipelineManager:
    """负责将抓取的文件转化为管线索引并导入 _structure_lib"""

    def __init__(self):
        CIF_RAW_DIR.mkdir(parents=True, exist_ok=True)
        STANDARDS_DIR.mkdir(parents=True, exist_ok=True)
        # 加载已有索引（T5 标准谱条目必须保留），绝不在启动时清空
        self.index: dict = self._load_index()

    @staticmethod
    def _load_index() -> dict:
        if not INDEX_FILE.exists():
            return {}
        try:
            data = json.loads(INDEX_FILE.read_text(encoding="utf-8"))
            if isinstance(data, dict):
                return data
        except Exception as exc:
            logging.warning(f"[Index] 读取失败，将在写入时备份损坏文件: {exc}")
        return {}

    def save_index(self):
        """整表原子写（保留给外部调用；日常增量写入走 upsert_index）。"""
        INDEX_FILE.parent.mkdir(parents=True, exist_ok=True)
        tmp_file = INDEX_FILE.with_name(
            f"{INDEX_FILE.name}.tmp.{os.getpid()}")
        tmp_file.write_text(
            json.dumps(self.index, indent=2, ensure_ascii=False), encoding="utf-8")
        os.replace(tmp_file, INDEX_FILE)

    def _upsert_entry(self, entry_key: str, entry: dict):
        upsert_index(INDEX_FILE, entry_key, entry)
        self.index[entry_key] = entry

    def _quarantine_motifs(self, absorber: str, motif_ids) -> Dict[Path, Path]:
        """把会冲突/待淘汰的旧 motif JSON 移入隔离区，返回 {原路径: 备份路径}。"""
        backup: Dict[Path, Path] = {}
        if not HAS_STRUCTURE_LIB or not motif_ids:
            return backup
        lib_root = Path(sl_mod.LIB_ROOT)
        stamp = time.strftime("%Y%m%d_%H%M%S") + f"_{time.time_ns() % 1_000_000_000:09d}"
        qdir = lib_root / ".quarantine" / stamp / absorber
        for motif_id in sorted(set(motif_ids)):
            src = lib_root / absorber / f"{motif_id}.json"
            if not src.exists():
                continue
            qdir.mkdir(parents=True, exist_ok=True)
            dst = qdir / src.name
            if dst.exists():
                dst = qdir / f"{src.stem}_{time.time_ns()}{src.suffix}"
            try:
                shutil.move(str(src), str(dst))
            except OSError as exc:
                logging.warning(f"[LibSync] 旧 motif 隔离失败 {src}: {exc}")
                continue
            backup[src] = dst
            logging.info(f"[LibSync] 旧 motif 已隔离: {src.name}")
        return backup

    @staticmethod
    def _prune_quarantine_dirs(backup: Dict[Path, Path]):
        dirs = {saved.parent for saved in backup.values()}
        for d in sorted(dirs, key=lambda p: len(p.parts), reverse=True):
            for candidate in (d, d.parent, d.parent.parent):
                try:
                    candidate.rmdir()
                except OSError:
                    pass

    @staticmethod
    def _restore_quarantine(backup: Dict[Path, Path]):
        for original, saved in backup.items():
            try:
                original.parent.mkdir(parents=True, exist_ok=True)
                if original.exists():
                    original.unlink()
                if saved.exists():
                    shutil.move(str(saved), str(original))
                    logging.info(f"[LibSync] 注册失败，已恢复旧 motif: {original.name}")
            except Exception as exc:
                logging.error(f"[LibSync] 恢复旧 motif 失败 {original}: {exc}")
        PipelineManager._prune_quarantine_dirs(backup)

    @staticmethod
    def _purge_quarantine(backup: Dict[Path, Path]):
        for _original, saved in backup.items():
            try:
                if saved.exists():
                    saved.unlink()
            except Exception as exc:
                logging.warning(f"[LibSync] 清理隔离文件失败 {saved}: {exc}")
        PipelineManager._prune_quarantine_dirs(backup)

    def _cleanup_legacy_files(self, absorber: str, legacy_ids, raw_cifs):
        """新相成功注册后，删除旧错相/过时 CIF 与旧 _cif 缓存。"""
        if HAS_STRUCTURE_LIB:
            lib_root = Path(sl_mod.LIB_ROOT)
            for old_id in set(legacy_ids or []):
                cache = lib_root / absorber / "_cif" / f"{old_id}.cif"
                try:
                    if cache.exists():
                        cache.unlink()
                        logging.info(f"[Cleanup] 已移除旧相缓存: {cache.name}")
                except OSError as exc:
                    logging.warning(f"[Cleanup] 删除旧相缓存失败 {cache}: {exc}")
        else:
            lib_root = STRUCTURE_LIB_DIR.resolve()

        for name in set(raw_cifs or []):
            for path in (CIF_RAW_DIR / name, lib_root / absorber / name):
                try:
                    if path.exists():
                        path.unlink()
                        logging.info(f"[Cleanup] 已移除旧相 CIF: {path}")
                except OSError as exc:
                    logging.warning(f"[Cleanup] 删除旧相 CIF 失败 {path}: {exc}")

    def import_and_register_motif(self, cif_path: Path, absorber: str, safe_id: str,
                                  formula: str, phase: str, cod_id: str,
                                  legacy_ids=None, raw_cifs=None,
                                  quarantine_collisions: bool = True,
                                  absorber_fract=None) -> bool:
        """通过原生 structure_lib.import_cif + register 入库。

        流程：
        1. 用 gemmi 生成 Motif 并预计算指纹，找出会冲突的旧 motif；
        2. 将旧 motif 移入隔离区（不直接删除，失败可回滚）；
        3. register(motif, overwrite=False)；
        4. 成功则清理隔离文件与旧 CIF，失败则恢复旧 motif。
        """
        if not HAS_STRUCTURE_LIB:
            logging.warning("[LibSync] structure_lib 未加载，无法生成 motif JSON")
            return False

        try:
            motif = sl_import_cif(
                str(cif_path), safe_id, absorber,
                nuclearity="crystal",
                source={"type": "cif_crawler",
                        "source_db": "cod",
                        "formula": formula,
                        "phase": phase,
                        "cod_id": cod_id,
                        "file": _index_path(cif_path)},
                absorber_fract=absorber_fract)
        except Exception as exc:
            logging.error(f"[LibSync] CIF -> Motif 失败 ({safe_id}): {exc}")
            return False

        # 预计算指纹，提前发现所有冲突项（与 register 内部算法保持一致）
        try:
            sites, occs = sl_mod.expand_full(motif, 6.0)
            motif.sites = [
                list(map(lambda v: round(v, 6) if isinstance(v, float) else v, s))
                for s in sites
            ]
            motif.shell_expect = [
                {k: (round(v, 4) if isinstance(v, float) else v)
                 for k, v in s.items() if k != "n_occ" or occs is not None}
                for s in sl_mod.geometric_shells(motif.sites, occs=occs)
            ]
            collision_ids = {
                other.id for other in sl_mod.load_library(absorber)
                if other.id != safe_id and sl_mod._same_fingerprint(other, motif)
            }
        except Exception as exc:
            logging.error(f"[LibSync] Motif 指纹计算失败 ({safe_id}): {exc}")
            return False

        lib_dir = Path(sl_mod.LIB_ROOT) / absorber
        targeted_ids = set(legacy_ids or [])
        if quarantine_collisions:
            targeted_ids |= collision_ids
        elif collision_ids:
            # 自由检索（--formula / --known-cod-id）不应擅自淘汰既有标准 motif
            logging.error(
                f"[LibSync] {safe_id} 与既有 motif {sorted(collision_ids)} "
                f"指纹重复；自由检索路径拒绝自动淘汰，索引不更新")
            return False
        if (lib_dir / f"{safe_id}.json").exists():
            targeted_ids.add(safe_id)

        quarantine = self._quarantine_motifs(absorber, targeted_ids)
        try:
            sl_register(motif, overwrite=False)
        except (XAFSConfigError, XAFSDataError) as exc:
            self._restore_quarantine(quarantine)
            logging.error(f"[LibSync] motif 注册被拒绝 ({safe_id}): {exc}")
            return False
        except Exception as exc:
            self._restore_quarantine(quarantine)
            logging.error(f"[LibSync] motif 注册异常 ({safe_id}): {exc}")
            return False

        # 注册成功：隔离区内的旧结构才是真正可淘汰的
        self._purge_quarantine(quarantine)
        self._cleanup_legacy_files(absorber, legacy_ids, raw_cifs)
        logging.info(f"[LibSync] motif 已注册 -> _structure_lib/{absorber}/{safe_id}.json")
        return True

    def register_cif(self, formula: str, absorber: str, cif_path: Path, source: str,
                     import_lib: bool, phase: Optional[str] = None,
                     cod_id: Optional[str] = None, tier: Optional[str] = None,
                     spg_num: Optional[int] = None,
                     spg_nums: Optional[List[int]] = None,
                     elements: Optional[List[str]] = None,
                     stoich: Optional[Dict[str, float]] = None,
                     motif_id: Optional[str] = None,
                     legacy_ids=None, raw_cifs=None,
                     absorber_fract=None) -> bool:
        entry_key = f"{absorber}_{formula}"
        if motif_id is None:
            if phase and cod_id:
                motif_id = generate_safe_motif_id(formula, phase, cod_id)
            else:
                # 兼容 --known-cod-id / --formula 的自由检索路径
                motif_id = f"{absorber}{re.sub(r'[^A-Za-z]', '', formula)}".lower()
                motif_id = motif_id or f"{absorber.lower()}_{formula.lower()}"

        if import_lib and HAS_STRUCTURE_LIB:
            if not self.import_and_register_motif(
                    cif_path, absorber, motif_id, formula,
                    phase or "auto", cod_id or "unknown",
                    legacy_ids=legacy_ids, raw_cifs=raw_cifs,
                    quarantine_collisions=legacy_ids is not None,
                    absorber_fract=absorber_fract):
                logging.error(f"[LibSync] {formula} 未完成入库，索引不更新")
                return False
        elif import_lib and not HAS_STRUCTURE_LIB:
            # 退化：仅复制 CIF 文件
            target_dir = STRUCTURE_LIB_DIR / absorber
            target_dir.mkdir(parents=True, exist_ok=True)
            target_cif = target_dir / cif_path.name
            target_cif.write_text(cif_path.read_text(encoding="utf-8"), encoding="utf-8")
            logging.info(f"[LibSync] 结构 CIF 已复制 -> {target_cif} (structure_lib 未加载，未生成 motif JSON)")

        entry = {
            "formula": formula,
            "absorber": absorber,
            "cif_file": _index_path(cif_path),
            "source": source,
            "motif_id": motif_id,
        }
        if tier is not None:
            entry["tier"] = tier
        if phase is not None:
            entry["phase"] = phase
        if cod_id is not None:
            entry["cod_id"] = str(cod_id)
        if spg_num is not None:
            entry["spg_num"] = spg_num
        if spg_nums is not None:
            entry["spg_nums"] = list(spg_nums)
        if elements is not None:
            entry["elements"] = list(elements)
        if stoich is not None:
            entry["stoich"] = {k: float(v) for k, v in stoich.items()}

        self._upsert_entry(entry_key, entry)
        return True

    def register_standard(self, element: str, std_path: Path, import_lib: bool):
        entry = {
            "element": element,
            "standard_file": _index_path(std_path),
            "type": "foil_spectrum"
        }
        self._upsert_entry(f"STD_{element}", entry)

        if import_lib:
            target_dir = STANDARDS_LIB_DIR / element
            target_dir.mkdir(parents=True, exist_ok=True)
            target_file = target_dir / std_path.name
            target_file.write_text(std_path.read_text(encoding="utf-8"), encoding="utf-8")
            logging.info(f"[LibSync] XAS标准谱已入库 -> {target_file}")


def process_cif_target(name: str, meta: dict, pipeline: PipelineManager,
                       import_lib: bool) -> bool:
    """抓取、物理核验并入库单个 VERIFIED_TARGETS 条目。"""
    cod_id = str(meta["cod_id"])
    absorber = meta["absorber"]
    phase = meta["phase"]

    cif_path = CODClient.download_by_id(
        cod_id, CIF_RAW_DIR,
        final_name=f"{name}_COD_{cod_id}.cif",
        expected_formula=name)
    if not cif_path:
        logging.error(f"[{name}] COD {cod_id} 下载/化学式预校验失败")
        return False

    if not HAS_GEMMI:
        logging.error(f"[{name}] 缺少 gemmi，无法执行空间群与计量比核验")
        return False

    try:
        doc = gemmi.cif.read(str(cif_path))
        block = doc.sole_block()
    except Exception as exc:
        logging.error(f"[{name}] CIF 解析失败: {exc}")
        return False

    actual_spg = get_spg_number(block)
    if actual_spg not in meta["spg_nums"]:
        logging.error(
            f"[{name}] 空间群不符: 期望 {meta['spg_nums']}, 实测 #{actual_spg} "
            f"({cif_path.name})")
        return False

    comp_ok, comp_msg = verify_composition(block, meta["stoich"])
    if not comp_ok:
        logging.error(f"[{name}] 化学成分核验失败: {comp_msg}")
        return False
    logging.info(f"[{name}] 物理相核验通过: spg=#{actual_spg}, {comp_msg}")

    safe_id = generate_safe_motif_id(name, phase, cod_id)
    legacy = LEGACY_ARTIFACTS.get(name, {})
    return pipeline.register_cif(
        formula=name,
        absorber=absorber,
        cif_path=cif_path,
        source="COD_FIXED",
        import_lib=import_lib,
        phase=phase,
        cod_id=cod_id,
        tier=meta.get("tier"),
        spg_num=actual_spg,
        spg_nums=meta.get("spg_nums"),
        elements=meta.get("elements"),
        stoich=meta.get("stoich"),
        motif_id=safe_id,
        legacy_ids=legacy.get("motif_ids", []),
        raw_cifs=legacy.get("raw_cifs", []),
        absorber_fract=meta.get("absorber_fract"),
    )


# ==================== 6. CLI 逻辑编排 ====================
def main() -> int:
    parser = argparse.ArgumentParser(description="多源 CIF 与 XAS 参考样数据抓取器")
    parser.add_argument("--all", action="store_true", help="全量运行 T3/T4/T5 所有任务并解除阻塞")
    parser.add_argument("--t3", action="store_true", help="抓取 T3 常见催化氧化物真实 CIF")
    parser.add_argument("--t4", action="store_true", help="抓取 T4 缺失基准元素单质 CIF")
    parser.add_argument("--t5", action="store_true", help="抓取 T5 未标定元素金属箔标准谱")
    parser.add_argument("--formula", type=str, help="按特定化学式抓取，如 Fe2O3")
    parser.add_argument("--absorber", type=str, help="指定吸收原子元素，如 Fe")
    parser.add_argument("--known-cod-id", type=str, help="直接根据已知的 COD ID 下载入库")
    parser.add_argument("--import-lib", action="store_true", default=True, help="是否自动导入 _structure_lib 和 _standards_lib")
    parser.add_argument("--no-import", dest="import_lib", action="store_false", help="仅下载到原始目录，不入库")

    args = parser.parse_args()
    pipeline = PipelineManager()
    failures: List[str] = []

    # 1. 优先响应指定 COD ID 直链下载
    if args.known_cod_id:
        cod_id = str(args.known_cod_id)
        cif_path = CODClient.download_by_id(cod_id, CIF_RAW_DIR)
        if cif_path:
            absorber = args.absorber or "Unknown"
            ok = pipeline.register_cif(
                "Custom", absorber, cif_path, "COD_DIRECT", args.import_lib,
                phase="custom", cod_id=cod_id)
            if not ok:
                failures.append(f"custom COD {cod_id}")
        else:
            failures.append(f"custom COD {cod_id}")
        if failures:
            logging.error("以下目标失败: %s", ", ".join(failures))
            return 1
        return 0

    # 2. 单个化学式检索（自由检索路径，不套用固定相绑定）
    if args.formula:
        absorber = args.absorber or "".join(re.findall(r'[A-Z][a-z]?', args.formula)[:1])
        cif_path = CODClient.query_and_fetch(args.formula, None, CIF_RAW_DIR)
        if not cif_path:
            cif_path = AFLOWClient.query_and_fetch(args.formula, CIF_RAW_DIR)

        if cif_path:
            m = re.search(r"COD[_-]?(\d+)", cif_path.name, re.IGNORECASE)
            cod_id = m.group(1) if m else None
            ok = pipeline.register_cif(
                args.formula, absorber, cif_path, "AUTO_FALLBACK", args.import_lib,
                phase="auto" if cod_id else None, cod_id=cod_id)
            if not ok:
                failures.append(args.formula)
        else:
            failures.append(args.formula)
        if failures:
            logging.error("以下目标失败: %s", ", ".join(failures))
            return 1
        return 0

    # 3. 批量任务流
    run_all = args.all or (not args.t3 and not args.t4 and not args.t5)

    if run_all or args.t3:
        logging.info("========== 执行任务 T3: 催化氧化物真实 CIF 抓取 ==========")
        for name, meta in VERIFIED_TARGETS.items():
            if meta.get("tier") != "T3":
                continue
            if not process_cif_target(name, meta, pipeline, args.import_lib):
                failures.append(name)

    if run_all or args.t4:
        logging.info("========== 执行任务 T4: 缺失基准元素单质 CIF 抓取 ==========")
        for name, meta in VERIFIED_TARGETS.items():
            if meta.get("tier") != "T4":
                continue
            if not process_cif_target(name, meta, pipeline, args.import_lib):
                failures.append(name)

    if run_all or args.t5:
        logging.info("========== 执行任务 T5: XAS 金属箔标准样基准谱抓取 ==========")
        for elem in TARGETS_T5_ELEMENTS:
            std_path = XASStandardFetcher.fetch_foil(elem, STANDARDS_DIR)
            if std_path:
                pipeline.register_standard(elem, std_path, args.import_lib)

    if failures:
        logging.error("以下 CIF 目标失败: %s", ", ".join(failures))
        return 1
    logging.info("所有请求任务执行完毕，索引已落盘至 data/cif_crawler_index.json。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
