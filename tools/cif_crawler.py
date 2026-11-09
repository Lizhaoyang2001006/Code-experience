# -*- coding: utf-8 -*-
"""cif_crawler.py — 构建期 CIF 抓取与 phase pack 导入

设计边界
========
本工具只服务于**构建期**：把真实、可追溯的 COD CIF 补充进
``_structure_lib/_phase_pack/``。它不接入 ``run_xafs.py`` 的运行路径，
因此现场分析仍然完全离线，不做在线检索。

硬约束
------
1. 只认 ``VERIFIED_TARGETS`` 中固定的 COD ID；不做化学式自由检索，
   不自动选择同质多象，也不使用 AFLOW/POSCAR 兜底。
2. 下载后必须同时通过：
   - gemmi 可解析；
   - IT 空间群号命中目标；
   - ``_chemical_formula_sum`` 的元素集合与化学计量命中目标；
   - 目标吸收元素存在至少一个位点；
   - 不是 ``structure_lib._cif_path`` 现场生成的小 CIF。
3. 原始 CIF 字节先落 ``data/cif_raw`` 缓存，再经
   ``phase_pack.import_cif`` 写入 phase pack；sha256 由 phase pack 记录。
4. 不写 motif JSON，不修改 ``_structure_lib/<元素>/`` 诊断库。
5. ``--online`` 必须显式给出；默认只读本地缓存（离线）。

用法
----
::

    python tools/cif_crawler.py --list
    python tools/cif_crawler.py --t3 --online
    python tools/cif_crawler.py --t4 --online
    python tools/cif_crawler.py --all --online
    python tools/cif_crawler.py --target Fe2O3 --offline

导入成功后默认重建 ``phase_index.json`` 与 ``phase_coverage.json``。
"""
from __future__ import annotations

import argparse
import json
import os
import re
import ssl
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import phase_pack                       # noqa: E402
from motif_provenance import is_generated_cif, sha256_file  # noqa: E402

COD_DOWNLOAD_URL = "https://www.crystallography.net/cod/{cod_id}.cif"
USER_AGENT = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "Catalysis-CIF-Crawler/3.0 (+phase-pack build tool)"
)
DEFAULT_CACHE = ROOT / "data" / "cif_raw"
DEFAULT_REPORT = ROOT / "_structure_lib" / "cif_crawl_report.json"
_COD_ID_RE = re.compile(r"^\d{6,9}$")

# ======================================================================
# 目标映射：COD ID 与目标相已逐条下载、gemmi 验真后固定。
# tier     任务分层（T3 催化常见氧化物/氢氧化物，T4 缺失基准元素单质）
# phase    相名（写入 tags / notes，避免同质多象混淆）
# formula  目标化学式，同时用于缓存命名与计量比核验
# stoich   目标化学计量（用于报告展示）
# spg_nums 该 COD 结构真实对应的 IT 空间群号；只放已验证号码
# tags     写入 phase pack manifest，供人工检索；不参与自动选相
# ======================================================================
VERIFIED_TARGETS = {
    # ---- T3: 催化常见氧化物与氢氧化物 ----
    "Fe2O3": {
        "tier": "T3", "phase": "hematite", "absorber": "Fe",
        "formula": "Fe2O3", "stoich": {"Fe": 2, "O": 3},
        "cod_id": "1011240", "spg_nums": [167],
        "tags": ["oxide", "hematite", "catalyst"],
    },
    "Fe3O4": {
        "tier": "T3", "phase": "magnetite", "absorber": "Fe",
        "formula": "Fe3O4", "stoich": {"Fe": 3, "O": 4},
        "cod_id": "1010369", "spg_nums": [227],
        "tags": ["oxide", "magnetite", "spinel", "catalyst"],
    },
    "CoO": {
        "tier": "T3", "phase": "rocksalt", "absorber": "Co",
        "formula": "CoO", "stoich": {"Co": 1, "O": 1},
        "cod_id": "1533087", "spg_nums": [225],
        "tags": ["oxide", "rocksalt", "catalyst"],
    },
    "Co3O4": {
        "tier": "T3", "phase": "spinel", "absorber": "Co",
        "formula": "Co3O4", "stoich": {"Co": 3, "O": 4},
        "cod_id": "9005888", "spg_nums": [227],
        "tags": ["oxide", "spinel", "catalyst"],
    },
    "NiO": {
        "tier": "T3", "phase": "bunsenite", "absorber": "Ni",
        "formula": "NiO", "stoich": {"Ni": 1, "O": 1},
        "cod_id": "1010093", "spg_nums": [225],
        "tags": ["oxide", "bunsenite", "catalyst"],
    },
    "Ni(OH)2": {
        "tier": "T3", "phase": "theophrastite", "absorber": "Ni",
        "formula": "Ni(OH)2", "stoich": {"Ni": 1, "O": 2, "H": 2},
        "cod_id": "1548811", "spg_nums": [164],
        "tags": ["hydroxide", "theophrastite", "catalyst"],
    },
    "CuO": {
        "tier": "T3", "phase": "tenorite", "absorber": "Cu",
        "formula": "CuO", "stoich": {"Cu": 1, "O": 1},
        "cod_id": "1011148", "spg_nums": [15],
        "tags": ["oxide", "tenorite", "catalyst"],
    },
    "Cu2O": {
        "tier": "T3", "phase": "cuprite", "absorber": "Cu",
        "formula": "Cu2O", "stoich": {"Cu": 2, "O": 1},
        "cod_id": "1010941", "spg_nums": [224],
        "tags": ["oxide", "cuprite", "catalyst"],
    },
    "CeO2": {
        "tier": "T3", "phase": "fluorite", "absorber": "Ce",
        "formula": "CeO2", "stoich": {"Ce": 1, "O": 2},
        "cod_id": "1562989", "spg_nums": [225],
        "tags": ["oxide", "fluorite", "catalyst"],
    },
    # ---- T4: 缺失基准元素单质 ----
    "Mn": {
        "tier": "T4", "phase": "alpha-Mn", "absorber": "Mn",
        "formula": "Mn", "stoich": {"Mn": 1},
        "cod_id": "9008589", "spg_nums": [217],
        "tags": ["metal", "alpha-Mn", "element_seed"],
    },
    "Ru": {
        "tier": "T4", "phase": "hcp", "absorber": "Ru",
        "formula": "Ru", "stoich": {"Ru": 1},
        "cod_id": "9008513", "spg_nums": [194],
        "tags": ["metal", "hcp", "element_seed"],
    },
    "Rh": {
        "tier": "T4", "phase": "fcc", "absorber": "Rh",
        "formula": "Rh", "stoich": {"Rh": 1},
        "cod_id": "9008482", "spg_nums": [225],
        "tags": ["metal", "fcc", "element_seed"],
    },
    "Pd": {
        "tier": "T4", "phase": "fcc", "absorber": "Pd",
        "formula": "Pd", "stoich": {"Pd": 1},
        "cod_id": "9008478", "spg_nums": [225],
        "tags": ["metal", "fcc", "element_seed"],
    },
    "Ag": {
        "tier": "T4", "phase": "fcc", "absorber": "Ag",
        "formula": "Ag", "stoich": {"Ag": 1},
        "cod_id": "9008459", "spg_nums": [225],
        "tags": ["metal", "fcc", "element_seed"],
    },
    "Ir": {
        "tier": "T4", "phase": "fcc", "absorber": "Ir",
        "formula": "Ir", "stoich": {"Ir": 1},
        "cod_id": "9008470", "spg_nums": [225],
        "tags": ["metal", "fcc", "element_seed"],
    },
    "Ce": {
        "tier": "T4", "phase": "fcc", "absorber": "Ce",
        "formula": "Ce", "stoich": {"Ce": 1},
        "cod_id": "9010987", "spg_nums": [225],
        "tags": ["metal", "fcc", "element_seed"],
    },
    "U": {
        "tier": "T4", "phase": "alpha-U", "absorber": "U",
        "formula": "U", "stoich": {"U": 1},
        "cod_id": "1527120", "spg_nums": [63],
        "tags": ["metal", "alpha-U", "element_seed"],
    },
    # ---- T4 续：补齐 Cr/Zn/Mo/Au，使 17 个吸收元素全部有 verified CIF ----
    "Cr": {
        "tier": "T4", "phase": "bcc", "absorber": "Cr",
        "formula": "Cr", "stoich": {"Cr": 1},
        "cod_id": "9008531", "spg_nums": [229],
        "tags": ["metal", "bcc", "element_seed"],
    },
    "Zn": {
        "tier": "T4", "phase": "hcp", "absorber": "Zn",
        "formula": "Zn", "stoich": {"Zn": 1},
        "cod_id": "9008522", "spg_nums": [194],
        "tags": ["metal", "hcp", "element_seed"],
    },
    "Mo": {
        "tier": "T4", "phase": "bcc", "absorber": "Mo",
        "formula": "Mo", "stoich": {"Mo": 1},
        "cod_id": "9008543", "spg_nums": [229],
        "tags": ["metal", "bcc", "element_seed"],
    },
    "Au": {
        "tier": "T4", "phase": "fcc", "absorber": "Au",
        "formula": "Au", "stoich": {"Au": 1},
        "cod_id": "9008463", "spg_nums": [225],
        "tags": ["metal", "fcc", "element_seed"],
    },
}

# 兼容旧报告/脚本命名。
TARGETS = VERIFIED_TARGETS


# ======================================================================
# 化学式
# ======================================================================
_ELEMENT_RE = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")


def parse_formula(formula: str) -> dict[str, float]:
    """把化学式解析为 ``{元素: 数量}``，支持括号与小数。

    例：``Fe2O3``、``Fe2 O3``、``Ni(OH)2``、``Al0.34 Ca0.66``。
    """
    text = str(formula or "").strip()
    if not text or text.lower() in ("?", "none", "n/a"):
        return {}
    out: dict[str, float] = {}

    def add_group(chunk: str, multiplier: float = 1.0) -> None:
        i, n = 0, len(chunk)
        while i < n:
            ch = chunk[i]
            if ch == "(":
                depth, j = 1, i + 1
                while j < n and depth:
                    if chunk[j] == "(":
                        depth += 1
                    elif chunk[j] == ")":
                        depth -= 1
                    j += 1
                inner = chunk[i + 1:j - 1]
                m = re.match(r"\d*\.?\d*", chunk[j:])
                num = m.group(0) if m else ""
                add_group(inner, multiplier * (float(num) if num else 1.0))
                i = j + len(num)
            elif ch.isupper():
                m = _ELEMENT_RE.match(chunk, i)
                if not m:
                    i += 1
                    continue
                el = m.group(1).capitalize()
                num = m.group(2)
                out[el] = out.get(el, 0.0) + (
                    float(num) if num else 1.0) * multiplier
                i = m.end()
            else:
                i += 1

    add_group(text)
    return out


def formula_matches(expected: str, actual: str,
                    tolerance: float = 0.05) -> tuple[bool, str]:
    """比较元素集合与归一化化学计量；返回 ``(ok, detail)``。"""
    tf = parse_formula(expected)
    cf = parse_formula(actual)
    if not tf:
        return False, f"目标化学式无法解析: {expected!r}"
    if not cf:
        return False, f"CIF 化学式无法解析: {actual!r}"
    if set(tf) != set(cf):
        return False, f"元素集合不同: 目标 {sorted(tf)} vs CIF {sorted(cf)}"
    t_sum = sum(tf.values())
    c_sum = sum(cf.values())
    if t_sum <= 0 or c_sum <= 0:
        return False, "化学式总数为 0"
    diffs = []
    for el in sorted(tf):
        target_ratio = tf[el] / t_sum
        actual_ratio = cf[el] / c_sum
        if abs(target_ratio - actual_ratio) > tolerance:
            diffs.append(
                f"{el}: 目标 {target_ratio:.4f} vs CIF {actual_ratio:.4f}")
    if diffs:
        return False, "化学计量不符: " + "; ".join(diffs)
    return True, "化学计量匹配"


# ======================================================================
# CIF 物理核验
# ======================================================================
def _cif_text(value) -> str:
    """去掉 CIF 词法层外引号，但保留 Hall 符号内部双引号。"""
    if value is None:
        return ""
    try:
        import gemmi
        return str(gemmi.cif.as_string(value)).strip()
    except Exception:
        return str(value).strip()


def get_spg_number(block) -> int | None:
    """从 CIF block 读取 IT 空间群号；优先 IT number，回退 H-M/Hall。"""
    try:
        import gemmi
    except ImportError:
        return None

    for tag in ("_space_group_IT_number",
                "_symmetry_Int_Tables_number"):
        raw = block.find_value(tag)
        if raw is not None:
            try:
                return int(float(_cif_text(raw)))
            except (TypeError, ValueError):
                pass

    for tag in ("_space_group_name_H-M_alt",
                "_symmetry_space_group_name_H-M",
                "_space_group_name_H-M_ref"):
        name = _cif_text(block.find_value(tag))
        if name:
            sg = gemmi.find_spacegroup_by_name(name)
            if sg is not None:
                return int(sg.number)

    for tag in ("_space_group_name_Hall",
                "_symmetry_space_group_name_Hall"):
        hall = _cif_text(block.find_value(tag))
        if hall:
            sg = gemmi.find_spacegroup_by_name(hall)
            if sg is not None:
                return int(sg.number)
            for candidate in gemmi.spacegroup_table():
                if str(candidate.hall).strip().lower() == hall.lower():
                    return int(candidate.number)
    return None


def verify_cif(path: str | Path, name: str,
               tolerance: float = 0.05) -> dict:
    """核验一个 CIF 是否可作为目标相的真实来源。

    返回 ``{ok, code, reason, sha256, spg, formula, parsed}``；
    不抛异常，便于批量报告。
    """
    target = VERIFIED_TARGETS[name]
    path = Path(path)
    result = {
        "ok": False, "code": "", "reason": "", "sha256": "",
        "spg": None, "formula": "", "path": str(path), "parsed": None,
    }
    if not path.is_file():
        result.update(code="cif_missing", reason=f"CIF 不存在: {path}")
        return result
    if is_generated_cif(path):
        result.update(code="generated_cif_not_allowed",
                      reason=f"{path.name}: 现场生成 CIF 不得进入 phase pack")
        return result

    try:
        parsed = phase_pack.parse_cif(path)
    except phase_pack.PhasePackError as exc:
        result.update(code=exc.code, reason=exc.detail)
        return result

    try:
        import gemmi
        doc = gemmi.cif.read(str(path))
        block = doc.sole_block() if len(doc) == 1 else doc[0]
    except Exception as exc:
        result.update(code="cif_parse_failed", reason=f"{path.name}: {exc}")
        return result

    formula = _cif_text(block.find_value("_chemical_formula_sum"))
    if not formula:
        result.update(code="missing_formula_sum",
                      reason=f"{path.name}: 缺 _chemical_formula_sum，"
                             f"无法做严格计量比核验")
        return result

    ok, detail = formula_matches(target["formula"], formula, tolerance)
    if not ok:
        result.update(code="composition_mismatch",
                      reason=f"{path.name}: {detail}")
        return result

    spg = get_spg_number(block)
    if spg is None:
        result.update(code="spg_unreadable",
                      reason=f"{path.name}: 无法读取 IT 空间群号")
        return result
    if spg not in target["spg_nums"]:
        result.update(code="spg_mismatch",
                      reason=f"{path.name}: 空间群 #{spg} 不在目标 "
                             f"{target['spg_nums']}")
        return result

    if target["absorber"] not in parsed["absorber_elements"]:
        result.update(code="absorber_missing",
                      reason=f"{path.name}: 缺目标吸收元素 "
                             f"{target['absorber']}")
        return result
    if not parsed["default_absorber_sites"]:
        result.update(code="absorber_site_missing",
                      reason=f"{path.name}: 目标吸收元素没有可解析位点")
        return result

    result.update(
        ok=True, code="ok",
        reason=f"spg=#{spg} {detail} sites="
               f"{','.join(parsed['default_absorber_sites'])}",
        sha256=sha256_file(path), spg=spg, formula=formula, parsed=parsed,
    )
    return result


# ======================================================================
# 缓存与下载
# ======================================================================
def _cache_candidates(cache_dir: Path, name: str, cod_id: str) -> list[Path]:
    return [
        cache_dir / f"{name}_COD_{cod_id}.cif",
        cache_dir / f"COD_{cod_id}.cif",
    ]


def _quarantine(cache_dir: Path, path: Path) -> Path:
    """把校验失败的缓存移到 ``_rejected/``，避免污染后续离线运行。"""
    rejected = cache_dir / "_rejected"
    rejected.mkdir(parents=True, exist_ok=True)
    digest = sha256_file(path)[:12] if path.is_file() else "missing"
    target = rejected / f"{path.name}.{digest}.rejected"
    try:
        os.replace(path, target)
    except OSError:
        pass
    return target


def _download_cif(cod_id: str, dest: Path, *,
                  timeout: float = 30.0, retries: int = 3) -> tuple[bool, str]:
    """把 COD 原始 CIF 下载到 ``dest``（原子写）；返回 ``(ok, url/reason)``。"""
    if not _COD_ID_RE.match(str(cod_id)):
        return False, f"非法 COD ID: {cod_id!r}"
    url = COD_DOWNLOAD_URL.format(cod_id=cod_id)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    context = ssl.create_default_context()
    headers = {"User-Agent": USER_AGENT, "Accept": "chemical/x-cif,text/plain,*/*"}

    last_error = ""
    for attempt in range(1, max(1, retries) + 1):
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(
                    request, timeout=timeout, context=context) as response:
                data = response.read()
            if not data or b"data_" not in data:
                last_error = "响应不是有效 CIF（缺 data_ 块）"
            else:
                tmp.write_bytes(data)
                os.replace(tmp, dest)
                return True, url
        except urllib.error.HTTPError as exc:
            last_error = f"HTTP {exc.code}"
            if exc.code in (400, 403, 404):
                break
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_error = str(exc)
        if attempt < retries:
            time.sleep(1.5 * attempt)
    try:
        tmp.unlink(missing_ok=True)
    except OSError:
        pass
    return False, f"{url}: {last_error}"


def acquire_target(name: str, *, cache_dir: Path, online: bool,
                   force: bool = False, timeout: float = 30.0,
                   retries: int = 3,
                   tolerance: float = 0.05) -> dict:
    """获取并核验目标 CIF；优先复用有效缓存，``online`` 才允许下载。"""
    target = VERIFIED_TARGETS[name]
    cod_id = target["cod_id"]
    cache_dir.mkdir(parents=True, exist_ok=True)

    existing = next(
        (p for p in _cache_candidates(cache_dir, name, cod_id) if p.is_file()),
        None)
    if existing is not None and not force:
        check = verify_cif(existing, name, tolerance)
        if check["ok"]:
            return {
                "ok": True, "status": "cached", "path": existing,
                "source_url": COD_DOWNLOAD_URL.format(cod_id=cod_id),
                "check": check, "rejected": None,
            }
        if not online:
            return {
                "ok": False, "status": "cache_invalid", "path": existing,
                "source_url": COD_DOWNLOAD_URL.format(cod_id=cod_id),
                "check": check, "rejected": None,
                "reason": f"缓存校验失败: {check['code']} — {check['reason']}",
            }
        rejected = _quarantine(cache_dir, existing)
    else:
        rejected = None

    if not online:
        return {
            "ok": False, "status": "offline_miss", "path": None,
            "source_url": COD_DOWNLOAD_URL.format(cod_id=cod_id),
            "check": None, "rejected": rejected,
            "reason": f"离线且没有有效缓存: {name} COD:{cod_id}",
        }

    final = cache_dir / f"{name}_COD_{cod_id}.cif"
    downloaded, detail = _download_cif(
        cod_id, final, timeout=timeout, retries=retries)
    if not downloaded:
        return {
            "ok": False, "status": "download_failed", "path": None,
            "source_url": COD_DOWNLOAD_URL.format(cod_id=cod_id),
            "check": None, "rejected": rejected, "reason": detail,
        }

    check = verify_cif(final, name, tolerance)
    if not check["ok"]:
        rejected = _quarantine(cache_dir, final)
        return {
            "ok": False, "status": "verify_failed", "path": final,
            "source_url": detail, "check": check, "rejected": rejected,
            "reason": f"{check['code']} — {check['reason']}",
        }
    return {
        "ok": True, "status": "downloaded", "path": final,
        "source_url": detail, "check": check, "rejected": rejected,
    }


# ======================================================================
# phase pack 导入
# ======================================================================
def _rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except (OSError, ValueError):
        return path.resolve().as_posix()


def _existing_entry(pack_root: Path, cod_id: str) -> dict | None:
    manifest_path = pack_root / phase_pack.MANIFEST_NAME
    if not manifest_path.is_file():
        return None
    try:
        manifest = phase_pack.load_manifest(pack_root)
    except phase_pack.PhasePackError:
        return None
    for entry in manifest.get("entries", []):
        if (str(entry.get("source", "")).lower() == "cod"
                and str(entry.get("id")) == str(cod_id)):
            return entry
    return None


def import_target(name: str, cif_path: Path, *, pack_root: Path,
                  source_url: str, replace: bool = False) -> dict:
    """把已核验 CIF 导入 phase pack；不写 motif。"""
    target = VERIFIED_TARGETS[name]
    meta = {
        "source": "COD",
        "id": target["cod_id"],
        "formula": target["formula"],
        "experimental": True,
        "model_level": "crystal_cif",
        # 目标名（如 Fe2O3）写入 tags，供 phase_coverage 的常见相
        # 匹配使用；不参与自动选相。
        "tags": list(dict.fromkeys([name] + list(target["tags"]))),
        "notes": (
            f"{target['tier']} {target['phase']}; "
            f"COD:{target['cod_id']}; source_url={source_url}"
        ),
    }
    return phase_pack.import_cif(
        cif_path, pack_root, meta=meta, replace=replace,
        added_by="tools/cif_crawler.py")


# ======================================================================
# CLI
# ======================================================================
def _select_targets(args) -> list[str]:
    selected: list[str] = []
    if args.all:
        selected.extend(VERIFIED_TARGETS)
    if args.t3:
        selected.extend(n for n, t in VERIFIED_TARGETS.items()
                        if t["tier"] == "T3")
    if args.t4:
        selected.extend(n for n, t in VERIFIED_TARGETS.items()
                        if t["tier"] == "T4")
    if args.target:
        lowered = {n.lower(): n for n in VERIFIED_TARGETS}
        for raw in args.target:
            name = lowered.get(str(raw).strip().lower())
            if name is None:
                raise SystemExit(
                    f"未知目标 {raw!r}；可用目标: "
                    f"{', '.join(sorted(VERIFIED_TARGETS))}")
            selected.append(name)
    # 去重并保持目标表顺序，报告稳定。
    order = {name: i for i, name in enumerate(VERIFIED_TARGETS)}
    return sorted(dict.fromkeys(selected), key=lambda n: order[n])


def _print_list() -> None:
    print(f"{'NAME':10s} {'TIER':5s} {'ABS':4s} {'COD_ID':9s} "
          f"{'SPG':5s} PHASE")
    for name, target in VERIFIED_TARGETS.items():
        spg = ",".join(f"#{n}" for n in target["spg_nums"])
        print(f"{name:10s} {target['tier']:5s} {target['absorber']:4s} "
              f"{target['cod_id']:9s} {spg:5s} {target['phase']}")


def _write_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(data, ensure_ascii=False, indent=2)
               + "\n").encode("utf-8")
    if path.is_file() and path.read_bytes() == payload:
        return
    # 统一 LF，保证构建产物跨平台字节一致。
    with path.open("wb") as fh:
        fh.write(payload)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="构建期 COD CIF 抓取与 phase pack 导入（运行期不联网）")
    parser.add_argument("--target", action="append", default=None,
                        help="目标名，可重复，例如 --target Fe2O3")
    parser.add_argument("--t3", action="store_true", help="全部 T3 目标")
    parser.add_argument("--t4", action="store_true", help="全部 T4 目标")
    parser.add_argument("--all", action="store_true", help="全部目标")
    parser.add_argument("--list", action="store_true", help="列出目标后退出")
    parser.add_argument("--online", action="store_true",
                        help="允许联网下载；默认只用本地缓存（离线）")
    parser.add_argument("--offline", action="store_true",
                        help="显式声明离线（默认行为；与 --online 互斥）")
    parser.add_argument("--cache", default=str(DEFAULT_CACHE),
                        help="原始 CIF 缓存目录（默认 data/cif_raw）")
    parser.add_argument("--pack", default=str(phase_pack.PACK_ROOT),
                        help="phase pack 根目录")
    parser.add_argument("--report", default=str(DEFAULT_REPORT),
                        help="构建期抓取报告 JSON 路径")
    parser.add_argument("--index", default=None,
                        help="phase_index.json 输出路径（默认主 pack 才写）")
    parser.add_argument("--coverage", default=None,
                        help="phase_coverage.json 输出路径（默认主 pack 才写）")
    parser.add_argument("--no-index", action="store_true",
                        help="不重建 phase_index/phase_coverage（测试用）")
    parser.add_argument("--force", action="store_true",
                        help="忽略缓存，强制重新下载（需 --online）")
    parser.add_argument("--replace", action="store_true",
                        help="同 COD ID 不同 sha256 时显式覆盖 phase pack")
    parser.add_argument("--dry-run", action="store_true",
                        help="只获取/核验，不写 phase pack")
    parser.add_argument("--json", action="store_true",
                        help="把报告 JSON 输出到 stdout")
    parser.add_argument("--timeout", type=float, default=30.0,
                        help="单次下载超时秒数（默认 30）")
    parser.add_argument("--retries", type=int, default=3,
                        help="下载重试次数（默认 3）")
    parser.add_argument("--composition-tolerance", type=float, default=0.05,
                        help="化学计量比容差（默认 0.05 = 5%%）")
    args = parser.parse_args(argv)

    if args.list:
        _print_list()
        return 0
    if args.online and args.offline:
        parser.error("--online 与 --offline 不能同时使用")
    if args.force and not args.online:
        parser.error("--force 需要 --online")

    selected = _select_targets(args)
    if not selected:
        parser.error("请指定 --target/--t3/--t4/--all，或使用 --list")

    cache_dir = Path(args.cache)
    pack_root = Path(args.pack)
    report_path = Path(args.report)
    is_default_pack = pack_root.resolve() == phase_pack.PACK_ROOT.resolve()
    index_path = (Path(args.index) if args.index
                  else (phase_pack.INDEX_PATH if is_default_pack else None))
    coverage_path = (Path(args.coverage) if args.coverage
                     else (phase_pack.COVERAGE_PATH
                           if is_default_pack else None))

    results: dict[str, dict] = {}
    failures = 0
    for name in selected:
        target = VERIFIED_TARGETS[name]
        acquired = acquire_target(
            name, cache_dir=cache_dir, online=args.online,
            force=args.force, timeout=args.timeout, retries=args.retries,
            tolerance=args.composition_tolerance)
        row = {
            "tier": target["tier"],
            "absorber": target["absorber"],
            "phase": target["phase"],
            "cod_id": target["cod_id"],
            "spg_nums": list(target["spg_nums"]),
            "status": acquired["status"],
            "source_url": acquired.get("source_url"),
            "cache_file": (_rel(acquired["path"])
                           if acquired.get("path") else None),
            "rejected_file": (_rel(acquired["rejected"])
                              if acquired.get("rejected") else None),
            "reason": acquired.get("reason", ""),
        }
        if acquired.get("check"):
            check = acquired["check"]
            row.update({
                "spg": check.get("spg"),
                "formula": check.get("formula"),
                "cif_sha256": check.get("sha256"),
                "absorber_sites": (check.get("parsed") or {}).get(
                    "default_absorber_sites", []),
            })

        if acquired["ok"] and not args.dry_run:
            existing = _existing_entry(pack_root, target["cod_id"])
            try:
                entry = import_target(
                    name, acquired["path"], pack_root=pack_root,
                    source_url=acquired["source_url"],
                    replace=args.replace)
            except phase_pack.PhasePackError as exc:
                row.update(status="import_failed",
                           reason=f"{exc.code}: {exc.detail}")
                failures += 1
            else:
                same = (existing is not None
                        and existing.get("cif_sha256") == entry["cif_sha256"])
                row.update(
                    status="already_present" if same else "imported",
                    pack_entry={
                        "source": entry["source"],
                        "id": entry["id"],
                        "cif_path": entry["cif_path"],
                        "cif_sha256": entry["cif_sha256"],
                    })
        elif acquired["ok"] and args.dry_run:
            row["status"] = "verified"
        else:
            failures += 1

        results[name] = row
        marker = "ok" if row["status"] in (
            "imported", "already_present", "verified", "cached") else "FAIL"
        print(f"[{marker}] {name:10s} COD:{target['cod_id']} "
              f"spg={row.get('spg')} status={row['status']} "
              f"sha={(row.get('cif_sha256') or '')[:12]}")
        if row.get("reason"):
            print(f"       {row['reason']}")

    if (not args.dry_run and not args.no_index
            and index_path is not None and coverage_path is not None
            and any(r["status"] in ("imported", "already_present")
                    for r in results.values())):
        try:
            index = phase_pack.build_index(pack_root, out_path=index_path)
            coverage = phase_pack.build_coverage(
                pack_root, index=index, out_path=coverage_path)
            print(f"[index] entries={index['summary']['entries']} "
                  f"verified={index['summary']['verified_cif']} "
                  f"elements={len(index['summary']['elements_covered'])}/"
                  f"{index['summary']['elements_total']}")
            print("[coverage] verified=" + ",".join(
                coverage["summary"]["verified_elements"]))
        except phase_pack.PhasePackError as exc:
            print(f"[index] FAILED: {exc}", file=sys.stderr)
            failures += 1

    report = {
        "schema_version": 1,
        "generated_by": "tools/cif_crawler.py",
        "online": bool(args.online),
        "dry_run": bool(args.dry_run),
        "pack_root": _rel(pack_root),
        "cache_root": _rel(cache_dir),
        "targets": results,
        "summary": {
            "selected": len(selected),
            "imported": sum(1 for r in results.values()
                            if r["status"] == "imported"),
            "already_present": sum(
                1 for r in results.values()
                if r["status"] == "already_present"),
            "failed": sum(1 for r in results.values()
                          if r["status"] not in (
                              "imported", "already_present", "verified",
                              "cached")),
        },
    }
    if not args.dry_run or args.report:
        _write_json(report_path, report)
    if args.json:
        print(json.dumps(report, ensure_ascii=False, indent=2))

    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
    raise SystemExit(main())
