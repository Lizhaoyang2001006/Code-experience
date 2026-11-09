# P6.1 执行指令：单晶位/多晶位正式 CN 端到端闭环修复（5 项）

> 本文件针对 `eb01675`（P6 完成态）做全链路穿透式体检后暴露出的 **5 个端到端问题**，提供可直接执行的精确替换方案与检验命令。
> 执行顺序：**先按 Step 1–5 直接替换代码，再运行 Step 6 检验**。

---

## 0. 硬性工程守则（最高优先级）

1. **只允许修改** `C:\Users\lizha\Desktop\Code-experience-cn`；原 `C:\Users\lizha\Desktop\Code-experience` 与 `C:\Users\lizha\Desktop\CN_CIF_执行方案` 严禁触碰。
2. **禁止破坏性 Git 命令**：严禁 `git checkout -- .`、`git reset --hard`、`git clean -fd`。
3. **严禁放宽 S₀² 标定门禁**：
   - `absorber.py` 静态表中 12 个扩展元素保持 `s02=None, calibrated=False`（守住 `#60`）；
   - `calibrate_s02.py` 门禁常量保持不变。
4. **禁止提交第三方原始谱与运行日志**：`_standards_lib/`、`data/standards_raw/`、`logs/` 严禁 `git add`。
5. **固定 Python 解释器**：
   `C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe`

---

## 1. 全链路体检发现的 5 个问题

| 编号 | 严重度 | 位置 | 现象与根因 |
|---|---|---|---|
| **问题 1** | **P0 阻塞** | `xafs_fit.py` `fit_identified()`（第 1138–1237、1733–1756 行） | **单晶位 `--formal-cn` 被诊断库 `_structure_lib` 劫持**：多位点 `run_multi_site_formal_cn`（第 2883 行）已用 `_site_motif_from_model(model)` 从已验证 CIF 构建 FEFF 路径；但单晶位 `fit_identified` 仍走 `_structure_lib` 识别分支。当用户对 16 个单晶位 P5-C 氧化物（`RuO2`、`Rh2O3`、`PdO`、`IrO2`、`PtO2`、`Au2O3`、`Ag2O`、`Cr2O3`、`MnO`、`MnO2`、`ZnO`、`MoO2`、`MoO3`、`UO2`、`FeO`、`FeOOH`）或自定义 `--phase-cif` 运行 `--fit --formal-cn --phase-id COD:...` 时，因 `_structure_lib/<El>` 中无该氧化物 motif（或 S6.5 返回 `partial`，如 `Ce-L3-CeO2.xdi` 配 `COD:1562989`），`fit_identified` 在第 1143–1191 行直接退回 `l1_first_shell` 并强制覆写 `formal_cn.status = "l1_diagnostic_only"`，根本不执行 P3 条件 CN；即使加 `--allow-ambiguous` 也会误取金属单质 motif（如 `ru_hcp`）导致 `map_n0_to_plan` 元素不匹配失败。 |
| **问题 2** | **P1 状态不一致** | `run_xafs.py`（第 636 行） | **多位点 `--absorber-site all` 未同步顶层 `cn_status`**：`multi_site_formal` 分支只写了 `result.fit["formal_cn"] = multi`，未写 `result.fit["cn_status"] = multi.get("status", "multi_site_failed")`，导致 `_meta.json` 中 `fit.cn_status` 仍残留 `"diagnostic_only"`。 |
| **问题 3** | **P1 溯源失真** | `phase_resolver.py`（第 377 行）、`run_xafs.py`（第 565–586 行） | **`StructureModel.occupancy_policy` 硬编码 `"explicit"`**：`_candidate_to_structure_model` 硬编码 `occupancy_policy="explicit"`，且 `run_xafs.py` 在 `cif_gate` 通过后未把 `gate_result.occupancy_source`（如 `"implicit_1.0"`）与 `n0_shells`/`warnings` 回填给 `structure_model`，导致 11 个无 `_atom_site_occupancy` 的有序氧化物与 12 个单质 CIF 在 `formal_cn.model.occupancy_policy` 中被误标为 `"explicit"`。 |
| **问题 4** | **P2 提示缺失** | `xafs_fit.py`（第 1086、1889、2847 行） | **跳过 S6.5 时缺 `--s02-from-t5` 提示**：`structure_id.py` 在 P6 已加上 `--s02-from-t5` 提示，但当用户用 `--catalytic-model`（跳过 `structure_id.identify_result`）或直接调用 `fit_identified` / `run_operando_exafs_series` / `run_multi_site_formal_cn` 未传 S₀² 时，报错未提示 `--s02-from-t5`。 |
| **问题 5** | **P2 入口对齐** | `xafs_report.py`（第 1424、1820 行） | **`xafs_report.py` 缺 `--s02-from-t5` 与 `--phase-*` 透传**：`xafs_report.py` 虽提供 `--formal-cn`，但未接入 `run_xafs.load_t5_s02` 与 `phase_resolver`/`cif_gate`，无法直接从报告 CLI 跑 T5 标定元素或 `--phase-id`/`--phase-cif` 正式 CN。 |

---

## Step 1：修复 `xafs_fit.py` 单晶位正式 CN 执行链与 `--s02-from-t5` 提示（问题 1 & 问题 4）

文件：`C:\Users\lizha\Desktop\Code-experience-cn\xafs_fit.py`

### 1.1 新增 `_t5_s02_hint` 与 `_attach_single_site_formal_cn` 辅助函数

在 `def _formal_cn_meta(`（第 993 行）定义结束之后、`def fit_identified(`（第 1035 行）**之前**插入：

**原文**：
```python
    return {
        "diagnostic_only": True,
        "cn_status": status,
        "formal_cn": {
            "status": status,
            "rows": [],
            "reject_reason": reason,
            "model": (structure_model.to_json()
                      if structure_model is not None else None),
        },
    }


def fit_identified(res, motif_id: str | None = None,
```

**替换为**：
```python
    return {
        "diagnostic_only": True,
        "cn_status": status,
        "formal_cn": {
            "status": status,
            "rows": [],
            "reject_reason": reason,
            "model": (structure_model.to_json()
                      if structure_model is not None else None),
        },
    }


def _t5_s02_hint(symbol: str) -> str:
    """若该元素已有 T5 独立标定产物，返回 CLI 提示后缀。"""
    t5_cal = (Path(__file__).resolve().parent / "data"
              / "s02_calibration" / symbol / "calibration.json")
    return ("（该元素已有 T5 独立标定值，CLI 可加 --s02-from-t5 或显式传 --s02）"
            if t5_cal.is_file() else "")


def _attach_single_site_formal_cn(
    out: dict,
    res,
    *,
    absorber,
    s02_label: str,
    formal_cn: bool,
    structure_model: StructureModel | None,
    gate_result=None,
    rmin: float | None = None,
    rmax: float | None = None,
    include_ms: bool = False,
    isotropic_deltar: bool = False,
    hetero_cocupancy: bool = False,
    is_l1: bool = False,
    verbose: bool = True,
) -> None:
    """挂载单晶位正式 CN 结果（与诊断轨 _structure_lib 解耦）。

    - 未提供 ``structure_model`` 且诊断轨为 L1 时，保持 #76 的 ``l1_diagnostic_only``；
    - 未提供 ``structure_model`` 且诊断轨为 motif 时，保持 #75 的 ``no_cif``；
    - 提供合格 ``structure_model`` + ``gate_result.ok`` 时，直接由已验证 CIF 构建
      ``_site_motif_from_model(structure_model)`` -> ``select_shell_paths`` ->
      ``map_n0_to_plan`` -> ``run_conditional_cn`` -> ``build_formal_cn_block``。
    """
    if is_l1 and structure_model is None:
        out.update(_formal_cn_meta(
            formal_cn, None,
            status_override=("l1_diagnostic_only" if formal_cn else None),
            reason_override=("L1 合成团簇只做诊断，不进入正式 CN 表"
                             if formal_cn else None)))
        return

    out.update(_formal_cn_meta(formal_cn, structure_model))
    if not (formal_cn and structure_model is not None):
        return
    formal_block = out.get("formal_cn") or {}
    if formal_block.get("status") != "pending_consistency":
        return

    conditional_cn = None
    conditional_cn_error = ""
    n0_mapping = None
    formal_plan = None
    if gate_result is not None and getattr(gate_result, "ok", False):
        try:
            formal_motif = _site_motif_from_model(structure_model)
            formal_plan = select_shell_paths(
                formal_motif, edge=absorber.edge,
                rmin=rmin, rmax=rmax,
                kmin=res.ft.kmin, kmax=res.ft.kmax,
                include_ms=include_ms,
                isotropic_deltar=isotropic_deltar,
                hetero_cocupancy=hetero_cocupancy)
            n0_mapping = map_n0_to_plan(
                formal_plan, gate_result.n0_shells,
                kmin=res.ft.kmin, kmax=res.ft.kmax)
            if n0_mapping["mapping_ok"]:
                conditional_cn = run_conditional_cn(
                    res, formal_plan, n0_mapping,
                    absorber_symbol=absorber.symbol, s02=absorber.s02,
                    kmin=res.ft.kmin, kmax=res.ft.kmax,
                    rmin=formal_plan["rmin"], rmax=formal_plan["rmax"],
                    verbose=verbose)
            else:
                conditional_cn_error = (
                    "N₀ 映射失败: " + n0_mapping["mapping_reason"])
            if verbose and conditional_cn is not None:
                print(f"  [P3] 条件 CN: N/N₀={conditional_cn['n_ratio']} "
                      f"({conditional_cn['n_status']}), "
                      f"consistency_passed={conditional_cn['consistency_passed']}")
        except (CNError, XAFSConfigError, XAFSDataError) as exc:
            conditional_cn_error = (
                "条件 CN 跳过: " + getattr(exc, "message", str(exc)))
            if verbose:
                print(f"  [P3] {conditional_cn_error}")
    else:
        conditional_cn_error = "缺少 CIF gate 结果，未运行条件 CN"
        if verbose:
            print(f"  [P3] {conditional_cn_error}")

    if conditional_cn is not None:
        out["conditional_cn"] = conditional_cn
    if (conditional_cn is not None and n0_mapping is not None
            and formal_plan is not None):
        block = build_formal_cn_block(
            conditional_cn=conditional_cn, plan=formal_plan,
            n0_mapping=n0_mapping,
            structure_model=structure_model,
            absorber=absorber,
            cn_systematics=out.get("cn_systematics"),
            s02_source=s02_label,
            data_sha256=(getattr(res, "data_sha256", "") or ""))
    else:
        block = {
            "status": "consistency_failed",
            "rows": [],
            "reject_reason": (conditional_cn_error or "条件 CN 未运行"),
            "model": structure_model.to_json(),
        }
    out["formal_cn"] = block
    out["cn_status"] = block["status"]


def fit_identified(res, motif_id: str | None = None,
```

### 1.2 在 `fit_identified` 中使用 `_t5_s02_hint` 与 `_attach_single_site_formal_cn`

**原文 1（第 1086–1090 行）**：
```python
    if absorber.s02 is None:
        raise XAFSConfigError(
            "FIT", f"{absorber.symbol} 未标定 S₀²，拒绝 CN 定量；"
                   f"先在 absorber.py 用标准样标定")
```

**替换为**：
```python
    if absorber.s02 is None:
        raise XAFSConfigError(
            "FIT", f"{absorber.symbol} 未标定 S₀²，拒绝 CN 定量；"
                   f"先在 absorber.py 用标准样标定{_t5_s02_hint(absorber.symbol)}")
```

**原文 2（第 1186–1238 行）**：
```python
            out.update(_formal_cn_meta(
                formal_cn, structure_model,
                status_override=("l1_diagnostic_only" if formal_cn else None),
                reason_override=("L1 合成团簇只做诊断，不进入正式 CN 表"
                                 if formal_cn else None)))
            return out

    motif = (catalytic_motif if catalytic_motif is not None
             else _find_motif(absorber.symbol, motif_id))
    if fit_weights and motif.kind != "mixture":
        raise XAFSConfigError(
            "FIT", f"{motif.id} 不是混合相（kind={motif.kind}），"
                   f"--fit-weights 只用于多位点吸收的 mixture")
    plan = select_shell_paths(motif, edge=absorber.edge,
                              rmin=rmin, rmax=rmax,
                              kmin=res.ft.kmin, kmax=res.ft.kmax,
                              group_by_component=fit_weights,
                              include_ms=include_ms,
                              isotropic_deltar=isotropic_deltar,
                              hetero_cocupancy=hetero_cocupancy)

    # ---- P3 条件 CN：固定 N → 自由 N（首壳单一比例因子）→ Hamilton F-test ----
    conditional_cn = None
    conditional_cn_error = ""
    n0_mapping = None
    if formal_cn and gate_result is not None and gate_result.ok:
        try:
            n0_mapping = map_n0_to_plan(plan, gate_result.n0_shells,
                                        kmin=res.ft.kmin, kmax=res.ft.kmax)
            if n0_mapping["mapping_ok"]:
                conditional_cn = run_conditional_cn(
                    res, plan, n0_mapping,
                    absorber_symbol=absorber.symbol, s02=absorber.s02,
                    kmin=res.ft.kmin, kmax=res.ft.kmax,
                    rmin=plan["rmin"], rmax=plan["rmax"],
                    verbose=verbose)
            else:
                conditional_cn_error = (
                    "N₀ 映射失败: " + n0_mapping["mapping_reason"])
            if verbose and conditional_cn is not None:
                print(f"  [P3] 条件 CN: N/N₀={conditional_cn['n_ratio']} "
                      f"({conditional_cn['n_status']}), "
                      f"consistency_passed={conditional_cn['consistency_passed']}")
        except (XAFSConfigError, XAFSDataError) as exc:
            conditional_cn_error = (
                "条件 CN 跳过: " + getattr(exc, "message", str(exc)))
            if verbose:
                print(f"  [P3] {conditional_cn_error}")
    elif formal_cn and structure_model is not None:
        conditional_cn_error = "缺少 CIF gate 结果，未运行条件 CN"
        if verbose:
            print(f"  [P3] {conditional_cn_error}")
```

**替换为**：
```python
            _attach_single_site_formal_cn(
                out, res, absorber=absorber, s02_label=s02_label,
                formal_cn=formal_cn, structure_model=structure_model,
                gate_result=gate_result, rmin=rmin, rmax=rmax,
                include_ms=include_ms, isotropic_deltar=isotropic_deltar,
                hetero_cocupancy=hetero_cocupancy, is_l1=True,
                verbose=verbose)
            return out

    motif = (catalytic_motif if catalytic_motif is not None
             else _find_motif(absorber.symbol, motif_id))
    if fit_weights and motif.kind != "mixture":
        raise XAFSConfigError(
            "FIT", f"{motif.id} 不是混合相（kind={motif.kind}），"
                   f"--fit-weights 只用于多位点吸收的 mixture")
    plan = select_shell_paths(motif, edge=absorber.edge,
                              rmin=rmin, rmax=rmax,
                              kmin=res.ft.kmin, kmax=res.ft.kmax,
                              group_by_component=fit_weights,
                              include_ms=include_ms,
                              isotropic_deltar=isotropic_deltar,
                              hetero_cocupancy=hetero_cocupancy)
```

**原文 3（第 1608–1615 行）**：
```python
           "cn_systematics": systematic_budget(
               fit, flags=flags, self_abs=self_abs,
               truncated_shells=len(plan["truncated"]),
               fit_weights=fit_weights)}
    if conditional_cn is not None:
        out["conditional_cn"] = conditional_cn
    if fit_weights and fit.mixture_weights:
```

**替换为**：
```python
           "cn_systematics": systematic_budget(
               fit, flags=flags, self_abs=self_abs,
               truncated_shells=len(plan["truncated"]),
               fit_weights=fit_weights)}
    if fit_weights and fit.mixture_weights:
```

**原文 4（第 1733–1767 行）**：
```python
    out.update(_formal_cn_meta(formal_cn, structure_model))
    if formal_cn and structure_model is not None:
        formal_block = out.get("formal_cn") or {}
        if formal_block.get("status") == "pending_consistency":
            if conditional_cn is not None and n0_mapping is not None:
                block = build_formal_cn_block(
                    conditional_cn=conditional_cn, plan=plan,
                    n0_mapping=n0_mapping,
                    structure_model=structure_model,
                    absorber=absorber,
                    cn_systematics=out.get("cn_systematics"),
                    s02_source=s02_label,
                    data_sha256=(getattr(res, "data_sha256", "") or ""))
            else:
                block = {
                    "status": "consistency_failed",
                    "rows": [],
                    "reject_reason": (conditional_cn_error
                                      or "条件 CN 未运行"),
                    "model": structure_model.to_json(),
                }
            out["formal_cn"] = block
            out["cn_status"] = block["status"]
    if formal_cn and catalytic_model is not None:
        out["formal_cn"] = {
            "status": "hypothesis_model_only",
            "rows": [],
            "reject_reason": (
                "catalytic_model 属于无实验 CIF 的催化位点假设模型，"
                "仅输出诊断 CN，不产出正式 CN"
            ),
            "model": None,
        }
        out["cn_status"] = "hypothesis_model_only"
    return out
```

**替换为**：
```python
    if formal_cn and catalytic_model is not None:
        out["diagnostic_only"] = True
        out["formal_cn"] = {
            "status": "hypothesis_model_only",
            "rows": [],
            "reject_reason": (
                "catalytic_model 属于无实验 CIF 的催化位点假设模型，"
                "仅输出诊断 CN，不产出正式 CN"
            ),
            "model": None,
        }
        out["cn_status"] = "hypothesis_model_only"
    else:
        _attach_single_site_formal_cn(
            out, res, absorber=absorber, s02_label=s02_label,
            formal_cn=formal_cn, structure_model=structure_model,
            gate_result=gate_result, rmin=rmin, rmax=rmax,
            include_ms=include_ms, isotropic_deltar=isotropic_deltar,
            hetero_cocupancy=hetero_cocupancy, is_l1=False,
            verbose=verbose)
    return out
```

### 1.3 加固 `_site_fractional`（支持无 `asym_sites` 的手动 `StructureModel` 回退读取 `source_cif`）及补齐其余两处 `--s02-from-t5` 提示

**原文 5（`run_operando_exafs_series` 第 1889–1894 行）**：
```python
    if absorber.s02 is None:
        raise XAFSConfigError(
            "OPERANDO",
            f"{absorber.symbol} 未标定 S₀²；原位 EXAFS 序列联合精修需要 "
            f"显式 s02（0.4–1.3）")
```

**替换为**：
```python
    if absorber.s02 is None:
        raise XAFSConfigError(
            "OPERANDO",
            f"{absorber.symbol} 未标定 S₀²；原位 EXAFS 序列联合精修需要 "
            f"显式 s02（0.4–1.3）{_t5_s02_hint(absorber.symbol)}")
```

**原文 6（`_site_fractional` 第 2727–2739 行）**：
```python
def _site_fractional(model: StructureModel) -> list[float]:
    """从 StructureModel 取吸收位的分数坐标。"""
    label = str(model.absorber_site or "")
    m = re.match(r"frac:\(([^,]+),([^,]+),([^)]+)\)", label)
    if m:
        return [float(m.group(1)), float(m.group(2)), float(m.group(3))]
    for site in (model.asym_sites or []):
        if str(site.get("label")) == label:
            return [float(site["x"]), float(site["y"]), float(site["z"])]
    raise CNError(
        "absorber_site_missing",
        f"{model.model_id}: 找不到吸收位 {label!r} 的分数坐标")
```

**替换为**：
```python
def _site_fractional(model: StructureModel) -> list[float]:
    """从 StructureModel 取吸收位的分数坐标（asym_sites 缺失时回退解析 source_cif）。"""
    label = str(model.absorber_site or "")
    m = re.match(r"frac:\(([^,]+),([^,]+),([^)]+)\)", label)
    if m:
        return [float(m.group(1)), float(m.group(2)), float(m.group(3))]
    for site in (model.asym_sites or []):
        if str(site.get("label")) == label:
            return [float(site["x"]), float(site["y"]), float(site["z"])]
    if model.source_cif and Path(model.source_cif).is_file():
        import phase_pack
        parsed = phase_pack.parse_cif(model.source_cif)
        for site in parsed.get("sites", []):
            if str(site.get("label")) == label:
                fx, fy, fz = site["fract"]
                return [float(fx), float(fy), float(fz)]
    raise CNError(
        "absorber_site_missing",
        f"{model.model_id}: 找不到吸收位 {label!r} 的分数坐标")
```

**原文 7（`run_multi_site_formal_cn` 第 2847–2850 行）**：
```python
    if absorber.s02 is None:
        raise XAFSConfigError(
            "FIT", f"{absorber.symbol} 未标定 S₀²，拒绝多位点 CN 定量")
```

**替换为**：
```python
    if absorber.s02 is None:
        raise XAFSConfigError(
            "FIT", f"{absorber.symbol} 未标定 S₀²，拒绝多位点 CN 定量"
                   f"{_t5_s02_hint(absorber.symbol)}")
```

---

## Step 2：修复 `phase_resolver.py` 的 `occupancy_policy` 硬编码（问题 3）

文件：`C:\Users\lizha\Desktop\Code-experience-cn\phase_resolver.py`

**原文（第 364–387 行）**：
```python
    model = StructureModel(
        model_id=f"{candidate.source}:{candidate.id}",
        absorber=absorber,
        absorber_site=site,
        model_level=candidate.model_level,
        source_type=source_type,
        source_ref=source_ref,
        source_cif=candidate.cif_path,
        source_cif_sha256=candidate.cif_sha256,
        cell=candidate.cell,
        spacegroup=candidate.spacegroup,
        asym_sites=asym_sites or None,
        cif_path=candidate.cif_path,
        occupancy_policy="explicit",
        edits=[],
        assumptions={
            "p1b": "phase_resolver resolved; verified_source_cif; "
                   "formal CN gate applies",
        },
        n0_shells=[],
        warnings=[],
        provenance_status="verified_source_cif",
    )
    return model
```

**替换为**：
```python
    has_occ_tag = False
    if candidate.cif_path and Path(candidate.cif_path).is_file():
        has_occ_tag = "_atom_site_occupancy" in Path(
            candidate.cif_path).read_text(encoding="utf-8", errors="ignore")
    occ_policy = "explicit" if has_occ_tag else "implicit_1.0"

    model = StructureModel(
        model_id=f"{candidate.source}:{candidate.id}",
        absorber=absorber,
        absorber_site=site,
        model_level=candidate.model_level,
        source_type=source_type,
        source_ref=source_ref,
        source_cif=candidate.cif_path,
        source_cif_sha256=candidate.cif_sha256,
        cell=candidate.cell,
        spacegroup=candidate.spacegroup,
        asym_sites=asym_sites or None,
        cif_path=candidate.cif_path,
        occupancy_policy=occ_policy,
        edits=[],
        assumptions={
            "p1b": "phase_resolver resolved; verified_source_cif; "
                   "formal CN gate applies",
        },
        n0_shells=[],
        warnings=[],
        provenance_status="verified_source_cif",
    )
    return model
```

---

## Step 3：修复 `run_xafs.py` 门禁属性回填与多位点 `cn_status` 同步（问题 2 & 问题 3）

文件：`C:\Users\lizha\Desktop\Code-experience-cn\run_xafs.py`

### 3.1 将 `cif_gate` 结果同步回 `structure_model` / `structure_models`

**原文（第 565–586 行）**：
```python
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
        else:
            gate_result = cif_gate.gate_cif(
                resolution.selected.cif_path, absorber,
                absorber_site=resolution.absorber_site,
                **gate_kwargs)
            if not gate_result.ok:
                raise XAFSConfigError(
                    "FIT", f"CIF gate rejected: {gate_result.code} — {gate_result.reason}")
```

**替换为**：
```python
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
```

### 3.2 多位点正式 CN 完成后同步 `result.fit["cn_status"]`

**原文（第 627–637 行）**：
```python
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
```

**替换为**：
```python
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
```

---

## Step 4：补齐 `xafs_report.py` 的 `--s02-from-t5` 与 `--phase-*` 正式 CN 入口（问题 5）

文件：`C:\Users\lizha\Desktop\Code-experience-cn\xafs_report.py`

### 4.1 增加 CLI 参数（`--s02-from-t5`、`--phase-cif`、`--phase-id`、`--phase`、`--phase-index`、`--phase-pack`、`--absorber-site`、`--offline`）

**原文（第 1424–1427 行）**：
```python
    parser.add_argument("--s02", type=float, default=None,
                        help="B18-2 运行时覆盖 S₀²（0.4–1.3），解锁未标定元素；"
                             "不修改 absorber.py 静态标定表")
```

**替换为**：
```python
    parser.add_argument("--s02", type=float, default=None,
                        help="B18-2 运行时覆盖 S₀²（0.4–1.3），解锁未标定元素；"
                             "不修改 absorber.py 静态标定表")
    parser.add_argument("--s02-from-t5", action="store_true",
                        help="P5 从 data/s02_calibration/<El>/calibration.json "
                             "加载通过冻结门禁的独立 S₀²（显式 --s02 优先）")
    parser.add_argument("--phase-cif", action="append", default=None,
                        help="P1b 用户真实 CIF 路径（可重复，最高优先级）")
    parser.add_argument("--phase-id", default=None,
                        help="P1b 精确数据库编号 SOURCE:ID，如 COD:1100108")
    parser.add_argument("--phase", default=None,
                        help="P1b 相 SPEC 字符串")
    parser.add_argument("--phase-index", default=None,
                        help="P1b 本地 phase_index.json 路径")
    parser.add_argument("--phase-pack", default=None,
                        help="P1b phase pack 根目录")
    parser.add_argument("--absorber-site", default=None,
                        help="P1b 吸收位标签或 'all'")
    parser.add_argument("--offline", action="store_true",
                        help="P1b 离线模式")
```

### 4.2 在 `xafs_report.py` 调用 `xafs_fit.run` 前接入 `load_t5_s02` 与 `phase_resolver` + `cif_gate`

**原文（第 1819–1837 行）**：
```python
    print("== S6.5 结构识别 + S7 拟合 ==")
    out_fit = xafs_fit.run(
        st["result"], identify=True, fit=not args.no_fit,
        absorber_symbol=args.absorber, fit_motif=args.fit_motif,
        fit_rmin=args.fit_rmin, fit_rmax=args.fit_rmax,
        allow_ambiguous=args.allow_ambiguous,
        l1_ligand=args.l1_ligand, l1_geometry=args.l1_geometry,
        fit_weights=args.fit_weights,
        shared_sigma2=shared_sigma2, kweights=fit_kweights,
        include_ms=args.include_ms, vary_c3=vary_c3, vary_c4=vary_c4,
        isotropic_deltar=args.isotropic_deltar,
        fit_space=args.fit_space, bvs_restraint=bvs_restraint,
        catalytic_model=args.catalytic_model,
        hetero_cocupancy=args.hetero_cocupancy,
        constraints=constraints,
        s02=args.s02,
        formal_cn=args.formal_cn,
        verbose=True)
```

**替换为**：
```python
    s02_val = args.s02
    s02_source = "--s02 runtime override" if s02_val is not None else None
    if s02_val is None and args.s02_from_t5:
        import run_xafs
        s02_val, s02_source = run_xafs.load_t5_s02(args.absorber)

    structure_model = None
    structure_models: list = []
    gate_result = None
    gate_results: list = []
    if args.formal_cn and (args.phase_cif or args.phase_id
                           or args.phase or args.phase_index):
        import cif_gate
        import phase_resolver
        resolution = phase_resolver.resolve_phase(
            args.absorber,
            phase_cif=args.phase_cif,
            phase_id=args.phase_id,
            phase_spec=args.phase,
            phase_index_path=args.phase_index,
            phase_pack_path=args.phase_pack,
            absorber_site=args.absorber_site,
            offline=args.offline)
        if resolution.status != "resolved":
            raise XAFSConfigError(
                "FIT", f"--phase-* 解析失败 ({resolution.code}): "
                       f"{resolution.reason}")
        structure_model = resolution.structure_model
        structure_models = list(resolution.site_models)
        multi_site = len(structure_models) > 1
        gate_kwargs = dict(
            kmin=args.kmin, kmax=args.kmax, rmax=6.0,
            allow_computed=False, allow_unusual_coordination=False,
            is_computed=not resolution.selected.experimental)
        if multi_site:
            gate_results = cif_gate.gate_cif_all_sites(
                resolution.selected.cif_path, args.absorber, **gate_kwargs)
            failed = [g for g in gate_results if not g.ok]
            if failed:
                raise XAFSConfigError(
                    "FIT", f"CIF gate rejected multi-site: {failed[0].reason}")
            gate_result = gate_results[0]
            for m, g in zip(structure_models, gate_results):
                m.occupancy_policy = g.occupancy_source
                m.n0_shells = list(g.n0_shells)
                m.warnings = list(g.warnings)
        else:
            gate_result = cif_gate.gate_cif(
                resolution.selected.cif_path, args.absorber,
                absorber_site=resolution.absorber_site, **gate_kwargs)
            if not gate_result.ok:
                raise XAFSConfigError(
                    "FIT", f"CIF gate rejected: {gate_result.code} — "
                           f"{gate_result.reason}")
            if structure_model is not None:
                structure_model.occupancy_policy = gate_result.occupancy_source
                structure_model.n0_shells = list(gate_result.n0_shells)
                structure_model.warnings = list(gate_result.warnings)

    multi_site_formal = bool(args.formal_cn and len(structure_models) > 1)
    print("== S6.5 结构识别 + S7 拟合 ==")
    out_fit = xafs_fit.run(
        st["result"], identify=True, fit=not args.no_fit,
        absorber_symbol=args.absorber, fit_motif=args.fit_motif,
        fit_rmin=args.fit_rmin, fit_rmax=args.fit_rmax,
        allow_ambiguous=args.allow_ambiguous,
        l1_ligand=args.l1_ligand, l1_geometry=args.l1_geometry,
        fit_weights=args.fit_weights,
        shared_sigma2=shared_sigma2, kweights=fit_kweights,
        include_ms=args.include_ms, vary_c3=vary_c3, vary_c4=vary_c4,
        isotropic_deltar=args.isotropic_deltar,
        fit_space=args.fit_space, bvs_restraint=bvs_restraint,
        catalytic_model=args.catalytic_model,
        hetero_cocupancy=args.hetero_cocupancy,
        constraints=constraints,
        s02=s02_val,
        s02_source=s02_source,
        formal_cn=args.formal_cn and not multi_site_formal,
        structure_model=structure_model,
        gate_result=gate_result,
        verbose=True)
    if multi_site_formal and out_fit.get("fit") is not None:
        multi = xafs_fit.run_multi_site_formal_cn(
            st["result"],
            site_models=structure_models,
            gate_results=gate_results,
            s02=s02_val,
            s02_source=s02_source,
            data_sha256=(getattr(st["result"], "data_sha256", "") or ""),
            verbose=True)
        out_fit["fit"]["formal_cn"] = multi
        out_fit["fit"]["cn_status"] = multi.get("status", "multi_site_failed")
```

---

## Step 5：新增回归契约 `#129`（单晶位 P5-C 氧化物端到端正式 CN 不被 L1 劫持 + `occupancy_policy` 同步）

文件：`C:\Users\lizha\Desktop\Code-experience-cn\regression_test.py`

在 `test_p6_all_phase_pack_entries_pass_cif_gate`（`#128`）末尾追加 `#129`：

```python
    def test_p6_1_single_site_oxide_formal_cn_end_to_end(self):
        """#129 P6.1：单晶位氧化物（不在 _structure_lib 中）通过已验证 CIF 完成正式 CN，不被 L1 劫持。"""
        from larch.xafs import feffpath
        from types import SimpleNamespace

        s02, s02_src = run_xafs.load_t5_s02("Ru")
        pr = phase_resolver.resolve_phase(
            "Ru", phase_id="COD:9007541", offline=True)
        self.assertEqual(pr.status, "resolved")
        self.assertEqual(pr.structure_model.occupancy_policy, "implicit_1.0")
        gate = cif_gate.gate_cif(
            pr.selected.cif_path, "Ru",
            absorber_site=pr.absorber_site, kmin=3.0, kmax=12.0)
        self.assertTrue(gate.ok, gate.reason)

        motif = xafs_fit._site_motif_from_model(pr.structure_model)
        plan = xafs_fit.select_shell_paths(
            motif, edge="K", kmin=3.0, kmax=12.0)
        k = np.arange(3.0, 12.0, 0.05)
        plist = []
        for sh, f in zip(plan["shells"], plan["files"]):
            p = feffpath(str(f), s02=s02, sigma2=0.004, e0=1.0, deltar=0.0)
            p.degen = float(sh["degen"])
            p._calc_chi()
            plist.append(p)
        chi = artemis.exafs_model_k(plist, k)
        rng = np.random.default_rng(42)
        chi = chi + rng.normal(0.0, 5e-4, size=chi.shape)
        chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                     mu_norm=np.ones_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
        res = SimpleNamespace(
            chi_k=chi_k, ft=ft, absorber_symbol="Ru", structure=None,
            merged=None, norm=None, self_absorption=None,
            data_sha256="b" * 64)

        out = xafs_fit.run(
            res, identify=True, fit=True, absorber_symbol="Ru",
            s02=s02, s02_source=s02_src, formal_cn=True,
            structure_model=pr.structure_model, gate_result=gate,
            verbose=False)
        # 诊断轨因 _structure_lib/Ru 只有 ru_hcp 返回 l1_first_shell，
        # 但正式轨必须独立用 COD:9007541 (RuO2) 跑通 P3 条件 CN 并输出 ok！
        self.assertEqual(out["fit"]["mode"], "l1_first_shell")
        self.assertEqual(out["cn_status"], "ok")
        self.assertEqual(out["formal_cn"]["status"], "ok")
        rows = out["formal_cn"]["rows"]
        self.assertGreaterEqual(len(rows), 1)
        self.assertEqual(rows[0]["phase_id"], "COD:9007541")
        self.assertAlmostEqual(float(rows[0]["n0"]), 6.0, places=3)
        self.assertAlmostEqual(float(rows[0]["n_fit"]), 6.0, delta=0.35)
```

---

## Step 6：检验、提交与补丁命令

```powershell
$py = "C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe"

# 1. 语法编译检查
& $py -m py_compile xafs_fit.py phase_resolver.py run_xafs.py xafs_report.py regression_test.py

# 2. 定向契约检查（#75 / #76 / #127 / #128 / #129）
& $py -m unittest `
    regression_test.RegressionTest.test_formal_cn_refusal_without_model `
    regression_test.RegressionTest.test_formal_cn_l1_not_in_table `
    regression_test.RegressionTest.test_p5d_multi_site_formal_cn `
    regression_test.RegressionTest.test_p6_all_phase_pack_entries_pass_cif_gate `
    regression_test.RegressionTest.test_p6_1_single_site_oxide_formal_cn_end_to_end

# 3. 全量回归 + 4 项压力测试（必须 129 契约 + 4 压力 PASS）
& $py run_ci.py --stress
```

---

## 7. 执行结果（2026-11-09 实测）

### 7.1 代码替换完成

| 文件 | 改动 |
|---|---|
| `xafs_fit.py` | 新增 `_t5_s02_hint()` 与 `_attach_single_site_formal_cn()`；`fit_identified` 的 L1 早退、motif 主路径与末尾挂载统一改走正式轨辅助函数；移除已被辅助函数接管的局部 `conditional_cn` / `n0_mapping` 逻辑；`_site_fractional` 增加 `source_cif` 回退；`run_operando_exafs_series` / `run_multi_site_formal_cn` 补 `--s02-from-t5` 提示。 |
| `phase_resolver.py` | `_candidate_to_structure_model` 按 CIF 是否含 `_atom_site_occupancy` 设置 `"explicit"` / `"implicit_1.0"`。 |
| `run_xafs.py` | `cif_gate` 通过后把 `occupancy_source` / `n0_shells` / `warnings` 回填单晶位与逐位点模型；多位点正式 CN 后同步 `result.fit["cn_status"]`。 |
| `xafs_report.py` | 新增 `--s02-from-t5`、`--phase-cif/--phase-id/--phase/--phase-index/--phase-pack/--absorber-site/--offline`；接入 `load_t5_s02`、`phase_resolver`、`cif_gate`；多位点分支调用 `run_multi_site_formal_cn` 并同步 `cn_status`。 |
| `regression_test.py` | 新增 `#129` `test_p6_1_single_site_oxide_formal_cn_end_to_end`。 |

### 7.2 端到端实测

- **单晶位 RuO2（`COD:9007541`，`Ru-K-RuO2.xdi`）**：
  - 诊断轨仍为 `l1_first_shell`（`_structure_lib/Ru` 只有 `ru_hcp`）；
  - 正式轨不再被劫持，`formal_cn.status=consistency_failed`，拒绝原因变为**物理原因**（`N/N₀=0.521 out_of_range`；`Rf=0.1314 > 0.08`），而非 `l1_diagnostic_only`；
  - `formal_cn.model.occupancy_policy="implicit_1.0"`。
- **单晶位 NiO（`COD:1010093`，`Ni-K-NickelOxide.xdi`）**：
  - 诊断轨 `motif`，正式轨 `ok`；3 条 shell 行，首壳 `Ni-O@2.08, N₀=6.0, N_fit=6.327, fitted`；
  - `occupancy_policy="explicit"`，`n0_shells=7` —— 与 P6 前行为一致，无回归。
- **多晶位 Fe3O4（`COD:1010369`，`Fe-Magnetite.xdi`，`--absorber-site all`）**：
  - `result.fit["cn_status"]="multi_site_failed"`，`_meta.json` 中 `fit.cn_status` 同为 `multi_site_failed`（P6 时残留 `diagnostic_only`）；
  - `site_weights={Fe1:1/3, Fe2:1/3, Fe3:1/3}`。
- **`xafs_report.py` CLI**：
  `--absorber Ru --no-strict --s02-from-t5 --formal-cn --phase-id COD:9007541` 可完整跑通并写出 `cn_table.csv`；拒绝行 `phase_id=COD:9007541`、`s02=0.7002`（Ru T5 标定值）、`gate_status=consistency_failed`。

### 7.3 全量检验

| 检验项 | 命令 | 结果 |
|---|---|---|
| 语法编译 | `py_compile xafs_fit.py phase_resolver.py run_xafs.py xafs_report.py regression_test.py` | PASS |
| 定向契约 | `#75 / #76 / #127 / #128 / #129` | 5/5 OK |
| 全量回归 | `regression_test.py` | **129 项 OK**（127.3 s） |
| CI + 压力 | `run_ci.py --stress` | **129 契约 + 4 压力 PASS**（186.8 s） |

`ci_report.json`：`{"ok": true, "contract_tests": 129, "stress_tests": 4}`。

### 7.4 边界确认

- `absorber.py` 静态表与 `calibrate_s02.py` 门禁常量未改动（`#60` 守住）；
- 44 条 phase pack 与 `Rh2O3=COD:1544346` 未改动；
- `_standards_lib/`、`data/standards_raw/`、`logs/` 未加入暂存区；
- 根目录 `lepidocrocite.cif` 保持原位；
- P6 遗留的根目录 `CIF_CRAWLER_INTEGRATION.md` 重复副本已删除，权威内容统一保留在 `docs/CIF_CRAWLER_INTEGRATION.md`（两份文件字节完全一致，sha256 前缀 `3b0dfeb23bcd09c9`）。

