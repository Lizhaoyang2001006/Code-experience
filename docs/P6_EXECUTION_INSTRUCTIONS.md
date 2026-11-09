# P6 执行指令：催化分析门禁修复（44/44 相全通）与总项目目录整合

> 本文件是 P6 的唯一权威执行记录。执行顺序：**先执行 Step 1–5，再运行 Step 6 检验**。
> 代码修改一律采用“原文精确替换”，禁止自行改写语义或放宽门禁。

---

## 0. 硬性工程守则（最高优先级）

1. **只允许操作** `C:\Users\lizha\Desktop\Code-experience-cn`；严禁修改 `C:\Users\lizha\Desktop\Code-experience`。
2. **禁止破坏性 Git 命令**：严禁 `git checkout -- .`、`git reset --hard`、`git clean -fd`。
3. **严禁放宽 S₀² 标定门禁**：
   - `absorber.py` 静态表中 12 个扩展元素必须保持 `s02=None, calibrated=False`（守住 `#60`）；
   - `calibrate_s02.py` 的门禁常量不得修改。
4. **禁止提交第三方原始谱**：`_standards_lib/`、`data/standards_raw/`、`logs/` 严禁 `git add`。
5. **固定 Python 解释器**：
   `C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe`
6. **零路径破坏**：
   - 根目录 `lepidocrocite.cif` 保持原位（被 `regression_test.py`、`xafs_report.py`、`_structure_lib/build_seed.py` 引用）；
   - 核心 `.py`、`_structure_lib/`、`data/`、`tools/` 保持原位。

---

## 1. 基线

- 当前主仓 HEAD：`13fba39`（P5-S02 + P5-D 闭环）
- 回归契约：127 项 OK；压力测试：4 项 PASS
- phase pack：44 条真实 CIF 条目
- 实测门禁失败：44 条中 15 条被 `cif_gate.gate_cif` 拒绝

已知断层：
1. 单晶位多元素晶体（22 条）不传 `--absorber-site` 时被误报 `absorber_site_ambiguous`；
2. 11 条有序多元素氧化物因省略 `_atom_site_occupancy` 被误报 `occupancy_missing`；
3. 6 条高价/重金属氧化物因金属单质共价半径被判 `bond_length_implausible`，且缺 `U` 半径；
4. `COD:1010584 (Rh2O3)` 为 1928 年粗糙结构，第一壳 `Rh–O=1.509 Å` 失真；
5. CLI 未提示 `--s02-from-t5`；`--catalytic-model + --formal-cn` 缺少 `hypothesis_model_only` 边界。

---

## Step 1：修复 `cif_gate.py` 三大运行期门禁断层

文件：`C:\Users\lizha\Desktop\Code-experience-cn\cif_gate.py`

### 1.1 补齐 `U` 共价半径并新增 `PAIR_BOND_REF`

**原文**：
```python
    "Rn": 1.82,
}

# 第一壳 CN 范围（§5.1）
CN_NORMAL_LO, CN_NORMAL_HI = 2, 12
CN_WARN_LO, CN_WARN_HI = 12, 16
BOND_TOL = 0.25           # 键长容差：共价半径之和 ±0.25 Å
OCCUPANCY_TOL = 0.08      # 共占位合并距离容差（Å），与 artemis 一致
```

**替换为**：
```python
    "Rn": 1.82, "U": 1.70,
}

# 第一壳 CN 范围（§5.1）
CN_NORMAL_LO, CN_NORMAL_HI = 2, 12
CN_WARN_LO, CN_WARN_HI = 12, 16
BOND_TOL = 0.25           # 键长容差：共价半径之和 ±0.25 Å
OCCUPANCY_TOL = 0.08      # 共占位合并距离容差（Å），与 artemis 一致

# 常见催化金属–配体第一壳实验键长允许区间 (Å)，解决高价态/重金属氧化物
# 与金属单质共价半径之和系统性偏离 0.3–0.45 Å 的问题。
PAIR_BOND_REF: dict[tuple[str, str], tuple[float, float]] = {
    ("Cr", "O"): (1.85, 2.15),
    ("Mn", "O"): (1.80, 2.30),
    ("Fe", "O"): (1.78, 2.25),
    ("Co", "O"): (1.80, 2.22),
    ("Ni", "O"): (1.90, 2.20),
    ("Cu", "O"): (1.78, 2.10),
    ("Zn", "O"): (1.88, 2.12),
    ("Mo", "O"): (1.65, 2.15),
    ("Ru", "O"): (1.88, 2.10),
    ("Rh", "O"): (1.92, 2.15),
    ("Pd", "O"): (1.92, 2.15),
    ("Ag", "O"): (1.95, 2.20),
    ("Ir", "O"): (1.90, 2.12),
    ("Pt", "O"): (1.90, 2.12),
    ("Au", "O"): (1.90, 2.15),
    ("Ce", "O"): (2.20, 2.50),
    ("U", "O"):  (1.75, 2.50),
}
PAIR_BOND_SLACK = 0.08    # Å，实验参考区间外的测量容差
```

### 1.2 新增化学计量比核验函数

在 `def _absorber_positions(st, cell, absorber: str) -> list:` **之前**插入：

```python
_FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")


def _formula_matches_unit_cell(block, st) -> bool:
    """检查 CIF 声明的 _chemical_formula_sum 是否与隐式 occ=1.0 的晶胞计量比严格一致。

    注意：若 CIF 未声明 _chemical_formula_sum（如 #98 测试的无化学式缺陷 CIF），
    必须返回 False，继续拒绝为 occupancy_missing。
    """
    raw = block.find_value("_chemical_formula_sum")
    if not raw:
        return False
    text = str(raw).strip().strip("'").strip('"').strip()
    if not text or text in ("?", "."):
        return False
    declared: dict[str, float] = {}
    for el, num in _FORMULA_TOKEN_RE.findall(text):
        if el in ("H", "D"):
            continue
        val = float(num) if num else 1.0
        if val <= 0:
            return False
        declared[el] = declared.get(el, 0.0) + val
    if not declared:
        return False

    cell_counts: dict[str, float] = {}
    for site in st.get_all_unit_cell_sites():
        el = _symbol(site)
        if not el or el in ("H", "D"):
            continue
        occ = float(site.occ) if site.occ is not None else 1.0
        if abs(occ - 1.0) > 1e-3:
            return False
        cell_counts[el] = cell_counts.get(el, 0.0) + 1.0

    if set(declared.keys()) != set(cell_counts.keys()):
        return False
    d_sum = sum(declared.values())
    c_sum = sum(cell_counts.values())
    if d_sum <= 0 or c_sum <= 0:
        return False
    for el in declared:
        if abs(declared[el] / d_sum - cell_counts[el] / c_sum) > 0.01:
            return False
    return True
```

### 1.3 修复单晶位对称展开误报 `absorber_site_ambiguous`

**原文**：
```python
    # ---- 4. 吸收位：同元素多个晶位时必须 --absorber-site ----
    groups = _absorber_positions(st, st.cell, absorber)
    if len(groups) == 0:
        return CIFGateResult(
            ok=False, code="absorber_missing", absorber=absorber,
            reason=f"{cif_path.name} 中找不到元素 {absorber} 的位点",
            cif_path=result.cif_path, cif_sha256=result.cif_sha256)
    if len(groups) > 1 and absorber_site is None:
        # 单元素结构：所有吸收元素位点对称等价，自动选第一个，不报 ambiguous
        if len(elements) == 1:
            result.warnings.append(
                f"{cif_path.name} 中元素 {absorber} 有 {len(groups)} 个对称等价晶位，"
                f"自动选 {groups[0]['members'][0].label}")
            abs_site = groups[0]["members"][0]
            result.absorber_site = str(abs_site.label)
        else:
            site_labels = []
            for g in groups:
                rep = g["members"][0]
                site_labels.append(str(rep.label))
            return CIFGateResult(
                ok=False, code="absorber_site_ambiguous", absorber=absorber,
                reason=(f"{cif_path.name} 中元素 {absorber} 有 {len(groups)} 个晶位 "
                        f"{site_labels}，必须 --absorber-site 指定；禁止按 occ×多重度自动挑主位点"),
                cif_path=result.cif_path, cif_sha256=result.cif_sha256)
```

**替换为**：
```python
    # ---- 4. 吸收位：不对称单元中同元素存在多个不等价晶位时必须 --absorber-site ----
    groups = _absorber_positions(st, st.cell, absorber)
    asym_abs_sites = [s for s in st.sites if _symbol(s) == absorber]
    if len(groups) == 0 or len(asym_abs_sites) == 0:
        return CIFGateResult(
            ok=False, code="absorber_missing", absorber=absorber,
            reason=f"{cif_path.name} 中找不到元素 {absorber} 的位点",
            cif_path=result.cif_path, cif_sha256=result.cif_sha256)
    if len(asym_abs_sites) > 1 and absorber_site is None:
        # 单元素结构：所有吸收元素位点对称等价，自动选第一个，不报 ambiguous
        if len(elements) == 1:
            result.warnings.append(
                f"{cif_path.name} 中元素 {absorber} 有 {len(asym_abs_sites)} 个对称等价晶位，"
                f"自动选 {asym_abs_sites[0].label}")
            abs_site = asym_abs_sites[0]
            result.absorber_site = str(abs_site.label)
        else:
            site_labels = [str(s.label) for s in asym_abs_sites]
            return CIFGateResult(
                ok=False, code="absorber_site_ambiguous", absorber=absorber,
                reason=(f"{cif_path.name} 中元素 {absorber} 有 {len(asym_abs_sites)} 个晶位 "
                        f"{site_labels}，必须 --absorber-site 指定；禁止按 occ×多重度自动挑主位点"),
                cif_path=result.cif_path, cif_sha256=result.cif_sha256)
```

**并且**将紧随其后的 `else:` 分支：

**原文**：
```python
    else:
        # 只有一个晶位，自动选
        abs_site = groups[0]["members"][0]
        result.absorber_site = str(abs_site.label)
```

**替换为**：
```python
    else:
        # 只有一个不等价晶位，自动选不对称单元代表位点
        abs_site = asym_abs_sites[0]
        result.absorber_site = str(abs_site.label)
```

### 1.4 修复隐式占据度与第一壳键长检查

**原文（occupancy 分支）**：
```python
    if not has_occ_tag:
        # 缺失 occupancy：判断是否允许隐式 1.0
        single_elem = len(elements) == 1
        ordered = not cocupancy
        # 化学式一致：单元素且没有共占位，化学式自然一致
        formula_consistent = single_elem and ordered
        if single_elem and ordered and formula_consistent:
            result.occupancy_source = "implicit_1.0"
            result.warnings.append(
                f"{cif_path.name} 缺少 _atom_site_occupancy，但为单元素有序结构，"
                f"隐式 occupancy=1.0")
        else:
            return CIFGateResult(
                ok=False, code="occupancy_missing", absorber=absorber,
                reason=(f"{cif_path.name} 缺少 _atom_site_occupancy，且不是单元素有序结构；"
                        f"无序/共占位/缺陷模型/化学式不匹配 → occupancy_missing"),
                cif_path=result.cif_path, cif_sha256=result.cif_sha256)
```

**替换为**：
```python
    if not has_occ_tag:
        # 缺失 occupancy：单元素有序，或声明化学式与全占据晶胞计量比严格一致的有序晶体允许隐式 1.0
        single_elem = len(elements) == 1
        ordered = not cocupancy
        formula_consistent = ordered and (
            single_elem or _formula_matches_unit_cell(block, st)
        )
        if ordered and formula_consistent:
            result.occupancy_source = "implicit_1.0"
            result.warnings.append(
                f"{cif_path.name} 缺少 _atom_site_occupancy，但为有序化学计量结构，"
                f"隐式 occupancy=1.0")
        else:
            return CIFGateResult(
                ok=False, code="occupancy_missing", absorber=absorber,
                reason=(f"{cif_path.name} 缺少 _atom_site_occupancy，且不是有序化学计量结构；"
                        f"无序/共占位/缺陷模型/化学式不匹配 → occupancy_missing"),
                cif_path=result.cif_path, cif_sha256=result.cif_sha256)
```

**原文（键长分支）**：
```python
    # ---- 8. 键长合理性（第一壳） ----
    first_el = first_shell["el"]
    if first_el != absorber:
        r_sum = (COVALENT_RADIUS.get(absorber) and COVALENT_RADIUS.get(first_el))
        if COVALENT_RADIUS.get(absorber) and COVALENT_RADIUS.get(first_el):
            expected = COVALENT_RADIUS[absorber] + COVALENT_RADIUS[first_el]
            r_obs = first_shell["r"]
            if abs(r_obs - expected) > BOND_TOL:
                return CIFGateResult(
                    ok=False, code="bond_length_implausible", absorber=absorber,
                    reason=(f"第一壳 {absorber}–{first_el} 键长 {r_obs:.3f} Å 与共价半径之和 "
                            f"{expected:.3f} Å 偏差 > {BOND_TOL} Å；未知元素对请提供参考键长"),
                    cif_path=result.cif_path, cif_sha256=result.cif_sha256)
```

**替换为**：
```python
    # ---- 8. 键长合理性（第一壳） ----
    first_el = first_shell["el"]
    if first_el != absorber:
        r_obs = first_shell["r"]
        pair_ref = (PAIR_BOND_REF.get((absorber, first_el))
                    or PAIR_BOND_REF.get((first_el, absorber)))
        if pair_ref is not None:
            r_lo, r_hi = pair_ref
            if (r_obs < r_lo - PAIR_BOND_SLACK
                    or r_obs > r_hi + PAIR_BOND_SLACK):
                return CIFGateResult(
                    ok=False, code="bond_length_implausible", absorber=absorber,
                    reason=(f"第一壳 {absorber}–{first_el} 键长 {r_obs:.3f} Å "
                            f"超出实验参考范围 [{r_lo:.2f}, {r_hi:.2f}] Å"
                            f"（容差 ±{PAIR_BOND_SLACK} Å）"),
                    cif_path=result.cif_path, cif_sha256=result.cif_sha256)
        elif COVALENT_RADIUS.get(absorber) and COVALENT_RADIUS.get(first_el):
            r_expected = COVALENT_RADIUS[absorber] + COVALENT_RADIUS[first_el]
            if abs(r_obs - r_expected) > BOND_TOL:
                return CIFGateResult(
                    ok=False, code="bond_length_implausible", absorber=absorber,
                    reason=(f"第一壳 {absorber}–{first_el} 键长 {r_obs:.3f} Å "
                            f"与共价半径之和 {r_expected:.3f} Å 偏差 > {BOND_TOL} Å；"
                            f"未知元素对请提供参考键长"),
                    cif_path=result.cif_path, cif_sha256=result.cif_sha256)
```

---

## Step 2：替换失真的 `Rh2O3`（`COD:1010584` → `COD:1544346`）

依据：Shannon & Prewitt 1970, *J. Solid State Chem.* 2, 134–136；`P b n a` IT #60；
显式 `_atom_site_occupancy 1.0`；单晶位 `Rh`；第一壳 `6 × Rh–O @ 2.050 Å`。

### 2.1 修改目标表

文件：`tools/cif_crawler.py`

**原文**：
```python
    "Rh2O3": {
        "tier": "P5C", "phase": "corundum", "absorber": "Rh",
        "formula": "Rh2O3", "stoich": {"Rh": 2, "O": 3},
        "cod_id": "1010584", "spg_nums": [167],
        "tags": ["oxide", "corundum", "catalyst"],
    },
```

**替换为**：
```python
    "Rh2O3": {
        "tier": "P5C", "phase": "orthorhombic-II", "absorber": "Rh",
        "formula": "Rh2O3", "stoich": {"Rh": 2, "O": 3},
        "cod_id": "1544346", "spg_nums": [60],
        "tags": ["oxide", "corundum", "catalyst"],
    },
```

### 2.2 替换缓存与 phase pack 条目

```powershell
$py = "C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe"
$dst = "C:\Users\lizha\Desktop\Code-experience-cn\data\cif_raw\Rh2O3_COD_1544346.cif"
Copy-Item "C:\Users\lizha\AppData\Local\Temp\opencode\pool_probe\1544346.cif" $dst -Force
Remove-Item "C:\Users\lizha\Desktop\Code-experience-cn\_structure_lib\_phase_pack\COD\1010584" -Recurse -Force
& $py tools/cif_crawler.py --target Rh2O3 --offline --replace
```

`rebuild_manifest()` 会扫描 `*/*/*.meta.json`，旧 `1010584` 目录删除后自动从 manifest 消失，
新 `1544346` 条目自动写入；总条目数必须仍为 **44**。

### 2.3 同步文档与谱源池文字

- `tools/build_source_pool.py`：`Rh2O3 corundum（P5-C COD 1010584）配套谱`
  → `Rh2O3（P5-C COD 1544346）配套谱`
- `README.md` 与 `CIF_CRAWLER_INTEGRATION.md` 的 Rh2O3 行：
  `1010584 | 167 | corundum` → `1544346 | 60 | orthorhombic-II`
- 重新生成谱源池：
```powershell
& $py tools/build_source_pool.py
```

---

## Step 3：CLI 提示与催化假设模型边界

### 3.1 `structure_id.py` 增加 `--s02-from-t5` 提示

文件：`structure_id.py`

**原文**：
```python
    if not absorber.calibrated:
        raise XAFSDataError(
            "structure_id", f"{absorber.symbol} 未标定 S₀²，"
                            f"识别结果不能用于 CN 定量，拒绝运行")
```

**替换为**（保留异常类型与前导关键字）：
```python
    if not absorber.calibrated:
        t5_cal = (Path(__file__).resolve().parent / "data"
                  / "s02_calibration" / absorber.symbol / "calibration.json")
        t5_hint = (
            "（该元素已有 T5 独立标定值，CLI 可加 --s02-from-t5 或显式传 --s02）"
            if t5_cal.is_file() else "")
        raise XAFSDataError(
            "structure_id", f"{absorber.symbol} 未标定 S₀²，"
                            f"识别结果不能用于 CN 定量，拒绝运行{t5_hint}")
```

### 3.2 `xafs_fit.py` 增加 `hypothesis_model_only`

文件：`xafs_fit.py`，在 `fit_identified()` 末尾、`return out` **之前**插入：

```python
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
```

---

## Step 4：回归契约更新 `#127` 与新增 `#128`

文件：`regression_test.py`

### 4.1 更新 `#127` 的 `Co3O4` 断言

**原文**：
```python
        # Co3O4：8a:16d → 1/3 : 2/3；但该 COD CIF 缺少
        # _atom_site_occupancy，冻结门禁必须拒绝（无合格来源不出 CN）。
        co = phase_resolver.resolve_phase(
            "Co", phase_id="COD:9005888", absorber_site="all", offline=True)
        self.assertEqual([m.absorber_site for m in co.site_models],
                         ["Co1", "Co2"])
        co_sites = cif_gate.list_inequivalent_absorber_sites(
            co.selected.cif_path, "Co")
        self.assertEqual([s["label"] for s in co_sites], ["Co1", "Co2"])
        self.assertEqual([s["multiplicity"] for s in co_sites], [8, 16])
        self.assertAlmostEqual(co_sites[0]["weight"], 1.0 / 3.0, places=12)
        self.assertAlmostEqual(co_sites[1]["weight"], 2.0 / 3.0, places=12)
        co_gates = cif_gate.gate_cif_all_sites(
            co.selected.cif_path, "Co", kmin=3.0, kmax=12.0)
        self.assertTrue(all(not g.ok for g in co_gates))
        self.assertTrue(all(g.code == "occupancy_missing" for g in co_gates))
```

**替换为**：
```python
        # Co3O4：8a:16d → 1/3 : 2/3；声明化学式 Co3 O4 与全占据晶胞一致，
        # 按 implicit_1.0 放行，加权第一壳 N₀ = 1/3×4 + 2/3×6 = 16/3。
        co = phase_resolver.resolve_phase(
            "Co", phase_id="COD:9005888", absorber_site="all", offline=True)
        self.assertEqual([m.absorber_site for m in co.site_models],
                         ["Co1", "Co2"])
        co_sites = cif_gate.list_inequivalent_absorber_sites(
            co.selected.cif_path, "Co")
        self.assertEqual([s["label"] for s in co_sites], ["Co1", "Co2"])
        self.assertEqual([s["multiplicity"] for s in co_sites], [8, 16])
        self.assertAlmostEqual(co_sites[0]["weight"], 1.0 / 3.0, places=12)
        self.assertAlmostEqual(co_sites[1]["weight"], 2.0 / 3.0, places=12)
        co_gates = cif_gate.gate_cif_all_sites(
            co.selected.cif_path, "Co", kmin=3.0, kmax=12.0)
        self.assertTrue(all(g.ok for g in co_gates),
                        [g.reason for g in co_gates])
        self.assertTrue(all(g.occupancy_source == "implicit_1.0"
                            for g in co_gates))
        co_weighted_n0 = sum(g.site_weight * g.first_shell_cn
                             for g in co_gates)
        self.assertAlmostEqual(co_weighted_n0, 16.0 / 3.0, places=6)
```

### 4.2 新增 `#128`

在 `#127` 之后、`if __name__ == "__main__":` 之前插入：

```python
    def test_p6_all_phase_pack_entries_pass_cif_gate(self):
        """#128 P6：phase_pack 全部 44 个入库晶体相 100% 通过 cif_gate。"""
        pack_root = ROOT / "_structure_lib" / "_phase_pack"
        manifest = phase_pack.load_manifest(pack_root)
        self.assertEqual(len(manifest["entries"]), 44)

        single_count = 0
        multi_count = 0
        for entry in manifest["entries"]:
            el = entry["absorber_elements"][0]
            cif_path = pack_root / entry["cif_path"]
            sites = entry["absorber_sites"]
            is_comp = not entry.get("experimental", True)
            if len(sites) == 1:
                single_count += 1
                g = cif_gate.gate_cif(cif_path, el, absorber_site=None,
                                      is_computed=is_comp)
                self.assertTrue(
                    g.ok, f"{entry['id']} {entry['formula']} failed: "
                          f"{g.code} ({g.reason})")
                self.assertEqual(g.absorber_site, sites[0])
                self.assertGreaterEqual(g.first_shell_cn, 2.0)
                self.assertLessEqual(g.first_shell_cn, 16.0)
            else:
                multi_count += 1
                g_amb = cif_gate.gate_cif(cif_path, el, absorber_site=None,
                                          is_computed=is_comp)
                parsed = phase_pack.parse_cif(cif_path)
                if len(parsed["formula_elements"]) == 1:
                    # 单元素多晶位：所有位点对称等价，自动选第一个
                    self.assertTrue(g_amb.ok, g_amb.reason)
                else:
                    # 多元素多晶位：未传 --absorber-site 必须继续拒绝
                    self.assertFalse(g_amb.ok)
                    self.assertEqual(g_amb.code, "absorber_site_ambiguous")
                gs = cif_gate.gate_cif_all_sites(
                    cif_path, el, is_computed=is_comp)
                self.assertEqual(len(gs), len(sites))
                self.assertTrue(
                    all(g.ok for g in gs),
                    f"{entry['id']} {entry['formula']} multi-site failed: "
                    f"{[(g.absorber_site, g.code, g.reason) for g in gs]}")
                self.assertAlmostEqual(
                    sum(g.site_weight for g in gs), 1.0, places=9)
        self.assertEqual(single_count, 39)
        self.assertEqual(multi_count, 5)

        # 催化位点假设模型与 --formal-cn 同时开启：明确写 hypothesis_model_only
        motif = structure_lib.parse_catalytic_cluster_spec(
            "Fe", "sac:N4@2.00,axial=O@2.10,C=none")
        plan = xafs_fit.select_shell_paths(
            motif, edge="K", rmin=1.0, rmax=2.6, kmin=3.0, kmax=12.0)
        k = np.arange(3.0, 12.0, 0.05)
        plist = []
        from larch.xafs import feffpath
        for sh, f in zip(plan["shells"], plan["files"]):
            p = feffpath(str(f), s02=0.803, sigma2=0.004, e0=2.0,
                         deltar=0.0)
            p.degen = float(sh["degen"])
            p._calc_chi()
            plist.append(p)
        chi = artemis.exafs_model_k(plist, k)
        rng = np.random.default_rng(7)
        chi = chi + rng.normal(0.0, 1e-3, size=chi.shape)
        chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                     mu_norm=np.ones_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
        from types import SimpleNamespace
        res = SimpleNamespace(
            chi_k=chi_k, ft=ft, absorber_symbol="Fe", structure=None,
            merged=None, norm=None, self_absorption=None)
        out = xafs_fit.fit_identified(
            res, absorber_symbol="Fe",
            catalytic_model="sac:N4@2.00,axial=O@2.10,C=none",
            rmin=1.0, rmax=2.6, verbose=False, formal_cn=True)
        self.assertEqual(out["formal_cn"]["status"], "hypothesis_model_only")
        self.assertEqual(out["cn_status"], "hypothesis_model_only")
```

---

## Step 5：项目整合与根目录清理

### 5.1 收编文档到 `docs/`

```powershell
$docs = "C:\Users\lizha\Desktop\Code-experience-cn\docs"
New-Item -ItemType Directory -Force -Path $docs | Out-Null
$src = "C:\Users\lizha\Desktop\CN_CIF_执行方案"
foreach ($f in @("ACCEPTANCE_LOG.md","EXECUTION_PLAN.md",
                 "IMPLEMENTATION_CHECKLIST.md","P4_PLUS_EXECUTION_PLAN.md",
                 "PHASE_PACK_README.md","REVIEW_AND_REVISED_PLAN.md",
                 "phase_pack_manifest.schema.json")) {
    Copy-Item (Join-Path $src $f) (Join-Path $docs $f) -Force
}
Copy-Item (Join-Path $src "README.md") (Join-Path $docs "CN_CIF_PLAN_README.md") -Force
Copy-Item "C:\Users\lizha\Desktop\Code-experience-cn\CIF_CRAWLER_INTEGRATION.md" $docs -Force
```

### 5.2 归档根目录日志与临时输出

```powershell
$root = "C:\Users\lizha\Desktop\Code-experience-cn"
New-Item -ItemType Directory -Force -Path (Join-Path $root "logs\archive") | Out-Null
Get-ChildItem (Join-Path $root "*.log") -File |
    Move-Item -Destination (Join-Path $root "logs\archive") -Force
if (Test-Path (Join-Path $root "xafs_out_batch")) {
    Move-Item (Join-Path $root "xafs_out_batch") (Join-Path $root "logs\archive\xafs_out_batch") -Force
}
```

禁止移动：
- 根目录 `lepidocrocite.cif`；
- `_structure_lib/_dynamic_cache/**/feffrun_feff.log` 等已跟踪缓存日志。

### 5.3 `.gitignore`

在末尾追加：
```gitignore
logs/
```

### 5.4 `README.md` 导航

在文件顶部增加“项目目录结构与文档导航”小节，列明：
- 核心入口：`run_xafs.py`、`calibrate_t5.py`、`run_ci.py`
- 构建期工具：`tools/cif_crawler.py`、`tools/build_source_pool.py`
- 文档：`docs/REVIEW_AND_REVISED_PLAN.md`、`docs/P4_PLUS_EXECUTION_PLAN.md`、
  `docs/ACCEPTANCE_LOG.md`、`docs/CIF_CRAWLER_INTEGRATION.md` 等

### 5.5 验收记录

在 `docs/ACCEPTANCE_LOG.md` 追加 P6 章节（范围 / 结果 / 测试 / 提交 / 结论），
在 `docs/P4_PLUS_EXECUTION_PLAN.md` 追加 P6 小节。外部
`C:\Users\lizha\Desktop\CN_CIF_执行方案` 保持只读，不修改。

---

## Step 6：检验、提交与补丁

按顺序执行，任一步失败即停止并修复，不得跳过。

### 6.1 回归测试（必须 `Ran 128 tests ... OK`）

```powershell
$py = "C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe"
& $py regression_test.py
```

### 6.2 CI + 压力测试（必须 128 + 4 PASS）

```powershell
& $py run_ci.py --stress
```

### 6.3 提交

```powershell
git add cif_gate.py structure_id.py xafs_fit.py regression_test.py .gitignore README.md CIF_CRAWLER_INTEGRATION.md tools/cif_crawler.py tools/build_source_pool.py data/cif_raw/Rh2O3_COD_1544346.cif data/s02_calibration/source_pool/index.json _structure_lib/_phase_pack/ _structure_lib/phase_index.json _structure_lib/phase_coverage.json docs/
git status --short
git commit -m "P6: fix phase-pack CIF gates (44/44), update Rh2O3 COD 1544346, and consolidate docs"
```

`git status --short` 中不得出现 `_standards_lib/`、`data/standards_raw/`、`logs/`。

### 6.4 导出补丁并反向校验

```powershell
git format-patch -1 HEAD --stdout | Set-Content -Encoding utf8 "patches/P6-catalytic-readiness-and-cleanup.patch"
git apply --check --reverse "patches/P6-catalytic-readiness-and-cleanup.patch"
```

### 6.5 最终验收清单

- [ ] `cif_gate.py` 含 `PAIR_BOND_REF`、`_formula_matches_unit_cell`；
- [ ] 单晶位 39 条不传 `--absorber-site` 全部 `ok=True`；
- [ ] 多晶位 5 条：4 条多元素继续 `absorber_site_ambiguous`，Mn 单元素自动选位；
- [ ] `Rh2O3` 条目为 `COD:1544346`，第一壳 `6 × Rh–O @ 2.050 Å`；
- [ ] phase pack 总数仍为 44；
- [ ] `regression_test.py`：128 OK；
- [ ] `run_ci.py --stress`：128 + 4 PASS；
- [ ] `docs/` 收录 9 份文档，根目录无 `*.log` 散落；
- [ ] `_standards_lib/`、`data/standards_raw/`、`logs/` 未被提交。
