# CN_CIF 执行验收日志

> 配套文档：`REVIEW_AND_REVISED_PLAN.md` §12
> 用法：每个批次开始前记录基线，完成后追加验收记录；不要覆盖历史记录。
> 规则：任何一项基线不匹配，先停下排查，不得继续改代码。

---

## 0. 冻结基线（P0a 开始前填写一次）

- 日期：2026-11-07
- 执行目录：`C:\Users\lizha\Desktop\Code-experience-cn`
- 副本基线 commit：`3b777c7 baseline: Phase V working tree before CN formalization`
- Python：`C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe`（报告 3.13.14）
- 操作系统：Windows / win32

### 0.1 CI

```text
命令：& $PY run_ci.py --stress
ci_report.json:
  ok=true
  contract_tests=74
  stress_tests=4
退出码：0（CI PASS，122.8s）
```

### 0.2 默认 LEPIDO 三扫描

```text
命令：& $PY xafs_report.py data/练手数据/02_Fe_standards/Fe_lepidocrocite.000 data/练手数据/02_Fe_standards/Fe_lepidocrocite.100 data/练手数据/02_Fe_standards/Fe_lepidocrocite.200 --absorber Fe --no-strict -o xafs_out/_baseline
r_factor=0.027281059410161693
第一壳 N=6.348637726145664
第一壳 R=2.0019558444559276
7 列拟合表路径：`xafs_out/_baseline/Fe_lepidocrocite_fit_table.csv`
7 列拟合表 SHA256=E098D0EDC24E49296710E44518F8A8F010DCD03D7D7EB8CB2E26831A7532E760
```

### 0.3 #60 静态 s02

```text
命令：& $PY -c "import absorber; ..."
结果：OK（12 个未标定元素均为 s02=None、calibrated=False）
```

### 0.4 工作区

```text
git status --short 摘要：基线提交后工作区干净
未跟踪产物摘要：已随基线提交纳入副本（_structure_lib/Au、Cr、Mo、Zn、氧化物 CIF/JSON、_dynamic_cache）
备份位置：原目录 `C:\Users\lizha\Desktop\Code-experience` 保持未改动
```

---

## 1. P0a：正式/诊断双轨止血

### 1.1 改动

- 新增文件：`structure_model.py`（`CNError` / `ShellN0` / `StructureEdit` / `StructureModel` / `sha256_file` / `model_from_motif` / `validate_formal`）、`cn_table.py`（38 列 `CN_TABLE_COLUMNS` / `refusal_row` / `write_cn_table` / `read_cn_table`）
- 修改文件：`xafs_fit.py`（占据度失败硬报错；`_formal_cn_meta`；`fit_identified` / `run` 增加 `formal_cn`、`structure_model`；L1 标 `l1_diagnostic_only`；无模型标 `no_cif`；合格模型标 `pending_consistency`）、`run_xafs.py`（`--formal-cn`、`process_sample(formal_cn=...)`、单样品与批处理写 `cn_table.csv`）、`xafs_report.py`（`--formal-cn`、透传 `xafs_fit.run`、写 `cn_table.csv`、与 `--no-fit` 冲突报错）、`xafs_config.py`（`DEFAULTS["formal_cn"]=False`）、`regression_test.py`（#75–#78）
- 新增接口：`--formal-cn`（`run_xafs.py` / `xafs_report.py`）；`xafs_fit.run(..., formal_cn=False, structure_model=None)`；`xafs_fit.fit_identified(..., formal_cn=False, structure_model=None)`；`StructureModel` / `validate_formal`；`cn_table.write_cn_table` / `refusal_row`
- 回滚 commit：`3b777c7 baseline: Phase V working tree before CN formalization`
- P0a commit：`4801140 P0a: formal CN / diagnostic dual track`
- 补丁：`patches/P0a.patch`（SHA256=`1E7510683CC1DFEFDF3F0CAA416883AF09541EBAAD22744E93B8A57A4209C1D1`；`git apply --check --reverse` 通过）

### 1.2 测试

```text
新增测试：
  - test_formal_cn_refusal_without_model: OK（无模型 → no_cif；诊断 fit 保留；拒绝行 n_fit 为空）
  - test_formal_cn_l1_not_in_table: OK（L1 → l1_diagnostic_only；不进入正式表）
  - test_occupancy_scale_failure_hard_fail: OK（占据度展开失败 → CNError(occupancy_scale_failed)）
  - test_cn_table_schema_and_refusal_row: OK（38 列顺序 + 拒绝行写读往返）

单独运行：#75–#78 → Ran 4 tests in 20.063s — OK
& $PY run_ci.py --stress
  ok=true
  contract_tests=78
  stress_tests=4
  退出码=0（CI PASS，134.3s）
```

### 1.3 零退化

```text
LEPIDO r_factor=0.027281059410161693
LEPIDO N=6.348637726145664
LEPIDO R=2.0019558444559276
fit_table SHA256=E098D0EDC24E49296710E44518F8A8F010DCD03D7D7EB8CB2E26831A7532E760
#60 静态 s02：OK（12 个未标定元素仍为 s02=None、calibrated=False）
诊断标记：diagnostic_only=true、cn_status=diagnostic_only
```

### 1.4 正式 CN 拒绝路径

```text
命令：& $PY run_xafs.py data/练手数据/02_Fe_standards/Fe_lepidocrocite.000 --absorber Fe --no-strict --identify --fit --formal-cn -o xafs_out/_p0a_formal
cn_table.csv 行数（CN 行）=0（仅 1 行拒绝记录）
gate_status=no_cif
reject_reason=未提供 StructureModel（无合格结构来源）
n_fit=（空）
诊断拟合是否仍生成：是（fit.csv / fit.json 均生成）

对照：
- 不带 --formal-cn：不生成 cn_table.csv
- xafs_report.py ... --formal-cn：生成 1 行 no_cif 拒绝记录
- xafs_report.py ... --formal-cn --no-fit：parser.error
```

### 1.5 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：正式 CN 数、CIF 门禁、N₀ 与固定 N 一致性按计划留到 P0b/P2/P3；P0a 只保证拒绝路径与表骨架
- 下一步：P0b（motif provenance 升级），预期 `contract_tests=82`、`stress_tests=4`

---

## 2. P0b：motif provenance 升级

### 2.1 改动

- 新增文件：`motif_provenance.py`（`PROVENANCE_STATUSES` / `MODEL_LEVELS` / `OCCUPANCY_SOURCES` / `sha256_file` / `is_generated_cif` / `absorber_site_label` / `classify_motif`）、`tools/backfill_motif_provenance.py`（幂等回填 + `--check`）
- 修改文件：`structure_lib.py`（`Motif` 新增 `source_ref`、`source_cif_sha256`、`absorber_site`、`occupancy_source`、`model_level`、`provenance_status`、`edits`；`register()` 覆盖注册时重算并保留 provenance）、`structure_model.py`（`model_from_motif` 接入分类；`validate_formal` 新增 `generated_cif_not_formal`；`sha256_file` 改由 `motif_provenance` 提供）、`regression_test.py`（#79–#82）
- 新增数据：`_structure_lib/provenance_report.json`、`_structure_lib/motif_coverage.json`；32 条 motif JSON 已回填
- 新增接口：`motif_provenance.classify_motif(motif)`；`structure_model.validate_formal()` 的 `generated_cif_not_formal` 拒绝码；`tools/backfill_motif_provenance.py --check`
- 回滚 commit：`4801140 P0a: formal CN / diagnostic dual track`
- P0b commit：`d34a44e P0b: classify motif provenance and reject generated CIFs`
- 补丁：`patches/P0b.patch`（SHA256=`EFD926630FFCF5C54146C22AA2928F3B5C3BBCB45C4E5822D21DDB4CFAC64825`；`git apply --check --reverse` 通过）

### 2.2 来源分层结果

```text
motif 总数=32
by_status:
  generated_only=21
  source_ref_only=2   （goethite=COD:2211652、lepidocrocite=COD:9002763）
  constructed=9       （cluster/SAC/DAC/alloy/mixture）
  verified_source_cif=0
  unverified=0
unclassified=0
真实源 CIF=0
覆盖吸收元素=9/17；完全缺失=8（Ag、Ce、Ir、Mn、Pd、Rh、Ru、U）
```

说明：P0b 不把现场生成 CIF 伪装成真实来源；`verified_source_cif` 要到 P1a phase pack / 用户 `--phase-cif` 才出现。`tools/backfill_motif_provenance.py --check` 在完整 `run_ci.py --stress` 之后退出码为 0，证明压力测试重注册不会丢 provenance。

### 2.3 测试

```text
新增测试：
  - test_motif_provenance_schema: OK（32/32 可加载；新字段齐全；状态/层级合法；报告一致）
  - test_generated_cif_not_formal_source: OK（fe_bcc=generated_only；lepidocrocite=source_ref_only；生成 CIF 即使哈希正确也被 generated_cif_not_formal 拒绝）
  - test_user_cif_is_verified_source: OK（临时用户 CIF + sha256 + 吸收位 → verified_source_cif / crystal_edited；validate_formal 通过）
  - test_b20_motifs_not_degraded: OK（B17–B20 motif 保留；Co3O4/CoN4 指纹不变；constructed 不升级；register(overwrite) 保留 provenance）

单独运行：#79–#82 → Ran 4 tests in 7.466s — OK
& $PY run_ci.py --stress
  ok=true
  contract_tests=82
  stress_tests=4
  退出码=0（CI PASS，152.8s）
```

### 2.4 零退化

```text
LEPIDO r_factor=0.027281059410161693
LEPIDO N=6.348637726145664
LEPIDO R=2.0019558444559276
fit_table SHA256=E098D0EDC24E49296710E44518F8A8F010DCD03D7D7EB8CB2E26831A7532E760
#60 静态 s02：OK（5 个已标定元素保持数值；12 个未标定元素仍为 s02=None、calibrated=False）
诊断标记：diagnostic_only=true、cn_status=diagnostic_only
```

### 2.5 正式来源门禁

```text
生成 CIF（_structure_lib/*/_cif/*.cif）→ generated_only → validate_formal=unverified_source
仅有 COD/AMCSD 编号、无源文件 → source_ref_only → validate_formal=unverified_source
用户真实 CIF + sha256 + 吸收位 → verified_source_cif → validate_formal=ok
伪造生成 CIF 为 verified_source_cif → generated_cif_not_formal
```

### 2.6 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：motif 库仍无真实源 CIF；P1a 通过 5 条 COD 金属种子包补 Fe/Co/Ni/Cu/Pt 的 verified 覆盖。
- 下一步：P1a（覆盖基线 + phase pack）

---

## 3. P1a：覆盖基线与 phase pack

### 3.1 改动

- 新增文件：`phase_pack.py`（`PhasePackError` / `parse_cif` / `import_cif` / `rebuild_manifest` / `validate_entry` / `build_index` / `build_coverage` / `verified_entries` / `COMMON_PHASES`）、`tools/import_phase_pack.py`、`tools/build_phase_index.py`
- 新增数据：`_structure_lib/_phase_pack/`（`manifest.json` + 5 条 COD CIF + 5 个 meta）、`_structure_lib/phase_index.json`、`_structure_lib/phase_coverage.json`
- 修改文件：`regression_test.py`（#83–#86；导入 `phase_pack`；契约数 86）、`.gitattributes`（`_structure_lib/_phase_pack/** -text`，保证 CIF 逐字节 sha256 不被换行转换破坏）
- 仓库外同步：`phase_pack_manifest.schema.json` 新增 `absorber_elements`、`original_filename`；`PHASE_PACK_README.md` 同步索引路径与导入要求
- 新增接口：`phase_pack.import_cif()` / `validate_entry()` / `build_index()` / `build_coverage()`；CLI `tools/import_phase_pack.py`、`tools/build_phase_index.py`
- 回滚 commit：`d34a44e P0b: classify motif provenance and reject generated CIFs`
- P1a commit：`b13f9db P1a: add phase pack, coverage audit and real CIF seeds`
- 补丁：`patches/P1a.patch`（SHA256=`CCC43CCE458EF38EB305729F6E9572903C99C17D91343915EAD5B7C57EE17FD5`；`git apply --check --reverse` 通过）

### 3.2 覆盖审计

```text
motif 总数=32
覆盖吸收元素=9/17（Fe、Co、Ni、Cu、Pt、Cr、Zn、Mo、Au）
完全缺失元素=8（Mn、Ru、Rh、Pd、Ag、Ir、Ce、U）
晶体=23；cluster/mixture=9
现场生成 CIF=21
带来源编号=2（goethite=COD:2211652、lepidocrocite=COD:9002763）
真实源 CIF=0
motif_coverage.json 未分类条数=0
```

### 3.3 phase pack

```text
pack 路径=_structure_lib/_phase_pack
导入条目数=5
verified_source_cif 条目数=5
sha256 冲突/拒绝记录=生产导入 0 拒绝；测试显式覆盖 sha256_mismatch、id_sha_conflict、replace=True
phase_coverage.json 元素数（应为 17）=17
Fe/Co/Ni/Cu/Pt verified 条数=1 / 1 / 1 / 1 / 1
```

| 元素 | source:id | CIF sha256（前 12 位） | 空间群 | 吸收位 |
|---|---|---|---|---|
| Fe | COD:1100108 | `41085e65debf` | Im-3m | Fe1 |
| Co | COD:9008492 | `bd9d57408d90` | P 63/m m c | Co |
| Ni | COD:1534892 | `adb86b38b531` | F m -3 m | Ni1 |
| Cu | COD:4105040 | `39726e1afd93` | F m -3 m | Cu |
| Pt | COD:1011103 | `e44ec4c99848` | F m -3 m | Pt1 |

### 3.4 测试

```text
新增测试：
  - test_motif_coverage_audit: OK（17/17 元素；9 覆盖 + 8 缺失；32 条 motif；未分类 0；motif verified 仍为 0）
  - test_phase_pack_import_and_index: OK（5 条 manifest 逐条校验通过；index/coverage 一致；重复导入字节级幂等；临时 pack 可重新索引）
  - test_phase_pack_sha256_mismatch_rejected: OK（meta 假 sha 拒绝；同 ID 不同 sha 冲突；replace=True 覆盖并删旧文件；篡改 manifest 被 sha256_mismatch 拒绝）
  - test_uncovered_element_no_formal_cn: OK（Mn not_covered、Zn motif_only；motif 模型 validate_formal 拒绝；CN 表 0 值 + no_cif 往返）

& $PY run_ci.py --stress
  ok=true
  contract_tests=86
  stress_tests=4
  退出码=0（CI PASS，146.5s）
& $PY tools/backfill_motif_provenance.py --check
  退出码=0（motifs=32，motif_json_updates=0，reports_changed=0）
```

### 3.5 零退化

```text
LEPIDO r_factor=0.027281059410161693
LEPIDO N=6.348637726145664
LEPIDO R=2.0019558444559276
fit_table SHA256=E098D0EDC24E49296710E44518F8A8F010DCD03D7D7EB8CB2E26831A7532E760
#60 静态 s02：OK（5 个已标定元素 Fe/Co/Ni/Cu/Pt；12 个未标定元素仍为 s02=None、calibrated=False）
诊断标记：diagnostic_only=true、cn_status=diagnostic_only
```

### 3.6 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：phase pack 首批只有 5 条金属；Fe/Co/Ni/Cu/Pt 的氧化物、氢氧化物等常见相仍缺。8 个元素（Mn、Ru、Rh、Pd、Ag、Ir、Ce、U）既无 motif 也无 phase pack 条目，正式 CN 只能等用户 `--phase-cif`。motif 库仍无真实源 CIF。
- 下一步：P1b（`phase_resolver.py` + cache），把 `--phase-cif` / `--phase-id` / `--phase` / `--offline` 接到正式 CN；预期 `contract_tests=93`

---

## 4. P1b：相解析与缓存

### 4.1 改动

- 新增文件：`phase_resolver.py`（`PhaseResolutionError` / `PhaseAmbiguousError` / `PhaseOfflineMiss` / `PhaseCandidate` / `PhaseResolution` / `parse_phase_spec` / `resolve_phase` / 缓存读写）、`_structure_lib/_phase_cache/`（P1b 缓存根目录，由测试临时创建并清理）
- 修改文件：`run_xafs.py`（`process_sample` 新增 `phase_cif/phase_id/phase_spec/phase_index_path/phase_pack_path/absorber_site/offline` 参数；在 `xafs_fit.run` 前解析 phase 得到 `StructureModel`；`main()` 新增 `--phase-cif/--phase-id/--phase/--phase-index/--phase-pack/--absorber-site/--offline` CLI 参数）、`regression_test.py`（#87–#93；导入 `phase_resolver`；契约数 93）
- 新增接口：`phase_resolver.resolve_phase(...)`；CLI `--phase-cif`（可重复）/`--phase-id SOURCE:ID`/`--phase SPEC`/`--phase-index`/`--phase-pack`/`--absorber-site`/`--offline`
- 回滚 commit：`b13f9db P1a: add phase pack, coverage audit and real CIF seeds`
- P1b commit：`591447e P1b: phase resolver + cache for formal CN structure sources`
- 补丁：`patches/P1b.patch`（SHA256=`41CE2DBDE75BCDA4BEF52AFC63BEE4C613EF109E0D2D95136F34BD7629803757`；`git apply --check --reverse` 通过）

### 4.2 测试

```text
新增测试：
  - test_phase_spec_parser_no_eval: OK（SPEC 只做 split(";")/split("=",1)；非法 key/格式拒绝；代码注入字符串只当普通字符串，不执行）
  - test_phase_cif_user_source: OK（用户真实 CIF → verified_source_cif；source_type=user_cif；validate_formal 通过）
  - test_phase_id_cache_hit: OK（COD:1100108 缓存命中；幂等；缓存文件正确写入）
  - test_phase_id_cache_miss_offline: OK（未知 ID + 离线 → offline_miss，不退回 motif，不产生 CN 行）
  - test_phase_name_ambiguous: OK（同 pack 下 2 条含同一吸收元素 → ambiguous，不自动 Top-1，selected=None）
  - test_phase_index_schema: OK（schema_version/pack_id/pack_root/manifest_sha256/summary；5 条 verified；sha256 与文件一致）
  - test_cache_sha256_reuse: OK（首次写入缓存；篡改 CIF 字节后复用时 sha256 不匹配，回退 phase pack 重新获取）

& $PY run_ci.py --stress
  ok=true
  contract_tests=93
  stress_tests=4
  退出码=0（CI PASS，150.2s）
```

### 4.3 零退化

```text
LEPIDO r_factor=0.027281059410161693
LEPIDO N=6.348637726145664
LEPIDO R=2.0019558444559276
fit_table SHA256=E098D0EDC24E49296710E44518F8A8F010DCD03D7D7EB8CB2E26831A7532E760
#60 静态 s02：OK（5 个已标定元素 Fe/Co/Ni/Cu/Pt；12 个未标定元素仍为 s02=None、calibrated=False）
诊断标记：diagnostic_only=true、cn_status=diagnostic_only
```

### 4.4 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：phase_resolver 当前只产出 StructureModel 供 `--formal-cn` 门禁使用；P2a 仍需把 CIF 门禁、N₀ 与固定 N 一致性接入正式 CN 行生成。8 个元素（Mn、Ru、Rh、Pd、Ag、Ir、Ce、U）仍无 phase pack 条目，正式 CN 只能等用户 `--phase-cif`。
- 下一步：P2a（CIF 门禁 + N₀ + 固定 N 一致性），预期 `contract_tests=101`

---

## 5. P2a：CIF 门禁 + N₀

### 5.1 改动

- 新增文件：`cif_gate.py`（`CIFGateError` / `CIFGateResult` / `gate_cif` / `gate_structure_source` / Rayleigh CN shell 分组 / `_rayleigh_delta_r` / `_group_shells_rayleigh`）
- 修改文件：`regression_test.py`（#94–#101；导入 `cif_gate`、`math`；契约数 101）
- 新增接口：`cif_gate.gate_cif(cif_path, absorber, *, absorber_site, kmin, kmax, rmax, allow_computed, allow_unusual_coordination, is_computed) -> CIFGateResult`
- 回滚 commit：`591447e P1b: phase resolver + cache for formal CN structure sources`
- P2a commit：`db0c8de P2a: CIF gate + N0 with Rayleigh shell grouping`
- 补丁：`patches/P2a.patch`（SHA256=`F8F55F706DB5EE4834459ED9E20CF58765A596614E5B1CC93844287EB5988BDB`；`git apply --check --reverse` 通过）

### 5.2 门禁检查结果

```text
测试用 CIF（phase pack 5 条金属种子）：
  Fe COD:1100108 (bcc, CN=8):   ok, occ_src=explicit, first_cn=8.0
  Co COD:9008492 (hcp, CN=12):  ok, occ_src=implicit_1.0, first_cn=12.0
  Ni COD:1534892 (fcc, CN=12):  ok, occ_src=explicit, first_cn=12.0
  Cu COD:4105040 (fcc, CN=12):  ok, occ_src=explicit, first_cn=12.0
  Pt COD:1011103 (fcc, CN=12):  ok, occ_src=explicit, first_cn=12.0

拒绝路径：
  parse_failed: 缺 _atom_site_* 循环的 CIF → 拒绝
  absorber_missing: CIF 中无吸收元素 → 拒绝
  absorber_site_ambiguous: 多元素结构中同元素多晶位，未 --absorber-site → 拒绝
  occupancy_missing: 多元素结构缺 _atom_site_occupancy → 拒绝
  cocupancy_unresolved: 同一坐标多元素共占位，未拆分 → 拒绝

允许路径：
  单元素有序结构缺 occupancy → implicit_1.0，放行
  fcc/hcp CN=12 金属 → 正常通过（不报警）
  bcc Fe CN=8 → 正常通过（在 2–12 范围内）
```

### 5.3 N₀ 与 Rayleigh shell 分组

```text
ΔR_Rayleigh = π/(2·(kmax−kmin))：
  k=3–6:   ΔR≈0.524 Å，第一壳 2.484 与第二壳 2.868 间距 0.384 < 0.524，合并（first_cn=14）
  k=3–12:  ΔR≈0.175 Å，不合并（first_cn=8，bcc Fe 第一壳 CN=8）
  k=8–14:  ΔR≈0.262 Å，不合并（first_cn=8）

N₀ = Σ(简并度 × 占据度) 验证通过：
  每个 shell 的 n0 等于 n0_members 中 occ 之和
```

### 5.4 测试

```text
新增测试：
  - test_cif_gate_parse_failed: OK（缺 _atom_site_* → parse_failed）
  - test_cif_gate_absorber_missing: OK（无吸收元素 → absorber_missing）
  - test_cif_gate_absorber_site_ambiguous: OK（多晶位未指定 → absorber_site_ambiguous；指定后放行）
  - test_cif_gate_occupancy_implicit_ok: OK（COD:9008492 Co，缺 occupancy 但单元素有序 → implicit_1.0，CN=12）
  - test_cif_gate_occupancy_disorder_rejected: OK（多元素缺 occupancy → occupancy_missing）
  - test_cif_gate_cocupancy_unresolved: OK（共占位未拆分 → cocupancy_unresolved）
  - test_cn_range_12_metals: OK（Fe/Co/Ni/Cu/Pt 各金属 CN 正常通过）
  - test_n0_shell_grouping_rayleigh: OK（Rayleigh delta_r 计算正确；窄 k 范围合并更多；N₀=Σocc）

& $PY run_ci.py --stress
  ok=true
  contract_tests=101
  stress_tests=4
  退出码=0（CI PASS，185.8s）
```

### 5.5 零退化

```text
LEPIDO r_factor=0.027281059410161693
LEPIDO N=6.348637726145664
LEPIDO R=2.0019558444559276
fit_table SHA256=E098D0EDC24E49296710E44518F8A8F010DCD03D7D7EB8CB2E26831A7532E760
#60 静态 s02：OK（5 个已标定元素 Fe/Co/Ni/Cu/Pt；12 个未标定元素仍为 s02=None、calibrated=False）
诊断标记：diagnostic_only=true、cn_status=diagnostic_only
```

### 5.6 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：cif_gate 目前只做“拒绝或放行”判断，不生成 cn_table.csv 行（正式 CN 行生成留到 P2b/P3）。计算结构检测依赖调用方传入 `is_computed` 参数（phase pack entry 的 `experimental` 字段可用于此）。键长合理性检查对未知元素对只警告不硬失败（需用户提供参考键长）。
- 下一步：P2b（固定 N 一致性 + 自由 N 对照），预期 `contract_tests=103`

## 6. P2b：固定 N 一致性

### 6.1 改动

- 修改文件：`xafs_fit.py`（新增 `map_n0_to_plan(plan, n0_shells, ...)`：N₀ shell → FEFF 路径映射，返回 n0_per_path / path_links / n_var_fixed_n / n_var_free_n_first；`select_shell_paths` 新增 `n_mode` 参数：`free`=3n+1、`fixed_n`=2n+1、`conditional`=2n+2）、`regression_test.py`（#102–#103；契约数 103）
- 新增接口：`xafs_fit.map_n0_to_plan()`、`xafs_fit.select_shell_paths(..., n_mode=...)`
- 回滚 commit：`db0c8de P2a: CIF gate + N0 with Rayleigh shell grouping`
- P2b commit：`ac24373 P2b: fixed-N consistency + free-N comparison`
- 补丁：`patches/P2b.patch`（SHA256=`894FB634C7CAA4744CC63B4BECF6A81CD35386943AB9B1C62450D927784A62C7`；`git apply --check --reverse` 通过）

### 6.2 固定 N 信息量公式

```text
select_shell_paths(n_mode=...) 的 n_var：
  free:      3n+1  （每壳 N/ΔR/σ²，加全局 ΔE₀）
  fixed_n:   2n+1  （每壳 ΔR/σ²，加全局 ΔE₀；N 固定为 N₀）
  conditional: 2n+2 （首壳多放开一个 N，其余固定）

map_n0_to_plan 返回：
  n_var_fixed_n   = 2n+1
  n_var_free_n_first = 2n+2
  n0_per_path     = 每条路径的 N₀（用于 path_degen）
  path_links      = 同 CN 组子壳的 degen_ratio 绑定（只放开参考路径 N）

验证（Fe bcc，n_paths=2）：
  fixed_n:   n_var = 5 = 2×2+1 ✅
  conditional: n_var = 6 = 2×2+2 ✅
  n0_per_path = [8.0, 6.0]（与 CIFGateResult n0_shells 前两壳一致）
  path_links  = [{target:1, source:0, degen_ratio:0.75}]（第二壳 N₀/第一壳 N₀）
```

### 6.3 测试

```text
新增测试：
  - test_fixed_n_nidp_formula: OK（fixed_n=2n+1，conditional=2n+2，map_n0_to_plan n_var 一致，n0_per_path 与前 N 个壳层 N₀ 一致）
  - test_wrong_phase_consistency_failed: OK（正确相 → mapping_ok=True；错误相元素不匹配 → mapping_ok=False，n0_per_path 为空，相当于 0 CN）

& $PY run_ci.py --stress
  ok=true
  contract_tests=103
  stress_tests=4
  退出码=0（CI PASS，160.0s）
```

### 6.4 零退化

```text
LEPIDO r_factor=0.027281059410161693
LEPIDO N=6.348637726145664
LEPIDO R=2.0019558444559276
fit_table SHA256=E098D0EDC24E49296710E44518F8A8F010DCD03D7D7EB8CB2E26831A7532E760
#60 静态 s02：OK（5 个已标定元素 Fe/Co/Ni/Cu/Pt；12 个未标定元素仍为 s02=None、calibrated=False）
诊断标记：diagnostic_only=true、cn_status=diagnostic_only
```

### 6.5 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：P2b 当前只实现了 N₀→路径映射和信息量公式，未接入实际的 artemis.fit_r_space 固定 N 拟合和 Hamilton F-test（留到 P3 或后续版本）。map_n0_to_plan 的 path_links 仅处理同元素相近壳层的 degen_ratio 绑定，未处理异核共占位分支绑定。
- 下一步：P3（条件 CN、多相与缺陷模型），预期 `contract_tests=110`

## 7. P3：条件 CN

### 7.1 改动

- 修改文件：`xafs_fit.py`（新增 `build_conditional_fit_args()`：首壳单一比例因子，仅放开参考路径 N；`compute_n_ratio()`：计算 N_fit_total/N0_total；`check_n_ratio()`：0.7–1.3 ok / 0.6–1.4 warning / 超出 out_of_range；`parse_defect_model()`：nanoparticle/vacancy/custom 结构化 JSON 解析；`check_defect_model()`：N/N₀ 与模型预期比对）、`regression_test.py`（#104–#110；契约数 110）
- 新增接口：`xafs_fit.build_conditional_fit_args()`、`xafs_fit.compute_n_ratio()`、`xafs_fit.check_n_ratio()`、`xafs_fit.parse_defect_model()`、`xafs_fit.check_defect_model()`
- 回滚 commit：`ac24373 P2b: fixed-N consistency + free-N comparison`
- P3 commit：`5900069 P3: conditional CN + defect models + multiphase no-auto-pick`
- 补丁：`patches/P3.patch`（SHA256=`8D3AF4855C13C4D97E02231324FBDBC2F482B2111DEF58E8B083E82B57C79A3A`；`git apply --check --reverse` 通过）

### 7.2 首壳单一比例因子

```text
build_conditional_fit_args(plan, n0_mapping):
  vary_n = [ref_idx]           # 仅放开首壳参考路径 N（N₀ 最大的首壳路径）
  path_degen = n0_per_path     # 所有路径 N 固定为 N₀ 初值
  path_links = 非参考首壳路径 → 参考路径，degen_ratio = N0_j/N0_ref
               share_deltar=False, share_sigma2=False
  n0_total = Σ 首壳 N₀
```

### 7.3 N/N₀ 范围检查

```text
check_n_ratio(n_ratio):
  0.7–1.3:   ok（通过）
  0.6–0.7 / 1.3–1.4: warning
  <0.6 或 >1.4: out_of_range（需定量缺陷模型）
```

### 7.4 缺陷模型

```text
parse_defect_model() 支持：
  nanoparticle: {"type": "nanoparticle", "diameter_nm": 2.5, "lattice": "fcc"}
    → expected_n_ratio = 1 − 1.5u + 0.5u³（u = r1/d）
    → 仅 fcc/hcp/bcc 晶格
  vacancy:    {"type": "vacancy", "site": "O1", "fraction": 0.1}
    → expected_n_ratio = 1 − fraction
  custom:     {"type": "custom", "expected_n_ratio": 0.75, "uncertainty": 0.10}

check_defect_model(n_ratio, defect_model):
  |n_ratio − expected| ≤ uncertainty → defect_ok
  否则 → defect_model_mismatch
  无模型 → defect_model_missing
```

### 7.5 测试

```text
新增测试：
  - test_first_shell_single_scale_factor: OK（vary_n 仅含 ref_idx；path_degen=n0_per_path；path_links degen_ratio 正确；n0_total 一致）
  - test_conditional_cn_out_of_range_requires_defect: OK（0.5/1.5 → out_of_range；无模型 → defect_model_missing；模型不匹配 → defect_model_mismatch；0.9/1.2 → ok）
  - test_conditional_cn_nanoparticle_model: OK（nanoparticle 解析正确；ratio 在预期内 → defect_ok；偏离 → mismatch；非 fcc/hcp/bcc 拒绝）
  - test_multiphase_no_auto_pick: OK（2 条含 Fe 相候选 → ambiguous，selected=None，不自动 Top-1）
  - test_cn_table_provenance_complete: OK（正式 CN 行含 model_id/absorber/absorber_site/model_level/source_type/source_ref/source_cif_sha256/provenance_status 等全部字段）
  - test_prototype_cluster_separate: OK（prototype 和 cluster 模型 validate_formal 均拒绝，不进正式 CN 表）
  - test_reproducible_cn_hash: OK（cn_table.csv 写读往返，sha256 逐字节一致）

& $PY run_ci.py --stress
  ok=true
  contract_tests=110
  stress_tests=4
  退出码=0（CI PASS，146.9s）
```

### 7.6 零退化

```text
LEPIDO r_factor=0.027281059410161693
LEPIDO N=6.348637726145664
LEPIDO R=2.0019558444559276
fit_table SHA256=E098D0EDC24E49296710E44518F8A8F010DCD03D7D7EB8CB2E26831A7532E760
#60 静态 s02：OK（5 个已标定元素 Fe/Co/Ni/Cu/Pt；12 个未标定元素仍为 s02=None、calibrated=False）
诊断标记：diagnostic_only=true、cn_status=diagnostic_only
```

### 7.7 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：P3 当前实现了条件 CN 的参数构建和 N/N₀ 检查，但未接入实际的 artemis.fit_r_space 固定 N 拟合 + 自由 N 对照 + Hamilton F-test（需要真实谱数据，留到后续版本）。多相不自动选已在 phase_resolver 的 ambiguous 状态中体现。cn_table 分表（cn_table.csv / cn_model_dependent.csv）的逻辑已在 P0a 中定义，P3 未修改。
- 下一步：P4（验证与迁移），预期 `contract_tests≥114`

## 8. P4：验证与迁移

### 8.1 改动

- 修改文件：`run_xafs.py`（`process_sample` 新增 `legacy_auto_fit` 参数；`main()` 新增 `--legacy-auto-fit` CLI 参数；formal_cn 处理时 legacy_auto_fit=True 标记 `legacy_unverified` 状态，不进正式表）、`regression_test.py`（#111–#114；契约数 114）
- 新增接口：CLI `--legacy-auto-fit`
- 回滚 commit：`5900069 P3: conditional CN + defect models + multiphase no-auto-pick`
- P4 commit：`21849bf P4: verification, migration, legacy-auto-fit marking`
- 补丁：`patches/P4.patch`（SHA256=`5E7E4EC50241F6BCD6250B0761C618B5C5671AED2D5ED43FF95909B709FC2EA1`；`git apply --check --reverse` 通过）

### 8.2 离线全链与可复现性

```text
test_offline_full_chain: OK（同一 CIF 两次 gate_cif 调用，cif_sha256/first_shell_cn/occupancy_source 完全一致；结果 dict 哈希一致）
```

### 8.3 标准箔 formal CN

```text
test_standard_foils_formal_cn: OK（Fe/Co/Ni/Cu 四种标准箔 CIF，gate_cif 放行；occupancy_source 为 explicit 或 implicit_1.0；StructureModel + validate_formal 均通过）
```

### 8.4 无 CIF 报告空 CN

```text
test_no_cif_report_empty: OK（cn_table.refusal_row("no_cif", ...) → gate_status=no_cif, n_fit="", cn_label="none", flags 包含 refused；写读往返后 n_fit 仍为空）
```

### 8.5 legacy-auto-fit 标记

```text
test_legacy_auto_fit_marked: OK（legacy_auto_fit=True → formal 状态=legacy_unverified，reject_reason 包含 "legacy_auto_fit"，legacy_auto_fit 标志=True）
```

### 8.6 测试

```text
新增测试：
  - test_offline_full_chain: OK（离线全链确定性）
  - test_standard_foils_formal_cn: OK（标准箔 formal CN 放行）
  - test_no_cif_report_empty: OK（无 CIF → 拒绝行，n_fit 为空）
  - test_legacy_auto_fit_marked: OK（legacy_auto_fit 标记 legacy_unverified）

& $PY run_ci.py --stress
  ok=true
  contract_tests=114
  stress_tests=4
  退出码=0（CI PASS，139.5s）
```

### 8.7 零退化

```text
LEPIDO r_factor=0.027281059410161693
LEPIDO N=6.348637726145664
LEPIDO R=2.0019558444559276
fit_table SHA256=E098D0EDC24E49296710E44518F8A8F010DCD03D7D7EB8CB2E26831A7532E760
#60 静态 s02：OK（5 个已标定元素 Fe/Co/Ni/Cu/Pt；12 个未标定元素仍为 s02=None、calibrated=False）
诊断标记：diagnostic_only=true、cn_status=diagnostic_only
```

### 8.8 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：P4 是本修订方案的最后一批。所有 8 个阶段（P0a→P4）已完毕，累计契约回归 114 项，压力测试 4 项，零退化。后续如需 P4+（真实氧化物 CIF、错误相一致性拟合、多相完整流程），可在本基线上继续扩展。
- 下一步：无（P0a–P4 全部完成）

---

## 9. P4+：真实 CIF 覆盖、构建期爬取整合与条件 CN 正式行

### 9.1 改动

- 新增 `tools/cif_crawler.py`：构建期 CIF 抓取/核验/phase pack 导入工具。
  - 固定 21 个 COD 目标，不做化学式自由检索；
  - 只允许 `--online` 联网，默认只读 `data/cif_raw/` 缓存；
  - 下载后必须通过 gemmi 解析、`_chemical_formula_sum` 化学计量、IT 空间群号、吸收元素位点核验；
  - 现场生成 CIF 拒绝导入；
  - 不写 motif、不修改 `_structure_lib/<元素>/` 诊断库。
- phase pack 从 5 条金属种子扩展到 **26 条**真实 COD CIF（5 种子 + 21 抓取目标），覆盖 **17/17** 吸收元素。
- `phase_resolver.py`：多吸收位相（Fe3O4/Co3O4/α-Mn）必须显式 `--absorber-site`；缺失/错误位点返回明确 reason，不静默取第一个位点。
- `xafs_fit.py`：新增 `build_formal_cn_block()`，把条件 CN 结果转为正式 CN 行；新增 `cn_table.data_sha256()`，正式行带输入数据哈希；修正 `map_n0_to_plan()`，不同 N₀ 壳层不再互相绑定，远壳严格固定 N₀。
- `run_xafs.py` / `xafs_report.py`：记录输入数据哈希，正式表使用条件 CN 行。
- `phase_pack.py`：`_phase_tag_match` 改为化学式精确 + tag 精确 + 长度≥4 相名子串；生成 JSON 统一 LF。
- 契约新增 #115–#119；#84/#86/#92 更新为 26 条 / 17 元素 / 21 目标。
- 提交：`d07adde P4+ / CIF crawler: build-time COD motifs + T5 standards, migrate legacy IDs`
- 补丁：`patches/P4+.patch`（SHA256=`40A353C810C4037827D43ACE96810AF33B4430CE1BEF57F50908459EE62D6692`；`git apply --check --reverse` 通过）

### 9.2 覆盖结果

```text
phase pack entries: 26
verified_cif: 26
elements_covered: 17/17
verified_elements: Fe,Co,Ni,Cu,Pt,Cr,Mn,Zn,Mo,Ru,Rh,Pd,Ag,Ir,Au,Ce,U
```

T3 常见氧化物/氢氧化物 9 条：Fe2O3、Fe3O4、CoO、Co3O4、NiO、Ni(OH)2、CuO、Cu2O、CeO2。  
T4 金属基准相 12 条：Mn、Ru、Rh、Pd、Ag、Ir、Ce、U、Cr、Zn、Mo、Au。

说明：以上 26 条是 `tools/cif_crawler.py` 的 phase pack 条目；根目录 `cif_crawler.py` 另建 17 条 COD motif（29 条索引 = 17 CIF + 12 STD），见 9.7。

### 9.3 phase pack CIF crawler 验收

```text
python tools/cif_crawler.py --list
  21 个目标，COD ID / 空间群 / 相名齐全

python tools/cif_crawler.py --all --offline --dry-run --no-index
  17/17 缓存核验通过

python tools/cif_crawler.py --all --offline
  26 条 phase pack，17/17 元素覆盖；二次运行 already_present，幂等

错相测试：Fe3O4 CIF 冒充 CeO2 → composition_mismatch，拒绝导入
离线缺缓存：offline_miss，不退回 motif
```

### 9.4 真实 Fe 箔端到端正式 CN

```text
run_xafs.py data/练手数据/01_Fe_foil/fe.060 \
  --absorber Fe --fit --formal-cn \
  --phase-id COD:1100108 --absorber-site Fe1 --offline --no-strict
```

关键输出：

```text
[P1b] phase resolved: COD:1100108 sha256=41085e65debf... site=Fe1
[P2a] CIF gate passed: occ_src=explicit, first_shell_cn=8.0, n0_shells=6
[P3-1] 固定 N 拟合：n_paths=2, vary_n=[]
[P3-2] 条件 CN 拟合：vary_n=[0], path_links=0
[P3-3] 自由 N 对照拟合：vary_n=True
[P3-4] N/N₀ = 0.7261 (ok), F-test significant=False, consistency_passed=True
```

`cn_table.csv` 产生 2 行正式 CN：

| shell | n_status | n0 | n_fit | n_ratio | source_cif_sha256 |
|---|---|---|---|---|---|
| Fe–Fe@2.48 | fitted | 8.0 | 5.8086 | 0.7261 | 41085e65debf... |
| Fe–Fe@2.87 | fixed | 6.0 | 6.0 | 1.0 | 41085e65debf... |

每行均含 `phase_id / model_level / source_type / source_ref / source_cif_sha256 / absorber_site / n0 / n_fit / n_status / n_ratio / s02 / s02_source / s02_calibrated / data_sha256 / software_version`。  
数据哈希：`a16245c385c25145abc36862e3f64f036d1adda46e4be4c07d8cf11be2b6b678`。

### 9.5 测试

```text
& $PY run_ci.py --stress
  ok=true
  contract_tests=121
  stress_tests=4
  总耗时=154.7s
  退出码=0（CI PASS）
```

新增契约：

- #115 21 个抓取目标 ID 唯一、元素/化学式/空间群自洽；
- #116 离线缓存 → 核验 → phase pack 导入，sha256 可追溯；
- #117 离线缺缓存 → offline_miss；错相 CIF → 拒绝导入；
- #118 Fe3O4 多吸收位必须显式 `--absorber-site`；
- #119 条件 CN 通过 → 正式行完整可追溯；一致性失败 → 0 行；
- #120 T7 多元素自动分组、失败隔离、并列报告、`recommended=null`；
- #121 T8 异常谱警告、干净谱无警告、`cn_reliable` 不被修改。

### 9.6 零退化

```text
LEPIDO r_factor=0.027281059410161693
LEPIDO N=6.348637726145664
LEPIDO R=2.0019558444559276
fit_table SHA256=E098D0EDC24E49296710E44518F8A8F010DCD03D7D7EB8CB2E26831A7532E760
#60 静态 s02：OK（5 个已标定元素 Fe/Co/Ni/Cu/Pt；12 个未标定元素仍为 s02=None、calibrated=False）
```

crawler 入库后按同一三扫描命令重跑，上述 r_factor/N/R 与 fit_table SHA256 逐位一致。  
注：`regression_test` 的 LEPIDO 片段只喂入 `.000` 单扫描（Rf≈0.0270、N≈6.41），与三扫描基线不是同一输入，不是退化。

### 9.7 根目录 CIF crawler 与 T5 标准谱

- 根目录 `cif_crawler.py` 作为构建期权威工具并入；`tools/cif_crawler.py` 继续只负责 phase pack 导入，两套资产独立，均不在 S1–S7 运行路径上。
- 构建命令与结果：
  - `python cif_crawler.py --all --import-lib` → 17 条 COD motif；
  - `python cif_crawler.py --t5` → 12 条 IIT-Lytle 真实金属箔谱；
  - `data/cif_crawler_index.json` 共 **29 条**：17 CIF + 12 STD；
  - `T5_PLACEHOLDERS 0`、`LEGACY_MOTIFS_LEFT 0`、`QUARANTINE_FILES 0`。
- 结构库：**43 个 motif**、**17/17** 元素覆盖、0 missing；`provenance_report` 为 `constructed=9 / generated_only=32 / source_ref_only=2`，`verified_source_cif=0`。
- 旧 ID 迁移：
  - `nio_rocksalt` → `nio_bunsenite_cod_1010093`；
  - `cu2o_cuprite` → `cu2o_cuprite_cod_1010941`；
  - 同批 `coo_rocksalt`、`co3o4_tet`、`ni_oh2_brucite`、`cuo_tenorite` 迁移到 COD 命名；
  - `co3o4_mixed` 的 tet 组件改为 `co3o4_spinel_cod_9005888`，保留 1/3 tet + 2/3 oct 语义。
- 多吸收位：Fe3O4/Co3O4/α-Mn 通过显式 `absorber_fract` 入库；未合入 `EnhancedMotifBuilder`，无 `meta`/`center_site`/BVS 改动。
- T5 边界（当批）：12 条真实谱写入 `_standards_lib/<元素>/*.dat`，`data/standards_raw/` 保存抓取缓存；当批尚未接入标定，12 元素仍 `s02=None`、`calibrated=False`；后续 §9.8 已完成独立标定产物（静态表仍不变）。
- 运行期边界：只通过 `structure_lib.load_library(absorber)` 读取 `_structure_lib`；不直接读 `data/cif_crawler_index.json` 取几何；`run_xafs.py` 不导入根爬虫。
- 数据提交：`data/cif_raw` 17 个被索引 CIF、`data/cif_crawler_index.json`、`_structure_lib/**/*.json`、`_cif`、`_phase_pack` 随仓库提交；IIT-Lytle 许可未明确，`data/standards_raw/` 与 `_standards_lib/` 保留为本地构建产物并用 `--t5` 重建。
- 回归（当批）：`regression_test.py` 119 项 OK；`run_ci.py --stress` PASS（119 + 4，170.5s）；NiO 压力测试 top1=`nio_bunsenite_cod_1010093`，`ni_fcc` 垫底；LEPIDO 三扫描基线 SHA256 不变（见 9.6）。后续 §9.8 批次为 121 + 4、154.7s。
- 结论：本批次完成；T5 S₀² 标定、常见相补全、多位点联合模型列为后续任务（T5 独立标定与 T7/T8 已在 §9.8 完成）。


### 9.8 T5 S₀² 标定、T7 多元素、T8 XANES–EXAFS 联合校验（2026-11-09）

#### T5 独立标定产物（#60 静态表不变）

- `calibrate_s02.py` 新增 `--motif-file` / `--rmin`；`--apply` 对 #60 锁定的
  12 元素硬拒绝（静态 `s02` 必须保持 `None`）。
- 标定专用 motif 落 `data/s02_calibration/motifs/`：`cr2o3_corundum`、
  `mno2_rutile`、`zno_wurtzite`、`ruo2_rutile`、`moo2_rutile`、`mos2_2h`、
  `moo4_tetra`；JSON 带展开 `sites`，不写 `_structure_lib`，不影响 43 条
  生产 motif 的覆盖率断言。
- `calibrate_t5.py` 批量跑 12 元素，产物：
  - `data/s02_calibration/index.json`（`pass_count=6`、`fail_count=6`）；
  - `<元素>/calibration.json`（完整 runs、门禁、标准谱/motif SHA256）；
  - `data/s02_calibration/REPORT.md`。
- 通过 6/12：Cr 0.6482、Mn 0.4311、Zn 0.6071、Ru 0.7002、
  Ag 0.7513、Au 0.8602。
- 未通过 6/12：
  - Mo：Na₂MoO₄ 四面体标准 s02≈0.38 低于 0.4 下界；Tugarinovite
    s02_err>0.06；MoS₂ 仅单窗口内部点；Mo 箔撞 1.3 上界；
  - Rh/Pd/Ir：IIT-Lytle 谱仅约 170–195 eV above E0，S5/S6 WT 无有效网格；
  - Ce：rmin=1.6 下前两窗口 ΔE₀>10 eV，三窗口门禁未通过；
  - U：α-U 拟合 s02 撞 0.3 下界、Rf 差。
- #60：`absorber.py` 未改；12 元素仍 `s02=None`、`calibrated=False`；
  运行时用 `run_xafs.py --s02` 显式覆盖。
- 许可边界：`_standards_lib/` 与 `data/standards_raw/` 不提交；Ag 标定依赖
  本地 IIT-Lytle 谱，产物记录源文件 SHA256，但不保证无本地标准时可复现。

#### T7 多元素联合判据

- `run_xafs.py` 新增 `--multi-element` 与 `--multi-element-s02
  'Fe=0.803,Co=0.7308'`；自动从 XDI `Element.symbol` / 纯文本
  `Absorber:` 头判定吸收元素；每个元素独立子目录、独立 `--s02`。
- 输出 `multi_element_summary.json` 与 `multi_element_summary.md`；
  `recommended=null`、`best_element=null`，各元素结果并列，不自动 Top-1。
- #120 回归：自动分组、单元素失败隔离、并列报告、无推荐字段、
  `absorber_source` 溯源。

#### T8 XANES–EXAFS 联合校验

- `xafs_fit.xanes_exafs_cross_check()`：归一化平台、原始边步、E0 导数双峰、
  白线拟合可靠性/位置、K 边极端 CN–白线对比、荧光自吸收风险。
- `run_xafs.py --xanes-cross-check` 可选开关（可独立使用，也可与
  `--identify`/`--fit` 同用）；结果写入 `<stem>_meta.json` 的
  `xanes_cross_check`；`affects_cn_reliable=False`，不改变正式 CN 门禁。
- #121 回归：异常谱产生 warning、干净谱无 warning、L3 边不做 CN–白线
  极端对比、`cn_reliable` 不被修改。

#### 提交与补丁

- 提交：`d7fa8c8`（P4+ / T5-T8: independent S02 calibration,
  multi-element batch, XANES-EXAFS cross-check）。
- `patches/P4+.patch` SHA256：`9B9779EBD39D3C493EF3D877797ADB9E1E0D33EA1982DAD302F4239524B5FF10`；
  `git apply --check --reverse` 通过。
- `#60` 复核：`absorber.py` 静态表仍为 Fe/Co/Ni/Cu/Pt 已标定，
  其余 12 元素 `s02=None`、`calibrated=False`。

#### 测试

- `regression_test.py` 121 项契约（新增 #120/#121）。
- `run_ci.py --stress`：`ok=true`、`contract_tests=121`、`stress_tests=4`、
  总耗时 154.7s、退出码 0（`ci_report.json/.md` 已更新）。

### 9.9 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：
  - T5 已产出独立标定产物（6/12 通过）；按 #60 静态表仍为 `s02=None`，正式 CN 需运行时 `--s02`；Mo/Rh/Pd/Ir/Ce/U 仍需更高质量的宽能区标准；
  - 常见氧化物/氢氧化物仍缺：FeO、FeOOH、PtO2、Cr2O3、MnO/MnO2/Mn3O4、ZnO、MoO2/MoO3、RuO2、Rh2O3、PdO、Ag2O、IrO2、Au2O3、UO2/U3O8；
  - 多位点联合模型尚未自动构建，当前只能逐位点显式指定；
  - T7/T8 已完成并通过 #120/#121；Mo/Rh/Pd/Ir/Ce/U 的 S₀² 标定仍待补。
- 下一步：补 Mo/Rh/Pd/Ir/Ce/U 的宽能区标准与标定；常见氧化物/氢氧化物补相（继续按固定 COD ID 流程，不自由检索）。

### 9.10 P5-A：T5 工程加固与运行期溯源

#### 9.10.1 改动

- `calibrate_t5.py`：
  - `--only <元素...>` 只重跑所选元素，其余元素从既有
    `<元素>/calibration.json` 读取并合并回 `index.json`/`REPORT.md`，
    不再把汇总覆盖成子集；损坏或元素不符的产物跳过并警告，不伪造。
  - `write_index` 兼容既有 payload 无 `attempts` 字段的情况。
- `run_xafs.py`：
  - 新增 `--s02-from-t5`、`load_t5_s02()`、`_resolve_s02_for_paths()`；
  - 只接受 `status="calibrated"` 且
    `static_absorber_py_unchanged=True` 的逐元素产物；未通过元素抛
    `XAFSConfigError`，错误含“未通过”；
  - 优先级 `--multi-element-s02` > 全局 `--s02` > T5 产物；显式
    `--absorber` 优先于文件头元素判定；
  - 单文件、批处理、operando、多元素分支均传递 `s02_source`；
  - T7 综合报告新增 `S₀²来源` 列。
- `xafs_fit.py`：
  - `fit_identified()`/`run()` 新增 `s02_source`；
  - L1/motif 输出记录 `s02` 与 `s02_source`；
  - 正式 CN 块用 `_s02_calibrated_flag()`：T5 产物来源记
    `s02_calibrated=True`，普通 `--s02` 手填值仍为 `False`。
- `regression_test.py`：新增 #122（T5 产物完整性/静态锁/`--only` 合并）
  与 #123（`load_t5_s02` 与来源优先级），契约总数 121 → 123。
- 文档：`README.md`、`CIF_CRAWLER_INTEGRATION.md`、
  `P4_PLUS_EXECUTION_PLAN.md` 同步 `--s02-from-t5`、溯源串与 123 项账本。

#### 9.10.2 验收

- `--fit` 冒烟（`Cr-K-Cr2O3.xdi --absorber Cr --s02-from-t5 --fit`）：
  `fit.s02=0.6482`，`fit.s02_source=
  s02_calibration:data\s02_calibration\Cr\calibration.json#sha256:b730941b751d;standard_sha256:b68293e5e588`；
  L1 `N=6.34±0.21`、`R=1.975 Å`、`σ²=0.0027`、`Rf=0.0216`。
- `--only Cr` 增量合并：`index.json` 保留既有 Mn/Au，`pass_count=3`；
  #122 覆盖该行为。
- #123：`load_t5_s02("Cr")=0.6482`；`load_t5_s02("Mo")` 抛
  `XAFSConfigError` 且含“未通过”；显式 `--s02` 覆盖 T5 产物。
- `#60` 复核：`absorber.py` 静态表仍为 Fe/Co/Ni/Cu/Pt 已标定，
  其余 12 元素 `s02=None`、`calibrated=False`。

#### 9.10.3 测试

- `regression_test.py`：123 项 OK（180.177s）。
- `run_ci.py --stress`：`ok=true`、`contract_tests=123`、
  `stress_tests=4`、PASS（183.1s）、退出码 0。
- 4 个压力测试（Co/Ni/Cu/NiO）全部通过。

#### 9.10.4 提交与补丁

- 提交：`c856ee9`（P5-A / T5: incremental calibration merge +
  runtime --s02-from-t5 provenance）。
- 补丁：`patches/P5-A.patch`（增量 diff `d7fa8c8..c856ee9`）；
  SHA256=`B81C9DDAF59C57DFE4780561796058573CD7CED57FF23D719582FBCA54D42BAF`；
  `git apply --check --reverse` 通过。
- `patches/P4+.patch` 保持 P4+ 历史补丁不变。

#### 9.10.5 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：
  - T5 通过仍为 6/12；Mo/Rh/Pd/Ir/Ce/U 的物理排查（P5-B）尚未开始，
    门禁不得放宽；
  - P5-C（固定 COD ID 扩 phase pack）与 P5-D（`--absorber-site all`
    多位点正式 CN）暂缓，待 P5-A 文档收尾后启动。
- 下一步：P5-B 先排查 Mo（可选双壳 bcc 模型）与 Ce（E₀ 方法），
  仍按现有能量范围/信噪/稳定性门禁判定。

### 9.11 P5-B Mo 双壳 bcc 与 Ce 物理排查（2026-11-09）

#### 9.11.1 改动

- `calibrate_s02.py`：
  - `calibrate_s02()` 新增 `e0_method`/`manual_e0`/`fit_shells`；
    CLI 新增 `--e0-method`、`--manual-e0`、`--fit-shells`；
  - 多壳拟合共享 S₀²/ΔE₀，各壳独立 ΔR/σ²；`runs` 记录
    `sigma2_at_bound` 与逐路径 `paths`；
  - 新增门禁 `sigma2_not_at_bound`：任一窗口 σ² 压到 1e-4 下界即失败。
    这是更严格的物理约束，不放松任何原有阈值；
  - 产物新增 `n_shells`/`shell_paths`/`e0`/`e0_method`，`n_variables`
    按 `2 + 2×n_shells` 计算。
- `calibrate_t5.py`：Mo 箔 case 改为 `half_height` + `fit_shells=2`；
  Ce 增加第二个 `strict=True` 的透射谱 case；attempt 记录
  E₀ 方法与壳层数；汇总层补写 `n_shells`/`shell_paths`/`e0`/`e0_method`。
- `regression_test.py`：#122 更新为 Mo 已标定、`pass_count>=7`；
  #123 改用 Ce 验证拒绝；新增 #124 锁定 Mo 双壳产物与 Ce 拒绝路径。
  契约总数 123 → 124。
- 文档：`README.md`、`CIF_CRAWLER_INTEGRATION.md`、
  `P4_PLUS_EXECUTION_PLAN.md` 同步 7/12、P5-B 小节与 124 项账本。

#### 9.11.2 验收

- Mo 箔：导数 E₀=20018.8 eV 落在近边肩峰，归一化平台 0.842
  （低于 0.85）；改用 `half_height` E₀=20012.8624 eV 后平台 0.8631。
  bcc 第一壳 8×2.7254 Å + 第二壳 6×3.147 Å 同属 Mo–Mo；单壳模型
  把第二壳尾巴折进第一壳振幅、S₀² 撞 1.3 上界；双壳模型共享
  S₀²/ΔE₀、各壳独立 ΔR/σ²，得到：
  - `S₀²=0.8322±0.0334`、`Rf=0.0040`、窗口极差 `0.0017`、
    `ΔE₀=−0.73 eV`、`nvar=6 < N_idp=6.88`；
  - 三个窗口全部 `interior`，门禁 5/5 通过。
- Mo CLI 冒烟：
  `run_xafs.py <Mo foil> --absorber Mo --e0-method half_height --s02-from-t5 --identify --fit`
  → `fit.s02=0.8322`，
  `fit.s02_source=s02_calibration:data\s02_calibration\Mo\calibration.json#sha256:3096ea30460a;standard_sha256:bfd033f23198`；
  L1 `Mo–Mo@2.73 N=8.03±0.41`、`Mo–Mo@3.15 N=6.31±0.66`、
  `Rf=0.0064`。不带 `--e0-method half_height` 时导数 E₀ 使平台落回
  0.842，被严格归一化门禁拒绝，属预期物理判定。
- Ce：保持 `not_calibrated`，不放宽门禁：
  - 自由 σ² 在三个窗口均塌到 1e-4 下界（`sigma2_at_bound`）；
  - 固定 σ² 扫描：S₀² 随 σ²=0.002→0.008 从 0.538 单调漂到 0.798，
    无唯一解；
  - 加 Ce–Ce 3.82 Å 第二壳后 Rf=0.127–0.130、ΔR 撞界；
  - `half_height` 把 ΔE₀ 推到 12.9–13.5 eV（导数法 9.8–10.2 eV）；
  - 备选 `Ce-L3-CeriumDioxide.xdi` 透射谱平台 2.271，严格归一化
    直接拒绝，属归一化伪迹，不得用于标定。
- `#60` 复核：`absorber.py` 静态表仍为 Fe/Co/Ni/Cu/Pt 已标定，
  其余 12 元素 `s02=None`、`calibrated=False`。
- T5 全量重跑：通过 7/12（Cr/Mn/Zn/Mo/Ru/Ag/Au），失败 5/12
  （Rh/Pd/Ir/Ce/U）。

#### 9.11.3 测试

- `regression_test.py`：124 项 OK（167.989s）。
- `run_ci.py --stress`：`ok=true`、`contract_tests=124`、
  `stress_tests=4`、PASS（202.46s）、退出码 0。
- 4 个压力测试（Co/Ni/Cu/NiO）全部通过。

#### 9.11.4 提交与补丁

- 提交：`d3e178c`（P5-B: Mo two-shell bcc calibration and Ce rejection）。
- 补丁：`patches/P5-B.patch`（增量 diff `c856ee9..d3e178c`）；
  大小 118862 字节；
  SHA256=`A46ADE770327A800B60569C4523B552E9F967CD7BC5D6659C6E20B9980866673`；
  `git apply --check --reverse` 通过。
- `patches/P4+.patch` 与 `patches/P5-A.patch` 保持不变。

#### 9.11.5 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：
  - Rh/Pd/Ir 的 IIT-Lytle 标准谱仅约 170–195 eV above E0，
    k_max 不足；Ce 的 σ² 与 S₀² 强相关、无唯一解；U 撞界；
  - P5-C（固定 COD ID 扩 phase pack）与 P5-D（`--absorber-site all`
    多位点正式 CN）尚未启动。
- 下一步：评估 P5-C/P5-D；T5 未通过元素不强行给值，门禁不得放宽。

## 9.12 P5-C：固定 COD ID 补齐 18 个常见相

#### 9.12.1 改动

- `tools/cif_crawler.py` 新增 tier=`P5C` 与 `--p5c`；固定目标由
  21 个增至 **39 个**，按 18 个固定 COD ID 补齐
  `phase_pack.COMMON_PHASES` 的 `missing_common_phases`：
  - FeO 1011169、FeOOH 9002158、PtO2 1008935、Cr2O3 9014036、
    MnO 9006658、MnO2 9009081、Mn3O4 1514120、ZnO 9004178、
    MoO2 9009090、MoO3 9009669、RuO2 9007541、Rh2O3 1010584、
    PdO 1009031、Ag2O 1010486、IrO2 2101854、Au2O3 2105390、
    UO2 9009049、U3O8 2310519。
- phase pack 由 26 条增至 **44 条**；`phase_index.json` 44 条；
  `phase_coverage.json` 的 `missing_common_phases` 全部为空。
- 选择规则：优先常压/室温实验结构；高压、高温或 DFT 条目排除
  （Cr2O3 的 63.9 GPa 条目、β-U3O8、高温 U3O8 均未采用）。
  PtO2 在 COD 中没有常压结构，只能采用高压合成 β-PtO2，已在
  phase/notes 明确标注为“非环境相”。
- Mn3O4（Mn1/Mn2）与 U3O8（U1/U2）为多位点相，正式 CN 仍必须
  显式 `--absorber-site`，不自动 Top-1。
- `data/cif_raw/` 保持 17 个根索引 CIF；P5-C 的 18 个下载缓存只在
  构建期使用，导入后清理。
- 测试：#84/#92 改为 44 条；#115 改为 39 个目标并允许 `P5C`；
  新增 #125 锁定 18 个目标、coverage 无缺口、多位点拒绝。
- 账本 #84 的旧断言 `elements_with_motif=9` 与实际
  `motif_coverage.json`（P4+ 后 17/17）不一致；P5-C 重建 coverage
  后暴露，改为 17/0，并清空 `without_motif_elements`。

#### 9.12.2 验收

- 18 条候选全部通过本地 `verify_cif`：gemmi 解析、IT 空间群号、
  化学计量比、吸收位点；全部为 experimental `crystal_cif`。
- 44 条 phase pack 条目逐条 `phase_pack.validate_entry` 通过。
- `manifest=44`、`phase_index=44`、`phase_coverage.verified_cif=44`；
  17/17 元素 `formal_cn_ready=True`，全部 `missing_common_phases=[]`。
- 多次位点仍拒绝自动选择：`COD:1514120`（Mn3O4）与
  `COD:2310519`（U3O8）返回“多个吸收位点，不自动 Top-1”。

#### 9.12.3 测试

- `regression_test.py`：125 项 OK（142.474s）。
- `run_ci.py --stress`：`ok=true`、`contract_tests=125`、
  `stress_tests=4`、PASS（215.1s）；4 个压力测试全部通过。

#### 9.12.4 提交与补丁

- 提交：`c8b8236`（P5-C: add 18 common-phase COD entries）。
- 补丁：`patches/P5-C.patch`（增量 diff `d3e178c..c8b8236`）；
  大小 164019 字节；
  SHA256=`5F712B52E796D1CFD0388D202BF66E5AB01372EF2B52DDF25CF8B2FE0FE19630`；
  `git apply --check --reverse` 通过。
- `patches/P4+.patch`、`patches/P5-A.patch`、`patches/P5-B.patch`
  保持不变。

#### 9.12.5 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：
  - PtO2 只有高压合成 β 相可用，非环境相边界已标注；
  - P5-D（`--absorber-site all` 多位点正式 CN）尚未启动。
- 下一步：P5 遗留项谱源隔离池（Rh/Pd/Ir/Ce/U）。

## 9.13 P5-prep：标准谱源隔离池（P5 遗留项准备）

#### 9.13.1 改动

- 新增 `tools/build_source_pool.py`：从
  `data/s02_calibration/{Rh,Pd,Ir,Ce,U}/calibration.json` 与
  `calibrate_s02.py` 的模块级门禁常量**确定性重建**
  `data/s02_calibration/source_pool/index.json`；门禁常量用 AST 读取，
  不导入/不执行标定模块。
- 新增 `data/s02_calibration/source_pool/README.md`：5 个未通过元素的
  状态分类、失败证据、Gemini 建议逐项评估、下一步动作与替换准入清单。
- **逻辑隔离**：不移动/不删除 IIT-Lytle 原始谱；池内元素
  `s02=null`、`usable_for_s02=false`，运行期不得用于 `--s02-from-t5`；
  池索引不进入生产 `_structure_lib`。
- 分类结果：
  - Rh/Pd/Ir → `short_spectrum_replace`：谱源边后仅约 122 eV，
    `k_max≈5.66–5.68 Å⁻¹`，S5/S6 WT 无有效网格；
  - Ce → `oxide_baseline_reprocess`：主谱已是 BMM
    `Ce-L3-CeO2.xdi`，motif 为 `ceo2_fluorite_cod_1562989`；问题是
    三窗口 ΔE₀=9.78–10.20 eV 贴 10 eV、σ² 压到 1e-4 下界、无
    interior 窗口；
  - U → `oxide_phase_verification`：谱长足够（`k_max≈20.47 Å⁻¹`），
    但 α-U 第一壳模型不稳（S₀²=0.3、Rf≈0.55–0.60、N_idp=2.48<4）；
    UO2 结构 CIF 已由 P5-C 入库（COD 9009049），但本地无 UO2 标准谱，
    切换前必须先核实可追溯谱源。
- 对 Gemini 建议的评估：
  - Rh/Pd“XANES 截断待换谱”**同意**；
  - Ce“换 CeO2”**不适用**：已经是 CeO2 基准；
  - U“核实是否换 UO2”**同意核实、不能直接换**；
  - Ir“提高 k_min 避开白线”**不足以救当前谱**：k_max≈5.66 Å⁻¹，
    k_min≥2.5 Å⁻¹ 时 Δk≈3.16 Å⁻¹，N_idp 仍小于 4 个自由参数，
    且失败发生在 WT 网格为空；仅保留为换谱后的诊断项，
    不得绕过门禁。
- 门禁冻结：`S02_PASS=[0.4,1.15]`、`s02_err<0.06`、`Rf≤0.12`、
  `|ΔE₀|<10 eV`、三窗口 interior、跨窗极差≤0.05、σ² 不贴界
  （优化下界 1e-4、≤5e-4 判 at-bound）、`vars<N_idp`；禁止改门禁/
  手调窗口/自动 Top-1。
- `absorber.py` 静态表不变：Fe/Co/Ni/Cu/Pt 已标定，其余 12 元素
  `s02=None`、`calibrated=False`。
- 新增契约 #126：锁定 5 元素分类、`s02=None`、T5 `not_calibrated`、
  门禁常量与索引可重建性。

#### 9.13.2 验收

- `python tools/build_source_pool.py --check` 通过：索引与标定产物/
  门禁常量一致；`py_compile` 通过。
- 隔离池 5 元素全部 `pool_state=quarantined`、`s02=None`、
  `usable_for_s02=false`；`classification_counts` 为
  `short_spectrum_replace=3`、`oxide_baseline_reprocess=1`、
  `oxide_phase_verification=1`。
- 索引中的数值证据固定为 2026-11-09 观测值，不依赖未提交的
  IIT-Lytle 文件；本机存在源文件时会校验 sha256。

#### 9.13.3 测试

- `regression_test.py`：126 项 OK（176.600s）。
- `run_ci.py --stress`：`ok=true`、`contract_tests=126`、
  `stress_tests=4`、PASS（200.4s；`total_seconds=200.43`）。

#### 9.13.4 提交与补丁

- 提交：`e785bdd`（P5-prep: add standard spectra source isolation pool）。
- 补丁：`patches/P5-prep.patch`（增量 diff `c8b8236..e785bdd`）；
  大小 57002 字节；
  SHA256=`9E13DD666BA4580DC791AD5F8BE9C9F2D8E5BD2532356C49E717698D7983C9D3`；
  `git apply --check --reverse` 通过。
- `patches/P4+.patch`、`patches/P5-A.patch`、`patches/P5-B.patch`、
  `patches/P5-C.patch` 保持不变。

#### 9.13.5 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：
  - Rh/Pd/Ir 需换 `k_max≥10–12 Å⁻¹` 的实验谱；
  - Ce 需复核归一化/E0 方法，仍无唯一解则保持 `not_calibrated`；
  - U 需先核实 UO2 L3 标准谱源许可与可追溯性；
  - P5-D（`--absorber-site all` 多位点正式 CN）尚未启动。
- 下一步：评估 P5-D 是否启动；未启动前保持 P5-D 暂缓，池内元素
  不强行给 S₀²。

## 9.14 P5-S02 隔离池换谱闭环 + P5-D 多位点正式 CN（2026-11-09）

#### 9.14.1 P5-S02：5 元素换谱并通过冻结门禁

- `cif_crawler.py` 新增 IIT-Lytle RAW 解析器与谱源：
  - Rh `RAW/Rh/rhfoil.003`（514 点，22694–24807 eV，透射）；
  - Pd `RAW/Pd/pdfl.404`（520 点，24156–26353 eV，NSLS 透射，
    数据列扣头部 `OFFSETS` I0=5172、I1=5557）；
  - Ir `RAW/Ir/irfoil.146`（645 点，11015–12820 eV，步进扫描；
    第 3 列随 L3 吸收升高，解析器取 `-ln(I0/I1)` 恢复 μ 边后升高
    约定，provenance 记录符号约定）；
  - Ce `RAW/Ce/ceo2.004`（327 点，5523–6158 eV，在 Ce L2 6164 eV
    前截断）；
  - U `RAW/U/uo2.379`（381 点，16966–18170 eV）。
- `calibrate_t5.py` 更新 5 元素 CASES；U 使用标定专用
  `uo2_fluorite`（COD 9009049，8×U–O 2.3678 Å），仅写
  `data/s02_calibration/motifs/`，不进生产 `_structure_lib`。
- 结果：T5 **12/12**，冻结门禁零放宽：

  | 元素 | S₀² | err | Rf | ΔE₀ (eV) | spread | plateau |
  |------|-----|-----|----|-----------|--------|---------|
  | Rh | 0.7022 | 0.0191 | 0.0041 | −0.71 | 0.0037 | 0.8590 |
  | Pd | 0.6559 | 0.0141 | 0.0035 | +1.14 | 0.0133 | 0.8138 |
  | Ir | 0.4752 | 0.0149 | 0.0043 | +9.30 | 0.0068 | 0.8604 |
  | Ce | 0.6279 | 0.0268 | 0.0223 | +5.94 | 0.0096 | 1.0645 |
  | U  | 0.7800 | 0.0229 | 0.0137 | +2.86 | 0.0024 | 0.8608 |

- 隔离池 `source_pool/index.json` 升级 schema 2：
  `pool_state="resolved"`、`usable_for_s02=true`、
  `pool_state_counts={"resolved":5}`；旧谱 sha256/点数/失败门禁保留在
  `quarantined_legacy_source`；`resolved_by.standard` 与
  `current_source` 一致。索引数值全部来自已提交 `calibration.json`
  + 固定观测常量，不依赖本机 `_standards_lib`。
- Ce 的 `attempts` 审计链保持：
  `Ce-L3-CeO2.xdi`（gated_out，σ² 贴界）→
  `Ce-L3-CeriumDioxide.xdi`（error，平台 2.271）→
  `Ce_metal_foil.dat`（passed）。
- `absorber.py` 静态表不变；`#60` 的 12 元素仍 `s02=None`/
  `calibrated=False`。

#### 9.14.2 P5-D：`--absorber-site all` 多位点正式 CN

- `phase_resolver.PhaseResolution` 新增 `site_models`；
  `resolve_phase(..., absorber_site="all")` 把多位点候选展开成每个
  不等价吸收位一个 `StructureModel`。未传 `--absorber-site` 时仍
  返回 `structure_model=None` + `absorber_site_ambiguous` 原因，
  不自动 Top-1。
- `cif_gate` 新增：
  - `list_inequivalent_absorber_sites`：读取每位的 multiplicity/occ，
    计算 `w_i = m_i·occ_i / Σ(m_j·occ_j)`；
  - `gate_cif_all_sites`：每位点独立过同一套冻结门禁，并附
    `site_multiplicity/site_occupancy/site_weight`。
- `xafs_fit.run_multi_site_formal_cn`：
  - 每位点构建 CIF 位点团簇 motif（不写生产库）→ 条件 CN →
    `build_formal_cn_block`；
  - 成员行带 `multi_site_member` 与 `site_weight=<w>`；
  - **仅当全部位点 ok** 才追加 `absorber_site="ALL_WEIGHTED"` 汇总行；
    任一位点失败则该位点写拒绝行，且不输出加权汇总；
  - 汇总 N₀/N_fit/σ²/ΔR/ΔE₀ 按 `w_i·N₀_i` 加权，Rf 取最大。
- `run_xafs.py` 集成：`--absorber-site all` 先对全部位点过门禁，
  任一位点失败即抛 `XAFSConfigError`；`xafs_fit.run` 后用多位点结果
  覆盖 `formal_cn`。`--absorber-site` 帮助文本已更新。
- 权重与 N₀：
  - Fe3O4（COD 1010369）：Fe1/Fe2 mult=16、occ=0.5，Fe3 mult=8、
    occ=1 → 权重各 1/3；第一壳 N₀=6/6/4 → 加权 N₀=16/3=5.3333；
  - Co3O4（COD 9005888）：Co1(8a):Co2(16d)=1/3:2/3；但该 COD CIF
    缺少 `_atom_site_occupancy`，冻结门禁按 `occupancy_missing` 拒绝，
    不出正式 CN（无合格来源不出 CN）。
- 真实 Fe3O4 谱端到端（BMM `Fe-Magnetite.xdi`，
  `--phase-id COD:1010369 --absorber-site all --offline`）：
  3 位点门禁全部通过（各 1/3）；Fe1/Fe2 条件 CN 因 N/N₀=0.553
  超出范围被拒绝，Fe3 通过（N/N₀=0.998）；
  最终 `status=multi_site_failed`，`cn_table.csv` 写拒绝行且无
  `ALL_WEIGHTED` 汇总行，符合“失败位点不出汇总”的契约。
- 契约 #127 锁定：
  未传 site 拒绝、Fe3O4 权重/加权 N₀、Co3O4 1/3:2/3、
  Co3O4 `occupancy_missing` 拒绝、全部位点成功才出汇总行。

#### 9.14.3 测试

- 定向测试：`#122`、`#123`、`#124`、`#126`、`#127` 全部 OK。
- 全量 `regression_test.py`：**127 项 OK**（181.507s）。
- `run_ci.py --stress`：`ok=true`、`contract_tests=127`、
  `stress_tests=4`、PASS（`total_seconds=202.3`，最终提交态）。
- `py_compile`：`cif_crawler.py`、`calibrate_t5.py`、
  `tools/build_source_pool.py`、`phase_resolver.py`、`cif_gate.py`、
  `xafs_fit.py`、`run_xafs.py`、`regression_test.py` 全部通过。
- `python tools/build_source_pool.py --check`：通过。

#### 9.14.4 提交与补丁

- 提交：`13fba39`（P5-S02 + P5-D: close 5-element S02 pool and
  multi-site formal CN）。
- 补丁：`patches/P5-S02-and-P5D.patch`（增量 diff `e785bdd..13fba39`）；
  大小 209770 字节；
  SHA256=`5E86EC05257919DF088AD4E0930DA90EED002E2E27FDC4CF03AAB5BFCCF21567`；
  `git apply --check --reverse` 通过。
- `_standards_lib/` 与 `data/standards_raw/` 中 5 个新 `.dat` 未
  `git add`；只提交代码、`data/s02_calibration/` 产物与文档。
- `patches/P4+.patch`、`patches/P5-A.patch`、`patches/P5-B.patch`、
  `patches/P5-C.patch`、`patches/P5-prep.patch` 保持不变。

#### 9.14.5 结论

- [x] 通过
- [ ] 不通过，原因：
- 遗留问题：
  - IIT-Lytle 谱许可未明确，原始谱不提交；标定产物记录 SHA256，
    换机后需重新获取同源谱才能复跑标定（运行期只读已提交 JSON）；
  - Co3O4 COD 9005888 缺少 occupancy 标签，P5-D 权重可计算但门禁
    拒绝；若需 Co3O4 正式 CN，必须换带 `_atom_site_occupancy` 的
    合格 CIF，不得放宽门禁；
  - P5-D 真实 Fe3O4 示例因 Fe1/Fe2 条件 CN 的 N/N₀=0.553 被拒绝，
    属于样本/模型一致性问题，不是门禁放宽的理由。
- 下一步：保持 12/12 冻结结果；后续新相/新谱按同一流程先过门禁，
  不自动 Top-1。


## 10. P6：催化分析门禁修复与项目整合（2026-11-09）

### 10.1 范围

- `cif_gate.py`：单晶位对称展开判定、隐式占据度、M–O 实验键长参考；
- `Rh2O3`：`COD:1010584` → `COD:1544346`（Shannon & Prewitt 1970）；
- CLI：`--s02-from-t5` 提示、`--catalytic-model + --formal-cn` 的
  `hypothesis_model_only` 边界；
- 回归：`#127` 更新、`#128` 新增（44/44 phase pack 全量门禁）；
- 项目整理：`docs/` 收编、根目录日志归档、`.gitignore` 增 `logs/`。

### 10.2 结果

- 44 条 phase pack：39 条单晶位全部 `ok=True`；5 条多晶位中
  4 条多元素继续 `absorber_site_ambiguous`（Fe3O4 / Mn3O4 / U3O8 /
  Co3O4），Mn 单元素自动选位；`gate_cif_all_sites` 全部位点通过，
  权重和 = 1.0。
- `PtO2 / IrO2 / Au2O3 / CeO2 / MoO3 / Cr2O3` 等高价氧化物键长门禁通过。
- `Co3O4` 按 `implicit_1.0` 放行，加权 N₀ = 16/3。
- `Rh2O3 (COD:1544346)` 第一壳 `6 × Rh–O @ 2.050 Å`。

### 10.3 测试

- `regression_test.py`：**128 项 OK**（144.7 s）；
- `run_ci.py --stress`：128 契约 + 4 压力 PASS（173.0 s）；
- 定向：#96 / #97 / #98 / #127 / #128 全部 OK。

### 10.4 提交与补丁

- 提交：`P6: fix phase-pack CIF gates (44/44), update Rh2O3 COD 1544346,
  and consolidate docs`
- 补丁：`patches/P6-catalytic-readiness-and-cleanup.patch`
- `_standards_lib/`、`data/standards_raw/`、`logs/` 未提交。

### 10.5 结论

- [x] 通过
- 遗留问题：
  - IIT-Lytle 谱许可未明确，原始谱仍不提交；
  - 12 个 T5 元素静态表保持 `s02=None`，运行期只经 `--s02-from-t5` 解锁；
  - 多元素多晶位相仍必须显式 `--absorber-site` 或 `all`，不自动 Top-1。

## 11. P6.1：单晶位/多晶位正式 CN 端到端闭环（2026-11-09）

### 11.1 范围

- `xafs_fit.py`：新增 `_attach_single_site_formal_cn`，让单晶位正式 CN
  与诊断轨 `_structure_lib` 解耦；单晶位 P5-C 氧化物与用户 `--phase-cif`
  不再因 S6.5 未命中而退回 `l1_diagnostic_only`；
- `phase_resolver.py`：`occupancy_policy` 按 CIF 是否含
  `_atom_site_occupancy` 初始化为 `explicit` / `implicit_1.0`；
- `run_xafs.py`：门禁结果回填 `occupancy_source` / `n0_shells` /
  `warnings`；多位点正式 CN 后同步 `fit.cn_status`；
- `xafs_report.py`：新增 `--s02-from-t5` 与 `--phase-*` /
  `--absorber-site` / `--offline`，接入 T5 标定、相解析与 CIF 门禁；
- `xafs_fit.py` 三处 S₀² 报错补 `--s02-from-t5` 提示；
- 契约 `#129` 锁定 RuO2 单晶位端到端正式 CN；
- 契约 `#130` 锁定 `frac:(x,y,z)` 分数坐标吸收位、`cif_gate` 晶胞正交化、P1b/P5 配置文件键镜像、XDI 文件头元素回退与 `xafs_report.py` 异常分支对齐。

### 11.2 结果

- RuO2（`COD:9007541`）：诊断轨仍 `l1_first_shell`，正式轨不再被劫持，
  输出物理拒绝 `N/N₀=0.521 out_of_range`；`occupancy_policy=implicit_1.0`。
- NiO（`COD:1010093`）：诊断 `motif`、正式 `ok`；首壳
  `Ni-O@2.08, N₀=6.0, N_fit=6.327`；无回归。
- Fe3O4（`COD:1010369`，`--absorber-site all`）：`fit.cn_status` 与
  `_meta.json` 均为 `multi_site_failed`；位点权重 1/3 : 1/3 : 1/3。
- `xafs_report.py --s02-from-t5 --formal-cn --phase-id COD:9007541` 可完整
  跑通并写出 `cn_table.csv`。

### 11.3 测试

- `regression_test.py`：**130 项 OK**；
- `run_ci.py --stress`：130 契约 + 4 压力 PASS；
- 定向 `#75 / #76 / #127 / #128 / #129 / #130`：6/6 OK。

### 11.4 提交与补丁

- 提交：`P6.1: decouple single-site formal CN from diagnostic library`
- 补丁：`patches/P6.1-single-site-formal-cn.patch`
- `_standards_lib/`、`data/standards_raw/`、`logs/` 未提交。

### 11.5 结论

- [x] 通过
- 遗留问题：
  - 稀释样品（如 Ru-K-RuO2.xdi）边后平台 0.106，仍须 `--no-strict`；
  - 条件 CN 的物理拒绝（N/N₀ 越界、Rf 超限）保持硬拒绝，不放宽。
