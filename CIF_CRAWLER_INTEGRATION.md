# CIF Crawler 整合报告（主项目）

日期：2026-11-09  
来源工作区：`C:\Users\lizha\.local\share\opencode\worktree\d5e1a5\witty-rocket`  
主项目：`C:\Users\lizha\Desktop\Code-experience-cn`

## 1. 结论

根目录 `cif_crawler.py` 已作为**构建期**的权威 motif / 标准谱构建工具并入主项目；它不在 S1–S7 运行路径上。

本批次完成：

- 29 条构建索引：**17 条 COD CIF + 12 条 T5 金属箔标准谱**；
- 17 个 COD motif 入库，结构库共 **43 个 motif**，覆盖 **17/17** 吸收元素；
- T5 12 条 IIT-Lytle 真实谱全部落盘，**0 个 placeholder**；
- 旧 ID 迁移完成：`LEGACY_MOTIFS_LEFT 0`、`QUARANTINE_FILES 0`；
- `regression_test.py` **125 项** + 4 个压力测试全部通过（`run_ci.py --stress` PASS；P5-C 后）。

运行期仍然完全离线：`run_xafs.py` 不导入根爬虫，不读取 `data/cif_crawler_index.json` 取几何，只通过
`structure_lib.load_library(absorber)` 读取 `_structure_lib`。

## 2. 两套构建资产，边界独立

| 工具 | 职责 | 产物 | 运行期依赖 |
|---|---|---|---|
| 根目录 `cif_crawler.py` | COD motif 抓取/核验/入库 + T5 标准谱抓取 | `data/cif_raw/`、`data/standards_raw/`、`data/cif_crawler_index.json`、`_structure_lib/<元素>/*.json`、`_standards_lib/<元素>/*.dat` | **否** |
| `tools/cif_crawler.py` | 仅做 phase pack 的固定 ID CIF 导入与索引重建（T3/T4/P5-C） | `_structure_lib/_phase_pack/`、`phase_index.json`、`phase_coverage.json` | **否** |
| `structure_lib.load_library()` | 运行期读取诊断 motif 库 | `_structure_lib/<元素>/*.json` | **是** |

硬边界：

1. 运行期不得直接读 `data/cif_crawler_index.json` 取几何；
2. 运行期不得调用任一爬虫、不得联网；
3. 爬虫不得成为 S1–S7 的隐式数据源；
4. T5 标准谱只作为**构建期产物**：`calibrate_t5.py` 读取本地标准谱生成独立标定产物，运行期不读标准谱、也不成为 S1–S7 的隐式数据源。

## 3. 适配与旧 ID 迁移

### 3.1 参考版适配

- `ROOT_DIR` 固定为脚本所在目录，避免从不同工作目录运行时把产物写到别处；
- 索引中的 `cif_file` / `standard_file` 统一写成项目内相对路径，不固化构建机绝对路径；
- 为多吸收位目标显式传 `absorber_fract`：
  - Fe3O4 `[0.625, 0.625, 0.625]`；
  - Co3O4 `[0.125, 0.125, 0.125]`；
  - α-Mn `[0.0, 0.0, 0.0]`；
- `import_cif(..., absorber_fract=...)` 只增加这一个参数；没有合入 `EnhancedMotifBuilder`，没有给 `Motif` 增加 `meta` 字段，也没有 `center_site` 参数或 BVS 质量评估改动。

### 3.2 旧 ID 迁移

爬虫按 `LEGACY_ARTIFACTS` 与指纹碰撞清理了 6 个旧 motif；其中用户点名的两条为：

| 旧 ID | 新 ID | 说明 |
|---|---|---|
| `nio_rocksalt` | `nio_bunsenite_cod_1010093` | COD 1010093，a=4.1684；Ni-O 2.0842×6、Ni-Ni 2.9475×12 |
| `cu2o_cuprite` | `cu2o_cuprite_cod_1010941` | COD 1010941，a=4.26；O 1.8446×2、Cu 3.0123×12 |

同批因指纹碰撞淘汰的还有：

- `coo_rocksalt` → `coo_rocksalt_cod_1533087`
- `co3o4_tet` → `co3o4_spinel_cod_9005888`
- `ni_oh2_brucite` → `nioh2_theophrastite_cod_1548811`
- `cuo_tenorite` → `cuo_tenorite_cod_1011148`

未碰撞的旧 Fe 相 `hematite`、`magnetite_oct/tet/mixed`、`wustite`、`goethite`、`lepidocrocite` 保留，避免影响磁铁矿与 LEPIDO 契约。

`co3o4_mixed` 保持“1/3 tet + 2/3 oct”语义，但 tet 组件改为 `co3o4_spinel_cod_9005888`，oct 组件仍为 `co3o4_oct`。

同步修改：

- `regression_test.py`：NiO/Cu2O/CoO/Co3O4/Ni(OH)2/CuO 的 ID、壳层断言、计数（43 motifs / 17 元素覆盖 / 0 missing）；
- `stress_test_nio.py`：删除手工注册，直接加载 `nio_bunsenite_cod_1010093`；
- `_structure_lib/build_seed.py`：Cu2O seed 改为 COD 1010941、a=4.26、source COD 1010941；
- `_structure_lib/Co/co3o4_mixed.json`：tet 组件迁移。

## 4. 构建结果

### 4.1 命令

```text
python cif_crawler.py --all --import-lib
python cif_crawler.py --t5
```

两次命令均 exit 0。

### 4.2 索引

`data/cif_crawler_index.json` 共 **29 条**：

- **17 条 CIF**：Fe2O3、Fe3O4、CoO、Co3O4、NiO、Ni(OH)2、CuO、Cu2O、CeO2、Mn、Ru、Rh、Pd、Ag、Ir、Ce、U；
- **12 条 STD**：Cr、Mn、Zn、Mo、Ru、Rh、Pd、Ag、Ir、Au、Ce、U。

索引中的路径全部为 `data/cif_raw/...`、`data/standards_raw/...` 形式的相对路径。

### 4.3 结构库

- 17 个 COD motif：`fe2o3_hematite_cod_1011240`、`fe3o4_magnetite_cod_1010369`、`coo_rocksalt_cod_1533087`、`co3o4_spinel_cod_9005888`、`nio_bunsenite_cod_1010093`、`nioh2_theophrastite_cod_1548811`、`cuo_tenorite_cod_1011148`、`cu2o_cuprite_cod_1010941`、`ceo2_fluorite_cod_1562989`、`mn_alpha_mn_cod_9008589`、`ru_hcp_cod_9008513`、`rh_fcc_cod_9008482`、`pd_fcc_cod_9008478`、`ag_fcc_cod_9008459`、`ir_fcc_cod_9008470`、`ce_fcc_cod_9010987`、`u_alpha_u_cod_1527120`；
- 结构库共 **43 个 motif**（原 32 − 6 + 17）；
- `motif_coverage.json`：17/17 元素覆盖、`missing_elements=[]`、`formal_cn_ready_elements=[]`；
- `provenance_report.json`：`by_status={"constructed":9,"generated_only":32,"source_ref_only":2}`、`unclassified=0`、`verified_source_cif=0`。

注意：COD motif 是 `generated_only` 诊断候选，不是正式 CN 来源。正式 CN 仍只接受 `verified_source_cif` + `crystal_cif/crystal_edited`，并通过 `validate_formal`。

### 4.4 T5 标准谱

- 12 个元素均有真实谱，写入 `_standards_lib/<元素>/<元素>_metal_foil.dat`；
- `data/standards_raw/` 保存同一份抓取缓存；
- `T5_PLACEHOLDERS 0`；
- 文件头保留 IIT-Lytle 来源、源文件名、E0、能量范围与点数；
- T5 标定由 `calibrate_t5.py` 批量执行，产物落 `data/s02_calibration/`（`index.json`、逐元素 `calibration.json`、`REPORT.md`），**不改 `absorber.py`**；#60 静态表仍为 Fe/Co/Ni/Cu/Pt 已标定，其余 12 元素 `s02=None`、`calibrated=False`；
- 运行时通过 `run_xafs.py --s02 <值>` 显式覆盖，或对已通过元素用
  `run_xafs.py --s02-from-t5 --identify/--fit` 从独立产物读取（只接受
  `status="calibrated"` 且 `static_absorber_py_unchanged=True`；未通过
  元素明确拒绝）；通过 7/12（Cr/Mn/Zn/Mo/Ru/Ag/Au），Rh/Pd/Ir/Ce/U
  记录未通过原因、不强行给值；`s02_source` 溯源串写入 L1/motif 输出与
  `<stem>_meta.json.fit`；
- Mo 于 P5-B 通过：BMM Mo 箔用 `half_height` E₀（20012.86 eV，平台
  0.863）+ 双壳 bcc（8×2.725 + 6×3.147 Å，共享 S₀²/ΔE₀、各壳独立
  ΔR/σ²），S₀²=0.8322±0.0334、Rf=0.0040、窗口极差 0.0017、ΔE₀=−0.73 eV；
  门禁含新增 `sigma2_not_at_bound`（σ² 压到 1e-4 下界即失败）。
  Ce 仍不通过：自由 σ² 压界、固定 σ² 无唯一解、加第二壳 Rf=0.127；
  备选透射谱平台 2.27 属归一化伪迹，严格门禁拒绝；
- `_standards_lib/` 与 `data/standards_raw/` 因许可未明确不提交；Ag 标定依赖本地 IIT-Lytle 谱，产物记录源文件 SHA256，但不保证无本地标准时可复现。

### 4.5 P5-C 常见相补齐

`tools/cif_crawler.py --p5c --online` 按固定 COD ID 补齐 18 个
`missing_common_phases`，phase pack 由 26 条增至 **44 条**，固定抓取目标
由 21 个增至 **39 个**。18 条全部为实验 `crystal_cif`，逐条通过 gemmi
解析、IT 空间群号与化学计量核验：

| 元素 | 相 | COD ID | 空间群 | 备注 |
|------|----|--------|--------|------|
| Fe | FeO | 1011169 | 225 | wüstite |
| Fe | FeOOH | 9002158 | 62 | goethite，298 K |
| Pt | PtO2 | 1008935 | 58 | β-PtO2；COD 无常压结构，已标注高压合成 |
| Cr | Cr2O3 | 9014036 | 167 | eskolaite，常压 |
| Mn | MnO | 9006658 | 225 | manganosite，P=0 |
| Mn | MnO2 | 9009081 | 136 | pyrolusite |
| Mn | Mn3O4 | 1514120 | 141 | hausmannite，Mn1/Mn2 |
| Zn | ZnO | 9004178 | 186 | zincite，293 K |
| Mo | MoO2 | 9009090 | 14 | tugarinovite |
| Mo | MoO3 | 9009669 | 62 | molybdite |
| Ru | RuO2 | 9007541 | 136 | rutile |
| Rh | Rh2O3 | 1010584 | 167 | corundum |
| Pd | PdO | 1009031 | 131 | palladinite |
| Ag | Ag2O | 1010486 | 224 | cuprite |
| Ir | IrO2 | 2101854 | 136 | rutile |
| Au | Au2O3 | 2105390 | 43 | gold(III) oxide |
| U | UO2 | 9009049 | 225 | uraninite |
| U | U3O8 | 2310519 | 21 | α-U3O8，中子衍射，U1/U2 |

高压/高温/DFT 条目显式排除（Cr2O3 的 63.9 GPa 条目、β-U3O8、高温
U3O8 均未采用）。`phase_coverage.json` 的 `missing_common_phases` 全部
为空；Mn3O4 与 U3O8 多位点仍需显式 `--absorber-site`。P5-C 下载缓存
只用于构建期导入，导入后清理，`data/cif_raw/` 仍保持 17 个根索引 CIF。

## 5. 验收

### 5.1 回归与压力

```text
regression_test.py: 125 tests OK (142.474s)
run_ci.py --stress: PASS (215.1s)
  stress_test_co.py   PASS
  stress_test_ni.py   PASS
  stress_test_cu.py   PASS
  stress_test_nio.py  PASS
```

NiO 压力测试：top1=`nio_bunsenite_cod_1010093`、score=0.008，`ni_fcc=0.832` 垫底；P5 Ni-O N=7.62/R=2.0786、Ni-Ni N=11.92/R=2.9368、Rf=0.0059。

### 5.2 LEPIDO 零退化

三扫描基线命令重跑后：

```text
r_factor = 0.027281059410161693
N(Fe-O)  = 6.348637726145664
R(Fe-O)  = 2.0019558444559276
fit_table SHA256 = E098D0EDC24E49296710E44518F8A8F010DCD03D7D7EB8CB2E26831A7532E760
```

结构库从 11 条 Fe 候选扩到 13 条后，`lepidocrocite` 仍为 top1（score=0.149），新增的
`fe2o3_hematite_cod_1011240` 与 `fe3o4_magnetite_cod_1010369` 只进入候选排序（第 4/12 位），S7 拟合值逐位不变。

回归测试中的 LEPIDO 片段只喂入 `.000` 单扫描，因此显示 `Rf=0.0270`、`N=6.41`；这不是退化，基线是 `.000/.100/.200` 三扫描合并。

## 6. 数据提交与许可

### 6.1 随仓库提交

- `data/cif_raw/` 中 **17 个被索引引用**的 COD CIF；
- `data/cif_crawler_index.json`；
- `_structure_lib/**/*.json`、`_structure_lib/*/_cif/*.cif`、`_structure_lib/_phase_pack/`；
- 整合文档与测试修改。

未索引的 9 个重复/旁路 CIF 已清理（5 个大小写重复副本 + Cr/Zn/Mo/Au 4 个仅 phase pack 使用的旁路缓存）。P5-C 新增的 18 个下载缓存同样只在构建期使用，导入后清理；phase pack 内的 `<sha256>.cif` 才是权威副本。

### 6.2 不随仓库提交

`data/standards_raw/` 与 `_standards_lib/` 是 IIT-Lytle 真实谱的本地构建产物。IIT-Lytle 站点只说明数据公开可用，未给出明确再分发许可证，因此本批次不把它们纳入版本库；文件头保留来源，任何有权限的环境可执行：

```text
python cif_crawler.py --t5
```

离线重建。`data/cif_crawler_index.json` 中的 `STD_*` 条目仍保留，作为该次构建的 manifest。

### 6.3 运行期几何来源

正式 CN 的几何只来自：

- `phase_pack` 的 `<sha256>.cif`（`verified_source_cif`）；
- 用户显式 `--phase-cif`；
- 或经 `structure_lib.load_library()` 读取的诊断 motif（不自动升级为正式来源）。

## 7. 仍然存在的缺口

1. **T5 S₀² 部分标定**：已产出独立标定产物并通过 7/12；Rh/Pd/Ir 缺足够宽能区标准、Ce 的 CeO₂ 标准 σ² 与 S₀² 强相关（无唯一解）、U 撞界；未标定元素仍需运行时 `--s02` 才能做正式 CN（已通过元素可用 `--s02-from-t5`）；
2. **常见相已补齐（P5-C）**：18 个 `missing_common_phases` 全部有 verified CIF；唯一边界是 PtO2 在 COD 中只有高压合成 β-PtO2、没有常压 α 相，已在条目 notes 标注；Mn3O4/U3O8 为多位点，正式 CN 仍需显式 `--absorber-site`；
3. **多位点联合模型未自动构建**：当前需显式 `--absorber-site` 或 `absorber_fract`；
4. **T7/T8 已实现**：`--multi-element` 逐元素并列报告、不推荐 Top-1；`--xanes-cross-check` 只报警告、不改 `cn_reliable`（#120/#121）。

## 8. 一句话状态

**根目录 `cif_crawler.py` 已成为主项目的构建期资产：29 条索引、17 个 COD motif、12 条 T5 真实谱、43 个 motif/17 元素覆盖，运行期离线且不依赖索引；旧 ID 迁移和 125+4 CI 全部通过；phase pack 已由 P5-C 扩到 44 条；T5 已产出 7/12 独立标定产物（静态表未改），运行期可用 `--s02-from-t5` 溯源，T7 多元素与 T8 联合校验已实现。**
