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
- `regression_test.py` **119 项** + 4 个压力测试全部通过（`run_ci.py --stress` PASS）。

运行期仍然完全离线：`run_xafs.py` 不导入根爬虫，不读取 `data/cif_crawler_index.json` 取几何，只通过
`structure_lib.load_library(absorber)` 读取 `_structure_lib`。

## 2. 两套构建资产，边界独立

| 工具 | 职责 | 产物 | 运行期依赖 |
|---|---|---|---|
| 根目录 `cif_crawler.py` | COD motif 抓取/核验/入库 + T5 标准谱抓取 | `data/cif_raw/`、`data/standards_raw/`、`data/cif_crawler_index.json`、`_structure_lib/<元素>/*.json`、`_standards_lib/<元素>/*.dat` | **否** |
| `tools/cif_crawler.py` | 仅做 phase pack 的固定 ID CIF 导入与索引重建 | `_structure_lib/_phase_pack/`、`phase_index.json`、`phase_coverage.json` | **否** |
| `structure_lib.load_library()` | 运行期读取诊断 motif 库 | `_structure_lib/<元素>/*.json` | **是** |

硬边界：

1. 运行期不得直接读 `data/cif_crawler_index.json` 取几何；
2. 运行期不得调用任一爬虫、不得联网；
3. 爬虫不得成为 S1–S7 的隐式数据源；
4. T5 标准谱目前只落盘，不参与 `calibrate_s02.py`，也不参与正式 CN。

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
- **尚未接入 `calibrate_s02.py` 或运行期**；#60 静态表仍为 Fe/Co/Ni/Cu/Pt 已标定，其余 12 元素 `s02=None`、`calibrated=False`。

## 5. 验收

### 5.1 回归与压力

```text
regression_test.py: 119 tests OK
run_ci.py --stress: PASS (170.5s)
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

未索引的 9 个重复/旁路 CIF 已清理（5 个大小写重复副本 + Cr/Zn/Mo/Au 4 个仅 phase pack 使用的旁路缓存）。

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

1. **T5 S₀² 标定未完成**：Cr/Mn/Zn/Mo/Ru/Rh/Pd/Ag/Ir/Au/Ce/U 12 元素 `s02=None`，正式 CN 仍会被拒绝；
2. **常见相仍缺**：FeO、FeOOH、PtO2、Cr2O3、MnO/MnO2/Mn3O4、ZnO、MoO2/MoO3、RuO2、Rh2O3、PdO、Ag2O、IrO2、Au2O3、UO2/U3O8；
3. **多位点联合模型未自动构建**：当前需显式 `--absorber-site` 或 `absorber_fract`；
4. **T7 多元素联合、T8 XANES-EXAFS 联合校验**仍待执行。

## 8. 一句话状态

**根目录 `cif_crawler.py` 已成为主项目的构建期资产：29 条索引、17 个 COD motif、12 条 T5 真实谱、43 个 motif/17 元素覆盖，运行期离线且不依赖索引；旧 ID 迁移和 119+4 CI 全部通过，T5 标定仍是后续独立任务。**
