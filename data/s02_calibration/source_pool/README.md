# S₀² 标准谱源池（P5-S02 已闭环）

日期：2026-11-09
范围：原 T5 未通过的 5 个元素——**Rh、Pd、Ir、Ce、U**
产物：`data/s02_calibration/source_pool/index.json`（schema 2）
生成/校验：`python tools/build_source_pool.py [--check]`

## 1. 当前状态

P5-S02 已完成 5 元素换谱闭环：**5/5 `pool_state="resolved"`、
`usable_for_s02=true`，T5 总索引 12/12 通过**。池索引同时保留
`schema_version=1` 时代的历史隔离证据：

- `resolved_by`：新谱路径、源文件 sha256、motif、窗口参数与门禁结果；
- `quarantined_legacy_source`：旧谱路径、sha256、点数、能量范围、
  失败门禁与错误原文（用于审计，不再作为 `current_source`）；
- `current_source` 与 `resolved_by.standard` 指向同一份新谱。

边界不变：

- `absorber.py` 静态表仍只有 Fe/Co/Ni/Cu/Pt 已标定，其余 12 元素
  `s02=None`、`calibrated=False`（#60）；
- T5 结果只落 `data/s02_calibration/`，运行期必须显式 `--s02` 或对
  已通过元素使用 `--s02-from-t5`；
- 池内谱源不进入生产 `_structure_lib`，不成为 S1–S7 的隐式数据源；
- IIT-Lytle 谱文件不提交，索引中的数值来自已提交的
  `calibration.json` 与固定观测常量，可在无本地谱的机器上重建。

## 2. 冻结门禁（未放宽）

来源：`calibrate_s02.py` 模块级常量，由 `tools/build_source_pool.py`
用 AST 读取，不导入/不执行标定模块。

| 门禁 | 阈值/规则 |
|------|-----------|
| S₀² 通过区间 | `[0.4, 1.15]` |
| S₀² 误差 | `< 0.06` |
| Rf | `≤ 0.12` |
| ΔE₀ | `|ΔE₀| < 10 eV` |
| 鲁棒性窗口 | `WINDOW_FRACS=(0.12, 0.18, 0.26)`，三窗口都必须 interior |
| 跨窗口极差 | `≤ 0.05` |
| σ² | 不得贴界；优化下界 `1e-4 Å²`，`≤5e-4 Å²` 判为 at-bound |
| 自由度 | 每个窗口 `n_variables < N_idp` |
| 其他禁止项 | 不改门禁、不手调窗口凑数、不自动 Top-1 |

**这些值在池建立后保持不变。** P5-S02 的通过值是原样跑冻结门禁的
结果，没有调低任何阈值。

## 3. 解析结果总表

| 元素 | 新谱（IIT-Lytle RAW） | edge | 点数 | 能量范围 (eV) | S₀² | err | Rf | ΔE₀ (eV) | spread | plateau |
|------|------------------------|------|------|----------------|-----|-----|----|-----------|--------|---------|
| Rh | `Rh/rhfoil.003` | K | 514 | 22694–24807 | 0.7022 | 0.0191 | 0.0041 | −0.71 | 0.0037 | 0.8590 |
| Pd | `Pd/pdfl.404` | K | 520 | 24156–26353 | 0.6559 | 0.0141 | 0.0035 | +1.14 | 0.0133 | 0.8138 |
| Ir | `Ir/irfoil.146` | L3 | 645 | 11015–12820 | 0.4752 | 0.0149 | 0.0043 | +9.30 | 0.0068 | 0.8604 |
| Ce | `Ce/ceo2.004` | L3 | 327 | 5523–6158 | 0.6279 | 0.0268 | 0.0223 | +5.94 | 0.0096 | 1.0645 |
| U  | `U/uo2.379` | L3 | 381 | 16966–18170 | 0.7800 | 0.0229 | 0.0137 | +2.86 | 0.0024 | 0.8608 |

逐元素要点：

- **Rh**：`half_height` E₀ + `strict=True`；平台 0.859 通过严格归一化。
- **Pd**：NSLS X11 透射谱，数据列扣 `OFFSETS`（I0=5172、I1=5557）；
  `strict=False` 记录 plateau=0.8138，五项拟合门禁仍全部通过。
- **Ir**：步进扫描（step→energy），第 3 列随 L3 吸收升高，解析器取
  `-ln(I0/I1)` 恢复 μ 边后升高约定；`kmin=3.0` 仅用于规避 L3 近边
  白线共振，门禁未放宽。
- **Ce**：CeO₂ L3 透射谱，在 Ce L2（6164 eV）前截断于 ≤6160 eV；
  BMM 的两条失败谱保留在 `calibration.json.attempts` 审计链：
  - `Ce-L3-CeO2.xdi`：σ² 贴界、无 interior 窗口（`gated_out`）；
  - `Ce-L3-CeriumDioxide.xdi`：平台 2.271、归一化伪迹（`error`）。
- **U**：UO₂ L3 透射谱 + 标定专用 `uo2_fluorite` motif
  （COD 9009049，8×U–O 2.3678 Å）；motif 只在
  `data/s02_calibration/motifs/`，不进生产 `_structure_lib`。

## 4. 旧谱隔离证据（保留，不再使用）

| 元素 | 旧谱 | 旧 sha256 前 12 位 | 失败类别 | 旧 k_max 估计 |
|------|------|--------------------|----------|----------------|
| Rh | `_standards_lib/Rh/Rh_metal_foil.dat`（`r/rhfoil.a03`） | `074792a1b5de` | WT 无有效网格 | 5.67 Å⁻¹ |
| Pd | `_standards_lib/Pd/Pd_metal_foil.dat`（`p/pdfoil.a71`） | `9f07525aaa52` | WT 无有效网格 | 5.68 Å⁻¹ |
| Ir | `_standards_lib/Ir/Ir_metal_foil.dat`（`i/irfoie.a46`） | `959da06b723b` | WT 无有效网格 | 5.66 Å⁻¹ |
| Ce | `data/bmm_standards/Data/Ce/Ce-L3-CeO2.xdi` | `1f679ebd0fb3` | ΔE₀ 贴界/σ² 贴界 | 10.58 Å⁻¹ |
| U | `_standards_lib/U/U_metal_foil.dat`（`U/umetal.351`） | `5b31461e25fe` | α-U 第一壳模型不稳 | 20.47 Å⁻¹ |

这些旧谱的完整 sha256、点数、能量范围、`gate_failures` 与
`failure_reason`/`failure_error` 都在 `index.json` 的
`quarantined_legacy_source` 中，替换不会丢失审计链。

## 5. 工作流

```text
T5 失败/拒绝
   │
   ├─→ quarantined：登记谱源、sha256、失败证据、分类、动作
   │
   ├─→ candidate：新谱进入候选（不产生 S₀²）
   │       │
   │       ├─ 未通过冻结门禁 → quarantined（记录原因）
   │       └─ 通过冻结门禁 → resolved（写 calibration.json + index.json）
   │
   └─ resolved：resolved_by 记录新谱；旧谱移入 quarantined_legacy_source
```

规则：

- `quarantined` 不等于“删除”；原始谱仍保留用于审计/复现；
- 旧谱即使被替换也继续留在池中，记录“已拒绝/已替换”原因；
- 任何状态变更都必须经过 `tools/build_source_pool.py --check`；
- 池内谱源不进入生产 `_structure_lib`，不进入运行期 `--s02-from-t5`
  之外的生产路径（`--s02-from-t5` 只读已提交的 `calibration.json`）。

## 6. 替换谱接收清单

新谱进入 `candidate` 前必须提供：

1. **来源**：实验谱，明确 beamline/文献/DOI；不得用 DFT 或模拟谱；
2. **许可**：再分发许可明确；许可不明时只保留本地，不提交；
3. **指纹**：源文件 sha256、能量范围、点数、E₀、边；
4. **长度**：EXAFS k_max ≥ 10 Å⁻¹（优先 ≥ 12 Å⁻¹）；三窗口必须有
   足够 `N_idp`；
5. **物相匹配**：标准谱物相与 motif 一致；换相必须同时换
   calibration-only motif，不能“谱是氧化物、模型是金属”；
6. **归一化**：边后平台 0.85–1.15（宽松模式需记录 plateau 且五项拟合
   门禁仍须通过）；边前无明显伪迹；
7. **门禁**：原样通过第 2 节全部条件；
8. **记录**：保留失败尝试与最终通过值，不覆盖历史。

## 7. 不变式与禁止项

- `absorber.py` 静态表：Fe/Co/Ni/Cu/Pt 已标定，其余 12 元素
  `s02=None`、`calibrated=False`；池不改变这一点。
- T5 结果只落 `data/s02_calibration/`，运行时必须显式
  `--s02` 或对已通过元素使用 `--s02-from-t5`。
- **禁止**：放宽 `S02_PASS`、`S02_ERR_MAX`、`RF_MAX`、`E0_MAX`、
  `WINDOW_SPREAD_MAX`；禁止手工裁剪窗口凑 interior；禁止自动 Top-1；
  禁止把隔离谱源写成通过。
- `tools/build_source_pool.py --check` 失败即表示索引与标定产物/门禁
  常量不一致，必须先修正再继续。

## 8. 复现与相关文件

```powershell
python tools/build_source_pool.py          # 生成/更新 index.json
python tools/build_source_pool.py --check  # 校验索引与标定产物一致
python tools/build_source_pool.py --json   # 打印完整索引
python regression_test.py                  # 契约 #126 锁定池；
                                           # #127 锁定 P5-D 多位点
```

相关文件：

- `data/s02_calibration/index.json`：T5 总索引（12/12 calibrated）
- `data/s02_calibration/{Rh,Pd,Ir,Ce,U}/calibration.json`：通过/失败尝试明细
- `data/s02_calibration/REPORT.md`：T5 报告
- `calibrate_s02.py`：门禁常量与生产等效拟合协议
- `calibrate_t5.py`：5 元素新谱 CASES
- `cif_crawler.py`：IIT-Lytle RAW 解析器（Rh/Pd/Ir/Ce/U）
- `data/s02_calibration/build_motifs.py`：calibration-only motif 生成器
- `_structure_lib/_phase_pack/`：P5-C 补齐的常见相 CIF（UO2 等）
