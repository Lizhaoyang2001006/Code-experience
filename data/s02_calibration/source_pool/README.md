# S₀² 标准谱源隔离池（P5 遗留项）

日期：2026-11-09  
范围：T5 未通过的 5 个元素——**Rh、Pd、Ir、Ce、U**  
产物：`data/s02_calibration/source_pool/index.json`（机器可读）  
生成/校验：`python tools/build_source_pool.py [--check]`

## 1. 目的与边界

这个池只做三件事：

1. 把未通过 T5 的谱源**逻辑隔离**，明确标记为不可用于 S₀²；
2. 记录每个元素的失败证据、建议分类、下一步动作与替换准入条件；
3. 冻结门禁，防止“为了通过而改阈值/手调窗口/自动选 Top-1”。

边界：

- 池内谱源**不产生 S₀²**，`s02` 永远为 `null`，`usable_for_s02=false`；
- 不移动、不删除原始谱，避免破坏既有标定溯源；隔离状态只存在于索引；
- 不修改 `absorber.py` 静态表；T5 结果仍只落独立产物；
- 运行期不能用池内谱源做 `--s02-from-t5`；
- 新谱只有通过**原样冻结**的全部门禁才可移出隔离池。

## 2. 冻结门禁（不得放宽）

来源：`calibrate_s02.py` 模块级常量，由 `tools/build_source_pool.py` 用 AST 读取，
不导入/不执行标定模块。

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

**这些值在隔离池建立后保持不变。** 任何“调 k_min / 换 rmin / 放宽 ΔE₀”的建议都必须
作为诊断项记录，不能作为通过判定的替代。

## 3. 状态分类总表

| 元素 | 当前谱源 | 边/E₀ | k_max 估计 | 失败类别 | 隔离状态 | 下一步 |
|------|----------|-------|-----------|----------|----------|--------|
| Rh | IIT-Lytle Rh 箔 `r/rhfoil.a03` | K / 23220 eV | 5.67 Å⁻¹ | 谱截断 | quarantined | 换长谱 |
| Pd | IIT-Lytle Pd 箔 `p/pdfoil.a71` | K / 24350 eV | 5.68 Å⁻¹ | 谱截断 | quarantined | 换长谱 |
| Ir | IIT-Lytle Ir 箔 `i/irfoie.a46` | L3 / 11215 eV | 5.66 Å⁻¹ | 谱截断 | quarantined | 换长谱，再做 k_min 诊断 |
| Ce | NIST BMM `Ce-L3-CeO2.xdi` | L3 / 5726 eV | 10.58 Å⁻¹ | 归一化/ΔE₀/σ² | quarantined | 保持 CeO2 基准重处理 |
| U | IIT-Lytle U 金属 `U/umetal.351` | L3 / 17163 eV | 20.47 Å⁻¹ | 物相/模型不适配 | quarantined | 核实 UO2 谱源后再切换 |

> 注：早期 `calibration.json` 对 Rh/Pd/Ir 写的“170–195 eV above E0”是近似描述；
> 按源文件实测，三者边后跨度均约 122 eV，`k_max≈5.66–5.68 Å⁻¹`。隔离池索引
> 采用实测值，`failure_reason` 保留标定产物原文以便审计。

### 与 Gemini 建议的逐项对照

| 元素 | Gemini 建议分类 | 本项目评估 |
|------|----------------|------------|
| Rh | XANES 截断待换谱 | **同意**；需换 k_max ≥ 10–12 Å⁻¹ 的长谱 |
| Pd | XANES 截断待换谱 | **同意**；需换 k_max ≥ 10–12 Å⁻¹ 的长谱 |
| Ir | 参数窗口调优（提高 k_min 避开白线） | **不同意仅靠调 k_min**；当前 k_max≈5.66 Å⁻¹，k_min≥2.5 时 N_idp 仍小于 4 个参数，且失败发生在 WT 网格为空；必须先换谱，换谱后 k_min 只能作为诊断项 |
| Ce | 物相再标定，转向 CeO2 | **物相已经是 CeO2**；主谱是 BMM `Ce-L3-CeO2.xdi`，motif 是 `ceo2_fluorite_cod_1562989`；真正问题是 ΔE₀ 贴 10 eV、σ² 压到 1e-4 下界、无 interior 窗口 |
| U | 核实是否换 UO2 | **同意核实，但不能直接换**；本地没有 UO2 标准谱，只有 RAW U 金属荧光谱；UO2 结构 CIF 已在 P5-C 入库，需先找到可追溯的 UO2 L3 谱 |

## 4. 逐元素状态卡

### 4.1 Rh：谱截断，换谱

- 谱源：`_standards_lib/Rh/Rh_metal_foil.dat`（IIT-Lytle，`r/rhfoil.a03`）；
  sha256 `074792a1…a9d85`；未提交。
- 数据：137 点，23169.73–23342.59 eV，E₀=23220 eV；
  边后仅 122.6 eV，k_max≈5.67 Å⁻¹。
- 失败：`ValueError: zero-size array to reduction operation maximum which has no identity`；
  WT 无有效网格，未进入拟合。
- 动作：获取 Rh K 边、k_max ≥ 10–12 Å⁻¹ 的 Rh 箔或 Rh2O3 谱；
  优先 Rh 箔（fcc，N=12）复核 S₀²；Rh2O3 需配套 corundum motif。
- 准入：原样通过第 2 节全部门禁。

### 4.2 Pd：谱截断，换谱

- 谱源：`_standards_lib/Pd/Pd_metal_foil.dat`（IIT-Lytle，`p/pdfoil.a71`）；
  sha256 `9f07525a…7a21e`；未提交。
- 数据：195 点，24275.21–24472.85 eV，E₀=24350 eV；
  边后仅 122.9 eV，k_max≈5.68 Å⁻¹。
- 失败：同 Rh，WT 无有效网格。
- 动作：获取 Pd K 边、k_max ≥ 10–12 Å⁻¹ 的 Pd 箔或 PdO 谱；
  优先 Pd 箔（fcc，N=12）；PdO 需配套 tetragonal motif。
- 准入：原样通过第 2 节全部门禁。

### 4.3 Ir：谱截断；k_min 调优不能替代换谱

- 谱源：`_standards_lib/Ir/Ir_metal_foil.dat`（IIT-Lytle，`i/irfoie.a46`）；
  sha256 `959da06b…d341b`；未提交。
- 数据：272 点，11142.70–11336.90 eV，E₀=11215 eV；
  边后仅 121.9 eV，k_max≈5.66 Å⁻¹。
- 失败：同 Rh/Pd，WT 无有效网格。
- 对“提高 k_min 避开白线”的评估：
  - 当前失败点是 WT 网格为空，不是白线污染；
  - 即使 `k_min=2.5 Å⁻¹`，`Δk≈3.16 Å⁻¹`；3 窗口的 `N_idp` 仍小于
    `s02 + ΔE₀ + ΔR + σ²` 四个自由参数；
  - 因此 k_min 调优**不能**救当前谱，也不能用来绕过门禁；
  - 换到 k_max ≥ 10–12 Å⁻¹ 的 Ir 谱后，可将提高 k_min 作为白线规避的
    诊断项，并记录 `kw=1/2/3` 敏感性。
- 动作：优先换 Ir L3 长谱（Ir 箔或 IrO2）；再做 k_min 诊断；
  三窗口 spread、N_idp、σ² 不贴界等门禁不变。

### 4.4 Ce：物相已是 CeO2，问题在归一化/ΔE₀/σ²

- 主谱：`data/bmm_standards/Data/Ce/Ce-L3-CeO2.xdi`（NIST BMM 公共数据，
  DOI 10.18434/mds2-4032）；sha256 `1f679ebd…c1e7e`；已提交。
- 结构：CeO2 fluorite，COD 1562989；第一壳 Ce–O 2.3414 Å，N=8。
- 数据：259 点，5523.405–6152.412 eV，E₀=5726 eV；k_max≈10.58 Å⁻¹。
  **谱长足够**，不是截断问题。
- 失败证据：
  - 三窗口 ΔE₀=9.78–10.20 eV，贴 `E0_MAX=10`；
  - σ² 三窗口均为 `1e-4`，`sigma2_at_bound`；
  - 无 interior 窗口，`s02_in_range`/`s02_err_ok`/`window_spread_ok` 失败；
  - 单窗口 S₀²≈0.46 只是数值解，不是通过值。
- 已拒绝的替代谱：`Ce-L3-CeriumDioxide.xdi`（sha256 `e4d10395…db6603`），
  边后平台 2.271，属归一化伪迹，严格门禁必须拒绝。
- 动作：保持 CeO2 基准，复核 rbkg/背景、E0 方法（derivative vs
  half_height）、可比较其它来源的 CeO2 L3 谱；不得放宽 ΔE₀/σ² 门禁；
  若仍无唯一解，保持 `not_calibrated`。
- 注意：`data/standards_raw/Ce_metal_foil.dat` 的文件头实际写的是
  “CeO2 L3 参考谱（IIT-Lytle 库中无 Ce 金属箔）”，不要被文件名误导。

### 4.5 U：先核实 UO2 谱源，再决定是否换相

- 当前谱：`_standards_lib/U/U_metal_foil.dat`（IIT-Lytle RAW 荧光，
  `U/umetal.351`）；sha256 `5b31461e…10dfe`；未提交。
- 数据：455 点，16966.00–18759.65 eV，E₀=17163.1 eV；k_max≈20.47 Å⁻¹。
  **谱长足够**，问题不是截断。
- 结构：α-U，COD 1527120；第一壳 U–U 2.7425 Å，固定 N=4。
- 失败证据：
  - S₀² 钉在 `0.3` 下界；
  - ΔE₀=9.41–10.94 eV；
  - Rf=0.55–0.60；
  - `N_idp=2.48 < 4` 个自由参数；
  - 5 项门禁全部失败。
- UO2 候选的物理依据：UO2 fluorite 第一壳为 8×U–O≈2.37 Å，
  比 α-U 金属的多壳层/强相关模型更适合作为 S₀² 标定基准。
- 现实约束：本地没有 UO2 标准谱；BMM 数据没有 U 目录；P5-C 只入库了
  UO2 结构 CIF（COD 9009049），没有配套谱。
- 动作：
  1. 先核实是否存在可追溯、许可明确的 UO2 L3 标准谱；
  2. 若有：用 COD 9009049 生成 **calibration-only** UO2 fluorite motif，
     以 8×U–O 第一壳重新标定；
  3. 若无：保持 `not_calibrated`，不用 RAW U 金属谱强行给值；
  4. 换相后仍原样通过第 2 节全部门禁。

## 5. 隔离池工作流

```text
T5 失败/拒绝
   │
   ├─→ quarantined：登记谱源、sha256、失败证据、分类、动作
   │
   ├─→ candidate：新谱进入候选（不产生 S₀²）
   │       │
   │       ├─ 未通过冻结门禁 → quarantined（记录原因）
   │       └─ 通过冻结门禁 → accepted
   │
   └─ accepted：更新 data/s02_calibration/<El>/calibration.json 与 index.json
```

规则：

- `quarantined` 不等于“删除”；原始谱仍保留用于审计/复现；
- 旧谱即使被替换也继续留在池中，记录“已拒绝/已替换”原因；
- 任何状态变更都必须经过 `tools/build_source_pool.py --check`；
- 池内谱源不进入生产 `_structure_lib`，不进入运行期 `--s02-from-t5`。

## 6. 替换谱接收清单

新谱进入 `candidate` 前必须提供：

1. **来源**：实验谱，明确 beamline/文献/DOI；不得用 DFT 或模拟谱；
2. **许可**：再分发许可明确；许可不明时只保留本地，不提交；
3. **指纹**：源文件 sha256、能量范围、点数、E₀、边；
4. **长度**：EXAFS k_max ≥ 10 Å⁻¹（优先 ≥ 12 Å⁻¹）；三窗口必须有
   足够 `N_idp`；
5. **物相匹配**：标准谱物相与 motif 一致；换相必须同时换
   calibration-only motif，不能“谱是氧化物、模型是金属”；
6. **归一化**：边后平台 0.85–1.15；边前无明显伪迹；
7. **门禁**：原样通过第 2 节全部条件；
8. **记录**：保留失败尝试与最终通过值，不覆盖历史。

## 7. 不变式与禁止项

- `absorber.py` 静态表：Fe/Co/Ni/Cu/Pt 已标定，其余 12 元素
  `s02=None`、`calibrated=False`；隔离池不改变这一点。
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
python regression_test.py                  # 契约 #126 锁定隔离池
```

相关文件：

- `data/s02_calibration/index.json`：T5 总索引（5 个 `not_calibrated`）
- `data/s02_calibration/{Rh,Pd,Ir,Ce,U}/calibration.json`：失败尝试明细
- `data/s02_calibration/REPORT.md`：T5 报告
- `calibrate_s02.py`：门禁常量与生产等效拟合协议
- `data/s02_calibration/build_motifs.py`：calibration-only motif 生成器
- `_structure_lib/_phase_pack/`：P5-C 补齐的常见相 CIF（UO2 等）
