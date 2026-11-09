# P4+ 执行计划：催化材料分析能力补齐

基于 P0a–P4 已完成（契约 114/4，零退化），本计划补齐催化材料实际定量分析所需的关键能力。

## 现状缺口

| 缺口 | 说明 | 对催化分析的影响 |
|------|------|----------------|
| G1 真实 CIF 稀缺 | 已整合构建期爬取器：phase pack 从 5 条金属种子扩展到 26 条真实 COD CIF（5 种子 + 21 抓取目标） | 氧化物/氢氧化物/金属基准相可产生正式 CN |
| G2 缺失元素覆盖 | Cr/Mn/Zn/Mo/Ru/Rh/Pd/Ag/Ir/Au/Ce/U 已全部补齐，phase pack 覆盖 17/17 元素 | 所有 17 个吸收元素都有 verified CIF |
| G3 S₀² 未标定 | Cr/Mn/Zn/Mo/Ru/Rh/Pd/Ag/Ir/Au/Ce/U 12 个元素 s02=None | 即使有 CIF，validate_formal 仍拒绝 CN 定量 |
| G4 P3 未端到端 | build_conditional_fit_args 仅构建参数，未接入 artemis.fit_r_space 固定 N 拟合 + 自由 N 对照 + Hamilton F-test | 正式条件 CN 无法自动触发 |
| G5 多元素联合 | 当前按单吸收元素独立处理 | 多金属催化剂需逐元素手动综合 |
| G6 XANES-EXAFS 联合 | 白线面积/边步与配位数无一致性校验 | 结构鉴别与定量断层 |

## 执行计划（按优先级）

| ID | 优先级 | 任务 | 预期产出 | 验收标准 | 状态 |
|----|--------|------|----------|----------|------|
| T1 | P0 | 整理现有 `_structure_lib/*/_cif/` 生成 CIF 目录，列出已有氧化物/标准相 | CIF 清单与归档 | 清单完整，各 CIF 可被 gemmi 解析 | 已完成 |
| T2 | P0 | 整理 `data/` 目录下标准样数据，整理为可用于 S₀² 标定的样本列表 | 标准样清单 | 清单含元素、边、文件路径、是否有 prj | 已完成 |
| T3 | P1 | 真实 CIF（催化常见氧化物/氢氧化物）导入 phase pack | 9 条 T3：Fe₂O₃ COD:1011240、Fe₃O₄ COD:1010369、CoO COD:1533087、Co₃O₄ COD:9005888、NiO COD:1010093、Ni(OH)₂ COD:1548811、CuO COD:1011148、Cu₂O COD:1010941、CeO₂ COD:1562989 | 每条 verified_source_cif，sha256 匹配，空间群+化学计量核验 | 已完成 |
| T4 | P1 | 扩展 phase pack 覆盖缺失元素 | 12 条 T4 金属：Mn COD:9008589、Ru COD:9008513、Rh COD:9008482、Pd COD:9008478、Ag COD:9008459、Ir COD:9008470、Ce COD:9010987、U COD:1527120、Cr COD:9008531、Zn COD:9008522、Mo COD:9008543、Au COD:9008463 | phase_coverage.json verified_elements=17/17 | 已完成 |
| T5 | P2 | S₀² 标定：对 12 个未标定元素跑批量标定；结果只落独立产物（`data/s02_calibration/`），不改 `absorber.py` 静态表 | `calibrate_t5.py` + 逐元素 `calibration.json` + `index.json` + `REPORT.md`；运行时用 `--s02` 显式覆盖，已通过元素可用 `--s02-from-t5` | 12 元素全部通过（12/12）：Cr/Mn/Zn/Mo/Ru/Ag/Au/Rh/Pd/Ir/Ce/U；`#60` 静态 12 元素仍为 `s02=None` | 已完成（P5-S02 换谱后 12/12；独立产物 + `--s02-from-t5`） |
| T6 | P2 | P3 端到端接入：把 `build_conditional_fit_args`/`compute_n_ratio`/`run_conditional_cn` 接入 `xafs_fit.fit_identified`，实现固定 N 拟合、自由 N 对照、Hamilton F-test | xafs_fit.py 新增 `run_conditional_cn()`；`fit_identified`/`run` 接受 `gate_result`；真实 Fe 箔 `N/N₀=0.726 (ok)`，正式表 2 行 | 固定 N 拟合成功；N/N₀ 计算正确；Hamilton F-test 可区分对错相 | 已完成（真实 Fe 箔 COD:1100108 端到端） |
| T7 | P3 | 多元素联合判据：多吸收元素样本自动分元素跑批 + 综合报告 | `run_xafs.py --multi-element` / `--multi-element-s02`；`multi_element_summary.json/.md`；逐元素独立子目录 | 多元素样本可一键处理；报告并列各元素 CN 对比；`recommended=null`，不自动 Top-1 | 已完成（#120） |
| T8 | P3 | XANES-EXAFS 联合校验：白线面积/边步与配位数一致性检查 | `xafs_fit.xanes_exafs_cross_check()`；`run_xafs.py --xanes-cross-check`；`<stem>_meta.json.xanes_cross_check` | 白线/边步/归一化异常时警告；`affects_cn_reliable=False`，不改正式 CN 门禁；可选开关 | 已完成（#121） |
| T9 | P4 | 更新 ACCEPTANCE_LOG.md §9（P4+），导出 patches/P4+.patch | 验收记录完整；补丁可反向应用 | run_ci.py --stress 全绿；契约数 ≥119 | 已完成（提交 `d7fa8c8`；补丁 SHA256 `9B9779EBD39D3C493EF3D877797ADB9E1E0D33EA1982DAD302F4239524B5FF10`；反向检查通过） |
| T10 | P1 | 整合构建期 CIF 爬取器到主项目 | `tools/cif_crawler.py`：固定 COD ID、IT 空间群+化学计量核验、原始 CIF 缓存、phase pack 导入、索引/覆盖重建；运行期不联网 | `--list` 21 目标；`--all --offline` 幂等；错相/离线缺失拒绝；契约 #115–#118 通过 | 已完成 |
| T11 | P5-A | T5 工程加固与运行期溯源 | `calibrate_t5.py --only` 增量合并；`run_xafs.py --s02-from-t5`；L1/motif/meta 记录 `s02_source`；契约 #122/#123 | `--only` 不丢其他元素；只接受 `status=calibrated` + `static_absorber_py_unchanged=True`；未标定元素明确拒绝；显式 `--s02`/`--multi-element-s02` 优先；T5 产物 `s02_calibrated=True` | 已完成 |
| T12 | P5-B | Mo/Ce 标定物理排查（不放宽门禁） | `calibrate_s02` 新增 `e0_method`/`manual_e0`/`fit_shells`；Mo 箔双壳 bcc；`sigma2_not_at_bound` 门禁；契约 #124 | Mo S₀²=0.8322、Rf=0.0040、ΔE₀=−0.73 eV、5/5 门禁；Ce 保持 `not_calibrated` 并记录 σ² 压界/无唯一解/伪迹拒绝原因；静态表仍 12 元素 `None` | 已完成 |
| T13 | P5-C | 固定 COD ID 补齐 18 个常见相 | `tools/cif_crawler.py` 新增 tier=`P5C` 与 `--p5c`；phase pack 26→44；39 个抓取目标；契约 #125 | 18 个 `missing_common_phases` 全部清空；每条 experimental `crystal_cif`、IT 空间群/化学计量核验；多位点仍需显式 `--absorber-site` | 已完成 |
| T14 | P5-prep | 标准谱源隔离池（为 P5 遗留项准备） | `tools/build_source_pool.py`；`data/s02_calibration/source_pool/{README.md,index.json}`；契约 #126 | Rh/Pd/Ir/Ce/U 状态分类、谱源 sha256/k_max/失败证据齐备；`s02=None`、`usable_for_s02=false`；门禁阈值冻结；Ir 的 k_min 建议仅作换谱后诊断 | 已完成 |
| T15 | P5-S02 | 隔离池 5 元素换谱闭环 | IIT-Lytle RAW 长谱解析（Rh/Pd/Ir/Ce/U）；`calibrate_t5.py` CASES 更新；`source_pool/index.json` schema 2（resolved） | 5 元素全部通过冻结门禁，T5 7/12 → 12/12；旧谱保留在 `quarantined_legacy_source`；门禁零放宽 | 已完成 |
| T16 | P5-D | 多位点正式 CN（`--absorber-site all`） | `phase_resolver` 展开每位点模型；`cif_gate.gate_cif_all_sites` 化学计量权重；`xafs_fit.run_multi_site_formal_cn`；CLI 集成；契约 #127 | Fe3O4 三位点各权重 1/3、加权 N₀=16/3；Co3O4 8a:16d=1/3:2/3；未传 site 仍拒绝；任一位点失败不出 ALL_WEIGHTED 汇总行 | 已完成 |

## 约束（继承 P0a–P4）

- 只改副本 `Code-experience-cn`；原目录 `Code-experience` 不动
- 禁止 `git checkout --.`、`git reset --hard`、`git clean -fd`
- 正式 CN 只认 `verified_source_cif`；generated_only 不得进正式表
- 禁止 `eval`/`exec`；`--phase SPEC` 只做字符串拆分
- Python 用 `C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe`
- 非 ASCII 输出先 `sys.stdout.reconfigure(encoding="utf-8")`
- 每批独立 `run_ci.py --stress` 全绿，独立回滚、独立验收

## 约定变更

- 测试账本继续累加：P4+ 契约从 114 增至 121（新增 #115–#121）；P5-A 增至 123（新增 #122/#123）；P5-B 增至 124（新增 #124）；P5-C 增至 125（新增 #125）；P5-prep 增至 126（新增 #126）；P5-S02/P5-D 增至 127（新增 #127）
- 合格结构来源：真实 CIF + sha256 + 吸收位；不自动 Top-1；L1 不进正式表
- 多吸收位相（Fe3O4/Co3O4/α-Mn）必须显式 `--absorber-site`，不自动取第一个位点；`--absorber-site all` 表示逐位点正式 CN + 化学计量加权汇总
- phase pack 现为 44 条：5 条金属种子 + 39 条爬取目标；覆盖 17/17 元素
- 条件 CN 通过才产生正式行；远壳严格固定 N₀，不与首壳绑定；正式行含数据哈希

## T5 口径修订（#60 优先，2026-11-09）

- 原验收“s02 计数从 5 增至 17”与 `EXECUTION_MASTER_PLAN.md` 的 #60 硬约束冲突：
  #60 要求 `absorber.py` 静态表中 Fe/Co/Ni/Cu/Pt 之外 12 元素保持
  `s02=None`/`calibrated=False`。
- 因此 T5 的实际产物改为**独立标定产物**：`data/s02_calibration/index.json`、
  逐元素 `calibration.json`、`REPORT.md`；不写 `absorber.py`。
- 运行时通过 `run_xafs.py --s02 <值>` 使用标定值，已通过元素也可用
  `run_xafs.py --s02-from-t5` 从独立产物读取；未给 S₀² 时未标定元素
  仍拒绝 CN 定量。
- 当前标定通过 12/12：Cr 0.6482、Mn 0.4311、Zn 0.6071、Mo 0.8322、
  Ru 0.7002、Rh 0.7022、Pd 0.6559、Ag 0.7513、Ir 0.4752、Au 0.8602、
  Ce 0.6279、U 0.7800。Mo 由 P5-B 通过（双壳 bcc + `half_height` E₀）；
  Rh/Pd/Ir/Ce/U 由 P5-S02 换 IIT-Lytle RAW 长谱后通过（冻结门禁零放宽）。
- IIT-Lytle 标准许可未明确，`_standards_lib/` 与 `data/standards_raw/` 不提交；
  相关结果记录 SHA256，但不保证无本地标准时可复现。
- 标定协议统一 larch 归一化；BMM RuO2 的 `xmu` 列符号为负，使用 `-xmu`
  派生谱并在 provenance 中记录。

## P5-A T5 运行期溯源（2026-11-09）

- `calibrate_t5.py --only <元素...>` 只重跑所选元素；其余元素从既有
  `<元素>/calibration.json` 合并回 `index.json`/`REPORT.md`，不把汇总
  覆盖成子集。全量运行仍重写全部 12 元素。
- `run_xafs.py --s02-from-t5 --identify/--fit` 读取逐元素产物：
  - 只接受 `status="calibrated"` 且 `static_absorber_py_unchanged=True`；
  - 未通过元素抛 `XAFSConfigError`（含“未通过”），不静默回退；
  - 优先级 `--multi-element-s02` > 全局 `--s02` > T5 产物；
  - 显式 `--absorber` 优先于文件头元素判定，避免 S₀² 与拟合元素错配。
- 溯源串格式：
  `s02_calibration:<相对路径>#sha256:<12位>;standard_sha256:<12位>`；
  写入 L1/motif 输出与 `<stem>_meta.json.fit.s02_source`。
- 正式 CN 表：T5 产物来源记 `s02_calibrated=True`；普通 `--s02` 手填值
  仍为 `False`；`absorber.py` 静态表不变（#60）。
- Cr₂O₃ 真实谱冒烟：`fit.s02=0.6482`，`s02_source` 以
  `s02_calibration:` 开头；L1 `N=6.34±0.21`、`R=1.975 Å`、`Rf=0.0216`。
- 契约 #122（产物完整性/静态锁/`--only` 合并）与 #123（加载器/来源
  优先级）通过；P5-A 时完整回归 123 项；P5-B 后 124 项；P5-C 后 125 项、
  CI 125+4 PASS（215.1s；完整回归 142.474s）；P5-prep 后 126 项、
  CI 126+4 PASS（200.4s；完整回归 176.600s）；P5-S02/P5-D 后 127 项、
  CI 127+4 PASS。

## P5-B Mo 双壳 bcc 与 Ce 物理排查（2026-11-09）

- `calibrate_s02.calibrate_s02` 新增 `e0_method`/`manual_e0`/
  `fit_shells`；CLI 增加 `--e0-method`、`--manual-e0`、`--fit-shells`。
- 新增门禁 `sigma2_not_at_bound`：任一窗口 σ² 压到 1e-4 下界即失败。
  这是更严格的物理约束，不放松原有阈值。
- **Mo 通过**：BMM Mo 箔导数峰 20018.8 eV 落在近边肩峰，平台 0.842；
  改用 `half_height` E₀=20012.86 eV 后平台 0.863。bcc 第一壳
  8×2.725 Å + 第二壳 6×3.147 Å 同属 Mo–Mo；单壳会把第二壳尾巴折进
  第一壳（S₀² 撞 1.3 上界），双壳共享 S₀²/ΔE₀、各壳独立 ΔR/σ²：
  S₀²=0.8322±0.0334、Rf=0.0040、窗口极差 0.0017、ΔE₀=−0.73 eV、
  n_var=6 < N_idp=6.88，5/5 门禁通过。
- **Ce 仍未通过**：自由 σ² 在窄窗内塌到 1e-4；固定 σ² 时 S₀² 随假设
  σ² 从 0.54 单调漂到 0.80，无唯一解；加 Ce–Ce 3.82 Å 第二壳后
  Rf=0.127（>0.12）、ΔR 撞界；备选 `Ce-L3-CeriumDioxide.xdi` 透射谱
  平台 2.27、边步被背景淹没，属归一化伪迹，严格门禁拒绝。
- 契约 #124 锁定 Mo 双壳产物与 Ce 拒绝路径。

## P5-C 固定 COD ID 补齐常见相（2026-11-09）

- `tools/cif_crawler.py` 新增 tier=`P5C` 与 `--p5c`，按固定 COD ID 补齐
  `phase_pack.COMMON_PHASES` 的 18 个 `missing_common_phases`：
  FeO 1011169、FeOOH 9002158、PtO2 1008935、Cr2O3 9014036、
  MnO 9006658、MnO2 9009081、Mn3O4 1514120、ZnO 9004178、
  MoO2 9009090、MoO3 9009669、RuO2 9007541、Rh2O3 1010584、
  PdO 1009031、Ag2O 1010486、IrO2 2101854、Au2O3 2105390、
  UO2 9009049、U3O8 2310519。
- phase pack 由 26 条增至 44 条；固定抓取目标由 21 个增至 39 个；
  `phase_coverage.json` 的 `missing_common_phases` 全部为空。
- 选择规则：优先常压/室温实验结构；高压、高温或 DFT 条目排除
  （Cr2O3 的 63.9 GPa 条目、β-U3O8、高温 U3O8 均未采用）。
  PtO2 在 COD 中没有常压结构，只能采用高压合成 β-PtO2，已在
  phase/notes 明确标注。
- Mn3O4（Mn1/Mn2）与 U3O8（U1/U2）有多个吸收位点，正式 CN 仍必须
  显式 `--absorber-site`，不自动 Top-1。
- 契约 #125 锁定 18 个目标、44 条索引与 coverage 无缺口。
- `data/cif_raw/` 仍只保留 17 个根索引 CIF；P5-C 的 18 个下载缓存
  只用于构建期导入，导入后清理，不作为运行期来源。

## P5 遗留项准备：标准谱源隔离池（2026-11-09）

- 新增 `tools/build_source_pool.py`：从 T5 标定产物与
  `calibrate_s02.py` 门禁常量确定性生成
  `data/s02_calibration/source_pool/index.json`。
- 新增 `data/s02_calibration/source_pool/README.md`：5 个未通过元素
  （Rh/Pd/Ir/Ce/U）的状态分类、失败证据、下一步动作与替换准入。
- 逻辑隔离：不移动/删除 IIT-Lytle 原始谱；池内元素 `s02=None`、
  `usable_for_s02=false`，运行期不得用于 `--s02-from-t5`。
- 分类判定：
  - Rh/Pd：短谱（k_max≈5.7 Å⁻¹），需换 k_max≥10–12 Å⁻¹ 的长谱；
  - Ir：同为短谱；Gemini 的“提高 k_min”建议**不足以**救当前谱，
    只保留为换谱后的诊断项；
  - Ce：已经是 CeO2 基准（BMM `Ce-L3-CeO2.xdi`），问题在 ΔE₀ 贴
    10 eV、σ² 压到 1e-4、无 interior 窗口；
  - U：谱长足够但 α-U 第一壳模型不稳；UO2 切换需先核实可追溯的
    UO2 L3 标准谱（本地暂无）。
- 门禁冻结：`S02_PASS=[0.4,1.15]`、`s02_err<0.06`、`Rf≤0.12`、
  `|ΔE₀|<10 eV`、三窗口 interior、跨窗极差≤0.05、σ² 不贴界、
  `vars<N_idp`；禁止改门禁/手调窗口/自动 Top-1。
- 契约 #126 锁定 5 元素分类、`s02=None`、门禁常量与 T5 失败状态。

## P5-S02 隔离池换谱闭环（2026-11-09）

- `cif_crawler.py` 新增 IIT-Lytle RAW 解析器与谱源：
  - Rh `RAW/Rh/rhfoil.003`（514 点，22694–24807 eV，透射）；
  - Pd `RAW/Pd/pdfl.404`（520 点，24156–26353 eV，NSLS 透射，扣 `OFFSETS`）；
  - Ir `RAW/Ir/irfoil.146`（645 点，11015–12820 eV，步进扫描；该谱
    `ln(I0/I1)` 为负边步，解析器按 `-ln(I0/I1)` 恢复 μ 边后升高约定，
    provenance 记录符号约定）；
  - Ce `RAW/Ce/ceo2.004`（327 点，5523–6158 eV，在 Ce L2 6164 eV 前截断）；
  - U `RAW/U/uo2.379`（381 点，16966–18170 eV）。
- `calibrate_t5.py` 更新 5 元素 CASES；U 使用标定专用 `uo2_fluorite`
  motif（COD 9009049，8×U–O 2.3678 Å），只写
  `data/s02_calibration/motifs/`，不进生产 `_structure_lib`。
- 结果为 T5 **12/12**（冻结门禁零放宽）：

  | 元素 | S₀² | err | Rf | ΔE₀ (eV) | spread | plateau |
  |------|-----|-----|----|-----------|--------|---------|
  | Rh | 0.7022 | 0.0191 | 0.0041 | −0.71 | 0.0037 | 0.8590 |
  | Pd | 0.6559 | 0.0141 | 0.0035 | +1.14 | 0.0133 | 0.8138 |
  | Ir | 0.4752 | 0.0149 | 0.0043 | +9.30 | 0.0068 | 0.8604 |
  | Ce | 0.6279 | 0.0268 | 0.0223 | +5.94 | 0.0096 | 1.0645 |
  | U  | 0.7800 | 0.0229 | 0.0137 | +2.86 | 0.0024 | 0.8608 |

- `source_pool/index.json` 升级 schema 2：5 元素
  `pool_state="resolved"`、`usable_for_s02=true`、`pool_state_counts=
  {"resolved":5}`；旧谱 sha256/点数/失败门禁保留在
  `quarantined_legacy_source`，`resolved_by` 记录新谱、motif 与门禁结果。
- 隔离池索引的数值全部来自已提交 `calibration.json` + 固定观测常量，
  不依赖本机 `_standards_lib`；`tools/build_source_pool.py --check` 通过。
- 契约 #122/#123/#124/#126 更新到 12/12 与新谱数值。

## P5-D 多位点正式 CN（2026-11-09）

- `phase_resolver.PhaseResolution` 新增 `site_models`；
  `resolve_phase(..., absorber_site="all")` 把多位点候选展开为每个
  不等价吸收位一个 `StructureModel`。未传 `--absorber-site` 时仍返回
  `structure_model=None` + `absorber_site_ambiguous` 原因，绝不自动 Top-1。
- `cif_gate` 新增 `list_inequivalent_absorber_sites` 与
  `gate_cif_all_sites`：每位点独立过同一套门禁，并计算化学计量权重
  `w_i = m_i·occ_i / Σ(m_j·occ_j)`。
  - Fe3O4（COD 1010369）：Fe1/Fe2 各 mult=16、occ=0.5，Fe3 mult=8、
    occ=1 → 权重各 1/3；第一壳 N₀ = 6、6、4，加权 N₀ = 16/3 = 5.3333。
  - Co3O4（COD 9005888）：Co1(8a):Co2(16d) = 1/3 : 2/3；但该 COD CIF
    缺少 `_atom_site_occupancy`，冻结门禁按 `occupancy_missing` 拒绝，
    不出正式 CN（无合格来源不出 CN）。
- `xafs_fit.run_multi_site_formal_cn`：逐位点构建团簇 motif → 条件 CN →
  `build_formal_cn_block`；每条成员行带 `multi_site_member` 与
  `site_weight`。**仅当全部位点 ok** 才追加
  `absorber_site="ALL_WEIGHTED"` 汇总行；任一位点失败则该位点写拒绝行，
  且不输出加权汇总。
- `run_xafs.py` 集成：`--absorber-site all` 时先对全部位点过门禁，
  任一位点失败即抛 `XAFSConfigError`；`xafs_fit.run` 后用多位点结果覆盖
  `formal_cn`。
- 真实 Fe3O4 谱（BMM `Fe-Magnetite.xdi`）端到端：3 位点门禁全部通过；
  Fe1/Fe2 的条件 CN 因 N/N₀=0.553 超出范围被拒绝，Fe3 通过；最终
  `status=multi_site_failed`，`cn_table.csv` 写拒绝行且无 ALL_WEIGHTED
  汇总行，符合“失败位点不出汇总”的契约。
- 契约 #127 锁定：未传 site 拒绝、Fe3O4 权重/加权 N₀、Co3O4 1/3:2/3、
  Co3O4 `occupancy_missing` 拒绝、全部位点成功才出汇总行。

## CIF 抓取整合（T10）

构建期工具：`tools/cif_crawler.py`

- `--list`：列出 39 个固定目标（COD ID / 空间群 / 相）
- `--t3` / `--t4` / `--p5c` / `--all`：选择目标
- `--online`：允许联网下载；默认只读 `data/cif_raw/` 缓存
- `--offline`：显式离线；缓存缺失返回 `offline_miss`
- `--dry-run`：只核验，不写 phase pack
- `--force`：忽略缓存重新下载（需 `--online`）
- `--replace`：同 COD ID 不同 sha256 时显式覆盖 phase pack
- `--no-index`：测试时不重建索引/覆盖报告

硬门禁：CIF 必须能被 gemmi 解析、`_chemical_formula_sum` 计量比匹配、IT 空间群号命中目标、目标吸收元素有位点，且不是 `structure_lib._cif_path` 现场生成 CIF。工具不写 motif、不做化学式自由检索、不使用 AFLOW/POSCAR 兜底。

运行期边界：`run_xafs.py` 不调用该工具；现场分析仍完全离线，只读 phase pack / phase index / 用户 `--phase-cif`。


## P6：催化分析可执行性修复（已完成）

- 单晶位多元素晶体：以不对称单元 `st.sites` 判定歧义，39 条单晶位相
  默认可直接过门禁；
- 有序化学计量晶体：`_chemical_formula_sum` 与全占据晶胞计量比一致时
  允许 `implicit_1.0`；
- 键长：新增 `PAIR_BOND_REF` + `PAIR_BOND_SLACK = 0.08 Å`，覆盖高价/
  镧系/锕系 M–O；
- `Rh2O3` 改用 `COD:1544346`（`P b n a` #60，`6 × Rh–O @ 2.050 Å`）；
- `--catalytic-model + --formal-cn` → `hypothesis_model_only`；
- 契约 `#128` 锁定 44/44 phase pack 全量门禁。
