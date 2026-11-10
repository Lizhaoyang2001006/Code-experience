# 催化剂 XAFS 分析：能力、结果、置信门禁与参数台账

> 适用版本：P6.6（图包分层与双入口产物对齐）
> 代码依据：所有阈值均来自当前代码常量，不是文档口径。
> 配套契约：`regression_test.py` 第 `#132` 项。

---

## 1. 这套代码的定位

本工具把催化剂 XAFS 分析分成两条**互不冒充**的轨道：

| 轨道 | 输入 | 产出 | 适用场景 |
|---|---|---|---|
| **诊断轨** | 任意 EXAFS 谱 + motif / 催化模型 / L1 第一壳 | `fit.json`、`fit_table.csv`、`fit.csv`、`rdf.csv` | 机制假设、模型筛选、论文 SI 的拟合表 |
| **正式轨** | `verified_source_cif` + `crystal_cif/crystal_edited` + 显式吸收位 | `cn_table.csv` | 可辩护的配位数、位点占据、相分数 |

两条轨道的结论**分开写、状态分开标**：

- 诊断轨永远可以出数，但 `cn_status` 可能是 `diagnostic_only` / `hypothesis_model_only`；
- 正式轨只在来源、CIF、键长、占据度、条件 CN 全部通过时签发 `gate_status=passed`；
- 正式轨拒签**不会删除**诊断轨结果，`reject_reason` 会写明物理原因。

一句话：**诊断轨回答“像什么”，正式轨回答“配位数能不能写进论文”。**

---

## 2. 执行逻辑（S1–S7 + 两条轨道）

```text
原始谱
  → S1 读取 / 多扫描合并 / 去死时间 / 去毛刺 / 去 glitch / 能量对齐
  → S2 元素判定 + E0 判定（导数 / 半高点 / 手动）
  → S3 归一化（边前线 + 边后窗口 + 平台门禁 + 自吸收判断）
  → S4 Rbkg 寻优（或手动指定）
  → S5 背景扣除 + χ(k) + k 权重 + FT → R 空间
  → S6 连续小波变换（CWT）
  → S6.5 结构识别 / 相来源解析 / 吸收位解析
  → S7 EXAFS 拟合
       ├─ 诊断轨：motif / L1 / SAC-DAC / 多壳层自由拟合
       └─ 正式轨：verified CIF → CIF 门禁 → N₀ 分壳 → 条件 CN
  → 报告层：HTML / CSV / PNG / 参数与门禁台账
```

| 阶段 | 做什么 | 关键判定 | 输出 |
|---|---|---|---|
| **S1** | 读取扫描、合并多扫描、剔除离群、死时间、去毛刺、去 glitch、对齐 | 扫描质量、离群比例 | 合并 μ(E) |
| **S2** | 元素与吸收边；E0 用导数、半高点或手动 | 元素必须在 17 元素表内 | `absorber`、`e0` |
| **S3** | 边前线拟合、归一化、平台检查、自吸收风险 | 平台 0.85–1.15；荧光 Δμ ≥1.0 必须校正 | `mu_norm`、`edge_step`、`plateau` |
| **S4** | Rbkg 自动寻优或手动 | 置信度、有效候选比例、崎岖度 | `rbkg` |
| **S5** | 背景扣除、χ(k)、k 加权、FT | kmin / kmax / kweight | `chi.dat`、`ft.dat` |
| **S6** | Cauchy 连续小波变换 | σ 默认 1.5；第一峰 R | `wt.dat`、`wavelet.csv` |
| **S6.5** | 结构识别：motif 库 / phase pack / 用户 CIF / 显式催化模型 | score、margin、未解释峰、coverage | `structure_id.json` |
| **S7** | 多壳层 EXAFS 拟合 | Rf、ΔE0、N_err、σ²/C3/C4 撞界、N_idp | `fit.json`、`fit.csv`、`fit_table.csv` |
| **正式 CN** | CIF 门禁 → N₀ 分壳 → 固定 N → 条件 N → Hamilton F | N/N₀、条件 Rf、F-test | `cn_table.csv` |
| **报告** | 汇总 HTML、XANES/EXAFS 图、拟合表、小波、LCF、置信椭圆 | — | `_report.html` 等 |

### 2.1 五个入口

| 入口 | 用途 | 关键限制 |
|---|---|---|
| `run_xafs.py` 单样品 | 常规单谱处理 + 可选识别/拟合/正式 CN | — |
| `run_xafs.py --batch-dir` | 批量同元素多样品 | 不逐文件识别元素 |
| `run_xafs.py --multi-element` | 多元素各自处理并汇总 | `--phase-id` 为单一值 |
| `run_xafs.py --operando` | 时间序列逐点 S1–S6 + PCA/MCR-ALS/LCF | 不执行 S6.5/S7 |
| `xafs_report.py` | 单样品完整报告 | 需要 `--absorber` |
| `xafs_figure_pack.py` | 对已有输出补 Figure 1–4 图包与数据包 | 只读、不重新拟合；operando 需对 `scan_XX` 单独运行 |

---

## 3. 催化剂体系能力矩阵

| 体系 | 能实现的部分 | 推荐用法 | 明确边界 |
|---|---|---|---|
| **Fe-N-C（SAC）** | N₄ 平面、轴向 O/N 配位、Fe–N 键长/配位数、Fe–C 远壳、小波判别金属团簇、Hamilton F-test、XANES 边前峰/白线/LCF | `--catalytic-model "sac:N4@2.00,C@2.90"`；轴向：`sac:N4@2.00,axial=O@2.10,C@2.90` | 无周期晶胞，正式 CN 标 `hypothesis_model_only`；诊断拟合表正常出数 |
| **Fe-N-C（DAC）** | Fe–Co / Fe–Ni 双原子模型、异核路径、双金属距离 | `--catalytic-model "dac:Co@2.50,N4@2.00"` | 同上；DAC 几何是假设模型 |
| **Fe-ZSM-5** | 第一壳 Fe–O、双核 Fe···Fe、低负载荧光谱自吸收校正、BVS 价态约束 | `--no-strict` + 自吸收参数 + `--bvs-target` | 第二壳 Fe–O–Si/Al 需用户提供含占据度的 CIF/DFT 团簇，再 `--phase-cif` + `--absorber-site` |
| **金属氧化物（单晶位）** | 44 条实验 COD 相包，17 元素；CIF 门禁、N₀ 分壳、条件 CN、正式 CN 表 | `--formal-cn --phase-id COD:xxxx` | 高缺陷纳米样可能被正式轨拒签；诊断轨仍出数 |
| **金属氧化物（多晶位/尖晶石）** | 不等价位点分别过门禁；`--absorber-site all` 按多重度×占据度加权；`--fit-weights` 放开相分数 | `--phase-id COD:xxxx --absorber-site all` | 多晶位不自动 Top-1；未显式给位点 → `absorber_site_ambiguous` |
| **掺杂氧化物 / 异核共占位** | `--hetero-cocupancy` 拆分共占位为 M–Ni / M–Fe 分支，共享 ΔR/σ² | `--hetero-cocupancy` | 需要 CIF 中真实共占位信息 |
| **双金属合金** | Cowley 短程有序参数 α_AB，判别随机固溶/有序/偏析 | 合金相 + `--fit` | 需要合金相 CIF 或 motif |
| **Operando** | 逐时间点 S1–S6，PCA（99% 方差）、LCF、MCR-ALS、组分浓度演化 | `--operando` | 不执行 S6.5/S7；`--fit/--formal-cn/--phase-*` 会显式报错 |
| **通用 XANES** | 边前峰质心/面积、白线面积、E0、导数谱、LCF 价态比例 | `--lcf-standards` | LCF 需要可靠标准谱 |
| **通用 EXAFS** | FT/反 FT、q 空间拟合、小波、相移校正、χ² 置信椭圆 | `--fit-space q`、`--phase-correct`、`--chi2-contour` | — |
| **标准谱 / T5 S₀²** | 12/12 元素独立标定产物；运行期 `--s02-from-t5` 加载 | `--s02-from-t5` | 12 元素静态表仍为 `None`，不写回 `absorber.py` |
| **批量/多元素** | 批量同元素；多元素分别处理并汇总 | `--batch-dir`、`--multi-element` | `--batch-dir` 不逐文件识别元素；`--multi-element` 的 `--phase-id` 为单一值 |

### 3.1 明确不会自动做的事

- 不会自动选多晶位 Top-1；必须显式 `--absorber-site`。
- 不会为缺占据度的无序 CIF 猜占据度；直接 `occupancy_missing`。
- 不会在 N/N₀ 越界或条件 Rf 超限时“放宽”正式 CN；只拒签，不篡改。
- 不会把 SAC/DAC 假设模型写成 `crystal_cif` 正式来源。
- 不会自动为 17 元素之外的吸收元素生成完整相库。
- 不会在 `--operando` 下静默忽略 S6.5/S7 参数；现在直接报错。
- 不会把 `data/standards_raw/` 的原始标准谱打包进仓库（许可原因）；运行期 T5 标定产物已随仓库。

---

## 4. 结果清单：每个有效输出对催化剂分析的价值

> 这一节是本文档的核心。每个输出都回答三个问题：
> **它是什么 → 能说明催化剂的什么 → 使用时要注意什么。**

### 4.1 归一化谱与 XANES

| 输出 | 是什么 | 对催化剂分析的价值 | 注意事项 |
|---|---|---|---|
| `<stem>_mu.dat`（run_xafs） | E、μ、μ_norm、σ_μ | 看吸收边位置、白线、边前峰；判断活性金属氧化态、电子结构、配位对称性 | 边移也受配位/共价性影响，不能只看 E0 定价态 |
| `<stem>_mu.csv` / `_xanes.csv`（report） | 归一化 μ(E) 与 XANES 区间 | 同一用途；可直接作为 XANES 图、LCF 输入、operando 时间序列 | 归一化平台必须可信；低负载荧光谱注意自吸收 |
| `<stem>_indicators.csv` | 边前峰质心/面积、白线面积、E0、边步等 | **氧化态**：边移 + 白线面积；**位点对称性**：边前峰强度/质心（四面体 vs 八面体、Fe–N₄ vs Fe–N₅）；**负载/浓度**：边步相对比较；**Operando**：这些指标随时间的突变对应活性相变 | 白线受自吸收、归一化、样品厚度影响；跨样品比较必须同一边、同一处理 |
| `<stem>_lcf.json` | XANES 线性组合拟合权重、Rf、权重误差 | **价态分布/相比例**：例如 Fe²⁺/Fe³⁺ 比例、金属/氧化物比例、反应过程中相分数；可直接画堆叠柱状图 | 标准谱必须与样品基体匹配、线性独立；看 Rf 与权重误差，不只看权重 |
| `<stem>_xas.png` | 汇总 XANES/EXAFS 图 | 论文正文图、组会/合作者快速审阅 | 图内标注的窗口/参数要与 `_meta.json` 一致 |

**催化解读示例：**

- Fe-N-C 的边前峰在 ~7112–7114 eV 出现强 1s→3d 特征，通常支持 Fe–N₄ 平面（近似 D4h）；若边前峰明显增强、白线向高能移动，可能说明轴向配位/氧化或结构畸变。
- Fe-ZSM-5 中 Fe²⁺ 与 Fe³⁺ 的边位置差约 2–3 eV；用 LCF 或边移跟踪还原/氧化处理后的价态变化。
- CeO₂ 的 Ce L₃ 边白线双峰与 Ce³⁺/Ce⁴⁺ 比例相关；LCF 比单看边移更稳。

### 4.2 χ(k)、FT、q 空间

| 输出 | 是什么 | 对催化剂分析的价值 | 注意事项 |
|---|---|---|---|
| `<stem>_chi.dat` / `_chik.csv` | k、χ、kⁿχ | **键长**：振荡频率；**配位数/无序**：振幅；**重散射体**：高 k 振荡；反应前后直接对比可看局域结构变化 | 振幅受 S₀²、自吸收、σ² 共同影响，不能单独把振幅当 N |
| `<stem>_ft.dat` / `_chir.csv` | R、Re、Im、\|χ(R)\| | **壳层结构**：第一壳 M–O/N、第二壳 M–M/M–C；峰位给出表观键长；峰高给出相对配位/无序；实部/虚部用于拟合对比 | FT 峰位未做相移校正，不等于真实键长；真实键长看拟合 R 或 `_phase_corrected` |
| `<stem>_chiq.csv` | q 空间反傅里叶 | 检查背景扣除/反变换质量；q 空间拟合可降低低 R 背景泄漏的影响 | 两个入口都会落盘；不是独立的结构信息，主要用于拟合质量与一致性 |
| `<stem>_rdf.csv` | 非高斯 g(R)（放开 C3/C4 时） | **非谐性/静态无序**：纳米颗粒表面弛豫、缺陷、界面应力；判断 σ² 是否被高斯模型低估 | 只有显式放开 C3/C4 才有；不要用它替代常规 CN |

**催化解读示例：**

- Fe-N-C 的第一壳 Fe–N 约 2.0 Å；若第二壳在 2.5–2.8 Å 出现明显 Fe–Fe，需高度怀疑金属团簇/纳米颗粒，而不是单原子。
- Fe-ZSM-5 的双核 [Fe₂(μ-O)] 常见 Fe···Fe ~3.0 Å；这个峰与 Fe–O–Si/Al 第二壳重叠时，必须用 WT + 多壳层拟合 + 约束区分。
- 金属氧化物中 M–O 第一壳缩短通常对应氧化态升高；M–M 第二壳缩短可能对应还原/氧空位。

### 4.3 小波与相移校正

| 输出 | 是什么 | 对催化剂分析的价值 | 注意事项 |
|---|---|---|---|
| `<stem>_wavelet.csv` / `_report` 内小波图 | k–R 二维强度 | **区分轻/重散射体**：N/O/C 在第一壳低 k 区，Fe/Co/Ni 在高 k 区；**单原子 vs 团簇**：金属–金属键在 2–3 Å、高 k 处出现强信号；**轴向配位**：不同 k 特征的壳层分离 | 小波峰位受 σ、k 权重影响；用 `--wt-sweep` 检查稳健性 |
| `<stem>_wt_sweep/*` | 多 kweight / σ 的 WT 对比 | 选择最优 k 权重、确认重金属壳不是伪峰；审稿人问“为什么用这个 k 权重”时直接给图 | 色标按行归一，不能跨行直接比绝对强度 |
| `<stem>_phase_corrected.*` | 相移校正后的 FT / 真实键长 | **精确 M–O/M–M 距离**：把表观峰位修到真实键长，便于和 CIF/DFT 对比；判断键长变化是否超出相移误差 | 需要第一壳路径；相移来自 FEFF，模型不对时校正也会偏 |

### 4.4 结构识别与拟合

| 输出 | 是什么 | 对催化剂分析的价值 | 注意事项 |
|---|---|---|---|
| `<stem>_structure_id.json` | 候选排名（全量 `ranked` + `top3`）、score、margin、未解释峰、coverage | **物相判断**：最像哪个相、是否有并列/混合相/库缺条目；`partial` 提示模型不完整；`no_match` 提示库没有正确相，应改用用户 CIF 或 L1 | 这是诊断轨，不是正式晶体学证明；`confident` 也不等于结构确证 |
| `<stem>_fit.json` | 拟合参数、误差、相关性、flags、`cn_status` | **配位数 N**：活性位饱和度/欠配位；**键长 R**：氧化态/吸附物；**σ²**：静态无序/缺陷；**ΔE0**：电子结构/模型偏差；**Rf**：模型质量；**flags**：哪些壳不可信 | N 与 σ² 强相关；看 `N_err`、`sigma2_at_bound`、`r_factor_bad` |
| `<stem>_fit.csv` / `_fit_table.csv` | 壳层、N、N_err、R、R_err、σ²、Rf、N_idp、N_var | **论文 SI 表 S1**：直接给出每个壳的拟合参数与误差；跨样品比较活性位结构 | 表中 N 是简并度/有效简并度；`motif_weights` 模式下 N 不是 CN |
| `<stem>_fit_table.png` | 拟合表图 | 组会/审稿人快速审阅；与 `_xas.png` 配套 | 只由报告入口生成；原 `_fit_table.txt` 与 CSV 重复，P6.6 起不再输出 |
| `<stem>_chi2_contour.*` | χ²(N, σ²) 置信椭圆 | **回应审稿人**：“N 和 σ² 是否强相关？” 直接展示 1σ/2σ/3σ 置信区；判断 CN 的不确定度是否被低估 | 只在有 S7 拟合时生成；椭圆形状依赖拟合窗口与模型 |

**催化解读示例：**

- Fe-N-C：N≈4、R(Fe–N)≈2.0 Å、σ²≈0.003–0.006 Å²，且无 Fe–Fe 壳 → 支持 Fe–N₄ 单原子位点；若 N 明显低于 4，可能是表面欠配位或模型缺失轴向配体。
- Fe-ZSM-5：N(Fe–O)≈4–6、R≈1.9–2.1 Å；若 N/N₀ 在 0.7–1.3 且无 Fe–Fe，支持孤立 Fe 位点；若出现 Fe–Fe，考虑双核/寡聚 FeOₓ。
- Co₃O₄：`--absorber-site all` 可分别得到四面体 Co²⁺ 与八面体 Co³⁺ 的 CN/键长；加权汇总对应整体化学计量。

### 4.5 正式 CN 表（论文级）

| 输出 | 是什么 | 对催化剂分析的价值 | 注意事项 |
|---|---|---|---|
| `cn_table.csv` | 39 列正式 CN：相来源、CIF sha256、吸收位、壳层、N₀、N_fit、N/N₀、形式/系统误差、ΔR、σ²、ΔE0、Rf、S₀² 来源、窗口、N_idp/N_var、gate_status、reject_reason、flags、assumptions、数据 sha256、软件版本 | **可辩护的配位数**：每个壳的 N、N₀、N/N₀ 与系统误差；**位点占据/相分数**：多位点加权行；**可复核性**：CIF 哈希、数据哈希、软件版本；**审稿人问“为什么 N 是这个值”**：直接给 N/N₀ 与条件 CN 门禁 | 只有 `gate_status=passed` 才能当正式 CN；`warning` 要说明 N/N₀ 偏差；`consistency_failed` 只能作为拒签记录 |
| `_meta.json` 的 `formal_cn` 块 | 正式 CN 完整状态、拒绝原因、模型来源 | 快速判断“为什么没签正式 CN”；区分 `no_cif` / `consistency_failed` / `hypothesis_model_only` | 与 `cn_table.csv` 的拒绝行一致 |

**`cn_table.csv` 的关键列怎么读：**

| 列 | 催化含义 |
|---|---|
| `phase_id` / `model_level` / `source_type` | 用的是哪个相、什么级别的模型（实验 CIF / 编辑 CIF） |
| `source_cif_sha256` | 相来源可复核；审稿人可追溯到具体 CIF |
| `absorber_site` | 哪个晶位（如 Fe1、Fe2、frac:(0,0,0)） |
| `shell_group` / `shell_label` | 第一壳 / 远壳，具体配位壳 |
| `n0` / `n0_members` | 晶体学 N₀ 与参与合并的路径数（Rayleigh 分壳） |
| `n_fit` | 条件 CN 拟合得到的 N |
| `n_status` | `fitted`（第一壳自由 N）/ `fixed`（远壳固定 N₀） |
| `n_ratio` | N/N₀；0.7–1.3 正常，0.6–0.7 / 1.3–1.4 警告，越界拒签 |
| `n_err_formal` | 形式误差（拟合协方差） |
| `n_err_systematic` | 系统误差预算（自吸收、参数相关、截断、σ² 撞界等） |
| `delta_r` | 相对晶体学键长的膨胀/收缩 |
| `sigma2` | 无序；与 N 强相关 |
| `delta_e0` | 条件拟合的 E0 偏移 |
| `rf` | 条件拟合 Rf |
| `s02` / `s02_source` / `s02_calibrated` | S₀² 值、来源（静态表 / T5 / 运行时覆盖）、是否标定 |
| `kmin/kmax/kweight/rmin/rmax` | 拟合窗口；跨样品比较必须一致 |
| `n_idp` / `n_variables` | 独立点数/变量数；判断过拟合 |
| `gate_status` | `passed` / `warning` / `consistency_failed` |
| `reject_reason` | 拒签的物理原因 |
| `flags` | 可靠性标记（如 `n_ratio_warning`） |
| `assumptions` | 条件 CN 的假设（首壳单一比例因子、远壳固定 N₀、Hamilton 不显著等） |
| `data_sha256` | 输入数据哈希 |
| `software_version` | 软件版本（`xafs-cn-p4plus`） |

### 4.6 报告、批量与多元素

| 输出 | 是什么 | 对催化剂分析的价值 | 注意事项 |
|---|---|---|---|
| `<stem>_report.html` | 单文件交互报告（图内联、门禁结论置顶） | 合作者/导师一键审阅；论文图与 SI 的草稿；把“数据质量、结构识别、拟合、门禁”放在同一页 | 不是正式发表文件；最终图仍需按期刊要求重绘 |
| `<stem>_window_tuner.html` | 交互式 k/R 窗口调参 | 检查窗口选择是否稳健；避免“挑窗口挑出想要的 N” | 调参后要重新跑并记录最终窗口 |
| `multi_element_summary.md/.json` | 多元素并列结果（不推荐 Top-1） | **主元素 vs 掺杂元素**：各自局域环境、配位、价态；判断掺杂是否进入晶格、是否形成分离相 | 不自动跨元素排序；科学问题由用户决定 |
| `batch_errors.json` | 批处理失败文件与原因 | 快速定位坏扫描/格式问题；保证批量分析可复现 | 单个失败不中断整批 |
| `_meta.json` 的 `params/gates/decisions/summary` | 本次运行的全部生效参数、来源、门禁命中、决策轨迹、总结 | **方法部分**：直接抄参数与窗口；**审稿人复核**：每个阈值/状态/拒绝原因都有记录；**失败诊断**：一眼看出是数据、CIF、拟合还是正式 CN 的问题 | 只增不改；现有字段保持原样 |

### 4.7 Operando 输出

| 输出 | 是什么 | 对催化剂分析的价值 | 注意事项 |
|---|---|---|---|
| `_operando_components.csv` | MCR-ALS 纯组分谱 | **反应中间体/活性相**：盲解谱得到随时间变化的纯组分；与标准谱/DFT 对比指认 | 组分顺序按特征能量排序；MCR 有旋转不确定性，需化学约束/标准谱验证 |
| `_operando_concentrations.csv` | 组分浓度 vs 时间 | **动力学/诱导期/失活**：组分分数随时间的变化；与催化活性数据做结构–活性关联 | 浓度是相对分数；closure 假设要记录 |
| `_operando_summary.png` | 四宫格：MCR 分数 / LCF 权重 / E0–白线描述符 / 纯组分谱 | 一张图讲清反应过程中结构与价态演变 | 图内标注的组分顺序要与 CSV 一致 |
| `_exafs_evolution.csv/.png` | FT 随时间的演化 | **键长/无序变化**：反应条件下活性位是否重构、是否有金属–金属键生成 | FT 峰位是表观值；关键时间点应做正式拟合 |
| `operando_meta.json` | 原位序列参数 + PCA/SPOIL/MCR 门禁 | 记录 `n_components`、接受/拒绝的标准谱、MCR Rf；保证原位分析可复现 | 与逐扫描的 `_meta.json` 分开 |

### 4.9 催化剂结构图包（P6.5 / P6.6）

`xafs_figure_pack.py` 只读已有输出，整理为可直接用于论文/组会/审稿回复
的图包；不重新拟合、不改变 `cn_table.csv`、`fit.json`、`_meta.json` 或
任何拒签逻辑。P6.6 起默认**不再重复渲染 28 张单面板图**，并取消
`figures/data/` 中与父目录/逐面板 CSV 重复的副本。

| 输出 | 是什么 | 对催化剂分析的价值 | 注意事项 |
|---|---|---|---|
| `figures/Figure1_master.png` | 8 面板结构总图：(a) XANES (b) χ(k) (c) FT (d) 小波 (e) 拟合壳层 (f) 结构识别 (g) 正式 CN 门禁 (h) 结构卡片 | 一张图解释“这个催化剂是什么结构、拟合到了什么、CN 是否过门禁” | 每个面板在 `panels/` 有 CSV；缺失数据写 `skipped`，不伪造 |
| `figures/Figure2_confidence.png` | 6 面板置信/过程仪表盘：门禁状态板、决策轨迹、参数来源、门禁表、决策表、关键参数表 | 回答“结果为什么可信/为什么被拒”，把 P6.4 台账直接变成图 | 阈值与 `_meta.json:gates` 一致；关键参数表限 14 行，完整值查 `all_parameters.csv` |
| `figures/Figure3_diagnostics.png` | 8 面板拟合与误差诊断：参数相关矩阵、路径 N–R（标注 ΔR/σ²）、CN 误差预算、Rbkg 寻优、归一化分解、拟合窗口、稳定性检查、扫描/数据质量 | 证明 N 的不确定度不是拍脑袋；定位 N–σ² 补偿、撞界、自吸收和过拟合 | 相关矩阵来自 `fit.json`；Rbkg 曲线需 `_meta.json:rbkg.table` |
| `figures/Figure4_features.png` | 6 面板数据特征与输入清单：dμ/dE、q 空间、FT 峰识别、CIF 配位壳层、特征总表、输入文件清单（行数/大小/SHA256） | 补齐“峰在哪、q 空间是否正常、CIF 有哪些壳层、结果基于哪些文件” | q 空间面板需两侧入口输出的 `_chiq.csv`；未识别峰提示模型缺壳层或多位点 |
| `figures/Figure2g_cn_table.png` / `Figure2h_fit_table.png` | Table S1/S2 表格直接出图 | 论文补充材料；无需再用 Excel 排版 | CSV 保留完整列 |
| `figures/panels/<panel>.csv` | 每个面板唯一的实际绘图数据 | 图与数据一一对应，可直接重绘/复查 | 面板 `status=skipped` 时无 CSV |
| `figures/panels/<panel>.png` | 可选的单面板独立图 | 单独插入 PPT/补充材料 | 仅 `--standalone-panels` 时生成；默认不重复渲染 |
| `figures/data/` | 跨面板汇总：`all_fit_paths.csv`、`all_parameters.csv`、`fit_variants.csv`、`gate_checks.csv`、`decisions.csv`、`panel_index.csv`、`summary.json`、`README.md` | 审稿回复可直接给长表；补图不必回挖输出目录 | 不再复制父目录源 JSON/CSV，也不重复逐面板表；权威原件仍在输出根目录 |
| `figures/figure_manifest.json` | Figure 1(a) ↔ 总图/独立图 ↔ CSV ↔ 源文件 ↔ 图注映射 | 机器可读索引，便于批量排版与追溯 | 默认 `image` 指向所属总图；`--standalone-panels` 后指向 `panels/<id>.png` |
| `figures/figure_pack.md` | 人可读索引 | 快速定位每个面板解释什么 | 与 manifest 同步生成 |

两个入口均支持 `--figure-pack` / `--figure-formats` / `--figure-dpi`；
独立 CLI 另有 `--standalone-panels`（默认关闭）。`--operando` 明确不支持
单结构总图，需对 `scan_XX` 目录单独运行 `xafs_figure_pack.py`。

### 4.10 诊断图

| 输出 | 是什么 | 对催化剂分析的价值 |
|---|---|---|
| `diagnostics/*` | 各阶段诊断图、反 FT 对比、相移校正、χ² 等高线、WT 切片 | 区分“真实结构变化”与“数据/背景/自吸收伪迹”；审稿人补充材料 |
| `diagnostics/*chi2_contour*` | N–σ² 置信椭圆 | 证明 CN 的不确定度不是拍脑袋 |
| `diagnostics/*phase_corrected*` | 相移校正真键长 | 与 CIF/DFT 键长直接对比 |

---

## 5. 关键数值指标怎么解读

| 指标 | 催化含义 | 可靠区间/门禁 | 常见误读 |
|---|---|---|---|
| **E0 边移** | 氧化态、电子结构 | 与同一边标准谱比较；ΔE0 拟合门禁 ±15 eV | 配位/共价性也会移动 E0；不能只看 E0 定价态 |
| **白线面积/高度** | d 带空位、氧化态 | 同一边、同一归一化比较 | 自吸收会压低白线；厚样/荧光样必须先评估 |
| **边前峰质心/面积** | 配位对称性、自旋态、共价性 | 与标准谱对比 | 双峰/多峰拟合方式会影响质心 |
| **N（配位数）** | 活性位饱和度、欠配位、吸附 | 看 `N_err`；`N_err >1` 或 `>0.5×简并度` 不可靠 | N 与 σ² 强相关；不能只看 N |
| **R（键长）** | 氧化态、吸附物、应力 | 与 CIF/DFT 比较；看 `R_err` | FT 峰位不是真实键长；用拟合 R 或相移校正 |
| **σ²** | 静态/热无序、缺陷、界面 | 撞界 → 不可靠；系统误差 +10% | σ² 很小不一定是“好”，可能是模型形状失配 |
| **ΔE0** | 电子结构/模型偏差 | `|ΔE0| ≥ 15 eV` 不可靠；`|ΔE0| > 5 eV` 记非阻断 advisory | 大 ΔE0 往往说明模型缺壳层/相不对；也可能是 FEFF6L 路径零点与实验 E0 的系统差，不能直接当边移 |
| **Rf** | 拟合质量 | `>0.08` 不可靠；`>0.12` 很差 | Rf 好不等于 N 可靠（N–σ² 补偿） |
| **N_idp / N_var** | 过拟合 | `N_var ≥ N_idp` 过拟合 | 加路径/加变量会虚降 Rf |
| **N/N₀** | 正式 CN 一致性 | 0.7–1.3 ok；0.6–0.7 / 1.3–1.4 warning；越界拒签 | N/N₀ 低不一定是缺陷，也可能是纳米/表面/模型不对 |
| **S₀²** | 振幅因子 | T5 通过 0.4–1.15；运行时覆盖 0.4–1.3 | 未标定元素不能拿 S₀² 硬凑 N |
| **小波 k–R** | 轻/重散射体、单原子/团簇 | 多 k 权重扫描稳健 | 小波峰不是定量配位数 |
| **LCF 权重** | 价态/相比例 | 看 Rf 与权重误差 | 标准谱不独立/基体不匹配会伪拟合 |
| **SPOIL** | 标准谱质量 | <3 good；3–6 fair；6–12 poor；>12 拒绝 | 只用于 LCF 标准筛选，不是样品质量 |
| **MCR Rf** | 盲解谱重构质量 | 越低越好；与组分合理性一起看 | 低 Rf 不代表组分化学意义唯一 |

---

## 6. 置信门禁卡

### A. 数据质量与归一化

| 门禁 | 阈值 | 触发动作 | 代码位置 |
|---|---|---|---|
| 元素白名单 | 17 元素 | 未知元素拒绝 | `absorber.py` |
| 边后归一化平台 | 0.85–1.15 | 超出 → strict 报错；`--no-strict` 降为警告 | `xafs_pipeline.py:196` |
| 边后窗口 | `0 < lo < hi`，默认 150–600 eV | 非法 → 参数错误 | `run_xafs.py` / `xafs_report.py` |
| 自吸收（边跃变） | ≥1.0 必须校正；0.3–1.0 警告 | 未校正且可能 → `cn_gate` → `cn_reliable=False` | `xafs_fit.py:56-57, 103-109` |
| XANES 交叉校验 | 白线/边步/归一化异常 | 仅警告；不改变 `cn_reliable` | T8 模块 |
| 结构覆盖 | `coverage_ok=False` | 拟合不可靠 | `xafs_fit.py:1511-1516` |

### B. 正式来源门禁（P0b）

| 门禁 | 要求 | 失败码 | 代码位置 |
|---|---|---|---|
| 来源状态 | `provenance_status == verified_source_cif` | `unverified_source` | `structure_model.py:202-208` |
| 模型级别 | `crystal_cif` 或 `crystal_edited` | `unverified_source` | `structure_model.py:209-215` |
| CIF 存在 | 路径存在 | `no_cif` / `source_cif_missing` | `structure_model.py:216-222` |
| 非生成 CIF | 不在 `_structure_lib/*/_cif/`，无生成标记 | `generated_cif_not_formal` | `structure_model.py:223-229` |
| SHA256 | 当前文件哈希 = 记录哈希 | `source_hash_missing` / `source_hash_mismatch` | `structure_model.py:230-244` |
| 吸收位 | 必须显式 `absorber_site` | `absorber_site_missing` | `structure_model.py:245-247` |

### C. CIF 晶体化学门禁（P2a）

| 门禁 | 阈值 | 触发动作 | 代码位置 |
|---|---|---|---|
| 第一壳 CN | 2–12 正常；12–16 警告；<2 或 >16 拒绝 | 拒绝，除非 `--allow-unusual-coordination` | `cif_gate.py:64-65, 616-627` |
| 键长 | `PAIR_BOND_REF` 区间 ±0.08 Å | `bond_length_implausible` | `cif_gate.py:90, 637-649` |
| 共占位合并 | 0.08 Å | 无法合并 → `cocupancy_unresolved` | `cif_gate.py:67, 517-575` |
| 占据度 | 缺 `_atom_site_occupancy` 且非有序化学计量 | `occupancy_missing` | `cif_gate.py:566-568` |
| 多面体闭合 | 配位多面体必须闭合 | `polyhedron_not_closed` | `cif_gate.py:660` |
| 吸收位解析 | 多晶位未指定 | `absorber_site_ambiguous` | `cif_gate.py:494-536` |

### D. 结构识别门禁（S6.5）

| 门禁 | 阈值 | 触发动作 | 代码位置 |
|---|---|---|---|
| 最优分数 | `score < 0.16` 才可能 confident | 超出 → `no_match` | `structure_id.py:66, 341-342` |
| 并列容差 | `Δscore < 0.05` | 并列报告 | `structure_id.py:54` |
| 最优-次优 margin | `< 0.03` | `ambiguous` | `structure_id.py:67, 344` |
| 未解释显著峰 | 存在 | `partial` | `structure_id.py:339-346` |
| 显著峰最小半高宽 | 0.08 Å | 窄于此视为噪声尖刺 | `structure_id.py:57` |
| 最小 prominence | 0.02 | 去除骑在峰肩的假峰 | `structure_id.py:59` |
| 指纹窗口 | R = 1.0–4.0 Å | — | `structure_id.py:63` |
| 指纹 σ² 网格 | 0.003 / 0.006 / 0.009 Å² | — | `structure_id.py:65` |

### E. 拟合可靠性门禁（S7）

| 门禁 | 阈值 | 触发动作 | 代码位置 |
|---|---|---|---|
| Rf | `> 0.08` | `cn_reliable=False`，flags `r_factor_bad` | `xafs_fit.py:1512, 1540` |
| ΔE0 | `|ΔE0| ≥ 15 eV` | 不可靠 | `xafs_fit.py:1513, 1541` |
| N 形式误差 | `N_err > 1.0` 或 `> 0.5×简并度` | 不可靠 | `xafs_fit.py:1500-1503` |
| σ² 撞界 | 任一壳撞界 | 不可靠；系统误差 +10% | `xafs_fit.py:1507, 1543` |
| C3/C4 撞界 | 任一壳撞界 | 不可靠 | `xafs_fit.py:1508-1509` |
| 过拟合 | `n_variables ≥ N_idp` | 不可靠 | `xafs_fit.py:1512` |
| 参数相关 | `|ρ| ≥ 0.9` 加 10%；0.8–0.9 加 5% | 系统误差预算增大 | `xafs_fit.py:164-172` |
| 远壳截断 | 因 N_idp 截断 | 系统误差 +5% | `xafs_fit.py:174-176` |
| 自吸收未校正 | `risk=possible` | 系统误差 +15% | `xafs_fit.py:159-162` |
| Booth 校正残留 | `corrected & reliable` | 系统误差 +10% | `xafs_fit.py:156-158` |

> `|ΔE0| > 5 eV` 只在 `cn_disclosure.py` / CN 表中给非阻断 advisory，
> 不改变 S7 的 15 eV 硬门禁，也不触发 `--strict-disclosure`。

### F. 条件 CN / 正式 CN 门禁（P3）

| 门禁 | 阈值 | 触发动作 | 代码位置 |
|---|---|---|---|
| N/N₀ | 0.7–1.3 → `ok`；0.6–0.7 / 1.3–1.4 → `warning`；<0.6 或 >1.4 → `out_of_range` | 越界拒签 | `xafs_fit.py:2204-2205, 2347-2349` |
| 条件拟合 Rf | `> 0.08` | 拒签 | `xafs_fit.py:2706-2707` |
| Hamilton F-test | 显著改善自由 N | 拒签（固定 N 不一致） | `xafs_fit.py:2698-2703` |
| 条件拟合过拟合 | `n_variables ≥ N_idp` | 拒签 | `xafs_fit.py:2704-2705` |
| 一致性 | `consistency_passed=False` | `consistency_failed`，`rows=[]` | `xafs_fit.py:2708-2714` |

### G. S₀² / T5 门禁

| 门禁 | 阈值 | 触发动作 | 代码位置 |
|---|---|---|---|
| T5 通过区间 | 0.4–1.15 | 超出不通过 | `calibrate_s02.py:53` |
| S₀² 误差 | ≤0.06 | 超出不通过 | `calibrate_s02.py:54` |
| 标定 Rf | ≤0.12 | 超出不通过 | `calibrate_s02.py:55` |
| ΔE0 | ≤10 eV | 超出不通过 | `calibrate_s02.py:56` |
| 跨窗口极差 | ≤0.05 | 超出不稳定 | `calibrate_s02.py:59` |
| ΔR 限 | −0.03–0.03 Å | 超出拒绝 | `calibrate_s02.py:60` |
| 静态表锁定 | 12 元素保持 `s02=None` | `--apply` 禁用 | `calibrate_s02.py:71-73` |

### H. 降级与拒绝行为

- `strict=True`：S3 平台越界直接中断；报错提示 `--no-strict`。
- `strict=False`：降级为警告，继续处理。
- 正式轨拒签：`cn_table.csv` 写 `consistency_failed` + `reject_reason`，不删除诊断拟合结果。
- 多晶位未指定：`absorber_site_ambiguous`，不自动 Top-1。
- SAC/DAC：`hypothesis_model_only`，不冒充晶体正式 CN。

---

## 7. 参数如何完全显示表露

### 7.1 三个入口、两套参数、一个台账

| 项 | 数量/内容 |
|---|---|
| `run_xafs.py` CLI 参数 | 68 个 |
| `xafs_report.py` CLI 参数 | 65 个 |
| `xafs_config.DEFAULTS` 配置键 | 45 个 |
| `_meta.json` 参数台账 | `params.effective` 记录全部生效值；`params.sources` 标记 `cli/config/default/auto/calibrated/runtime_override` |

### 7.2 来源优先级

```text
CLI 显式值 > 配置文件 > 默认值
自动推导值（E0 / Rbkg / S₀² / 相 ID / 吸收位）标记为 auto 或 calibrated
```

### 7.3 `_meta.json` 新增的四个块

```json
{
  "params": {
    "effective": { "kmin": 3.0, "kmax": 12.0, "...": "..." },
    "sources": { "kmin": "cli", "kmax": "config", "kweight": "default" },
    "units": { "kmin": "Å⁻¹", "rbkg": "Å" },
    "stages": { "kmin": "S5", "rbkg": "S4" },
    "affects_formal_cn": ["absorber", "phase_id", "..."],
    "cli_argv": ["run_xafs.py", "sample.xdi", "..."],
    "inputs": ["sample.xdi"],
    "config_file": "cfg.json",
    "config_keys_applied": ["phase_id", "s02_from_t5"]
  },
  "gates": {
    "normalization": {
      "status": "passed",
      "threshold": [0.85, 1.15],
      "observed": 0.978,
      "action": "none"
    },
    "self_absorption": { "...": "..." },
    "structure_id": { "...": "..." },
    "cif": { "...": "..." },
    "fit": { "...": "..." },
    "formal_cn": { "...": "..." },
    "xanes_cross_check": { "...": "..." },
    "summary": {
      "passed": ["normalization", "fit"],
      "warned": [],
      "failed": [],
      "not_applicable": ["cif", "formal_cn"]
    }
  },
  "decisions": [
    { "step": "S2", "decision": "absorber_and_e0", "detail": { "...": "..." }, "status": "info" },
    { "step": "S6.5", "decision": "phase_resolution", "detail": { "...": "..." }, "status": "info" },
    { "step": "P3", "decision": "formal_cn", "detail": { "...": "..." }, "status": "warning" }
  ],
  "summary": {
    "passed": ["normalization"],
    "warned": ["structure_id"],
    "failed": ["formal_cn"],
    "not_applicable": ["cif"],
    "cn_reliable": false,
    "cn_status": "consistency_failed",
    "formal_cn_status": "consistency_failed",
    "formal_cn_reject_reason": "N/N₀=0.521 (out_of_range) ..."
  }
}
```

### 7.4 终端打印

两个入口都支持：

```bash
python run_xafs.py sample.xdi --absorber Fe --print-gates --print-params
python xafs_report.py sample.xdi --absorber Fe --print-gates --print-params
```

- `--print-gates`：打印门禁表（门禁 / 状态 / 动作 / 观测 / 阈值 / 拒绝原因）。
- `--print-params`：打印全部生效参数（参数 / 生效值 / 来源 / 单位 / 阶段）。
- 默认关闭；即使不打印，`_meta.json` 里也有完整台账。

### 7.5 参数分组（按阶段）

| 阶段 | 参数 |
|---|---|
| **S1** | `align`、`batch_dir`、`deadtime`、`deglitch`、`despike`、`reject_outlier_scans`、`multi_element`、`multi_element_s02`、`operando` |
| **S2** | `absorber`、`edge_energy`、`e0_method`、`manual_e0`、`edge`（report） |
| **S3** | `post_window`、`strict`、`self_absorption_formula`、`self_absorption_thickness`、`self_absorption_density`、`self_absorption_angles`、`self_absorption_line` |
| **S4** | `rbkg` |
| **S5** | `kmin`、`kmax`、`kweight`、`mback`、`chi_std`、`clamp_lo`、`clamp_hi` |
| **S6** | `wt_sigma`（report） |
| **S6.5** | `identify`、`allow_ambiguous`、`fit_motif`、`l1_ligand`、`l1_geometry`、`catalytic_model`、`phase_cif`、`phase_id`、`phase`、`phase_index`、`phase_pack`、`absorber_site`、`offline`、`s02_from_t5` |
| **S7** | `fit`、`fit_rmin`、`fit_rmax`、`fit_weights`、`fit_kweights`、`fit_space`、`include_ms`、`vary_c3`、`vary_c4`、`isotropic_deltar`、`shared_sigma2`、`constraints`、`bvs_target`、`bvs_sigma`、`bvs_b`、`hetero_cocupancy`、`s02` |
| **T8 / P2a / P3** | `xanes_cross_check`、`allow_unusual_coordination`、`formal_cn`、`legacy_auto_fit` |
| **输出** | `diagnostics`、`phase_correct`、`chi2_contour`、`lcf_standards`、`demo`、`no_fit`、`no_html`、`no_window_tuner`、`sample`、`shells`、`stem`、`wt_sweep`、`print_gates`、`print_params` |

---

## 8. 典型工作流

### 8.1 Fe-N-C 单原子

```bash
python xafs_report.py Fe-N-C.xdi --absorber Fe --no-strict \
  --catalytic-model "sac:N4@2.00,C@2.90" --fit \
  --self-absorption-formula booth --self-absorption-thickness 100 \
  --self-absorption-density 0.8 --self-absorption-angles 45,45 \
  --chi2-contour --wt-sweep --print-gates -o out/FeNC
```

看什么：

1. `_wavelet.csv` / WT 图：有没有 Fe–Fe 高 k 信号；
2. `_fit_table.csv`：N(Fe–N)、R、σ²、Rf；
3. `_indicators.csv`：边前峰/白线；
4. `_meta.json`：`gates.fit`、`summary`、`params`。

### 8.2 Fe-ZSM-5

```bash
python xafs_report.py Fe-ZSM5.xdi --absorber Fe --no-strict \
  --self-absorption-formula booth --self-absorption-thickness 80 \
  --self-absorption-density 0.7 \
  --bvs-target 3.0 --fit --print-gates -o out/FeZSM5
```

第二壳 Fe–O–Si/Al 需要用户 CIF：

```bash
python xafs_report.py Fe-ZSM5.xdi --absorber Fe --no-strict \
  --phase-cif Fe_ZSM5_DFT.cif --absorber-site Fe1 \
  --formal-cn --fit --print-gates -o out/FeZSM5_formal
```

### 8.3 单晶位金属氧化物（正式 CN）

```bash
python xafs_report.py RuO2.xdi --absorber Ru --no-strict \
  --formal-cn --phase-id COD:9007541 --s02-from-t5 \
  --fit --print-gates -o out/RuO2
```

看 `cn_table.csv` 的 `gate_status`：
- `passed` → 可直接写论文；
- `warning` → N/N₀ 在 0.6–0.7 / 1.3–1.4，讨论里说明；
- `consistency_failed` → 只作为拒签记录，诊断轨仍可用。

### 8.4 多晶位氧化物（尖晶石）

```bash
python xafs_report.py Co3O4.xdi --absorber Co --no-strict \
  --formal-cn --phase-id COD:xxxx --absorber-site all \
  --s02-from-t5 --fit --print-gates -o out/Co3O4
```

`cn_table.csv` 会给出每个不等价位的行，并追加 `ALL_WEIGHTED` 化学计量加权汇总行（全部位点通过才追加）。

### 8.5 掺杂/双金属

```bash
python run_xafs.py Co-K.xdi Fe-K.xdi --multi-element \
  --absorber Co --fit --print-gates -o out/doped
```

看 `multi_element_summary.md`：主元素与掺杂元素各自的 N/R/σ²，不自动排序。

### 8.6 Operando

```bash
python xafs_report.py scan_*.xdi --operando --absorber Fe \
  --lcf-standards FeO.xdi,Fe2O3.xdi,Fe3O4.xdi \
  --print-gates -o out/operando
```

看 `_operando_concentrations.csv` 与 `_exafs_evolution.png`；关键时间点再单独跑正式拟合。

---

## 9. 结果解读示例

### 9.1 正式 CN 通过

```text
gate_status=passed
n_ratio=0.96
conditional_rf=0.031
hamilton_significant=false
```

解释：N 与晶体学 N₀ 一致，条件拟合质量好，自由 N 没有显著改善；可以写“EXAFS 给出 N=…±…，与体相 CIF 的 N₀ 一致”。

### 9.2 正式 CN 拒签

```text
gate_status=consistency_failed
reject_reason=N/N₀=0.521 (out_of_range) 超出范围，且未提供定量缺陷模型；条件 CN Rf=0.1314 > 0.08
```

解释：样品第一壳配位明显低于体相 CIF，且条件拟合质量差；不能把体相 N₀ 硬写成样品 CN。诊断轨的 L1/自由拟合仍给出 N；论文中应报告诊断 CN 并说明“体相 CIF 不适用于该缺陷/纳米样品”。

### 9.3 SAC/DAC 假设模型

```text
cn_status=hypothesis_model_only
formal_cn.status=hypothesis_model_only
```

解释：这是无实验 CIF 的催化位点假设模型；`fit_table.csv` 的 N/R/σ² 是“模型条件下的诊断值”，不是晶体学正式 CN。论文中应写成“EXAFS 拟合与 Fe–N₄ 模型一致”，而不是“正式 CN=4”。

### 9.4 结构识别 no_match

```text
structure_id.status=no_match
gates.structure_id.action=library_miss_l1_fallback
```

解释：诊断 motif 库没有匹配相；S7 退回 L1 第一壳。若要做正式 CN，请提供用户 CIF 或 `--phase-id`；不要因为 `no_match` 就认定数据无效。

---

## 10. 明确边界

1. 17 元素白名单；`Ti/V/W/Nb/Sn/In` 等需按 `absorber.py` 三步流程接入。
2. 12/17 元素静态 `s02=None`；须 `--s02-from-t5` 或显式 `--s02`。
3. 5 条多晶位相必须显式 `--absorber-site`，不自动 Top-1。
4. 用户 CIF 无 `_atom_site_occupancy` 且非有序化学计量 → `occupancy_missing`。
5. 条件 CN 的 `N/N₀` 越界、`Rf > 0.08`、Hamilton 显著 → 硬拒绝。
6. `data/standards_raw/` 因许可不入库；运行期 T5 标定产物已随仓库。
7. B19 原位 EXAFS 联合精修仍为 Python API。
8. `--batch-dir` 不做逐文件元素识别；混合元素用 `--multi-element`。
9. `--multi-element` 的 `--phase-id` 是单一值，逐元素相需分次运行。
10. `--operando` 不执行 S6.5/S7；相关参数会显式报错。

---

## 11. 版本与契约

- 参数/门禁/决策/摘要台账：`xafs_provenance.py`
- 催化剂结构图包：`xafs_figure_pack.py`
- 契约：`regression_test.py` 第 `#132` 项（台账）、`#133` 项（图包）、
  `#134` 项（双入口产物对齐/图包分层瘦身）与化学态 S₀² 多态契约
- 当前契约基线：135 项（化学态多态新增后）
- 软件版本：`xafs-cn-p4plus`
- 门禁常量冻结：`absorber.py` 静态 S₀² 表、`calibrate_s02.py` 门禁、44 条 phase pack 均未改动（化学态映射为显式 opt-in，不写静态表）
