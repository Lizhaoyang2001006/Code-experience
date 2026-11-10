# XAFS 数据处理 / 结构识别 / CN 拟合管线

> **项目最高级执行规范（Master Spec）**：当前对标 Athena / Artemis / Larch-WT 全功能补齐工程（Phase I：B1–B4 已验收 36 项回归基线；Phase II：B5–B8 进阶补齐至 45 项；Phase III：B9–B12 高阶能力已验收 60 项回归基线；Phase IV：B13–B16 面向催化化学全面分析已验收 68 项回归基线；**Phase V：B17–B20 催化剂精细结构深度解析与实战完备化已验收 74 项回归基线**）的完整数学契约、函数签名、小模型防呆清单与旗舰模型终审标准，一律以 [`EXECUTION_MASTER_PLAN.md`](EXECUTION_MASTER_PLAN.md) 为唯一最高优先级准则。

本仓库把「读数 → 合并 → E0 → 归一化 → Rbkg → χ(k) → FT → 小波 →
结构识别 → 多壳层 CN 拟合 → 报告」串成一条可复现的 Python 管线。
物理计算基于 xraylarch，结构展开基于 gemmi/FEFF6，不依赖 Athena /
Artemis 的图形界面。

## 项目目录结构与文档导航

| 类别 | 路径 | 职责 |
|---|---|---|
| 核心入口 | `run_xafs.py` | 单样品全流程（处理 / 识别 / 拟合 / 正式 CN） |
| 核心入口 | `calibrate_t5.py` | T5 独立 S₀² 标定（冻结门禁，不放宽） |
| 核心入口 | `run_ci.py` | 回归契约 + 压力测试总入口 |
| 构建期工具 | `tools/cif_crawler.py` | 固定 COD ID → 正式 phase pack（44 条真实 CIF） |
| 构建期工具 | `tools/build_source_pool.py` | T5 谱源隔离池生成与校验 |
| 结构库 | `_structure_lib/` | 诊断 motif + `_phase_pack/` 正式 CIF |
| 标定资产 | `data/s02_calibration/` | T5 标定产物、motifs、source_pool |
| 设计规范 | `docs/REVIEW_AND_REVISED_PLAN.md` | 正式 CN 物理门禁总纲（§1–§6） |
| 执行计划 | `docs/P4_PLUS_EXECUTION_PLAN.md` | P4+/P5/P6 催化能力扩展计划 |
| 验收账本 | `docs/ACCEPTANCE_LOG.md` | P0a–P6.5 完整验收记录 |
| 催化能力/门禁 | `docs/CATALYSIS_CAPABILITY_AND_GATES.md` | 催化剂体系能力、每个输出的分析价值、置信门禁卡、参数台账 |
| 爬虫说明 | `docs/CIF_CRAWLER_INTEGRATION.md` | COD 目标表、P5-C 入库与 T5 谱源整合 |
| 执行记录 | `docs/P6_EXECUTION_INSTRUCTIONS.md` | P6 门禁修复与项目整合的权威执行指令 |
| 执行记录 | `docs/P6_1_EXECUTION_INSTRUCTIONS.md` | P6.1 单晶位/多晶位正式 CN 端到端闭环修复 |

## 安装

```powershell
python -m pip install -r requirements.txt
python env_check.py          # 检查 numpy/scipy/matplotlib/larch/gemmi
```

## 快速开始

```powershell
# 常规处理（默认 strict：边后平台不在 0.85–1.15 就中断）
python run_xafs.py 文件1 文件2 --absorber Fe -o xafs_out/样品名

# 追加结构识别 + 自动拟合
python run_xafs.py 文件 --absorber Fe --fit -o xafs_out/样品名

# 强制指定结构；识别 ambiguous/partial 时也允许 Top-1
python run_xafs.py 文件 --absorber Fe --fit --fit-motif goethite `
    --allow-ambiguous -o xafs_out/样品名

# 多位点混合相：放开组分权重做相分数拟合
python run_xafs.py 磁铁矿.xdi --absorber Fe --fit --fit-weights -o xafs_out/磁铁矿

# 生成五联图 + 拟合表 + CSV
python xafs_report.py 文件 --absorber Fe -o xafs_out/报告名

# 内置演示（lepidocrocite + Fe 箔）
python xafs_report.py
```

`--fit` 的行为：

- 识别 `confident`：按 motif 自动选壳层，做多壳层 CN；
- `partial` / `ambiguous` / `no_match`：默认只做 L1 第一壳；
- 未在 absorber.py 标定 S₀² 的新元素：拒绝 CN 定量；可用运行时
  `--s02 <值>` 显式覆盖，或对 T5 已通过元素用 `--s02-from-t5`
  从独立产物读取（现有 Fe/Co/Ni/Cu/Pt 均已标定；T5 产物见下文，
  不写静态表）；Fe 金属相用 `--chemical-state metal`（0.6566），
  Fe 氧化物/羟基氧化物保持默认 0.803；
- `--fit-weights`（仅 `kind="mixture"`）：固定每路径晶体学简并度，
  用 softmax 放开各组分权重（恒正、和=1），输出 `mixture_weights`
  `[{group, weight, weight_err}]`。此时表里的 `N` 是**有效简并度
  N_crystal·w，不是配位数**；相分数可靠性由 `weights_reliable` 单独
  门禁（Rf/过拟合/ΔE0/权重撞界/权重误差/组分简并）。
- 所有 S7 拟合额外输出 `parameter_names` 与 `correlation`（自由变量
  相关矩阵）：误差棒不是独立的，`N–σ²–ΔE0` 常强相关。实测磁铁矿
  权重拟合最强相关为 ΔE0 vs oct Fe–Fe ΔR = 0.904。
- 荧光模式自吸收筛查：没有样品组成/厚度就不做不可验证的定量校正；
  `Δμ≥1.0` 判 `risk=high` 并拒绝 CN/相分数结论，`0.3–1.0` 只警告，
  `<0.3` 视为稀样。透射/文件自带 `xmu` 不受影响。
- 荧光模式可选 Booth 自吸收校正（S4.5）：给出**样品总化学式 + 厚度 +
  有效密度 + 入射/出射角**后，在 k 空间按 Booth & Bridges 2005 Eq.4
  反解振幅。缺任一参数即拒绝启动，不做假设。校正后仍由 `reliable`
  标记把关；化学式写错会直接把误差带进 N。

## 关键文件

| 文件 | 作用 |
|------|------|
| `xafs_pipeline.py` | S1–S6 归约编排（唯一上游） |
| `structure_id.py` | S6.5 结构识别、置信度、L1 第一壳 |
| `xafs_fit.py` | S7 自动/指定结构拟合、壳层合并、信息量截断 |
| `run_xafs.py` | 数据处理 + 拟合 CLI |
| `xafs_report.py` | 通用报告 CLI + 内置演示 |
| `artemis.py` | CIF → FEFF → R 空间拟合 |
| `athena.py` | 归一化、AUTOBK、FT |
| `larch_wt.py` | 连续小波变换 |
| `calibrate_s02.py` | 标定 S₀²（金属箔/氧化物标样，生产等效协议；`--motif-file` 用标定专用结构） |
| `calibrate_t5.py` | T5 批量标定 12 元素，产物落 `data/s02_calibration/`（不改静态表） |
| `absorber.py` | 吸收元素配置表（边能量/S₀²/kmax），元素泛化唯一入口 |
| `regression_test.py` | 135 项契约回归测试（B1–B4 基线、B5 复数小波/元素 k 标尺、B6 多扫描离群剔除/XANES 自吸收/边前分峰、B7 各向同性晶格约束/部分占据 CIF/C4 与 g(R)、B8 多数据集联合拟合与爱因斯坦模型、B9 配置全键/C4 降级/σ_μ/SNR/MCA/.prj/白线面积、B10 MBACK/chi_std/PCA SPOIL/MCR-ALS、B11 相移 FT-WT/异核共占位/关联 Debye/BVS、B12 三空间拟合/AIC-BIC/双吸收边 ΔE₀/SAC-DAC、B13 SAC/DAC 动态 Motif 与 Hamilton F 检验/异核双分支 path_links、B14 Mini-GDS 化学计量约束/相移校正与 χ² 等高线/关联 Debye 联合精修、B15 Operando 时间序列、B16 Athena Perl .prj 与 Co/Ni/Cu/Pt 催化结构库、B17 纳米粒径反演/双金属合金 Cowley SRO、B18 Rayleigh 分辨护栏/运行时 --s02、B19 原位 EXAFS 序列联合精修、B20 电催化氧化物/羟基氧化物库与真实 BMM 标样盲测、化学态 S₀² 多态（Fe metal 0.6566）、P4+ CIF 抓取/条件 CN/T7 多元素/T8 联合校验、P5-A T5 产物完整性/--only 增量合并/--s02-from-t5 溯源、P5-B Mo 双壳 bcc/Ce σ² 压界拒绝、P5-C 18 个常见相 COD 补齐/coverage 无缺口、P5-prep 标准谱源隔离池/5 元素门禁冻结、P5-D 多位点正式 CN、P6 phase pack 44/44 门禁、P6.1 单晶位正式 CN 与诊断库解耦、P6.2 边界入口对齐、P6.3 双入口配置/operando 守卫、P6.4 参数/门禁/决策/摘要台账与打印开关、P6.5 催化剂结构图包/逐面板 CSV/数据包、P6.6 双入口产物对齐/图包分层瘦身/q 空间与全量结构候选落盘） |
| `xafs_provenance.py` | P6.4 输出台账：参数来源、门禁命中、决策轨迹、结果摘要（只记录，不判定） |
| `xafs_figure_pack.py` | P6.5/P6.6 图包后处理：Figure 1–4 总图、表格图、`figures/panels/` 逐面板 CSV、`figures/data/` 跨面板汇总；默认不重复渲染单面板（只读已有输出，不重新拟合） |
| `run_ci.py` | 跨平台 CI：编译 + 契约回归 + 可选压力测试，输出 ci_report.json/md |
| `stress_test_{co,ni,cu,nio}.py` | 元素泛化/护栏/CN 回收压力测试 |
| `xafs_config.py` | 配置文件支持（.yaml/.json）、批处理入口（--batch-dir） |

## 当前验证结果

**S₀² 标定状态（真实箔，BMM 标样，2026-11）：**

| 元素 | S₀² | 箔相假设 | 数据 | 窗口极差 |
|------|-------|----------|------|----------|
| Fe | 0.803 | 金属 bcc / 氧化物折中 | 生产值：氧化物价态锚定；BMM 箔交叉标定 0.6566（IIT 0.6335） | 0.0205（箔） |
| Co | 0.7308 | hcp | Co-K-CobaltFoil | 0.0005 |
| Ni | 0.7795 | fcc | Ni-K-NickelFoil | 0.0026 |
| Cu | 0.8670 | fcc | Cu-K-CopperFoil | 0.0092 |
| Pt | 0.8075 | fcc (L3) | Pt-L3-PlatinumFoil | 0.0026 |

标定工具：`python calibrate_s02.py <箔数据> --absorber Co --no-strict [--apply]`
（协议与门禁详见 `calibrate_s02.py` 模块注释）。

**Fe 箔重标定（2026-11-10，交叉检查）**：BMM `Fe-K-IronFoil.xdi`，
`fe_bcc` 第一壳 N=8 固定，E0=half_height，3 窗口
S₀²=0.6471/0.6566/0.6676（极差 0.0205），Rf≤0.0422，ΔE₀≈+8.9 eV，
门禁 5/5；IIT-Lytle 5 µm Fe 箔独立交叉检查 0.6335（极差 0.02）。
生产单壳拟合（`--fit-rmax 2.60 --chemical-state metal`）回收 N=7.93±0.49
（晶体学 8）、R=2.517 Å、Rf=0.0219。**该值只对单壳 bcc 模型成立**：
默认双壳拟合会把 N 抬到 11.5/13.4。它不能直接替代生产值 0.803——
0.803 是氧化物价态验收锚定的跨材料折中值（原始流程不可复原），改用
0.6566 会使 Fe 氧化物/羟基氧化物 CN 高约 22%（lepidocrocite 第一壳
N≈7.8）。金属 Fe 研究用 `--chemical-state metal`（等价 `--s02 0.6566`，
来源记为 `chemical_state:Fe=metal`）；氧化物 CN 保持默认 0.803。

金属 Fe 生产链示例：

```powershell
python xafs_report.py foil.xdi --absorber Fe --fit-motif fe_bcc \
  --chemical-state metal --fit-rmax 2.60
```

**T5 独立标定产物（12 个未标定元素，2026-11-09）：**

`python calibrate_t5.py` 对 Cr/Mn/Zn/Mo/Ru/Rh/Pd/Ag/Ir/Au/Ce/U 逐元素
跑批量标定，产物在 `data/s02_calibration/`（`index.json`、逐元素
`calibration.json`、`REPORT.md`）。**不写 `absorber.py`**——#60 要求
静态表中这 12 个元素保持 `s02=None`/`calibrated=False`；运行时用
`run_xafs.py --s02 <值>` 显式覆盖，或对已通过元素用
`run_xafs.py --s02-from-t5 --fit` 读取独立产物（见下文 P5-A/P5-B）。

| 元素 | S₀² | 标样 / motif | Rf | 窗口极差 |
|------|------|--------------|----|----------|
| Cr | 0.6482 | BMM Cr₂O₃ / `cr2o3_corundum` | 0.0074 | 0.0053 |
| Mn | 0.4311 | BMM Pyrolusite / `mno2_rutile` | 0.0402 | 0.0213 |
| Zn | 0.6071 | BMM Zincite / `zno_wurtzite` | 0.0090 | 0.0070 |
| Mo | 0.8322 | BMM Mo foil / `mo_bcc`（双壳 + `half_height` E₀） | 0.0040 | 0.0017 |
| Ru | 0.7002 | BMM RuO₂（`-xmu` 派生）/ `ruo2_rutile` | 0.0154 | 0.0123 |
| Ag | 0.7513 | IIT-Lytle Ag foil / `ag_fcc_cod_9008459` | 0.0030 | 0.0141 |
| Au | 0.8602 | BMM Au foil / `au_fcc` | 0.0003 | 0.0073 |

未通过 5/12（不强行给值）：Rh/Pd/Ir（IIT-Lytle 谱仅约 170–195 eV
above E0，k_max 不足）、Ce（自由 σ² 塌到 1e-4 下界；固定 σ² 时
S₀² 随假设值从 0.54 漂到 0.80；加 Ce–Ce 第二壳 Rf=0.127；备选透射谱
平台 2.27、属归一化伪迹）、U（α-U 撞下界、Rf 差）。Mo 已在 P5-B
通过：双壳 bcc（8×2.725 + 6×3.147 Å）+ `half_height` E₀，
S₀²=0.8322，门禁 5/5；详见下文 P5-B。
`_standards_lib/` 与 `data/standards_raw/` 因 IIT-Lytle 许可未明确
**不提交**；Ag 的标定依赖本地标准，产物记录源文件 SHA256，但不保证
在无本地标准的机器上可复现。

**T7 多元素与 T8 XANES–EXAFS 联合校验（可选开关）：**

```text
python run_xafs.py a.dat b.dat --multi-element --multi-element-s02 "Fe=0.803,Co=0.7308"
python run_xafs.py sample.xdi --absorber Fe --xanes-cross-check
python run_xafs.py Cr-K-Cr2O3.xdi --absorber Cr --s02-from-t5 --identify --fit
```

- `--multi-element`：逐文件从 XDI `Element.symbol` / `Absorber:` 头自动判定
  吸收元素，分元素独立子目录跑批，输出 `multi_element_summary.json/.md`；
  各元素结果**并列**，`recommended=null`、`best_element=null`，不自动
  选择/推荐 Top-1。
- `--s02-from-t5`：从 `data/s02_calibration/<元素>/calibration.json` 读取
  T5 已通过元素的 S₀²；只接受 `status="calibrated"` 且
  `static_absorber_py_unchanged=True` 的产物，未通过元素明确拒绝。
  优先级 `--multi-element-s02` > 全局 `--s02` > `--chemical-state`
  > T5 产物 > 静态表。
- `--chemical-state`：化学态 S₀² 多态，支持 `metal` 或
  `Fe=metal,Co=metal`。当前仅 Fe 金属相有交叉标定 0.6566；默认
  生产值仍为 0.803，静态表只读。未列出的逐元素映射继续走 T5/静态表。
- `--xanes-cross-check`：归一化平台、边步、E0 导数双峰、白线可靠性与位置、
  K 边极端 CN–白线对比、荧光自吸收风险；只写警告到
  `<stem>_meta.json.xanes_cross_check`，`affects_cn_reliable=False`，
  不改变正式 CN 门禁。

真实 BMM Fe 标准（K 边）：

| 标准 | 判定 | 说明 |
|------|------|------|
| Hematite | confident | 正确 |
| Goethite | confident | 正确（COD 2211652 单晶结构） |
| Wustite | confident | 正确（理想岩盐 FeO） |
| Magnetite | confident | 正确——补 8a 四面体位 + `magnetite_mixed` 联合模型后由 partial 升级 |
| Ferrihydrite | no_match | 弱结晶，无可靠体相 motif → L1 第一壳 |

真实 BMM Co/Ni/Cu 标准（K 边）：

| 标准 | 判定 | 说明 |
|------|------|------|
| CobaltFoil | confident | co_hcp，第一壳 N=11.85 vs 晶体学 12（−1.3%） |
| NickelFoil | confident | ni_fcc，第一壳 N=11.80 vs 12（−1.7%） |
| CopperFoil | confident | cu_fcc，第一壳 N=11.60 vs 12（−3.3%）；第二壳 N 撞界 → `cn_reliable=False` 正确拒标 |
| NickelOxide | confident | nio_rocksalt（score 0.025，与 ni_fcc 隔 0.80）；Ni 标定在第二种结构类型上的独立验证 |

真实 BMM Pt 标准（L3 边）：

| 标准 | 判定 | 说明 |
|------|------|------|
| PlatinumFoil | confident | pt_fcc，第一壳 N=11.69±0.22 vs 12（−2.6%），R=2.764 Å，Rf=0.0038 |
| TetraamminePtCl₂ | no_match → L1 | 库仅金属箔；L1 默认 O 拟合 Rf=0.031（N=4.33），`--l1-ligand N` Rf=0.024 |
| K₂PtCl₄ | no_match → L1 | `--l1-ligand Cl` 后 N=3.98±0.14、R=2.302 Å、Rf=0.023（平方平面 CN=4 正确恢复）；默认 O 会 Rf=0.32 并正确判 `cn_reliable=False` |

合成数据压力测试全部通过：

- Co / Ni / Cu 金属箔：N≈11.3–11.9（真值 12），R 误差 <0.01 Å；
- NiO：Ni–O N≈7.6（真值 6）、Ni–Ni N≈11.9（真值 12），Rf≈0.006；
- 护栏：未标定 S₀²、无元素/无文件头、无 s02 均按预期报错。

回归测试：

```powershell
python regression_test.py
python stress_test_co.py
python stress_test_ni.py
python stress_test_cu.py
python stress_test_nio.py
```

## P1 工程化功能（2026-11）

### 配置文件（`xafs_config.py`）
支持 `.yaml` / `.json` 配置文件，CLI 参数优先级高于配置文件：
```powershell
python run_xafs.py <谱文件> --absorber Fe --config xafs.yaml
```
配置项：`kmin/kmax/kweight/post_window/rbkg/strict/wt_sigma/despike`，
以及自吸收子表（formula/thickness/density/angles/line）和批处理（pattern/outdir）。

### 批处理（`--batch-dir`）
对目录下所有匹配文件循环处理，单个失败不中断整批，失败记录写入
`<outdir>/batch_errors.json`。实测 BMM Fe 标准 38 个文件 38/38 成功。

### 去毛刺（`--despike`）
S1.5–S3 之间插入中值滑窗（窗口=5）+ MAD 门禁（阈值=3×MAD），默认关闭。
实测 lepidocrocite 替换 79 个离群点，E0/edge_step/Rbkg 结果与未去毛刺一致。

### FT–WT 联合诊断图（`diagnostics/s7_ft_wt_joint.png`）
该图属于**显式诊断产物**，默认不生成，`--figure-pack` 也不包含它。
需要在命令中加 `--diagnostics`：

```powershell
# 报告入口
python xafs_report.py sample.xdi --absorber Fe --fit-motif lepidocrocite --diagnostics -o xafs_out/demo

# 通用入口
python run_xafs.py sample.xdi --absorber Fe --identify --fit --diagnostics -o xafs_out/demo
```

输出目录为 `<outdir>/diagnostics/`。两个入口都生成 9 张 PNG：
`s1_raw_scans.png`、`s2_e0.png`、`s3_normalization.png`、`s4_rbkg.png`、
`s5_chik_chir.png`、`s5b_backft.png`、`s6_wavelet.png`、
`s7_ft_wt_joint.png`、`s8_wt_sweep.png`，以及配套 `wt_sweep.csv`
（9 行：σ × kweight 的峰位/幅度/相对 FT 偏移）。

`s7_ft_wt_joint.png` 左图 `|χ(R)|` 标注 FT 第一峰 r1；右图 WT contour
叠加水平线 R=r1 与垂直线 k=k_wt，交点标出，直接回答「这个 FT 峰主要
来自哪个 k」。对标 Larch-WT 的交互诊断，避免 WT 在多壳分支间 argmax
跳壳（Fe3O4 实测：Booth 校正把低 k 抬高 7% 后 WT 峰从 2.60 跳到 1.42，
FT 峰纹丝不动；联合图可肉眼核对）。

### 多 kweight / σ 的 WT 对比（S8，`--wt-sweep`）
`xafs_diagnostics.wt_sweep` 在 FT 第一峰 ±0.5 Å 锚定下扫描 σ∈{1.0,1.5,2.0}
× kweight∈{1,2,3}，输出 `s8_wt_sweep.png` + `wt_sweep.csv`（9 行：
k_peak / r_peak / amplitude / offset_from_ft）。`--wt-sweep` 另存
`<stem>_wt_sweep.png/csv` 并内联进 HTML 报告。

色标**按 kweight 行归一**，不是全局归一：kweight=3 的幅度约是 kweight=1 的
k² 倍（k≈6 Å⁻¹ 实测约 36×），全局色标会把低 kweight 整块压成深蓝。
因此同一行内 σ 可比，跨行幅度不可比——这是 k 加权的定义，不是数据处理。
σ=1.0 用于暴露相邻壳粘连；σ=1.5 是产线默认。k 轴用全部面板的并集范围
（`sharex` 下逐面板 `set_xlim` 会被最后一个面板覆盖，曾把 σ=1 的 k<5 裁掉）。

### 交互式 k 窗口（`<stem>_window_tuner.html`）
`xafs_window_tuner.write_window_tuner` 生成单文件 HTML：χ(k) 与 Hanning 窗
FT 内核嵌进页面，浏览器里拖 kmin/kmax/kweight/dk，实时看 `|χ(R)|` 峰位与形状。
后端先用 `larch.xftf` 扫描候选窗口（kmin 2.5–4.0、kmax 9–15），只保留峰位
落在产线 r1±0.15 Å 的候选，按
`score = amplitude/(1+20·stability+5·|R−r1|)` 给出**参考**推荐，不自动改产线，
页面给出对应 CLI 命令。浏览器 FT 与产线同用 Hanning + e^{−2ikR}，幅度按各自
最大值归一，只比较峰位/形状；不改变 Rbkg/归一化。`--no-window-tuner` 关闭。

### σ² 共享组（对标 Artemis 路径约束）
`--shared-sigma2 all` 让所有路径共享一个 σ² 变量，减少自由参数、避免
σ² 与 ΔR 过度相关。或用 JSON 列表显式分组，如 `--shared-sigma2 '[[0,1],[2]]'`。
共享组输出相同的 σ² 和 σ²_err。与 `--fit-weights` 互斥。

### 单文件 HTML 报告（默认生成）
每次 `xafs_report.py` 运行都会生成 `<stem>_report.html`：图片以 base64
内联，可脱离目录单独发送；关键指标与门禁结论置顶，强相关
（`|ρ|≥0.8`）与不可靠 CN 标红；CSV 数据文件保留链接。`--no-html` 关闭。

### 动态 CN 系统误差预算
`cn_systematics` 不再恒为 ±30%，而是按本次拟合的实际条件合成：
基础 0.30 + Booth 校正残留 0.10 / 未校正自吸收 0.15 + 参数强相关
0.05–0.10 + 远壳层截断 0.05 + σ² 撞界 0.10（平方和）。每项来源写入
`extra_components` 与 `note`，报告层直接显示合成后的预算。
相分数另有 `weights_systematics`（结构模型绑定 0.20 + 拟合权重误差）。

### CI 自动化回归
```powershell
python run_ci.py            # 编译 + 135 项契约回归
python run_ci.py --stress   # 额外跑 Co/Ni/Cu/NiO 压力测试
```
产物 `ci_report.json` / `ci_report.md`；任一步失败即非零退出。
`.github/workflows/xafs-ci.yml` 在 push/PR 时自动执行。

### 催化剂结构图包（P6.5 / P6.6）
```powershell
# 已有输出目录，独立补图包（默认只出总图 + 表格图）
python xafs_figure_pack.py xafs_out/RuO2 --formats png,pdf --dpi 300

# 需要逐面板独立图片时显式打开（输出到 figures/panels/）
python xafs_figure_pack.py xafs_out/RuO2 --standalone-panels

# 或在两个入口直接生成
python run_xafs.py sample.xdi --absorber Fe --identify --fit --figure-pack
python xafs_report.py sample.xdi --absorber Fe --fit --figure-pack
```
输出 `figures/`：Figure 1 结构总图（a–h）、Figure 2 置信/过程仪表盘
（a–f）、Figure 3 拟合与误差诊断（a–h）、Figure 4 数据特征与输入清单
（a–f）、Table S1/S2 表格图。目录分层固定为：

- `figures/` 根目录：4 张总图 + `Figure2g_cn_table.png` /
  `Figure2h_fit_table.png` + `figure_manifest.json` / `figure_pack.md`；
- `figures/panels/`：每个面板唯一的 `.csv` 实际绘图数据；加
  `--standalone-panels` 时才额外生成 `panels/<id>.png`；
- `figures/data/`：只放跨面板汇总（`summary.json`、`panel_index.csv`、
  `all_fit_paths.csv`、`fit_variants.csv`、`all_parameters.csv`、
  `gate_checks.csv`、`decisions.csv`、`README.md`），不再复制父目录
  源 JSON/CSV，也不重复 `fit_correlation.csv` / `cif_shells.csv` /
  `structure_candidates.csv`。

`figure_manifest.json` 记录每个面板的图、CSV、源文件与图注；
`figure_pack.md` 为人可读索引。图 1(f) 的 CSV 保留全部结构候选
（图面只画前 6 个），并带 `nuclearity/score_mag` 列。图包只读已有
输出，不重新拟合、不改变任何判定；`--operando` 不支持单结构总图。

### 通用论文主图（`paper_figures.py`）
`xafs_figure_pack.py` 面向“已有输出目录的完整图包”；`paper_figures.py`
只做论文主图 + CN 拟合表，适合跨体系复用：

```powershell
python paper_figures.py --sample sample_mu.dat --absorber Fe --edge K `
  --ref "Fe foil=foil_mu.dat" `
  --kmin 3 --kmax 12 --kweight 3 --rbkg 1.3 `
  --guide "Fe-O=<本图峰位>" --guide "Fe-Fe=<本图峰位>" `
  --fit-json sample_fit.json --out paper_out
```

有 `--fit-json` 时输出六面板 (a) XANES (b) kⁿχ(k) (c) |FT| (d) R 空间
拟合 (e) k 空间拟合 (f) CN 表，并另出 `*_cn_table.{png,csv,md}` 与
`*_ft.csv`；无拟合时只出 (a)–(c)。模块不绑定元素/边/壳层：标线由
`--guide LABEL=X` 或 `LABEL=X..Y` 传入，S₀²、E0、拟合窗口原样读取，
不重新拟合。CN 表逐行标注 formal / diagnostic / upper limit /
σ² at bound，诊断级配位数不会被写成正式值。图注、CN 表、CSV 与 fit JSON
会同时披露本次拟合采用的**结构与 CIF 来源**（结构类型、来源描述、phase ID、
吸收位点、CIF 路径与 SHA256、FEFF 输入与 SHA256）。若拟合用的是手工构造
团簇而非 CIF 派生相，输出会明确写 `CIF: none (no CIF used)`，不会拿无关
CIF 顶替。可用 `--structure-json` 传入或覆盖 provenance；模块 docstring
内是最小调用示例；本机 Pt 案例的完整产物和可复现脚本在
`xafs_out/Pt1CeO2-flu_paper/`（分析产物不入库），通用检查表见
`xafs_out/Pt1CeO2-flu_analysis/lessons_learned.md`。

当前目录缺少 `data/standards_raw/`、`_standards_lib/` 或需要补齐 COD CIF 时，
先运行构建期爬虫：

```powershell
python cif_crawler.py --all    # T3/T4 CIF + T5 金属箔标准谱
python cif_crawler.py --t5     # 只补标准谱
```

T3/T4 会优先复用 `data/cif_raw/` 的已验证 CIF，缺失才联网取 COD；T5 从
IIT-Lytle 下载 12 个元素的标准谱并写入 `_standards_lib/`。这些是本机构建
资产，已在 `.gitignore` 中，不随仓库提交。运行期相解析仍以本地 cache /
phase pack 为主；`--offline` 可显式禁止外部 fetcher。

### 批处理输出目录
`--batch-dir` 的输出优先级：显式 `-o` > 配置文件 `batch.outdir` >
默认 `xafs_out_batch`。

## 四批次对标补齐（B1–B4，2026-11）

> 数学契约与验收标准以 `EXECUTION_MASTER_PLAN.md` 为准；本节只列
> 已落地接口与实测数字。所有增强均为**显式开启**，默认路径与历史
> 行为逐位一致（默认拟合 Rf=0.027281、第一壳 N=6.3486、R=2.00196 Å）。

### B1-1 多 kweight 联合拟合
`--fit-kweights 1,2,3`（API：`fit_identified(..., kweights=(1,2,3))`）
把各 kweight 分支的 χ(R) 实部/虚部同时放入残差。每个分支先除以该
分支实验谱 RMS 范数 `s_w`，避免高 kweight 因数值尺度大而垄断代价
函数；`N_idp` **不乘**权重个数（同一批数据不因加权方式产生新信息）。
单权路径保持历史逐位一致。

### B1-2 Back-FT χ(q)
`athena.back_ft_to_q_space(r, chi_r, rmin, rmax, dr)` 对选定 R 窗做
反傅里叶变换，两个入口都输出 `<stem>_chiq.csv`；用于检查第一壳在 q 空间
的相位/振幅，不引入新的拟合变量。

### B1-3 受约束多散射（MS）
`--include-ms`（API：`fit_identified(..., include_ms=True)`）只纳入
`nleg∈{3,4}`、`reff≤rmax+0.15 Å` 的近共线路径；ΔR/σ² 按
`reff_ratio`、N 按 `degen_ratio × N_parent` 由母单散射壳层驱动，
**不新增自由变量**。所有中间原子（含闭合腿）满足 `|cos|≥0.85`，
并加 `reff ≥ 1.5×第一壳距离` 的聚焦判据。MS 标签带完整原子序列与
FEFF 路径号，避免同名路径互相覆盖。NiO 共线 MS 回归已通过。

### B2-1 XANES flattened 与颈部特征
`athena.extract_xanes_features` 在 `[E0+100, E0+400]` 平台均值归一后
给出 flattened 谱；larch 分支用 `g.flat`，native 分支对 `mu_norm`
拟合二次背景后同样归一。输出平台值、白线强度/位置与颈部特征
（双台阶位置/高度）。这是**补充功能**，不参与 EXAFS 定量门禁。

### B2-2 LCF / PCA
`xafs_lcf.run_lcf` 对 XANES 或 χ(k) 做线性组合拟合（权重和=1 可选），
`xafs_lcf.run_pca` 给出主成分方差与目标谱重构残差。CLI：
`--lcf-standards a.xdi,b.xdi`。结果写入报告 JSON 与 HTML；LCF 是
端元线性组合，不是独立物相定量，权重解释必须附端元选择依据。

### B2-3 死时间与显式 deglitch
`--deadtime` 仅对荧光模式生效：`dtime≥1` 视为倍乘因子，
`0<dtime<1` 视为死时间分数；透射模式直接跳过。
`--deglitch '7250:7262,7410:7415'` 直接删除能量区间内的点，**不插值**；
删除比例 >15% 或边核心剩余 <5 点抛 `XAFSConfigError("S1.2", ...)`。
去毛刺（`--despike`）仍是独立的中值滑窗 + MAD 门禁，默认关闭。

### B3-1 WT 壳层一维切片与实验/模型残差
- `larch_wt.extract_k_profile(wave_map, r_center, half_width=0.25)`：
  返回 `ridge_profile`（R 带内嵴线）、`slice_profile`（最近 R 单行）、
  `k_peak`、`r_track`、`bimodal`/`k_secondary`/`peak_separation`。
- `larch_wt.compute_model_wt_residual(chi_k, model_k, ...)`：实验与模型
  用**完全相同**的 WT 参数（`normalize=False`），共享同一物理标度后
  相减，给出二维 WT R-factor。
- `xafs_diagnostics.write_wt_shell_profiles`：输出
  `<stem>_wt_profiles.png`（各壳层切片 + 二维残差）与 `.csv`；
  `--wt-sweep` 时自动内联进 HTML 报告。

实测 LEPIDO（kmin=2.5 分析窗）：Fe–O `k_peak=5.5 Å⁻¹`、Fe–Fe
`k_peak=6.4 Å⁻¹`（差 0.9），二维 WT R-factor=0.046。产线默认
kmin=3.0 时差约 0.7，排序不变。

### B4-1 ε_R、χ²_ν 与增减壳层 F-test
- `artemis.estimate_r_noise`：对 `[15,25] Å` 高 R 区实部/虚部取均方根
  得 ε_R（IFEFFIT 口径），显式要求 `xftf(rmax_out=25)`，否则 larch
  默认 R 网格只到 ~10 Å，区间为空。
- `artemis.shell_f_test(fit_simple, fit_complex)`：
  `F = [(χ²₁−χ²₂)/(p₂−p₁)] / [χ²₂/(N_idp−p₂)]`，
  `confidence = F.cdf(F; dfn, dfd)`；`confidence≥0.95` 且 `dfd≥1`
  才判新增壳层显著。自动 S7 在同一 R 窗生成单壳对照并调用。
- 报告指标新增 `epsilon_r` / `reduced_chi_square` /
  `outer_shells_f_test_confidence`，拟合 JSON 同步
  `epsilon_r/chi_square/reduced_chi_square` 与 `shell_f_test`。

实测默认 LEPIDO（3 扫描）：ε_R=0.009899、χ²_ν=314.2、
F=32.34、confidence=0.9991。单扫描对照：ε_R=0.01064、
χ²=1405.43、χ²_ν=273.32、F=32.67、confidence=0.99909。

### B4-2 三阶累积量 C3（非谐）
`--vary-c3 first|all`（API：`fit_identified(..., vary_c3=True/[0])`）
放开路径三阶累积量 C₃（Å³）。边界固定为 `±artemis.C3_LIMIT=±0.01 Å³`；
`|C₃|≥0.008` 标 `c3_at_bound` 并令 `cn_reliable=False`（非谐项在
吞噬相位失配）。若放开后 `n_var ≥ N_idp`，先自动降级为仅第一壳，
仍超限则彻底关闭并记录 `flags.c3_dropped_by_nidp=True`。
C3 与 mixture 权重拟合互斥（权重模式的有效简并度不是配位数）。
`PathFit.c3/c3_err/c3_at_bound` 全程序列化；拟合表只在 C3 真的
放开时增加 `C3 (Å³)` 扩展列，默认表仍为 7 列。

合成验证：C3 真值 `+1.5×10⁻⁴ Å³`，rmax=2.6 时精确回收；
rmax=3.5 回收 `1.808×10⁻⁴`（误差 3.08×10⁻⁵）。真实 LEPIDO
`--vary-c3 first` 得 C3=+1.03×10⁻³ Å³、R=2.046 Å、CN=5.5（仅在
显式开启时改变；默认仍为 R=2.002 Å、N=6.35）。

## Phase II 进阶补齐（B5–B8，2026-11）

> 全部显式开启，默认路径不变。Phase II–V 合计 **74/74 回归 +
> 4 压力测试** 全绿（`ci_report.md`）。

### B5-1 复数小波与相位盲区
`larch_wt.continuous_wavelet_transform` 新增 `complex_map` 与 `phase`；
`WaveletMap` 同步归一化。`larch_wt.compute_complex_wt_residual` 对
实部/虚部同时比较，反相模型 `χ_model=−χ_data` 的
`complex_r_factor=4.0`（幅值残差 0.0），证明幅值相等、相位反 π 的
退化不再漏检。诊断层 `write_wt_shell_profiles` 输出
`data_phase_ridge` 与 `complex_r_factor`。

### B5-2 候选元素理论 k_peak 标尺
`structure_id.reference_element_k_peaks("Fe", elements=("O","Cl","Fe"), ...)`
用缓存 FEFF 路径生成理论 WT 峰位；实测 `k_peak(O/Cl/Fe)=4.9/5.6/7.0 Å⁻¹`，
Fe 的 σ² 从 0.003 增至 0.015 时峰位左移 1.4 Å⁻¹。该标尺用于给 WT
壳层切片提供元素判据，不替代拟合。

### B6-1 多扫描离群剔除
`athena.merge_scans(..., reject_outliers=False, z_threshold=3.5)`：
默认关闭且与历史均值逐位一致；开启后按中位谱 RMS+MAD 剔除离群
扫描，保底保留 2 条。`RawScan.sigma_mu` 传播逐点标准差；CLI
`--reject-outlier-scans`。

### B6-2 XANES 荧光自吸收（Tröger 反解）
`athena.fluo_xanes_self_absorption(norm, formula, elem, edge, line,
thickness_um, density, anginp=45, angout=45)`：
`μ_true = α·μ_exp/(α+1−μ_exp)`，分母下界 `0.25(α+1)`；只改 XANES，
不动 `chi_k`。Fe₂O₃/Fe K/45°–20° 实测 α≈0.673、白线误差 0.0、
LCF 权重回到 0.5/0.5；默认 45°/45° 的 α≈0.42 触发 74 点钳制并
`reliable=False`（诚实降级）。

### B6-3 边前 Pseudo-Voigt 分峰
`athena.fit_pre_edge_peaks` / `_pseudo_voigt_area`：主边拖尾基线 +
单/双 Pseudo-Voigt。合成双峰（7113.0/7114.6 eV）质心回收
7114.1025 eV（真值 7114.0667，误差 0.036 eV），Rf=0.00123。
指标表新增 `pre_edge_centroid / pre_edge_n_peaks /
pre_edge_voigt_r_factor`。

### B7-1 各向同性晶格膨胀约束
`fit_identified(..., isotropic_deltar=True)` / CLI `--isotropic-deltar`：
所有 SS 路径共用唯一 `(-5,"alpha")` 槽位，`ΔR_i=α·Reff_i`，边界 ±0.08。
合成 NiO 四壳真值 α=0.004：变量 13→10、α 回收 0.00399999999995、
Rf=4.2e-23，四壳 `ΔR/Reff` 与 α 严格一致（<1e-12）。

**物理边界**：实验 NiO 的各壳自由位移为 Ni-O −0.0098、Ni-Ni +0.0030、
Ni-O −0.0794、Ni-Ni +0.0983 Å，不服从统一晶格缩放；四壳约束最优
α=0.00091、Rf=0.135（自由拟合 0.024）。因此实验数据必须
`cn_reliable=False` 并给出“各向同性晶格约束 Rf>0.08”警告，不能用
该约束下的 N 作配位数结论。

### B7-2 CIF 部分占据与共占位
`artemis.read_cif_sites(..., return_occupancy=True)` 返回 5 元组
`(x,y,z,symbol,occ)`，默认仍 4 元组；距离 <0.08 Å 共占位合并，
元素取占据度较大者、有效占据度 `min(1,Σocc)`；吸收中心按
`occ×多重度` 选主位点。`xafs_fit._ss_rows` 用
`cif_shell_occupancy_scale` 按距离/元素匹配平均占据度，把 FEFF
几何简并度缩放为 `N_cryst=N_geom·mean(occ)`。合成 Pm-3m 验证：
第一壳 O occ=0.75 → `degens[0]` 简并度 6×0.75=4.5；第二壳 Co/Ni 共占位合并为
Co occ=1.0（简并度 8×1.0=8.0），FEFF 不崩溃。

### B7-3 C4 累积量与非高斯 g(R)
`PathFit.c4/c4_err/c4_at_bound`、`artemis.C4_LIMIT=2e-4 Å⁴`、
`set_path_parameters(..., fourth=...)`、`fit_r_space(..., vary_c4=True)`。
合成 C4=6e-5 精确回收（误差 6e-17），Rf≈5e-25。
`artemis.reconstruct_pair_distribution(r, n, r0, sigma2, c3, c4)` 用
Gram-Charlier/Edgeworth 展开、非负截断与梯形重归一；C3=2.5e-4 时
面积=6.0、`r_mode=2.032 < r_mean=2.05`、`skewness=0.538`。

**判据说明**：规范原文“C3>0 时右半侧积分面积大于左半侧”与归一化
正偏态分布数学矛盾——在均值处切分必然左面积>右面积（实测
3.159 vs 2.779）。回归 #44 改用等价尾长判据：
`q95−r0=0.142 > r0−q05=0.118`，并断言左面积>右面积以锁定归一化
正确性；不伪造面积。

### B8-1 多数据集/变温全局联合拟合
`artemis.einstein_sigma2(T, theta_e, reduced_mass_amu, sigma2_static)`
（`HBAR2_OVER_2KB=24.254316 amu·Å²·K`，T≤0 抛 `XAFSConfigError`）；
`artemis.fit_multi_datasets(datasets, shared_bonds, einstein_models)`：
每个数据集残差除以自身 RMS 范数 s_d，跨数据集绑定 ΔR/σ²/N，
σ²(T) 由 σ²_static+Θ_E 两个全局参数驱动。合成 LEPIDO 第一壳
（真值 Θ_E=420 K、σ²_static=0.0015 Å²、ΔR=−0.02 Å、N=6）双温度
80/350 K：全局变量 8→5、Θ_E=419.99999999 K、σ²_static=0.0015、
两数据集 Rf≈1e-22。`ΔE₀` 在所有数据集共享一个全局变量；双吸收边
如需独立 ΔE₀，应拆成两次单边拟合后用 R/σ² 约束对齐。

## Phase III 高阶能力补齐（B9–B12，2026-11）

> 全部显式开启，默认路径零退化：LEPIDO 三扫描默认仍
> `r_factor=0.027281059410161693`、第一壳 `N=6.348637726`、
> `R=2.0019558445 Å`、7 列拟合表逐位不变。CI 实测 **74/74 回归 +
> 4 压力测试** 全绿。

### B9-1 配置全键同步 / C4 的 N_idp 降级 / g(R) 导出
`xafs_config.DEFAULTS` 覆盖 27 键（含 `vary_c4`、`fit_space`、
`bvs_restraint`、`mback`、`chi_std`、`clamp_lo/hi`、`target_transform`、
`mcr_als`），`apply_config_to_kwargs` 优先级为 CLI 显式 dest >
配置文件 > 代码默认值。`fit_identified(..., vary_c4=True|"all")`
在变量数触顶时先降 C4、再降 C3，并分别记录
`flags.c4_dropped_by_nidp` / `flags.c3_dropped_by_nidp`。

实测：合成 Fe–O 单壳真值 `C4=6×10⁻⁵ Å⁴` 回收 `5.96×10⁻⁵`；
LEPIDO 4 壳层 `vary_c4="all"` → `n_variables=14 < N_idp=15.007`，
但仍先降为仅 `Fe-O@2.04` 一层并置 `c4_dropped_by_nidp=True`。
`xafs_fit.rdf_table` 的 g(R) 面积严格等于拟合 N（5.995/5.995、
7.436/7.436），`write_rdf_csv` 输出 7 列 `_rdf.csv`。

### B9-2 σ_μ 误差带 / 高 k SNR / WT 标尺 / L3 白线
`RawScan.sigma_mu` 随三扫描合并传播并写入 `_mu.dat` 第 4 列；
`athena.high_k_snr` 用 `σ_χ=σ_μ/edge_step` 与滚动 RMS 给
`k_snr2/reliable_kmax`。注入噪声标尺 1e-3/3e-3/1e-2/3e-2 →
`k_snr2=13.0/9.5/4.5/4.15`；无 `sigma_mu` 时返回 `None` 并说明原因。
`xafs_report.auto_ref_k_peaks` 从拟合路径 ∪ {O, 吸收元素} 自动生成
WT 标尺（LEPIDO 实测 `{'O':4.9,'Fe':7.0}`，同参数命中缓存）。
`athena.extract_xanes_features` 按边分流：K 边保留
`pre_edge_n_peaks`；L3 边不再拟合 1s→3d 边前峰，改给
`white_line_area / white_line_area_norm / whiteline_peak_height`
（Pt 箔实测 4.5073 / 5.3414 / 0.5427）。

### B9-3 .prj 容器 / MCA 坏道剔除
`xafs_io.write_prj/read_prj` 使用 gzip-JSON（魔数 `1f 8b`），
energy/mu 逐位 round-trip；兼容未压缩 JSON，对 Athena 原生 Perl
`.prj` 干净拒绝。`athena.merge_mca_channels` 用通道间 MAD 鲁棒
z 分数剔除 Bragg spike，再按通道离散度加权合并：8 道中第 2、6 道
注入 spike → `rejected=[2,6]`、`n_kept=6`、权重和=1，合并谱与
干净道加权和逐位一致。

### B10-1 MBACK 归一化 / AUTOBK 先验与钳制
`athena.normalize_mback` 用 Cromer–Liberman 理论截面拟合比例因子
与低阶基线，适用于边后仅 60–150 eV 的短窗口。Fe K 合成谱实测
`edge_step` 相对误差 `5.8×10⁻¹⁶`、`plateau=1.1353`；同窗口传统
二次多项式归一化误差 `2.29%` 且平台 `NaN`。`a≤0`（理论与实测边跳
反向）显式拒绝。
`athena.extract_chi_k(..., chi_std, k_std, clamp_lo, clamp_hi)`
把标准谱作为 AUTOBK 先验：合成谱低 R（<1.2 Å）“提取−真值”泄漏
从 `0.316` 降到 `0.0618`；`clamp_lo/hi` 改变样条端点权重且保持
有限；长度不一致的 `chi_std/k_std` 抛 `XAFSConfigError`。

### B10-2 PCA 目标转换 / MCR-ALS 盲解谱
`xafs_lcf.run_pca` 返回 `components/scores/loadings`；
`target_transform` 按 Malinowski SPOIL 分级（good/fair/poor/
unacceptable）。3 端元 + 9 条混合谱（0.5% 噪声）实测：真实端元
SPOIL≈0.017–0.018，随机谱≈17。`mcr_als` 用 SVD/SIMPLISMA 初始化 +
精确 NNLS 交替最小二乘，`closure=True` 时 621 次迭代收敛，
组分相关 >0.9997、浓度 RMSE≈0.0069、`r_factor=0.01953`。

### B11-1 相移校正 FT/WT / 异核共占位双分支
`artemis.path_total_phase(datfile, k)` 从 FEFF `feffNNNN.dat` 的
`pha` 段取总相移（失败时 Hilbert 兜底并标记来源）。
`athena.ft_to_r_space(..., phase_correct=φ)` 与
`larch_wt.continuous_wavelet_transform(..., phase_correct=φ)` 把表观
峰拉回真实 Reff：Fe–O 实测 FT 峰 1.534→1.9942 Å（Reff=1.9808，
误差 0.0134），WT 峰 1.48→2.0 Å（误差 0.019）。
`artemis.cif_hetero_cocupancy_branches` 按 Reff 距离聚壳，同一晶位
按元素拆支、几何简并度按不同晶位数计：Ti(0.7)/Ru(0.3) 共占 1b 位
→ `n_geom=8`、Ti degen=5.6、Ru degen=2.4、占据度和=1，
`single_element_approx=False`（不再把展开条目重复乘占据度）。

### B11-2 关联 Debye / BVS 软约束
`artemis.correlated_debye_sigma2(T, Θ_D, Reff, m₁, m₂,
k_debye=0.35 Å⁻¹, sigma2_static=0)`：
`σ²=σ²_static+[1−sin(2k_D R)/(2k_D R)]·σ²_un`。实测 300 K 近壳
0.003406、R=4 Å 远壳 0.010127；1000 K 远壳/Einstein=1.309；
T→0 收敛到零点平台 0.0010828，`k_debye_source` 显式记录默认来源。
`artemis.bond_valence_sum(paths, absorber=...)` 计算
`BVS=Σ N_i·exp[(R₀−R_i)/b]`（b=0.37，缺 R₀ 时抛错）。
Fe–O 6×2.0 Å → 3.1280；拟合层 `bvs_restraint={"target":3.0,
"sigma":0.05}` 为软约束：实测 N=4.8567、R=1.9371、
`bvs=3.00097`、`bvs_chi2=3.79×10⁻⁴`、Rf=0.0109，且 `fit.bvs`
与用最终路径重算的手工 BVS 逐位一致。

### B11-3 / B12-1 三空间拟合 / AIC-BIC / 双吸收边 ΔE₀
`fit_r_space(..., fit_space="r"|"q"|"wt")`、`FitResult.aic/bic/
fit_space`；`fit_multi_datasets(..., share_e0=True|False|"auto")`
按 `edge_group` 分组，`"auto"` 在两组 ΔE₀ 差异显著时拆分。
`artemis.profile_chi2_contour` 输出二维 Δχ² 等高线
（1σ/2σ/3σ=2.30/6.18/11.83）。合成单壳三空间 N/R/σ² 最大差
`1.2×10⁻⁹ / 1.1×10⁻⁹ / 1×10⁻¹⁷`；双吸收边真值 +2.0/−6.0 eV：
`"auto"` 精确回收且 `global_r_factor≈2.1×10⁻²⁹`，`True` 只剩
单组 ΔE₀ 且 Rf≈0.0203；`aic_gain > bic_gain`（BIC 惩罚更重）。

### B12-2 SAC/DAC 团簇 / 元素扩充
`artemis.sac_cluster_sites` 生成 M–N₄ 平面（可加轴向配体与
R≈ring_radius 的 8 个 C 载体壳）；`diatomic_m1m2_sites` 生成
M1–M2 双原子中心 + 双配体壳。实测 M–N₄+C 共 12 位点
（4×N@2.0、8×C@2.9），FEFF 生成 9 条路径、第一壳 Reff=2.0；
DAC 9 位点、M2 在 2.5 Å、8×N。
`absorber.py` 新增 Cr/Mn/Zn/Mo/Ru/Rh/Pd/Ag/Ir/Au/Ce/U 共 12 元素，
全部 `s02=None, calibrated=False`（未标定即拒绝 CN 定量）；
Fe/Co/Ni/Cu/Pt 的标定 S₀² 逐位未改。

### Phase III 新增 CLI
`--vary-c4 first|all`、`--fit-space r|q|wt`、
`--bvs-target/--bvs-sigma/--bvs-b`、`--mback`、`--chi-std`、
`--clamp-lo/--clamp-hi` 已同时接入 `run_xafs.py` 与
`xafs_report.py`；B9–B12 其余增强通过 API 调用。
`python run_xafs.py --help` / `python xafs_report.py --help` 可核对
Phase I–III 全部开关。

## Phase IV 催化化学全面分析收官（B13–B16，2026-11）

> 全部新能力显式开启；不传新开关时默认路径零退化：LEPIDO 三扫描仍
> `r_factor=0.027281059410161693`、第一壳 `N=6.348637726`、
> `R=2.0019558445 Å`，7 列拟合表逐位不变。CI 实测 **74/74 回归 +
> 4 压力测试** 全绿（`ci_report.md`）。

### B13-1 SAC/DAC 动态 Motif 与 M₁–M₂ 金属键 F 检验
`structure_lib.build_sac_motif` / `build_dac_motif` /
`parse_catalytic_cluster_spec` 把 B12-2 的团簇构造器升级为可拟合
Motif，FEFF 缓存固定在 `_structure_lib/_dynamic_cache/<hash>/`
（不随临时目录退出而失效）。CLI `--catalytic-model
"dac:Co@2.50,N4@2.00"` / `"sac:N4@2.00,axial=O@2.10,C@2.90"` 绕过
晶相库搜索直接进入 S7。DAC 拟合完成后自动剔除 M₁–M₂ 路径、在同一
窗口重跑纯配体壳层，用 Hamilton F-test 输出 `mm_bond_ftest`
（`f_stat` / `p_value` / `significant` / `delta_r_factor`）：真
Fe–Co 谱回收 R=2.50±0.03 Å、N≈1 且 p<0.01；纯 SAC 谱 p>0.05 不
误报；SAC 轴向氧谱回收 N_N≈4、N_O≈1（误差 <0.5）。

### B13-2 异核共占位双分支与 `path_links`
`artemis.fit_r_space(..., path_links=[...])` 支持同一数据集内跨元素
分支共享 ΔR/σ²（`share_deltar` / `share_sigma2`）并按 `degen_ratio`
绑定 N；`xafs_fit.select_shell_paths(..., hetero_cocupancy=True)`
（CLI `--hetero-cocupancy`）自动把共占位 CIF 拆成双元素分支并生成
链接。Ti0.7Ru0.3 共占 N_geom=8：简并度 5.6/2.4，自由变量 7→5；
合成谱回收两分支 N，ΔR/σ² 及其误差逐位相等，双分支 Rf 比单元素
近似小 5 倍以上。

### B14-1 Mini-GDS 化学计量约束
`artemis.normalize_constraints` + `fit_r_space(..., constraints=...)`
（CLI `--constraints`）支持 `N(Fe-O)+N(Fe-Cl)=6.0`、
`N(Fe-Cl)=1.5*N(Fe-O)`、`dR(Fe-Cl)=dR(Fe-O)+0.02` 与逐参数边界
覆盖；严格禁止 `eval`/`exec`，循环依赖或未知路径标签抛
`XAFSConfigError`。n_sum 使自由变量真实减 1，N 之和逐位守恒，
σ(N) 按系数严格传递。

### B14-2 相移校正 FT/WT 与 χ² 置信椭圆
`--phase-correct` 生成 `<stem>_ft_phase_corrected.csv/.png`（未校正
vs 相移校正 |χ(R)|，并标注拟合 R₁ 竖线）；`--chi2-contour` 生成
`<stem>_chi2_contour.csv/.png`（15×15 Δχ² 曲面 + 1σ/2σ/3σ 等高线
及最优解星标）；`xafs_report.py` 把两张图内联进 HTML。LEPIDO 实测
表观峰 1.4–1.7 Å，校正后与 S7 键长差 <0.06 Å；默认关闭时产物与
耗时不变。

### B14-3 关联 Debye 多数据集联合精修
`artemis.fit_multi_datasets(..., debye_models=[...])` 用单个 Θ_D
同时约束近壳与远壳、跨温度数据集，返回 `debye_results`。合成
80/350 K 双温度数据回收 Θ_D 误差 <15 K、σ²_static 误差
<5×10⁻⁴ Å²，同一温度下远壳 σ² 严格大于近壳，
`global_r_factor < 1e-4`。

### B15-1 Operando 时间序列一键解谱
`xafs_lcf.run_operando_series`（≥3 条谱）串联 PCA 定组分 →
SPOIL 筛标样 → LCF → MCR-ALS 盲解 → E₀/半高/白线/FT r₁ 轨迹；
MCR 组分按主峰能量升序消除排列随机性。
`xafs_diagnostics.write_operando_diagnostics` 输出
`<stem>_operando_concentrations.csv`、`_components.csv` 与四宫格
`_summary.png`。CLI `--operando` 已接入 `run_xafs.py` /
`xafs_report.py`。实测零均值噪声伪标样 SPOIL≈16.95 判
unacceptable 并被 `rejected_standards` 剔除；LCF 权重 RMSE<0.02，
MCR 浓度随序号单调（|r|>0.98）。

### B16-1 Athena 原生 Perl `.prj` 只读导入
`xafs_io.read_athena_perl_prj(path_or_text, as_scans=False,
include_derived=False)` 按行首 `[record]` 分块，解析 `@args` /
`@x` / `@y` / `@i0` / `@signal`；默认过滤 `Ref ` 参比通道与
`is_merge` 派生记录。`read_prj` 仅当文本同时含 `@args =` 与非空
`@x`/`@y` 时才转发，任意 `$VAR1` Perl 片段仍抛
`XAFSConfigError`（保护 #50）。仓库自带 `Fe_standards.prj` 默认
解析 4 条主扫描、`include_derived=True` 共 10 条；`.100` 与原始
文件 462 点逐点一致。

### B16-2 Co/Ni/Cu/Pt 催化结构库扩充
`_structure_lib/build_seed.py` 新增 `con4_planar`、`nin4_planar`、
`cun4_planar`（平面 M–N₄ + 4×C 第二壳）、`cu2o_cuprite`
（`P n -3 m:1`：Cu–O 2×1.849 + Cu–Cu 12×3.019）与 `ptcl4_square`
（Pt–Cl 4×2.31 + 4×Cl@3.45）。全部可 `ensure_feff`。真实
`Pt-L3-PlatinumChloride.xdi`（Sample.name=PtCl2）识别
`status=confident`、top1=`ptcl4_square`（score=0.0864）、S7
`Rf=0.0054`、第一壳 N=3.54、R=2.3099 Å。新增 Motif 第一壳
1.85–2.31 Å 与金属箔 2.49–2.77 Å / NiO 2.95 Å 完全错开，
Co/Ni/Cu/NiO 压力测试 top1 不变。

## Phase V 催化剂精细结构深度解析（B17–B20，2026-11）

> 全部新能力显式开启；不传新开关时默认路径零退化：LEPIDO 三扫描仍
> `r_factor=0.027281059410161693`、第一壳 `N=6.348637726`、
> `R=2.0019558445 Å`，7 列拟合表 SHA256
> `E098D0ED…7532E760` 逐字节不变。CI 实测 **74/74 回归 + 4 压力测试**
> 全绿（`ci_report.md`）。

### B17-1 纳米颗粒尺寸效应反演
`artemis.estimate_nanoparticle_size(n1, r1, n1_err=None, n2=None,
lattice="fcc", shape="spherical")` 用 Calvin/Beale 球截断方程
`N(r,d)=N_bulk[1−1.5u+0.5u³]`（`u=r/d`）经 `brentq` 反解粒径：
`N₁≥0.99·N_bulk` 返回 `regime="bulk_limit"`（`diameter=inf`、
`dispersion=0`）；否则输出 `diameter_nm/diameter_A/n_atoms/dispersion/
n2_spherical_expected/shape_hint`。fcc/hcp 取 `N_bulk=12`、bcc 取
`N_bulk=8`；第二壳显著低于球形预期（`n2<0.65·N₂,sphere`）判
`raft_or_flattened`。S7 仅在拟合壳层全部为同核金属散射且
`N₁<0.96·N_bulk` 时自动挂载 `out["nanoparticle_size"]`，体相箔保持
`None`。实测 `N₁=8.0, r₁=2.774 Å` → `d=1.2270 nm`、`n_atoms=64.1`、
`dispersion=0.9995`。

### B17-2 双金属纳米合金与 Cowley 短程有序
`structure_lib.build_alloy_motif(absorber, partner, x_partner=0.5,
r_aa=None, r_ab=None, n_total=12, lattice="fcc")` 在 fcc 立方八面体
12 方向（bcc 为 8 顶点方向）上按 `x_partner` 放置 A–A / A–B 两条
路径，缓存固定为 `_structure_lib/_dynamic_cache/alloy_<hash12>`，
不进入主库、默认识别零影响。`parse_catalytic_cluster_spec` 新增
`alloy:` 字符串（`"alloy:Co@2.64,Pt@2.76,x=0.25,n=12"`）与
`{"type":"alloy",...}` 字典。`artemis.compute_cowley_sro` 输出
`alpha_cowley = 1 − p_AB/x_partner`、`alpha_err`、`alloying_extent_pct`
与 `tendency`（`|α|≤0.15` 随机固溶体、`α<−0.15` 异核有序、
`α>0.15` 同核团聚/偏析）。合成 PtCo 谱回收
`N_Pt-Pt=5.915`、`N_Pt-Co=2.142`、`α=−0.0636`、`d=1.224 nm`。

### B18-1 Rayleigh 分辨极限与轻元素 Z 简并护栏
`artemis.check_shell_resolution_and_z_degeneracy(shells, kmin=3.0,
kmax=12.0)` 按 `ΔR_Rayleigh = π/[2(kmax−kmin)]` 检查任意两壳；
当 `|ΔZ|≤2` 且 `ΔR<ΔR_Rayleigh` 时给出
`light_element_z_degeneracy`（重元素为 `sub_rayleigh_split`）告警并
建议绑定 `shared_sigma2`/`constraints` 或 BVS 消歧。k∈[3,12] 时
`ΔR_Rayleigh=0.17453 Å`；Fe–N@1.98 vs Fe–O@2.06 命中，Z 差大或
ΔR 超限不误报。`fit_identified` 返回 `resolution_check` 与
`flags.sub_rayleigh_z_degeneracy`，未绑定约束时追加物理告警。

### B18-2 运行时 `--s02` 与 Au/Mo/Cr/Zn 金属箔解锁
`structure_id.validate_s02_override` / `absorber_with_s02` 以
`dataclasses.replace` 构造局部 `AbsorberConfig` 副本，合法区间
`[0.4, 1.3]`，绝不修改 `absorber.py` 静态表（12 个未标定元素仍
`s02=None`，#60 契约不变）。`--s02` 已同步接入 `structure_id`、
`xafs_fit`、`xafs_config.DEFAULTS`、`run_xafs.py` 与
`xafs_report.py`；`_structure_lib/build_seed.py` 新增 `au_fcc`、
`mo_bcc`、`cr_bcc`、`zn_hcp` 与 `calibrate_s02.DEFAULT_MOTIF`。
真实 Au 箔 `s02=0.82`（`rbkg=1.2`）：top1=`au_fcc` confident，
N=12.259、R=2.8591 Å、Rf=0.0022；真实 Mo 箔 `s02=0.88`：
top1=`mo_bcc`，N₁=7.47/R₁=2.7196 Å、N₂=6.20/R₂=3.1355 Å、
Rf=0.0082。

### P5-A T5 运行期溯源（`--s02-from-t5`，2026-11-09）
`run_xafs.py --s02-from-t5 --identify/--fit` 从逐元素
`data/s02_calibration/<元素>/calibration.json` 读取 T5 S₀²，不改
`absorber.py` 静态表。读取规则：

- 只接受 `status="calibrated"` 且 `static_absorber_py_unchanged=True`；
  未通过元素（如 Mo/Rh/Pd/Ir/Ce/U）抛 `XAFSConfigError`，错误含
  “未通过”，不静默回退到静态表或猜测值；
- 显式 `--s02` / `--multi-element-s02` 优先；显式 `--absorber` 优先于
  文件头元素判定，避免 S₀² 与拟合元素错配；
- 溯源串写入 L1/motif 输出与 `<stem>_meta.json.fit.s02_source`：
  `s02_calibration:<相对路径>#sha256:<12位>;standard_sha256:<12位>`；
- 正式 CN 表把 T5 产物来源记为 `s02_calibrated=True`（普通 `--s02`
  手填值仍为 `False`）；
- `calibrate_t5.py --only <元素>` 只重跑所选元素，其余元素从既有
  `<元素>/calibration.json` 合并回 `index.json`/`REPORT.md`，不丢结果。

Cr₂O₃ 真实谱冒烟：`fit.s02=0.6482`，`s02_source` 以
`s02_calibration:` 开头；L1 `N=6.34±0.21`、`R=1.975 Å`、`Rf=0.0216`。
契约 #122 锁定产物完整性与 `--only` 增量合并，#123 锁定加载器与
来源优先级。

### P5-B Mo 双壳 bcc 与 Ce 物理排查（2026-11-09）
**Mo：已通过。** BMM Mo 箔的导数峰 20018.8 eV 落在近边肩峰上，
归一化平台 0.842（低于 0.85）；改用半高点 E₀=20012.86 eV 后平台
0.863。bcc Mo 第一壳 8×2.725 Å 与第二壳 6×3.147 Å 同属 Mo–Mo，
单壳模型把第二壳尾巴折进第一壳振幅（S₀² 撞 1.3 上界）；双壳模型
共享 S₀²/ΔE₀、各壳独立 ΔR/σ²，得到：

| 指标 | 值 |
|------|-----|
| S₀² | 0.8322 ± 0.0334 |
| Rf | 0.0040 |
| 窗口极差 | 0.0017 |
| ΔE₀ | −0.73 eV |
| n_variables / N_idp | 6 / 6.88 |

三个窗口全部 `interior`，门禁 5/5（含新增 `sigma2_not_at_bound`）。
`python calibrate_t5.py --only Mo` 重跑，产物 `data/s02_calibration/Mo/calibration.json`；
运行期需同时指定 E₀ 方法：
`run_xafs.py <Mo 谱> --absorber Mo --e0-method half_height --s02-from-t5 --fit`
（不带 `--e0-method half_height` 时导数 E₀ 会让平台落回 0.842，被严格归一化
门禁拒绝——这是物理门禁，不是缺陷）。

**Ce：仍未通过，不放宽门禁。** 唯一可用的 CeO₂ 荧光谱在导数 E₀ 下
ΔE₀=9.8–10.2 eV，贴 `E0_MAX=10`；`half_height` 把 ΔE₀ 推到
12.9–13.5 eV。自由 σ² 在窄窗内塌到 1e-4 下界（`sigma2_at_bound`），
S₀²=0.46 与这个非物理 σ² 绑定；固定 σ² 时 S₀² 随假设值从 0.54
单调漂到 0.80，无唯一解；加入 Ce–Ce 3.82 Å 第二壳后 Rf=0.127
（超 0.12）、ΔR 撞界。备选 `Ce-L3-CeriumDioxide.xdi` 是透射谱，
边步被背景淹没（归一化平台 2.27），数值上虽可通过宽松拟合，但属于
归一化伪迹，严格门禁必须拒绝。以上均记录在 `Ce/calibration.json`
的 `attempts` 中，静态表仍 `s02=None`。

新增门禁 `sigma2_not_at_bound`：任一窗口 σ² 压到 1e-4 下界即判失败。
这是更严格的物理约束，不放松任何原有阈值；契约 #124 锁定。

### P5-C 固定 COD ID 补齐常见相（2026-11-09）
`tools/cif_crawler.py` 新增 `--p5c` 目标组（tier=`P5C`），按固定 COD ID
补齐 `phase_pack.COMMON_PHASES` 中的 18 个 `missing_common_phases`：

| 元素 | 常见相 | COD ID | 空间群 | 说明 |
|------|--------|--------|--------|------|
| Fe | FeO | 1011169 | 225 | wüstite |
| Fe | FeOOH | 9002158 | 62 | goethite，298 K |
| Pt | PtO2 | 1008935 | 58 | β-PtO2；COD 无常压结构，已标注高压合成 |
| Cr | Cr2O3 | 9014036 | 167 | eskolaite，常压 |
| Mn | MnO | 9006658 | 225 | manganosite，P=0 |
| Mn | MnO2 | 9009081 | 136 | pyrolusite |
| Mn | Mn3O4 | 1514120 | 141 | hausmannite（Mn1/Mn2 两位点） |
| Zn | ZnO | 9004178 | 186 | zincite，293 K |
| Mo | MoO2 | 9009090 | 14 | tugarinovite |
| Mo | MoO3 | 9009669 | 62 | molybdite |
| Ru | RuO2 | 9007541 | 136 | rutile |
| Rh | Rh2O3 | 1544346 | 60 | orthorhombic-II（Shannon & Prewitt 1970） |
| Pd | PdO | 1009031 | 131 | palladinite |
| Ag | Ag2O | 1010486 | 224 | cuprite |
| Ir | IrO2 | 2101854 | 136 | rutile |
| Au | Au2O3 | 2105390 | 43 | gold(III) oxide |
| U | UO2 | 9009049 | 225 | uraninite |
| U | U3O8 | 2310519 | 21 | α-U3O8（中子衍射，U1/U2 两位点） |

结果：phase pack 由 26 条增至 **44 条**；`tools/cif_crawler.py` 固定目标
由 21 个增至 **39 个**；`phase_coverage.json` 的 `missing_common_phases`
全部为空。18 条全部为实验 `crystal_cif`，每条都通过 gemmi 解析、IT 空间群
与化学计量核验。高压/高温/DFT 条目显式排除（如 Cr2O3 的 63.9 GPa 条目、
β-U3O8、高温 U3O8）。Mn3O4 与 U3O8 有多个吸收位点，正式 CN 仍必须显式
`--absorber-site`，不自动 Top-1。

运行期仍只读 `_structure_lib/_phase_pack/`；`data/cif_raw/` 保持 17 个
根索引 CIF，P5-C 下载缓存仅用于构建，导入后清理。契约 #125 锁定
18 个目标、44 条索引与全覆盖。

### P5-prep 标准谱源隔离池（2026-11-09）
T5 未通过的 5 个元素进入**逻辑隔离池**：不移动/删除原始谱，不产生
S₀²，运行期不得用于 `--s02-from-t5`。索引
`data/s02_calibration/source_pool/index.json` 由
`tools/build_source_pool.py` 从标定产物与门禁常量确定性重建：
Rh/Pd/Ir 为短谱（k_max≈5.7 Å⁻¹，需换 k_max≥10–12 Å⁻¹ 的谱）；
Ce 已经是 CeO2 基准，问题是 ΔE₀ 贴 10 eV、σ² 压到 1e-4、无 interior
窗口；U 谱长足够但 α-U 第一壳模型不稳，UO2 切换需先核实 UO2 标准谱。
Ir 的“提高 k_min”建议只保留为换谱后的诊断项——当前 k_max 下
N_idp 不足以支撑 4 个参数，不能靠调窗绕过门禁。全部门禁阈值冻结，
`absorber.py` 静态表不变。契约 #126 锁定 5 元素分类、`s02=None`
与门禁常量；详见 `data/s02_calibration/source_pool/README.md`。

### B19-1 原位序列 EXAFS 联合精修
`xafs_fit.run_operando_exafs_series(results, absorber_symbol, motif=None,
catalytic_model=None, ..., share_sigma2_across_series=False)` 对 ≥3 个
时间点共享单一全局 ΔE₀（`share_e0=True`），逐点放开 N/ΔR；
`share_sigma2_across_series=True` 时同一壳层跨点共享 σ²。输出
`trajectories[label] = {reff,N,N_err,R,R_err,sigma2,sigma2_err}`、
`global_r_factor`、`shared_delta_e0` 与 `r_factors`。
`xafs_diagnostics.write_operando_exafs_diagnostics` 落盘
`<stem>_exafs_evolution.csv`（每点 N/R/σ²）与 `_exafs_evolution.png`
（N(t)±σ_N、R(t)±σ_R）。合成 R=[2.02,1.99,1.96,1.93] Å、
N=[4.0,4.6,5.3,6.0] 精确回收，`global_r_factor≈6e-21`。

### B20-1 电催化氧化物/羟基氧化物库与真实标样盲测
`_structure_lib/build_seed.py` 落盘 7 个 Motif：`coo_rocksalt`
（Co–O 6×2.1306、Co–Co 12×3.0131）、`co3o4_tet`（Co–O 4×1.9323）、
`co3o4_oct`（Co–O 6×1.9217、Co–Co 6×2.8581）、`co3o4_mixed`
（`kind="mixture"`，tet:oct=1:2）、`ni_oh2_brucite`（Ni–O 6×2.0697、
Ni–Ni 6×3.126）、`cuo_tenorite`（Cu–O 4×1.9558）、`pto2_beta`
（Pt–O 6×2.0495、Pt–Pt 6×3.100）。真实 BMM 标样盲测：
CoO top1=`coo_rocksalt`（score 0.033，confident，S7 Rf=0.0076）；
Co₃O₄ top1=`co3o4_mixed`（score 0.0249，confident，无未解释峰）；
Ni(OH)₂ top1=`ni_oh2_brucite`（score 0.1159，confident，S7 Rf=0.0387，
Ni–O 6.67/2.0698 Å、Ni–Ni 7.64/3.1255 Å）；CuO top1=`cuo_tenorite`
（score 0.0747；`uncovered_peaks=[2.056]` 故 status=partial，按 #74
契约只要求 top1 与 score<0.10）。新增 Motif 第一壳 1.92–2.13 Å
与金属箔 ≥2.49 Å 错开，Co/Ni/Cu/NiO 压力测试 top1/score 零退化。

## 科学边界（必须读）

- **S₀² 是本管线归一化约定下的相对值**。`absorber.py` 已有 Fe 0.803、
  Co 0.7308、Ni 0.7795、Cu 0.8670、Pt 0.8075（真实箔标定，协议见
  `calibrate_s02.py` 模块注释：N 固定晶体学、生产等效拟合、3 窗口鲁棒性
  门禁、kw 敏感性只记录不门禁）。S₀² 跨管线不可移植；Fe 的 0.803 是
  氧化物价态锚定的跨材料折中值，BMM 箔单壳重标定 0.6566 只适用于
  金属 Fe（氧化物 CN 用 0.6566 会高约 22%）；金属链显式
  `--chemical-state metal`。
- **形式误差 ≠ 绝对精度**。`N_err` 只来自最小二乘协方差；S₀² 传递、
  代表路径/壳层合并、N–σ²–ΔE₀ 相关和归一化偏差会带来约 ±20–40%
  的系统误差。输出里的 `cn_reliable` 和 `cn_systematic_budget` 必须一起看。
- **代表路径近似**。同一壳层的多条 SS 路径合并成一条代表路径 + 总简并度；
  换一条代表路径会改变 N 与 Rf，`cn_reliable=False` 时不要拿 N 做结论。
  （Cu 箔即例：第一壳 N=11.6 可靠，第二壳形式误差 1.62 → 整体拒标。）
- **简并簇不可辨识**（针铁矿实测）。3.0–3.44 Å 的 Fe/O 简并簇导致 S7
  Rf≈0.17、`cn_reliable=False`——这是正确判定而非 bug：换 Fe@3.31 能把
  Rf 压到 0.11，但 N 全撞上界 6.0（简并吞噬，失去定量意义）；中距离壳
  σ² 被压到 0.0001–0.0007 Å²（≤5e-4 触发 `flags.sigma2_at_bound`）
  是形状失配证据。第一壳仍可靠（N=5.89±0.51 vs 晶体学 6）。详见
  `references/structure-identification.md` §4。
- **覆盖清单**。`partial` 表示最优候选没有解释全部显著峰（多位点残缺/
  混合相/库缺条目），此时默认只做 L1。**mixture（多位点联合）条目例外**：
  χ 是各位点 χ 的加权求和，复合谱受干涉位移/压制，覆盖判定按
  "复合谱 ∪ 各分量分别解释"（磁铁矿 tet 分量解释 1.35/3.712 Å，
  oct 分量解释 2.608 Å）；mixture 拟合时各位点路径简并度按权重缩放。
- **mixture 相分数**（`--fit-weights`）。识别用固定晶体学权重；定量
  时改为固定 `N_crystal`、放开 softmax 组分权重（和=1），输出
  `mixture_weights` 与 `weights_reliable`。合成磁铁矿 tet/oct 第一壳
  0.25/0.75 可回收至 ≤0.03（回归测试 13）；真实 `Fe-Magnetite.xdi`
  默认窗 Rf=0.097、w_tet=0.17±0.07，因 Rf>0.08 被正确判为
  `weights_reliable=False`——相分数是模型组分 χ 的线性组合权重，
  不是独立物相定量，N 列也不是配位数。
- **弱结晶/非晶**（水铁矿等）没有可靠体相 motif，识别会 `no_match`，
  L1 第一壳是唯一诚实的定量路径。
- **荧光自吸收校正（S4.5，可选）**。默认只做风险筛查；要真正校正必须
  同时提供以下五项，缺一即报错，管线不做假设：

  | 参数 | CLI | 说明 |
  |------|-----|------|
  | 样品总化学式 | `--self-absorption-formula` | 含基体/稀释剂（如 `Fe3O4`、`((NH4)2SO4)0.1(H2O)54.8`），不能只写活性相 |
  | 厚度 | `--self-absorption-thickness` | μm；厚样可取大于信息深度的值 |
  | 有效密度 | `--self-absorption-density` | g/cm³；粉末压片用几何密度，不是晶体密度 |
  | 入射/出射角 | `--self-absorption-angles 49.4,41` | 相对样品表面，度（90° = 垂直入射/出射；常见 45,45） |
  | 荧光线 | `--self-absorption-line` | 默认 K→Kα、L3→Lα；探测器收 Kβ 等需显式覆盖 |

  数学：Booth & Bridges 2005 Eq.4，在 k 空间反解
  `χ_exp → χ`；只改振幅、不改相位。μ_T/μ_a 用 xraydb 光滑截面，
  μ_a 扣边前基线、能量锚定表值边 `E_edge+0.1`（与 Athena/Demeter
  的 Booth 实现一致，避免样品 E0 被自吸收拉低后 k=0 落到边以下）。
  薄/厚样由厚度自动落入对应极限，d→0 时严格回到 χ_exp。
  合成自洽回收：max|Δχ|=0.0015、rms=0.00050；**全链合成验证**
  （真实磁铁矿 μ → 前向 Eq.2 生成荧光谱 → 归一化/autobk → Booth）：
  χ 振幅比从 0.723（未校正）回到 1.054，rms=0.020。厚样与薄样极限
  均有回归测试（第 15 项）。
  边界：假设均匀平表面；颗粒/孔洞/不均匀会让有效 d、ρ 失真；化学式
  偏差直接进入校正因子。**自吸收会压低并移动 E0**（合成厚样实测
  7111.3 vs 表值 7112 eV），校正不修 E0——拟合前必须用参考箔或已知
  边能量校准 E0。校正后 `self_absorption.reliable=False` 时
  CN 与相分数仍被拒标。**它不能替代透射谱**——能提供同样品透射谱时，
  应优先用透射做独立对照。
- **WT 峰位定位**。WT 的 R 窗以 FT 已选峰为中心（±0.5 Å），不是固定
  `[1.0, 2.6]`。多壳体系在宽窗内会同时出现两条分支（Fe3O4 实测：
  Fe–O R≈1.42/k≈4.5 与 Fe–Fe R≈2.72/k≈5.8 只差 3%），任何轻微振幅
  变化都会让 argmax 跳壳；锚定后 S4.5 校正前后 WT 峰位移 ≤0.02 Å
  （第 16 项回归）。WT 回答的是"这个 FT 峰主要来自哪个 k"，不是独立
  重新选壳。
- **严格门禁**。`strict=True` 下边后平台不在 0.85–1.15 会中断；金属箔
  合并扫描/厚样可用 `--no-strict`（标定门禁在拟合层），但 FT 峰位/幅度
  不可直接引用。
- **L 边**。`artemis.build_feff_input(..., edge="L3")` 已生成正确 HOLE 卡；
  Pt L3 全流程已验收（`Pt/pt_fcc` 种子库 + `Pt-L3-PlatinumFoil.xdi` 标定
  S₀²=0.8075，第一壳 N=11.69±0.22、R=2.764 Å、Rf=0.0038）。
  Pt 库目前只有金属箔一条，分子配合物会正确 `no_match` → L1（用
  `--l1-ligand Cl/N` 指定配体；默认 O 只适用于氧化物）。
