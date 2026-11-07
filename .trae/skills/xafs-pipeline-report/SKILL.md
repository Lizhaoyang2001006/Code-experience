---
name: xafs-pipeline-report
description: Run the XAFS processing pipeline (athena/artemis/larch_wt modules) and produce the five-panel XAS figure plus EXAFS CN fitting table. Use when the user asks to process XAFS/EXAFS data, fit coordination numbers from a CIF, or generate XAS report figures and tables. Do not use for general matplotlib plotting or non-XAFS data.
---

# XAFS Pipeline & Report

> **最高优先级工程规范（必读）**：涉及补齐 Athena / Artemis / Larch-WT 全功能（Phase I：B1–B4 已完成 36 项回归基线；Phase II：B5–B8 进阶补齐至 45 项；Phase III：B9–B12 高阶能力已完成 60 项回归基线；Phase IV：B13–B16 面向催化化学全面分析已验收 68 项；**Phase V：B17–B20 催化剂精细结构深度解析与实战完备化已验收 74 项回归基线**）的小模型执行与旗舰模型批改验收任务，一律以仓库根目录 [`C:\Users\lizha\Desktop\Code-experience\EXECUTION_MASTER_PLAN.md`](../../../EXECUTION_MASTER_PLAN.md) 为唯一最高级准则。

目标：任何模型按本 skill 操作，都能从原始 XAFS 数据稳定产出
①五联图（μ(E)/XANES/k²χ(k)/|χ(R)|/WT）②EXAFS CN 拟合参数表（Table S 风格）。
所有函数签名、陷阱、数值基准均已在真实数据上验证（2026-10）。

## 0. 环境（先读，错了全盘崩）

- **解释器**：必须用 WorkBuddy Python（PATH 上的 `python` 是 Microsoft Store
  stub，`Get-Command python` 指向 WindowsApps，无任何依赖库）：
  `C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe`
  （目录名 3.13.12 是历史遗留，实际 `sys.version` = **3.13.14**）
- **依赖库**（env_check.py 实测版本，缺一不可跑管线）：
  numpy 2.5.4 / scipy 1.18.1 / matplotlib 3.11.2 / gemmi 0.7.5 /
  xraylarch 2026.3.1
- **环境自检**（agent 接手任何任务前必跑一次）：
  1. `Get-Command python` → 确认 PATH python 是 Store stub，不得裸用
  2. `& "C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe" "C:\Users\lizha\Desktop\Code-experience\env_check.py"`
     → 逐库打印版本，全部可导入 exit 0；缺库 exit 1 并列出缺哪些
- **源码根**：`C:\Users\lizha\Desktop\Code-experience\`（xafs_core.py / xafs_io.py /
  athena.py / artemis.py / larch_wt.py / xafs_pipeline.py /
  xafs_diagnostics.py / xafs_fit.py / xafs_report.py / run_xafs.py /
  absorber.py / structure_lib.py / structure_id.py / _structure_lib/ /
  lepidocrocite.cif）
- **测试数据**：`C:\Users\lizha\Desktop\Code-experience\data\练手数据\`
  （01_Fe_foil\fe.060/150/300 三列 XDAC；02_Fe_standards\Fe_lepidocrocite.000/100/200 五列 XDAC）
- **依赖方向**（两链一契约，全部单向无环）：
  core ← io ← athena ← {artemis, larch_wt} ← xafs_pipeline ← {run_xafs, xafs_report}；
  xafs_pipeline ← xafs_diagnostics；
  artemis ← structure_lib ← structure_id ← xafs_fit ← {run_xafs, xafs_report}；
  absorber 独立（仅依赖 core；pipeline 顶层引用其配置查表）。
  归约链（S1–S6）与拟合链（S6.5/S7）以 PipelineResult 为唯一契约；
  S6.5 只经 xafs_fit.run 触发，入口不得直调 structure_id。
  每个模块可 `python xxx.py` 独立自测。
- **PowerShell 铁律**：禁止 `python -c "..."` 内联多行/含花括号代码（转义必崩）。
  一律写临时脚本 `_tmp.py` 执行，跑完删除。
- 模块自测命令（改动后回归必跑；先过 env_check.py）：
  `& "C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe" "C:\Users\lizha\Desktop\Code-experience\<module>.py"`
- **元素普适性压力测试**（动过 absorber / structure_* / athena 窗口后必跑；
  新元素接入后按同模板各建一条）：
  `& "C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe" "C:\Users\lizha\Desktop\Code-experience\stress_test_co.py"`
  `& "C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe" "C:\Users\lizha\Desktop\Code-experience\stress_test_ni.py"`
  `& "C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe" "C:\Users\lizha\Desktop\Code-experience\stress_test_cu.py"`
  `& "C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe" "C:\Users\lizha\Desktop\Code-experience\stress_test_nio.py"`
  → 五阶段：合成真值数据 → S1–S6 归约 → 三护栏逐条触发 → 内存接入
    → CN 真值回收（金属箔 Rf 均 <0.005；NiO 双壳层
    N(Ni-O)=7.64/N(Ni-Ni)=11.92，Rf=0.006）。顺产 `_structure_lib/{Co,Ni,Cu}`
    金属参考相 + `Ni/nio_rocksalt` 氧化物种子与 `xafs_out/_*_stress/` 留档。
    stress_test_nio 兼测氧化物 vs 金属鉴别力（nio score 0.003 vs ni_fcc 0.678）。

## 0.5 用户输入清单（agent 先问后跑：只问化学，不问物理）

服务对象是催化方向的科研人员：用户只回答**化学问题**；一切物理处理
参数（Rbkg、k 窗口、kweight、E0 定义、σ…）由管线自动决定，并在
indicators.csv 中**完整披露**，不向用户征询物理参数。

**开跑前必答（化学四问，缺一不跑）**：

1. **元素与边**：测的哪个元素、哪条边（Fe K？Pt L₃？）。S1–S6 算法
   元素无关，但 S2 定边窗口与 S6.5/S7（选库、S₀²）必须知道元素；
   文件头无 e0_nominal 又不给元素时 S2 直接报错而不是猜。
   CLI：`--absorber Fe`（无默认）；absorber.py 未标定（s02=None）
   的元素拒绝 CN 定量（见"新元素接入"）。
2. **材料是什么**：有没有 XRD 谱 / 精修 CIF / 文献结构？已知晶相
   → 库匹配、strict 门禁可降级；未知/弱结晶 → 预期走 L1 只报第一壳。
   XRD 物相结论直接交给 agent：缩小假设库、裁决并列候选。
3. **化学组成**：化学式 / 负载量 / 载体（ICP、XRF 或投料比）——
   决定第一壳散射体可能是 O/N/C/S/Cl 中哪些、哪些物相可被排除。
4. **合成与后处理**：气氛（N₂/H₂/空气）、温度、酸洗与否——直接改变
   候选物相概率（如 N₂ 碳化 + 未酸洗 → 金属/碳化物颗粒共存概率大）。

多扫描关系 agent 从数据自动判断（独立参考边→对齐；同一样品→直接合并），
仅歧义时补问一句"这几个 scan 是同一样品吗"。

**文献对比原则**：用户要"和某文献直接对比"时，让用户给出**文献本身**，
agent 从中提取处理参数（Rbkg/k 窗口…），不让用户翻译物理参数。

**触发式追问**（数据说话时才问，且用化学语言）：
- S1 嗅探失败/通道歧义 → "哪几列？透射还是荧光？"
  （荧光默认只做自吸收风险筛查：Δμ≥1.0 拒绝对 CN 下结论；
  要真正校正须主动向用户索要**样品总化学式、厚度 μm、有效密度
  g/cm³、入射/出射角**，缺一不可，不能代猜）
- S3 平台越界 → "样品压得比较厚，还是稀释在低 Z 基体里？"
  （agent 据此决定门禁降级或窗口调整，不问 post_window）
- S6.5 并列候选 / 未解释峰 / L1 降级 → 拿候选列表问"合成上更可能是哪个？"
- S7 前 → 确认结构来源（用户 CIF / 识别 top1 / L1）与报告口径

**新元素接入三步**（用户提供化学信息，agent 不得代猜）：
1. 元素与边 → absorber.py 加条目（s02=None 起步；L 边配 next_edge）
2. 金属箔标样 → `python calibrate_s02.py <箔数据> --absorber X --no-strict [--apply]`
   （生产等效协议 + 3 窗口鲁棒性门禁，详见 references/pipeline-reference.md。
   当前已标定：Fe=0.803 / Co=0.7308 / Ni=0.7795 / Cu=0.8670 / Pt=0.8075（L3）；
   **Fe 箔数据在现行协议下被门禁拒绝**（厚样谐波嫌疑），Fe 值由
   氧化物价态验收结果锚定保留）
3. 该元素结构假设（CIF / 文献坐标 / 构造团簇）
   → 见 [references/structure-identification.md](references/structure-identification.md)

接入后用该元素标样重跑并记录数值基准（Fe 基准见 §1.4 / §2.7，
是回归锚点，不代表管线只支持 Fe）。

## 1. 执行层管线（S1–S6）

整条链已封装为一次调用（run_xafs 与 xafs_report 共用）：

```python
import sys; sys.path.insert(0, r"C:\Users\lizha\Desktop\Code-experience")
import xafs_pipeline

res = xafs_pipeline.process_sample(files, align=True, absorber="Fe", strict=False)
# res: PipelineResult — scans/merged/e0/norm/rbkg/chi_k/ft/wt/r1/r_wt/
#      absorber_symbol/rbkg_table（S6.5 见 xafs_fit）
# strict=True（默认）：边后平台 ∉[0.85,1.15] 直接抛错（未知样品用）；
# strict=False：降级为警告（已知标样如厚 Fe 箔平台 0.79 用）。
```

逐阶段诊断图（每阶段一张 PNG，只读 PipelineResult 不重算）：

```python
import xafs_diagnostics
xafs_diagnostics.write_diagnostics(res, outdir / "diagnostics")
# s1_raw_scans / s2_e0 / s3_normalization / s4_rbkg / s5_chik_chir / s6_wavelet
```

需要逐步调用时的底层签名与陷阱见
[references/pipeline-reference.md](references/pipeline-reference.md)。

S6.5 结构识别（拟合链；指纹匹配排序 + 覆盖清单，不做 CN 定量）：

```python
import xafs_fit
sid = xafs_fit.run(res)["structure"]   # 唯一触发点；回填 res.structure
# sid: ranked（按 score 升序）/ indistinguishable（并列候选）/
#      uncovered_peaks（未解释实验峰 → 应扩库）
```

motif 库（`_structure_lib/<元素>/`）、注册/入库/扩展流程与设计哲学见
[references/structure-identification.md](references/structure-identification.md)。

**执行层不变量**（违反即结果不可信）：
1. `ft_to_r_space` 的第 3 参 `r` 被忽略（R 网格由 xftf 决定）；**绝不自写 FFT**
   （相位因子须为 e^{−2ikR}，自研版实测差 0.05 Å）。
2. 小波 `sigma >= 1.2`（更小则相邻壳层粘连）；`kmax` 已在入口钳制。
3. XDAC 有 3 列（无参考）与 5 列（带 Ir 参考）两种，读取均已兼容；
   2 列纯文本的 μ 中位数须在 0.001–50，否则按 I0 冒充 μ 拦截。
4. 数值基准（Fe_lepidocrocite 3 扫描合并）：
   E0=7128.60 / edge_step≈0.4695 / Rbkg=1.20 / FT 第一峰 1.473 Å。
   偏离即查上游步骤。

## 2. CN 拟合层（CIF → FEFF → R 空间拟合）

```python
sites  = artemis.read_cif_sites(cif_path, absorber="Fe", rmax=5.0)
inp    = artemis.build_feff_input(sites, absorber="Fe", rmax=5.0, title="...")
files  = artemis.run_feff6l(inp, workdir)        # 返回绝对路径列表
rows   = artemis.describe_paths(files, 20)       # rows[i]["file"] 是相对名！
from absorber import get_absorber
fit    = artemis.fit_r_space(chi_k, ft, datfiles, rmin, rmax, kmin=3.0,
                             kmax=12.0, s02=get_absorber("Fe").s02,
                             vary_n=True)   # s02 无默认：未标定元素直接拒绝
```

**拟合层铁律**：
1. **S₀² 固定、N 放开**（`vary_n=True` 或索引列表）——两者在振幅上只以
   乘积出现，同时放开必简并。S₀² 取 absorber.py 标定表（
   `get_absorber(symbol).s02`；Fe=0.803 / Co=0.7308 / Ni=0.7795 /
   Cu=0.8670 / Pt=0.8075），函数无默认——未标定元素拒绝定量。
2. **不可分辨路径合并**：ΔR < π/(2Δk)（k=3–12 时 ≈0.26 Å）的 SS 路径并为
   一条有效壳层，N 语义变为该壳层总配位数（如 Fe-O 1.98×2+2.07×4 →
   单壳层 N≈6）。
3. **MS 路径慎加**：与主壳层简并时 N 误差从 ±0.6 爆到 ±13；加前对比
   变量数与 N_idp（`fit.n_variables < fit.n_independent`）。
4. **kmax 默认 12**：数据再好也别贸然到 14——模型缺路径时 kmax 越大
   R-factor 越差（实测 0.074→0.092）。
5. **ΔE₀ 判据用 |ΔE₀| < 15 eV**（参考 Table S4 脚注同款；该表自己的
   Fe-foil 行就是 6.5 eV）。ΔE₀ 是实验 E₀ 标记与 FEFF6 阈值的定标差，
   不是待最小化的物理量：E₀ 定义平移会 1:1 传给 ΔE₀ 而 R 纹丝不动
   （实测 derivative 7128.6→ΔE₀=−5.2；half_height 7125.5→−2.2；
   两者的 R 都=1.994）。γ-FeOOH 在 derivative 下 −5.2 eV 正常
   （Fe³⁺ + FEFF6 非自洽势；同管线 Fe⁰ 箔仅 +0.9）。**绝不为压
   ΔE₀ 而硬设 E₀**——实测 E₀=7122 会让 N(Fe-O) 6.5→4.9、
   Rf 0.074→0.083；ΔE₀≡0 更糟：Rf→0.112、R 虚涨到 2.03。
   需要小 ΔE₀ 就用 half_height 定义（Athena 内置，合法）。
6. **CIF 可信度自检**：展开后配位多面体必须闭合（Fe-O 得 6 键）且键长
   对照文献；从 PDF/表格抄的坐标常串行，不闭合即弃用并查别的来源。
7. 验收（γ-FeOOH, k=3–12, S₀²=0.803）：
   Fe-O N=6.5±0.4 R=1.994；Fe-Fe(1) N=7.3±0.6 R=3.052（合并 6+2 边共享）；
   Fe-Fe(2) N=1.0±0.6 R=3.979；R-factor 0.074（简化 3 壳层模型）。
   Fe 箔标定（2026-11 并列裁决后复核值）：Rbkg=0.80，N=6.8/7.6,
   R=2.446/2.821, R-factor 0.0036, ΔE0=+0.64 eV。
   （注：optimize_rbkg 的 5% 并列裁决将箔 Rbkg 从 0.90 改为 0.80，
   拟合 R/Rf 不变、N 在误差棒内移动，S₀²=0.803 标定仍成立。）

## 3. 绘图与导出层（xafs_report.py）

报告层已模块化：`process_scans`（S1–S6 一层封装）→ `plot_five_panel` /
`plot_fit_table`（图）→ `export_panel_csvs` / `export_fit_table` /
`export_indicators`（数据）。**优先复用该模块，不要从零写图。**
逐面板规格见 [references/figure-spec.md](references/figure-spec.md)。

一次运行产出（`xafs_out\<sample>_report\`，共 11 个文件 + 诊断图 6 张）：

| 文件 | 内容 |
|------|------|
| `<stem>_xas.png` | 五联图；(a) 标注**吸收边 + E₀**（如 Fe K-edge）；(b) 标注**白线**；(d)(e) 标注**壳层归属**（Fe–O / Fe–Fe，未相位修正的 R） |
| `<stem>_fit_table.png` | CN 拟合参数表（Table S 风格） |
| `<stem>_fit_table.csv` / `.txt` | 同表：utf-8-sig CSV + 制表符 TSV（粘贴 Word/Excel 自动成表） |
| `<stem>_indicators.csv` | **关键物理指标汇总**（元素/边/E₀/边跳变/白线/Rbkg/k 窗口/相位修正状态/WT 壳层位置/S₀²/ΔE₀/R-factor/各壳层 N,R,σ²）——读图猜测的终结者 |
| `<stem>_mu.csv` `_xanes.csv` `_chik.csv` `_chir.csv` | (a)–(d) 面板逐点数据 |
| `<stem>_wavelet.csv` | (e) W(R,k) 矩阵（行=R，列=k） |
| `structure_id.json` | S6.5 结构识别：候选排序 score/ΔE0、并列候选、未解释峰清单 |
| `diagnostics/` | S1–S6 各阶段诊断图 6 张（只读 PipelineResult，不重算） |

**导出三铁律**：①每张图必有对应 CSV；②WT 面板必标壳层线+峰标记
（`larch_wt.locate_shell` 定峰）；③indicators.csv 必含"元素+边"、
"not phase-corrected"两条防误读声明。

## 4. 交付前检查清单

- [ ] 用户化学四问已确认（元素与边/材料与晶相/化学组成/合成处理）；
      非 Fe 元素已完成新元素接入——缺答禁止开跑
- [ ] 全部用 WorkBuddy Python 跑通，无 traceback
- [ ] 数值基准 §1.4 与验收 §2.7 在容差内
- [ ] 五联图逐面板对照 figure-spec.md：E₀ 线、白线标注、(d)(e) 壳层标注
- [ ] 每张图有对应 CSV；(e) 的 wavelet.csv 行=R 列=k
- [ ] indicators.csv 已生成且含 absorbing_edge / phase_correction 两行
- [ ] 表格脚注含 S₀²、结构来源、拟合窗口、路径合并说明
- [ ] 11 个输出文件落盘（图 2 + 表 3 + 面板数据 5 + structure_id.json）
- [ ] S6.5 已跑：top1 与样品一致；indistinguishable/uncovered_peaks 已人工过目
- [ ] 元素普适性：动过元素路径/窗口/Rbkg 寻优时 stress_test_{co,ni,cu,nio} 全绿
- [ ] 临时 `_*.py` 脚本已删除
