# XAFS 当前出图清单

本文件说明当前代码在各入口**默认**和**显式开关**下会生成哪些图，以及每张图的
文件名、面板组成和触发条件。以代码为准；本文件最后核对于
`2918b4f`（2026-11-10）。

`xafs_out/` 下的历史目录可能包含旧脚本或手工脚本生成的图（例如
`*_compare.png`）。它们**不是**当前入口的产物，不能用来推断运行行为。

## 0. 快速总表

| 入口 | 默认出图 | 显式开关 | 主要图文件 |
|---|---|---|---|
| `paper_figures.py` | 3 面板论文主图 | `--fit-json` 增加 3 个拟合面板和 CN 表 | `<stem>.png`、`<stem>_ft.csv`；有拟合时再加 `<stem>_cn_table.{png,csv,md}` |
| `xafs_report.py` | 5 联报告图；默认跑 S7，有拟合行时加拟合表 | `--diagnostics` 增加 S1–S8；`--figure-pack` 增加图包 | `<stem>_xas.png`、`<stem>_fit_table.png` |
| `run_xafs.py` | 不出图，只落数据 | `--diagnostics`、`--phase-correct`、`--chi2-contour`、`--figure-pack` | `diagnostics/*.png`、`figures/*` |
| `xafs_figure_pack.py` | 4 张总图 + 2 张表图 | `--standalone-panels` 再出逐面板 PNG | `figures/Figure*.png`、`figures/panels/*.png` |

默认文件名中的 `stem`：

- `paper_figures.py`：`--stem`，未给时为 `<样品文件 stem>_paper`。
- `xafs_report.py` / `run_xafs.py`：输入文件 stem，多文件时为合并后的 stem。
- `xafs_figure_pack.py`：从输出目录的 `*_meta.json` 自动识别，或 `--stem` 指定。

---

## 1. `paper_figures.py`：通用论文主图

这是跨体系复用的论文版式，只画图，不重新拟合；S₀²、E0、k/R 窗口、kweight
都从输入或 fit JSON 原样读取。

### 1.1 无 `--fit-json`：3 面板

输出：

- `<stem>.png`（默认 `<sample>_paper.png`）
- `<stem>_ft.csv`（|FT| 面板实际绘制数据）

| 面板 | 内容 |
|---|---|
| (a) | XANES，横轴相对 E0；标注吸收元素与边 |
| (b) | `k^n χ(k)`，`n = --kweight`（默认 3） |
| (c) | `\|FT(k^n χ(k))\|`，带 `--guide` 标线、Rbkg 以下阴影 |

### 1.2 有 `--fit-json`：6 面板 + CN 表

在上面的 (a)–(c) 之外增加：

| 面板 | 内容 |
|---|---|
| (d) | R 空间拟合：数据 vs 模型、拟合窗口、Rf/ΔE₀ |
| (e) | k 空间拟合：数据 vs 模型，kweight 来自 fit JSON |
| (f) | CN 汇总表（逐壳层 N、R、σ²、状态） |

同时输出：

- `<stem>_cn_table.png`：独立 CN 表图
- `<stem>_cn_table.csv`：CN 表数据，含 `disclosure_status` / `disclosure_issues`
- `<stem>_cn_table.md`：CN 表 Markdown，含结构/CIF provenance 与披露警告

触发与限制：

- fit JSON 必须含足够的 `path_files`，否则直接报错（无法重建 d/e 面板）。
- (b)(c) 用命令行 `--kweight`；(d)(e) 用 fit JSON 的 `kweights`，可用
  `--fit-kweight` 只改显示。不同 kweight 的 |FT| 峰高不可直接比较。
- CN 披露检查默认只打印警告；`--strict-disclosure` 会升级为错误。
- 常用命令：

```powershell
python paper_figures.py --sample sample_mu.dat --absorber Fe --edge K `
  --fit-json sample_fit.json --out paper_out
```

---

## 2. `xafs_report.py`：五联报告图 + 拟合表

### 2.1 默认产物：五联图

输出 `<stem>_xas.png`，300 dpi：

| 面板 | 内容 |
|---|---|
| (a) | μ(E) 多扫描 + E0 标注 |
| (b) | 归一化 XANES + 白线 |
| (c) | `k²χ(k)` |
| (d) | `\|χ(R)\|`，未做相位修正，标出壳层归属 |
| (e) | 小波 `\|W(R,k)\|`，标出壳层定峰位置 |

注意：(c) 是固定的 `k²χ(k)`；与 `paper_figures.py` 的 `--kweight` 不是同一套
显示参数。

### 2.2 有拟合结果时：拟合表

运行 S7 并得到拟合行时，额外输出：

- `<stem>_fit_table.png`：Table S 风格拟合参数表
- `<stem>_fit_table.csv`：同表数据

`--no-fit` 或没有得到拟合行时，跳过并打印“未运行 S7 或无拟合结果”。

### 2.3 `--diagnostics`：S1–S8 阶段诊断

输出到 `<outdir>/diagnostics/`，共 **9 张 PNG**：

| 文件 | 内容 |
|---|---|
| `s1_raw_scans.png` | 原始扫描一致性、通道角色 |
| `s2_e0.png` | 导数峰 / 参考边判定 E0 |
| `s3_normalization.png` | 边前/边后基线外推与归一化 |
| `s4_rbkg.png` | Rbkg 寻优曲线 |
| `s5_chik_chir.png` | χ(k) 与 \|χ(R)\| 的窗口关系 |
| `s5b_backft.png` | 反 FT：实验 vs 模型 q 空间 |
| `s6_wavelet.png` | k–R 二维小波响应 |
| `s7_ft_wt_joint.png` | FT 峰与小波峰的联合定位 |
| `s8_wt_sweep.png` | 多 σ / kweight 的小波对比 |

配套数据：`diagnostics/wt_sweep.csv`（9 行，σ × kweight 的峰位、幅度和
相对 FT 偏移）。

### 2.4 其他显式图

| 开关 | 输出 | 条件 |
|---|---|---|
| `--phase-correct` | `diagnostics/<stem>_phase_corrected.png` + `<stem>_ft_phase_corrected.csv` | 需要第一壳 FEFF 路径 |
| `--chi2-contour` | `diagnostics/<stem>_chi2_contour.png` + `<stem>_chi2_contour.csv` | 需要 S7 拟合/路径；缺条件时跳过 |
| `--wt-sweep` | `<stem>_wt_sweep.png` + `<stem>_wt_sweep.csv` | 锚定 FT 第一峰 |
| `--wt-sweep` + S7 模型 | `<stem>_wt_profiles.png` + `<stem>_wt_profiles.csv` | 拟合对象含 `model_k` |
| `--figure-pack` | `figures/` 图包（见第 4 节） | 图包失败只警告，不影响已落盘结果 |

### 2.5 HTML 报告与交互页

默认还会生成：

- `<stem>_report.html`：单文件 HTML，图片 base64 内联；`--no-html` 关闭
- `<stem>_window_tuner.html`：交互式 k 窗口调参页；`--no-window-tuner` 关闭

它们不是静态 PNG，但属于当前入口的默认可视化产物。

---

## 3. `run_xafs.py`：数据优先，图按开关生成

默认运行**不生成任何 PNG**，只写 CSV / JSON / 数据文件和进度日志。

| 开关 | 输出 |
|---|---|
| `--diagnostics` | `<outdir>/diagnostics/` 下与第 2.3 节相同的 9 张 PNG + `wt_sweep.csv` |
| `--phase-correct` | `diagnostics/<stem>_phase_corrected.png` + CSV |
| `--chi2-contour` | `diagnostics/<stem>_chi2_contour.png` + CSV |
| `--figure-pack` | `<outdir>/figures/` 图包（见第 4 节） |
| `--operando` | `<outdir>/operando/` 下的演化图与 CSV（见下） |

`--operando` 的图文件名由 `write_operando_diagnostics(..., stem="operando")`
生成，实际为：

- `operando_operando_summary.png`：四宫格演化图
- `operando_operando_concentrations.csv`：MCR 组分分数、LCF 权重、E0/白线/FT 描述符
- `operando_operando_components.csv`：MCR 纯组分谱

该模式只做 S1–S6 + PCA/MCR/LCF，不执行 S6.5/S7，也不支持 `--figure-pack`。

---

## 4. `xafs_figure_pack.py`：出版级图包

只读取已有输出，不重新拟合、不改物理判定。默认输出到
`<outdir>/figures/`。

### 4.1 四张总图

| 文件 | 面板 | 内容 |
|---|---|---|
| `Figure1_master.png` | (a)–(h) | XANES、χ(k)、FT、小波、拟合壳层、结构识别、正式 CN 门禁、结构卡片 |
| `Figure2_confidence.png` | (a)–(f) | 门禁状态板、决策轨迹、参数来源、门禁表、决策表、关键参数表 |
| `Figure3_diagnostics.png` | (a)–(h) | 参数相关矩阵、路径 N–R 图、CN 误差预算、Rbkg 寻优、归一化分解、拟合窗口表、稳定性检查、扫描质量 |
| `Figure4_features.png` | (a)–(f) | dμ/dE、q 空间、FT 峰识别、CIF 配位壳层、特征总表、输入文件清单 |

`--formats png,pdf` 时，每张总图同时生成同名 `.pdf`。

### 4.2 两张表图

| 文件 | 内容 | 条件 |
|---|---|---|
| `Figure2g_cn_table.png` | Table S1，正式 CN 表 | 存在 `cn_table.csv` |
| `Figure2h_fit_table.png` | Table S2，EXAFS 拟合表 | 存在拟合行 |

缺数据时面板标记为 `skipped`，不生成对应 PNG，也不会伪造数据。

### 4.3 逐面板数据与清单

- `figures/panels/<PanelID>.csv`：每个面板实际绘制的数据，例如
  `Figure1a_xanes.csv`、`Figure3a_fitcorr.csv`。
- `figures/data/`：跨面板汇总，通常包括 `summary.json`、
  `panel_index.csv`、`all_fit_paths.csv`、`fit_variants.csv`、
  `all_parameters.csv`、`gate_checks.csv`、`decisions.csv`、`README.md`。
- `figures/figure_manifest.json`：面板 ↔ 图片 ↔ CSV ↔ 源文件 ↔ 图注映射。
- `figures/figure_pack.md`：人可读索引。

### 4.4 独立面板 PNG

默认**不**逐面板渲染 PNG，避免根目录被淹没。需要时：

```powershell
python xafs_figure_pack.py xafs_out/RuO2 --standalone-panels
```

输出到 `figures/panels/<PanelID>.png`，例如 `Figure1a_xanes.png`。
总图仍然保留。

---

## 5. 当前 Pt 案例的实际产物（仅示例）

这些文件在 `xafs_out/` 下，已被 `.gitignore` 忽略，只用于核对版式：

| 路径 | 说明 |
|---|---|
| `xafs_out/Pt1CeO2-flu_paper/Pt1CeO2_paper.png` | 六面板论文主图（a–f） |
| `xafs_out/Pt1CeO2-flu_paper/Pt1CeO2_paper_cn_table.{png,csv,md}` | CN 表图、CSV、Markdown |
| `xafs_out/Pt1CeO2-flu_paper/Pt1CeO2_paper_ft.csv` | \|FT\| 面板数据 |
| `xafs_out/Pt1CeO2-flu_report/Pt1CeO2-flu_mu_xas.png` | 五联报告图（a–e） |
| `xafs_out/Pt1CeO2-flu_report/diagnostics/s1–s8` | 9 张阶段诊断图 + `wt_sweep.csv` |

当前 `xafs_out/Pt1CeO2-flu_report/Pt1CeO2-flu_mu_fit.json` 内容为
`{"mode": "skipped"}`，因此该报告目录**没有** `*_fit_table.png`。这是
“无拟合结果时跳过拟合表”的正常行为，不是绘图失败。

`xafs_out/Pt1CeO2-flu_report/Pt1CeO2-flu_mu_compare.png` 不由当前代码生成；
它来自历史本地脚本，不能作为当前出图入口的产物。

---

## 6. 出图后检查

- 面板字母、坐标轴单位、图例、标线不压峰。
- 拟合面板同时有数据与模型，并标注拟合窗口、Rf、ΔE₀。
- 不同面板的 kweight 明确，峰高不跨 kweight 直接比较。
- CN 表逐壳层标状态；diagnostic 不写成正式 CN。
- 结构/CIF provenance 与 `disclosure_status` 已写入 CN 表；
  交付前运行 `python cn_disclosure.py <fit.json>`。
- 出图后目视检查，不只看脚本退出码。
