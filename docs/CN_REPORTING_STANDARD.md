# CN 报告与结构/CIF 披露标准

版本：1.0（2026-11-10）

适用范围：所有从 EXAFS 拟合输出配位数（CN）的图、表、CSV、JSON 和论文 /
Supplementary Information 草稿。本文件是项目级标准，位于 `docs/`，**不依赖
`xafs_out/` 下的任何案例目录**；案例只用于说明，不是运行依赖。

## 0. 一句话规则

每个 CN 数字必须同时回答四个问题：

1. 用的是哪个结构 / CIF？
2. 拟合参数是什么？
3. 误差和可信度是什么级别？
4. 哪些系统误差还没有校正？

缺任何一项，都不能把该数字写成正式 CN。

## 1. 必须披露的字段

| 字段 | 要求 |
|---|---|
| 结构来源类型 | `source_type`：CIF、constructed cluster、motif 等 |
| 结构引用 | `source_ref`：CIF 文件名 / COD ID，或构造方法 |
| CIF 路径 | 有 CIF 时写路径；没有时写 `cif=None`，输出 `CIF: none (no CIF used)` |
| CIF SHA256 | 有 CIF 时必须写；用于确认拟合时用的就是这份 CIF |
| 吸收位点 | `absorber_site`：晶位标签或 `frac:(x,y,z)`；构造团簇写中心原子 |
| Phase / model | `phase_id`、`model_level`（`crystal_cif` / `cluster` / ...） |
| FEFF 输入 | 每个 `feff.inp` 的路径 + SHA256；多壳层要全部列出 |
| S₀² | 数值 + 来源（静态表 / 显式 `--s02` / T5 / 化学态） |
| ΔE₀ | 数值；有误差时一并给出 |
| k 窗口 / kweight | `kmin`、`kmax`、`kweight`；每个面板/每次拟合分别注明 |
| R 窗口 | `rmin`、`rmax` |
| 自由度 | `N_idp`、自由变量数、R-factor |
| CN 状态 | 逐壳层标注，不能只给一个全局“fit result” |
| 系统误差预算 | `cn_systematics.relative_1sigma` 及主要分量 |

## 2. CN 状态分级

| 状态 | 含义 | 能否写进论文正文 |
|---|---|---|
| `formal` | 合格结构/CIF + 门禁通过 + 壳层可分辨 + 参数未压界 | 可以，需带误差与预算 |
| `diagnostic` | 构造/候选结构、第二壳或 MS 简并、相关强、仅存在性 | 只能作诊断，不能当正式 CN |
| `upper limit` | 未检出或只能给上限 | 写成“< X”或“未检出”，不能写 N=0 |
| `σ² at bound` | σ² 压到边界，N 与 σ² 不可分 | 该壳层降级，不报精确 N |

强制规则：

- 逐壳层状态必须可见；一张表不能只写一个全局状态。
- `diagnostic_only=True` 时，所有相关壳层都必须标成 diagnostic 或更低级别。
- 以下任一情况必须降级：手工构造团簇、第二壳/MS 分解不唯一、
  N–σ² 强相关、σ² 撞界、自吸收未校正。
- “没检测到”只能写成上限，不能写成 N=0。

## 3. 结构与 CIF 披露

### 3.1 CIF 派生结构

- 记录 COD ID、CIF 路径、SHA256、吸收位点。
- 必须经过 `cif_gate.py` / phase pack 的空间群、化学计量、吸收位点检查。
- 换 CIF、换相或换吸收位点后必须重新拟合，不能沿用旧 CN。

### 3.2 手工构造团簇

- 写清构造函数与参数，例如：
  `artemis.sac_cluster_sites(1.995 Å, O, CN=4, carbon_shell=False)`
  + 手动 Ce 壳层 @3.10 Å。
- `model_level="cluster"`，`provenance_status="constructed"`。
- 明确写“not derived from a CIF phase”，不能暗示有晶体学来源。
- 构造团簇的 CN 默认只能算 diagnostic，除非另有独立结构依据和门禁。

### 3.3 没有 CIF

- 用 `cif=None` 显式声明，输出必须出现 `CIF: none (no CIF used)`。
- 不能省略字段，也不能拿同元素的其他相 CIF 顶替。
- 没有合格 CIF/晶体模型时，不产生正式 CN；只能保留诊断级或拒绝行。

### 3.4 FEFF 输入

- 每个 `feff.inp` 都要记录路径和 SHA256；多壳层、多路径全部列出。
- FEFF 输入变了，原来的 CN provenance 就失效，必须重新拟合。

## 4. 拟合参数与误差

- **S₀² 固定并写来源**；不要为了“让 CN 好看”而放开或替换 S₀²。
- **形式误差与系统误差分开**：形式误差来自协方差矩阵；系统误差预算
  至少覆盖 S₀² 传递、壳层合并、N–σ²–ΔE₀ 相关、归一化。
- **报告 N_idp 和变量数**；变量数接近或超过 N_idp 时视为过拟合。
- **kweight 必须逐面板注明**；不同 kweight 的 |FT| 峰高不能直接比较。
- **R < Rbkg 不解释**；图上报出的峰位必须来自可解释 R 区间。
- **自吸收未校正**时，绝对 CN 只能给预算，不能给“精确值”。
- **第二壳与 MS 简并**时，只报存在性或上限，不报精确 N。

## 5. 运行时检查

### 5.1 独立检查任意 fit JSON

```powershell
python cn_disclosure.py path/to/fit.json
python cn_disclosure.py fit1.json fit2.json --warn-only
```

- 返回码 1 表示存在披露缺口；`--warn-only` 只打印、不返回失败。
- 只读 JSON，不重新拟合、不修改数值、不猜结构。
- `|Delta E0| > 5 eV` 只打印 `[NOTE]` 非阻断提示，不改变返回码；
  项目硬门禁仍是 `|Delta E0| < 15 eV`。该提示要求检查 FEFF 路径
  能量零点、实验 E0 约定和模型完整性，不能直接写成吸收边位移。

### 5.2 在 Python 中检查

```python
import paper_figures as pf

fit = pf.load_fit("sample_fit.json")
for issue in pf.validate_cn_disclosure(fit):
    print("WARN:", issue)
for note in pf.validate_cn_advisories(fit):
    print("NOTE:", note)   # 非阻断；当前为大 |Delta E0|
```

### 5.3 出图时强制门禁

```powershell
python paper_figures.py --sample sample_mu.dat --absorber M `
  --fit-json sample_fit.json --strict-disclosure --out paper_out
```

- 默认只打印警告，保持历史行为。
- `--strict-disclosure` 在披露不完整时直接报错；交付论文图建议开启。
- `*_cn_table.md` 会写出 `Disclosure warnings` 段；非阻断提示标为
  `(advisory, non-blocking)`。
- `*_cn_table.csv` 末尾有 `disclosure_status`、`disclosure_issues` 和
  `disclosure_advisories` 三列；`disclosure_status` 只反映阻断缺口。

## 6. 交付清单

- [ ] fit JSON 有 `structure` 块，或已用 `--structure-json` 补录
- [ ] CIF 路径 + SHA256，或显式 `cif=None`
- [ ] 构造团簇写清构造方法和所有壳层
- [ ] FEFF 输入路径 + SHA256 完整
- [ ] S₀² 数值 + 来源
- [ ] ΔE₀、kmin/kmax/kweight、rmin/rmax
- [ ] N_idp、变量数、R-factor
- [ ] 逐壳层 CN 状态
- [ ] CN 系统误差预算
- [ ] 限制说明：自吸收、第二壳/MS 简并、k 范围
- [ ] `cn_disclosure.py` 返回 0（或明确接受 `--warn-only` 的缺口）
- [ ] CN 表 PNG / CSV / MD / fit JSON 数值一致
- [ ] 出图后目视检查，不只看脚本退出码

## 7. 与现有模块的关系

| 模块 | 职责 |
|---|---|
| `cn_table.py` | 正式 CN 表结构与拒绝行；正式 CN 的权威 CSV |
| `motif_provenance.py` | motif 来源分类、CIF/SHA256 与 model level |
| `cif_gate.py` | CIF 空间群、化学计量、吸收位点门禁 |
| `paper_figures.py` | 论文图、CN 表显示、`validate_cn_disclosure()` |
| `cn_disclosure.py` | 本标准的独立检查器；也被 `paper_figures.py` 调用 |
| `xafs_figure_pack.py` | 已有输出目录的图包；不替代正式 CN 表 |

## 8. 参考案例（非依赖）

Pt₁/CeO₂@CeO₂ 案例使用手工构造的两壳层团簇：

- Pt–O4 方形平面，1.995 Å；
- 四个 Ce 邻居，3.10 Å；
- 没有 CIF，输出明确写 `CIF: none (no CIF used)`；
- CN 全部为 `diagnostic_only`，不写成正式配位数；
- 记录两个 `feff.inp` 的 SHA256。

示例产物在 `xafs_out/Pt1CeO2-flu_paper/`（分析产物，已被 `.gitignore`
忽略）。该目录只用于复现示例；本标准由 `docs/CN_REPORTING_STANDARD.md`
和 `cn_disclosure.py` 执行，不依赖它是否存在。
