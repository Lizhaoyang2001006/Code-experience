# CN_CIF 方案评审、修订版设计与执行交接

> 评审对象：`EXECUTION_PLAN.md`、`IMPLEMENTATION_CHECKLIST.md`、`README.md`  
> 评审仓库：`C:\Users\lizha\Desktop\Code-experience`  
> 日期：2026-11-07  
> 状态：评审与设计完成；执行交接见 §12；本文档不包含仓库代码改动

---

## 0. 结论摘要

### 0.1 总判断

原方案的方向是正确的：

- CN 必须有结构来源；
- 固定 N=N₀ 先检验结构，再放开 N；
- 没有合格来源时拒绝出 CN；
- 每个 CN 必须可追溯；
- XRD 只圈定候选相，EXAFS 只在候选集内排除。

但**不能按原方案原样落地**。主要问题是：

1. 把“没有用户 CIF”直接等同于“不能出 CN”，会误伤现有已验证的金属箔、标准样、B17–B20 催化模型和 74 项回归。
2. 把 motif 库整体降级为“不能产 CN”，属于过度反应。正确做法是升级 motif 的 provenance，而不是废弃整个库。
3. 方案里有若干科学/工程细节会导致错误拒绝或错误通过，包括第一壳 CN 范围、缺失 occupancy 的处理、固定 N 的 N_idp 公式、子壳合并与 Rayleigh 极限、S₀² 与系统误差缺失等。

### 0.2 推荐路线

采用“**正式 CN / 诊断拟合**”双轨：

- 现有 `--fit` 保持为诊断拟合，输出 `diagnostic_only`；
- 新增 `--formal-cn`（或 `--cn-mode strict`），只有它要求合格结构来源并写 `cn_table.csv`；
- 正式 CN 表与诊断拟合表分开，报告只从 `cn_table.csv` 读正式 CN；
- 先做非破坏性的 P0a，再做 motif provenance 升级、相解析、CIF 门禁、固定 N 一致性、条件 CN；
- motif 库实测仅 32 条、覆盖 9/17 个吸收元素，且晶体 CIF 均为现场生成文件；它只做候选索引，不作正式来源。P1a 增加 phase pack + coverage report，正式覆盖由用户 CIF、精确编号和本地种子包决定，不靠扩充 motif 数量。

这样既能立刻停止“正式 CN 被污染”，又不会破坏现有 74 项回归和 Phase V 验收链。

---

## 1. 当前代码中已确认的风险

| 风险 | 代码位置 | 现状 | 结论 |
|---|---|---|---|
| 自动 Top-1 进 CN | `xafs_fit.py:815-818` | `confident` 或 `--allow-ambiguous` 时取 `ranked[0]` | 真实风险：`confident` 只代表库内相对最优 |
| L1 合成 O/oct 团簇 | `structure_id.py` 的 `l1_first_shell_fit` | 默认 `ligand="O", geometry="oct"`，N 仍放开 | 真实风险：可能“可靠但化学错误” |
| 现场重建 CIF | `structure_lib.py:115-154` | 由 `asym_sites` 拼最小 CIF，`_atom_site_occupancy` 默认 1.0 | 真实问题：不是独立实验 CIF |
| 占据度缩放失败被吞 | `xafs_fit.py:220-234` | `except XAFSError: scale = []`，匹配不到就保持满占据 | 真实 bug：N 可能被高估 1/occ |
| 全壳 `vary_n=True` | `xafs_fit.py:945` | `vary_n=not fit_weights`，所有壳层 N 自由 | 真实风险：N 可补偿错误结构 |
| 报告直接读 fit 字典 | `xafs_report.py:351-354`、`607-616` | 从 `fit_out["fit"].paths` 取 N | 迁移时必须保留诊断表兼容 |
| 壳层合并固定 0.20 Å | `xafs_fit.select_shell_paths` | 按 Reff + 元素合并 | 与 CIF 子壳/ Rayleigh 极限需要重新协调 |
| N_idp 截断按 3n+1 | `xafs_fit.select_shell_paths` | 自由 N 模式 | 固定 N 模式需要单独公式 |

同时，当前代码已经具备一些可直接复用的能力：

- `artemis.fit_r_space(vary_n=[...])` 已支持只放开指定路径的 N（`artemis.py:1850-1857`）；
- `path_links` 已支持 `degen_ratio` 绑定 N 比例（`artemis.py:1804-1813`）；
- `xafs_fit.systematic_budget` 已实现 CN 系统误差预算（`xafs_fit.py:130-189`）；
- `structure_lib` 已有 CIF 展开、占据度保留和 FEFF 缓存；
- B18-1 已有 Rayleigh / Z 简并护栏，可直接用于 CN shell 分组；
- B17 已有 `estimate_nanoparticle_size`，可作为定量缺陷/纳米模型。

---

## 2. 设计原则（修订版）

### 2.1 数据与假设分离

- 数据：E0、χ(k)、FT、小波、XRD 峰位。
- 假设：相名、化学式、吸收位、配位几何、缺陷模型、CIF。
- CN 是“数据 × 假设”的条件量，不是数据本身。

### 2.2 正式 CN 与诊断拟合分离

| 输出 | 来源 | 是否进入正式 CN 表 |
|---|---|---|
| `formal_cn` | 合格结构模型 + 固定 N 一致性 + 条件 CN | 是 |
| `diagnostic_fit` | 现有自动 Top-1 / L1 / motif 拟合 | 否，仅诊断 |
| `legacy_auto_fit` | 旧行为，显式开关 | 否，标 `legacy_unverified` |
| `cluster_hypothesis` | SAC/DAC/合金/显式团簇 | 单独表，标 `cn_hypothesis_only` |
| `prototype` | 原型结构 | 单独表，标 `cn_model_dependent` |

### 2.3 结构来源分级

| 级别 | 来源 | CN 标签 | 正式表 |
|---|---|---|---|
| `crystal_cif` | 用户 CIF 或数据库实验 CIF + sha256 | `cn_conditional` | 是 |
| `crystal_edited` | CIF + 显式位点编辑 | `cn_conditional` | 是，附编辑记录 |
| `prototype` | 原型结构 + XRD/文献晶胞 | `cn_model_dependent` | 单独表 |
| `cluster` | 显式团簇/位点模型 | `cn_hypothesis_only` | 单独表 |
| `none` | 无合格来源 | 无 | 否，写拒绝原因 |

### 2.4 固定 N 优先

先锁 N=N₀，只放 ΔR、σ²、ΔE₀，回答“这份结构能不能解释这张谱”。  
通过后才允许第一壳 N 浮动。远壳 N 始终固定 N₀。

### 2.5 拒绝优先

程序可以拒绝出数，不可以编数。拒绝原因必须写入 JSON 和报告。

### 2.6 分辨率感知的 CN shell

CN 的“壳层”不是 CIF 的每一个距离，而是**在当前 k 窗口下可分辨的配位组**：

- `ΔR_Rayleigh = π / (2·(kmax − kmin))`；
- 低于 Rayleigh 极限的同元素子壳归为一个 CN 组；
- N₀ 求和，记录被合并的距离和 `sub_rayleigh_grouped=true`；
- 不把不可分辨的子壳分别放开 N。

这解决原方案“不允许合并”与 B18-1 的直接冲突。

### 2.7 可复现

- CIF 原始文件 + sha256 + 来源 URL/编号 + 检索日期；
- 结构编辑记录；
- 数据文件 sha256 + 处理参数 + S₀² 来源；
- 相同输入 + 相同缓存 → 相同输出 JSON 哈希（哈希不含时间戳/绝对路径）。

### 2.8 不做全库本地化

- 不镜像 COD/ICSD；
- 只缓存本次样品实际用到的候选 CIF；
- ICSD 由用户导出单条 CIF，走 `user_cif` 通道。

---

## 3. 修订版架构

### 3.1 模块图

```
run_xafs.py
  ├─ phase_resolver.py      相名/编号/CIF → 候选表（0 或 1）
  ├─ cif_gate.py            CIF 资格门禁 + N₀ 计算
  ├─ structure_model.py     结构模型：CIF / 编辑 / 原型 / 团簇
  ├─ structure_id.py        候选排序与诊断，不产正式 CN
  ├─ xafs_fit.py            诊断拟合 + 固定 N 一致性 + 条件 CN
  └─ xafs_report.py         正式 CN 表 + 诊断表 + 来源面板
```

### 3.2 数据流

```
用户输入：XRD 候选相 + 结构假设卡
  → phase_resolution.json（候选表 / 唯一选择 / 拒绝原因）
  → cif_gate.json（逐项门禁 + N₀ shell 组）
  → structure_model.json + model.cif / model.xyz
  → consistency_fit.json（固定 N=N₀ + 自由 N 对照）
  → cn_conditional.json（仅通过者）
  → cn_table.csv（正式 CN 唯一来源）
  → diagnostic_fit.json / fit_table.csv（诊断拟合，不进正式表）
```

### 3.3 输出目录

```
<outdir>/
  structures/
    raw/
      COD-2211652_ab12cd34.cif
      COD-2211652_ab12cd34.meta.json
    models/
      hematite_Fe1.cif
      hematite_Fe1.model.json
    phase_resolution.json
    cif_gate.json
    consistency_fit.json
    cn_conditional.json
  cn_table.csv
  diagnostic_fit.json
  fit_table.csv
  report.md / report.html
```

---

## 4. 接口契约

### 4.1 `StructureModel` 数据类

```python
@dataclass
class ShellN0:
    group_id: str
    shell_label: str
    scatterer: str
    r_cif: float
    n0: float
    n0_members: list[dict]       # [{r, n, occ, site_label}, ...]
    sub_rayleigh_grouped: bool
    merge_reason: str | None


@dataclass
class StructureEdit:
    op: str                      # substitute / vacancy / dope / split_cocupancy
    site: str
    old_value: object
    new_value: object
    reason: str
    approved_by: str
    timestamp: str


@dataclass
class StructureModel:
    model_id: str
    absorber: str
    absorber_site: str
    model_level: str             # crystal_cif / crystal_edited / prototype / cluster
    source_type: str             # user_cif / database_cif / prototype / cluster
    source_ref: str              # COD:xxxx / AMCSD:xxxx / user path / prototype id
    source_cif: str | None
    source_cif_sha256: str | None
    cell: dict | None
    spacegroup: str | None
    asym_sites: list | None
    cif_path: str | None
    occupancy_policy: str        # explicit / implicit_1.0 / rejected
    edits: list[StructureEdit]
    assumptions: dict
    n0_shells: list[ShellN0]
    warnings: list[str]
```

### 4.2 CLI 参数

| 参数 | 说明 | 默认 |
|---|---|---|
| `--phase-cif PATH` | 用户 CIF，可重复 | 无 |
| `--phase-id SOURCE:ID` | 精确编号，如 `COD:2211652` | 无 |
| `--phase SPEC` | `name=...;formula=...;sg=...;cell=...` | 无 |
| `--absorber-site SITE` | 多吸收位时必须指定 | 无 |
| `--build-model none\|prototype\|cluster` | 无 CIF 时是否允许构建模型 | `none` |
| `--prototype TYPE` | `rocksalt` / `spinel` / `corundum` / `rutile` / `perovskite` / `layered` | 无 |
| `--lattice SPEC` | 原型模型晶胞参数 | 无 |
| `--cn-mode off\|consistency-only\|conditional` | 正式 CN 模式 | `off` |
| `--formal-cn` | 写正式 `cn_table.csv` | 关闭 |
| `--legacy-auto-fit` | 旧行为，标 `legacy_unverified` | 关闭 |
| `--allow-computed` | 允许计算结构，降级为 hypothesis | 关闭 |
| `--allow-unusual-coordination` | 允许 CN 范围外，写警告 | 关闭 |
| `--cn-shell-tol FLOAT` | CN shell 分组容差（默认 Rayleigh 极限） | 自动 |
| `--defect-model JSON` | 定量缺陷/纳米模型 | 无 |
| `--phase-index PATH` | 本地相索引 | 内置 |
| `--offline` | 禁止网络 | 关闭 |

### 4.3 CLI 行为矩阵

| 用户命令 | 正式 CN | 诊断输出 | 说明 |
|---|---|---|---|
| `--fit` | 否 | 是 | 保持旧行为，标 `diagnostic_only` |
| `--fit --formal-cn` | 需结构来源 | 是 | 无来源 → `cn_table.csv` 0 行 + `no_cif` |
| `--fit --formal-cn --cn-mode conditional` | 是 | 是 | 固定 N 通过后只放开第一壳 N |
| `--fit --legacy-auto-fit` | 否 | 是 | 旧 Top-1/L1，标 `legacy_unverified` |
| `--identify` | 否 | 候选排序 | JSON 中不得出现 N |
| `--catalytic-model` | 否 | 是 | 标 `cn_hypothesis_only` |
| `--fit-weights` | 否 | 相分数 | 不产生 CN |

### 4.4 `cn_table.csv` 列

正式 CN 表必须包含：

```text
phase_id, model_level, source_type, source_ref, source_cif_sha256,
absorber_site, shell_group, shell_label, scatterer, r_cif,
n0, n0_members, n_fit, n_status, n_ratio,
n_err_formal, n_err_systematic, delta_r, sigma2, delta_e0, rf,
s02, s02_source, s02_calibrated,
kmin, kmax, kweight, rmin, rmax, n_idp, n_variables,
cn_label, gate_status, reject_reason, flags, assumptions,
data_sha256, software_version
```

说明：

- `n_status`: `fitted` / `fixed` / `diagnostic`；
- `n_err_systematic`: 相对 1σ，绝对值 = `n_err_systematic × n_fit`；
- 无 CN 时输出表头 + 一行 `gate_status=no_cif`，`reject_reason` 写原因。

### 4.5 JSON 输出

- `phase_resolution.json`：输入、候选表、选中条目、拒绝原因、缓存路径、sha256；
- `cif_gate.json`：逐项门禁、失败码、吸收位、忽略元素、occupancy 摘要、N₀ shell 组；
- `consistency_fit.json`：固定 N 结果、自由 N 对照、F-test、Rf、ΔR、σ²、ΔE₀、N_idp、变量数、通过/失败；
- `cn_conditional.json`：每个通过相的每个 CN shell 的完整字段；
- `diagnostic_fit.json`：现有拟合输出，不进入正式 CN 表。

哈希规则：

- JSON 排序键、固定浮点格式；
- 哈希不包含 `timestamp`、绝对路径、随机临时目录；
- `created_at` 单独字段，不参与哈希。

---

## 5. 关键算法

### 5.1 CIF 门禁

#### 解析

- 复用 `artemis.read_cif_sites` / `structure_lib.expand_full`，不重造 gemmi 展开；
- 必须有 `_atom_site_*` 循环；
- 记录原始 CIF 路径与 sha256。

#### 吸收位

- 吸收元素必须存在；
- 同元素多个晶位时，必须 `--absorber-site`；
- 禁止按 `occ × 多重度` 自动挑主位点；
- 未指定：`absorber_site_ambiguous`，不出正式 CN。

#### occupancy

- 显式 `_atom_site_occupancy`：使用；
- 缺失时：
  - 有序、单元素、化学式与占位一致：允许隐式 1.0，记录 `occupancy_source=implicit_1.0`；
  - 存在共占位、无序、化学式不匹配、或用户缺陷模型：`occupancy_missing`，不出正式 CN；
- 共占位：按元素拆分为不同分支，未拆分：`cocupancy_unresolved`；
- 占据度缩放失败：硬失败 `occupancy_scale_failed`，不允许按满占据继续。

#### 配位与键长

- 第一壳 CN 正常范围 **2–12**；12–16 警告；>16 或 <2 需 `--allow-unusual-coordination`；
- 金属箔 fcc/hcp 的 CN=12 必须正常通过；
- 键长合理区间用元素对共价/离子半径和 ± 容差；未知元素对要求用户提供参考键长；
- “多面体闭合”定义为：第一壳近邻数在合理范围、键长在合理区间、无异常空缺；不再使用模糊的“闭合”一词。

#### 计算结构

- `computed` 结构不得作为实验 N₀；
- 用户显式 `--allow-computed` 时，CN 标签 `cn_hypothesis_only`，写 `structure_is_computed=true`。

#### 失败码

```text
parse_failed
absorber_missing
absorber_site_ambiguous
polyhedron_not_closed
bond_length_implausible
occupancy_missing
occupancy_scale_failed
cocupancy_unresolved
cif_computed
cif_ambiguous
offline_miss
consistency_failed
n_not_identifiable
n_out_of_range
correlation_too_high
defect_model_missing
no_cif
legacy_unverified
```

### 5.2 N₀ 计算与 CN shell 分组

- 晶体：展开到 `rmax`（默认 6.0 Å），按 `(元素, 距离容差)` 分壳；
- `N₀ = Σ(简并度 × 占据度)`；
- 团簇：`N₀ = Σ(同壳位点数 × 占据度)`；
- CN shell 分组：
  - `ΔR_Rayleigh = π/(2·(kmax−kmin))`；
  - 同元素子壳间距 < `ΔR_Rayleigh` 时归为一个 CN 组；
  - 记录 `n0_members`、`sub_rayleigh_grouped=true`、`merge_reason`；
  - `N₀` 为组内求和；
- 不允许用“代表路径 + 0.20 Å 聚类”重新定义 CIF 的配位；0.20 Å 仅可作为路径选择容差，不能替代 CN shell 分组。

### 5.3 固定 N 一致性

```python
fit_fixed = artemis.fit_r_space(
    chi_k, ft, files,
    rmin=rmin, rmax=rmax, kmin=kmin, kmax=kmax,
    s02=s02,
    vary_delta_e0=True,
    vary_n=[],                 # 全部固定 N0
    path_degen=n0_list,
    ...
)
```

变量数：

- 每条壳层：ΔR、σ² = 2 个；
- 全局 ΔE₀ = 1 个；
- **`n_var = 2n + 1`**；
- 若启用 C3/C4、额外约束，按实际变量增减。

门禁：

- `n_var < n_idp`，否则先截断远壳或拒绝；
- `Rf ≤ 0.08` 通过，`0.08–0.10` 警告，`>0.10` 失败；
- `|ΔE₀| < 10 eV` 通过，10–15 警告，>15 失败；
- ΔR、σ² 在合理范围，不撞界；
- 最大参数相关 `|ρ| < 0.9`；
- 自吸收/荧光风险门禁；
- S₀² 必须已标定或显式给出。

自由 N 对照：

```python
fit_free = artemis.fit_r_space(..., vary_n=first_shell_indices, ...)
```

- 用 Hamilton F-test / AIC-BIC 比较固定 N 与自由 N；
- 若自由 N 显著改善（p<0.05）且无定量缺陷模型：`consistency_failed` 或 `structure_incomplete`；
- 若改善不显著：固定 N 通过，进入条件 CN。

### 5.4 条件 CN

只有固定 N 通过且 `--cn-mode conditional`：

```python
fit_cn = artemis.fit_r_space(
    ...,
    vary_n=first_shell_indices,   # 只放开第一壳
    path_degen=n0_list,
    path_links=first_shell_links, # 子壳 N 比例绑定
    ...
)
```

- 第一壳总 N：`N_fit = Σ first_shell_N`；
- `N/N₀ = N_fit / N₀_total`；
- 远壳 N 固定 N₀，`n_status=fixed`；
- 第一壳 `n_status=fitted`；
- `N/N₀` 0.7–1.3 通过；0.6–0.7 / 1.3–1.4 警告；<0.6 或 >1.4 需定量缺陷模型，否则 `n_out_of_range`；
- `N_err_formal > 25% N` 或系统误差过大：不给点值，给区间或 `n_not_identifiable`；
- `|corr(N, σ²)| > 0.9` 或最大相关 >0.9：给二维置信区间，不给点值。

### 5.5 多相处理

- 每个候选相独立解析、独立门禁、独立一致性检验；
- 输出 `phase_comparison` 表；
- 不按 Rf 最低自动选相；
- 只有一相通过：报告该相条件 CN，其他相写淘汰原因；
- 多相通过：分列 CN，不合并；
- 全部失败：不出 CN，输出诊断；
- 输出 `phase_discrimination_power`：若候选间只有第一壳可比，明确写“无法唯一指认相”。

### 5.6 定量缺陷/纳米模型

缺陷模型必须是结构化 JSON，例如：

```json
{"type": "nanoparticle", "diameter_nm": 2.5, "shape": "spherical"}
{"type": "vacancy", "site": "O1", "fraction": 0.1}
{"type": "surface", "particle_size_nm": 3.0, "shell": "first"}
```

要求：

- 用 B17 `estimate_nanoparticle_size` 或其他模型给出预期 `N/N₀` 与不确定度；
- 拟合 `N/N₀` 必须与模型预期一致；
- 没有定量预期 → `defect_model_missing`；
- 不接受自由文本“缺陷”声明。

---

## 6. 分阶段实施（修订版）

### P0a：正式 CN / 诊断双轨（非破坏）

**目标**：立刻停止“正式 CN 被污染”，但不翻默认、不破坏 74 项回归。

改动：

- `xafs_fit.fit_identified` 增加 `formal_cn: bool = False`；
- `formal_cn=False`：保持现有诊断输出，标 `diagnostic_only`；
- `formal_cn=True`：要求 `StructureModel`，无模型抛 `no_cif`；
- L1 输出加 `cn_status="diagnostic_only"`，正式表排除；
- occupancy 缩放失败改为硬报错 `occupancy_scale_failed`；
- 新增 `cn_table.csv` 写入函数与 provenance 字段；
- `xafs_report.py` 同时支持 `fit_table.csv`（诊断）和 `cn_table.csv`（正式）；
- 现有 `--fit` 行为不变。

验收：

- 无 CIF 时 `--fit --formal-cn` → 0 行 CN + `no_cif`；
- `--fit` 仍输出诊断拟合；
- `--identify` 无 N；
- 74 项回归全绿；
- 4 项压力测试全绿。

### P0b：motif provenance 升级

**目标**：让有来源的 motif 能作为正式结构来源，而不是整体废弃。

改动：

- `Motif`/JSON 增加 `source_ref`、`source_cif_sha256`、`absorber_site`、`occupancy_source`、`model_level`；
- 对每个晶体 motif：
  - 有源 CIF：校验晶胞/空间群/位点/距离与生成 CIF 一致；
  - 只有来源编号无源文件：标 `prototype` / `cn_model_dependent`；
  - constructed/SAC/DAC/alloy：标 `cluster` / `cn_hypothesis_only`；
- 保留现有 motif 库，不删除 B17–B20 条目。

验收：

- 至少 Fe/Co/Ni/Cu/Pt 金属箔、lepidocrocite、hematite、goethite、magnetite、CoO、Co₃O₄、Ni(OH)₂、CuO 的 provenance 可追溯；
- B17–B20 回归不降级；
- motif 生成的 CIF 与源 CIF 几何一致（有源时）。

### P1：覆盖基线、phase pack、相解析与缓存（本地优先）

**目标**：先量化“缺什么”，再解决“怎么取”；不依赖在线数据库也能工作。

P1a（覆盖基线 + phase pack）：

- 审计 motif 库，生成 `motif_coverage.json`；
- 新增 `_structure_lib/_phase_pack/`、`manifest.json`、`phase_index.json`；
- `tools/import_phase_pack.py`：导入用户/离线 CIF，校验 sha256，写索引；
- `tools/build_phase_index.py`：生成 `phase_index.json` 与 `phase_coverage.json`；
- 覆盖报告按吸收元素列出 verified / prototype / cluster / not_covered；
- 没有 verified CIF 的元素，正式 CN 只能走 `--phase-cif`，不得回退 motif。

P1b（相解析 + 缓存）：

- `phase_resolver.py`；
- `--phase-cif`、`--phase-id`、`--phase`、`--phase-index`、`--phase-pack`、`--offline`；
- 候选表、唯一性规则、缓存元数据；
- 在线 COD/AMCSD 仅作为可选 fetcher，默认关闭。

验收：

- 精确 ID 只取一条；
- 同名多条目返回候选表，0 CN；
- 离线未命中 `offline_miss`，不退回 motif；
- 缓存重复使用不重复下载；
- CLI 解析无 `eval/exec`；
- Fe/Co/Ni/Cu/Pt 每种至少 1 条 `verified_source_cif`；
- 未覆盖元素正式 CN 为 0 行 + `no_cif`。

### P2：CIF 门禁与固定 N 一致性

**目标**：有结构也不乱出 CN。

改动：

- `cif_gate.py`；
- `structure_model.py` Level A/C；
- 门禁：解析、吸收位、CN 范围、键长、occupancy、共占位、计算结构；
- N₀ 计算与 Rayleigh-aware CN shell 分组；
- 固定 N 一致性拟合；
- 自由 N 对照 + F-test/AIC-BIC；
- 失败码与 JSON；
- `xafs_fit.fit_identified` 接收 `StructureModel`。

验收：

- 错误相 `consistency_failed`，0 CN；
- 正确相通过，输出一致性指标；
- 缺 occupancy 按规则处理；
- 多晶位未指定 `absorber_site_ambiguous`；
- 共占位未拆分 `cocupancy_unresolved`；
- 计算结构降级/拒绝；
- fcc/hcp 金属 CN=12 正常通过。

### P3：条件 CN 与多相报告

**目标**：只放开第一壳 N，输出可复核的条件 CN。

改动：

- `cn_conditional.py`；
- `--cn-mode conditional`；
- 第一壳单一比例因子（`vary_n` + `path_links`）；
- `N₀/N_fit/N/N₀`；
- 多相对照表；
- 定量缺陷/纳米模型；
- `xafs_report.py` 来源面板 + 正式 CN 表。

验收：

- 通过者输出条件 CN；未通过者 0 CN；
- `N/N₀` 超范围且无定量缺陷模型：拒绝；
- 多相通过：分列，不自动选；
- `prototype` / `cluster` 单独表；
- 报告显示 CIF 来源、sha256、吸收位、N₀、S₀²、门禁状态。

### P4：验证与迁移

改动：

- 标准样、合成数据、离线、可复现测试；
- 文档、技能、README；
- `--legacy-auto-fit` 标记 deprecated；
- 移除对旧 `--fit-motif` 的正式 CN 依赖。

验收：

- 测试矩阵全部通过；
- 正式 CN 表 100% 行含结构/数据/S₀² provenance；
- 无 CIF 样品报告为空 CN + 原因；
- 现有 74 项回归在诊断模式下保持全绿。

---

## 7. 测试矩阵

### 7.1 单元测试

| 用例 | 断言 |
|---|---|
| `test_phase_spec_parser` | `name=...;formula=...;sg=...;cell=...` 解析正确；非法输入抛 `XAFSConfigError`；无 `eval/exec` |
| `test_cif_gate_parse_failed` | 坏 CIF → `parse_failed` |
| `test_cif_gate_absorber_missing` | 无吸收元素 → `absorber_missing` |
| `test_cif_gate_multiple_absorber_sites` | 多晶位未指定 → `absorber_site_ambiguous` |
| `test_cif_gate_occupancy_explicit` | 显式 occupancy → 正确 N₀ |
| `test_cif_gate_occupancy_implicit_ok` | 有序 CIF 缺 occupancy → 隐式 1.0，记录来源 |
| `test_cif_gate_occupancy_missing_disorder` | 无序/共占位缺 occupancy → `occupancy_missing` |
| `test_cif_gate_cocupancy_unresolved` | 共占位未拆分 → `cocupancy_unresolved` |
| `test_cif_gate_computed` | 计算结构 → `cif_computed` 或 hypothesis |
| `test_cn_range_metals` | fcc/hcp CN=12 通过；>16 需开关 |
| `test_n0_shell_grouping_rayleigh` | 1.95/2.10 Å 子壳在 k 3–12 下归一组，N₀ 求和 |
| `test_n0_occupancy_weighted` | N₀ = Σ( degeneracy × occupancy ) |
| `test_fixed_n_nidp_formula` | 固定 N：`2n+1`；条件 CN：`2n+2` |
| `test_provenance_hash_canonical` | 相同输入+缓存 → 相同 JSON 哈希；时间戳不影响 |
| `test_defect_model_schema` | 非法缺陷模型 → `XAFSConfigError` |

### 7.2 集成测试

| 用例 | 断言 |
|---|---|
| `test_no_phase_refuses_formal_cn` | `--formal-cn` 无来源 → 0 CN + `no_cif` |
| `test_motif_with_provenance_formal_cn` | 有源 CIF + sha 的 motif 可进正式表 |
| `test_motif_without_provenance_diagnostic_only` | constructed motif → 诊断表，0 正式 CN |
| `test_phase_id_unique_fetch_and_cache` | 精确 ID 只取一条；缓存命中不重复下载 |
| `test_phase_name_ambiguous_returns_candidates` | 多候选 → 0 CN + 候选表 |
| `test_phase_name_unique_with_spacegroup` | 名称+化学式+空间群唯一 → 进入门禁 |
| `test_wrong_phase_consistency_failed` | 赤铁矿 CIF 拟合针铁矿数据 → `consistency_failed` |
| `test_fixed_n_consistency_pass_then_conditional` | 正确相 → 固定 N 通过 → 条件 CN |
| `test_free_n_improves_significantly` | 自由 N 显著改善 → 结构不完整/拒绝 |
| `test_conditional_cn_out_of_range_requires_defect` | N/N₀ 超范围无定量缺陷模型 → 拒绝 |
| `test_conditional_cn_with_nanoparticle_model` | B17 粒径模型 → N/N₀ 预期一致 → 通过 |
| `test_multiphase_no_auto_pick` | 两相通过 → 分列，不自动选 |
| `test_prototype_labeled_model_dependent` | prototype → `cn_model_dependent` |
| `test_cluster_labeled_hypothesis` | cluster → `cn_hypothesis_only` |
| `test_offline_cache_hit_and_miss` | 命中可用；未命中 `offline_miss` |
| `test_legacy_auto_fit_marked_unverified` | legacy → `legacy_unverified`，不进正式表 |
| `test_l1_not_in_formal_table` | L1 诊断不出现于 `cn_table.csv` |
| `test_occupancy_scale_failure_hard_fail` | 缩放失败 → 硬报错，不按满占据 |
| `test_reproducible_output_hash` | 相同输入+缓存 → 输出哈希一致 |
| `test_report_reads_cn_table` | 报告只从 `cn_table.csv` 读正式 CN；无表显示原因 |

### 7.3 标准样与真实数据

| 样品 | 期望 |
|---|---|
| Fe 箔（bcc） | N₀=8；固定 N 检验通过；条件 CN 在误差内 |
| Co 箔（hcp） | N₀=12；固定 N 通过；fcc/hcp CN=12 不被误判 |
| Ni 箔（fcc） | N₀=12；固定 N 通过 |
| Cu 箔（fcc） | N₀=12；固定 N 通过 |
| 赤铁矿 | 合格 CIF 固定 N 检验；已知不可定量壳层保持拒绝 |
| 针铁矿 | 子壳 Rayleigh 分组；第一壳总 N₀=6 |
| 磁铁矿 | 多位点/混合相只用于相分数，不产正式 CN |
| CoO / Co₃O₄ / Ni(OH)₂ / CuO | B20 真实标准谱；有 CIF 时进入正式链 |
| 错误相 | 赤铁矿 CIF 拟合针铁矿 → `consistency_failed` |
| 缺 occupancy CIF | 按隐式/硬失败规则处理 |
| 纳米/团簇 | `cn_hypothesis_only`，不进正式表 |

### 7.4 回归与压力

- 现有 74 项契约回归：诊断模式下必须全绿；
- 4 项压力测试：`stress_test_{co,ni,cu,nio}.py` 全绿；
- 默认 LEPIDO 三扫描：`r_factor≈0.0273`、`N≈6.35`、`R≈2.002 Å`、7 列拟合表哈希不变；
- 正式 CN 新测试全部通过。

---

## 8. 迁移与兼容

### 8.1 保留旧行为

- `--fit` 不带 `--formal-cn`：现有诊断行为不变；
- `--fit-motif`：继续可用于诊断；
- `--allow-ambiguous`：继续影响诊断排序；
- `--fit-weights`：继续输出相分数，不产 CN；
- `--catalytic-model`：继续输出 `cn_hypothesis_only`。

### 8.2 新增正式链

- `--fit --formal-cn`：要求结构来源；
- `--cn-mode consistency-only`：只做固定 N 检验，不输出 N_fit；
- `--cn-mode conditional`：固定 N 通过后输出条件 CN；
- `--legacy-auto-fit`：旧行为，标 `legacy_unverified`。

### 8.3 报告

- `fit_table.csv`：诊断表，保留现有 7 列格式；
- `cn_table.csv`：正式表，新增 provenance 列；
- 报告中的“CN”只读正式表；诊断表单独标注；
- 无正式 CN 时显示拒绝原因，不显示空误差棒。

### 8.4 motif 库

- 保留为候选索引和结构来源；
- 有源 CIF + sha 的 motif 可进正式表；
- constructed/SAC/DAC/alloy 标 hypothesis；
- `_structure_lib/<元素>/_cif/*.cif` 若为生成文件，不得单独作为 `source_cif`，必须附源编号/哈希。

---

## 9. 需要拍板的产品决策

1. **正式 CN 是否默认严格？**  
   推荐：不翻默认，`--fit` 诊断，`--formal-cn` 严格。

2. **motif 库是否允许在 provenance 合格时作为正式来源？**  
   推荐：允许。只有 constructed/SAC/DAC/alloy 降级。

3. **第一版是否做 COD/AMCSD 在线检索？**  
   推荐：先不做。用户 CIF + 精确 ID + 本地索引优先；在线作为可选增强。

4. **缺失 occupancy 的隐式 1.0 是否允许？**  
   推荐：有序、单元素、化学式一致时允许，并记录；无序/共占位/缺陷模型时硬失败。

5. **Rf 阈值以哪套为准？**  
   推荐：一致性筛选 `≤0.08` 通过、`0.08–0.10` 警告、`>0.10` 失败；正式 CN 报告 `≤0.05` 通过、`0.05–0.08` 警告、`>0.08` 拒绝。

---

## 10. Definition of Done

一个批次完成必须同时满足：

- [ ] 代码改动可独立回滚。
- [ ] 新增/更新测试通过。
- [ ] 输出 JSON/CSV 字段有文档。
- [ ] 失败路径有明确拒绝码。
- [ ] 正式报告不出现无来源 CN。
- [ ] 诊断表与正式表分离。
- [ ] 现有 74 项回归在诊断模式下保持全绿。
- [ ] 4 项压力测试全绿。
- [ ] 不引入全库镜像。
- [ ] 不自动 Top-1 进正式 CN。
- [ ] L1 不进入正式 CN 表。
- [ ] 每个正式 CN 行含 CIF 来源 + sha256 + 吸收位 + N₀ + S₀² + 模型级别 + 数据哈希。
- [ ] 相同输入+缓存输出哈希一致。

---

## 11. 与原方案的差异总表

| 原方案 | 修订版 | 原因 |
|---|---|---|
| 无 CIF → 0 CN，默认 `--fit` 拒绝 | 正式/诊断双轨，`--formal-cn` 才严格 | 避免破坏 74 项回归和现有标准样 |
| motif 库降级为候选索引，不产 CN | motif provenance 升级，合格者可产 CN | 有源 CIF 的 motif 是合法结构来源 |
| CN 范围 2–8 | 2–12 正常，12–16 警告 | fcc/hcp 金属 CN=12 |
| 缺 occupancy 一律失败 | 有序可隐式 1.0，无序/共占位失败 | 避免大量假阴性 |
| 固定 N 仍用 3n+1 | 固定 N 用 2n+1；条件 CN 用 2n+2 | 变量数正确 |
| 不允许合并壳层 | Rayleigh-aware CN shell 分组 | 与 B18-1 一致，避免病态拟合 |
| Rf 单指标 | Rf + ΔR + σ² + ΔE₀ + F-test/AIC-BIC | 防止错误结构靠参数补偿 |
| 无 S₀² / 系统误差 | CN 表含 S₀²、系统误差、数据哈希 | CN 依赖 S₀² 和处理参数 |
| 缺陷模型自由声明 | 定量 JSON + B17 粒径模型 | 防止“声明即通过” |
| COD/AMCSD 在线检索 | 用户 CIF + 精确 ID + 本地索引优先 | 降低脆弱性和工程量 |
| 只追 CIF 来源 | CIF + 数据 + S₀² + 软件版本 | 完整可复现 |

---

## 12. 执行交接手册（给后续模型）

> 定位：§0–§11 是设计与验收标准；本节是执行顺序、命令、文件级改动与交接约束。两者冲突时，以本节的操作性约定为准，并在 `ACCEPTANCE_LOG.md` 记录原因。
>
> 本手册默认执行者只做当前批次，不提前实现后续批次。

### 12.1 启动指令（可直接复制给执行模型）

```text
你在 C:\Users\lizha\Desktop\Code-experience-cn 执行 CN_CIF 修订方案（先按 §12.2 从原目录复制并提交基线；原目录 C:\Users\lizha\Desktop\Code-experience 只读参考，不得直接修改）。
先完整阅读 C:\Users\lizha\Desktop\CN_CIF_执行方案\REVIEW_AND_REVISED_PLAN.md 的 §0–§5、§7–§12。
当前只执行 P0a；不要提前做 P0b/P1a/P1b/P2a/P2b/P3/P4。

硬约束：
1. 不改变现有 --fit 的数值行为。
2. 74 项契约回归 + 4 项压力测试必须全绿。
3. 默认 LEPIDO 三扫描的 r_factor/N/R 和 7 列拟合表 SHA256 不得变化。
4. 不改 absorber.py 的静态 s02；12 个未标定元素仍必须是 s02=None、calibrated=False。
5. 不用 eval/exec。
6. 不自动 Top-1 进正式 CN；L1 不进正式 CN 表。
7. 每批独立验收、独立回滚。
8. 工作区当前有未提交的 Phase V 改动；禁止 git checkout -- . / git reset --hard / git clean -fd。

完成后输出：改动文件清单、测试结果、LEPIDO 基线哈希、遗留问题、回滚点。
```

### 12.2 工作区隔离（第一步，必须做）

当前仓库不是干净状态：Phase V 改动未提交，`_structure_lib/{Au,Mo,Cr,Zn}`、部分氧化物 CIF/JSON、`_dynamic_cache` 是未跟踪产物。执行模型不得在未备份的情况下直接改原目录。

推荐做法：复制一份工作副本，在副本里建立基线提交，所有批次都在副本里完成；原目录保持不动。

```powershell
$src = "C:\Users\lizha\Desktop\Code-experience"
$dst = "C:\Users\lizha\Desktop\Code-experience-cn"
if (Test-Path $dst) { throw "目标已存在，先确认是否要复用" }
robocopy $src $dst /E /XD __pycache__ .pytest_cache /XF *.pyc
if ($LASTEXITCODE -gt 7) { throw "robocopy 失败: $LASTEXITCODE" }
Set-Location $dst
git add -A
git commit -m "baseline: Phase V working tree before CN formalization"
git status --short
```

说明：

- `robocopy` 退出码 0–7 都表示成功，只有 `>7` 才是失败。
- 基线提交只在副本里做，不推送、不改原目录。
- 如果用户明确要求在原始目录执行，必须先做完整备份或基线提交；仍然禁止 `git checkout -- .`、`git reset --hard`、`git clean -fd`。
- 每批完成后在副本里单独提交，回滚用 `git revert <commit>`，不要 reset 到基线。

### 12.3 冻结基线（改代码前必须记录）

```powershell
$PY = "C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe"

# 1) 全量 CI：74 项契约回归 + 4 项压力测试
& $PY run_ci.py --stress

# 2) 检查 CI 报告
& $PY -c "import json,sys; d=json.load(open('ci_report.json',encoding='utf-8')); print({'ok':d['ok'],'contract':d['contract_tests'],'stress':d['stress_tests']}); sys.exit(0 if d['ok'] and d['contract_tests']==74 and d['stress_tests']==4 else 1)"

# 3) 默认 LEPIDO 三扫描基线
& $PY xafs_report.py data/练手数据/02_Fe_standards/Fe_lepidocrocite.000 data/练手数据/02_Fe_standards/Fe_lepidocrocite.100 data/练手数据/02_Fe_standards/Fe_lepidocrocite.200 --absorber Fe --no-strict -o xafs_out/_baseline

# 4) 7 列拟合表 SHA256
Get-FileHash "xafs_out\_baseline\Fe_lepidocrocite_fit_table.csv" -Algorithm SHA256

# 5) #60 静态 s02 检查：12 个元素必须仍为 None/False
& $PY -c "import absorber; syms=['Cr','Mn','Zn','Mo','Ru','Rh','Pd','Ag','Ir','Au','Ce','U']; bad=[(s,absorber.get_absorber(s).s02,absorber.get_absorber(s).calibrated) for s in syms if absorber.get_absorber(s).s02 is not None or absorber.get_absorber(s).calibrated]; print('OK' if not bad else bad); raise SystemExit(1 if bad else 0)"
```

必须记录到 `ACCEPTANCE_LOG.md` 的基线值：

- `ci_report.json`：`ok=true`、`contract_tests=74`、`stress_tests=4`。
- LEPIDO：`r_factor=0.027281059410161693`。
- 第一壳：`N=6.348637726`、`R=2.0019558445 Å`。
- 7 列拟合表 SHA256：`E098D0EDC24E49296710E44518F8A8F010DCD03D7D7EB8CB2E26831A7532E760`。
- #60：12 个未标定元素的静态 `s02 is None`、`calibrated is False`。

任何一项不匹配，先停下排查，不要继续改代码。

### 12.4 每批通用流程

1. 只读当前批次对应的设计段落，确认不提前实现后续批次。
2. 先加失败测试：追加到 `regression_test.py` 的 `RegressionTest` 末尾，继续编号（#75 起），不改旧测试编号。
3. 做最小实现，不改无关文件。
4. 单测：
   ```powershell
   & $PY -m unittest regression_test.RegressionTest.test_name -v 2>&1 | Out-String
   ```
   注意：`unittest` 默认写 stderr，PowerShell 可能显示 NativeCommandError；以 `Ran 1 test` + `OK` 为准，最终门禁以 `run_ci.py` 退出码为准。
5. 全量门禁：
   ```powershell
   & $PY run_ci.py --stress
   ```
6. 重跑 §12.3 的 LEPIDO 基线命令与 SHA256。
7. 把结果写入 `C:\Users\lizha\Desktop\CN_CIF_执行方案\ACCEPTANCE_LOG.md`。
8. 在副本里提交：
   ```powershell
   git add -A
   git commit -m "P0a: formal CN / diagnostic dual track"
   git format-patch -1 --stdout > "C:\Users\lizha\Desktop\CN_CIF_执行方案\patches\P0a.patch"
   ```
9. 回滚：`git revert <commit>`，禁止 `git reset --hard`。

### 12.5 P0a 执行卡：正式/诊断双轨止血

**目标**

- 立即停止“正式 CN 被污染”，但**不翻默认**。
- `--fit` 继续输出原有诊断拟合，数值、CSV 列、HTML 默认输出不变。
- 新增 `--formal-cn`；没有合格结构模型时正式表 0 行 CN，写拒绝原因。
- 占据度缩放失败改为硬报错，不再按满占据继续。
- L1 只做诊断，永不进入正式 CN 表。

**非目标**

- 不做相解析、不做 CIF 门禁、不做固定 N 一致性、不做条件 CN。
- 不删除 motif、不改 `structure_id.identify_result` 的 Top-1/L1 自动回退行为。
- 不改 `absorber.py`、不改 `vary_n=not fit_weights` 的诊断路径、不改任何默认阈值。

**文件级改动**

| 文件 | 改动 |
|---|---|
| `structure_model.py`（新增） | `CNError`、`ShellN0`、`StructureEdit`、`StructureModel`、`sha256_file`、`model_from_motif`、`validate_formal` |
| `cn_table.py`（新增） | 正式表列定义、`write_cn_table`、`read_cn_table`、`refusal_row` |
| `xafs_fit.py` | `_ss_rows` 占据度失败硬报错；`fit_identified`/`run` 增加 `formal_cn` 参数；输出 `formal_cn` 块与诊断标记 |
| `run_xafs.py` | 新增 `--formal-cn`；`process_sample` 透传；单样品/批处理都写 `cn_table.csv` |
| `xafs_report.py` | 新增 `--formal-cn` 透传；生成/读取 `cn_table.csv`；仅在请求正式 CN 时增加正式段 |
| `xafs_config.py` | 建议同步 `formal_cn` 配置键；若不做，必须在文档中注明 CLI-only |
| `regression_test.py` | 新增 4 项测试，追加在 #74 之后 |

**接口契约**

```python
# structure_model.py
class CNError(XAFSError):
    def __init__(self, code: str, message: str):
        self.code = code
        super().__init__("CN", f"{code}: {message}")


@dataclass
class ShellN0:
    group_id: str
    shell_label: str
    scatterer: str
    r_cif: float
    n0: float
    n0_members: list[dict]
    sub_rayleigh_grouped: bool
    merge_reason: str | None


@dataclass
class StructureEdit:
    op: str
    site: str
    old_value: object
    new_value: object
    reason: str
    approved_by: str
    timestamp: str


@dataclass
class StructureModel:
    model_id: str
    absorber: str
    absorber_site: str | None
    model_level: str          # crystal_cif / crystal_edited / prototype / cluster
    source_type: str          # user_cif / database_cif / prototype / cluster / none
    source_ref: str | None
    source_cif: str | None
    source_cif_sha256: str | None
    cell: dict | None
    spacegroup: str | None
    asym_sites: list | None
    cif_path: str | None
    occupancy_policy: str     # explicit / implicit_1.0 / rejected
    edits: list[StructureEdit]
    assumptions: dict
    n0_shells: list[ShellN0]
    warnings: list[str]
    provenance_status: str    # verified_source_cif / source_ref_only / generated_only / constructed / unverified


def sha256_file(path: str | Path) -> str: ...
def model_from_motif(motif, *, source_type=None, source_ref=None,
                     source_cif_sha256=None, absorber_site=None,
                     occupancy_policy=None, model_level=None) -> StructureModel: ...

# cn_table.py
CN_TABLE_COLUMNS = [...]
def write_cn_table(path: str | Path, rows: list[dict]) -> None: ...
def read_cn_table(path: str | Path) -> list[dict]: ...
def refusal_row(gate_status: str, reject_reason: str, **meta) -> dict: ...

# xafs_fit.py
def fit_identified(..., formal_cn: bool = False,
                   structure_model: StructureModel | None = None) -> dict: ...
def run(..., formal_cn: bool = False,
        structure_model: StructureModel | None = None) -> dict: ...
```

**行为约定**

- `formal_cn=False`：完全保持现状。只允许在返回 JSON 中增加 `diagnostic_only=true`、`cn_status="diagnostic_only"`；不得改 `fit`、`fit.paths`、7 列 CSV、HTML 默认输出。
- `formal_cn=True`：
  - 无 `structure_model`：`formal_cn.status="no_cif"`，诊断拟合照常完成；`cn_table.csv` 写一行拒绝记录，CN 数值列全空。
  - 有模型但 P0a 未实现一致性：`formal_cn.status="pending_consistency"`，仍然不写任何 CN 数值。
  - L1：`formal_cn.status="l1_diagnostic_only"`，L1 结果只留在诊断 JSON/表，不写正式表。
  - P0a 绝不输出 `n_fit` 点值；正式 CN 从 P3 才开始。
- `formal_cn` 块建议结构：

```python
out["formal_cn"] = {
    "status": "no_cif" | "pending_consistency" | "l1_diagnostic_only" | "ok",
    "rows": [],                       # P0a 永远为空
    "reject_reason": "..." | None,
    "model": model.to_json() if model else None,
}
out["diagnostic_only"] = True
out["cn_status"] = "diagnostic_only"  # L1/诊断路径
```

- `xafs_fit.run` 不因正式链拒绝而抛异常；拒绝用 `formal_cn.status` 表达，保证诊断输出始终存在。`CNError` 只用于真正的配置/数据错误（例如 `occupancy_scale_failed`）。
- `cn_table.csv` 只在 `--formal-cn` 开启时写；没有 CN 时写表头 + 一行 `gate_status=no_cif` 的拒绝记录。不要把拒绝记录当 CN 行统计。
- 报告中的“正式 CN”段只读 `cn_table.csv`；诊断拟合表继续标注为诊断，不得让读者把诊断 N 当正式 CN。

**占据度硬报错**

位置：`xafs_fit._ss_rows` 的 `needs_occ` 分支（约 216–235 行）。

- `cif_shell_occupancy_scale` 抛异常 → 抛 `CNError("occupancy_scale_failed", ...)`，不要 `scale=[]`。
- 返回空列表同样视为失败：`raise CNError("occupancy_scale_failed", ...)`。
- 不得在 `select_shell_paths` 或上层静默吞掉该错误。
- 单条路径匹配不到占据度时，保留原简并度但记录 `occupancy_match_missing=true`；整表失败才硬报错。正式链的 N₀ 由 P2 的 `cif_gate` 独立计算，不依赖这里的近似。

**P0a 新增测试**

| 测试 | 断言 |
|---|---|
| `test_formal_cn_refusal_without_model` | `xafs_fit.run(..., formal_cn=True)` 返回 `formal_cn.status=="no_cif"`；诊断拟合仍存在；`cn_table.csv` 只有拒绝行、`n_fit` 为空 |
| `test_formal_cn_l1_not_in_table` | 用 `FERRIHYDRITE` 触发 L1；`mode=="l1_first_shell"`；`formal_cn.status=="l1_diagnostic_only"`；正式表 0 个 CN 行 |
| `test_occupancy_scale_failure_hard_fail` | 对触发 `needs_occ` 的 motif 打补丁使 `cif_shell_occupancy_scale` 失败；`_ss_rows` 抛 `occupancy_scale_failed`，不再按满占据返回 |
| `test_cn_table_schema_and_refusal_row` | `CN_TABLE_COLUMNS` 顺序正确；拒绝行 `gate_status=no_cif`、`n_fit==""`、`reject_reason` 非空；写读往返一致 |

**P0a 验收**

- `run_ci.py --stress`：`ok=true`、`contract_tests=78`、`stress_tests=4`。
- LEPIDO 基线：`r_factor`、`N`、`R`、7 列拟合表 SHA256 全部不变。
- `--fit --formal-cn` 无来源：`cn_table.csv` 存在且只有拒绝行；`fit.csv`/`fit.json` 仍生成。
- 不带 `--formal-cn` 的 `--fit`：不生成 `cn_table.csv`，7 列拟合表哈希不变。
- `--identify` 仍不输出 N。
- #60 静态 `s02` 检查通过。

**回滚**

`git revert <P0a commit>`。P0a 不写任何持久化结构缓存，不需要数据迁移。

### 12.6 P0b 执行卡：motif provenance 升级

**目标**：让有真实来源的 motif 能作为正式结构来源；不删除现有 motif，不影响 B17–B20。

**字段**

在 `structure_lib.Motif` 增加并持久化：

```python
source_ref: str | None = None          # COD:xxxx / AMCSD:xxxx / 用户文件标识
source_cif_sha256: str | None = None   # 仅当 source_cif 是真实来源文件时填写
absorber_site: str | None = None       # 吸收位点标签/分数坐标标识
occupancy_source: str | None = None    # explicit / implicit_1.0 / generated
model_level: str | None = None         # crystal_cif / crystal_edited / prototype / cluster
provenance_status: str = "unverified" # verified_source_cif / source_ref_only / generated_only / constructed
edits: list | None = None
```

**分类规则**

| 条件 | `provenance_status` | `model_level` | 是否进主正式表 |
|---|---|---|---|
| 有真实 CIF + sha256 匹配 + 吸收位明确 | `verified_source_cif` | `crystal_cif` / `crystal_edited` | 是 |
| 只有 COD/AMCSD 编号，没有源 CIF | `source_ref_only` | `prototype` | 否，单独 prototype 表 |
| `_structure_lib/*/_cif/*.cif` 现场生成的小 CIF | `generated_only` | `prototype` | 否 |
| `kind=="cluster"` / SAC / DAC / alloy / constructed | `constructed` | `cluster` | 否，`cn_hypothesis_only` |
| 其他 | `unverified` | `prototype` 或 `cluster` | 否 |

关键点：

- `source.ref="COD-9002763"` 但 `cif_path=None` 的 motif（例如当前 lepidocrocite）只能是 `source_ref_only`/`prototype`，不能仅凭编号宣称 `crystal_cif`。
- `_structure_lib/Fe/_cif/lepidocrocite.cif` 这类 300–500 字节的现场生成 CIF 是 `generated_only`，不得作为 `source_cif_sha256` 的来源。
- 只有 P1 缓存到真实 COD/AMCSD CIF，或用户提供 `--phase-cif` 后，才能升级为 `verified_source_cif`。

**实现**

- 新增 `tools/backfill_motif_provenance.py`（或等价脚本）：遍历 `_structure_lib/**/*.json`，补齐字段，写出 `provenance_report.json`。
- 同时生成 `motif_coverage.json`：按吸收元素统计 motif 数、真实 CIF 数、provenance 分布、完全缺失元素；该报告只用于覆盖审计，不参与正式 CN。
- 脚本必须可重复执行；已存在且哈希匹配的字段不重复改写。
- `model_from_motif` 按上表分类，不要根据 `kind=="crystal"` 就直接当正式来源。
- `validate_formal()` 只允许 `provenance_status=="verified_source_cif"` 且 `source_cif_sha256` 与文件当前哈希一致；否则返回 `no_cif`/`unverified_source`/`source_hash_mismatch`。

**P0b 新增测试**

| 测试 | 断言 |
|---|---|
| `test_motif_provenance_schema` | 全部 motif 可加载；新字段存在；`provenance_status` 在允许集合内 |
| `test_generated_cif_not_formal_source` | lepidocrocite 等生成 CIF 的 motif 为 `generated_only`；`validate_formal()` 拒绝 |
| `test_user_cif_is_verified_source` | 临时真实 CIF + sha256 + 吸收位 → `verified_source_cif`；`validate_formal()` 通过 |
| `test_b20_motifs_not_degraded` | Co/Ni/Cu/Pt 的 B17–B20 motif 仍存在、几何指纹不变；constructed/SAC/DAC/alloy 标 `constructed` |

**P0b 验收**

- `run_ci.py --stress`：`ok=true`、`contract_tests=82`、`stress_tests=4`。
- LEPIDO 基线与 #60 不变。
- `provenance_report.json` 中每一条 motif 都有明确状态；没有“未分类”。
- #74 与 4 项压力测试 top1 零退化。

### 12.7 P1 执行卡：覆盖基线、phase pack、相解析与缓存

**目标**：先量化 motif/相覆盖，再解析 CIF；不依赖在线数据库也能工作。

#### P1a：motif 覆盖审计 + phase pack

**背景**：当前 motif 库 32 条、覆盖 9/17 个吸收元素；23 条晶体 motif 的 CIF 均为 339–415 B 现场生成文件，不能作为正式来源。详见 §13。

**文件**

- 新增 `_structure_lib/_phase_pack/` 目录与 `manifest.json`。
- 新增 `tools/import_phase_pack.py`。
- 新增 `tools/build_phase_index.py`。
- 新增 `motif_coverage.json`、`phase_coverage.json`、`phase_index.json`。

**Phase pack 规范**

- 目录：`_structure_lib/_phase_pack/<source>/<id>/<sha256>.cif` + `.meta.json`。
- 每条必须有 `source`、`id`、`formula`、`spacegroup`、`cell`、`absorber_sites`、`cif_sha256`、`experimental`、`model_level`。
- `import_phase_pack.py` 必须校验 CIF 可解析、sha256 与 meta 一致、吸收元素/位点存在。
- 幂等：同一 `(source, id, sha256)` 重复导入不重复写。
- 不同 sha256 的同 ID 冲突必须报错，除非显式 `--replace`。
- 不下载、不联网；只导入用户给定目录。

**覆盖规则**

- 5 个已验证元素 Fe/Co/Ni/Cu/Pt：每种至少 1 条 `verified_source_cif`。
- 其余 12 个吸收元素允许 `not_covered`，但 `phase_coverage.json` 必须明确列出。
- 没有 verified CIF 的元素，正式 CN 只能走 `--phase-cif`；不得回退 motif。
- generated CIF 只能算 `prototype`，不得计入 verified 覆盖。

**P1a 新增测试（4 项）**

`test_motif_coverage_audit`、`test_phase_pack_import_and_index`、`test_phase_pack_sha256_mismatch_rejected`、`test_uncovered_element_no_formal_cn`。

**P1a 验收**

- `contract_tests=86`、`stress_tests=4`。
- `motif_coverage.json` 覆盖全部 32 条 motif，无“未分类”。
- `phase_coverage.json` 明确 9 个有 motif 元素与 8 个无 motif 元素。
- 至少 Fe/Co/Ni/Cu/Pt 各 1 条 verified CIF（由用户提供或 P1a 种子包导入）。
- 未覆盖元素正式 CN 为 0 行 + `no_cif`。
- LEPIDO 基线与 #60 不变。

#### P1b：相解析与缓存

**文件**

- 新增 `phase_resolver.py`。
- `phase_index.json`（由 P1a 生成）。
- `_structure_lib/_phase_cache/<source>/<id>/<sha256>.cif` + `meta.json`。

**CLI**

`--phase-cif PATH`（可重复）、`--phase-id SOURCE:ID`、`--phase SPEC`、`--phase-index PATH`、`--phase-pack PATH`、`--absorber-site SITE`、`--offline`。

**解析优先级**

1. 用户 `--phase-cif`。
2. 精确 `--phase-id`（本地缓存/phase pack 优先；未命中且非离线才允许 fetcher）。
3. 本地 `phase_index.json` 唯一命中。
4. `--phase SPEC` 唯一命中。
5. 多候选：返回候选表，`cn_table.csv` 0 行 + `cif_ambiguous`，不自动选。
6. 离线未命中：`offline_miss`，不退回 motif。

**安全与可复现**

- `--phase SPEC` 只用 `split(";")` / `split("=", 1)` 解析 `name/formula/sg/cell`；禁止 `eval`/`exec`。
- `cell` 必须解析为 6 个有限浮点数；非法输入抛 `XAFSConfigError`。
- 缓存复用前校验 sha256；不匹配则拒绝并重新获取。
- 网络 fetcher 必须依赖注入；CI 默认 `fetcher=None` + `--offline`，测试用 stub，不发真实网络请求。
- 不做 COD/ICSD 全库镜像；只缓存实际用到的条目。

**P1b 新增测试（7 项）**

`test_phase_spec_parser_no_eval`、`test_phase_cif_user_source`、`test_phase_id_cache_hit`、`test_phase_id_cache_miss_offline`、`test_phase_name_ambiguous`、`test_phase_index_schema`、`test_cache_sha256_reuse`。

**P1b 验收**

- `contract_tests=93`、`stress_tests=4`。
- 精确 ID 只取一条；多候选 0 CN；离线可复现；缓存不重复下载。
- LEPIDO 基线与 #60 不变。

### 12.8 P2 执行卡：CIF 门禁与固定 N 一致性

建议拆成 P2a（门禁 + N₀）和 P2b（一致性拟合），分别验收。

#### P2a：CIF 门禁 + N₀ + Rayleigh 分组

**文件**：新增 `cif_gate.py`；扩展 `structure_model.py`。

**门禁检查**

- CIF 可解析；必须有 `_atom_site_*` 循环。
- 吸收元素存在。
- 同元素多个晶位时必须 `--absorber-site`；不得按 `occ × 多重度` 自动挑主位点。
- occupancy：
  - 显式 `_atom_site_occupancy`：使用。
  - 缺失 occupancy：有序、单元素、化学式一致 → 隐式 1.0，记录 `occupancy_source="implicit_1.0"`。
  - 无序/共占位/缺陷模型/化学式不匹配 → `occupancy_missing`，拒绝。
  - 共占位未拆分 → `cocupancy_unresolved`。
- 计算结构（computed）不得作为实验 N₀；用户显式 `--allow-computed` 时降级为 `cn_hypothesis_only`。
- 第一壳 CN 正常范围 2–12；12–16 警告；>16 或 <2 需 `--allow-unusual-coordination`。fcc/hcp 的 CN=12 必须正常通过。
- 键长合理区间用元素对共价/离子半径和 ± 容差；未知元素对要求用户提供参考键长。
- 失败码按 §5.1 的列表返回。

**N₀ 与 CN shell 分组**

- 用 `structure_lib.expand_full`（或 `artemis.read_cif_sites(..., with_occupancy=True)`）拿到近邻 + 占据度；**不要**用默认 4 元组模式做定量。
- `N₀ = Σ(简并度 × 占据度)`。
- `ΔR_Rayleigh = π/(2·(kmax−kmin))`；同元素子壳间距低于 Rayleigh 时归为一个 CN 组，N₀ 求和，记录 `n0_members`、`sub_rayleigh_grouped=true`、`merge_reason`。
- 不允许用 0.20 Å 的路径聚类直接定义 CN；0.20 Å 只用于 FEFF 路径选择。
- 读取 CIF 时要显式检查 `_atom_site_occupancy` 标签是否存在；`artemis.read_cif_sites` 对缺失 occupancy 默认 1.0，不能用来判断“缺失”。

**P2a 测试（至少 8 项）**

`test_cif_gate_parse_failed`、`test_cif_gate_absorber_missing`、`test_cif_gate_absorber_site_ambiguous`、`test_cif_gate_occupancy_implicit_ok`、`test_cif_gate_occupancy_disorder_rejected`、`test_cif_gate_cocupancy_unresolved`、`test_cn_range_12_metals`、`test_n0_shell_grouping_rayleigh`。

#### P2b：固定 N 一致性 + 自由 N 对照

**文件**：`xafs_fit.py`、`artemis.py`（尽量少改）、`cn_table.py`。

**路径映射**

- 新增 `map_n0_to_plan(plan, model)`：
  - 对 `plan["files"]` 的每条 FEFF 路径，用 `(scatterer, reff)` 匹配 `ShellN0` 的成员。
  - 一个 CN 组可对应多条路径（子壳）；每条路径的 `path_degen` 用该路径自己的 N₀，不要用组总和。
  - 同一 CN 组的多条路径用 `path_links.degen_ratio = N0_j / N0_ref` 绑定，只放开参考路径的 N。
  - 匹配不到 / 多义 → `n0_mapping_failed`，不出正式 CN。

**固定 N 拟合**

```python
fit_fixed = artemis.fit_r_space(
    res.chi_k, res.ft, plan["files"],
    rmin=..., rmax=..., kmin=..., kmax=..., kweight=...,
    s02=s02, vary_delta_e0=True,
    vary_n=[],                 # 全部固定 N0
    path_degen=n0_per_path,
    ...)
```

- `path_degen` 长度必须等于路径数；`vary_n=[]`。
- 变量数：每条壳层 ΔR + σ² = 2，加全局 ΔE₀ = 1，即 `n_var = 2n + 1`。
- `select_shell_paths` 增加信息量模式：
  - 诊断/自由 N（默认）：`3n+1`，保持现状。
  - 固定 N：`2n+1`。
  - 条件 CN：`2n+2`（第一壳多一个 N）。
- 不要直接复用现有 `3n+1` 截断，否则固定 N 会不必要地丢掉远壳。

**自由 N 对照**

```python
fit_free = artemis.fit_r_space(
    ..., vary_n=first_shell_indices,
    path_degen=n0_per_path,
    path_links=first_shell_links, ...)
```

- 用 `artemis.shell_f_test(fit_fixed, fit_free)` 做 Hamilton F-test；它已实现 `F = [(χ²₁−χ²₂)/(p₂−p₁)] / [χ²₂/(N_idp−p₂)]`。
- `confidence ≥ 0.95` 且 `significant=True`：自由 N 显著改善，说明结构不完整/错误或存在缺陷；无定量缺陷模型时 `consistency_failed` / `structure_incomplete`。
- 不显著：固定 N 通过，进入条件 CN。
- 同时门禁 Rf、ΔE₀、ΔR、σ² 撞界、最大参数相关、自吸收、S₀² 已标定。

**P2b 测试（至少 2 项）**

`test_fixed_n_nidp_formula`（`2n+1` / `2n+2`）、`test_wrong_phase_consistency_failed`（错误相 0 CN；正确相通过）。可再加 `test_free_n_improves_significantly`。

**P2 验收**

- P2a 完成：`contract_tests=101`；P2b 完成：`contract_tests=103`。
- 错误相 `consistency_failed`，0 CN；正确相通过并输出一致性指标。
- fcc/hcp CN=12 正常通过；缺 occupancy 按规则处理；多晶位/共占位按失败码拒绝。
- LEPIDO 基线与 #60 不变。

### 12.9 P3 执行卡：条件 CN、多相与缺陷模型

**目标**：固定 N 通过后，只放开第一壳总 N，输出可复核的 `N₀ / N_fit / N/N₀`。

**第一壳单一比例因子**

- 选第一壳参考路径 `ref`（建议 N₀ 最大的路径）。
- `path_degen` 用每条路径自己的 N₀。
- `vary_n=[ref_idx]`，只放开参考路径 N。
- 其他第一壳路径用 `path_links` 绑定：

```python
{"target": j, "source": ref_idx,
 "degen_ratio": N0_j / N0_ref,
 "share_deltar": False,      # 关键：默认 True 会把 ΔR 也绑死
 "share_sigma2": False}      # 关键：默认 True 会把 σ² 也绑死
```

- `N_fit_total = Σ first_shell_paths.N`；`N0_total = Σ N0_i`；`N/N0 = N_fit_total / N0_total`。
- 远壳 N 固定 N₀，`n_status="fixed"`；第一壳 `n_status="fitted"`。

**门禁**

- `N/N₀` 0.7–1.3 通过；0.6–0.7 / 1.3–1.4 警告；<0.6 或 >1.4 需要定量缺陷模型，否则 `n_out_of_range`。
- `N_err_formal > 25% N` 或系统误差过大：不给点值，给区间或 `n_not_identifiable`。
- 最大参数相关 >0.9：给二维置信区间，不给点值。
- 正式 CN 行必须包含 §4.4 的全部 provenance 字段；缺一列即验收失败。

**缺陷模型**

只接受结构化 JSON，不接受自由文本。P3 至少实现：

```json
{"type": "nanoparticle", "diameter_nm": 2.5, "shape": "spherical", "lattice": "fcc"}
{"type": "vacancy", "site": "O1", "fraction": 0.1}
{"type": "custom", "expected_n_ratio": 0.75, "uncertainty": 0.10, "justification": "..."}
```

- `nanoparticle`：用 B17 的球截断公式 `N/N_bulk = 1 − 1.5u + 0.5u³`（`u = r1/d`）计算预期比例；只对 fcc/hcp/bcc 金属开放，其他结构拒绝或要求 `custom`。
- `vacancy`：预期比例 = 1 − fraction，前提是 vacancy 位点在当前 CN 组内；位点映射失败则拒绝。
- `custom`：必须有数值 `expected_n_ratio` 和 `uncertainty`；只标 `cn_model_dependent`，不能进主正式表。
- 拟合 `N/N₀` 与模型预期的差异必须在声明的不确定度内；否则 `defect_model_mismatch`。
- 没有定量预期 → `defect_model_missing`。

**多相**

- 每个候选相独立解析、独立门禁、独立一致性检验。
- 输出 `phase_comparison`，不按 Rf 最低自动选相。
- 只有一相通过：报告该相，其他相写淘汰原因。
- 多相通过：分列，不合并；输出 `phase_discrimination_power`，若只有第一壳可比，明确写“无法唯一指认相”。
- 全部失败：0 CN，输出诊断与拒绝原因。

**表拆分**

- `cn_table.csv`：只放 `verified_source_cif` 的 `crystal_cif` / `crystal_edited`。
- `cn_model_dependent.csv`：`prototype` / `cluster`，标签 `cn_model_dependent` / `cn_hypothesis_only`。
- 报告不得把两张表混成一个“CN”列。

**P3 测试（至少 7 项）**

`test_first_shell_single_scale_factor`、`test_conditional_cn_out_of_range_requires_defect`、`test_conditional_cn_nanoparticle_model`、`test_multiphase_no_auto_pick`、`test_cn_table_provenance_complete`、`test_prototype_cluster_separate`、`test_reproducible_cn_hash`。

**P3 验收**

- `contract_tests=110`、`stress_tests=4`。
- 通过者输出条件 CN；未通过者 0 CN；超范围无定量缺陷模型拒绝。
- 多相不自动选；prototype/cluster 分表；报告显示 CIF 来源、sha256、吸收位、N₀、S₀²、门禁状态。
- LEPIDO 基线与 #60 不变。

### 12.10 P4 执行卡：验证、迁移与文档

**内容**

- 标准样：Fe 箔（bcc）、Co 箔（hcp）、Ni 箔（fcc）、Cu 箔（fcc）；优先补真实来源 CIF，否则只能进 `cn_model_dependent.csv`。
- 真实氧化物：赤铁矿、针铁矿、磁铁矿、CoO、Co₃O₄、Ni(OH)₂、CuO；有真实 CIF 的进正式链。
- 错误相测试：赤铁矿 CIF 拟合针铁矿数据 → `consistency_failed`。
- 离线与可复现：相同输入 + 相同缓存 → 相同输出哈希；哈希不含时间戳/绝对路径。
- 文档：README、技能文档、`EXECUTION_MASTER_PLAN.md` 更新正式 CN 说明与测试数。
- `--legacy-auto-fit` 标 `legacy_unverified`，不进正式表；旧 `--fit-motif` 只保留诊断用途。

**P4 测试（至少 4 项）**

`test_offline_full_chain`、`test_standard_foils_formal_cn`、`test_no_cif_report_empty`、`test_legacy_auto_fit_marked`。

**P4 验收**

- `contract_tests≥114`、`stress_tests=4`。
- 正式 CN 表 100% 行含结构/数据/S₀² provenance。
- 无 CIF 样品报告为空 CN + 原因。
- 现有 74 项回归在诊断模式下保持全绿；LEPIDO 与 #60 不变。

### 12.11 测试账本（最低要求）

| 批次 | 新增测试数 | 累计契约回归 | 关键测试 |
|---|---:|---:|---|
| P0a | 4 | 78 | `test_formal_cn_refusal_without_model`、`test_formal_cn_l1_not_in_table`、`test_occupancy_scale_failure_hard_fail`、`test_cn_table_schema_and_refusal_row` |
| P0b | 4 | 82 | `test_motif_provenance_schema`、`test_generated_cif_not_formal_source`、`test_user_cif_is_verified_source`、`test_b20_motifs_not_degraded` |
| P1a | 4 | 86 | `test_motif_coverage_audit`、`test_phase_pack_import_and_index`、`test_phase_pack_sha256_mismatch_rejected`、`test_uncovered_element_no_formal_cn` |
| P1b | 7 | 93 | `test_phase_spec_parser_no_eval`、`test_phase_cif_user_source`、`test_phase_id_cache_hit`、`test_phase_id_cache_miss_offline`、`test_phase_name_ambiguous`、`test_phase_index_schema`、`test_cache_sha256_reuse` |
| P2a | 8 | 101 | `test_cif_gate_parse_failed`、`test_cif_gate_absorber_missing`、`test_cif_gate_absorber_site_ambiguous`、`test_cif_gate_occupancy_implicit_ok`、`test_cif_gate_occupancy_disorder_rejected`、`test_cif_gate_cocupancy_unresolved`、`test_cn_range_12_metals`、`test_n0_shell_grouping_rayleigh` |
| P2b | 2 | 103 | `test_fixed_n_nidp_formula`、`test_wrong_phase_consistency_failed` |
| P3 | 7 | 110 | `test_first_shell_single_scale_factor`、`test_conditional_cn_out_of_range_requires_defect`、`test_conditional_cn_nanoparticle_model`、`test_multiphase_no_auto_pick`、`test_cn_table_provenance_complete`、`test_prototype_cluster_separate`、`test_reproducible_cn_hash` |
| P4 | 4 | 114 | `test_offline_full_chain`、`test_standard_foils_formal_cn`、`test_no_cif_report_empty`、`test_legacy_auto_fit_marked` |

如果执行时新增了额外测试，更新本表；不要为了凑数删测试。

### 12.12 已知陷阱（执行前逐条读）

1. **PowerShell 不支持 `&&`**：用分行、`;` 或 `if ($LASTEXITCODE -ne 0) { ... }`。
2. **Python 路径有空格**：统一用 `& $PY ...`，不要直接写 `python`。
3. **unittest 写 stderr**：PowerShell 可能把通过显示成 NativeCommandError；以 `Ran N tests` + `OK` 或 `run_ci.py` 退出码为准。
4. **非 ASCII 输出**：临时 Python 脚本先 `sys.stdout.reconfigure(encoding="utf-8")`。
5. **robocopy 退出码**：0–7 成功，>7 失败。
6. **LEPIDO 拟合表是 UTF-8 BOM**：哈希包含 BOM，不要用会改写编码的工具。
7. **`path_links` 默认值**：`share_deltar=True`、`share_sigma2=True`；N-only 绑定必须显式设 false。
8. **`select_shell_paths` 默认 `merge_tol=0.20` + `3n+1`**：正式固定 N 必须用 `2n+1` 模式；CN 分组用 Rayleigh，不要直接用 0.20 定义配位。
9. **`read_cif_sites` 默认缺失 occupancy=1.0**：正式门禁必须直接检查 CIF 标签，不能靠这个默认值判断“缺失”。
10. **`read_cif_sites` 会自动按 occ×多重度选吸收位**：正式链必须显式传吸收位/分数坐标，多晶位未指定直接拒绝。
11. **`fit_summary` 会序列化 `fit_out` 中除 `structure`/`fit` 外的字段**：`formal_cn` 块必须是纯 dict/list/str/float，不能塞 dataclass 对象。
12. **默认 HTML/CSV 输出零退化**：正式段只在 `--formal-cn` 或存在 `cn_table.csv` 时出现；不要改默认 7 列表。
13. **CI 计数来自 `Ran N tests`**：不要 skip 测试；skip 会让 `contract_tests` 与账本不一致。
14. **工作区有未提交改动**：禁止破坏性 git 命令；在副本里做基线提交。
15. **不要用 `eval`/`exec` 解析 `--phase`、`--defect-model`**：只用 `json.loads` 和白名单 `split`。
16. **不要为了让测试通过而放宽门禁**：拒绝优先；失败码写进 JSON/报告。

### 12.13 交接报告模板

每批完成后，把下面内容写入 `C:\Users\lizha\Desktop\CN_CIF_执行方案\ACCEPTANCE_LOG.md`：

```text
批次：P0a
日期：
执行目录：
改动文件：
新增接口：
新增/修改测试：
  - Ran N tests — OK
  - ci_report.json: ok=true, contract_tests=78, stress_tests=4
LEPIDO 基线：
  - r_factor=
  - N=
  - R=
  - fit_table SHA256=
#60 静态 s02：
  - OK / 失败详情
正式 CN 拒绝路径：
  - 命令：
  - cn_table.csv 行数 / gate_status / reject_reason
回滚点：
  - commit / patch 路径
遗留问题：
  - ...
```

---

## 13. motif 库稀疏性与覆盖策略（2026-11-08 补充）

> 本节回应“motif 库是否太稀少”的担忧。结论：**确实稀疏，但它不是正式 CN 的数据源；正式覆盖由 CIF / phase pack 决定。**

### 13.1 审计结果

对 `C:\Users\lizha\Desktop\Code-experience\_structure_lib` 的实测：

| 项目 | 实测 |
|---|---|
| motif 总数 | 32 |
| 覆盖吸收元素 | 9 / 17：Au、Co、Cr、Cu、Fe、Mo、Ni、Pt、Zn |
| 完全缺失元素 | 8：Ag、Ce、Ir、Mn、Pd、Rh、Ru、U |
| 晶体 motif | 23 |
| 混合物 motif | 2 |
| 团簇 motif | 7 |
| 现场生成 CIF | 23 / 23 个晶体 CIF，339–415 B |
| 带来源编号 | 2：goethite（COD 2211652）、lepidocrocite（COD-9002763） |
| 真实源 CIF | 0 |

关键点：**带来源编号不等于有源文件**。当前这两个 CIF 仍是现场生成的小文件，只能标 `generated_only`，不能作为 `source_cif_sha256` 的来源。

### 13.2 为什么稀疏不阻塞正式 CN

正式 CN 链不消费 motif，而消费 CIF：

1. `--phase-cif`：用户 CIF，覆盖无上限；
2. `--phase-id`：精确编号 + 本地缓存，只取实际用到的条目；
3. `phase_pack`：本地 curated CIF 包，可增量扩展；
4. `prototype` / `cluster`：显式模型，单独表。

motif 只影响候选排序与诊断。稀疏最多导致“没有候选”，不会导致“用错结构出正式 CN”。

### 13.3 覆盖分层

| 层 | 来源 | 可进正式 CN | 数量策略 |
|---|---|---|---|
| A | 用户 CIF（`--phase-cif`） | 是 | 无上限，用户负责 |
| B | 精确数据库 CIF（`--phase-id` + cache） | 是 | 只缓存用到的条目 |
| C | phase pack（本地 curated） | 是 | 小而有据，首批 15–25 条 |
| D | prototype（生成 CIF） | 否，`cn_model_dependent` | 按需 |
| E | cluster（SAC/DAC/alloy） | 否，`cn_hypothesis_only` | 按需 |
| F | motif only（无 CIF） | 否 | 仅诊断 |

### 13.4 phase pack 规范

- 目录：`_structure_lib/_phase_pack/<source>/<id>/<sha256>.cif` + `.meta.json`。
- 索引：`_structure_lib/_phase_index.json`，由 pack 与 cache 生成。
- 每条必须含：`source`、`id`、`formula`、`spacegroup`、`cell`、`absorber_sites`、`cif_sha256`、`experimental`、`model_level`。
- 导入：`tools/import_phase_pack.py`，只读本地目录，不下载。
- 幂等：同一 `(source, id, sha256)` 重复导入不重复写。
- 冲突：同 ID 不同 sha256 必须报错，除非显式 `--replace`。
- 校验：CIF 可解析、sha256 匹配、吸收元素/位点存在；任一失败拒绝导入。

schema 示例见 `phase_pack_manifest.schema.json`。

### 13.5 覆盖 KPI

- 5 个已验证元素 Fe/Co/Ni/Cu/Pt：每种至少 1 条 A/B/C 级 `verified_source_cif`。
- 其余 12 个吸收元素：允许 `not_covered`，但 `phase_coverage.json` 必须明确列出。
- 未覆盖元素要出正式 CN，只能走 `--phase-cif`；否则 0 行 + `no_cif`。
- generated CIF 只能计入 `prototype`，不得计入 verified 覆盖。
- `phase_coverage.json` 每个元素至少包含：`verified_cif`、`prototype`、`cluster`、`formal_cn_ready`、`missing_common_phases`。

### 13.6 种子包建议

- 首批目标 15–25 条，优先 5 个已验证元素与 B17–B20 涉及的相。
- 来源优先级：用户 XRD 精修 CIF > COD/AMCSD 单条导出（开发期人工导入）> 原型结构。
- 运行时保持离线；不在第一版做在线检索。
- 不把现有 23 条生成 CIF 伪装成真实来源。
- 暂时没有真实 CIF 时：先出 prototype 表，正式 CN 留空；不要把 prototype 数字写进正式表。
- 可选的“诊断覆盖”扩展：为缺失元素补 prototype motif，只提升候选排序，不改变正式 CN 规则。

### 13.7 与批次的关系

- P0b：motif provenance 分类，生成 `motif_coverage.json`。
- P1a：phase pack + coverage report；新增 4 项测试。
- P1b：resolver/cache；新增 7 项测试。
- P2a/P2b：CIF 门禁/一致性；未覆盖元素继续 `no_cif`。
- P3：正式 CN 只从 verified CIF 出。
- P4：补种子包、标准样、文档。

### 13.8 验收底线

- `motif_coverage.json` 覆盖 32/32，无“未分类”。
- `phase_coverage.json` 覆盖 17/17 吸收元素，明确 verified / prototype / cluster / not_covered。
- Fe/Co/Ni/Cu/Pt 每种至少 1 条 verified CIF。
- 未覆盖元素正式 CN 为 0 行 + `no_cif`。
- 任何正式 CN 行都能追到真实 CIF + sha256 + 吸收位，不依赖 motif。

---

*本文档只做评审、设计与执行交接，不包含仓库代码改动。后续按 P0a → P0b → P1a → P1b → P2a → P2b → P3 → P4 分批推进，每批独立回滚、独立验收。*
