# CN 拟合结构溯源执行方案

> 版本：v1.0  
> 日期：2026-11-07  
> 适用仓库：`C:\Users\lizha\.local\share\opencode\worktree\d5e1a5\brave-sailor`  
> 目标：把 CN 从“模型输出”改为“有 CIF/结构模型依据的条件测量”；没有合格结构依据时，拒绝输出 CN。

---

## 0. 一页摘要

### 0.1 核心转变

当前逻辑是：

```
--fit → 结构识别 → 取 Top-1 motif 或 L1 合成 O/oct 团簇 → 放开 N → 输出 CN
```

新逻辑是：

```
XRD 候选相 + 用户结构假设
→ 解析到唯一 CIF / 显式结构模型
→ CIF 资格门禁
→ 固定 N=N₀ 的一致性检验
→ 通过后才放开第一壳 N
→ 输出带来源、N₀、N/N₀、假设标签的条件 CN
```

### 0.2 四条铁律

1. **没有合格结构来源，不出 CN。** 谱、FT、小波可以照常出；CN 列空着并写拒绝原因。
2. **XRD 只圈定候选相，不用来定量。** EXAFS 只能在候选集内做排除，不能发现候选集之外的相。
3. **先固定 N 检验结构，再放开 N。** 不允许一上来就对全部壳层 `vary_n=True`。
4. **每个 CN 必须可追溯。** 至少包含 CIF 来源、文件哈希、吸收位、N₀、模型级别、门禁状态。

### 0.3 五个批次

| 批次 | 目标 | 关键结果 |
|---|---|---|
| P0 | 止血 | 禁止无来源 CN；L1 移出自动链；占据度失败硬报错 |
| P1 | 相解析 | 按相名/编号/化学式解析到 0 或 1 份 CIF；本地只缓存用到的条目 |
| P2 | 门禁与一致性 | CIF 资格门禁；固定 N=N₀ 检验；错误相被淘汰 |
| P3 | 条件 CN | 只放开第一壳 N；输出 N₀/N_fit/N/N₀；多相对照 |
| P4 | 验证与迁移 | 标准样、合成数据、离线、回归；文档与旧接口迁移 |

### 0.4 最终效果

- 无 CIF、CIF 歧义、门禁失败、一致性失败：**0 个 CN**。
- 有合格结构：CN 带 `source_ref + sha256 + absorber_site + N0 + model_level`。
- 同名多 CIF：输出候选表，不自动取第一条。
- 错误晶相：在固定 N 一致性阶段被淘汰。
- 非晶/纳米/单原子：用显式团簇模型，CN 标为 `hypothesis_only`，不与晶体 CN 混表。
- 不需要本地镜像 COD/ICSD；缓存只随实际用到的相增长。

---

## 1. 现状与问题

### 1.1 当前链路

1. `run_xafs.py --fit` 调用 `xafs_fit.run()`。
2. `xafs_fit.run()` 调用 `structure_id.identify_result()` 做库内指纹匹配。
3. `xafs_fit.fit_identified()` 在 `confident` 或 `--allow-ambiguous` 时取 `ranked[0]` 做多壳层拟合。
4. 非 `confident` 时，`structure_id.identify_result()` 自动调用 `l1_first_shell_fit()`，默认 O/oct 合成团簇。
5. 拟合时 `vary_n=True`，所有壳层的 N 都作为自由变量。

### 1.2 四个虚构路径

| 路径 | 位置 | 问题 |
|---|---|---|
| 库内 Top-1 直接进 CN | `xafs_fit.py:812-815` | 库外相可能被误判 `confident`，整表 CN 绑错结构 |
| L1 合成 O/oct 团簇 | `structure_id.py:378-393`、`456-565` | 配体/几何是程序假设，N 仍被放开，可产出“可靠但化学错误”的 CN |
| 现场重建 CIF 当结构来源 | `structure_lib.py:115-154` | 用 motif JSON 拼 CIF，运行时“有 CIF”不等于用户提供了结构 |
| 占据度缩放失败被吞 | `xafs_fit.py:216-225` | CIF 缺失或骨架不足时按满占据算 N，N 被高估 1/occ |

### 1.3 本方案要解决什么

- 用户只有 XRD 的“大概组分/可能相”，没有精修 CIF。
- 用户希望按自己提出的结构去匹配 CIF，或构建结构文件。
- CN 拟合不能因为“需要一个模型”而虚构模型。
- 结构来源、模型级别、N₀、拟合偏差必须能被复核。

---

## 2. 设计原则

### 2.1 数据与假设分离

- 数据：E0、χ(k)、FT、小波、XRD 峰位。
- 假设：相名、化学式、吸收位、配位几何、缺陷模型、CIF。
- CN 是“数据 × 假设”的条件量，不是数据本身。

### 2.2 结构来源分级

| 级别 | 来源 | 允许的 CN 标签 |
|---|---|---|
| `crystal_cif` | 用户 CIF 或数据库实验 CIF | `cn_conditional` |
| `crystal_edited` | CIF + 显式位点编辑 | `cn_conditional`，附编辑记录 |
| `prototype` | 原型结构 + XRD/文献晶胞 | `cn_model_dependent` |
| `cluster` | 显式团簇/位点模型 | `cn_hypothesis_only` |
| `none` | 无结构 | 不出 CN |

### 2.3 CIF 资格门禁

“有 .cif 文件”不等于“可用于 CN”。必须通过解析、吸收位、配位多面体、占据度、元素一致性等门禁。

### 2.4 固定 N 优先

先锁 N=N₀，只放 ΔR、σ²、ΔE₀，回答“这份结构能不能解释这张谱”。通过后才允许第一壳 N 浮动。

### 2.5 拒绝优先

程序可以拒绝出数，不可以编数。拒绝原因必须写入 JSON 和报告。

### 2.6 可复现

- CIF 原始文件 + sha256 + 来源 URL/编号 + 检索日期。
- 结构编辑记录。
- 相同输入 + 相同缓存 → 相同输出。

### 2.7 不做全库本地化

- 不镜像 COD/ICSD。
- 只缓存本次样品实际用到的候选 CIF。
- ICSD 由用户导出单条 CIF，走 `user_cif` 通道。

---

## 3. 目标架构

### 3.1 模块图

```
run_xafs.py
  ├─ phase_resolver.py      相名/编号/化学式 → 候选 CIF（0 或 1）
  ├─ cif_gate.py            CIF 资格门禁
  ├─ structure_model.py     结构模型：CIF / 编辑 / 原型 / 团簇
  ├─ structure_id.py        仅候选排序，不产 CN
  ├─ xafs_fit.py            固定 N 一致性 + 条件 CN
  └─ xafs_report.py         只读 cn_table.csv，附结构来源面板
```

### 3.2 数据流

```
用户输入：XRD 候选相 + 结构假设卡
  → phase_resolution.json（候选表 / 唯一选择 / 拒绝原因）
  → cif_gate.json（逐项门禁）
  → structure_model.json + model.cif / model.xyz
  → consistency_fit.json（固定 N=N₀）
  → cn_conditional.json（仅通过者）
  → cn_table.csv（报告唯一 CN 来源）
```

### 3.3 输出目录

```
<outdir>/
  structures/
    raw/
      COD-2211652.cif
      COD-2211652.meta.json
    models/
      hematite_Fe1.cif
      hematite_Fe1.model.json
    phase_resolution.json
    cif_gate.json
    consistency_fit.json
    cn_conditional.json
  cn_table.csv
  report.md / report.html
```

---

## 4. 输入契约

### 4.1 CLI 参数

| 参数 | 说明 |
|---|---|
| `--phase-cif PATH` | 用户提供的 CIF，可重复 |
| `--phase-id SOURCE:ID` | 精确编号，如 `COD:2211652`、`AMCSD:0001234` |
| `--phase SPEC` | 相名/化学式/空间群/晶胞，如 `name=hematite;formula=Fe2O3;sg=R-3c;cell=a=5.0356,c=13.7489` |
| `--absorber-site SITE` | 多吸收位时必须指定，如 `Fe1` |
| `--build-model none\|prototype\|cluster` | 没有 CIF 时是否允许构建模型；默认 `none` |
| `--prototype TYPE` | `rocksalt`、`spinel`、`corundum`、`rutile`、`perovskite`、`layered` 等 |
| `--lattice SPEC` | 原型模型的晶胞参数 |
| `--cn-mode consistency-only\|conditional` | 默认 `consistency-only`；只有 `conditional` 才允许放开第一壳 N |
| `--legacy-auto-fit` | 仅用于回归对照，输出标 `legacy_unverified`，不得进入正式 CN 表 |
| `--identify` | 保留，仅候选排序；JSON 中不得出现 N |

### 4.2 结构假设卡

建议用 YAML/JSON 保存到样品目录，作为可复核输入：

```yaml
phase: hematite
formula: Fe2O3
absorber: Fe
site: Fe1
coordination: octahedral
expected_shells:
  - {el: O,  r: 1.95, n: 3}
  - {el: O,  r: 2.11, n: 3}
  - {el: Fe, r: 2.97, n: 1}
evidence: XRD peak set / literature / DFT
defects: none
notes: ""
```

字段说明：

- `phase` / `formula`：XRD 或文献给出的候选相。
- `absorber` / `site`：吸收元素与晶位；多晶位必须指定。
- `coordination`：预期配位几何，用于门禁和原型构建。
- `expected_shells`：可留空；有则用于局域指纹比对。
- `evidence`：结构假设来源，必须是独立于 EXAFS 拟合的证据。
- `defects`：空位、掺杂、表面、纳米等；没有就写 `none`。

### 4.3 输入优先级

1. `--phase-cif`：用户精修 CIF，最高优先级。
2. `--phase-id`：精确数据库编号。
3. `--phase`：名称/化学式/空间群/晶胞，解析后必须唯一。
4. `--build-model prototype|cluster`：显式授权构建模型。
5. 都没有：`no_cif`，不出 CN。

### 4.4 与现有参数的关系

- `--fit-motif`：废弃为旧接口；仅 `--legacy-auto-fit` 下可用。
- `--allow-ambiguous`：不再放行 `no_match`；只影响候选排序展示。
- `--fit-weights`：保留给显式 mixture 相分数；不产生 CN。
- `--hetero-cocupancy`：保留，但必须基于合格 CIF。
- `--catalytic-model`：保留为显式团簇模型；CN 标签 `cn_hypothesis_only`。

---

## 5. 相解析与 CIF 获取

### 5.1 来源

| 来源 | 用法 | 备注 |
|---|---|---|
| 用户 CIF | 直接使用 | 最优先；同批 XRD Rietveld 结果最佳 |
| COD | 按编号或检索 | 开放；只取单条 |
| AMCSD | 按编号或矿物名 | 矿物学补充 |
| ICSD | 用户导出单条 CIF | 不镜像、不批量抓取 |
| 原型结构 | 用户显式 `--build-model prototype` | 标 `prototype` |
| 团簇模型 | 用户显式 `--build-model cluster` | 标 `cluster` |

### 5.2 硬过滤

候选 CIF 必须同时满足：

1. 化学式含吸收元素，且与 XRD 推断的组分一致。
2. 若 XRD 给出空间群/晶胞，必须在容差内相符。
3. 吸收元素存在；多晶位可由 `--absorber-site` 唯一指认。
4. 条目类型为实验结构；计算/预测结构标 `computed`，不得作为实验 N₀。
5. 占据度字段完整；缺失则门禁失败。

### 5.3 唯一性规则

- 有 `--phase-id`：按编号取，不搜索。
- 只有 `--phase`：检索后必须恰好一条通过硬过滤，才自动采用。
- 多条通过：输出候选表，停止，不出 CN，等待用户指定编号。
- 零条通过：输出 `no_cif` 或 `cif_ambiguous`，不退回 motif 库。

### 5.4 候选表字段

`phase_resolution.json` 至少包含：

```json
{
  "input": {"phase": "hematite", "formula": "Fe2O3", "absorber": "Fe"},
  "status": "unique|ambiguous|not_found|offline_miss",
  "candidates": [
    {
      "source": "COD",
      "ref": "COD:2211652",
      "formula": "Fe2O3",
      "spacegroup": "R-3c",
      "cell": {"a": 5.0356, "b": 5.0356, "c": 13.7489},
      "experimental": true,
      "absorber_sites": ["Fe1"],
      "local_fingerprint": {
        "shells": [
          {"el": "O", "r": 1.95, "n0": 3},
          {"el": "O", "r": 2.11, "n0": 3},
          {"el": "Fe", "r": 2.97, "n0": 1}
        ]
      },
      "match_notes": ["spacegroup matches XRD", "first-shell R within 0.05 Å"],
      "sha256": "...",
      "cache_path": "structures/raw/COD-2211652.cif"
    }
  ],
  "selected": null,
  "reason": "multiple candidates; specify --phase-id"
}
```

### 5.5 缓存

- 目录：`<outdir>/structures/raw/` 或用户配置的共享缓存。
- 文件名：`<SOURCE>-<ID>_<sha8>.cif`。
- 元数据：`<SOURCE>-<ID>_<sha8>.meta.json`。
- 元数据字段：URL、检索日期、来源、编号、化学式、空间群、晶胞、sha256、实验/计算、用户确认人/时间。
- 同一编号重复使用直接读缓存；不重复下载。

### 5.6 离线行为

- 缓存命中：正常使用。
- 缓存未命中：该相 `offline_miss`，不出 CN。
- 不退回 motif 库，不静默使用相似相。

### 5.7 ICSD 与计算结构

- ICSD：不镜像；用户导出单条 CIF，走 `user_cif`。
- Materials Project 等计算结构：标 `computed`；只能用于候选排序或原型参考，不能作为实验 N₀。
- 若用户明确要用计算结构做假设检验，必须 `--allow-computed`，CN 标签降级为 `cn_hypothesis_only`。

---

## 6. CIF 资格门禁

### 6.1 解析

- gemmi 能读取晶胞、空间群、原子位。
- 必须有 `_atom_site_*` 循环；缺失则 `parse_failed`。
- 记录 CIF 原始文件路径与 sha256。

### 6.2 吸收位

- 吸收元素必须存在。
- 同元素多个晶位时，必须 `--absorber-site` 指定。
- 禁止再按 `occ × 多重度` 自动挑主位点来出 CN。
- 未指定多晶位：`absorber_site_ambiguous`，不出 CN。

### 6.3 配位多面体

- 展开后第一壳配位多面体必须闭合。
- 第一壳配位数常见范围 2–8；超出需 `--allow-unusual-coordination` 并写警告。
- 键长必须落在该元素对的化学合理范围；未知元素对需用户提供参考键长。
- 不闭合或键长明显异常：`polyhedron_not_closed` / `bond_length_implausible`。

### 6.4 占据度

- 所有相关位点必须有 `_atom_site_occupancy`。
- 缺失：`occupancy_missing`，不出 CN。
- 共占位：按元素拆分为不同分支；不能合并成单一元素近似。
- 未拆分：`cocupancy_unresolved`，不出 CN。
- 占据度缩放失败：硬失败，不允许按满占据继续。

### 6.5 元素与化学式

- CIF 化学式必须与假设卡/XRD 候选一致。
- 吸收元素价态/配位环境异常时写警告；不自动改结构。
- H/D 默认忽略；必须在输出中记录 `ignored_elements`。

### 6.6 计算结构

- `computed` 结构不得作为实验 N₀。
- 用户显式 `--allow-computed` 时，CN 标签 `cn_hypothesis_only`，并写 `structure_is_computed=true`。

### 6.7 失败码表

| 失败码 | 含义 | 动作 |
|---|---|---|
| `parse_failed` | CIF 解析失败 | 不出 CN |
| `absorber_missing` | CIF 中无吸收元素 | 不出 CN |
| `absorber_site_ambiguous` | 多晶位未指定 | 不出 CN，要求 `--absorber-site` |
| `polyhedron_not_closed` | 第一壳多面体不闭合 | 不出 CN |
| `bond_length_implausible` | 键长异常 | 不出 CN |
| `occupancy_missing` | 缺占据度 | 不出 CN |
| `cocupancy_unresolved` | 共占位未拆分 | 不出 CN |
| `cif_computed` | 计算结构 | 降级或拒绝 |
| `cif_ambiguous` | 多条候选 | 输出候选表，停止 |
| `offline_miss` | 离线且无缓存 | 不出 CN |
| `consistency_failed` | 固定 N 检验失败 | 不出 CN |
| `n_not_identifiable` | 信息量不足/相关过高 | 不出点值 CN |
| `no_cif` | 无结构来源 | 不出 CN |

---

## 7. 结构模型构建

### 7.1 级别 A：CIF 导入与位点编辑

允许的编辑：

- 替位：A/B 位元素替换。
- 空位：删除位点或降低占据度。
- 掺杂：指定晶位放入掺杂元素。
- 超胞/切团簇：为 FEFF 扩到 5–6 Å。
- 共占位拆分：按元素分支。

要求：

- 每次编辑写入 `structure_edit`：操作、目标位点、原值、新值、理由、操作人/时间。
- 编辑后重新计算 N₀。
- 不允许静默改结构。

### 7.2 级别 B：原型结构构建

适用：XRD 能判断结构类型，但没有同相 CIF。

支持原型：

- `rocksalt`、`spinel`、`corundum`、`rutile`、`perovskite`、`layered` 等。
- 晶胞参数来自 XRD 峰位、文献或用户 `--lattice`。
- 吸收位与配位几何由假设卡指定。

输出：

- `model_level=prototype`
- `cn_label=cn_model_dependent`
- N₀ 是理想配位数，不是实验 CN。

### 7.3 级别 C：团簇/位点模型

适用：非晶、纳米、表面、单原子催化剂、界面。

输入：

- 吸收原子在原点；
- 近邻 `(element, r, n, occ)` 显式列出；
- 或从体相 CIF 切团簇后加表面/空位。

输出：

- `model_level=cluster`
- `cn_label=cn_hypothesis_only`
- 必须写清配体、几何、缺陷假设。

### 7.4 编辑记录

```json
{
  "model_id": "hematite_Fe1_vacancy_10pct",
  "base_cif": "COD-2211652",
  "edits": [
    {"op": "scale_occupancy", "site": "O1", "from": 1.0, "to": 0.9,
     "reason": "XRD/文献提示氧空位"}
  ],
  "n0_recomputed": true,
  "approved_by": "user",
  "timestamp": "2026-11-07T12:00:00"
}
```

### 7.5 N₀ 计算

- 晶体：展开后按 `(元素, 距离容差)` 分壳，`N₀ = Σ(简并度 × 占据度)`。
- 团簇：`N₀ = Σ(同壳位点数 × 占据度)`。
- 合并壳层：只有 CIF 本身允许或用户显式 `--merge-shells` 时才合并；必须记录被合并的壳层和 N₀ 求和依据。
- 不允许用“代表路径 + 0.20 Å 容差聚类”重新定义 CIF 的配位。

### 7.6 模型级别与 CN 标签

| model_level | 说明 | CN 标签 | 是否进入正式 CN 表 |
|---|---|---|---|
| `crystal_cif` | 实验 CIF | `cn_conditional` | 是 |
| `crystal_edited` | CIF + 显式编辑 | `cn_conditional` | 是，附编辑记录 |
| `prototype` | 原型结构 | `cn_model_dependent` | 单独表，不与 CIF CN 混排 |
| `cluster` | 团簇/位点模型 | `cn_hypothesis_only` | 单独表，不与 CIF CN 混排 |
| `none` | 无结构 | 无 | 否 |

---

## 8. 拟合流程

### 8.1 Step 1：固定 N=N₀ 的一致性检验

- N 全部锁在 N₀。
- 只放 ΔR、σ²、ΔE₀。
- 输出 Rf、各壳 ΔR/σ²、ΔE₀、N_idp、变量数。
- 不通过：`consistency_failed`，不出 CN。

### 8.2 Step 2：条件 N

只有 Step 1 通过，且 `--cn-mode conditional`：

- 只放开第一壳 N；远壳 N 保持 N₀。
- 输出 `N₀`、`N_fit`、`N/N₀`、`N_err`、`ΔR`、`σ²`。
- `N/N₀` 超出 0.7–1.3 时，必须声明缺陷模型（表面/空位/纳米），否则该壳不进入结论。
- `N_err > 25% N`：不给点值，只给区间或 `n_not_identifiable`。
- `|corr(N, σ²)| > 0.9`：不给点值，给二维置信区间。

### 8.3 Step 3：多相处理

- 每个候选相独立解析、独立门禁、独立一致性检验。
- 输出 `phase_comparison` 表。
- 不按 Rf 最低自动选相。
- 只有一相通过：报告该相条件 CN，其他相写淘汰原因。
- 多相通过：分列 CN，不合并。
- 全部失败：不出 CN，输出诊断。

### 8.4 门禁阈值（初始值，按标准样校准）

| 项目 | 通过 | 警告 | 失败 |
|---|---|---|---|
| Rf（一致性） | ≤ 0.08 | 0.05–0.08 | > 0.08 |
| Rf（CN 报告） | ≤ 0.05 | 0.05–0.08 且 N_err ≤ 20% | > 0.08 |
| ΔR（相对 CIF 壳层） | ≤ 0.10 Å | 0.10–0.15 Å | > 0.15 Å 或撞界 |
| ΔE₀ | < 10 eV | 10–15 eV | > 15 eV |
| σ² | 1e-4–0.03 Å² | 接近边界 | 撞界或负值 |
| N/N₀（无缺陷声明） | 0.7–1.3 | 0.6–0.7 / 1.3–1.4 | < 0.6 或 > 1.4 |
| N_err | ≤ 25% N | 25–40% | > 40% |
| N_idp | 3n+1 < N_idp | 接近上限 | 超限 |

### 8.5 拒绝码

- `no_cif`
- `cif_ambiguous`
- `offline_miss`
- `gate_failed:<code>`
- `consistency_failed`
- `n_not_identifiable`
- `n_out_of_range`
- `correlation_too_high`
- `defect_model_missing`
- `legacy_unverified`

### 8.6 与现有 `fit_r_space` 的关系

- `artemis.fit_r_space` 保留，继续作为数值引擎。
- `xafs_fit.fit_identified` 不再自动决定结构；必须接收 `structure_model`。
- `vary_n` 只在 `cn_mode=conditional` 且 Step 1 通过后，对第一壳启用。
- `fit_weights` 继续用于显式 mixture 相分数，不产生 CN。
- `include_ms` 仅在 CIF 几何支持时启用；MS 路径必须由母壳层约束，不新增自由变量。

---

## 9. 输出契约

### 9.1 `phase_resolution.json`

- 输入假设卡
- 候选表
- 选中条目
- 拒绝原因
- 缓存路径与 sha256

### 9.2 `cif_gate.json`

- 逐项门禁结果
- 失败码
- 吸收位
- 忽略元素
- 占据度摘要
- 第一壳 N₀

### 9.3 `consistency_fit.json`

- 固定 N=N₀ 的拟合结果
- Rf、ΔR、σ²、ΔE₀、N_idp、变量数
- 通过/失败
- 失败原因

### 9.4 `cn_conditional.json`

- 每个通过相的每个壳层：
  - `phase_id`
  - `model_level`
  - `source_ref`
  - `cif_sha256`
  - `absorber_site`
  - `shell_label`
  - `scatterer`
  - `r_cif`
  - `n0`
  - `n_fit`
  - `n_ratio`
  - `n_err`
  - `delta_r`
  - `sigma2`
  - `rf`
  - `flags`
  - `assumptions`
  - `defect_model`

### 9.5 `cn_table.csv`

报告唯一 CN 来源。列：

```
phase_id, model_level, source_ref, cif_sha256, absorber_site,
shell_label, scatterer, r_cif, n0, n_fit, n_ratio, n_err,
delta_r, sigma2, rf, cn_label, gate_status, flags, assumptions
```

无 CN 时输出表头 + 一行拒绝原因，不输出空表。

### 9.6 报告结构来源面板

报告必须显示：

- 结构来源：`user_cif` / `COD:xxxx` / `AMCSD:xxxx` / `prototype` / `cluster`
- CIF sha256 前 8 位
- 吸收位
- 模型级别
- N₀
- 门禁状态
- 拒绝原因（如有）

---

## 10. 代码改动清单

### 10.1 新增模块

| 文件 | 职责 |
|---|---|
| `phase_resolver.py` | 相名/编号/化学式解析；COD/AMCSD 查询；缓存；候选表 |
| `cif_gate.py` | CIF 资格门禁；失败码；N₀ 计算 |
| `structure_model.py` | 结构模型抽象；CIF/编辑/原型/团簇；模型级别 |
| `cn_conditional.py` | 固定 N 一致性 + 条件 CN + 多相对照（也可并入 `xafs_fit.py`） |

### 10.2 修改模块

| 文件 | 改动 |
|---|---|
| `run_xafs.py` | 新增 `--phase-*`、`--absorber-site`、`--build-model`、`--cn-mode`、`--legacy-auto-fit`；编排新链路 |
| `xafs_fit.py` | `fit_identified` 接收 `structure_model`；移除自动 Top-1/L1；`vary_n` 仅条件模式；占据度失败硬报错；输出 provenance |
| `structure_id.py` | `identify_result` 只做候选排序，不调用 L1；L1 改为显式诊断函数；不得写 N |
| `structure_lib.py` | 保留 motif 作为候选索引；`_cif_path` 仅供内部参考；新增 `local_fingerprint` / `find_candidates`；不得作为 CN 结构来源 |
| `xafs_report.py` | 只读 `cn_table.csv`；新增结构来源面板；无 CN 时写拒绝原因；移除从 `fit` 字典直接取 N |
| `artemis.py` | 保留 `fit_r_space`；如需失败码，增加 `CIFError` 细分；不改数值算法 |
| `regression_test.py` | 新增第 12 节用例 |

### 10.3 废弃/降级

- `--fit-motif`：仅 `--legacy-auto-fit` 下可用。
- `--allow-ambiguous`：不再放行 `no_match`；只影响排序展示。
- `structure_lib` 的 motif 自动识别：降级为“候选建议”，不产 CN。
- `xafs_report.py:1308` 的 `lepidocrocite.cif` 直接读取：仅保留为 demo，不进入正式 CN 链。

### 10.4 文档与技能

- `README.md`：更新 `--fit` 行为、无 CIF 行为、CN 标签。
- `.trae/skills/xafs-pipeline-report/SKILL.md`：更新 S7 流程与拒绝规则。
- `.trae/skills/xafs-pipeline-report/references/pipeline-reference.md`：更新 CIF/CN 章节。
- `.trae/skills/xafs-pipeline-report/references/structure-identification.md`：明确 motif 仅候选排序。

---

## 11. 分阶段实施

### P0：止血（1 个 PR）

改动：

- `--fit` 必须有 `--phase-cif` / `--phase-id` / `--phase` / `--build-model`。
- 移除自动 Top-1 全壳 CN。
- L1 移出自动链，改为显式 `l1_first_shell_diagnostic`。
- `--allow-ambiguous` 不再放行 `no_match`。
- 占据度缩放失败硬报错。
- 输出 `structure_source`；`none` 时 CN 列为空。

验收：

- 无 CIF 测试：0 个 CN，原因 `no_cif`。
- `--identify` 输出无 N。
- `--legacy-auto-fit` 输出标 `legacy_unverified`，不进正式表。
- 现有回归测试通过或按预期更新。

回滚：保留 `--legacy-auto-fit` 一个版本。

### P1：相解析与缓存（1–2 个 PR）

改动：

- `phase_resolver.py`。
- `--phase-cif` / `--phase-id` / `--phase`。
- COD/AMCSD 查询与缓存。
- 候选表输出；歧义停止。
- 离线缓存命中/未命中。

验收：

- 精确编号只取一条。
- 同名多条目返回候选表，0 CN。
- 离线未命中拒绝，不退回 motif。
- 缓存元数据含 URL、日期、sha256。

### P2：门禁与一致性（1–2 个 PR）

改动：

- `cif_gate.py`。
- `structure_model.py` Level A。
- 固定 N=N₀ 一致性拟合。
- 失败码与输出 JSON。
- `xafs_fit.fit_identified` 接收模型。

验收：

- 错误相 `consistency_failed`，0 CN。
- 正确相通过，输出一致性指标。
- 缺占据度、多晶位未指定：门禁失败。
- CIF 计算结构：降级或拒绝。

### P3：条件 CN 与多相报告（1–2 个 PR）

改动：

- `cn_conditional.py`。
- `--cn-mode conditional`。
- 第一壳 N 放开；N₀/N_fit/N/N₀。
- 多相对照表。
- `xafs_report.py` 结构来源面板。

验收：

- 通过者输出条件 CN；未通过者 0 CN。
- `N/N₀` 超范围且无缺陷模型：拒绝。
- 多相通过：分列，不自动选。
- `prototype` / `cluster` 单独表。

### P4：验证与迁移（1 个 PR）

改动：

- 标准样与合成数据回归。
- 离线与缓存测试。
- 文档、技能、README。
- 移除 `--legacy-auto-fit`（或标记 deprecated）。

验收：

- 第 12 节全部用例通过。
- 报告可复核每个 CN。
- 无 CIF 样品报告为空 CN + 原因。

---

## 12. 验收测试

### 12.1 必测用例

1. `test_no_phase_refuses_cn`
   - 无 `--phase-*`，`--fit` 返回 0 CN，原因 `no_cif`。
2. `test_phase_id_unique_fetch_and_cache`
   - `--phase-id COD:xxxx` 只取一条，缓存命中不重复下载。
3. `test_phase_name_ambiguous_returns_candidates`
   - 同名多条，返回候选表，0 CN。
4. `test_phase_name_unique_with_spacegroup`
   - 名称+化学式+空间群唯一，进入门禁。
5. `test_cif_gate_missing_occupancy`
   - 缺 `_atom_site_occupancy`，门禁失败，0 CN。
6. `test_cif_gate_multiple_absorber_sites_requires_site`
   - 多晶位未指定，失败，0 CN。
7. `test_cif_gate_cocupancy_unresolved`
   - 共占位未拆分，失败，0 CN。
8. `test_wrong_phase_consistency_failed`
   - 错误相固定 N 检验失败，0 CN。
9. `test_fixed_n_consistency_pass_then_conditional_cn`
   - 正确相通过，输出 N₀/N_fit/N/N₀。
10. `test_conditional_cn_out_of_range_requires_defect_model`
    - N/N₀ 超范围且无缺陷声明，拒绝。
11. `test_multiphase_comparison_no_auto_pick`
    - 两相都通过，分列，不自动选。
12. `test_cluster_model_labeled_hypothesis`
    - 团簇模型 CN 标 `cn_hypothesis_only`，不进正式表。
13. `test_prototype_model_labeled_dependent`
    - 原型模型 CN 标 `cn_model_dependent`。
14. `test_offline_cache_hit_and_miss`
    - 命中可用；未命中 `offline_miss`。
15. `test_legacy_auto_fit_marked_unverified`
    - 旧路径输出 `legacy_unverified`，不进正式表。
16. `test_l1_diagnostic_not_in_cn_table`
    - L1 显式诊断不出现在 `cn_table.csv`。
17. `test_occupancy_scale_failure_hard_fail`
    - 占据度缩放失败，硬报错，不按满占据。
18. `test_reproducible_output_hash`
    - 相同输入+缓存，输出 JSON 哈希一致。

### 12.2 标准样与合成数据

- Fe 箔（bcc）：第一壳 N₀=8，固定 N 检验应通过；条件 CN 回收在误差内。
- 赤铁矿、针铁矿、磁铁矿：合格 CIF 下固定 N 检验；已知不可定量壳层保持拒绝。
- 合成磁铁矿 mixture：仅用于相分数，不产 CN。
- 错误相：用赤铁矿 CIF 拟合针铁矿数据，应 `consistency_failed`。
- 缺占据度 CIF：门禁失败。
- 纳米/团簇模型：CN 标 `cn_hypothesis_only`。

### 12.3 验收指标

- 无 CIF/歧义/门禁失败/一致性失败：**0 个 CN**。
- 正式 CN 表：**100%** 行含 `source_ref + sha256 + absorber_site + N0 + model_level`。
- 多相：**0 次**自动 Top-1。
- 离线：缓存命中可用；未命中拒绝。
- 可复现：相同输入输出一致。
- 缓存规模：随候选相数量增长，不随数据库总量增长。

---

## 13. 迁移与兼容

### 13.1 破坏性变更

- `--fit` 不再自动识别并出 CN。
- `--fit` 无结构来源 → 拒绝。
- `--fit-motif` 降级为旧接口。
- `--allow-ambiguous` 不再放行 `no_match`。

### 13.2 旧路径

- `--legacy-auto-fit`：保留一个版本，输出 `legacy_unverified=true`，报告不纳入正式 CN 表。
- 旧 JSON 字段尽量保留，新增字段不破坏读取。
- 旧报告若直接读 `fit` 字典取 N，改为读 `cn_table.csv`；无表则显示拒绝原因。

### 13.3 motif 库

- 保留为候选索引和参考结构。
- 不再作为 CN 结构来源。
- `_structure_lib/<元素>/_cif/*.cif` 仅供内部参考，不进入 `source_ref`。
- 若用户想用某个 motif，必须提供对应 CIF 或显式 `--build-model` 并接受模型级别标签。

### 13.4 文档

- README 增加“无 CIF 时会发生什么”。
- 技能文档更新 S7：先相解析、门禁、固定 N、条件 N。
- 增加拒绝码字典。

---

## 14. 风险与边界

| 风险 | 处理 |
|---|---|
| CIF 质量差 | 门禁 + 固定 N 检验；不合格拒绝 |
| 同名多 CIF | 候选表，人工指定编号 |
| XRD 候选集漏掉真实相 | 明确声明：EXAFS 不能发现候选集之外的相 |
| ICSD 授权 | 不镜像；用户导出单条 CIF |
| 计算结构 | 标 `computed`；不作为实验 N₀ |
| 纳米/缺陷 | N₀ 是理想值；报告 N/N₀ 与缺陷模型 |
| 用户强行指定错误相 | 记录 `user_confirmed`；一致性失败仍拒绝 CN |
| 模型编辑偏置 | 记录 `structure_edit`；重新计算 N₀；报告可见 |
| CN 数量变少 | 设计目标；用拒绝原因替代虚构数字 |
| 旧流程依赖 | `--legacy-auto-fit` 过渡；明确标记不可用于正式结论 |

---

## 15. 效果

### 15.1 对你当前场景

你只根据 XRD 判断“大概是 Fe₂O₃ / FeOOH / NiO”也可以开始：

1. 输入相名、化学式、吸收元素、可能的吸收位。
2. 系统给出候选 CIF 表：编号、空间群、晶胞、吸收位、第一壳 N₀/R、与假设的差异。
3. 没有合适 CIF 时，你可以显式选择 `--build-model prototype` 或 `cluster`，系统构建结构文件并标明模型级别。
4. 每个候选做固定 N 一致性检验。
5. 通过者输出条件 CN：`N₀ / N_fit / N/N₀`；未通过者写原因。

### 15.2 对 CN 质量

- 散射体身份由 CIF/显式模型固定，不再由程序猜。
- N₀ 由 CIF 几何和占据度计算，不再用“八面体=6”这种默认。
- 先固定 N 检验，再放开 N，避免用 N 去补偿错误结构。
- N/N₀ 和缺陷模型让“体相理想配位数”与“样品有效配位数”分开。
- 每个 CN 可追溯、可复现、可复核。

### 15.3 对工程

- 不需要本地镜像 COD/ICSD。
- 只缓存实际用到的候选 CIF。
- 离线可复现。
- 新增相只需加编号或 CIF，不需要重建全库。
- 旧流程有过渡开关，但正式报告只认新链路。

---

## 16. 附：最小落地顺序

如果只能先做一件事：

1. **P0**：让 `--fit` 无 CIF 时拒绝出 CN，并把 L1 移出 CN 表。
2. **P2**：加固定 N 一致性检验。
3. **P1**：再加相解析和缓存。
4. **P3**：最后加条件 N 和多相报告。

理由：P0 立刻停止虚构；P2 保证“有结构也不乱出 CN”；P1/P3 提升可用性和覆盖面。

---

## 17. 附：术语

- **N₀**：结构模型给出的参考配位数（几何简并度 × 占据度），不是实验 CN。
- **N_fit**：固定结构后放开第一壳 N 得到的拟合值。
- **N/N₀**：有效配位比；偏离 1 需要缺陷/表面/纳米模型解释。
- **model_level**：`crystal_cif` / `crystal_edited` / `prototype` / `cluster` / `none`。
- **cn_label**：`cn_conditional` / `cn_model_dependent` / `cn_hypothesis_only`。
- **consistency-only**：只做固定 N 检验，不出 N_fit。
- **conditional**：固定 N 检验通过后，放开第一壳 N。

---

*本方案只定义执行逻辑与验收标准，不包含可直接运行的代码。后续落地时按 P0→P4 分批修改当前仓库，并保持每批可独立回滚。*
