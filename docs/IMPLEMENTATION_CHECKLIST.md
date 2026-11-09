# CN 拟合结构溯源实施清单

> 配套文档：`EXECUTION_PLAN.md`  
> 用法：每完成一项勾选；每个批次结束必须跑对应验收用例，未通过不得进入下一批。

---

## P0：止血

### 代码

- [ ] `run_xafs.py`：`--fit` 增加结构来源检查；无 `--phase-cif` / `--phase-id` / `--phase` / `--build-model` 时拒绝 CN。
- [ ] `run_xafs.py`：新增 `--legacy-auto-fit`，仅回归对照。
- [ ] `xafs_fit.py`：移除 `confident`/`--allow-ambiguous` 自动取 Top-1 出 CN 的路径。
- [ ] `xafs_fit.py`：`fit_identified` 必须接收结构模型；无模型时抛 `no_cif`。
- [ ] `xafs_fit.py`：`vary_n` 默认关闭；仅 `--cn-mode conditional` 且一致性通过后放开第一壳。
- [ ] `structure_id.py`：`identify_result` 不再自动调用 L1。
- [ ] `structure_id.py`：L1 改为显式 `l1_first_shell_diagnostic`，输出不进入 CN 表。
- [ ] `xafs_fit.py:216-225`：占据度缩放失败改为硬报错，不允许按满占据继续。
- [ ] `xafs_fit.py:812-815`：`--allow-ambiguous` 不再放行 `no_match`。
- [ ] 输出新增 `structure_source`；`none` 时 CN 列为空。
- [ ] `xafs_report.py`：无 CN 时显示拒绝原因，不显示空误差棒。

### 验收

- [ ] 无 CIF：`cn_table.csv` 0 行 CN，原因 `no_cif`。
- [ ] `--identify`：只输出候选排序，JSON 无 N。
- [ ] L1 诊断：输出存在，但不进入 `cn_table.csv`。
- [ ] `--legacy-auto-fit`：输出 `legacy_unverified=true`，正式报告排除。
- [ ] 现有回归测试按预期更新并通过。

---

## P1：相解析与缓存

### 代码

- [ ] 新增 `phase_resolver.py`。
- [ ] 支持 `--phase-cif`。
- [ ] 支持 `--phase-id SOURCE:ID`（COD/AMCSD）。
- [ ] 支持 `--phase "name=...;formula=...;sg=...;cell=..."`。
- [ ] 硬过滤：化学式含吸收元素、与 XRD 组分一致、空间群/晶胞匹配。
- [ ] 唯一性：恰好一条才自动采用；多条输出候选表并停止。
- [ ] 候选表输出 `phase_resolution.json`。
- [ ] 缓存：`structures/raw/<SOURCE>-<ID>_<sha8>.cif`。
- [ ] 元数据：URL、检索日期、来源、编号、化学式、空间群、晶胞、sha256、实验/计算。
- [ ] 离线：缓存命中可用；未命中 `offline_miss`，不退回 motif。
- [ ] `computed` 结构标记；默认不作为实验 N₀。

### 验收

- [ ] 精确编号只取一条。
- [ ] 同名多条目返回候选表，0 CN。
- [ ] 名称+化学式+空间群唯一时进入门禁。
- [ ] 离线未命中拒绝。
- [ ] 缓存重复使用不重复下载。
- [ ] 不本地化全库；缓存只随候选相增长。

---

## P2：CIF 门禁与固定 N 一致性

### 代码

- [ ] 新增 `cif_gate.py`。
- [ ] 新增 `structure_model.py` 的 Level A：CIF 导入 + 位点编辑。
- [ ] 门禁：解析、吸收元素、吸收位、配位多面体、键长、占据度、共占位、计算结构。
- [ ] 多晶位未指定：`absorber_site_ambiguous`，不出 CN。
- [ ] 缺占据度：`occupancy_missing`，不出 CN。
- [ ] 共占位未拆分：`cocupancy_unresolved`，不出 CN。
- [ ] N₀ 计算：Σ(简并度 × 占据度)。
- [ ] 固定 N=N₀ 一致性拟合：只放 ΔR、σ²、ΔE₀。
- [ ] 输出 `cif_gate.json`、`consistency_fit.json`。
- [ ] 失败码字典与输出。
- [ ] `xafs_fit.fit_identified` 接入 `structure_model`。

### 验收

- [ ] 错误相：`consistency_failed`，0 CN。
- [ ] 正确相：一致性通过，输出 Rf/ΔR/σ²/ΔE₀。
- [ ] 缺占据度：门禁失败。
- [ ] 多晶位未指定：门禁失败。
- [ ] 计算结构：降级或拒绝。
- [ ] CIF 解析失败：`parse_failed`。

---

## P3：条件 CN 与多相报告

### 代码

- [ ] 新增 `cn_conditional.py`（或并入 `xafs_fit.py`）。
- [ ] `--cn-mode consistency-only|conditional`；默认 `consistency-only`。
- [ ] `conditional` 仅在第一壳放开 N；远壳 N 固定 N₀。
- [ ] 输出 `cn_conditional.json`、`cn_table.csv`。
- [ ] 字段：`phase_id`、`model_level`、`source_ref`、`cif_sha256`、`absorber_site`、`shell_label`、`scatterer`、`r_cif`、`n0`、`n_fit`、`n_ratio`、`n_err`、`delta_r`、`sigma2`、`rf`、`cn_label`、`flags`、`assumptions`。
- [ ] `N/N₀` 超范围且无缺陷模型：拒绝。
- [ ] `N_err > 25% N`：不给点值。
- [ ] `|corr(N, σ²)| > 0.9`：给区间，不给点值。
- [ ] 多相：独立检验，输出对照表，不自动选。
- [ ] `prototype` / `cluster` 单独表，不与 CIF CN 混排。
- [ ] `xafs_report.py`：结构来源面板 + CN 表只读 `cn_table.csv`。

### 验收

- [ ] 通过者输出条件 CN；未通过者 0 CN。
- [ ] `N/N₀` 超范围且无缺陷声明：拒绝。
- [ ] 多相通过：分列，不自动选。
- [ ] 团簇模型：`cn_hypothesis_only`。
- [ ] 原型模型：`cn_model_dependent`。
- [ ] 报告可显示 CIF 来源、sha256、吸收位、N₀、门禁状态。

---

## P4：验证与迁移

### 代码/文档

- [ ] 新增第 12 节全部回归用例。
- [ ] 标准样：Fe 箔、赤铁矿、针铁矿、磁铁矿。
- [ ] 错误相：赤铁矿 CIF 拟合针铁矿数据。
- [ ] 离线与缓存测试。
- [ ] 可复现性测试：相同输入+缓存，输出哈希一致。
- [ ] 更新 `README.md`。
- [ ] 更新 `.trae/skills/xafs-pipeline-report/SKILL.md`。
- [ ] 更新 `pipeline-reference.md`、`structure-identification.md`。
- [ ] 标记 `--legacy-auto-fit` deprecated 或移除。
- [ ] 更新拒绝码字典。

### 验收

- [ ] 无 CIF/歧义/门禁失败/一致性失败：0 CN。
- [ ] 正式 CN 表 100% 含 `source_ref + sha256 + absorber_site + N0 + model_level`。
- [ ] 多相 0 次自动 Top-1。
- [ ] 离线缓存命中/未命中行为正确。
- [ ] 报告无 CN 时给出原因。
- [ ] 文档与技能说明与新流程一致。

---

## 测试矩阵

| 用例 | P0 | P1 | P2 | P3 | P4 |
|---|---|---|---|---|---|
| 无 CIF 拒绝 CN | ✓ |  |  |  | ✓ |
| `--identify` 无 N | ✓ |  |  |  | ✓ |
| L1 不进 CN 表 | ✓ |  |  |  | ✓ |
| legacy 标记 | ✓ |  |  |  | ✓ |
| 编号唯一取 CIF |  | ✓ |  |  | ✓ |
| 同名多条目候选表 |  | ✓ |  |  | ✓ |
| 离线未命中拒绝 |  | ✓ |  |  | ✓ |
| 缺占据度门禁 |  |  | ✓ |  | ✓ |
| 多晶位未指定 |  |  | ✓ |  | ✓ |
| 错误相一致性失败 |  |  | ✓ |  | ✓ |
| 固定 N 通过 |  |  | ✓ |  | ✓ |
| 条件 CN 输出 |  |  |  | ✓ | ✓ |
| N/N₀ 超范围拒绝 |  |  |  | ✓ | ✓ |
| 多相不自动选 |  |  |  | ✓ | ✓ |
| 团簇/原型标签 |  |  |  | ✓ | ✓ |
| 可复现哈希 |  |  |  |  | ✓ |

---

## Definition of Done

一个批次完成必须同时满足：

- [ ] 代码改动可独立回滚。
- [ ] 新增/更新测试通过。
- [ ] 输出 JSON/CSV 字段有文档。
- [ ] 失败路径有明确拒绝码。
- [ ] 报告不出现无来源 CN。
- [ ] 旧接口有迁移说明或过渡开关。
- [ ] 不引入全库镜像、不自动 Top-1、不让 L1 进入正式 CN 表。

---

## 落地建议

1. 先做 P0，哪怕只改门控和输出，也能立刻停止虚构 CN。
2. 再做 P2，先把“有结构也不乱出 CN”锁住。
3. P1 解决你“XRD 只给大概组分”的输入问题。
4. P3 才放开条件 CN，避免过早追求数字。
5. P4 用标准样和错误相把边界钉死。
