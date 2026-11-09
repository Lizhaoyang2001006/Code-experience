# CN_CIF 执行方案

本目录是后续加入当前项目执行的方案包，不包含可直接运行的代码。

## 文件

- `EXECUTION_PLAN.md`：原主方案。背景、设计原则、目标架构、输入/输出契约、代码改动清单、分阶段实施、验收测试、迁移与风险。
- `IMPLEMENTATION_CHECKLIST.md`：原清单，按 P0–P4 拆分的可勾选实施清单和测试矩阵；已被修订版取代。
- `REVIEW_AND_REVISED_PLAN.md`：**唯一修订方案**。§0–§13 包含评审、双轨设计、接口、算法、测试矩阵、迁移、执行交接手册、motif 覆盖策略。
- `ACCEPTANCE_LOG.md`：基线、每批验收、零退化、回滚点记录。
- `PHASE_PACK_README.md`：phase pack 的目录、导入、索引和覆盖报告说明。
- `phase_pack_manifest.schema.json`：phase pack `manifest.json` 的 JSON Schema。
- `patches/README.md`：批次补丁导出约定。

## 建议阅读顺序

1. `EXECUTION_PLAN.md` 第 0 节：一页摘要。
2. 第 2 节：设计原则。
3. 第 11 节：分阶段实施。
4. `IMPLEMENTATION_CHECKLIST.md`：实际落地时逐项勾选。

## 核心决策

- CN 必须有结构来源：用户 CIF、精确数据库编号，或显式构建的模型。
- XRD 只圈定候选相；EXAFS 只在候选集内做排除。
- 先固定 N=N₀ 做一致性检验，通过后才放开第一壳 N。
- 无 CIF、CIF 歧义、门禁失败、一致性失败：不出 CN。
- 不本地化 COD/ICSD 全库；只缓存实际用到的候选 CIF。
- motif 库降级为候选索引，不再作为 CN 结构来源。
- L1 合成团簇只做显式诊断，不进入正式 CN 表。

## 落地顺序

P0（止血）→ P2（门禁与一致性）→ P1（相解析与缓存）→ P3（条件 CN）→ P4（验证与迁移）。

先停止虚构，再提升覆盖面；不要反过来先追求“每个样品都有 CN”。

## 与当前仓库的关系

- 现有数值引擎 `artemis.fit_r_space` 保留。
- `structure_id.py` 保留候选排序，但不再自动产 CN。
- `structure_lib.py` 保留 motif 作为参考，但 `_cif_path` 生成的 CIF 不作为正式结构来源。
- 新链路从 `run_xafs.py` 编排，新增 `phase_resolver.py`、`cif_gate.py`、`structure_model.py`、`cn_conditional.py`。

## 验收底线

- 正式 CN 表每一行都能追到：CIF 来源 + sha256 + 吸收位 + N₀ + 模型级别。
- 没有合格结构来源时，报告显示空 CN 和拒绝原因，而不是一个看起来合理的数。
