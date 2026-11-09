# Phase pack 使用说明

`phase_pack` 是正式 CN 的**本地真实 CIF 包**。它不是 motif 库，也不是 COD/ICSD 镜像。

## 为什么需要它

当前 motif 库只有 32 条、覆盖 9/17 个吸收元素，且 23 条晶体 motif 的 CIF
都是 339–415 字节的现场生成文件。它能做候选排序和诊断，但不能作为正式
CN 的结构来源。正式 CN 必须绑定真实 CIF + sha256 + 明确吸收位。

phase pack 解决的是“离线、可复现、可增量扩展”的正式来源问题：

- 用户 CIF 永远优先（`--phase-cif`）；
- 精确数据库编号按需缓存（`--phase-id`）；
- phase pack 是本地 curated 小包，首批目标 15–25 条；
- 不做全库镜像，不在运行时联网。

## 目录结构

```text
_structure_lib/_phase_pack/
  manifest.json
  COD/9002763/<sha256>.cif
  COD/9002763/<sha256>.meta.json
  user/my_sample_Fe1/<sha256>.cif
  user/my_sample_Fe1/<sha256>.meta.json
```

`manifest.json` 的 schema 见 `phase_pack_manifest.schema.json`。

## 导入流程

```powershell
& $PY tools/import_phase_pack.py `
  --src "D:\cifs\seed" `
  --pack "_structure_lib\_phase_pack" `
  --source COD `
  --offline
```

要求：

- 每个 CIF 必须有对应 `meta.json`，或由 `--source` + 文件名提供 `source/id`；
- 导入时计算 sha256，与 meta 不一致则拒绝；
- CIF 必须能被 gemmi 解析；
- 吸收元素和 `absorber_sites` 必须存在；
- `absorber_elements` 必须与 `absorber_sites` 的元素一致；
- 同一 `(source, id, sha256)` 重复导入是幂等的；
- 同 ID 不同 sha256 必须显式 `--replace` 才能覆盖。

## 生成索引与覆盖报告

```powershell
& $PY tools/build_phase_index.py `
  --pack "_structure_lib\_phase_pack" `
  --cache "_structure_lib\_phase_cache" `
  --out "_structure_lib\phase_index.json" `
  --coverage "_structure_lib\phase_coverage.json"
```

覆盖报告按吸收元素列出：

- `verified_cif`：真实实验 CIF 数量；
- `prototype`：原型/生成结构数量；
- `cluster`：团簇模型数量；
- `formal_cn_ready`：是否至少有一条 verified CIF；
- `missing_common_phases`：常见但尚未覆盖的相。

## 覆盖底线

- Fe / Co / Ni / Cu / Pt：每种至少 1 条 `verified_source_cif`；
- 其余 12 个吸收元素允许 `not_covered`，但必须明确报告；
- 未覆盖元素要出正式 CN，只能走 `--phase-cif`；
- generated CIF 只能算 prototype，不得计入 verified；
- 任何正式 CN 行都必须能追到真实 CIF + sha256 + 吸收位。

## 首批已导入（P1a）

| 元素 | source:id | 空间群 | 吸收位 | 说明 |
|---|---|---|---|---|
| Fe | COD:1100108 | Im-3m | Fe1 | α-Fe，a=2.8680 Å |
| Co | COD:9008492 | P 63/m m c | Co | hcp Co，a=2.5071/c=4.0686 Å |
| Ni | COD:1534892 | F m -3 m | Ni1 | fcc Ni，a=3.45 Å |
| Cu | COD:4105040 | F m -3 m | Cu | fcc Cu，a=3.5819 Å |
| Pt | COD:1011103 | F m -3 m | Pt1 | fcc Pt，a=3.944 Å |

当前 `phase_coverage.json`：Fe/Co/Ni/Cu/Pt 各 1 条 verified；Cr/Zn/Mo/Au 只有 motif（`motif_only`）；Mn/Ru/Rh/Pd/Ag/Ir/Ce/U 为 `not_covered`。

## 种子包建议

- 优先使用用户自己的 XRD 精修 CIF；
- 其次使用 COD/AMCSD 单条导出，在开发期人工导入；
- 运行时保持离线，不做在线检索；
- 不把现有 23 条生成 CIF 伪装成真实来源；
- 暂时没有真实 CIF 时，正式 CN 留空并写 `no_cif`，不要用 prototype 数字顶替。
