# T5 S₀² 标定报告（独立产物）

- 生成时间：2026-11-09T17:27:37
- 静态表：`absorber.py` **未改动**（#60：12 元素 s02=None）
- 运行时使用：`run_xafs.py --s02 <值>`
- 协议：larch 归一化 + 生产等效拟合（N 固定、S₀²/ΔR/σ²/ΔE₀，kw=2，3 窗口鲁棒性）

通过 12 / 12，未通过 0。

| 元素 | 状态 | S₀² | err | Rf | spread | motif | 标准谱 | 未通过原因 |
|------|------|-----|-----|----|--------|-------|--------|------------|
| Cr | calibrated | 0.6482 | 0.0198 | 0.0074 | 0.0053 | cr2o3_corundum | data\bmm_standards\Data\Cr\Cr-K-Cr2O3.xdi | — |
| Mn | calibrated | 0.4311 | 0.0313 | 0.0402 | 0.0213 | mno2_rutile | data\bmm_standards\Data\Mn\Mn-Pyrolucite.xdi | — |
| Zn | calibrated | 0.6071 | 0.0213 | 0.0090 | 0.0070 | zno_wurtzite | data\bmm_standards\Data\Zn\Zn-Zincite.xdi | — |
| Mo | calibrated | 0.8322 | 0.0334 | 0.0040 | 0.0017 | mo_bcc | data\bmm_standards\Data\Mo\Mo-K-MolybdenumFoil.xdi | — |
| Ru | calibrated | 0.7002 | 0.0275 | 0.0154 | 0.0123 | ruo2_rutile | data\bmm_standards\Data\Ru\Ru-K-RuO2.xdi | — |
| Rh | calibrated | 0.7022 | 0.0191 | 0.0041 | 0.0037 | rh_fcc_cod_9008482 | _standards_lib\Rh\Rh_metal_foil.dat | — |
| Pd | calibrated | 0.6559 | 0.0141 | 0.0035 | 0.0133 | pd_fcc_cod_9008478 | _standards_lib\Pd\Pd_metal_foil.dat | — |
| Ag | calibrated | 0.7513 | 0.0115 | 0.0030 | 0.0141 | ag_fcc_cod_9008459 | _standards_lib\Ag\Ag_metal_foil.dat | — |
| Ir | calibrated | 0.4752 | 0.0149 | 0.0043 | 0.0068 | ir_fcc_cod_9008470 | _standards_lib\Ir\Ir_metal_foil.dat | — |
| Au | calibrated | 0.8602 | 0.0060 | 0.0003 | 0.0073 | au_fcc | data\bmm_standards\Data\Au\Au-L3-GoldFoil.xdi | — |
| Ce | calibrated | 0.6279 | 0.0268 | 0.0223 | 0.0096 | ceo2_fluorite_cod_1562989 | _standards_lib\Ce\Ce_metal_foil.dat | — |
| U | calibrated | 0.7800 | 0.0229 | 0.0137 | 0.0024 | uo2_fluorite | _standards_lib\U\U_metal_foil.dat | — |

## 诚实边界

- 未通过元素不产生 S₀²；不以放宽门禁/手工调窗口凑数。P5-S02 后 12/12 通过，门禁阈值零放宽。
- `_standards_lib/` 与 `data/standards_raw/` 为 IIT-Lytle 标准，许可未明确，**不提交**；相关标定结果记录源文件 SHA256，但不保证在无本地标准的机器上可复现。
- Ce 的 BMM 两条失败谱保留在 attempts 审计链：`Ce-L3-CeO2.xdi`（σ² 贴界、无 interior 窗口）与 `Ce-L3-CeriumDioxide.xdi`（平台 2.271，归一化伪迹）；RAW `Ce/ceo2.004` 在 Ce L2 前截断后通过严格门禁。
- Mo 由 P5-B 双壳 bcc + `half_height` E0 通过；历史失败候选（Na₂MoO₄ 四面体、Tugarinovite、MoS₂）保留在 attempts。
- Rh/Pd/Ir/Ce/U 由 P5-S02 换 IIT-Lytle RAW 长谱后通过；旧短谱证据保留在 `source_pool/quarantined_legacy_source`，门禁未放宽。
