# T5 S₀² 标定报告（独立产物）

- 生成时间：2026-11-09T12:18:56
- 静态表：`absorber.py` **未改动**（#60：12 元素 s02=None）
- 运行时使用：`run_xafs.py --s02 <值>`
- 协议：larch 归一化 + 生产等效拟合（N 固定、S₀²/ΔR/σ²/ΔE₀，kw=2，3 窗口鲁棒性）

通过 6 / 12，未通过 6。

| 元素 | 状态 | S₀² | err | Rf | spread | motif | 标准谱 | 未通过原因 |
|------|------|-----|-----|----|--------|-------|--------|------------|
| Cr | calibrated | 0.6482 | 0.0198 | 0.0074 | 0.0053 | cr2o3_corundum | data\bmm_standards\Data\Cr\Cr-K-Cr2O3.xdi | — |
| Mn | calibrated | 0.4311 | 0.0313 | 0.0402 | 0.0213 | mno2_rutile | data\bmm_standards\Data\Mn\Mn-Pyrolucite.xdi | — |
| Zn | calibrated | 0.6071 | 0.0213 | 0.0090 | 0.0070 | zno_wurtzite | data\bmm_standards\Data\Zn\Zn-Zincite.xdi | — |
| Mo | not_calibrated | — | — | — | — | — | — | 门禁未通过: window_runs_interior, s02_in_range, s02_err_ok, window_spread_ok |
| Ru | calibrated | 0.7002 | 0.0275 | 0.0154 | 0.0123 | ruo2_rutile | data\bmm_standards\Data\Ru\Ru-K-RuO2.xdi | — |
| Rh | not_calibrated | — | — | — | — | — | — | 标准谱过短：k_max 不足（IIT-Lytle 仅约 170–195 eV above E0），S5/S6 WT 无有效网格 |
| Pd | not_calibrated | — | — | — | — | — | — | 标准谱过短：k_max 不足（IIT-Lytle 仅约 170–195 eV above E0），S5/S6 WT 无有效网格 |
| Ag | calibrated | 0.7513 | 0.0115 | 0.0030 | 0.0141 | ag_fcc_cod_9008459 | _standards_lib\Ag\Ag_metal_foil.dat | — |
| Ir | not_calibrated | — | — | — | — | — | — | 标准谱过短：k_max 不足（IIT-Lytle 仅约 170–195 eV above E0），S5/S6 WT 无有效网格 |
| Au | calibrated | 0.8602 | 0.0060 | 0.0003 | 0.0073 | au_fcc | data\bmm_standards\Data\Au\Au-L3-GoldFoil.xdi | — |
| Ce | not_calibrated | — | — | — | — | — | — | 门禁未通过: window_runs_interior, window_spread_ok |
| U | not_calibrated | — | — | — | — | — | — | 门禁未通过: window_runs_interior, s02_in_range, s02_err_ok, window_spread_ok |

## 诚实边界

- 未通过元素不产生 S₀²；不以放宽门禁/手工调窗口凑数。
- `_standards_lib/` 与 `data/standards_raw/` 为 IIT-Lytle 标准，许可未明确，**不提交**；相关标定结果记录源文件 SHA256，但不保证在无本地标准的机器上可复现。
- Ce 在 rmin=1.6 下前两窗口 ΔE₀>10 eV，三窗口门禁未通过；该结果按失败记录。
- Mo 的 Na₂MoO₄ 四面体标准给出 s02≈0.38（低于 0.4 下界），Tugarinovite 的 s02_err>0.06，MoS₂ 仅单窗口内部点，Mo 箔撞 1.3 上界——Mo 记为未标定。
- Rh/Pd/Ir/U 的 IIT-Lytle 谱仅约 170–195 eV above E0，归一化平台/独立点数不足，记为未标定。
