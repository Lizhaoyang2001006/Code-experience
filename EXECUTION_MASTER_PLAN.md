# XAFS 管线最高优先级执行总纲（Athena / Artemis / Larch-WT 全功能对标规范）

> **文档级别**：本项目**唯一最高优先级工程规范与判分标准（Master Execution & Grading Spec）**。
> **用途**：
> 1. 供**高性价比执行模型（小模型 / Execution Model）**按步骤、按文件行号、按数学公式机械无歧义地落地代码，杜绝自由发挥与踩坑；
> 2. 供**旗舰审查模型（Grading Model）**按文末 Checklist、函数签名、物理护栏与数值容差逐条严苛验收。

---

## 0. 全局不可违背铁律（小模型动工前必读）

1. **单向无环依赖铁律**：
   - 归约链：`xafs_core` ← `xafs_io` ← `athena` ← `{artemis, larch_wt, xafs_lcf}` ← `xafs_pipeline` ← `{run_xafs, xafs_report}`
   - 拟合链：`artemis` ← `structure_lib` ← `structure_id` ← `xafs_fit` ← `{run_xafs, xafs_report}`
   - 诊断层：`xafs_pipeline` ← `xafs_diagnostics`
   - **严禁**在 `larch_wt.py` 中顶层 `import artemis`；**严禁**在 `athena.py` 中 `import xafs_pipeline`。
2. **零退化铁律**：
   - 第一至第四阶段（Phase I–IV：B1–B16）已完成 **68 项契约回归 + 4 项压力测试**（全部锁死为基线）。
   - 第五阶段（Phase V：B17–B20，**当前最高优先级执行任务：催化剂精细结构深度解析与电催化/多相催化实战完备化**）每完成一个子项，必须立即运行 `python run_ci.py --stress`，确保既有 68 项回归 + 4 项压力测试 **100% exit=0**，且新增回归全部通过（最终达到 **74 项契约回归**）。
3. **向后兼容铁律**：
   - 所有新增函数参数与 CLI 开关**默认关闭或取向后兼容默认值**；
   - 不传新标志跑默认 `LEPIDO` 三扫描时，必须保持 `r_factor ≈ 0.0273`、第一壳 `N ≈ 6.35`、`R ≈ 2.002 Å` 逐位不变；默认拟合表仍为 7 列。
4. **物理诚实铁律**：
   - **严禁无标定地修改** `absorber.py` 中已标定的 $S_0^2$ 值（Fe=0.803, Co=0.7308, Ni=0.7795, Cu=0.8670, Pt=0.8075）；任何重标定必须通过 `calibrate_s02.py` 门禁并记录数据、命令与窗口极差；
   - Fe=0.803 是氧化物价态验收锚定的跨材料生产值（原始流程不可复原）；2026-11-10 BMM Fe-K-IronFoil.xdi 单壳 bcc 交叉标定 0.6566（IIT-Lytle 5 µm 箔 0.6335）**只适用于金属 Fe**，不替换 0.803（替换会使 Fe 氧化物 CN 高约 22%，违反第 2 条零退化铁律）；金属 Fe 定量用 `--s02 0.6566 --fit-rmax 2.60`；
   - **严禁**为了降低 R-factor 放宽 `cn_reliable` 门禁或插值伪造实验数据点。
5. **Windows PowerShell 环境防坑铁律**：
   - PowerShell 不支持 `&&`（必须用 `;` 或分行），不支持 `python -c` 多行含花括号脚本（一律写临时 `.py` 文件执行后删除）；
   - 打印含 `Å` / `σ²` 等非 ASCII 字符时，脚本入口必须加 `sys.stdout.reconfigure(encoding="utf-8", errors="replace")`。

---

## 第一阶段（Phase I：B1–B4）已验收锁死基线（36 项回归，严禁退化）

> 状态：**100% 已完成并通过旗舰终审**。小模型在执行 Phase II 时不得破坏以下任何一项：

| 批次 | 编号 | 已落地模块 | 核心文件与关键函数 | 对应回归测试（#28–#36） |
|---|---|---|---|---|
| 第一批 | **B1-1** | 多 k-weight `(1,2,3)` 联合最小二乘拟合（分支 RMS $s_w$ 归一，$N_{\text{idp}}$ 不变） | `artemis.py:fit_r_space`, `xafs_fit.py:fit_identified` | `test_multi_kweight_joint_fit` |
| 第一批 | **B1-2** | R 窗反傅里叶变换（Back-FT / $\chi(q)$）滤波与包络诊断 | `athena.py:back_ft_to_q_space`, `xafs_diagnostics.py:write_back_ft_comparison` | `test_back_ft_q_space_phase_and_envelope` |
| 第一批 | **B1-3** | 共线多散射（MS，`nleg ≤ 4`）受母壳层约束纳入 S7（自由变量数不变） | `xafs_fit.py:_collinear_ms_rows`, `artemis.py:fit_r_space(ms_links=...)` | `test_collinear_ms_constrained_fit_nio` |
| 第二批 | **B2-1** | XANES 边后展平（Flattened）+ $E_{1/2}$ / 白线 / 边前峰面积 / 颈部导数双峰 `e0_jump_risk` | `athena.py:normalize_absorption`, `extract_xanes_features` | `test_xanes_flatten_and_neck_features` |
| 第二批 | **B2-2** | XANES / $\chi(k)$ 线性组合拟合（LCF）与 PCA 主成分判定 | `xafs_lcf.py:run_lcf`, `run_pca`, `lcf_from_scans` | `test_lcf_and_pca_synthetic_mixture` |
| 第二批 | **B2-3** | 荧光死时间校正与显式能量区间 Deglitch（删点 >15% 或边核心 <5 点报错） | `athena.py:apply_deadtime_correction`, `deglitch_ranges`, `despike_scan` | `test_deadtime_and_explicit_deglitch` |
| 第三批 | **B3-1** | 壳层 $R_i$ 一维小波包络切片 $|W(k; R_i)|$ + 实验/模型二维 WT 幅值残差 | `larch_wt.py:extract_k_profile`, `compute_model_wt_residual`, `xafs_diagnostics.py:write_wt_shell_profiles` | `test_wt_shell_k_profiles_and_model_residual` |
| 第四批 | **B4-1** | 高 R 噪声估计 $\varepsilon_R$（`rmax_out=25`）、$\chi^2_\nu$ 与增减壳层 Hamilton F-test | `artemis.py:estimate_r_noise`, `shell_f_test` | `test_reduced_chi2_and_shell_ftest` |
| 第四批 | **B4-2** | 三阶累积量 $C_3$ 全链打通 + $N_{\text{idp}}$ 自动降级 `c3_dropped_by_nidp` + `c3_at_bound` 门禁 | `artemis.py:C3_LIMIT`, `xafs_fit.py:_normalize_c3_paths`, `structure_id.py:l1_first_shell_fit` | `test_c3_cumulant_cli_and_nidp_gate` |

---

## 第二阶段（Phase II：B5–B8）最高优先级执行计划（待执行，回归 36 → 45 项）

### 全局执行顺序总览（Phase II）

| 批次 | 编号 | 功能模块 | 对标软件 | 修改文件 | 新增回归测试（`regression_test.py`） |
|---|---|---|---|---|---|
| **第五批** | **B5-1** | 复数小波矩阵 `complex_map`、相位图 $\arg W(R,k)$ 与二维复数域残差 $|W_{\text{data}} - W_{\text{model}}|$ | Larch-WT / EvAX | `xafs_core.py`, `larch_wt.py`, `xafs_diagnostics.py` | `test_wt_complex_residual_catches_phase_inversion` (#37) |
| **第五批** | **B5-2** | 候选散射元素理论小波包络峰位标尺 `reference_element_k_peaks`（含 $\sigma^2$ 调制） | Larch-WT / Hama | `structure_id.py`, `xafs_diagnostics.py` | `test_wt_element_reference_ruler_sigma2_shift` (#38) |
| **第六批** | **B6-1** | 多扫描相异度筛查、异常扫描自动剔除（Outlier Rejection）与标准差传播 $\sigma_\mu(E)$ | Athena | `xafs_core.py`, `athena.py`, `xafs_pipeline.py`, `run_xafs.py`, `xafs_report.py` | `test_merge_scans_outlier_rejection_and_sigma` (#39) |
| **第六批** | **B6-2** | XANES 区域荧光自吸收校正（Tröger/FLUO 能量域反解，与 EXAFS Booth 严格隔离） | Athena | `athena.py`, `xafs_pipeline.py`, `xafs_report.py` | `test_xanes_fluo_self_absorption_restores_whiteline_and_lcf` (#40) |
| **第六批** | **B6-3** | 边前峰（Pre-edge）主边拖尾基线扣除与多 Pseudo-Voigt 分峰拟合（Wilke 质心能量） | Athena | `athena.py`, `xafs_report.py` | `test_pre_edge_voigt_deconvolution_and_centroid` (#41) |
| **第七批** | **B7-1** | 各向同性晶格膨胀约束（`isotropic_deltar=True`，$\Delta R_i = \alpha \cdot R_{\text{eff},i}$） | Artemis | `artemis.py`, `xafs_fit.py`, `run_xafs.py`, `xafs_report.py` | `test_isotropic_lattice_expansion_constraint` (#42) |
| **第七批** | **B7-2** | CIF 部分占据（Partial Occupancy）与空位/共占位简并度自动缩放 | Artemis | `artemis.py`, `structure_lib.py`, `xafs_fit.py` | `test_cif_partial_occupancy_and_vacancy_scaling` (#43) |
| **第七批** | **B7-3** | 四阶累积量 $C_4$ 与非高斯非对称径向分布函数 $g(R)$ 重构 | Artemis | `xafs_core.py`, `artemis.py`, `xafs_fit.py` | `test_c4_cumulant_and_asymmetric_rdf_reconstruction` (#44) |
| **第八批** | **B8-1** | 多数据集/双吸收边/变温爱因斯坦模型全局受约束联合拟合 `fit_multi_datasets` | Artemis | `artemis.py` | `test_multi_dataset_corefinement_and_einstein_model` (#45) |

---

## 第五批（B5）：Larch-WT 复数域与元素散射标尺深度补齐

### B5-1 复数小波矩阵 `complex_map`、相位图与二维复数域残差

#### 1. 物理背景与数学契约
- **当前缺陷**：`WaveletMap` 只保存了模值 `magnitude = |W(R, k)|`；`compute_model_wt_residual` 只计算幅值差 $|W_{\text{data}}| - |W_{\text{model}}|$。若模型在某 $(R, k)$ 处包络高度与实验相同但**相位反了 $\pi$**（如 $\Delta R$ 错半个波长或散射体选反），幅值差 $|W_{\text{data}}| - |W_{\text{model}}| = 0$ 会完全漏报失配！
- **数学定义**：
  1. 在 `larch_wt.continuous_wavelet_transform`（`larch_wt.py:90–102`）中，原本已算出复数向量 `transform * dk`（长度为 `len(r_grid)`）。将其存入复数矩阵 `complex_map`（形状 `(len(k_centers), len(r_grid))`，`dtype=complex`）。
  2. 当 `normalize=True` 且 `peak > 0` 时，`complex_map = complex_map / peak`，保证 `np.abs(complex_map)` 与 `magnitude` 逐点相等。
  3. 在 `extract_k_profile`（`larch_wt.py:123`）中：若 `wave_map.complex_map is not None`，额外提取沿嵴线轨迹和中心单行的相位 `phase_ridge: np.ndarray`、`phase_slice: np.ndarray`（`np.angle(...)`），并计算双峰之间的最小幅值谷处是否存在相位跳变 `beat_node_phase_jump: float`（单位 rad，用于区分“双元素叠加”与“同元素亚壳层相消拍频节点”）。
  4. 在 `compute_model_wt_residual`（`larch_wt.py:220`）中，保持原有 `residual_map`（幅值差）和 `r_factor`（幅值 R-factor）不变以兼容 B3-1，**新增**：
     $$\text{complex\_residual\_map} = \frac{|W_{\text{data}}(R, k) - W_{\text{model}}(R, k)|}{\text{scale}}$$
     $$\text{complex\_r\_factor} = \frac{\sum_{k, R} |W_{\text{data}}(R, k) - W_{\text{model}}(R, k)|^2}{\sum_{k, R} |W_{\text{data}}(R, k)|^2}$$
     （其中 `scale = float(wt_data.magnitude.max())`）。

#### 2. 小模型逐步修改指南
- **Step 1**：打开 `xafs_core.py` 第 173 行 `class WaveletMap`，在 `sigma: float` 后新增默认字段與属性：
  ```python
  complex_map: np.ndarray | None = None  # 复数 W(R,k)，形状 (len(k_centers), len(r))

  @property
  def phase(self) -> np.ndarray | None:
      return None if self.complex_map is None else np.angle(self.complex_map)
  ```
- **Step 2**：打开 `larch_wt.py` 第 90–102 行 `continuous_wavelet_transform`：
  - 初始化 `c_map = np.zeros((len(k_centers), len(r_grid)), dtype=complex)`；
  - 循环内 `c_row = transform * dk`，`c_map[i, :] = c_row`，`wave_map[i, :] = np.abs(c_row)`；
  - 若 `normalize and peak > 0`，同步 `c_map = c_map / peak`；
  - 返回 `WaveletMap(k_centers=k_centers, r=r_grid, magnitude=wave_map, sigma=float(sigma), complex_map=c_map)`。
- **Step 3**：在 `larch_wt.py:extract_k_profile` 返回字典中追加 `"phase_slice"` 与 `"phase_ridge"`；在 `compute_model_wt_residual` 返回字典中追加 `"complex_residual_map"` 与 `"complex_r_factor"`。

#### 3. 旗舰模型验收标准（新增回归 #37：`test_wt_complex_residual_catches_phase_inversion`）
- 取 `self.lep.chi_k`，令反相模型 `model_inverted = -self.lep.chi_k.chi`（包络模值 100% 相同，但相位差 $\pi$）：
  1. 调用 `res = larch_wt.compute_model_wt_residual(self.lep.chi_k, model_inverted)`；
  2. 验证纯幅值残差 `res["r_factor"] < 1e-12`（证明仅看模值完全失效），而 `res["complex_r_factor"]` 严格等于 `4.0`（容差 `1e-6`，因为 $|W - (-W)|^2 / |W|^2 = 4$）；
  3. 对真实 S7 拟合模型 `fit.model_k`，验证 `res_fit["complex_r_factor"] < 0.25` 且 `wt.complex_map.shape == wt.magnitude.shape`。

---

### B5-2 候选散射元素理论小波包络标尺 `reference_element_k_peaks`（含 $\sigma^2$ 调制）

#### 1. 物理背景与数学契约
- **当前缺陷**：未知样品（L1 状态）做小波切片时，用户不知道当前窗口 `(kmin, kmax, kweight)` 和热无序 $\sigma^2$ 下，C/N/O、S/Cl、3d 金属、5d 重金属的理论 `k_peak` 到底在多少 $\text{Å}^{-1}$（高 $\sigma^2$ 会把 3d 金属的峰拉低 $1\text{–}2\text{ Å}^{-1}$）。
- **数学契约**：
  - 在 `structure_id.py` 新增函数：
    ```python
    def reference_element_k_peaks(
        absorber_symbol: str,
        elements: tuple[str, ...] = ("O", "Cl", "Fe"),
        r_bond: float = 2.10,
        sigma2: float = 0.005,
        kmin: float = 3.0,
        kmax: float = 12.0,
        kweight: int = 2,
        wt_sigma: float = 1.5,
        edge: str | None = None,
    ) -> dict[str, float]:
    ```
  - 对每个候选元素 `el in elements`：
    1. 用 `artemis.octahedron_sites(r_bond, symbol=el)` 生成单壳层模型，调 `artemis.build_feff_input` + `run_feff6l`（使用缓存目录 `_structure_lib/_wt_ref_cache/{absorber}_{edge}_{el}_{r_bond:.2f}` 避免重复跑 FEFF）；
    2. 取第一条单散射路径 `feff0001.dat`，用 `artemis.set_path_parameters(path, s02=1.0, delta_r=0.0, sigma2=sigma2, delta_e0=0.0)` 生成理论 $\chi(k)$（在均匀网格 `k = np.arange(0.0, kmax + 2.0, 0.05)` 上插值）；
    3. 包装为 `ChiK` 后调用 `larch_wt.continuous_wavelet_transform(..., sigma=wt_sigma, kmin=kmin, kmax=kmax, kweight=kweight, rmax=max(4.0, r_bond + 1.0))`；
    4. 在表观峰附近调用 `larch_wt.extract_k_profile` 取出理论 `k_peak`，返回 `{el: float(k_peak)}`。
  - **注意依赖铁律**：函数放在 `structure_id.py`（需要同时调 `artemis` 和 `larch_wt` 时，可在函数体内局部 `import larch_wt`，保持模块顶层依赖图简洁无环）。

#### 2. 小模型逐步修改指南
- **Step 1**：在 `structure_id.py` 末尾（`if __name__ == "__main__":` 之前）实现 `reference_element_k_peaks`，缓存路径放在 `structure_lib.CACHE_ROOT / f"wt_ref_{absorber_symbol}_{edge}_{el}_{int(round(r_bond*100))}"`。
- **Step 2**：在 `xafs_diagnostics.py:write_wt_shell_profiles` 中，新增可选参数 `ref_k_peaks: dict[str, float] | None = None`；若传入，则在切片子图顶部用灰色点划线画出各参考元素的理论 `k_peak` 标尺。

#### 3. 旗舰模型验收标准（新增回归 #38：`test_wt_element_reference_ruler_sigma2_shift`）
- 调用 `peaks_low = structure_id.reference_element_k_peaks("Fe", elements=("O", "Cl", "Fe"), r_bond=2.2, sigma2=0.003)` 与 `peaks_high = structure_id.reference_element_k_peaks("Fe", elements=("Fe",), r_bond=2.2, sigma2=0.015)`：
  1. 原子序数排序单调性：`peaks_low["O"] < peaks_low["Cl"] < peaks_low["Fe"]`，且 `peaks_low["Fe"] - peaks_low["O"] >= 1.2`；
  2. 热无序偏移效应：高无序下 Fe 的理论峰位显著向低 $k$ 移动，`peaks_low["Fe"] - peaks_high["Fe"] >= 0.5`。

---

## 第六批（B6）：Athena 进阶预处理与 XANES 精细定量补齐

### B6-1 多扫描相异度筛查、异常扫描剔除（Outlier Rejection）与方差传播 $\sigma_\mu(E)$

#### 1. 物理背景与数学契约
- **当前缺陷**：`athena.merge_scans(scans)`（`athena.py:91`）做无条件等权平均。多条扫描中若有 1 条因气泡/束流丢失出现局部跳变，会直接污染合并谱；且未输出逐点标准差 $\sigma_\mu(E)$。
- **数学契约**：
  1. `RawScan`（`xafs_core.py:67`）新增可选字段 `sigma_mu: Optional[np.ndarray] = None`。
  2. `athena.merge_scans` 签名扩展为（保持默认参数 100% 向后兼容）：
     ```python
     def merge_scans(scans, reject_outliers: bool = False,
                     z_threshold: float = 3.5) -> RawScan:
     ```
  3. 当 `len(scans) >= 2` 时，在公共能量网格上构建矩阵 $M$（形状 `(n_scans, n_energy)`）：
     - 计算中位谱 $\mu_{\text{med}}(E) = \text{median}(M, \text{axis}=0)$；
     - 每条扫描 $i$ 的均方根偏离度 $d_i = \sqrt{\text{mean}_E \left[(\mu_i(E) - \mu_{\text{med}}(E))^2\right]}$；
     - 当 `reject_outliers=True` 且 `len(scans) >= 3` 时：计算 $d$ 的中位数 $d_{\text{med}} = \text{median}(d)$ 与绝对中位差 $\text{MAD} = \text{median}(|d - d_{\text{med}}|)$。若某条扫描满足 $d_i > d_{\text{med}} + z_{\text{threshold}} \cdot \max(1.4826 \cdot \text{MAD}, 10^{-4})$ 且 $d_i > 0.02 \cdot \text{ptp}(\mu_{\text{med}})$，则将其标记为异常扫描剔除（**保底护栏**：至少保留 2 条最接近中位谱的扫描，绝不全部剔光）。
     - 在保留的扫描集合上计算算术平均 $\bar{\mu}(E)$ 与样本标准差 $\sigma_\mu(E) = \text{std}(M_{\text{kept}}, \text{axis}=0, \text{ddof}=1)$，填入返回的 `RawScan.sigma_mu`，并在 `merged.meta["merge_info"]` 中记录 `{"n_input": len(scans), "n_kept": int(kept.sum()), "rejected_indices": [...], "scan_rms_deviations": [...]}`。

#### 2. 小模型逐步修改指南
- **Step 1**：`xafs_core.py:RawScan` 增加 `sigma_mu: Optional[np.ndarray] = None`；`athena._copy_raw_scan` 同步拷贝 `sigma_mu`。
- **Step 2**：修改 `athena.py:merge_scans` 实现上述逻辑（注意：**当 `reject_outliers=False` 时，合并出的 `energy` 与 `mu` 必须与原代码逐位完全一致**！只是额外填充了 `sigma_mu` 和 `meta["merge_info"]`）。
- **Step 3**：`xafs_pipeline.process_sample`、`run_xafs.py`、`xafs_report.py` 新增可选参数 `reject_outlier_scans: bool = False` 及 CLI `--reject-outlier-scans`。

#### 3. 旗舰模型验收标准（新增回归 #39：`test_merge_scans_outlier_rejection_and_sigma`）
- 读取 `LEPIDO` 的 3 条原始扫描（`.000`, `.100`, `.200`）并复制出第 4 条扫描，在第 4 条扫描上注入 `+0.15` 的吸收跳变畸变：
  1. `m_clean = athena.merge_scans(scans[:3])`：`m_clean.sigma_mu` 非 None，长度等于 `len(m_clean.energy)`，且 `np.all(m_clean.sigma_mu >= 0)`；
  2. `m_rej = athena.merge_scans(scans_4_with_bad, reject_outliers=True)`：`m_rej.meta["merge_info"]["rejected_indices"] == [3]`，且 `np.max(np.abs(m_rej.mu - m_clean.mu)) < 1e-12`。

---

### B6-2 XANES 区域荧光自吸收校正（Tröger / FLUO 能量域反解）

#### 1. 物理背景与数学契约
- **当前缺陷**：S4.5 的 `booth_self_absorption` 只修 $k$ 空间 $\chi(k)$。浓样荧光 XANES 的白线峰被自吸收压扁，导致 XANES LCF 定量失真。而 Larch 内置 `fluo_corr` 若用于全谱会把高能 EXAFS 放大 7.28 倍。因此必须实现**严格限制在 XANES 能量域、带阻尼分母门禁的 Tröger/FLUO 校正**，且绝不污染 EXAFS $\chi(k)$。
- **数学契约**：
  - 在 `athena.py` 新增函数：
    ```python
    def fluo_xanes_self_absorption(
        norm: NormalizedScan,
        formula: str, elem: str, edge: str, line: str,
        thickness_um: float, density: float,
        anginp: float = 45.0, angout: float = 45.0,
    ) -> tuple[NormalizedScan, dict]:
    ```
  - 复用 `athena._booth_mu_terms(formula, elem, edge, line, norm.e0, np.maximum(norm.energy, norm.e0 + 0.1), anginp, angout)`（`athena.py:756`）获取边后基准点 $E_{\text{ref}} = E_0 + 50\text{ eV}$ 处的吸收比参数：
    $$\alpha = \frac{\mu_{\text{bg}}(E_{\text{ref}}) + g \cdot \mu_T(E_f)}{\mu_a(E_{\text{ref}})}$$
    对于厚样极限（或含厚度修正因子 $\eta = 1 - \exp[-(\mu_T/\sin\theta_{\text{in}} + \mu_f/\sin\theta_{\text{out}})\rho d]$），XANES 归一化谱 $\mu_{\text{flat}}(E)$ 满足正演压缩关系：
    $$\mu_{\text{exp}}(E) = \frac{(\alpha + 1)\,\mu_{\text{true}}(E)}{\alpha + \mu_{\text{true}}(E)}$$
    其闭式反解（FLUO / Tröger 方程）为：
    $$\mu_{\text{true}}(E) = \frac{\alpha \, \mu_{\text{exp}}(E)}{\alpha + 1 - \mu_{\text{exp}}(E)}$$
  - **关键物理护栏**：
    1. 仅对 $E \ge E_0 - 15\text{ eV}$ 且 $\mu_{\text{exp}} > 0$ 的点应用反解；
    2. **分母防爆门禁**：令分母 $D(E) = \alpha + 1 - \mu_{\text{exp}}(E)$。若 $\min D(E) < 0.25 \cdot (\alpha + 1)$，将分母钳制在 $0.25 \cdot (\alpha + 1)$ 并标记 `info["reliable"] = False`（防止极端白线尖峰处分母趋零炸出无穷大）；
    3. 返回新的 `NormalizedScan`（更新 `mu_norm` 与 `flattened`），原输入对象不原地修改。

#### 2. 小模型逐步修改指南
- **Step 1**：在 `athena.py` 的 `booth_self_absorption` 下方（约第 950 行）新增 `fluo_xanes_self_absorption`。
- **Step 2**：在 `xafs_pipeline.process_sample` 中，当荧光模式提供了五项自吸收参数时：
  - 先用原始 `norm` 提取 `chi_k` 并走 `booth_self_absorption`（**保持 S4/S4.5 EXAFS 路径 100% 不变**！）；
  - 再对 `norm` 调 `fluo_xanes_self_absorption` 更新 `norm.flattened` 并重算 `extract_xanes_features`，在 `self_absorption["xanes_corrected"] = True` 中记录。

#### 3. 旗舰模型验收标准（新增回归 #40：`test_xanes_fluo_self_absorption_restores_whiteline_and_lcf`）
- 用真实 `LEPIDO` 与 `GOETHITE` 的 flattened XANES 构造真值混合谱 $y_{\text{true}} = 0.5 y_{\text{lep}} + 0.5 y_{\text{goe}}$，按 $\text{Fe}_2\text{O}_3$ 浓样正演公式压缩生成自吸收畸变谱 $y_{\text{atten}}$（白线高度被压低 $>15\%$）：
  1. 调 `fluo_xanes_self_absorption` 反解得到 $y_{\text{corr}}$；
  2. 验证白线高度恢复误差 `|max(y_corr) - max(y_true)| < 0.02`；
  3. 投入 `xafs_lcf.run_lcf`：校正后回收权重与 `0.50` 偏差 `< 0.02`，且 `r_factor` 比未校正时下降 $> 10$ 倍。

---

### B6-3 边前峰（Pre-edge）主边拖尾扣除与多 Pseudo-Voigt 分峰拟合（Wilke 质心）

#### 1. 物理背景与数学契约
- **当前缺陷**：直线基线无法扣除主吸收边的洛伦兹上升拖尾，也无法拆分 $t_{2g}/e_g$ 双峰并给出面积加权质心能量 `pre_edge_centroid`。
- **数学契约**：
  - 在 `athena.py` 新增函数：
    ```python
    def fit_pre_edge_peaks(
        norm: NormalizedScan,
        window: tuple[float, float] = (-14.0, -2.0),
        max_peaks: int = 2,
    ) -> dict:
    ```
  - 定义边前窗口 $[E_{\min}, E_{\max}] = [E_0 + \text{window}[0], E_0 + \text{window}[1]]$：
    1. **基线模型（扣主边拖尾）**：取窗口两端各 $1.5\text{ eV}$ 的邊界带，用线性 + 洛伦兹尾基线 $B(E) = a + b(E - E_0) + \frac{c}{(E - E_{\text{edge}})^2 + \Gamma_0^2}$（$c \ge 0, E_{\text{edge}} = E_0 + 2.0, \Gamma_0 = 3.0$）拟合两端边界点，得到净边前谱 $y_{\text{net}}(E) = \max(\mu_{\text{flat}}(E) - B(E), 0)$。
    2. **Pseudo-Voigt 峰函数**：单峰定义为面积参数化形式：
       $$V(E; A, E_c, \gamma, \eta) = A \left[ \eta \frac{\gamma / \pi}{(E - E_c)^2 + \gamma^2} + (1 - \eta) \frac{\sqrt{\ln 2}}{\gamma \sqrt{\pi}} \exp\left(-\ln 2 \frac{(E - E_c)^2}{\gamma^2}\right) \right]$$
       其中 $A \ge 0$ 为峰积分面积，$E_c \in [E_{\min}+0.5, E_{\max}-0.5]$ 为峰中心，半高半宽 $\gamma \in [0.4, 2.5]\text{ eV}$，洛伦兹权重 $\eta \in [0, 1]$。
    3. 分别用 1 个峰和 2 个峰（当 `max_peaks >= 2`）通过 `scipy.optimize.least_squares` 拟合 $y_{\text{net}}(E)$；当双峰模型的残差平方和 $\text{RSS}_2 < 0.70 \cdot \text{RSS}_1$ 且两峰间距 $|E_{c,1} - E_{c,2}| \ge 0.8\text{ eV}$、次峰面积 $\ge 0.10 \times$ 主峰面积时采用双峰，否则采用单峰。
    4. **面积加权质心能量（Wilke Centroid）**：
       $$E_{\text{centroid}} = \frac{\sum_k A_k E_{c, k}}{\sum_k A_k}, \qquad A_{\text{total}} = \sum_k A_k$$
    5. 返回 `{"centroid_energy": float, "total_area": float, "n_peaks": int, "peaks": [{"area": A, "center": Ec, "hwhm": gamma, "eta": eta}, ...], "r_factor": float}`。
  - 在 `extract_xanes_features`（`athena.py:580`）中调用 `fit_pre_edge_peaks`，将 `pre_edge_centroid`、`pre_edge_voigt_area`、`pre_edge_n_peaks` 写入返回字典。

#### 2. 小模型逐步修改指南
- **Step 1**：在 `athena.py`（约第 685 行，`extract_xanes_features` 旁）实现 `fit_pre_edge_peaks`。
- **Step 2**：在 `extract_xanes_features` 中调用并合并字段；若边前窗口点数 `< 6`，安全降级返回 `pre_edge_centroid = pre_edge_peak_energy`。
- **Step 3**：在 `xafs_report.py:export_indicators_single` 中新增指标行 `pre_edge_centroid`。

#### 3. 旗舰模型验收标准（新增回归 #41：`test_pre_edge_voigt_deconvolution_and_centroid`）
- 在均匀能量网格 $E \in [7090, 7160]\text{ eV}$（$E_0 = 7125.0\text{ eV}$）上构造含主边洛伦兹拖尾 + 两个已知 Pseudo-Voigt 边前峰（峰 1：$E_{c1} = 7113.0\text{ eV}, A_1 = 0.08$；峰 2：$E_{c2} = 7114.6\text{ eV}, A_2 = 0.16$；真值加权质心 $E_{\text{true}} = \frac{0.08 \times 7113.0 + 0.16 \times 7114.6}{0.24} = 7114.067\text{ eV}$）：
  1. 调 `res = athena.fit_pre_edge_peaks(norm)`；
  2. 验证 `res["n_peaks"] == 2`；
  3. 质心能量误差 `abs(res["centroid_energy"] - 7114.067) < 0.15`；两个子峰中心误差均 `< 0.20 eV`。

---

## 第七批（B7）：Artemis 晶体学约束、部分占据与累积量进阶补齐

### B7-1 各向同性晶格膨胀约束（`isotropic_deltar=True`，$\Delta R_i = \alpha \cdot R_{\text{eff}, i}$）

#### 1. 物理背景与数学契约
- **当前缺陷**：多壳层晶体拟合时，每条单散射路径独立分配一个 `(i, "deltar")` 自由变量。对于立方/高对称晶相（如金属箔、NiO 岩盐矿），各壳层距离变化受同一晶格常数缩放支配：$\Delta R_i = \alpha \cdot R_{\text{eff}, i}$。若把 $M$ 个 $\Delta R_i$ 合并为 1 个全局膨胀系数 $\alpha$，立省 $M-1$ 个自由变量，避免远壳层因 $N_{\text{idp}}$ 不足被截断！
- **数学契约**：
  1. `artemis.fit_r_space`（`artemis.py:543`）新增可选参数 `isotropic_deltar: bool = False`（默认 `False`）。
  2. 当 `isotropic_deltar=True` 且非单路径时：
     - 不再为每条单散射路径建立 `(i, "deltar")` 槽位，而是建立唯一一个全局槽位 `(-5, "alpha")`，初值 `0.0`，上下界 `(-0.08, +0.08)`（对应 $\pm 8\%$ 晶格应变）；
     - 在 `apply_params(x)` 中，对所有单散射路径 $i$ 写入 `by_path[i]["deltar"] = alpha * float(paths[i].reff)`；
     - 在构造返回 `PathFit` 时：第 $i$ 条路径的 `delta_r = alpha * float(path.reff)`，误差按线性传播 `delta_r_err = alpha_err * float(path.reff)`；
     - `FitResult` 新增可选字段 `lattice_strain: float | None = None`、`lattice_strain_err: float | None = None`（记录 $\alpha$ 及其标准误）。
  3. 在 `xafs_fit.select_shell_paths`（`xafs_fit.py:318`）与 `fit_identified` 中新增 `isotropic_deltar: bool = False`：
     - 当 `isotropic_deltar=True` 时，信息量上限公式从每壳 3 个变量（$3n + 1 < N_{\text{idp}}$）放宽为每壳 2 个变量加 2 个全局变量（$2n + 2 < N_{\text{idp}}$，即 $N_i, \sigma^2_i$ 加全局 $\alpha, \Delta E_0$）！
  4. `run_xafs.py` 与 `xafs_report.py` 新增 CLI `--isotropic-deltar`。

#### 2. 小模型逐步修改指南
- **Step 1**：`xafs_core.py:FitResult` 新增 `lattice_strain: float | None = None`、`lattice_strain_err: float | None = None`。
- **Step 2**：`artemis.py:fit_r_space` 接入 `isotropic_deltar: bool = False`：在构建 `slots`、`param_names`（命名为 `"lattice_alpha"`）、`apply_params` 及 `PathFit` 生成处分别处理 `(-5, "alpha")`。
- **Step 3**：`xafs_fit.py`、`run_xafs.py`、`xafs_report.py` 透传 `isotropic_deltar`。

#### 3. 旗舰模型验收标准（新增回归 #42：`test_isotropic_lattice_expansion_constraint`）
- 对 `NiO`（`rmax=4.3`，含 4 个单散射壳层）：
  1. 分别运行 `fit_free = xafs_fit.fit_identified(res_nio, absorber_symbol="Ni", rmax=4.3, isotropic_deltar=False)` 与 `fit_iso = xafs_fit.fit_identified(res_nio, absorber_symbol="Ni", rmax=4.3, isotropic_deltar=True)`；
  2. 验证变量数精确减少：`fit_iso["fit"].n_variables == fit_free["fit"].n_variables - (len(fit_iso["fit"].paths) - 1)`（4 壳层从 13 个变量降为 **10 个变量**）；
  3. 验证所有壳层的应变比值严格一致：对任意两条路径 $p_a, p_b$，`abs(p_a.delta_r / p_a.reff - fit_iso["fit"].lattice_strain) < 1e-12`；
  4. 拟合收敛且 `fit_iso["fit"].r_factor < 0.08`。

---

### B7-2 CIF 部分占据（Partial Occupancy）与空位/共占位简并度自动缩放

#### 1. 物理背景与数学契约
- **当前缺陷**：`artemis.read_cif_sites`（`artemis.py:175`）忽略 CIF 中的 `_atom_site_occupancy`，把所有位点按满占据 `1.0` 展开；若某位点是掺杂共占位（如同一坐标 `0.5 Fe + 0.5 Co`）或含氧空位（`occ = 0.75`），FEFF 会因两个原子距离为 0 崩溃，或把理论简并度高估。
- **数学契约**：
  1. 在 `artemis.read_cif_sites` 中读取 `site.occ`（`gemmi` 的 `site.occ` 属性，默认 `1.0`，钳制在 `[0.0, 1.0]`）：
     - 新增可选参数 `return_occupancy: bool = False`（默认 `False` 返回原有的 `(x, y, z, symbol)` 4 元组，**保证 100% 向后兼容**；当 `return_occupancy=True` 时返回 `(x, y, z, symbol, occ)` 5 元组）。
     - **共占位去重（防 FEFF 崩溃）**：对展开后的散射原子，若两个散射原子坐标距离 `< 0.08 Å`（同一晶位共占据），合并为一个位点：元素符号取占据度较大者，总有效占据度取 $\min(1.0, \sum \text{occ})$（或按元素分别记录）。
     - 吸收中心原子若有多个同元素不等价位点或不同占据度，按 `occ * multiplicity` 选择主位点。
  2. 在 `structure_lib.Motif` 或 `xafs_fit._ss_rows` 中：当 CIF 位点存在 `occ < 0.999` 时，计算该距离壳层的平均占据因子 $\overline{\text{occ}}(R_{\text{shell}})$，将 FEFF 几何简并度缩放为 $N_{\text{cryst}} = N_{\text{geom}} \cdot \overline{\text{occ}}$。

#### 2. 小模型逐步修改指南
- **Step 1**：修改 `artemis.py:read_cif_sites`（第 175–283 行），读取 `float( getattr(site, "occ", 1.0) )`，在周期像展开去重时合并 `< 0.08 Å` 的共占位点，并根据 `return_occupancy` 返回 4 元组或 5 元组。
- **Step 2**：在 `artemis.py` 新增辅助函数 `cif_shell_occupancy_scale(cif_path, absorber: str, rmax: float = 6.0) -> list[tuple[float, str, float]]`，返回每个 `(distance, element, occ)` 供壳层简并度缩放。
- **Step 3**：在 `xafs_fit.py:_ss_rows` 中，若 `motif.cif_path` 存在，按距离匹配（容差 `0.05 Å`）乘上平均占据度 $\overline{\text{occ}}$。

#### 3. 旗舰模型验收标准（新增回归 #43：`test_cif_partial_occupancy_and_vacancy_scaling`）
- 构造一个含氧空位（O 位点 `_atom_site_occupancy 0.75`）和共占位（同一坐标放两个原子，occ 分别为 `0.6` 和 `0.4`）的临时 CIF：
  1. `sites4 = artemis.read_cif_sites(cif, absorber="Fe")` 返回 4 元组且不抛异常；
  2. `sites5 = artemis.read_cif_sites(cif, absorber="Fe", return_occupancy=True)` 返回 5 元组，第一壳 O 的 `occ == 0.75`，且共占位点已合并去重（无距离 `< 0.08 Å` 的重叠散射体），丢进 `build_feff_input` + `run_feff6l` 正常生成路径；
  3. 经 `xafs_fit.select_shell_paths` 选出的第一壳简并度 `degens[0]` 精确等于满占据的 `0.75` 倍（如从 `6.0` 变为 `4.5`）。

---

### B7-3 四阶累积量 $C_4$ 与非高斯径向分布函数 $g(R)$ 重构

#### 1. 物理背景与数学契约
- **数学契约**：
  1. **四阶累积量 $C_4$**：在 EXAFS 累积量展开中，振幅项含 $\exp(-2 k^2 \sigma^2 + \frac{2}{3} k^4 C_4)$。Larch 的 `FeffPathGroup` 参数名为 `"fourth"`。
     - `PathFit`（`xafs_core.py:182`）新增 `c4: Optional[float] = None`、`c4_err: Optional[float] = None`、`c4_at_bound: bool = False`。
     - `artemis.py` 定义 `C4_LIMIT = 2e-4`（单位 $\text{Å}^4$）；`set_path_parameters` 新增 `fourth=None`（映射 `"fourth": fourth`）；`fit_r_space` 新增 `vary_c4=False`（槽位 `(i, "c4")`，初值 `0.0`，边界 `[-C4_LIMIT, +C4_LIMIT]`，`|c4| >= 0.8 * C4_LIMIT` 置 `c4_at_bound = True`）。
  2. **非高斯径向分布函数 $g(R)$ 重构**：
     - 在 `artemis.py` 新增纯数学函数：
       ```python
       def reconstruct_pair_distribution(
           r_grid: np.ndarray,
           n: float,
           r0: float,
           sigma2: float,
           c3: float = 0.0,
           c4: float = 0.0,
       ) -> dict:
       ```
     - 采用 Gauss-Hermite（Gram-Charlier / Edgeworth）累积量展开：令归一化无量纲坐标 $z = \frac{R - r_0}{\sigma}$（其中 $\sigma = \sqrt{\max(\sigma^2, 10^{-6})}$），标准高斯分布 $\phi(z) = \frac{1}{\sqrt{2\pi}\sigma} \exp(-z^2 / 2)$，厄米多项式 $H_3(z) = z^3 - 3z$，$H_4(z) = z^4 - 6z^2 + 3$，无量纲偏度 $\gamma_1 = \frac{C_3}{\sigma^3}$，超额峰度 $\gamma_2 = \frac{C_4}{\sigma^4}$：
       $$g_{\text{raw}}(R) = N \cdot \phi(z) \cdot \max\left(0.0, \; 1 + \frac{\gamma_1}{6} H_3(z) + \frac{\gamma_2}{24} H_4(z)\right)$$
     - 对 $g_{\text{raw}}(R)$ 按梯形积分重归一使其总面积严格等于配位数 $N$，返回 `{"r": r_grid, "g_r": g_norm, "r_mode": float(r_grid[np.argmax(g_norm)]), "r_mean": float(r0), "skewness": float(gamma_1), "area": float(np.trapezoid(g_norm, r_grid))}`。

#### 2. 小模型逐步修改指南
- **Step 1**：`xafs_core.py:PathFit` 增加 `c4`, `c4_err`, `c4_at_bound`。
- **Step 2**：`artemis.py` 增加 `C4_LIMIT = 2e-4`、`set_path_parameters(..., fourth=None)`、`fit_r_space(..., vary_c4=False)` 及 `reconstruct_pair_distribution`。
- **Step 3**：`xafs_fit.py` 的 `fit_identified` 支持可选 `vary_c4=False`，并在 `flags` 中同步 `c4_at_bound`（若撞界则 `cn_reliable = False`）。

#### 3. 旗舰模型验收标准（新增回归 #44：`test_c4_cumulant_and_asymmetric_rdf_reconstruction`）
- 1) 用 `set_path_parameters(..., fourth=6e-5)` 生成合成单壳谱，验证 `fit_r_space(..., vary_c4=True)` 回收的 `c4` 与真值 `6e-5` 误差 `< 2.5e-5`；
- 2) 调 `rdf = artemis.reconstruct_pair_distribution(np.linspace(1.5, 2.7, 601), n=6.0, r0=2.05, sigma2=0.006, c3=2.5e-4)`：验证积分面积 `abs(rdf["area"] - 6.0) < 1e-3`，且当 $C_3 > 0$ 时右半侧长尾积分面积大于左半侧（峰位 `r_mode < r_mean` 体现正偏态分布物理特性）。

---

## 第八批（B8）：多数据集/双吸收边/变温全局联合拟合（Co-refinement）

### B8-1 多数据集全局受约束联合拟合 `fit_multi_datasets`（含双边守恒与变温爱因斯坦模型）

#### 1. 物理背景与数学契约
- **解决痛点**：
  1. **双金属双吸收边联合拟合**（如 Pt-Co 同时拟合 Pt $L_3$ 和 Co $K$ 边）：强制异核键满足 $R_{\text{A-B}} = R_{\text{B-A}}$、$\sigma^2_{\text{A-B}} = \sigma^2_{\text{B-A}}$；
  2. **变温多谱联合拟合**（如 100 K 与 300 K）：用爱因斯坦温度 $\Theta_E$ 与静态无序 $\sigma^2_{\text{static}}$ 两个全局参数统一驱动所有温度下的热无序 $\sigma^2(T)$，彻底解耦 $N$ 与 $\sigma^2$。
- **数学契约**：
  - 在 `artemis.py` 新增函数：
    ```python
    def fit_multi_datasets(
        datasets: list[dict],
        shared_bonds: list[dict] | None = None,
        einstein_models: list[dict] | None = None,
    ) -> dict:
    ```
  - **输入格式**：
    - `datasets`：列表，每项为 `{"name": str, "chi_k": ChiK, "ft": FourierTransform, "datfiles": list, "rmin": float, "rmax": float, "kmin": float, "kmax": float, "kweight": int, "s02": float, "vary_n": bool | list, "path_degen": list | None, "path_labels": list | None, "temperature_K": float | None}`。
    - `shared_bonds`：跨数据集路径绑定列表，每项为 `{"target": (ds_idx_b, path_idx_b), "source": (ds_idx_a, path_idx_a), "share_deltar": bool, "share_sigma2": bool, "degen_ratio": float | None}`。被绑定的 `target` 路径不建独立槽位，直接复用 `source` 路径的 $\Delta R$、$\sigma^2$（及可选按化学计量比缩放的 $N_b = \text{degen\_ratio} \cdot N_a$）。
    - `einstein_models`：变温热无序模型列表，每项为 `{"paths": [(ds_idx, path_idx), ...], "reduced_mass_amu": float, "theta_e_init": float = 350.0}`。
      - 对每个爱因斯坦组，建立两个全局自由变量：静态无序 $\sigma^2_{\text{static}} \in [0.0, 0.015]\text{ Å}^2$ 与爱因斯坦温度 $\Theta_E \in [80.0, 1200.0]\text{ K}$；
      - 组内所有指定 `(ds_idx, path_idx)` 不建独立 `sigma2` 槽位，其 $\sigma^2$ 由对应数据集温度 $T = \text{datasets}[\text{ds\_idx}][\text{"temperature\_K"}]$ 按量子统计谐振子精确公式计算（内置辅助函数 `einstein_sigma2(T, theta_e, reduced_mass_amu, sigma2_static=0.0)`）：
        $$\sigma^2(T) = \sigma^2_{\text{static}} + \frac{\hbar^2}{2 \mu k_B \Theta_E} \coth\left(\frac{\Theta_E}{2 T}\right), \qquad \frac{\hbar^2}{2 k_B} \approx 24.2543 \text{ amu}\cdot\text{Å}^2\cdot\text{K}$$
        即：
        $$\sigma^2_{\text{thermal}}(T) = \frac{24.2543}{\mu_{\text{amu}} \cdot \Theta_E} \cdot \frac{1}{\tanh\left(\frac{\Theta_E}{2 T}\right)}$$
  - **代价函数与统计量**：
    - 每个数据集 $d$ 在自身 `mask_d` 内的实/虚部残差向量除以该数据集实验谱在 `mask_d` 内的均方根范数 $s_d = \sqrt{\text{mean}(|\chi_{\text{data}, d}(R)|^2)}$，再拼接为全局残差向量；
    - 总独立点数 $N_{\text{idp, total}} = \sum_d N_{\text{idp}, d}$，全局自由变量数 $N_{\text{var, global}}$ 为实际优化槽位数；
    - 返回 `{"datasets": [FitResult, ...], "n_independent_total": float, "n_variables_total": int, "global_r_factor": float, "einstein_results": [{"theta_e": float, "theta_e_err": float, "sigma2_static": float, "sigma2_static_err": float}, ...], "parameter_names": list, "correlation": np.ndarray}`。

#### 2. 小模型逐步修改指南
- **Step 1**：在 `artemis.py` 中先实现纯函数 `einstein_sigma2(temperature_K: float, theta_e: float, reduced_mass_amu: float, sigma2_static: float = 0.0) -> float` 并校验常量 `HBAR2_OVER_2KB = 24.254316`（对 Fe–O 键：$\mu = \frac{55.845 \times 16.0}{55.845 + 16.0} = 12.436\text{ amu}$，当 $\Theta_E = 500\text{ K}, T = 300\text{ K}$ 时，$\sigma^2_{\text{thermal}} \approx 0.00442\text{ Å}^2$，量级完全正确）。
- **Step 2**：在 `artemis.py` 末尾实现 `fit_multi_datasets`，复用 `load_feff_path`、`set_path_parameters`、`exafs_model_k` 与 `scipy.optimize.least_squares`。

#### 3. 旗舰模型验收标准（新增回归 #45：`test_multi_dataset_corefinement_and_einstein_model`）
- 用 `LEPIDO` 第一壳 Fe–O 路径（折合质量 $\mu = 12.436\text{ amu}$）在真值 $\Theta_E = 420.0\text{ K}, \sigma^2_{\text{static}} = 0.0015\text{ Å}^2, \Delta R = -0.02\text{ Å}, N = 6.0$ 下合成 $T_1 = 80\text{ K}$ 与 $T_2 = 350\text{ K}$ 两条实验谱：
  1. 调用 `artemis.fit_multi_datasets` 同时绑定两条谱的 $\Delta R$、$N$、$\Delta E_0$ 并启用 `einstein_models`；
  2. 验证全局自由变量数从独立拟合的 `2 × 4 = 8` 个降为 **5 个**（`N, deltar, e0, sigma2_static, theta_e`）；
  3. 回收的 $\Theta_E$ 与真值 `420.0 K` 误差 `< 15.0 K`，回收的 $\sigma^2_{\text{static}}$ 与 `0.0015` 误差 `< 5e-4`，两条数据集各自的 `r_factor < 1e-4`。

---

## 小模型（Execution Model）防呆自检清单（每完成一步打勾）

小模型在编写代码时，必须严格按以下顺序逐批推进，**禁止跨批跳步**：

- [ ] **准备检查**：确认使用环境 `C:\Users\lizha\.workbuddy\binaries\python\versions\3.13.12\python.exe`，先跑一次 `python run_ci.py`，确认起点 36/36 测试全绿。
- [ ] **B5 批次检查**：
  - [ ] `WaveletMap.complex_map` 形状是否与 `magnitude` 完全一致 `(len(k_centers), len(r))`？
  - [ ] `normalize=True` 时是否同步对 `complex_map` 除以了同一个 `peak`？
  - [ ] `reference_element_k_peaks` 是否使用了缓存目录避免每次重复跑 FEFF6？
  - [ ] 新增 #37、#38 两条回归测试后，`python regression_test.py` 是否显示 **Ran 38 tests — OK**？
- [ ] **B6 批次检查**：
  - [ ] `merge_scans` 在默认 `reject_outliers=False` 下输出的 `energy` 和 `mu` 是否与改动前逐位一致？
  - [ ] `fluo_xanes_self_absorption` 是否有分母下界保护 `>= 0.25 * (alpha + 1)`？是否严格只改 XANES 而不动 `chi_k`？
  - [ ] `fit_pre_edge_peaks` 的面积加权质心能量公式分母为 0 时是否有安全回退？
  - [ ] 新增 #39、#40、#41 三条回归测试后，`python regression_test.py` 是否显示 **Ran 41 tests — OK**？
- [ ] **B7 批次检查**：
  - [ ] `isotropic_deltar=True` 时，单散射路径是否不再创建独立 `(i, "deltar")` 槽位，而是共用唯一 `(-5, "alpha")`？
  - [ ] `read_cif_sites` 默认 `return_occupancy=False` 时是否仍返回 4 元组 `(x, y, z, symbol)`（不破坏既有 CIF 调用）？
  - [ ] `reconstruct_pair_distribution` 输出的 $g(R)$ 是否经过了非负截断与梯形积分面积归一化？
  - [ ] 新增 #42、#43、#44 三条回归测试后，`python regression_test.py` 是否显示 **Ran 44 tests — OK**？
- [ ] **B8 批次检查**：
  - [ ] `einstein_sigma2` 的物理常数 `24.254316` 是否准确？温度 $T \le 0$ 时是否抛 `XAFSConfigError`？
  - [ ] `fit_multi_datasets` 每个数据集的残差是否除以了该数据集自身的 RMS 范数 $s_d$？
  - [ ] 新增 #45 回归测试后，`python regression_test.py` 是否显示 **Ran 45 tests — OK**？
- [ ] **收尾同步**：更新 `README.md` 与 `.trae/skills/xafs-pipeline-report/references/pipeline-reference.md` 的回归计数（36 → 45 项）及 B5–B8 新函数说明。

---

## 旗舰模型（Grading Model）终审验收 Checklist（Phase I + Phase II 全量判分）

当小模型提交全部代码后，旗舰模型将按以下 **7 条硬指标**逐项执行命令核查，任何一条不达标即判退回重做：

1. **全量 CI 与压力测试一键通关**：
   ```powershell
   python run_ci.py --stress
   ```
   - 核查 `ci_report.md`：`py_compile` 退出码 `0`；`regression_test.py` 退出码 `0`，且测试总数达到 **45 项**（#1–#27 初始基线 + #28–#36 Phase I + #37–#45 Phase II）；`stress_test_{co,ni,cu,nio}.py` 全部退出码 `0`。
2. **默认路径向后兼容零退化核查**：
   ```powershell
   python xafs_report.py data/练手数据/02_Fe_standards/Fe_lepidocrocite.000 data/练手数据/02_Fe_standards/Fe_lepidocrocite.100 data/练手数据/02_Fe_standards/Fe_lepidocrocite.200 --absorber Fe --no-strict -o xafs_out/_grade_default
   ```
   - 核查 `r_factor`（`0.0273`）、第一壳 `N`（`6.35`）、`R`（`2.002 Å`）及默认 7 列拟合表与历史基线逐位一致。
3. **CLI 完整性核查**：
   - `python run_xafs.py --help` 与 `python xafs_report.py --help` 均包含 Phase I 与 Phase II 全部开关：`--fit-kweights`、`--include-ms`、`--vary-c3`、`--deglitch`、`--deadtime`、`--lcf-standards`、`--reject-outlier-scans`、`--isotropic-deltar`。
4. **Phase II 关键数学与物理契约核查**：
   - **B5-1**：反相信号 $\chi_{\text{model}} = -\chi_{\text{data}}$ 的二维复数小波残差 `complex_r_factor == 4.0`（幅值残差为 `0.0`），证明复数相位盲区已被消除；
   - **B5-2**：理论小波标尺满足 $k_{\text{peak}}(\text{O}) < k_{\text{peak}}(\text{Cl}) < k_{\text{peak}}(\text{Fe})$ 且随 $\sigma^2$ 增大向低 $k$ 漂移；
   - **B6-1 / B6-2 / B6-3**：异常扫描被 MAD 门禁剔除且至少保底保留 2 条；XANES FLUO 自吸收校正带分母防爆且不改 `chi_k`；边前峰 Pseudo-Voigt 分峰质心误差 `< 0.15 eV`；
   - **B7-1 / B7-2 / B7-3**：`isotropic_deltar=True` 严格将 $M$ 条路径的 $\Delta R$ 变量缩减为 1 个全局 $\alpha$；部分占据 CIF 正确缩放简并度且共占位不崩溃；非高斯 $g(R)$ 积分面积等于 $N$；
   - **B8-1**：`fit_multi_datasets` 正确执行数据集间 RMS 均衡与变温爱因斯坦模型 $\sigma^2(T)$ 约束，回收 $\Theta_E$ 误差 `< 15 K`。
5. **模块单向依赖无环核查**：
   - 依次独立运行 `python xafs_core.py`、`python athena.py`、`python artemis.py`、`python larch_wt.py`、`python xafs_lcf.py`，全部退出码 `0`，无循环依赖或编码报错。
6. **S₀² 与门禁神圣不可侵犯核查**：
   - 核查 `absorber.py` 中 Fe/Co/Ni/Cu/Pt 的 $S_0^2$ 值未被篡改，`xafs_fit.py` 中 `cn_reliable` 门禁阈值（`r_factor <= 0.08`, `|delta_e0| < 15.0`, `sigma2_at_bound`, `c3_at_bound`, `c4_at_bound`）严密生效。
7. **文档同步核查**：
   - `README.md` 与 `.trae/skills/xafs-pipeline-report/references/pipeline-reference.md` 已同步更新 **45 项契约回归**与 B5–B8 新接口文档。

---

## 第三阶段（Phase III：B9–B12）工程闭环、原位解谱、物理深度与统计推断（回归 45 → 60 项）

> **目标**：缝合 Phase I–II 留下的 8 个“半打通断层/边界漏洞”，并补齐对标 Athena / Artemis / Larch-WT 的 10 项高阶能力。
> **铁律不变**：默认参数零退化（`r_factor≈0.0273`、第一壳 `N≈6.35`、`R≈2.002 Å`、默认 7 列拟合表）；S₀² 表神圣不可侵犯；无法定量时显式降级或拒绝。

### 全局执行顺序总览（Phase III）

| 批次 | 编号 | 功能模块 | 修改文件 | 新增回归（`regression_test.py`） |
|---|---|---|---|---|
| **第九批** | **B9-1** | 全参数配置文件同步 + `--vary-c4` CLI + $N_{\text{idp}}$ 降级 + 非高斯 $g(R)$ 自动导出 | `xafs_config.py`, `run_xafs.py`, `xafs_report.py`, `xafs_fit.py` | `test_config_file_covers_all_cli_switches` (#46), `test_vary_c4_cli_nidp_gate_and_rdf_export` (#47) |
| | **B9-2** | $\sigma_\mu(E)$ 误差带/高 $k$ 信噪比 + WT 元素标尺自动挂载 + $K/L$ 边 XANES 智能分流（$L_3$ 白线净面积） | `athena.py`, `xafs_pipeline.py`, `run_xafs.py`, `xafs_report.py`, `xafs_diagnostics.py` | `test_sigma_mu_export_and_high_k_snr_gate` (#48), `test_wt_ruler_auto_mount_and_l3_whiteline_area` (#49) |
| | **B9-3** | `.prj` 工程容器读写 + 多通道荧光探测器（MCA）坏道剔除与加权合并 | `xafs_io.py`, `athena.py` | `test_prj_roundtrip_and_multichannel_mca_merge` (#50) |
| **第十批** | **B10-1** | 短边后窗口 Cromer–Liberman（MBACK）归一化 + AUTOBK `chi_std` 先验与 `clamp_lo/hi` | `athena.py` | `test_mback_short_window_normalization` (#51), `test_autobk_chi_std_prior_and_clamp` (#52) |
| | **B10-2** | PCA Malinowski Target Transformation（`SPOIL`）+ 原位序列谱 MCR-ALS 盲解谱 | `xafs_lcf.py` | `test_pca_target_transform_spoil_and_mcr_als` (#53) |
| **第十一批** | **B11-1** | FEFF 相移校正 FT/WT（Phase-Corrected）+ 异核共占位（$\Delta Z\ge3$）双分支路径 | `athena.py`, `larch_wt.py`, `artemis.py` | `test_phase_corrected_ft_and_wt` (#54), `test_hetero_cocupancy_dual_branch_paths` (#55) |
| | **B11-2** | 关联 Debye 声子模型 + 键价和（BVS）计算与软约束 | `artemis.py`, `xafs_core.py`, `xafs_fit.py` | `test_correlated_debye_sigma2_contract` (#56), `test_bond_valence_sum_restraint` (#57) |
| | **B11-3** | 多空间拟合切换（$R$ / $q$ / 二维复数小波 `wt`） | `artemis.py` | `test_multispace_q_and_wt_fit` (#58) |
| **第十二批** | **B12-1** | 双吸收边独立 $\Delta E_0$ 分组 + $\text{AIC}_c/\text{BIC}$ 模型选择 | `artemis.py`, `xafs_core.py` | `test_multi_edge_e0_groups_and_aic_bic` (#59) |
| | **B12-2** | 单原子/双原子催化剂（SAC/DAC）非周期团簇生成器 + 催化元素扩充 | `artemis.py`, `absorber.py` | `test_sac_dac_cluster_generator_and_element_expansion` (#60) |

---

### B9-1 工程闭环：配置同步、`--vary-c4` 与非高斯 $g(R)$ 落盘

#### 数学与接口契约
1. `xafs_config.DEFAULTS` 必须覆盖 Phase I–III 全部 CLI 开关：`deglitch`、`deadtime`、`reject_outlier_scans`、`fit_kweights`、`include_ms`、`vary_c3`、`vary_c4`、`isotropic_deltar`、`lcf_standards`、`shared_sigma2`、`mback`、`chi_std`、`clamp_lo`、`clamp_hi`、`fit_space`、`target_transform`、`mcr_als`。
2. `xafs_config.apply_config_to_kwargs(config, kwargs, explicit=None)`：布尔开关只在 CLI 未显式给出时由配置打开；`explicit` 为 CLI 显式 dest 集合时优先级最高。
3. `xafs_fit.fit_identified(..., vary_c4=...)` 增加与 C3 完全对称的 $N_{\text{idp}}$ 防护：若放开 C4 后 `n_var >= N_idp`，先降级为仅第一壳，仍超限则关闭并写 `flags.c4_dropped_by_nidp=True`。
4. 报告层：只要拟合中存在 `c3` 或 `c4`，必须对每条路径调用 `artemis.reconstruct_pair_distribution`，输出 `<stem>_rdf.csv`（列：`label, r_A, g_r, r_mode, r_mean, skewness, area`）与 `_fit.json` 中的 `rdf` 摘要。

#### 验收标准
- **#46**：写一个含全部 17 个键的临时 YAML/JSON，`load_config` + `apply_config_to_kwargs` 后逐键生效；未出现在配置中的键保持默认。
- **#47**：合成单壳谱 `vary_c4=True` 回收 $C_4$；`vary_c4=True` 且路径数过多时触发 `c4_dropped_by_nidp=True`；真实 LEPIDO 运行时生成 `_rdf.csv` 且 `area≈N`。

---

### B9-2 误差传播落地、元素标尺自动挂载与边类型分流

#### 数学与接口契约
1. `athena.high_k_snr(chi_k, sigma_mu=None, edge_step=None, kweight=2)`：把 $\sigma_\mu(E)$ 换算为 $\sigma_\chi(k)=\sigma_\mu/|\text{edge\_step}|$，输出 `snr_k`（逐点）、`k_snr2`（滚动中位数 SNR 首次跌破 2 的 $k$）、`reliable_kmax`。无 `sigma_mu` 时返回 `None` 并说明原因，不伪造。
2. `xafs_pipeline` 把 `merged.sigma_mu` 保留到 `PipelineResult`；`run_xafs` 导出的 `<stem>_mu.dat` 增加 `sigma_mu` 列；`xafs_diagnostics._s1` 绘制 $\mu\pm\sigma_\mu$ 误差带。
3. `xafs_report.py` 在调用 `write_wt_shell_profiles` 时自动计算并传入 `ref_k_peaks`：元素集合 = 拟合路径散射元素 ∪ `("O", absorber)`，带缓存。
4. `athena.extract_xanes_features(norm, edge="K", absorber=None)`：
   - $K$ 边：保留 `fit_pre_edge_peaks`；
   - $L_2/L_3$ 边：不拟合边前峰，改为 `fit_whiteline_area`：对 $[E_0-15, E_0+40]$ 用 `arctan` 阶跃 + 伪 Voigt 拟合并积分扣除阶跃后的净白线面积，输出 `white_line_area`、`white_line_area_norm`、`whiteline_peak_height`。

#### 验收标准
- **#48**：三扫描合并后 `_mu.dat` 含 `sigma_mu` 列；注入 $\sigma_\mu$ 后 `k_snr2` 随噪声增大单调下降；无 `sigma_mu` 时返回 `None` 且 `reason` 明确。
- **#49**：`extract_xanes_features(..., edge="L3")` 不产生 `pre_edge_n_peaks`，而给出有限 `white_line_area`；`edge="K"` 保持原字段；WT 诊断自动带 `ref_k_peaks`。

---

### B9-3 `.prj` 工程容器与多通道 MCA

#### 接口契约
1. `xafs_io.write_prj(path, payload)` / `read_prj(path)`：gzip 压缩的 JSON 容器（`{"xafs_prj": 1, "scans": [...], "standards": {...}, "meta": {...}}`）；兼容未压缩 JSON；对 Athena 原生 Perl `.prj` 给出明确 `XAFSConfigError` 提示“请先在 Athena 导出为文本/JSON”。
2. `xafs_io.read_xdi` 识别 `mca1..mcaN` / `if1..ifN` 与 `dt1..dtN` 多通道列，生成 `scan.meta["mca"]={"channels": ndarray(n, ch), "dtime": ndarray|None}`，`mu` 先用通道求和。
3. `athena.merge_mca_channels(channels, dtime=None, z_threshold=3.5, min_channels=2)`：逐道死时间校正 → 以通道间中位谱 RMS+MAD 剔除坏道（保底 2 道）→ 按边阶（或均值）加权求和；返回 `(summed, info)`，`info` 含 `rejected`, `n_kept`, `weights`。

#### 验收标准
**#50**：`.prj` 写入/读取 round-trip 的 `energy/mu` 逐位相同；合成 8 通道（2 道含 Bragg spike）正确剔除 2 道且保留 6 道；求和谱与干净通道加权结果一致。

---

### B10-1 MBACK 归一化与 AUTOBK 先验

#### 数学与接口契约
1. `athena.normalize_mback(scan, e0, elem, edge="K", pre_window=(-200,-30), post_window=(60,150), fit_order=2)`：
   - 理论截面 $\mu_{\text{th}}(E)$ 用 `xraydb.mu_elam(elem, E)`；
   - 在边前与边后窗口联合最小二乘拟合 $\mu_{\text{meas}} \approx a\,\mu_{\text{th}} + b$（边后窗口可短至 60 eV）；
   - 归一化 $\mu_{\text{norm}} = (\mu_{\text{meas}}-b)/(a\,\Delta\mu_{\text{th}})$，其中 $\Delta\mu_{\text{th}}=\mu_{\text{th}}(E_0+50)-\mu_{\text{th}}(E_0-50)$；
   - 返回 `NormalizedScan` 并带 `meta["mback"]={a,b,theory_edge_step,post_window}`；`a<=0` 时干净拒绝。
2. `athena.extract_chi_k(norm, rbkg=1.2, kweight=2, chi_std=None, clamp_lo=None, clamp_hi=None)`：把参数透传 larch `autobk`（不支持的版本显式降级并记录）；`chi_std` 作为先验标准谱保护低 $R$ 真实信号，`clamp_lo/hi` 抑制样条末端甩尾。

#### 验收标准
- **#51**：用 Fe K 理论截面合成短边后窗口（$E_0+60\sim+150$ eV）的 $\mu$，`normalize_mback` 回收的 `edge_step` 相对误差 `< 3%`，而传统 `normalize_absorption(post_window=(60,150))` 的平台偏离更远。
- **#52**：`chi_std` 传入后第一壳 $R<1.2$ Å 的泄漏降低；`clamp_lo/hi` 生效且 `ChiK` 有限；非法 `chi_std` 长度干净报错。

---

### B10-2 PCA 目标转换与 MCR-ALS

#### 数学与接口契约
1. `xafs_lcf.run_pca` 增加返回 `components`（$U$ 矩阵，`(n_energy, n_samples)`）与 `n_components_ind`。
2. `xafs_lcf.target_transform(pca, candidates, n_components=None, spoil_good=3.0, spoil_poor=6.0)`：
   - 对每个候选谱 $c$（插值到 PCA 网格并扣均值），投影 $p=U_k^T c$，重建 $r=U_k p$；
   - $\text{RET}=\sqrt{\sum(c-r)^2/(n-k)}$，$\text{REP}=\sqrt{\sum r^2/n}$，$\text{SPOIL}=\text{RET}/\text{REP}$；
   - 返回 `{name: {spoils, ret, rep, rating}}`，`rating ∈ good/fair/poor/unacceptable`；`REP=0` 时 `SPOIL=inf` 并标记 `degenerate`。
3. `xafs_lcf.mcr_als(spectra, n_components, max_iter=200, tol=1e-6, nonneg=True, closure=False, seed=0)`：
   - 输入 `[(name, x, y)]`，构造 $Y$；PCA 初始化 $C,S$；
   - 交替非负最小二乘：$C \leftarrow \arg\min_{C\ge0}\|Y-SC\|$，$S \leftarrow \arg\min_{S\ge0}\|Y-SC\|$；
   - 返回 `concentrations`（$k\times n$）、`components`（$n\times k$）、`r_factor`、`iterations`、`converged`。

#### 验收标准
**#53**：3 个纯端元的线性混合（浓度和为 1）+ 0.5% 噪声：`target_transform` 对真实端元 `SPOIL < 3`、对随机谱 `SPOIL > 6`；MCR-ALS 回收的组分谱与真值相关系数 `> 0.99`、浓度 RMSE `< 0.02`。

---

### B11-1 相移校正 FT/WT 与异核共占位双分支

#### 接口契约
1. `artemis.path_total_phase(datfile, k)`：返回该 FEFF 路径的总相移 $\phi_{\text{tot}}(k)=2\delta_c+\phi_s$（优先读路径 `pha`；否则由 `p.chi` 的复数相位解缠提取），长度与 $k$ 对齐。
2. `athena.ft_to_r_space(..., phase_correct=None)` 与 `larch_wt.continuous_wavelet_transform(..., phase_correct=None)`：在变换核中对 $\chi(k)$ 乘 $e^{-i\phi_{\text{tot}}(k)}$（相位校正后峰位应回到 $R_{\text{eff}}$）。
3. `artemis.cif_hetero_cocupancy_branches(cif_path, absorber, rmax=6.0, delta_z=3, occ_min=0.10)`：对共占位点按元素返回分支 `[{"distance", "element", "occ", "degen"}]`；同壳层各分支共享 $R$、$\sigma^2$（在拟合层由 `shared_bonds` 绑定），简并度按 `N_geom*occ` 分配。

#### 验收标准
- **#54**：合成 Fe–O 单壳：未校正 FT 峰位 $\approx1.47$ Å；相位校正后峰位与 `Reff` 误差 `<0.05 Å`；WT 同理。
- **#55**：合成 `0.7 Ti + 0.3 Ru` 共占位 CIF：`cif_hetero_cocupancy_branches` 返回 Ti/Ru 两支，占据度和为 1，简并度合计等于该壳层几何简并度；单元素近似被显式标记 `single_element_approx=False`。

---

### B11-2 关联 Debye 模型与 BVS 软约束

#### 数学与接口契约
1. `artemis.correlated_debye_sigma2(temperature_K, theta_d, reff, mass_absorber, mass_scatterer, k_debye=None, sigma2_static=0.0)`：
   - 约化质量 $\mu=\frac{m_1m_2}{m_1+m_2}$；
   - 非关联项 $\sigma^2_{\text{un}}(T)=\frac{3\hbar^2}{2\mu k_B\Theta_D}\left[\frac14+\left(\frac{T}{\Theta_D}\right)^2\int_0^{\Theta_D/T}\frac{x}{e^x-1}dx\right]$；
   - 关联因子 $C(R)=1-\frac{\sin(2k_DR)}{2k_DR}$（$k_D$ 缺省由最近邻密度给出，显式记录），$\sigma^2=\sigma^2_{\text{static}}+C(R)\sigma^2_{\text{un}}$；
   - $T\le0$、$\Theta_D\le0$、$\mu\le0$ 抛 `XAFSConfigError`。
2. `artemis.bond_valence_sum(paths, r0_table=None, b=0.37, target=None, sigma=None)`：
   - $\text{BVS}=\sum_i N_i\exp[(R_{0,i}-R_i)/b]$；
   - `r0_table` 默认内置 `BVS_R0`（Fe–O 1.759、Co–O 1.692、Ni–O 1.654、Cu–O 1.679、Pt–Cl 2.19 等）；
   - `fit_r_space(..., bvs_restraint={"target":..., "sigma":..., "r0_table":...})` 把 $(\text{BVS}-\text{target})/\sigma$ 追加到残差；`FitResult.bvs`、`bvs_target`、`bvs_chi2` 输出。
3. `FitResult` 增加 `aic`, `bic`, `bvs`, `bvs_target`, `bvs_chi2`, `fit_space`。

#### 验收标准
- **#56**：$\Theta_D=400$ K、Fe–O：$\sigma^2(300K)>\sigma^2(80K)$；远壳 $R=4$ Å 的 $\sigma^2$ 大于近壳 $R=2$ Å；$T\to0$ 时收敛到静态值；与 Einstein 模型在高温极限同量级（比值 0.5–2）。
- **#57**：Fe–O $N=6,R=2.0$ 的 BVS $\approx3.13$（Fe³⁺）；`bvs_restraint` 开启后拟合输出 `bvs` 且 `|bvs-target|` 小于无约束拟合；`bvs_chi2` 有限。

---

### B11-3 多空间拟合（$R$ / $q$ / WT）

#### 接口契约
`artemis.fit_r_space(..., fit_space="r"|"q"|"wt", wt_sigma=1.5, wt_rmin=None, wt_rmax=None, wt_kmin=None, wt_kmax=None)`：
1. `"r"`：保持历史逐位一致（默认）；
2. `"q"`：数据与模型各自经 `athena.back_ft_to_q_space(r, chi_r, rmin, rmax)` 后，在公共 $q$ 网格比较实/虚部；
3. `"wt"`：数据与模型各自做同参数 CWT，在指定 $(k,R)$ 矩形窗内比较复数 $W$ 的实/虚部；
4. 返回 `FitResult.fit_space`；`r_factor` 的口径随空间切换并在文档中注明（不同空间 R-factor 不可横向比较）。

#### 验收标准
**#58**：同一条合成单壳谱分别用 `"r"`、`"q"`、`"wt"` 拟合，三种空间回收的 $N$、$R$、$\sigma^2$ 相互一致（$|\Delta N|<0.3$、$|\Delta R|<0.02$ Å、$|\Delta\sigma^2|<0.001$ Å²），且默认 `"r"` 的 `r_factor` 与历史基线逐位一致。

---

### B12-1 双吸收边独立 $\Delta E_0$ 与 AIC/BIC

#### 接口契约
1. `fit_multi_datasets(..., share_e0="auto")`：
   - `share_e0=True`：所有数据集共享一个 `("global","e0")`（历史行为）；
   - `share_e0=False`：每个数据集独立 `(di,"e0")`；
   - `"auto"`：按 `datasets[di]["edge_group"]` 分组（缺省用 `edge` 字符串，再缺省全共享）；同组共享一个 $\Delta E_0$。
2. `FitResult.aic = n\ln(\text{RSS}/n)+2k$、`bic = n\ln(\text{RSS}/n)+k\ln n`（$n$ 为残差点数，$k$ 为自由变量数）；`fit_multi_datasets` 返回全局 `aic`/`bic`。
3. `artemis.profile_chi2_contour(fit_func, x_name, y_name, ...)`：二维网格剖面扫描，返回 `{x, y, chi2_grid, levels}`（$1\sigma/2\sigma/3\sigma$），供诊断图使用。

#### 验收标准
**#59**：两条同边谱真值 $\Delta E_0=+2.0$ eV、不同边谱真值 $-6.0$ eV：`share_e0="auto"` 按 `edge_group` 分出 2 个 $\Delta E_0$ 槽位并各自回收误差 `<0.5$ eV`；`share_e0=True` 时只有 1 个槽位且 $R_f$ 显著变差；AIC/BIC 有限且 BIC 对多参数模型惩罚更重。

---

### B12-2 SAC/DAC 团簇生成器与元素扩充

#### 接口契约
1. `artemis.sac_cluster_sites(bond_m_ligand=2.00, ligand="N", n_coord=4, axial=None, axial_bond=2.10, carbon_shell=True, ring_radius=2.90)`：
   - 平面 $M$–$N_4$（正方形，键长 `bond_m_ligand`）；
   - `carbon_shell=True` 时在 $R\approx ring_radius$ 添加 8–12 个 C 原子；
   - `axial` 可给 `"O"/"Cl"` 等，生成轴向配位（oct 五/六配位）；
   - 返回 `[(x,y,z,element), ...]`，中心在原点，可直接送 `build_feff_input`。
2. `artemis.diatomic_m1m2_sites(...)`：双原子位点 $M_1$–$M_2$ 沿 $z$ 轴，各自带 $N$ 配位壳。
3. `absorber.py` 扩充常用催化/环境元素（Cr, Mn, Zn, Mo, Ru, Rh, Pd, Ag, Ir, Au, Ce, U），全部 `s02=None`（未标定，CN 定量拒绝），`e0_edge`/`kmax_ceiling` 填表；**Fe/Co/Ni/Cu/Pt 的 S₀² 逐位不变**。
4. `artemis.COVALENT_RADII` 同步扩充上述元素。

#### 验收标准
**#60**：`sac_cluster_sites` 的 $M$–$N$ 恰好 4 条、键长误差 `<1e-9`；C 壳层存在且距离在 `ring_radius±0.2`；`build_feff_input` 成功生成且 FEFF 可跑；新元素 `calibrated=False`，五元素标定值不变。

---

## 第三阶段防呆自检清单（每完成一批打勾）

- [ ] **B9 批次检查**：
  - [ ] 配置文件 17 键全部生效，未配置键保持默认？
  - [ ] `vary_c4` 有 `c4_dropped_by_nidp` 降级路径？
  - [ ] `_mu.dat` 含 `sigma_mu` 列，`_rdf.csv` 在 C3/C4 时必出？
  - [ ] `.prj` round-trip 与 MCA 坏道剔除保底 2 道？
  - [ ] 新增 #46–#50 后 `Ran 50 tests — OK`？
- [ ] **B10 批次检查**：
  - [ ] MBACK 短窗口 edge_step 误差 <3%？
  - [ ] `chi_std`/`clamp_lo/hi` 透传 larch 且失败时显式降级？
  - [ ] SPOIL 对真实端元 <3、随机谱 >6？
  - [ ] MCR-ALS 相关系数 >0.99、浓度 RMSE <0.02？
  - [ ] 新增 #51–#53 后 `Ran 53 tests — OK`？
- [ ] **B11 批次检查**：
  - [ ] 相移校正后 FT/WT 峰位回到 Reff（误差 <0.05 Å）？
  - [ ] 异核共占位双分支简并度合计等于几何简并度？
  - [ ] 关联 Debye 的 T→0 收敛、远壳 > 近壳、高温与 Einstein 同量级？
  - [ ] BVS 软约束生效且 `bvs_chi2` 有限？
  - [ ] 三种 fit_space 参数一致，默认 `"r"` 逐位不变？
  - [ ] 新增 #54–#58 后 `Ran 58 tests — OK`？
- [ ] **B12 批次检查**：
  - [ ] `share_e0="auto"` 按 `edge_group` 正确分组？
  - [ ] AIC/BIC 有限且 BIC 惩罚更重？
  - [ ] SAC 团簇 M–N 四条、C 壳层在位、FEFF 可跑？
  - [ ] 新元素未标定且五元素 S₀² 不变？
  - [ ] 新增 #59–#60 后 `Ran 60 tests — OK`？
- [ ] **收尾同步**：`README.md` 与 `.trae/skills/xafs-pipeline-report/references/pipeline-reference.md` 更新到 **60 项契约回归**与 B9–B12 接口。

---

## 第三阶段终审验收 Checklist（Phase I + II + III 全量判分，8 条硬指标）

1. **全量 CI**：`python run_ci.py --stress` 中 `py_compile` 退出码 0、`regression_test.py` **60 项全绿**、4 项压力测试退出码 0。
2. **零退化**：默认 `LEPIDO` 三扫描 `r_factor≈0.0273`、第一壳 `N≈6.35`、`R≈2.002 Å`、7 列拟合表逐位不变。
3. **CLI/配置完整**：`--help` 含 Phase I–III 全部开关；YAML/JSON 17 键与 CLI 100% 镜像。
4. **B9 契约**：`vary_c4` 降级、`sigma_mu` 导出、WT 自动标尺、$L_3$ 白线面积、`.prj` round-trip、MCA 坏道剔除全部通过 #46–#50。
5. **B10 契约**：MBACK 短窗口、`chi_std`/`clamp`、SPOIL、MCR-ALS 全部通过 #51–#53。
6. **B11 契约**：相移校正、异核共占位、关联 Debye、BVS、三空间拟合全部通过 #54–#58。
7. **B12 契约**：多边 $\Delta E_0$、AIC/BIC、SAC/DAC、元素扩充全部通过 #59–#60。
8. **物理诚实**：S₀² 五元素值未改；未标定新元素 `calibrated=False`；无法定量时 `cn_reliable=False` 或显式拒绝；不伪造面积/谱点。

---

## 第四阶段（Phase IV：B13–B16）最高优先级执行计划：面向催化化学全面分析收官（待执行，回归 60 → 68 项）

> **战略定位（严禁盲目堆砌）**：
> - 前三阶段（Phase I–III，60 项回归）已完成了 Athena / Artemis / Larch-WT 的底层数学内核建设。
> - 第四阶段（Phase IV）**不堆砌任何与催化化学无关的边缘功能**（不做光学单色器机械标定、不做掠入射薄膜 GIXAFS、不做从头算量子化学接口、不做重型桌面 GUI），**只聚焦催化化学顶刊表征（SAC 单原子、DAC 双原子、双金属固溶体/合金、混合配位化学计量约束、原位 *Operando* 时间序列解谱、存量 Athena 标样库导入）的最后一公里全链闭环**，彻底消灭 Phase III 留下的 5 个“API 已备但未进主产线/CLI”的半连接点。

### 全局执行顺序总览（Phase IV）

| 批次 | 编号 | 催化分析核心模块 | 修改文件 | 新增回归测试（`regression_test.py`） |
|---|---|---|---|---|
| **第十三批** | **B13-1** | SAC / DAC 团簇动态 Motif 构建、CLI 一键建模（`--sac-model` / `--dac-model`）与 $M_1\text{–}M_2$ 金属键存在性 Hamilton F-test 判据 | `structure_lib.py`, `xafs_fit.py`, `xafs_config.py`, `run_xafs.py`, `xafs_report.py` | `test_sac_axial_and_dac_mm_bond_ftest_pipeline` (#61) |
| **第十三批** | **B13-2** | 异核共占位双分支自动拆壳与单数据集内 `path_links` 跨元素共享 $\Delta R/\sigma^2$（`--hetero-cocupancy`） | `artemis.py`, `xafs_fit.py`, `xafs_config.py`, `run_xafs.py`, `xafs_report.py` | `test_hetero_cocupancy_end_to_end_dual_branch_fit` (#62) |
| **第十四批** | **B14-1** | 催化配位化学计量与线性代数约束引擎（Mini-GDS：`--constraints`，支持配位数守恒 $N_A+N_B=C$、比例 $N_B=rN_A$、键长差偏移及误差传播） | `artemis.py`, `xafs_fit.py`, `xafs_config.py`, `run_xafs.py`, `xafs_report.py` | `test_mini_gds_stoichiometric_sum_and_ratio_constraints` (#63) |
| **第十四批** | **B14-2** | 相移校正真键长 FT/WT（`--phase-correct`）与二维 $\chi^2$ 相关性置信椭圆（`--contour`）全链接入诊断图与 HTML 报告 | `xafs_pipeline.py`, `xafs_diagnostics.py`, `xafs_report.py`, `run_xafs.py`, `xafs_config.py` | `test_phase_corrected_diagnostics_and_chi2_contour_report` (#64) |
| **第十四批** | **B14-3** | 关联 Debye 模型（`debye_models`）纳入 `fit_multi_datasets` 全局变温/多壳层联合精修（单 $\Theta_D$ 同时约束近壳与远壳） | `artemis.py` | `test_multi_dataset_correlated_debye_corefinement` (#65) |
| **第十五批** | **B15-1** | 原位（*Operando* / *In-situ*）催化时间序列一键解谱编排器 `run_operando_series` 与 `--operando` 报告模式（PCA + SPOIL 筛标样 + LCF + MCR-ALS + $E_{1/2}/r_1$ 轨迹） | `xafs_lcf.py`, `xafs_diagnostics.py`, `xafs_report.py`, `xafs_config.py` | `test_operando_series_pca_spoil_lcf_and_mcr_tracking` (#66) |
| **第十六批** | **B16-1** | 存量 Athena 原生 Perl `.prj` 只读导入器 `read_athena_perl_prj`（直接读通仓库 `Fe_standards.prj`，兼容 #50） | `xafs_io.py` | `test_read_real_athena_perl_prj_fe_standards` (#67) |
| **第十六批** | **B16-2** | Co / Ni / Cu / Pt 四大催化元素内置结构库扩充（$\text{M–N}_4$ 单原子、$\text{Co}_3\text{O}_4$ 尖晶石、$\text{Cu}_2\text{O}$、$\text{PtCl}_4$ 平方平面） | `_structure_lib/build_seed.py`, `_structure_lib/{Co,Ni,Cu,Pt}/` | `test_catalysis_library_motifs_for_co_ni_cu_pt` (#68) |

---

## 第十三批（B13）：催化活性位点（SAC / DAC / 异核双金属固溶体）一键建模与拟合闭环

### B13-1 SAC / DAC 团簇动态 Motif 构建与 $M_1\text{–}M_2$ 金属键 Hamilton F-test 判据

#### 1. 催化痛点与数学契约
- **痛点**：B12-2 在 `artemis.py:823, 878` 实现了 `sac_cluster_sites` 与 `diatomic_m1m2_sites`，但未包装成 `structure_lib.Motif`，导致用户无法在 CLI（`run_xafs.py` / `xafs_report.py`）对未知的 $\text{M–N}_4\text{O}_{\text{axial}}\text{/C}$ 或 $\text{M}_1\text{–M}_2\text{/N–C}$ 样品直接跑 S7 多壳层拟合并判断“究竟是纯单原子还是存在双原子金属键”。
- **接口与行为契约**：
  1. 在 `structure_lib.py` 新增两个动态 Motif 构造函数（缓存落盘到 `_structure_lib/_dynamic_cache/<hash>/`）：
     ```python
     def build_sac_motif(
         absorber: str,
         ligand: str = "N",
         n_coord: int = 4,
         bond_m_ligand: float = 2.00,
         axial: str | None = None,
         axial_bond: float = 2.10,
         carbon_shell: bool = True,
         ring_radius: float = 2.90,
     ) -> Motif

     def build_dac_motif(
         absorber: str,
         partner: str,
         bond_length: float = 2.50,
         ligand: str = "N",
         n_coord: int = 4,
         ligand_bond: float = 2.00,
     ) -> Motif
     ```
     - 内部调用 `artemis.sac_cluster_sites` / `artemis.diatomic_m1m2_sites` 生成 `sites`（`[(x, y, z, element), ...]`），构造 `Motif(id=..., absorber=absorber, kind="cluster", nuclearity=("mononuclear" if sac else "binuclear"), sites=sites, shells=..., r_applicable=(1.0, 3.6))`，并确保 `ensure_feff(motif, edge=...)` 可直接在其缓存目录跑通 `feff6l`。
     - 支持精简字符串规格解析函数 `structure_lib.parse_catalytic_cluster_spec(absorber: str, spec: str | dict) -> Motif`：
       - SAC 示例：`"sac:N4@2.00,axial=O@2.10,C@2.90"` 或 `{"type": "sac", "ligand": "N", "n_coord": 4, "bond_m_ligand": 2.00, "axial": "O", "axial_bond": 2.10, "carbon_shell": True, "ring_radius": 2.90}`；
       - DAC 示例：`"dac:Co@2.50,N4@2.00"` 或 `{"type": "dac", "partner": "Co", "bond_length": 2.50, "ligand": "N", "n_coord": 4, "ligand_bond": 2.00}`。
  2. 在 `xafs_fit.fit_identified` 与 `xafs_fit.run` 新增参数 `catalytic_model: str | dict | Motif | None = None`（CLI 对应 `--catalytic-model`，`xafs_config.DEFAULTS["catalytic_model"] = None`）：
     - 当 `catalytic_model` 非空时，直接用 `parse_catalytic_cluster_spec` 构造目标 `Motif`，绕过晶相库搜索，直接进入 S7 多壳层精修；
     - **DAC 金属键显著性自动判据（Hamilton F-test）**：若目标 Motif 为 `nuclearity == "binuclear"`（含 $M_1\text{–}M_2$ 金属键路径），`fit_identified` 在完成含 $M_1\text{–}M_2$ 路径的完整拟合 `fit_complex` 后，自动剔除 $M_1\text{–}M_2$ 路径重跑一次纯配体壳层拟合 `fit_simple`（同一 $k/R$ 窗口），调用 `artemis.shell_f_test(fit_simple, fit_complex)`，在返回字典顶层写入：
       ```python
       "mm_bond_ftest": {
           "partner": partner_symbol,
           "f_stat": float,
           "p_value": float,
           "significant": bool,   # p_value < 0.05 且 M-M 的 N > N_err
           "delta_r_factor": float,
       }
       ```

#### 2. 验收标准（新增回归 #61：`test_sac_axial_and_dac_mm_bond_ftest_pipeline`）
- 用 `build_dac_motif("Fe", "Co", bond_length=2.50, ligand="N", n_coord=4, ligand_bond=2.00)` 的 FEFF 路径合成两条谱（均加 $\varepsilon \approx 10^{-3}$ 噪声）：
  1. **纯 SAC 谱**（仅含 $\text{Fe–N}\times 4$ @ 2.00 Å）：以 `catalytic_model="dac:Co@2.50,N4@2.00"` 跑 `xafs_fit.fit_identified`，验证 `out["mm_bond_ftest"]["significant"] is False`（`p_value > 0.05`）；以 `"sac:N4@2.00,axial=O@2.10,C=none"` 拟合含轴向氧合成谱，回收 $N_{\text{N}} \approx 4, N_{\text{O}} \approx 1$（误差 `<0.5`）。
  2. **真 DAC 谱**（$\text{Fe–N}\times 4$ @ 2.00 Å + $\text{Fe–Co}\times 1$ @ 2.50 Å）：以 `catalytic_model="dac:Co@2.50,N4@2.00"` 跑 `fit_identified`，验证回收 $R_{\text{Fe-Co}} \in [2.47, 2.53]$ Å、$N_{\text{Fe-Co}} \in [0.7, 1.3]$，且 `out["mm_bond_ftest"]["significant"] is True`（`p_value < 0.01`）。

---

### B13-2 异核共占位双分支自动拆壳与跨元素受约束精修（`--hetero-cocupancy`）

#### 1. 催化痛点与数学契约
- **痛点**：B11-1 的 `artemis.cif_hetero_cocupancy_branches`（`artemis.py:384`）已经能从共占位 CIF 算准双元素分支与简并度，但 `artemis.fit_r_space` 缺少单数据集内的路径间参数绑定 `path_links`，`xafs_fit.select_shell_paths` 也还未调用该函数拆出双分支 FEFF 路径。
- **数学与接口契约**：
  1. 在 `artemis.fit_r_space`（`artemis.py:1125`）新增可选参数 `path_links: list[dict] | None = None`（默认 `None`）：
     - 每项格式：`{"target": int, "source": int, "share_deltar": bool = True, "share_sigma2": bool = True, "degen_ratio": float | None = None}`；
     - 当 `target` 路径绑定到 `source` 路径时：
       - 若 `share_deltar=True`，`target` 不分配 `(target, "deltar")` 自由变量槽位，在 `apply_params` 中令 $\Delta R_{\text{target}} = \Delta R_{\text{source}}$，误差 $\sigma(\Delta R_{\text{target}}) = \sigma(\Delta R_{\text{source}})$；
       - 若 `share_sigma2=True`，`target` 不分配 `(target, "sigma2")` 自由变量槽位，令 $\sigma^2_{\text{target}} = \sigma^2_{\text{source}}$；
       - 若 `degen_ratio is not None`，`target` 不分配 `(target, "n")` 槽位，令 $N_{\text{target}} = \text{degen\_ratio} \cdot N_{\text{source}}$；若 `degen_ratio is None` 且 `vary_n=True`，两分支各自独立拟合 $N_A, N_B$（用于测定固溶体局部真实掺杂比！）。
  2. 在 `xafs_fit.select_shell_paths` 与 `fit_identified` 新增 `hetero_cocupancy: bool = False`（CLI `--hetero-cocupancy`，`xafs_config.DEFAULTS["hetero_cocupancy"] = False`）：
     - 当 `hetero_cocupancy=True` 且 `motif.cif_file` 存在时，调用 `artemis.cif_hetero_cocupancy_branches`；
     - 对 `single_element_approx == False` 的共占位壳层，为每个分支元素生成对应的单散射 `feffNNNN.dat`（在 `_structure_lib/_hetero_cache/` 下将共占位位点纯化为该分支元素跑一次 FEFF 缓存，取该距离处的 SS 路径），在 `plan` 中输出两条路径（如 `Fe-Ti@3.55`、`Fe-Ru@3.55`），初始简并度分别取 `branch["degen"]`，并自动生成 `path_links=[{"target": idx_b, "source": idx_a, "share_deltar": True, "share_sigma2": True, "degen_ratio": None}]` 传入 `fit_r_space`。

#### 2. 验收标准（新增回归 #62：`test_hetero_cocupancy_end_to_end_dual_branch_fit`）
- 构造吸收原子 Fe、近邻壳层为 $0.7\text{ Ti} + 0.3\text{ Ru}$ 共占位的 CIF（$N_{\text{geom}} = 8$，真值 $N_{\text{Ti}} = 5.6, N_{\text{Ru}} = 2.4, \Delta R = -0.02\text{ Å}, \sigma^2 = 0.004\text{ Å}^2$）并合成 $\chi(k)$：
  1. 验证 `fit_r_space(..., path_links=...)` 将双分支自由变量数从 $2\times 3 + 1 = 7$ 压到 $4 + 1 = 5$（共享 $\Delta R$ 与 $\sigma^2$，各自放开 $N_{\text{Ti}}, N_{\text{Ru}}$）；
  2. 回收 $N_{\text{Ti}} \in [5.2, 6.0]$、$N_{\text{Ru}} \in [2.0, 2.8]$、两路径 $\Delta R$ 与 $\sigma^2$ 逐位相等，且双分支拟合 `r_factor < 1e-4`，比单元素近似（仅用 Ti 一条路径拟合同一谱）的 `r_factor` 小 5 倍以上。

---

## 第十四批（B14）：催化配位化学计量约束（Mini-GDS）与审稿级相移/相关性诊断闭环

### B14-1 催化配位化学计量与线性代数约束引擎（Mini-GDS：`constraints`）

#### 1. 催化痛点与数学契约
- **痛点**：催化样品常为混合配体（如氧氮共配位 $\text{M–N}_x\text{O}_{6-x}$、氯氧共配位 $\text{Pt–O}_x\text{Cl}_{4-x}$）。若同时放开 $N_A, N_B, \sigma^2_A, \sigma^2_B$，振幅完全简并；必须像 Artemis 的 GDS 一样施加**总配位数守恒** $N_A + N_B = N_{\text{tot}}$ 或**固定化学计量比** $N_B = r \cdot N_A$，并严格传播协方差形式误差。
- **接口与数学契约**：
  1. 在 `artemis.py` 新增解析与校验函数：
     ```python
     def normalize_constraints(
         constraints: list[dict | str] | str | None,
         path_labels: list[str],
     ) -> list[dict]
     ```
     支持两种等价输入格式（绝对禁止使用 `eval` / `exec`）：
     - **结构化字典列表**（供 API 与 YAML/JSON 配置使用）：
       - 配位数守恒：`{"kind": "n_sum", "target": "Fe-Cl", "source": "Fe-O", "total": 6.0, "coeff": 1.0}` $\implies N_{\text{target}} = \text{total} - \text{coeff} \cdot N_{\text{source}}$；
       - 配位数比例：`{"kind": "n_ratio", "target": "Fe-Cl", "source": "Fe-O", "ratio": 1.5}` $\implies N_{\text{target}} = \text{ratio} \cdot N_{\text{source}}$；
       - 键长差偏移：`{"kind": "deltar_offset", "target": "Fe-Cl", "source": "Fe-O", "offset": 0.0}` $\implies \Delta R_{\text{target}} = \Delta R_{\text{source}} + \text{offset}$；
       - 参数边界覆盖：`{"kind": "bounds", "path": "Fe-O", "param": "sigma2", "min": 0.002, "max": 0.010}`。
     - **简洁字符串表达式**（供 CLI `--constraints` 使用，分号或逗号分隔）：
       - `"N(Fe-O)+N(Fe-Cl)=6.0"` $\to$ `kind="n_sum"`；
       - `"N(Fe-Cl)=1.5*N(Fe-O)"` $\to$ `kind="n_ratio"`；
       - `"dR(Fe-Cl)=dR(Fe-O)+0.02"` $\to$ `kind="deltar_offset"`。
  2. 在 `artemis.fit_r_space`（`artemis.py:1125`）与 `xafs_fit.fit_identified` / `run` 新增参数 `constraints=None`（`xafs_config.DEFAULTS["constraints"] = None`）：
     - 被 `n_sum` / `n_ratio` 约束的 `target` 路径**不分配独立 `(target, "n")` 自由变量槽位**（`n_variables` 真实减少 1）；
     - 为保证物理非负性（$N_{\text{target}} \ge 0$），当 `kind == "n_sum"` 时，自动将 `source` 路径的 $N_{\text{source}}$ 上界收紧为 `min(upper[source], total / coeff)`；
     - **误差传播**：
       - `n_sum`：$\sigma(N_{\text{target}}) = |\text{coeff}| \cdot \sigma(N_{\text{source}})$；
       - `n_ratio`：$\sigma(N_{\text{target}}) = |\text{ratio}| \cdot \sigma(N_{\text{source}})$；
       - `deltar_offset`：$\sigma(\Delta R_{\text{target}}) = \sigma(\Delta R_{\text{source}})$。
     - 检测循环依赖（如 $A \to B \to A$）或未知路径标签时，显式抛 `XAFSConfigError`。

#### 2. 验收标准（新增回归 #63：`test_mini_gds_stoichiometric_sum_and_ratio_constraints`）
- 用 Fe–O（$R_{\text{eff}} \approx 1.98$ Å）与 Fe–Cl（或第二条独立路径）合成总配位数为 6.0 的双路径谱（真值 $N_1 = 2.4, N_2 = 3.6$）：
  1. 传入 `constraints="N(Fe-O)+N(Fe-Cl)=6.0"`（或字典形式），验证 `fit.n_variables` 比无约束减少 1；
  2. 验证拟合结果 $N_1 + N_2 = 6.0$（误差 `<1e-12`），$|N_1 - 2.4| < 0.25$，且 `fit.paths[1].n_err == fit.paths[0].n_err`（误差严格传递）；
  3. 验证非法循环约束 `"N(Fe-O)=1.0*N(Fe-Cl);N(Fe-Cl)=1.0*N(Fe-O)"` 抛出 `XAFSConfigError`。

---

### B14-2 相移校正真键长 FT/WT 与二维 $\chi^2$ 置信椭圆全链接入诊断与 HTML 报告

#### 1. 催化痛点与接口契约
- **痛点**：B11-1 的 `phase_correct` 和 B12-1 的 `profile_chi2_contour` 已在底层跑通，但 `xafs_diagnostics.py`、`run_xafs.py`、`xafs_report.py` 尚未接入，用户跑报告看不到相移校正后的真键长图和参数相关性等高线图。
- **接口与落盘契约**：
  1. 在 `xafs_config.DEFAULTS` 新增 `"phase_correct": False` 与 `"chi2_contour": False`；CLI `run_xafs.py` 与 `xafs_report.py` 同步新增 `--phase-correct`（布尔开关）与 `--chi2-contour`（布尔开关）。
  2. 在 `xafs_diagnostics.py` 新增两个标准诊断产物函数：
     - `write_phase_corrected_diagnostics(res, fit_out, outdir, stem=None) -> dict`：
       - 从 `fit_out` 的第一壳路径（若无 S7 则从 L1 第一壳缓存路径）调 `artemis.path_total_phase(datfile, res.chi_k.k)` 获取 $\phi_{\text{tot}}(k)$；
       - 调用 `athena.ft_to_r_space(..., phase_correct=phi)` 与 `larch_wt.continuous_wavelet_transform(..., phase_correct=phi)`；
       - 导出 `<stem>_ft_phase_corrected.csv`（列：`R_A, mag_raw, real_raw, mag_phase_corrected, real_phase_corrected`）与 `<stem>_phase_corrected.png`（左图：未校正 vs 相移校正 $|\chi(R)|$ 并标注拟合真键长 $R_1$ 竖线；右图：相移校正二维 WT 图），返回 `{"csv": ..., "png": ..., "r_peak_raw": float, "r_peak_corrected": float, "r_fit": float}`。
     - `write_chi2_contour_diagnostics(res, fit_out, outdir, stem=None, n_grid: int = 15) -> dict | None`：
       - 从 `fit_out["fit"]["fit"]` 取第一壳的 `(N, σ²)`（或相关矩阵中非对角 $|\rho|$ 最大的一对自由变量），围绕最优解 $\pm 3\sigma$（若 $\sigma$ 过小则取 $\pm 20\%$）构造 `n_grid × n_grid` 网格，调用 `artemis.profile_chi2_contour` 计算二维 $\chi^2$ 曲面；
       - 导出 `<stem>_chi2_contour.csv` 与 `<stem>_chi2_contour.png`（含 $1\sigma/2\sigma/3\sigma$ 即 $\Delta\chi^2 = 2.30, 6.18, 11.83$ 等高线及最优解星标），返回 `{"csv": ..., "png": ..., "x_name": ..., "y_name": ..., "min_chi2": float}`。
  3. 当 `--phase-correct` 或 `--chi2-contour` 开启时，`xafs_report.py` 自动调用上述函数并将生成的 PNG base64 内联到 `<stem>_report.html` 的诊断章节中；默认关闭时不生成，保持默认耗时与产物不变。

#### 2. 验收标准（新增回归 #64：`test_phase_corrected_diagnostics_and_chi2_contour_report`）
- 对 `self.lep` 与 `self.lep_fit` 调用 `write_phase_corrected_diagnostics` 与 `write_chi2_contour_diagnostics`：
  1. 验证相移校正前 `r_peak_raw` 在 $1.4\text{–}1.7$ Å，校正后 `|r_peak_corrected - r_fit| < 0.06` Å，且 CSV 与 PNG 文件存在且非空；
  2. 验证 `write_chi2_contour_diagnostics` 返回的 `min_chi2` 有限、网格形状为 `(15, 15)`、CSV 与 PNG 存在；通过 `xafs_report.py` 开启两开关生成的 HTML 报告中包含对应图像标题。

---

### B14-3 关联 Debye 模型纳入多数据集/变温全局联合拟合（`debye_models`）

#### 1. 数学与接口契约
- 在 `artemis.fit_multi_datasets`（`artemis.py:2085`）新增参数 `debye_models: list[dict] | None = None`：
  - 每项格式：`{"paths": [(ds_idx, path_idx), ...], "mass_absorber": float, "mass_scatterer": float, "theta_d_init": 400.0, "k_debye": None}`；
  - 每个 Debye 组分配两个全局变量槽位 `("debye", gi, "sigma2_static")`（初值 `0.001`，边界 `[0.0, 0.015]`）与 `("debye", gi, "theta_d")`（初值 `theta_d_init`，边界 `[80.0, 1200.0]`）；
  - 被纳入 `debye_models` 的 `(ds_idx, path_idx)` 不再分配独立 `(di, pi, "sigma2")` 槽位（且不可与 `einstein_models` 或 `shared_bonds.share_sigma2` 重复绑定同一路径，冲突时抛 `XAFSConfigError`）；
  - 在 `resolved_sigma2(x, key)` 中：根据该数据集的 `temperature_K` 与该路径的晶体学 `reff = float(ds_paths[di][pi].reff)`，调用 `artemis.correlated_debye_sigma2(T, theta_d, reff, mass_absorber, mass_scatterer, k_debye=k_debye, sigma2_static=s2_static)` 实时计算 $\sigma^2(T_d, R_{\text{eff}})$，并用有限差分（步长 `1e-4` 和 `0.5 K`）传递 `sigma2_err`；
  - 返回字典新增 `"debye_results": [{"theta_d", "theta_d_err", "sigma2_static", "sigma2_static_err", "paths": ...}]`。

#### 2. 验收标准（新增回归 #65：`test_multi_dataset_correlated_debye_corefinement`）
- 取 `lepidocrocite` 第一壳 Fe–O（$R_{\text{eff}} \approx 1.98$ Å）与一远壳路径（或人工设定 $R_{\text{eff}} \approx 3.68$ Å 的 Fe–O 路径），用真值 $\Theta_D = 420$ K、$\sigma^2_{\text{static}} = 0.0012\text{ Å}^2$ 合成 80 K 与 350 K 两个温度下的数据集：
  1. 用单个 `debye_models` 组同时绑定两温度下的路径，验证 `theta_d` 回收误差 `<15` K、`sigma2_static` 误差 `<5e-4`；
  2. 验证同一温度下远壳拟合出的 $\sigma^2$ 严格大于近壳 $\sigma^2$，且 `global_r_factor < 1e-4`。

---

## 第十五批（B15）：原位（*Operando* / *In-situ*）催化时间序列一键解谱与演变报告

### B15-1 原位序列谱编排器 `xafs_lcf.run_operando_series` 与 `--operando` 报告模式

#### 1. 催化痛点与数学契约
- **痛点**：B10-2 实现了 `run_pca`、`target_transform`（SPOIL）与 `mcr_als`，但催化研究者面对一组随电位/时间变化的原位谱（如 `scan_01.xdi ... scan_10.xdi`）时，需要一次运行自动完成“主成分定物种数 $\to$ SPOIL 筛除伪标样 $\to$ LCF 标样拟合 + MCR-ALS 无标样盲解 $\to$ 价态/键长与物种分数同步轨迹图”。
- **接口与落盘契约**：
  1. 在 `xafs_lcf.py` 新增函数：
     ```python
     def run_operando_series(
         items: list,
         standards: dict | None = None,
         n_components: int | None = None,
         mode: str = "xanes",
         xmin: float = -20.0,
         xmax: float = 60.0,
         closure: bool = True,
         max_iter: int = 800,
         spoil_good: float = 3.0,
         spoil_poor: float = 6.0,
     ) -> dict
     ```
     - `items`：可为 `[(name, NormalizedScan), ...]` 或 `[PipelineResult, ...]` 或 `[(name, x, y), ...]`（至少 3 条谱，否则抛 `XAFSConfigError`）；
     - 步骤 1（PCA）：调用 `run_pca(spectra, xmin, xmax)`；若 `n_components is None`，取 `k = max(2, min(int(pca["n_components"]), len(items) - 1))`；
     - 步骤 2（Target Transform + LCF，当 `standards` 非空时）：对所有候选标样调 `target_transform(pca, cand_list, n_components=k, spoil_good=spoil_good, spoil_poor=spoil_poor)`；仅保留 `rating != "unacceptable"`（即 `spoils < 2 * spoil_poor`）的标样作为有效标样集 `accepted_standards`（若全部被拒则退回全部标样并记警告），对序列每条谱调 `run_lcf` 得到浓度矩阵 `lcf_weights`（形状 `(len(accepted_standards), len(items))`）；
     - 步骤 3（MCR-ALS 盲解谱）：调用 `mcr_als(spectra, n_components=k, max_iter=max_iter, closure=closure)`，并且**对 MCR-ALS 解出的组分按特征能量（主峰位置或重心）从低到高排序**，消除组分排列随机性；
     - 步骤 4（描述符同步轨迹）：若输入为 `PipelineResult` 或 `NormalizedScan`，提取每条谱的 `e0`、`e_half_height`、`white_line_height`、`ft_r1`（若有 `res.r1`）组装为 `descriptors: list[dict]`；
     - 返回 `{"names": [...], "n_components": k, "pca": pca, "target_transform": tt_out, "accepted_standards": [...], "rejected_standards": [...], "lcf_series": [...], "mcr_als": mcr_out, "descriptors": [...]}`。
  2. 在 `xafs_diagnostics.py` 新增 `write_operando_diagnostics(operando_out: dict, outdir, stem: str = "operando") -> dict`：
     - 写出 `<stem>_operando_concentrations.csv`（含每条谱的 MCR 组分分数、LCF 标样权重与结构/价态描述符）、`<stem>_operando_components.csv`（能量网格与 MCR 纯组分谱）及四宫格演化图 `<stem>_operando_summary.png`；
  3. 在 `xafs_report.py` 与 `xafs_config.DEFAULTS` 新增 `--operando`（布尔开关，`DEFAULTS["operando"] = False`）：当开启且输入样品数 $\ge 3$ 时，在生成各单样品报告后自动调用 `run_operando_series` 与 `write_operando_diagnostics`，并将原位演变总表与四宫格图写入汇总报告。

#### 2. 验收标准（新增回归 #66：`test_operando_series_pca_spoil_lcf_and_mcr_tracking`）
- 用真实 `LEPIDO` 与 `GOETHITE` 的归一化 XANES 谱按权重 $w_i \in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]$ 线性混合构造 6 条原位相变序列（加 $10^{-3}$ 噪声），候选标样传入 `{"lepido": lep, "goethite": goe, "fake_noise": random_scan}`：
  1. 验证 `run_operando_series` 的 `rejected_standards == ["fake_noise"]`、`accepted_standards == ["lepido", "goethite"]`；
  2. 验证 `lcf_series` 回收的 6 个点 `goethite` 权重与真值 `[0.0..1.0]` 的 RMSE `< 0.02`；
  3. 验证 `mcr_als` 回收的 2 条组分浓度随序号呈单调变化（相关系数 `> 0.98`），且 `write_operando_diagnostics` 成功导出 CSV 与 PNG。

---

## 第十六批（B16）：存量 Athena 原生 `.prj` 只读导入与多元素催化结构库扩充

### B16-1 存量 Athena 原生 Perl `.prj` 只读导入器（`read_athena_perl_prj`）

#### 1. 痛点、防呆陷阱与接口契约
- **痛点**：课题组历史标样（包括本仓库自带的 `data/练手数据/02_Fe_standards/Fe_standards.prj`）是 Athena 原生 Perl 格式（gzip 压缩的文本，每个数据组以 `[record]` 分隔，内含 `@args = (...)`、`@x = (...)`、`@y = (...)`、`@i0 = (...)`、`@signal = (...)`）。B9-3 的 `read_prj` 遇到非 JSON 直接抛错，导致连仓库自带的 `Fe_standards.prj` 都无法导入。
- **核心防呆陷阱（必读，否则破坏 #50 回归）**：
  - 既有回归 #50（`regression_test.py:2169`）向临时文件写入了 `"$VAR1 = { 'x' => 1 };"` 并断言 `with self.assertRaises(XAFSConfigError): xafs_io.read_prj(perl)`！
  - 因此，**判定为合法 Athena 原生 `.prj` 的充要条件**必须是：文本中同时包含 `@args =` 与至少一组非空的 `@x =` 和 `@y =` 数组记录！对于像 `"$VAR1 = { 'x' => 1 };"` 这种不含 `@x =` / `@y =` 扫描记录的任意 Perl 片段，**必须继续抛出 `XAFSConfigError`**！
- **解析算法与接口契约**：
  1. 在 `xafs_io.py` 新增 `read_athena_perl_prj(path_or_text, as_scans: bool = False, include_derived: bool = False) -> dict`：
     - 按 `[record]` 切分记录块；用正则提取 `@args = (...)` 中的键值序列（偶数下标为键、奇数下标为字符串值，遇到嵌套 `{...}` / `[...]` 时安全跳过或按平衡括号提取字符串键值），以及 `@x = (...)`、`@y = (...)`、`@i0 = (...)`、`@signal = (...)` 中的浮点数列表；
     - 对每条记录：
       - 读取 `label`（若无则取 `tag` 或 `old_group`）、`bkg_z`（吸收元素）、`bkg_e0`、`is_merge`、`is_nor`、`ln`；
       - 默认（`include_derived=False`）跳过参比通道派生记录（`label` 去空格后以 `"Ref "` 开头或 `tag` 以 `"d_Ref_"` 开头）以及合并派生记录（`str(is_merge) not in ("", "0")`），只保留原始实验扫描；当 `include_derived=True` 时全部保留；
       - 吸收系数 $\mu(E)$ 还原规则：
         - 在 Athena `.prj` 中，当 `@i0` 与 `@signal` 长度均等于 `@x` 且非全零时，Athena 的 `@y` 存的是已做过内部预处理或对数符号相反的数组（实测 `Fe_standards.prj` 第二条记录 `@y` 为负值 $-\ln(I_0/I_t)$）；此时若 `str(args.get("ln", "1")) == "1"` 且 $I_0/I_{\text{signal}} > 0$，优先按标准透射公式 $\mu = \ln(I_0 / I_{\text{signal}})$ 计算（并同时把 `i0`、`it=signal` 填入 `RawScan`）；若 `@i0` 缺失或 `is_nor == "1"`，则直接取 `@y` 作为 $\mu(E)$；
     - 返回与 `read_prj` 完全同构的字典：`{"xafs_prj": "athena_perl_v1", "scans": [...], "standards": {}, "meta": {"source_format": "athena_perl", "n_records": ...}}`。
  2. 在 `xafs_io.read_prj(path, as_scans=False, include_derived=False)` 中：当文本不是 JSON 时，若满足上述合法 Athena 记录条件（含 `@args` 且含 `@x`、`@y`），转发给 `read_athena_perl_prj`；否则抛 `XAFSConfigError`。

#### 2. 验收标准（新增回归 #67：`test_read_real_athena_perl_prj_fe_standards`）
- 直接调用 `out = xafs_io.read_prj( Path("data/练手数据/02_Fe_standards/Fe_standards.prj"), as_scans=True)`：
  1. 验证 `out["meta"]["source_format"] == "athena_perl"`，且解析出的原始扫描包含 `Fe_lepidocrocite.000`、`Fe_lepidocrocite.100`、`Fe_lepidocrocite.200`、`Fe_pyrite.000`（共 4 条主扫描，`Ref` 与 `merge` 默认被过滤；传 `include_derived=True` 时为 10 条）；
  2. 取其中 `Fe_lepidocrocite.100`（或 `.000`）的 `energy` 与 `mu`，与直接用 `xafs_io.read_scan` 读取同目录下 `Fe_lepidocrocite.100` 原始文件的 `energy` 与 `mu` 对比，长度完全相等（462 点）且 `np.max(np.abs(s_prj.energy - s_raw.energy)) < 1e-4`、`np.corrcoef(s_prj.mu, s_raw.mu)[0, 1] > 0.999999`；
  3. 同时确认 #50 的 `"$VAR1 = { 'x' => 1 };"` 依然抛 `XAFSConfigError`。

---

### B16-2 Co / Ni / Cu / Pt 四大催化元素内置结构库扩充

#### 1. 催化痛点与接口契约
- **痛点**：目前 `_structure_lib/` 中只有 Fe 拥有丰富的催化 Motif（氧化物、单原子 $\text{Fe–N}_4$、双核 $\text{Fe}_2\text{N}_6$），而另外 4 个已完成 $S_0^2$ 标定的催化核心元素（Co, Ni, Cu, Pt）仅有金属单质（及 NiO），遇到单原子催化剂或典型配合物（如 `Pt-L3-PlatinumChloride.xdi`）只能退化走 L1。
- **扩充清单与防呆契约**：
  1. 在 `_structure_lib/build_seed.py`（并运行落盘到 `_structure_lib/{Co,Ni,Cu,Pt}/`）新增以下催化基准 Motif：
     - **Co**：`con4_planar`（平面四配位单原子 $\text{Co–N}_4$ + 碳壳层，键长 $1.96$ Å / $2.88$ Å，`kind="cluster"`）；
     - **Ni**：`nin4_planar`（平面四配位单原子 $\text{Ni–N}_4$ + 碳壳层，键长 $1.92$ Å / $2.85$ Å，`kind="cluster"`）；
     - **Cu**：`cun4_planar`（平面四配位单原子 $\text{Cu–N}_4$ + 碳壳层，键长 $1.95$ Å / $2.88$ Å，`kind="cluster"`）与 `cu2o_cuprite`（立方 $\text{Cu}_2\text{O}$，空间群 $Pn\bar{3}m, a=4.2696$ Å，第一壳 $\text{Cu–O}\times 2 @ 1.849$ Å、第二壳 $\text{Cu–Cu}\times 12 @ 3.019$ Å）；
     - **Pt**：`ptcl4_square`（平面四配位 $[\text{PtCl}_4]^{2-}$ 团簇，第一壳 $\text{Pt–Cl}\times 4 @ 2.31$ Å，`kind="cluster"`，对应仓库真实标样 `Pt-L3-PlatinumChloride.xdi`）。
  2. **压力测试零退化铁律**：上述新增 Motif 的第一壳距离（$1.85\text{–}2.31$ Å）与 Co/Ni/Cu/Pt 金属箔的第一壳（$2.49\text{–}2.77$ Å）及 NiO 的主导峰（$2.95$ Å）完全错开。落盘后必须立即运行 `stress_test_{co,ni,cu,nio}.py`，确保金属箔与 NiO 在 P4 结构识别中依然稳居 `top1`（`score < 0.02`）且与新增 Motif 的 `margin > 0.40`！

#### 2. 验收标准（新增回归 #68：`test_catalysis_library_motifs_for_co_ni_cu_pt`）
1. 验证 `structure_lib.load_library` 对 `"Co"`, `"Ni"`, `"Cu"`, `"Pt"` 分别包含新增的 `con4_planar`、`nin4_planar`、`cun4_planar`、`cu2o_cuprite`、`ptcl4_square`，且每个 Motif 调用 `structure_lib.ensure_feff(m, edge=get_absorber(el).edge)` 均能生成有效 `feff0001.dat`；
2. 对真实实验谱 `PT_CL`（`Pt-L3-PlatinumChloride.xdi`）运行 `xafs_fit.run(res_ptcl, identify=True, fit=True, absorber_symbol="Pt")`，验证 `structure["ranked"][0]["id"] == "ptcl4_square"`，且 S7 拟合回收第一壳 $N \in [3.5, 4.5]$、$R \in [2.26, 2.35]$ Å、`r_factor < 0.05`！

---

## 第四阶段防呆自检清单（小模型每完成一批必须逐项打勾）

- [ ] **B13 批次检查（SAC/DAC 与异核共占位）**：
  - [ ] `build_sac_motif` / `build_dac_motif` 缓存目录固定在 `_structure_lib/_dynamic_cache/`，不会在临时目录退出后失效？
  - [ ] `--catalytic-model` 同时接入 `run_xafs.py`、`xafs_report.py` 与 `xafs_config.DEFAULTS`？
  - [ ] DAC 拟合自动输出 `mm_bond_ftest`，且纯 SAC 上 `significant=False`、真 DAC 上 `significant=True`？
  - [ ] `fit_r_space(..., path_links=...)` 正确减少自由变量数并同步传递 `delta_r_err`/`sigma2_err`？
  - [ ] 新增 #61–#62 后运行 `python regression_test.py` 显示 `Ran 62 tests — OK`？
- [ ] **B14 批次检查（Mini-GDS 约束、相移/等高线诊断、关联 Debye 联合拟合）**：
  - [ ] `normalize_constraints` 绝不使用 `eval`/`exec`，循环约束或未知路径抛 `XAFSConfigError`？
  - [ ] `n_sum` 约束后 $N_A + N_B = C$ 精确成立且 $\sigma(N_B) = |c|\cdot\sigma(N_A)$？
  - [ ] `--phase-correct` 与 `--chi2-contour` 默认关闭，开启时生成 CSV + PNG 并内联进 HTML 报告？
  - [ ] `fit_multi_datasets(..., debye_models=...)` 在单 $\Theta_D$ 下同时约束近壳与远壳 $\sigma^2(T, R)$？
  - [ ] 新增 #63–#65 后运行 `python regression_test.py` 显示 `Ran 65 tests — OK`？
- [ ] **B15 批次检查（原位 Operando 时间序列）**：
  - [ ] `run_operando_series` 自动用 SPOIL 过滤 `unacceptable` 伪标样？
  - [ ] MCR-ALS 解出的纯组分谱已按特征能量排序，不会因组分顺序翻转导致测试抖动？
  - [ ] `--operando` 模式生成 `_operando_concentrations.csv`、`_operando_components.csv` 与 `_operando_summary.png`？
  - [ ] 新增 #66 后运行 `python regression_test.py` 显示 `Ran 66 tests — OK`？
- [ ] **B16 批次检查（原生 Athena `.prj` 与多元素催化库）**：
  - [ ] `read_prj` 对仓库真实 `Fe_standards.prj` 解析出 4 条主扫描且与原始文件一致，同时对 `#50` 的 `"$VAR1 = { 'x' => 1 };"` 依然抛 `XAFSConfigError`？
  - [ ] 新增 Co/Ni/Cu/Pt 催化 Motif 后，`stress_test_{co,ni,cu,nio}.py` 4 项压力测试全部 `exit=0`？
  - [ ] 真实实验谱 `Pt-L3-PlatinumChloride.xdi` 直接命中 `ptcl4_square` 并完成 S7 拟合？
  - [ ] 新增 #67–#68 后运行 `python run_ci.py --stress` 显示 **68 项契约回归 + 4 项压力测试 100% 全绿**？
- [ ] **文档同步**：`README.md` 与 `.trae/skills/xafs-pipeline-report/references/pipeline-reference.md` 更新至 **68 项契约回归**与 B13–B16 接口说明。

---

## 第四阶段终审验收 Checklist（Phase I–IV 全量判分，8 条硬指标）

1. **全量 CI 门禁**：`python run_ci.py --stress` 中 `py_compile` 退出码 0、`regression_test.py` **68 项全绿**、4 项压力测试退出码 0。
2. **默认零退化**：不传新开关跑默认 `LEPIDO` 三扫描，`r_factor ≈ 0.0273`、第一壳 `N ≈ 6.35`、`R ≈ 2.002 Å`、7 列拟合表逐字节不变。
3. **B13 催化位点建模**：`build_sac_motif` / `build_dac_motif`、CLI `--catalytic-model`、DAC 金属键 Hamilton F-test、异核共占位 `path_links` 双分支精修全部通过 #61–#62。
4. **B14 约束与审稿诊断**：Mini-GDS（`n_sum` / `n_ratio` / `deltar_offset` / `bounds`）、`--phase-correct` 真键长诊断、`--chi2-contour` 置信椭圆、`fit_multi_datasets(debye_models=...)` 全部通过 #63–#65。
5. **B15 原位序列分析**：`run_operando_series`（PCA + SPOIL 筛标样 + LCF + MCR-ALS + 描述符轨迹）与 `--operando` 报告产物通过 #66。
6. **B16 存量生态与多元素库**：`read_prj` 直接读通真实 `Fe_standards.prj` 且不破坏 #50；Co/Ni/Cu/Pt 催化 Motif 全部入库，`Pt-L3-PlatinumChloride.xdi` 命中 `ptcl4_square`，通过 #67–#68。
7. **CLI 与配置镜像**：`run_xafs.py --help`、`xafs_report.py --help` 与 `xafs_config.DEFAULTS` 完整包含 `--catalytic-model`、`--hetero-cocupancy`、`--constraints`、`--phase-correct`、`--chi2-contour`、`--operando`。
8. **物理诚实铁律**：`absorber.py` 五元素标定 $S_0^2$ 值未改动；一切新增约束与扩展均真实减少/核算自由度 $N_{\text{var}}$ 与 $N_{\text{idp}}$，绝不伪造数据或放宽 `cn_reliable` 门禁。

---

## 第五阶段（Phase V：B17–B20）最高优先级执行计划：催化剂精细结构深度解析与实战完备化（回归 68 → 74 项）

> **战略定位**：
> - 第四阶段（Phase IV，68 项回归）已打通 SAC/DAC、共占位固溶体、Mini-GDS、相移/等高线诊断与原位 XANES 时间序列。
> - 第五阶段（Phase V：B17–B20）聚焦真实催化剂结构解析中最核心的 **6 个深水区痛点**：
>   1. **负载型纳米颗粒尺寸效应**：从第一/二壳金属配位数 $(N_1, N_2)$ 定量反演团簇粒径 $d\ (\text{nm})$、原子总数、分散度 $D$ 及筏状平铺形貌判据；
>   2. **双金属纳米合金（如 PtCo / PtRu / CuAg）**：动态构建无需 CIF 的双金属合金 Motif（`alloy:...`），并定量计算 **Cowley 短程有序参数 $\alpha_{\text{Cowley}}$ 与合金度 $J$**（区分均匀固溶体、有序金属间化合物与核壳相分离）；
>   3. **轻元素（C/N/O）近简并壳层 Rayleigh 分辨极限护栏**：自动拦截 $\Delta R < \frac{\pi}{2\Delta k}$ 且 $|\Delta Z| \le 2$ 的过拟合拆壳并给出物理告警；
>   4. **多元素运行时 `--s02` 覆盖与真实 Au/Mo/Cr/Zn 金属箔解锁**：在不改动 `absorber.py` 静态表（严守 #60 契约）的前提下，让全部 17 个元素支持显式 `--s02` 跑通 S6.5/S7；
>   5. **原位（*Operando*）序列 EXAFS 壳层参数联合追踪**：全序列共享 $\Delta E_0$ 消除键长漂移伪影，批量导出 $N(t), R(t), \sigma^2(t)$ 演化曲线；
>   6. **电催化核心氧化物/羟基氧化物库完备化**：补齐 $\text{CoO}$、$\text{Co}_3\text{O}_4$（尖晶石双位点 `mixture`）、$\beta\text{-Ni(OH)}_2$、$\text{CuO}$、$\text{PtO}_2$，并在仓库真实实验标样上全盲测通过。

### 全局执行顺序总览（Phase V）

| 批次 | 编号 | 催化分析核心模块 | 修改文件 | 新增回归测试（`regression_test.py`） |
|---|---|---|---|---|
| **第十七批** | **B17-1** | 负载型金属纳米颗粒尺寸效应反演 `estimate_nanoparticle_size`（粒径 $d_{\text{nm}}$、分散度 $D$、球形 vs 筏状平铺判据）与 S7 自动挂载 | `artemis.py`, `xafs_fit.py` | `test_nanoparticle_size_and_dispersion_inversion` (#69) |
| **第十七批** | **B17-2** | 双金属纳米合金动态建模 `build_alloy_motif`（`alloy:...`）与 Cowley 短程有序参数 `compute_cowley_sro` | `structure_lib.py`, `artemis.py`, `xafs_fit.py` | `test_bimetallic_alloy_motif_and_cowley_sro_pipeline` (#70) |
| **第十八批** | **B18-1** | 轻元素（C/N/O）近简并壳层与 Rayleigh 径向分辨极限门禁 `check_shell_resolution_and_z_degeneracy` | `artemis.py`, `xafs_fit.py` | `test_rayleigh_resolution_and_light_element_z_degeneracy_guard` (#71) |
| **第十八批** | **B18-2** | 运行时 `--s02` 覆盖解锁未标定催化元素 + `au_fcc`/`mo_bcc`/`cr_bcc`/`zn_hcp` 入库与真实 Au/Mo 箔 S7 验证（严守 #60） | `structure_id.py`, `xafs_fit.py`, `xafs_config.py`, `run_xafs.py`, `xafs_report.py`, `_structure_lib/build_seed.py`, `calibrate_s02.py` | `test_s02_runtime_override_and_real_au_mo_foil_calibration` (#72) |
| **第十九批** | **B19-1** | 原位（*Operando*）序列 EXAFS 联合精修与 $N(t)/R(t)/\sigma^2(t)$ 轨迹导出 `run_operando_exafs_series` + `write_operando_exafs_diagnostics` | `xafs_fit.py`, `xafs_diagnostics.py` | `test_operando_exafs_series_joint_e0_bond_and_cn_tracking` (#73) |
| **第二十批** | **B20-1** | 电催化氧化物/羟基氧化物库完备化（`coo_rocksalt`, `co3o4_tet/oct/mixed`, `ni_oh2_brucite`, `cuo_tenorite`, `pto2_beta`）与真实标准谱盲测 | `_structure_lib/build_seed.py` | `test_electrocatalysis_oxide_hydroxide_library_real_standards` (#74) |

---

## 第十七批（B17）：负载型金属纳米颗粒尺寸反演与双金属纳米合金短程有序（CSRO）

### B17-1 负载型金属纳米颗粒尺寸效应与形貌/分散度反演（`estimate_nanoparticle_size`）

#### 1. 催化痛点与数学契约
- **痛点**：负载型金属催化剂（Pt/C、Cu/ZnO、Ni/Al₂O₃ 等）因表面原子配位不饱和，第一壳金属配位数 $N_1 < N_{1,\text{bulk}}$（fcc/hcp 为 12，bcc 为 8）。S7 拟合目前仅给出 $N_1$，缺少从 $(N_1, R_1, N_2)$ 反推纳米颗粒等效直径 $d\ (\text{nm})$、表面分散度 $D$ 以及区分“三维球形颗粒”与“载体强相互作用下的二维筏状/扁平团簇（raft）”的定量模块。
- **数学公式与接口契约**：
  1. 在 `artemis.py` 新增函数：
     ```python
     def estimate_nanoparticle_size(
         n1: float,
         r1: float,
         n1_err: float | None = None,
         n2: float | None = None,
         lattice: str = "fcc",
         shape: str = "spherical",
     ) -> dict
     ```
     - 支持 `lattice in ("fcc", "hcp", "bcc")`（其余抛 `XAFSConfigError`）：
       - `fcc` / `hcp`：$N_{1,\text{bulk}} = 12.0, N_{2,\text{bulk}} = 6.0, r_2 = \sqrt{2} r_1$，原子体积填充因子 $\eta = \frac{\pi}{3\sqrt{2}} \approx 0.74048$；
       - `bcc`：$N_{1,\text{bulk}} = 8.0, N_{2,\text{bulk}} = 6.0, r_2 = \frac{2}{\sqrt{3}} r_1$，填充因子 $\eta = \frac{\sqrt{3}\pi}{8} \approx 0.68017$；
     - **球形团簇解析截断方程**（Calvin / Beale 对关联球截断公式）：
       对于直径 $d$（单位 Å）的球形纳米颗粒，距离 $r$ 壳层的平均配位数为：
       $$N(r, d) = N_{\text{bulk}} \left[ 1 - \frac{3}{2} u + \frac{1}{2} u^3 \right], \quad u = \frac{r}{d} \in (0, 1]$$
       - 当 $N_1 \ge 0.99 \cdot N_{1,\text{bulk}}$ 时：返回体相极限 `regime = "bulk_limit"`，`diameter_A = float("inf")`，`diameter_nm = float("inf")`，`diameter_nm_err = None`，`n_atoms = float("inf")`，`dispersion = 0.0`；
       - 当 $0 < N_1 < 0.99 \cdot N_{1,\text{bulk}}$ 时：令 $f = \text{clip}(N_1 / N_{1,\text{bulk}}, 0.05, 0.99)$，在 $u \in (0, 1]$ 上用 `scipy.optimize.brentq` 精确求解 $1 - \frac{3}{2}u + \frac{1}{2}u^3 - f = 0$，得到球形直径 $d_{\text{Å}} = r_1 / u$（若 `shape == "hemispherical"`，基底直径乘以 $2^{1/3}$），并换算为纳米：
         $$\text{diameter\_nm} = \frac{d_{\text{Å}}}{10.0}$$
       - **形式误差传播**：由隐函数求导 $\frac{\partial f}{\partial u} = -\frac{3}{2}(1 - u^2)$，$\frac{\partial d_{\text{Å}}}{\partial N_1} = -\frac{r_1}{u^2 N_{1,\text{bulk}} (\partial f / \partial u)} = \frac{2 r_1}{3 u^2 (1 - u^2) N_{1,\text{bulk}}}$，当 `n1_err` 非空时计算：
         $$\text{diameter\_nm\_err} = \frac{1}{10.0} \left|\frac{\partial d_{\text{Å}}}{\partial N_1}\right| \cdot \sigma(N_1)$$
       - **原子总数与表面分散度**：
         - 单原子体积 $V_{\text{atom}} = \frac{4\pi}{3} \left(\frac{r_1}{2}\right)^3 / \eta$，颗粒内原子总数 $N_{\text{atoms}} = \max\left(2.0, \eta \left(\frac{d_{\text{Å}}}{r_1}\right)^3\right)$；
         - 表面分散度估计 $D = \text{clip}\left(\frac{4.0}{N_{\text{atoms}}^{1/3}}, 0.0, 1.0\right)$；
       - **第二壳形貌自洽检验（球形 vs 二维筏状）**：
         - 计算直径 $d_{\text{Å}}$ 下的球形理论第二壳配位数：若 $d_{\text{Å}} \le r_2$ 则 $N_{2,\text{sphere}} = 0.0$，否则 $N_{2,\text{sphere}} = N_{2,\text{bulk}} \left[1 - \frac{3}{2}(r_2/d_{\text{Å}}) + \frac{1}{2}(r_2/d_{\text{Å}})^3\right]$；
         - 若传入了 `n2 is not None`：当 $N_{2,\text{sphere}} > 0.5$ 且 $n_2 < 0.65 \cdot N_{2,\text{sphere}}$ 时，置 `"shape_hint": "raft_or_flattened"`（第二壳显著低于球形预期，提示二维平铺/筏状团簇），否则 `"spherical_consistent"`；未传 `n2` 时置 `"single_shell_only"`。
     - 返回字典必须包含：`{"diameter_nm": float, "diameter_A": float, "diameter_nm_err": float | None, "n_atoms": float, "dispersion": float, "lattice": str, "n1_bulk": float, "n2_spherical_expected": float, "shape_hint": str, "regime": str}`。
  2. 在 `xafs_fit.fit_identified` 中：当拟合出的壳层全部为同核金属散射体（`all(sh["scatterer"] == absorber.symbol for sh in plan["shells"])`）且第一壳拟合 $N_1 < 0.96 \cdot N_{1,\text{bulk}}$ 时，自动根据 Motif 晶格（`fe_bcc`/`mo_bcc`/`cr_bcc` 取 `"bcc"`，`co_hcp`/`zn_hcp` 取 `"hcp"`，其余取 `"fcc"`）调用 `artemis.estimate_nanoparticle_size` 并写入返回字典顶层 `out["nanoparticle_size"]`（体相金属箔 $N_1 \ge 0.96 \cdot N_{1,\text{bulk}}$ 时 `out["nanoparticle_size"] = None`，保持默认输出干净）。

#### 2. 验收标准（新增回归 #69：`test_nanoparticle_size_and_dispersion_inversion`）
1. 验证 `artemis.estimate_nanoparticle_size(n1=8.0, r1=2.774, n1_err=0.5, lattice="fcc")`：精确反解得 $u = r_1 / d_{\text{Å}}$ 满足 $N_{\text{bulk}}(1 - 1.5u + 0.5u^3) = 8.0$（误差 $< 10^{-9}$），`diameter_nm` 在 $[1.0, 1.5]\ \text{nm}$ 之间，`diameter_nm_err > 0`，`0.3 < dispersion <= 1.0`；
2. 验证当 $N_1$ 从 $6.0 \to 9.0 \to 11.0$ 增大时，`diameter_nm` 严格单调递增而 `dispersion` 严格单调递减；当 $N_1 = 12.0$ 时返回 `regime == "bulk_limit"` 且 `dispersion == 0.0`；
3. 验证传入显著偏低的第二壳 `n2=1.0`（球形预期 $>3.0$）时，`shape_hint == "raft_or_flattened"`；传入与球形一致的 `n2` 时为 `"spherical_consistent"`。

---

### B17-2 双金属纳米合金动态建模（`build_alloy_motif`）与 Cowley 短程有序参数（`compute_cowley_sro`）

#### 1. 催化痛点与数学契约
- **痛点**：双金属纳米合金催化剂（PtCo、PtRu、PdAu、CuNi）需要同时拟合同核金属路径 $A\text{–}A$ 与异核金属路径 $A\text{–}B$，并用 **Cowley 短程有序参数 $\alpha_{\text{Cowley}}$** 判断合金内部是均匀随机固溶体、有序金属间化合物还是核壳相分离。此前用户必须手工造共占位 CIF 才能跑双金属拟合。
- **数学公式与接口契约**：
  1. 在 `structure_lib.py` 新增动态双金属合金团簇构造函数：
     ```python
     def build_alloy_motif(
         absorber: str,
         partner: str,
         x_partner: float = 0.5,
         r_aa: float | None = None,
         r_ab: float | None = None,
         n_total: int = 12,
         lattice: str = "fcc",
     ) -> Motif
     ```
     - 参数校验：`absorber != partner`，$0.0 < x_{\text{partner}} < 1.0$，`n_total >= 2`（否则抛 `XAFSConfigError`）；
     - 缺省键长：若 `r_aa is None` 取 $2 \cdot R_{\text{cov}}(\text{absorber})$，若 `r_ab is None` 取 $R_{\text{cov}}(\text{absorber}) + R_{\text{cov}}(\text{partner})$；若 $|r_{aa} - r_{ab}| < 0.01\ \text{Å}$，将两距离设成精确相等以代表同壳层；
     - 在第一配位壳的 `n_total` 个立方八面体（fcc 12 方向 `(±1, ±1, 0) / √2` 等）方向上：前 $n_{ab} = \max(1, \min(n_{\text{total}} - 1, \text{ int}(\text{round}(n_{\text{total}} \cdot x_{\text{partner}}))))$ 个方向放距离为 `r_ab` 的 `partner` 原子，其余 $n_{aa} = n_{\text{total}} - n_{ab}$ 个方向放距离为 `r_aa` 的 `absorber` 原子；
     - 构造 `Motif(id=f"alloy_{absorber.lower()}_{partner.lower()}_{hash8}", absorber=absorber, kind="cluster", nuclearity="alloy", source={"type": "dynamic_alloy", "partner": partner, "x_partner": float(x_partner), "r_aa": float(r_aa), "r_ab": float(r_ab), "n_total": int(n_total)}, sites=..., shell_expect=..., cache_dir=str(DYNAMIC_CACHE_ROOT / f"alloy_{hash12}"))`；
     - 扩展 `structure_lib.parse_catalytic_cluster_spec(absorber, spec)`：当字符串以 `"alloy:"` 开头（如 `"alloy:Co@2.65,Pt@2.77,x=0.25,n=12"` 或 `"alloy:partner=Co,r_ab=2.65,r_aa=2.77,x=0.25"`）或字典 `spec.get("type") == "alloy"` 时，调用 `build_alloy_motif` 返回该 Motif。
  2. 在 `artemis.py` 新增 Cowley 短程有序参数计算函数：
     ```python
     def compute_cowley_sro(
         n_aa: float,
         n_ab: float,
         x_partner: float,
         n_aa_err: float | None = None,
         n_ab_err: float | None = None,
     ) -> dict
     ```
     - 校验 $0.0 < x_{\text{partner}} < 1.0$ 且 $N_{AA} \ge 0, N_{AB} \ge 0, N_{AA} + N_{AB} > 0$（否则抛 `XAFSConfigError`）；
     - 计算总金属配位数 $N_{\text{total}} = N_{AA} + N_{AB}$，异核配位分数 $p_{AB} = \frac{N_{AB}}{N_{\text{total}}}$，**Cowley 短程有序参数**：
       $$\alpha_{\text{Cowley}} = 1.0 - \frac{p_{AB}}{x_{\text{partner}}}$$
       以及合金化程度百分比 $J_{\text{alloy}} = \frac{p_{AB}}{x_{\text{partner}}} \times 100.0$；
     - 误差传播：当 `n_aa_err` 与 `n_ab_err` 非空时，按偏导数
       $$\sigma(\alpha) = \frac{1}{x_{\text{partner}} N_{\text{total}}^2} \sqrt{(N_{AB} \cdot \sigma_{N_{AA}})^2 + (N_{AA} \cdot \sigma_{N_{AB}})^2}$$
       计算 `alpha_err`；
     - 局域化学有序趋势 `tendency`：
       - $\alpha_{\text{Cowley}} < -0.15$：`"heteroatomic_ordering"`；
       - $|\alpha_{\text{Cowley}}| \le 0.15$：`"random_solid_solution"`；
       - $\alpha_{\text{Cowley}} > 0.15$：`"homoatomic_clustering_or_segregation"`；
     - 返回 `{"n_aa": float, "n_ab": float, "n_total": float, "p_ab": float, "x_partner": float, "alpha_cowley": float, "alpha_err": float | None, "alloying_extent_pct": float, "tendency": str}`。
  3. 在 `xafs_fit.fit_identified` 中：当 `motif.nuclearity == "alloy"` 时，从拟合出的路径中提取同核路径 $N_{AA}$ 与异核路径 $N_{AB}$，自动调用 `artemis.compute_cowley_sro` 写入返回字典顶层 `out["alloy_sro"]`（同时若 $N_{\text{total}} < 11.5$ 也自动计算 `out["nanoparticle_size"]`）。

#### 2. 验收标准（新增回归 #70：`test_bimetallic_alloy_motif_and_cowley_sro_pipeline`）
1. 验证 `compute_cowley_sro` 对三种典型情形的精确数学判别：
   - 随机固溶体（$x_B = 0.25, N_{AA} = 9.0, N_{AB} = 3.0$）：`alpha_cowley == 0.0`、`alloying_extent_pct == 100.0`、`tendency == "random_solid_solution"`；
   - 核壳偏析/同核团聚（$x_B = 0.50, N_{AA} = 8.0, N_{AB} = 2.0$）：`alpha_cowley == 0.60`、`tendency == "homoatomic_clustering_or_segregation"`；
   - 异核有序（$x_B = 0.25, N_{AA} = 6.0, N_{AB} = 6.0$）：`alpha_cowley == -1.0`、`tendency == "heteroatomic_ordering"`；
2. 用 `parse_catalytic_cluster_spec("Pt", "alloy:Co@2.64,Pt@2.76,x=0.25,n=12")` 生成动态合金 Motif，取其两条 FEFF 路径按真值 $N_{\text{Pt-Pt}} = 6.0, N_{\text{Pt-Co}} = 2.0$（对应总配位数 $N_{\text{total}} = 8.0$ 的约 $1.2\ \text{nm}$ 随机固溶体合金颗粒，$x_{\text{Co}} = 0.25, \alpha_{\text{true}} = 0.0$）合成 $\chi(k)$，调用 `xafs_fit.fit_identified(..., catalytic_model="alloy:Co@2.64,Pt@2.76,x=0.25,n=12", shared_sigma2="all")`：
   - 验证回收 $|N_{\text{Pt-Pt}} - 6.0| < 0.5$、$|N_{\text{Pt-Co}} - 2.0| < 0.5$；
   - 验证 `out["alloy_sro"]["tendency"] == "random_solid_solution"` 且 $|\alpha_{\text{Cowley}}| < 0.10$；
   - 验证 `out["nanoparticle_size"]["diameter_nm"]` 在 $[1.0, 1.5]\ \text{nm}$ 之间。

---

## 第十八批（B18）：轻元素近简并 Rayleigh 分辨极限门禁与 17 元素 `--s02` 运行时解锁

### B18-1 轻元素（C/N/O）近简并壳层与 Rayleigh 分辨率门禁（`check_shell_resolution_and_z_degeneracy`）

#### 1. 催化痛点与数学契约
- **痛点**：EXAFS 在有限 $k$ 窗口 $[k_{\min}, k_{\max}]$ 内对同元素/近邻元素双壳层的最小可分辨距离由 Rayleigh 判据给出：
  $$\Delta R_{\text{Rayleigh}} = \frac{\pi}{2 (k_{\max} - k_{\min})}$$
  当 $|\Delta Z| \le 2$（如 C/N/O，背散射振幅包络与相移函数几乎重合）且两壳层距离差 $|\Delta R| < \Delta R_{\text{Rayleigh}}$ 时，若不加约束地同时拟合两条路径，极易发生反相干涉抵消，产生毫无物理意义的配位数。
- **接口与数学契约**：
  1. 在 `artemis.py` 新增函数：
     ```python
     def check_shell_resolution_and_z_degeneracy(
         shells: list[dict],
         kmin: float = 3.0,
         kmax: float = 12.0,
     ) -> dict
     ```
     - 输入 `shells`：每个元素为含 `{"label": str, "scatterer": str, "r": float (或 "reff": float)}` 的字典（兼容 `select_shell_paths` 的 `plan["shells"]` 与 `FitResult.paths` 转换）；
     - 计算 $\Delta k = \max(float(k_{\max}) - float(k_{\min}), 10^{-6})$ 与 $\Delta R_{\text{rayleigh}} = \frac{\pi}{2 \Delta k}$；
     - 对所有路径对 $(i < j)$：
       - 用 `absorber.ATOMIC_NUMBERS`（或内置元素周期表映射）查散射元素原子序数 $Z_i, Z_j$；
       - 计算 $\Delta Z = |Z_i - Z_j|$ 与 $\Delta R = |R_i - R_j|$；
       - 若 $\Delta Z \le 2$ 且 $\Delta R < \Delta R_{\text{rayleigh}}$，加入 `pair_warnings`：
         ```python
         {
             "path_a": label_i,
             "path_b": label_j,
             "scatterer_a": el_i,
             "scatterer_b": el_j,
             "delta_z": int(delta_z),
             "delta_r": float(delta_r),
             "rayleigh_limit_A": float(rayleigh_limit),
             "kind": ("light_element_z_degeneracy" if max(z_i, z_j) <= 10 else "sub_rayleigh_split"),
         }
         ```
     - 返回 `{"rayleigh_limit_A": float, "delta_k": float, "resolvable": len(pair_warnings) == 0, "pair_warnings": pair_warnings, "recommendation": str}`。
  2. 在 `xafs_fit.fit_identified` 中：根据拟合后各路径的实际键长 `p.r` 与散射元素调用 `artemis.check_shell_resolution_and_z_degeneracy`，将结果存入 `out["resolution_check"]`，并在 `out["flags"]["sub_rayleigh_z_degeneracy"] = not res_check["resolvable"]`；若触发且 `shared_sigma2 is None and constraints is None`，在 `out["warning"]` 中追加 `"壳层间距低于 Rayleigh 极限且散射元素近简并，建议绑定 shared_sigma2/constraints 或用 BVS 消歧"`。

#### 2. 验收标准（新增回归 #71：`test_rayleigh_resolution_and_light_element_z_degeneracy_guard`）
1. 验证 $k \in [3.0, 12.0]\ \text{Å}^{-1}$ 时 `rayleigh_limit_A` 精确等于 $\pi / 18 \approx 0.17453\ \text{Å}$；
2. 对 `Fe-N@1.98` 与 `Fe-O@2.06`（$\Delta Z = 1, \Delta R = 0.08 < 0.1745\ \text{Å}$），验证 `resolvable is False`、`pair_warnings[0]["kind"] == "light_element_z_degeneracy"`；
3. 对 `Fe-N@1.98` 与 `Fe-Co@2.50`（$\Delta Z = 20, \Delta R = 0.52\ \text{Å}$）或 `Fe-O@1.98` 与 `Fe-O@2.25`（$\Delta R = 0.27 > 0.1745\ \text{Å}$），验证 `resolvable is True`、`pair_warnings == []`；
4. 验证 `fit_identified` 在拟合含近简并 `Fe-N@2.00` + `Fe-O@2.10` 的 SAC 模型时，`out["resolution_check"]["resolvable"] is False` 且 `out["flags"]["sub_rayleigh_z_degeneracy"] is True`。

---

### B18-2 17 元素全表运行时 `--s02` 覆盖与真实 Au / Mo 金属箔全链解锁（严守 #60）

#### 1. 核心防呆陷阱与接口契约
- **防呆红线（必读）**：
  - `regression_test.py:2726`（回归 #60）显式断言 `absorber.py` 中 `("Cr", "Mn", "Zn", "Mo", "Ru", "Rh", "Pd", "Ag", "Ir", "Au", "Ce", "U")` 共 12 个元素的静态配置 `get_absorber(sym).s02 is None` 且 `calibrated is False`！
  - **绝对禁止**直接把 `absorber.py` 静态字典里这 12 个元素的 `s02` 改成浮点数！
- **解锁机制与接口契约**：
  1. 在 `_structure_lib/build_seed.py` 新增并落盘 4 个高对称性催化金属单质种子 Motif（供 `calibrate_s02.py` 与 S6.5/S7 使用）：
     - **Au**：`au_fcc`（`absorber="Au"`, `kind="crystal"`, `spacegroup="Fm-3m"`, `cell={"a": 4.0782}`，`Au(0,0,0)`，第一壳 Au–Au $12 \times 2.8837\ \text{Å}$）；
     - **Mo**：`mo_bcc`（`absorber="Mo"`, `kind="crystal"`, `spacegroup="Im-3m"`, `cell={"a": 3.1470}`，`Mo(0,0,0)`，第一壳 Mo–Mo $8 \times 2.7254\ \text{Å}$，第二壳 $6 \times 3.1470\ \text{Å}$）；
     - **Cr**：`cr_bcc`（`absorber="Cr"`, `kind="crystal"`, `spacegroup="Im-3m"`, `cell={"a": 2.8846}`，`Cr(0,0,0)`，第一壳 Cr–Cr $8 \times 2.4981\ \text{Å}$）；
     - **Zn**：`zn_hcp`（`absorber="Zn"`, `kind="crystal"`, `spacegroup="P 63/m m c"`, `cell={"a": 2.6649, "b": 2.6649, "c": 4.9468, "alpha": 90, "beta": 90, "gamma": 120}`，`Zn(1/3, 2/3, 0.25)`）；
     - 在 `calibrate_s02.DEFAULT_MOTIF` 同步加入 `"Au": "au_fcc", "Mo": "mo_bcc", "Cr": "cr_bcc", "Zn": "zn_hcp"`。
  2. 在 `structure_id.identify_result`、`xafs_fit.fit_identified`、`xafs_fit.run` 新增可选参数 `s02: float | None = None`；在 `xafs_config.DEFAULTS` 新增 `"s02": None`；在 `run_xafs.py` 与 `xafs_report.py` 新增 CLI 参数 `--s02`（`type=float, default=None`）：
     - 当 `s02 is not None` 时：校验 $0.4 \le \text{float}(s02) \le 1.3$（越界抛 `XAFSConfigError`），用 `dataclasses.replace(get_absorber(absorber_symbol), s02=float(s02))` 构造运行时已标定 `AbsorberConfig` 实例（不修改全局 `absorber.ABSORBERS` 字典！），使 `identify_result` 与 `fit_identified` 直接放行并使用该 `s02` 完成识别与拟合；
     - 当 `s02 is None` 时：保持原有逻辑不变（未标定元素依然抛 `XAFSDataError`）。

#### 2. 验收标准（新增回归 #72：`test_s02_runtime_override_and_real_au_mo_foil_calibration`）
1. 验证不传 `s02` 时，对 `"Au"` 调用 `structure_id.identify_result` 仍抛 `XAFSDataError`（且 `get_absorber("Au").s02 is None`）；
2. 对仓库真实实验谱 `data/bmm_standards/Data/Au/Au-L3-GoldFoil.xdi` 传入 `s02=0.82` 运行 `xafs_fit.run(res_au, identify=True, fit=True, absorber_symbol="Au", s02=0.82)`：
   - 验证 `out["structure"]["ranked"][0]["id"] == "au_fcc"` 且 `status == "confident"`；
   - 验证 S7 第一壳拟合回收 $N \in [11.0, 13.0]$、$R \in [2.83, 2.90]\ \text{Å}$、`r_factor < 0.02`；
3. 对仓库真实实验谱 `data/bmm_standards/Data/Mo/Mo-K-MolybdenumFoil.xdi` 传入 `s02=0.88` 运行 `xafs_fit.run(res_mo, identify=True, fit=True, absorber_symbol="Mo", s02=0.88)`：
   - 验证 `out["structure"]["ranked"][0]["id"] == "mo_bcc"`；
   - 验证 S7 双壳拟合回收第一壳 $N_1 \in [6.8, 9.0]$、$R_1 \in [2.69, 2.75]\ \text{Å}$ 及第二壳 $N_2 \in [5.0, 7.2]$、$R_2 \in [3.11, 3.18]\ \text{Å}$、`r_factor < 0.02`。

---

## 第十九批（B19）：原位（*Operando*）序列 EXAFS 壳层参数联合追踪

### B19-1 原位时间序列 EXAFS 联合精修与键长/配位数轨迹导出（`run_operando_exafs_series`）

#### 1. 催化痛点与数学契约
- **痛点**：原位催化过程中，研究者不仅要看 XANES 的价态/组分变化（B15-1），更要定量追踪**活性位点键长 $R(t)$ 的伸缩（如电催化 OER 中高价金属–氧键缩短 $0.05\ \text{Å}$）与配位数 $N(t)$ 的演变**。若逐条谱独立跑 S7，$\Delta E_0$ 的随机涨落会通过 $\Delta E_0\text{–}\Delta R$ 强相关性转化为虚假的键长抖动。
- **数学与接口契约**：
  1. 在 `xafs_fit.py` 新增函数：
     ```python
     def run_operando_exafs_series(
         results: list,
         absorber_symbol: str,
         motif: str | Motif | None = None,
         catalytic_model: str | dict | Motif | None = None,
         rmin: float | None = None,
         rmax: float | None = None,
         kmin: float = 3.0,
         kmax: float = 12.0,
         kweight: int = 2,
         s02: float | None = None,
         share_sigma2_across_series: bool = False,
         names: list[str] | None = None,
     ) -> dict
     ```
     - 输入 `results`：至少 3 个 `PipelineResult`（或含 `chi_k, ft` 的对象，`<3` 抛 `XAFSConfigError`）；
     - 解析目标 `Motif`：优先用 `catalytic_model`（调 `parse_catalytic_cluster_spec`），其次用显式 `motif`（id 或 `Motif` 对象），否则对 `results[0]` 调用 `structure_id.identify_result` 取 `ranked[0]["id"]`；
     - 调用 `select_shell_paths(target_motif, edge=absorber.edge, rmin=rmin, rmax=rmax, kmin=kmin, kmax=kmax)` 获得统一的壳层路径文件 `plan["files"]`、标签 `plan["labels"]` 与初始简并度 `plan["degens"]`；
     - 构造 `datasets` 列表并调用 `artemis.fit_multi_datasets(datasets, shared_bonds=shared_bonds, share_e0=True)`：
       - **全序列共享单一全局 $\Delta E_0$**（`share_e0=True`，彻底消除时间点间 $\Delta E_0$ 漂移造成的键长伪波动）；
       - 每个时间点独立放开各壳层的 $N_{d,p}$ 与 $\Delta R_{d,p}$；
       - 当 `share_sigma2_across_series=True` 时，为每个壳层 $p$ 构造一条 `shared_bonds` 项：`{"paths": [(d, p) for d in range(len(results))], "share_deltar": False, "share_sigma2": True, "share_n": False}`，使同一壳层在全序列共享同一个 $\sigma^2_p$，显著压低单点 $N(t)$ 误差；
     - 组装返回字典：
       ```python
       {
           "names": names_list,
           "motif": target_motif.id,
           "shells": [sh["label"] for sh in plan["shells"]],
           "global_r_factor": float(multi_out["global_r_factor"]),
           "shared_delta_e0": float(multi_out["datasets"][0].delta_e0),
           "shared_delta_e0_err": multi_out["datasets"][0].delta_e0_err,
           "r_factors": [float(f.r_factor) for f in multi_out["datasets"]],
           "trajectories": {
               label: {
                   "reff": float(plan["shells"][pi]["reff"]),
                   "N": [float(f.paths[pi].degeneracy) for f in multi_out["datasets"]],
                   "N_err": [f.paths[pi].n_err for f in multi_out["datasets"]],
                   "R": [float(f.paths[pi].r) for f in multi_out["datasets"]],
                   "R_err": [f.paths[pi].delta_r_err for f in multi_out["datasets"]],
                   "sigma2": [float(f.paths[pi].sigma2) for f in multi_out["datasets"]],
                   "sigma2_err": [f.paths[pi].sigma2_err for f in multi_out["datasets"]],
               }
               for pi, label in enumerate(plan["labels"])
           },
           "multi_fit": multi_out,
       }
       ```
  2. 在 `xafs_diagnostics.py` 新增诊断落盘函数：
     ```python
     def write_operando_exafs_diagnostics(
         exafs_series_out: dict,
         outdir,
         stem: str = "operando_exafs",
     ) -> dict
     ```
     - 导出 `<stem>_exafs_evolution.csv`（表头：`index,name,r_factor,delta_e0` + 每个壳层的 `<label>_N,<label>_N_err,<label>_R,<label>_R_err,<label>_sigma2`）；
     - 导出双面板图 `<stem>_exafs_evolution.png`（上图：各壳层 $N(t) \pm \sigma_N$ 随时间点演化；下图：各壳层真实键长 $R(t) \pm \sigma_R$ 随时间点演化）；
     - 返回 `{"csv": str(csv_path), "png": str(png_path)}`。

#### 2. 验收标准（新增回归 #73：`test_operando_exafs_series_joint_e0_bond_and_cn_tracking`）
- 用 `Fe–O` 第一壳 FEFF 路径（$R_{\text{eff}} \approx 1.98\ \text{Å}$）合成 4 个原位电位点的序列谱（模拟 OER 氧化过程：真值键长从 $2.02 \to 1.99 \to 1.96 \to 1.93\ \text{Å}$ 线性缩短，配位数从 $4.0 \to 4.6 \to 5.3 \to 6.0$ 单调增加，全局 $\Delta E_0 = +2.5\ \text{eV}, \sigma^2 = 0.004\ \text{Å}^2$）：
  1. 调用 `xafs_fit.run_operando_exafs_series(results, absorber_symbol="Fe", catalytic_model="sac:O4@1.98,C=none", share_sigma2_across_series=True)`；
  2. 验证回收的 4 个时间点键长 `traj["R"]` 与真值 `[2.02, 1.99, 1.96, 1.93]` 的最大绝对误差 `< 0.005 Å`（严格单调递减），配位数 `traj["N"]` 与真值最大绝对误差 `< 0.15`（严格单调递增），`global_r_factor < 1e-3`；
  3. 验证 `write_operando_exafs_diagnostics` 生成非空的 `_exafs_evolution.csv`（5 行：表头 + 4 个时间点）与 `_exafs_evolution.png`。

---

## 第二十批（B20）：电催化核心氧化物/羟基氧化物库完备化与真实标准谱盲测

### B20-1 Co / Ni / Cu / Pt 电催化核心氧化物与尖晶石双位点入库

#### 1. 扩充清单与晶体学契约（已用仓库真实 BMM 标准谱预验证）
- 在 `_structure_lib/build_seed.py` 新增并运行落盘以下 7 个电催化核心 Motif：
  1. **Co 岩盐氧化物 `coo_rocksalt`**：
     - `absorber="Co"`, `kind="crystal"`, `spacegroup="Fm-3m"`, `cell={"a": 4.2612}`，`asym_sites`: $\text{Co}(0,0,0), \text{O}(0.5,0.5,0.5)$（第一壳 $\text{Co–O}\times 6 @ 2.1306\ \text{Å}$、第二壳 $\text{Co–Co}\times 12 @ 3.0131\ \text{Å}$）。
  2. **Co 尖晶石 $\text{Co}_3\text{O}_4$ 双位点体系（`co3o4_tet` + `co3o4_oct` + `co3o4_mixed`）**：
     - 晶体骨架：`spacegroup="F d -3 m:2"`, `cell={"a": 8.084}`，`asym_sites`: $\text{Co}_{8a}(0.125, 0.125, 0.125)$, $\text{Co}_{16d}(0.5, 0.5, 0.5)$, $\text{O}_{32e}(0.263, 0.263, 0.263)$；
     - `co3o4_tet`：`absorber_fract=[0.125, 0.125, 0.125]`（四面体位，占 $1/3$，第一壳 $\text{Co–O}\times 4 @ 1.932\ \text{Å}$）；
     - `co3o4_oct`：`absorber_fract=[0.5, 0.5, 0.5]`（八面体位，占 $2/3$，第一壳 $\text{Co–O}\times 6 @ 1.922\ \text{Å}$、第二壳 $\text{Co–Co}\times 6 @ 2.858\ \text{Å}$）；
     - `co3o4_mixed`：`kind="mixture"`, `nuclearity="mixed-site"`, `components=[{"id": "co3o4_tet", "weight": 1.0/3.0}, {"id": "co3o4_oct", "weight": 2.0/3.0}]`。
  3. **Ni 氢氧化镍 `ni_oh2_brucite`（$\beta\text{-Ni(OH)}_2$）**：
     - `absorber="Ni"`, `kind="crystal"`, `spacegroup="P -3 m 1"`, `cell={"a": 3.126, "b": 3.126, "c": 4.605, "alpha": 90, "beta": 90, "gamma": 120}`，`asym_sites`: $\text{Ni}(0,0,0), \text{O}(1/3, 2/3, 0.22)$（第一壳 $\text{Ni–O}\times 6 @ 2.070\ \text{Å}$、第二壳 $\text{Ni–Ni}\times 6 @ 3.126\ \text{Å}$）。
  4. **Cu 单斜氧化铜 `cuo_tenorite`（$\text{CuO}$）**：
     - `absorber="Cu"`, `kind="crystal"`, `spacegroup="C 1 2/c 1"`, `cell={"a": 4.6837, "b": 3.4226, "c": 5.1288, "alpha": 90, "beta": 99.54, "gamma": 90}`，`asym_sites`: $\text{Cu}(0.25, 0.25, 0.0), \text{O}(0.0, 0.4184, 0.25)$，`absorber_fract=[0.25, 0.25, 0.0]`（第一壳 $\text{Cu–O}\times 4 @ 1.956\ \text{Å}$）。
  5. **Pt 氧化物 `pto2_beta`（$\alpha/\beta\text{-PtO}_2$ 层状六方参照相）**：
     - `absorber="Pt"`, `kind="crystal"`, `spacegroup="P -3 m 1"`, `cell={"a": 3.100, "b": 3.100, "c": 4.161, "alpha": 90, "beta": 90, "gamma": 120}`，`asym_sites`: $\text{Pt}(0,0,0), \text{O}(1/3, 2/3, 0.24)$（第一壳 $\text{Pt–O}\times 6 @ 2.047\ \text{Å}$、第二壳 $\text{Pt–Pt}\times 6 @ 3.100\ \text{Å}$）。

#### 2. 验收标准（新增回归 #74：`test_electrocatalysis_oxide_hydroxide_library_real_standards`）
1. 验证仓库 4 条真实电催化标准实验谱在 `structure_id.identify_result` / `xafs_fit.run` 下的盲测表现：
   - `data/bmm_standards/Data/Co/Co-K-CoO.xdi`：`ranked[0]["id"] == "coo_rocksalt"`、`status == "confident"`（`score < 0.06`），S7 拟合 `r_factor < 0.02`；
   - `data/bmm_standards/Data/Co/Co-K-Co3O4.xdi`：`ranked[0]["id"] == "co3o4_mixed"`、`status == "confident"`（`score < 0.05`，`uncovered_peaks == []`）；
   - `data/bmm_standards/Data/Ni/Ni-K-NickelHydroxide.xdi`：`ranked[0]["id"] == "ni_oh2_brucite"`、`status == "confident"`（`score < 0.14`），S7 拟合 `r_factor < 0.05` 且第一壳 $R \in [2.03, 2.11]\ \text{Å}$、第二壳 $R \in [3.08, 3.17]\ \text{Å}$；
   - `data/bmm_standards/Data/Cu/Cu-CuO.xdi`：`ranked[0]["id"] == "cuo_tenorite"`（`score < 0.10`）；
2. 运行 `stress_test_{co,ni,cu,nio}.py`，验证 4 项压力测试 100% 通过（金属箔与 NiO 的 `top1` 零退化）。

---

## 第五阶段终审验收 Checklist（Phase I–V 全量判分，8 条硬指标）

1. **全量 CI 门禁**：`python run_ci.py --stress` 中 `py_compile` 退出码 0、`regression_test.py` **74 项全绿**、4 项压力测试退出码 0。
2. **默认零退化**：不传新开关跑默认 `LEPIDO` 三扫描，`r_factor ≈ 0.0273`、第一壳 `N ≈ 6.35`、`R ≈ 2.002 Å`、7 列拟合表逐字节不变。
3. **B17 纳米粒径与双金属合金**：`estimate_nanoparticle_size`、动态合金 `build_alloy_motif`（`alloy:...`）与 `compute_cowley_sro` 全部通过 #69–#70。
4. **B18 分辨率护栏与多元素 `--s02` 解锁**：`check_shell_resolution_and_z_degeneracy` 准确拦截亚 Rayleigh 近简并拆壳；`--s02` 解锁真实 `Au-L3-GoldFoil.xdi` 与 `Mo-K-MolybdenumFoil.xdi` 全链 S6.5/S7 且不破坏 #60，通过 #71–#72。
5. **B19 原位 EXAFS 序列追踪**：`run_operando_exafs_series` 全序列共享 $\Delta E_0$ 准确回收 $R(t), N(t)$ 单调演变并导出 CSV/PNG，通过 #73。
6. **B20 电催化氧化物/羟基氧化物库**：`CoO`、`Co3O4`（双位点 `mixture`）、`Ni(OH)2`、`CuO`、`PtO2` 入库并在 4 条真实 BMM 氧化物标准谱上盲测命中，通过 #74。
7. **CLI 与配置同步**：`--s02` 同步接入 `run_xafs.py`、`xafs_report.py` 与 `xafs_config.DEFAULTS`。
8. **文档与清理**：`README.md` 与 `.trae/skills/xafs-pipeline-report/references/pipeline-reference.md` 同步更新至 **74 项契约回归**与 Phase V 接口说明；工作区无残留 `_scratch_*.py`。
