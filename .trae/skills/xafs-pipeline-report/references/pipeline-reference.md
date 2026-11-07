# Pipeline Reference — 执行层完整参考

SKILL.md 之外的完整签名、陷阱与诊断流程。按阶段组织。

## S1 读数（xafs_io）

| 函数 | 说明 |
|------|------|
| `read_scan(path)` | 自动判格式：.000/.100 XDAC、XDI、纯文本 |
| `read_xdac` / `read_xdi` / `read_simple` | 底层解析器 |
| `_safe_float` | 容错转换（"7112.0 eV" 这类带单位值） |

陷阱：
- XDAC 数据行 ≥3 列即收（3 列=无参考通道变体，01_Fe_foil 就是）。
- 2 列纯文本视作 (E, μ)，μ 中位数必须在 0.001–50，否则抛
  XAFSDataError 提示改 3 列 (E, I0, It)。
- RawScan 字段：energy/i0/it/ifluor/iref/dtime/mu/meta。

## S1.5 对齐与合并（athena）

```python
align_by_reference(scans, edge_energy, tolerance=0.2)
    # 抛 XAFSDataError 的两种情况：无 Ir 参考通道 / 参考边位与样品差 <3 eV
    # （= 参考测的是同一样品，本数据集即此情况）
merge_scans(scans)   # 网格一致直接平均；不一致插值到 0.2 eV 公共网格
```

**手写脚本必须自己包 try/except 降级**（run_xafs.process_sample 已内置）：

```python
try:
    scans, info = athena.align_by_reference(scans, 7112.0)
except XAFSDataError as exc:
    print(f"参考通道不可用，跳过对齐直接合并（{exc.message}）")
merged = athena.merge_scans(scans)
```

## S2 E0（athena）

```python
determine_e0(scan, method="derivative")   # 可选 half_height / manual
# 返回 dict：e0_sample / e0_reference / e0_derivative_peak / chemical_shift
```

- derivative 默认；half_height 会低约 3 eV。
- e0_reference 仅用于判断参考通道独立性，不参与 k 轴。

## S3 归一化（athena）

```python
normalize_absorption(scan, e0, pre_window=(-150,-30), post_window=(150,600), engine="larch")
optimize_rbkg(norm, rbkg_range=(0.7,1.6), n_grid=10)   # 返回 (best行, 全表)
```

- engine="larch" 与 Athena 数字对齐（edge_step 用二次多项式外推）；
  "native" 边后平台更平但 edge_step 小 ~10%。
- best 行字段：rbkg / r1 / amplitude / leakage / ratio / score /
  confidence（high|medium|low）。

## S4–S5 χ(k) 与 FT（athena）

```python
extract_chi_k(norm, rbkg=1.2, kweight=2)     # larch AUTOBK
ft_to_r_space(k, chi, r=None, kmin=3.0, kmax=12.0, kweight=2, dk=1.0)
```

**r 参数被忽略**——R 网格由 xftf 内部决定，二次插值到更细网格会把
主峰抹平并压低峰位（实测 1.473→1.440）。自写 FFT 的相位因子若误用
e^{−2πikR}，R 轴压缩 π 倍，Fe–S 1.82 Å 会跑到 0.58 Å。

FourTransform 字段：r / chi_r（复数）/ magnitude / real / imag。

## S4.5 荧光自吸收校正（可选，athena）

```python
booth_self_absorption(chi_k, formula, elem, edge, line, e0,
                      thickness_um, density, anginp=45.0, angout=45.0)
# 返回 (新 ChiK, info)；info 含 corrected/reliable/correction_factor_median
```

- 仅荧光模式；透射谱没有可反推的自吸收。
- **五项输入缺一不可**：样品总化学式（含基体）、厚度 μm、有效密度
  g/cm³、入射/出射角（相对表面）、荧光线（默认 K→Ka / L3→La）。
- 数学：Booth & Bridges 2005 Eq.4，在 k 空间闭式反解；只改振幅，
  不改相位。`mu_T`/`mu_a` 用 xraydb 光滑截面；`mu_a` 扣边前基线
  （`cross_section(E)-cross_section(E_edge-200eV)`），能量锚定表值边
  `E_edge+0.1`，与 Athena/Demeter 实现一致。
- 厚度自动落极限：d→0 严格 χ_corr=χ_exp；d 大于信息深度进厚样分支。
- 合成自洽回收 max|Δχ|=0.0015、rms=0.00050；薄极限逐点恒等。
  全链合成（真实 μ → Eq.2 前向 → 归一化/autobk → Booth）：
  χ 振幅比 0.723 → 1.054（rms=0.020）。
- 边界：均匀平表面假设；颗粒/孔洞/化学式错误会让有效 d、ρ、μ 失真。
  自吸收还会压低/移动 E0（实测 7111.3 vs 表值 7112 eV），校正不修 E0，
  拟合前必须用参考箔或已知边能量校准。`info["reliable"]=False`
  （截断/非有限值）时 CN 和相分数仍拒标。
- 不替代透射谱；有同样品透射时优先做独立对照。

## S6 小波（larch_wt）

```python
continuous_wavelet_transform(chi_k, sigma=1.5, kmin=3.0, kmax=14.0,
                             kweight=2, rmax=4.0, dk=0.1,
                             normalize=True)   # False 时保留原始幅度，供跨图比较
locate_shell(wave_map, r_lo, r_hi)          # 宽窗 argmax
locate_shell_near(wave_map, r_center, half_width=0.5)  # 锚定 FT 峰
extract_k_profile(wave_map, r_center, half_width=0.25)  # B3-1 一维切片
compute_model_wt_residual(chi_k_data, model_k, sigma=1.5,
                          kmin=3.0, kmax=12.0, kweight=2,
                          rmax=4.0)          # B3-1 实验/模型二维残差
```

- sigma ≥ 1.2（=1.2 以下 Fe–O 与 Fe–S 粘连；1.5 是推荐值）。
- kmax 入口已钳制到 max(chi_k.k)。
- 零载频高斯包络（非经典 Morlet）——R 轴与 |χ(R)| 严格对齐。
- **管线用 `locate_shell_near`**：以 FT 第一峰为中心 ±0.5 Å。宽窗
  [1.0,2.6] 在多壳体系会同时装进两条分支，argmax 会跳壳（Fe3O4 实测：
  S4.5 把低 k 抬高 7% 后 WT 从 2.60 跳到 1.42，FT 峰不变）。锚定后
  只回答“这个 FT 峰来自哪个 k”，不独立重新选壳。

### B3-1 一维壳层切片与模型残差

`extract_k_profile` 在 `[r_center−half_width, r_center+half_width]`
柱带内取 `ridge_profile = mag[:, cols].max(axis=1)`（嵴线）和最近
R 单行的 `slice_profile`，返回 `k_peak` / `r_track` / `r_peak` /
`slice_r`。双峰判据：次峰 ≥20% 主峰且间距 ≥1.5 Å⁻¹ → `bimodal=True`
（同一 R 窗叠加轻元素 M–O 与 M–M 的证据）。

`compute_model_wt_residual` 用**完全相同**的
`(sigma, kmin, kmax, kweight, rmax, normalize=False)` 分别算实验与
模型 WT，统一除以实验峰值后相减；二维
`R_WT = Σ(|W_data|−|W_model|)² / Σ|W_data|²`。两侧禁止各自归一化，
否则残差图失去物理标度。

`xafs_diagnostics.write_wt_shell_profiles(res, fit_out, outdir, stem)`
输出 `<stem>_wt_profiles.png`（N 个壳层切片 + 二维残差面板）与
`<stem>_wt_profiles.csv`（长表：`data_ridge/model_ridge/data_slice/
model_slice`）。壳层 R 用拟合键长经 `_ft_anchor_near_ft` 映射到
FT 表观峰（Fe–O 2.00→1.47、Fe–Fe 3.07→2.70），不手填 R 窗。
实测 LEPIDO：kmin=2.5 时 Fe–O `k_peak=5.5`、Fe–Fe `k_peak=6.4`
（差 0.9 Å⁻¹），2D WT R-factor=0.046；产线 kmin=3.0 时差约 0.7，
排序不变。

## S7 CIF → FEFF → 拟合（artemis）

### CIF 读取
```python
read_cif_sites(cif_path, absorber="Fe", rmax=6.0, ignore=("H","D"))
    # 返回 [(x,y,z,元素), ...]（笛卡尔 Å，吸收原子在原点）
```
- gemmi 对称展开 + 周期像补全。
- **从文献 PDF 抄坐标必须自检**：展开后近邻键数与键长对照文献
  （例：γ-FeOOH Fe-O 必须 6 键 1.98×2+2.07×4；若只得 4 键说明坐标
  串行/错位，弃用）。
- 本仓库已有 `C:\Users\lizha\Desktop\Code-experience\lepidocrocite.cif`
  （Zhukhlistov ED 数据 + O1 几何补全，头注有推导记录）。
- 部分占据与共占位（B7-2）：`return_occupancy=True` 返回 `(x,y,z,元素,occ)` 5 元组（默认仍为 4 元组）；距离 `<0.08 Å` 的共占位自动合并，且 `xafs_fit._ss_rows` 通过 `cif_shell_occupancy_scale` 将 FEFF 几何简并度按壳层平均占据度缩放为 `N_geom·mean(occ)`。

### FEFF
```python
build_feff_input(sites, absorber="Fe", rmax=5.0, title="...", edge="K")
run_feff6l(inp_text, workdir)          # 返回绝对路径列表
describe_paths(files, n)               # [{file, reff, degeneracy, legs}, ...]
```
- `edge` 决定 HOLE 卡：K=1 / L1=2 / L2=3 / L3=4（FEFF8 约定）。
  非 K 边使用独立 FEFF 缓存键（`lib_<motif>_<hash>_l3`），避免复用 K 边路径。
- L 边接入（Pt L3 已验收）：absorber.py 配 `next_edge=("L2", 13273)`，
  kmax_ceiling 据此估算（Pt 预置 16.0）；L₃ 白线强、能量窗口宽，
  Pt 箔实测 S3 edge_step=0.844、平台 0.945（strict 可通过）。
  分子配合物若库内无条目会正确 no_match → L1，配体必须显式给
  （`--l1-ligand Cl/N`；默认 O 只适用于氧化物）。

## S6.5 识别与 S₀² 标定（structure_id / calibrate_s02）

完整方法、覆盖清单与诚实边界见 `structure-identification.md`；执行层要点：

```python
structure_id.identify(ft, k, absorber_symbol="Fe", kmin=3.0, kmax=12.0,
                      kweight=2)   # → ranked / indistinguishable / uncovered_peaks
structure_id.identify_result(res, absorber_symbol="Fe",
                             vary_c3=False)   # → status + L1 降级
# vary_c3=False/True/"first"/"all"/[idx]；L1 降级时同样生效
```

- 未标定 S₀² 的元素在 identify 入口即拒绝（只排序不定量的结果是伪定量）。
- 判定四态阈值：`CONFIDENT_SCORE=0.16` / `CONFIDENT_MARGIN=0.03`；
  非 confident 默认 L1 第一壳，`--allow-ambiguous` 才强制 Top-1。
- **mixture 候选**（`kind="mixture"` 多位点联合，χ = Σ w·χ_位点）：覆盖判定按
  "复合谱 ∪ 各分量谱分别解释"——复合谱中分量峰受干涉位移/压制（实测
  磁铁矿 tet 位 Fe-Fe@3.48 峰在复合谱中相对峰高 14.7% < 15% 被过滤、
  1.35 Å 只余斜坡），只看复合峰会把正确联合模型错标 partial。
  拟合侧 `_mixture_rows` 按位点权重缩放简并度，合并容差收紧 ≤0.10 Å
  （磁铁矿 Fe-O 1.886/2.060 Å 只差 0.17 Å，0.20 会错并成一层）。
- **mixture 相分数**（`--fit-weights`）：识别阶段权重固定；定量阶段
  `select_shell_paths(..., group_by_component=True)` 保留**未乘权重**的
  晶体学简并度并按组分分组，`fit_r_space(vary_n=False, vary_weights=True)`
  用 softmax 放开 K-1 个 logit（恒正、和=1），有效简并度 = N_crystal·w。
  禁止同时 `vary_n`（N·w 严格简并）；权重误差由协方差 + softmax 雅可比
  做 delta method 传播，`PathFit.n_err = N_crystal·weight_err`。
  输出 `FitResult.mixture_weights` / `weights_reliable`；不同组分在同一
  壳层不可分辨时标记 `flags.weight_degenerate` 并拒绝相分数结论。
  合成磁铁矿 0.25/0.75 回收误差 <1e-6（回归测试 13）；真实数据
  `Fe-Magnetite.xdi` 默认窗 Rf=0.097 → `weights_reliable=False`。
- **参数相关矩阵**：`FitResult.parameter_names` / `FitResult.correlation`
  随每次 S7 拟合输出（`<stem>_fit.json` 的 `fit.correlation`），报告层
  额外给出 `max_parameter_correlation`。协方差是 `pinv(JᵀJ)·(2cost/dof)`，
  相关矩阵按 `cov_ij/(σ_iσ_j)` 归一；强相关（实测磁铁矿 ΔE0 vs
  oct Fe–Fe ΔR = 0.904）说明单个误差棒不能独立解读。
- **荧光自吸收筛查 / 可选校正**：`xafs_fit._self_absorption_diag` 从
  `merged.meta["mode"]` 与 `norm.edge_step` 判定风险；默认**不做定量
  校正**（缺少组成/密度/厚度/几何）。契约：荧光 `Δμ≥1.0` →
  `risk=high`、`cn_gate=True`，CN/相分数结论被拒绝；`0.3–1.0` →
  `possible` 仅警告；`<0.3` → `low`；透射/`given` → `none`；荧光但
  Δμ 缺失 → `unknown` 且不放行。若用户提供总化学式+厚度+有效密度+
  角度（+荧光线），S4.5 按 Booth 2005 Eq.4 在 k 空间校正，`risk=
  corrected`；校正不可靠时 `risk=corrected_unreliable` 仍拒标。
  结果随 `<stem>_fit.json` 的 `self_absorption` 输出。

### S₀² 标定工具（calibrate_s02.py）

```bash
python calibrate_s02.py <箔数据> --absorber Co --no-strict -o xafs_out/_cal_co   # 干跑
python calibrate_s02.py <箔数据> --absorber Co --no-strict --apply              # 通过门禁才写入
python calibrate_s02.py ... --force     # 覆盖已有标定值（重标定是显式决定）
```

协议 = **生产等效拟合**：N 固定晶体学简并度，S₀² + 每路径 ΔR/σ² +
全局 ΔE₀ 自由、kw=2、k∈[3,12]、与生产拟合同一 FT/R 空间约定——标定值
写入后，生产拟合（变 N）的 CN 直接复现晶体学简并度。数学上
"定 N 变 S₀²" ≡ "定 S₀² 变 N"（振幅只以乘积出现），但数值上必须同构。

- **3 窗口鲁棒性**：rmax = 第一壳 Reff + 0.12/+0.18/+0.26 各拟一次，
  S₀² 取中位数；`window_spread`（极差）> 0.05 拒绝。历史教训：单一
  窗口的 4 自由变量窄窗拟合曾给出跨配置 8%–崩塌级（撞参数边界）
  的抖动，"恰好稳定"的窗口不可作为生产值。
- 门禁：s02∈[0.4,1.15]、s02_err<0.06、Rf≤0.12、|ΔE₀|<10 eV、
  vars<N_idp（逐窗口）、全窗口 interior。
- **kw∈{1,3} 拟合只做诊断记录（`kw_sensitivity` 字段），不门禁**：
  kw=2 才是生产配置；kw≠2 下实测箔样撞 s02 上下界（1.3/0.3）、
  Rf 数倍膨胀——那是数据-FEFF 形状失配的放大镜，不携带生产语义。
- 金属箔 plateau 0.79–0.83（强 EXAFS 振荡 + E⁻³ 慢衰减样条雏形）
  被 strict 拒绝是预期行为；标定与箔识别都必须 `--no-strict`，
  数值可靠性由拟合层门禁兜底。
- 已执行（2026-11）：Co=0.7308（极差 0.0005）、Ni=0.7795（0.0026）、
  Cu=0.8670（0.0092）、Pt=0.8075（L3，0.0026，Rf=0.0004，ΔE₀=+8.25）。
- **Fe 箔负例**：fe.060（XDAC）与 Fe-K-IronFoil.xdi（XDI）在标定协议下
  均撞 ΔE₀=+15 eV（开 C₃ 也不解）→ 厚样谐波污染嫌疑，门禁正确拒绝；
  Fe=0.803 保留历史值——原始推导不可复原，由氧化物价态验收结果锚定。
- 标定结果存 `xafs_out/_cal_<elem>/calibration.json`（runs/kw 诊断/门禁全留档）。
- **describe_paths 返回的 file 是相对文件名**（feff0001.dat），拼路径
  必须加 FEFF 工作目录，否则 feffpath() 抛 "not found"。

内置生成器：`bcc_lattice_sites(a, rmax)`（Fe 箔）、
`octahedron_sites(bond, symbol)`（无结构时第一壳层）。

### 拟合
```python
fit_r_space(chi_k, ft, datfiles, rmin, rmax, kmin=3.0, kmax=12.0,
            kweight=2,             # 或 (1, 2, 3) 多权联合
            s02=absorber.s02, vary_delta_e0=True, vary_n=False,
            vary_s02=False,        # 仅标定模式（calibrate_s02.py）
            vary_sigma2=True,      # False 时 σ² 固定在 sigma2_init
            vary_c3=False,         # True / [路径下标]；边界 ±C3_LIMIT
            sigma2_init=0.005, deltar_lim=(-0.3, 0.3),
            path_degen=None, path_labels=None,
            sigma2_groups=None, ms_links=None,
            vary_weights=False, weight_init=None)
```
- vary_n：False/None=全固定；True=全部路径放开；[0,2]=只放开指定下标。
- **B1-1 多 kweight**：`kweight=(1,2,3)` 时每个分支的残差先除以该
  分支实验谱 RMS 范数 `s_w`，再相加；`N_idp` 不乘权重个数。
  单权路径与历史逐位一致。
- **B4-2 C3**：`vary_c3=True/[0]` 放开三阶累积量，边界
  `±artemis.C3_LIMIT=±0.01 Å³`；`|C3|≥0.008` → `c3_at_bound=True`，
  上层 `cn_reliable=False`。C3 与 `vary_weights` 互斥。
- `ms_links`：B1-3 受约束 MS 的母壳层绑定；只选 `nleg∈{3,4}`、
  `reff≤rmax+0.15`、所有中间原子 `|cos|≥0.85`、且
  `reff≥1.5×第一壳距离` 的近共线路径，不新增自由变量。
- vary_s02：生产拟合**禁止**开启（S₀²·N 简并）；仅 calibrate_s02
  标定用。s02=None 且 vary_s02 时初值 1.0。FitResult 新增
  `s02`/`s02_err` 字段（s02_err 仅标定模式非 None）。
- vary_sigma2/sigma2_init/deltar_lim：σ² 可按路径固定（标定协议用，
  把窄窗过参数化的 4 变量拟合降为 (S₀²,ΔE₀) 近线性）；sigma2 误差
  在固定路径上为 None。
- `path_degen`/`path_labels`：合并壳层的代表路径可显式给总简并度与
  可读壳层名；拟合出的 N 是壳层总配位数。
- N 的实现是每轮迭代改 path.degen 属性（larch _calc_chi 直接读它）。
- 返回 FitResult：paths（PathFit 列表）/ s02(_err) / delta_e0(_err) /
  r_factor / n_independent / n_variables / overfit /
  **epsilon_r / chi_square / reduced_chi_square** /
  **model_k / model_k_grid**（B3-1 WT 残差与反 FT 复用）/
  **mixture_weights / parameter_names / correlation**。
- PathFit：reff / degeneracy（放开时即拟合 N）/ delta_r(_err) /
  sigma2(_err，固定路径为 None) / n_err（未放开为 None）/ r
  （=reff+delta_r）/ **c3 / c3_err / c3_at_bound**（B4-2）/
  **nleg / ms_parent**（B1-3）/ **sigma2_at_bound**（σ²≤5e-4·Å²
  视为撞下界——形状失配把 σ² 往下压而非吸收的诊断旗，见 goethite
  案例：Fe-O@3.21 σ²=0.0001 触发；Fe-Fe@3.02 的 0.0007 不触发，
  只作辅助观察 → cn_reliable=False 的证据链之一）。
- **B4-1 统计量**：`estimate_r_noise` 对
  `|χ(R)|` 在 `[15,25] Å` 的实部/虚部取均方根得 ε_R；必须显式
  `xftf(rmax_out=25)`（larch 默认网格只到 ~10 Å）。χ² 用主 kweight
  分支、`N_idp` 不乘权重个数：
  `χ² = N_idp·Σ|χ_data−χ_model|² / (2·N_grid·ε_R²)`，
  `χ²_ν = χ² / max(N_idp−N_var, 1)`。
- **B4-1 壳层 F-test**：`shell_f_test(fit_simple, fit_complex)` 要求
  两个模型同一 R 窗、同一数据；
  `F = [(χ²₁−χ²₂)/(p₂−p₁)] / [χ²₂/(N_idp−p₂)]`，
  `confidence = F.cdf(F; dfn, dfd)`；`confidence≥0.95 且 dfd≥1`
  才判新增壳层显著。自动 S7 用第 1 壳对照生成 `shell_f_test`，
  `delta_chi_square = χ²_simple − χ²_complex`（>0 表示多壳改善）。
- **n_err 只是形式误差**：不含 S₀² 传递、代表路径/壳层合并、
  N–σ²–ΔE₀ 相关与归一化偏差。按 Fe 基准实测，CN 绝对值系统误差
  约 ±20–40%（1σ 量级估计），报告必须注明，不能把 ±0.3 当绝对精度。
- `xafs_fit` 输出 `cn_reliable`：Rf>0.08、过拟合、|ΔE0|≥15 eV、
  N 形式误差过大（>1 或 >50%）或覆盖清单未通过时为 False；
  此时拟合表标题会标注 “CN unreliable”。实测 Cu 箔第二壳
  N=18.0±1.6（上界附近）→ 正确判 False。

### 自动 S7（xafs_fit）
`xafs_fit.run(res, fit=True)` 在识别 `confident` 时自动执行：
`select_shell_paths` 按 (散射元素, Reff) 分壳层（合并容差 0.20 Å，
小于 π/(2Δk)≈0.26 Å），每条壳层取简并度最大的 SS 路径作代表，
`path_degen` 传入壳层总简并度，`path_labels` 写成 `Fe-O@2.04` 形式；
默认 R 上限 3.5 Å，避免把 3.6–4.4 Å 的弱远壳全部塞进自动拟合。
自动选壳层还受信息量限制：`3n+1 < N_idp`，超出时从最远壳层截断，
并写 `truncated_shells`。非 confident 时返回 L1 第一壳。
`--fit-motif` 可强制指定结构，`--allow-ambiguous` 才允许 Top-1 全拟合。

### L1 第一壳（无结构时的诚实降级）

`structure_id.l1_first_shell_fit` 用合成团簇（oct 6 / tet 4 个配体）
拟 N/R/σ²/ΔE₀。两条实测规则：

- **配体感知选峰**：r0 不取"最低显著峰"也不取 argmax，而是用共价半径
  和（`artemis.COVALENT_RADII`）估预期键长、减 0.45 Å 相位偏移，在显著
  峰里选最近者。反例：Fe 氧化物里 Fe–Fe 第二壳（2.55）比 Fe–O 第一壳
  （1.35）更强；Pt L3 白线在 1.2–1.4 Å 产生显著伪峰——两种朴素规则
  都会选错。缺半径时退回最低显著峰。
- **R 窗上限随表观峰扩展**：`min(3.0, max(2.3, peak+0.15))`。固定 2.3 Å
  会把 Pt–Pt（2.77）/Fe–Fe（2.48）的模型峰截在窗外，曾给出 Rf=1.0 的
  伪失败。
- `ligand`/`geometry` 是**显式假设**，默认 O/oct 只适合氧化物；氯化物、
  硫化物、金属必须覆盖（`--l1-ligand Cl`、`--l1-geometry tet`）。输出
  保留假设，`cn_reliable` 按 Rf≤0.08 且 |ΔE₀|<15 判定。

### 模型选择决策树（实测经验）
1. 列全部 SS 路径（legs==2），按 Reff 排序。
2. ΔR < 0.26 Å 的相邻路径合并为一条有效壳层，N 语义变为该壳层总配位数。
   **代表路径一经基准验证即固定**（γ-FeOOH：Fe-O→2.072、Fe-Fe(1)→3.062、
   Fe-Fe(2)→3.873；Fe 箔→2.482/2.866）——同壳层内换代表路径会改变
   F(k)/φ(k) 形状，N 与 Rf 随之漂移（实测 Fe-Fe(2) 换 3.87→3.93，
   N 从 1.0 漂到 3.8）。xafs_report.pick_paths 会对每条路径校验 Reff。
3. 首选 3–4 条有效壳层 + k=3–12。变量数必须 < N_idp。
4. R-factor 仍 > 0.08 时才考虑加 MS/聚焦路径；加了必须复查各 N 误差
   是否爆炸（简并信号），爆炸即回退。
5. 拟合区间：rmin/rmax 盖住目标壳层。示例数值是 Fe 基准
   （lepido 1.0–3.9；Fe 箔 1.8–3.1）——**不是普适**：Pt L3 金属
   第一壳 Pt–Pt 2.77 Å → 表观 ~2.45，窗口应为 (1.9, 3.2)；区间
   按元素/结构的第一壳键长设置。

### 诊断脚本模板（临时文件 _diag.py，跑完删）
```python
from absorber import get_absorber        # s02 一律查标定表，不再硬编码
for name, kmax, files, vn in configs:
    ft = athena.ft_to_r_space(chi_k.k, chi_k.chi, None, 3.0, kmax, 2)
    fit = artemis.fit_r_space(chi_k, ft, files, 1.0, 3.9,
                              kmin=3.0, kmax=kmax,
                              s02=get_absorber("Fe").s02, vary_n=vn)
    print(name, fit.r_factor, fit.delta_e0,
          [(p.degeneracy, p.n_err) for p in fit.paths])
```

## 已知数值基准汇总

| 量 | 基准值 | 数据 |
|----|--------|------|
| E0 (derivative) | 7128.60 eV | lepido 3 扫描合并 |
| edge_step | 0.4695 | 同上 |
| Rbkg（optimize 选中） | 1.20, confidence high | 同上 |
| FT 第一峰 | 1.473 Å (k=3–12) | 同上 |
| WT 第一壳 | k≈5.7, R≈1.46 Å | 同上 |
| Fe 箔拟合 | N=6.8±0.7/7.4±0.8, R=2.447/2.822, Rf=0.0037 | fe.060/150/300 |
| lepido 拟合 | N=6.5/7.3/1.0, R=1.994/3.052/3.979, Rf=0.0737 | 3 壳层, k=3–12 |
| S₀²（absorber.py 标定表） | Fe 0.803※ / Co 0.7308 / Ni 0.7795 / Cu 0.8670 / Pt 0.8075(L3) | 真实箔 calibrate_s02.py；窗口极差 ≤0.0092 |
| Pt 箔 L3 拟合 | N=11.69±0.22, R=2.764 Å, Rf=0.0038, ΔE₀=+8.12 | Pt-L3-PlatinumFoil, pt_fcc |
| Pt–Cl L1（ligand=Cl） | N=3.98±0.14, R=2.302 Å, Rf=0.0228 | K₂PtCl₄（平方平面 CN=4） |
| LEPIDO 默认拟合（B4-1） | Rf=0.0273, N=6.35, R=2.002 Å；ε_R=0.009899, χ²_ν=314.2 | 3 扫描合并；F-test F=32.34, confidence=0.9991 |
| LEPIDO 单扫描拟合（B4-1） | ε_R=0.01064, χ²=1405.43, χ²_ν=273.32 | F=32.67, confidence=0.99909 |
| LEPIDO WT 切片（B3-1） | Fe–O k_peak=5.5, Fe–Fe k_peak=6.4 Å⁻¹（kmin=2.5） | 2D WT R-factor=0.046 |
| LEPIDO C3（B4-2） | C3=+1.0315e-3 Å³, R=2.046 Å, CN=5.5 | `--vary-c3 first`；默认 R=2.002/N=6.35 |

※ Fe=0.803 为历史标定（原始推导不可复原）；fe.060 与
Fe-K-IronFoil.xdi 在现行协议下被门禁正确拒绝（厚样谐波嫌疑）。
另：Fe 箔拟合基准行是"每条路径各自放开 N"的旧模型数值（壳层合并
铁律落定前），两路径 N=6.8/7.4 对 hcp/bcc 均无物理意义，仅供对照
拟合链的历史行为；现行基准见 Co/Ni/Cu 箔生产拟合（第一壳
N≈11.6–11.9 vs 晶体学 12）。

## 工程化与诊断增强（P1，2026-11）

### 配置文件（xafs_config.py）
支持 `.yaml` / `.json`，CLI 参数优先级高于配置文件。仅在 CLI 仍为代码
默认值时才被配置文件覆盖。

| 配置项 | CLI 参数 | 默认值 |
|--------|---------|--------|
| kmin/kmax/kweight | `--kmin/--kmax/--kweight` | 3.0/12.0/2 |
| post_window | — | [150, 600] |
| rbkg | `--rbkg` | None（自动寻优） |
| strict | `--strict/--no-strict` | True |
| wt_sigma | — | 1.5 |
| despike | `--despike` | False |
| self_absorption.* | `--self-absorption-*` | 各 null |
| batch.pattern/outdir | `--batch-dir` | `*.xdi` / `xafs_out_batch` |

用法：
```powershell
python run_xafs.py <谱文件> --absorber Fe --config xafs.yaml
python xafs_report.py <谱文件> --absorber Fe --config xafs.yaml --no-strict
```

### 批处理（`--batch-dir`）
对目录下所有匹配文件循环处理，单个失败不中断整批，失败记录写入
`<outdir>/batch_errors.json`。实测 BMM Fe 标准 38 个文件 38/38 成功。
与 `files` 位置参数互斥。

### 去毛刺（`--despike`）
S1.5–S3 之间插入中值滑窗（窗口=5）+ MAD 门禁（阈值=3×MAD），默认关闭。
实测 lepidocrocite 替换 79 个离群点，E0/edge_step/Rbkg 结果与未去毛刺一致。

### FT–WT 联合诊断图（`s7_ft_wt_joint.png`）
`--diagnostics` 现输出 7 张阶段图（原 6 张 + 联合图）。左图 `|χ(R)|` 标注
FT 第一峰 r1；右图 WT contour 叠加水平线 R=r1 与垂直线 k=k_wt，交点标出，
直接回答「这个 FT 峰主要来自哪个 k」。对标 Larch-WT 的交互诊断，避免
WT 在多壳分支间 argmax 跳壳（Fe3O4 实测：Booth 校正把低 k 抬高 7% 后
WT 峰从 2.60 跳到 1.42，FT 峰纹丝不动；联合图可肉眼核对）。

### σ² 共享组（对标 Artemis 路径约束，`fit_r_space` 新增）
`--shared-sigma2 all` 让所有路径共享一个 σ² 变量，减少自由参数、避免
σ² 与 ΔR 过度相关。或用 JSON 列表显式分组，如
`--shared-sigma2 '[[0,1],[2]]'`。共享组输出相同的 σ² 和 σ²_err。
与 `--fit-weights` 互斥（权重模式下 σ² 已按组分缩放，不再共享）。

FitResult 新增字段：
- `sigma2_groups`：[[0, 1], [2]] 形式的路径下标分组
- `sigma2_group_names`：["all_shells", "s2grp1"] 等组名称

xafs_fit 输出 `sigma2_shared` 标志，报告层据此标注哪些壳层共享 σ²。

### 单文件 HTML 报告（xafs_report.export_html_report）
每次报告默认生成 `<stem>_report.html`（`--no-html` 关闭）。图片以 base64
内联，HTML 可脱离目录独立发送；强相关（|ρ|≥0.8）与不可靠 CN/相分数标红；
关键指标、flags、拟合表、诊断图（折叠）与 CSV 文件链接同页呈现。

### 动态 CN 系统误差预算（xafs_fit.systematic_budget）
基础预算 0.30 不再写死。额外项按平方和合成：
`self_absorption_residual`（Booth 校正可靠 +0.10）、
`self_absorption_uncorrected`（risk=possible +0.15）、
`parameter_correlation`（|ρ|≥0.9 +0.10；0.8–0.9 +0.05）、
`truncated_shells`（+0.05）、`sigma2_at_bound`（+0.10）。
输出 `cn_systematics.relative_1sigma` / `extra_components` / `note`；
相分数另有 `weights_systematics`（模型绑定 0.20 + 拟合权重误差）。
实测：lepidocrocite ρ=0.913 → ±32%；磁铁矿 SA 合成谱
（Booth + ρ=0.999 + 截断 + σ² 撞界）→ ±35%。

### S8 WT 对比与交互式 k 窗口
`xafs_diagnostics.wt_sweep(res)` 在 FT 第一峰 ±0.5 Å 锚定下扫描
σ∈{1.0,1.5,2.0} × kweight∈{1,2,3}（`normalize=False`），返回
`rows`（9 行）/`maps`/`row_max`/`anchor_r`。`write_wt_sweep` 与
`--diagnostics` 的 `_s8_wt_sweep` 共用 `_render_wt_sweep`：色标按
**kweight 行**归一（kweight=3 幅度约是 kweight=1 的 k² 倍，全局色标
会把低 kweight 压成深蓝），k 轴用全部面板并集（`sharex` 下逐面板
`set_xlim` 会被最后一个面板覆盖）。`--wt-sweep` 另存
`<stem>_wt_sweep.png/csv`。

`xafs_window_tuner.write_window_tuner` 生成 `<stem>_window_tuner.html`：
浏览器端 FT 与产线同用 Hanning + e^{−2ikR}（`hann_window` 是契约测试用的
同一窗函数），幅度按各自最大值归一，只比较峰位/形状。
`recommend_window` 用 `larch.xftf` 扫描候选网格，只保留峰位落在
r1±0.15 Å 的候选，`score = amplitude/(1+20·stability+5·|R−r1|)`，
给出参考推荐，**不自动改产线**。`--no-window-tuner` 关闭。

### B2-1 XANES flattened / 颈部特征
`athena.extract_xanes_features(norm)` 在 `[E0+100, E0+400]` 平台均值
归一得到 flattened 谱；larch 分支用 `g.flat`，native 分支对
`mu_norm` 拟合二次背景后同样归一。返回平台值、白线强度/位置与
颈部双台阶特征。`run_xafs`/`xafs_report` 输出 `<stem>_xanes.csv`。
这是补充表征，不参与 EXAFS 定量门禁。

### B2-2 LCF / PCA（xafs_lcf）
```python
run_lcf(target_x, target_y, standards, xmin, xmax, mode="xanes",
        sum_to_one=True)   # → weights / fit / r_factor / residual
run_pca(spectra, xmin, xmax)  # → components / explained_variance / ...
lcf_from_scans(target, standards, mode="xanes")  # RawScan → LCF
```
端元权重是线性组合系数，不是独立物相定量；报告必须附端元选择
依据。CLI：`--lcf-standards a.xdi,b.xdi`。

### B2-3 死时间 / deglitch
```python
apply_deadtime_correction(scan)   # 仅荧光模式
deglitch_ranges(scan, [(7250.0, 7262.0), ...])  # 直接删点，不插值
```
`dtime≥1` 为倍乘因子，`0<dtime<1` 为死时间分数；透射模式跳过。
deglitch 删除比例 >15% 或边核心剩余 <5 点抛
`XAFSConfigError("S1.2", ...)`。CLI：`--deadtime`、
`--deglitch '7250:7262,7410:7415'`。

### B4-2 C3 全链
`fit_identified(..., vary_c3=False|True|"first"|"all"|[idx])` 与
`structure_id.l1_first_shell_fit(..., vary_c3=...)` 均已接入；
CLI `--vary-c3 first|all`。放开后若 `n_var ≥ N_idp`：先降级为
仅第一壳，仍超限则关闭 C3 并写 `c3_dropped_by_nidp=True`。
`|C3|≥0.008 Å³` → `c3_at_bound=True`、`cn_reliable=False`。
`PathFit.c3/c3_err/c3_at_bound` 写入 `_fit.json`；拟合表仅在有
C3 时追加 `C3 (Å³)` 列（默认仍 7 列）。合成 C3=1.5e-4 在
rmax=2.6 精确回收；rmax=3.5 回收 1.808e-4（误差 3.08e-5）。

### B5-1 复数小波 / 相位切片
`larch_wt.continuous_wavelet_transform` 返回 `complex_map`、`phase`；
`WaveletMap.magnitude` 与 `|complex_map|` 在 `normalize=True` 时逐位一致。
`larch_wt.compute_complex_wt_residual` 同时比较实部/虚部：
反相模型 `χ_model=−χ_data` 的 `complex_r_factor=4.0`、幅值残差=0.0。
`xafs_diagnostics.write_wt_shell_profiles` 的 CSV 增加
`data_phase_ridge`，并返回 `complex_r_factor`。

### B5-2 元素理论 k_peak 标尺
`structure_id.reference_element_k_peaks(absorber, elements=(...),
r_bond=..., sigma2=...)` → `{element: k_peak}`，FEFF 结果带缓存。
实测 O/Cl/Fe = 4.9/5.6/7.0 Å⁻¹；Fe 的 σ² 0.003→0.015 时峰位左移
1.4 Å⁻¹。仅作 WT 元素判据，不替代 S7 拟合。

### B6-1 多扫描离群剔除与 σ_μ
```python
athena.merge_scans(scans, reject_outliers=False, z_threshold=3.5)
# RawScan.sigma_mu: 保留扫描的逐点标准差
```
默认关闭，输出与历史均值逐位相同；开启后按中位谱 RMS+MAD 剔除，
至少保留 2 条。CLI：`--reject-outlier-scans`。

### B6-2 XANES 荧光自吸收
```python
athena.fluo_xanes_self_absorption(norm, formula, elem, edge, line,
                                  thickness_um, density,
                                  anginp=45.0, angout=45.0)
# μ_true = α·μ_exp / (α + 1 − μ_exp)，分母 ≥ 0.25(α+1)
```
只改 XANES，不动 `chi_k`。Fe₂O₃/Fe K/45°–20°：α≈0.673、
白线误差 0.0、LCF 权重 0.5/0.5；默认 45°/45° 的 α≈0.42
触发 74 点钳制并 `reliable=False`。

### B6-3 边前 Pseudo-Voigt 分峰
```python
athena.fit_pre_edge_peaks(norm)  # n_peaks/centroid_energy/r_factor
athena._pseudo_voigt_area(...)
```
合成双峰 7113.0/7114.6 eV → 质心 7114.1025 eV（真值 7114.0667，
误差 0.036 eV），Rf=0.00123。指标表新增 `pre_edge_centroid`、
`pre_edge_n_peaks`、`pre_edge_voigt_r_factor`。

### B7-1 各向同性晶格约束
```python
fit_identified(..., isotropic_deltar=True)
# 唯一全局槽位 (-5, "alpha")，ΔR_i = alpha·Reff_i，|alpha| ≤ 0.08
```
合成 NiO 四壳 α=0.004：变量 13→10、α 回收 0.00399999999995、
Rf=4.2e-23，四壳 `ΔR/Reff` 与 α 一致（<1e-12）。
实验 NiO 各壳自由位移不服从统一缩放（α=0.00091，Rf=0.135，
自由拟合 0.024），必须 `cn_reliable=False` + 警告；不能用该约束
下的 N 作配位数结论。CLI：`--isotropic-deltar`。

### B7-2 CIF 部分占据 / 共占位
```python
artemis.read_cif_sites(cif, absorber, rmax, return_occupancy=True)
# → [(x, y, z, symbol, occ), ...]；默认仍 4 元组
artemis.cif_shell_occupancy_scale(cif, absorber, rmax)
# → [(distance, element, occ), ...]
```
距离 <0.08 Å 共占位合并，元素取高占据者，有效占据度
`min(1, Σocc)`；吸收中心按 `occ×多重度` 选主位点。
`xafs_fit._ss_rows` 按距离/元素匹配平均占据度，把 FEFF 简并度
缩放为 `N_geom·mean(occ)`。合成验证：O occ=0.75 → 6×0.75=4.5；
Co/Ni 共占位 → Co occ=1.0。

### B7-3 C4 与 g(R)
```python
artemis.C4_LIMIT                      # 2e-4 Å⁴
artemis.set_path_parameters(..., fourth=c4)
artemis.fit_r_space(..., vary_c4=True)  # 槽位 (i, "c4")
artemis.reconstruct_pair_distribution(r_grid, n, r0, sigma2,
                                      c3=0.0, c4=0.0)
```
C4=6e-5 合成回收误差 6e-17。g(R) 面积严格=N（梯形重归一），
C3=2.5e-4 时 `r_mode=2.032 < r_mean=2.05`、`skewness=0.538`。
**注意**：规范原文“右半侧面积>左半侧”与归一化正偏态分布矛盾；
本实现用尾长判据（q95−r0=0.142 > r0−q05=0.118），并显式断言
左面积>右面积（3.159 vs 2.779）以锁定归一化正确性，不伪造面积。
`PathFit.c4/c4_err/c4_at_bound` 写入 `_fit.json`。

### B8-1 多数据集全局联合拟合
```python
artemis.HBAR2_OVER_2KB                 # 24.254316 amu·Å²·K
artemis.einstein_sigma2(T, theta_e, reduced_mass_amu, sigma2_static=0.0)
artemis.fit_multi_datasets(datasets, shared_bonds=None,
                           einstein_models=None)
# → datasets / n_independent_total / n_variables_total /
#   global_r_factor / einstein_results / parameter_names / correlation
```
每个数据集残差除以自身 RMS 范数 `s_d`；`shared_bonds` 绑定
ΔR/σ²/N（N 可按 degen_ratio 缩放）；`einstein_models` 用
σ²_static+Θ_E 驱动各温度 σ²(T)。合成 LEPIDO 第一壳双温度
80/350 K（真值 Θ_E=420 K、σ²_static=0.0015 Å²、ΔR=−0.02 Å、N=6）：
全局变量 8→5、Θ_E=419.99999999 K、σ²_static=0.0015、Rf≈1e-22。
`ΔE₀` 为所有数据集共享的全局变量；双吸收边需独立 ΔE₀ 时应拆分
为两次单边拟合，再用 R/σ² 约束对齐。Phase III 已补 `share_e0`，
见 B12-1。

## Phase III（B9–B12）签名与判据

### B9-1 配置全键 / C4 降级 / g(R)
```python
xafs_config.DEFAULTS                     # 27 键，含 B9–B12 新开关
xafs_config.apply_config_to_kwargs(config, kwargs, explicit=None)
xafs_fit.fit_identified(..., vary_c4=True|False|[0]|"all",
                        fit_space="r", bvs_restraint=None)
xafs_fit.rdf_table(fit)                  # {"rows", "summary"}
xafs_fit.write_rdf_csv(path, rdf)        # 7 列 long-form
```
优先级：CLI 显式 dest > 配置文件 > 代码默认。C4 因信息量上限被降级
时记 `flags.c4_dropped_by_nidp=True`（先降 C4、再降 C3，各自独立
标记）。实测合成 C4=6e-5 → 5.96e-5；LEPIDO 4 壳 `vary_c4="all"`
→ n_var=14 < N_idp=15.007，仍只保留第一壳 C4。g(R) 面积=N
（5.995/5.995、7.436/7.436）。

### B9-2 误差带 / SNR / WT 标尺 / L3 白线
```python
athena.high_k_snr(chi_k, sigma_mu=None, edge_step=None, kweight=2)
# → {"snr_k", "k_snr2", "reliable_kmax", "sigma_chi", "reason", ...}
xafs_report.auto_ref_k_peaks(res, fit, absorber)   # 带缓存的 dict
athena.extract_xanes_features(norm, edge="K"|"L3", absorber="Pt")
athena.fit_whiteline_area(norm, window=(-6, 40))
```
`σ_χ=σ_μ/|edge_step|`；无 `sigma_mu` 返回 None+reason，不编曲线。
噪声 1e-3/3e-3/1e-2/3e-2 → `k_snr2=13.0/9.5/4.5/4.15`。K 边有
`pre_edge_n_peaks`、无 `white_line_area`；L3 反之（Pt 箔
4.5073/5.3414/0.5427）。WT 标尺实测 `{'O':4.9,'Fe':7.0}`。

### B9-3 .prj / MCA
```python
xafs_io.write_prj(path, {"scans": [...], "standards": ..., "meta": ...})
xafs_io.read_prj(path, as_scans=False)
athena.merge_mca_channels(channels, dtime=None, z_threshold=3.5,
                          min_channels=2)
# → (summed, {"rejected", "n_kept", "weights", "scores", ...})
```
`.prj` 为 gzip-JSON（非 Athena Perl）；Athena 原生工程显式拒绝。
MCA 用通道间 MAD-z 剔除坏道，再按离散度加权；8 道 2 道 spike
→ `rejected=[2,6]`，合并谱=干净道加权和（逐位）。

### B10-1 MBACK / chi_std
```python
athena.normalize_mback(scan, e0, elem, edge="K",
                       pre_window=(-200,-30), post_window=(60,150),
                       fit_order=2)
athena.extract_chi_k(norm, rbkg=1.2, kweight=2,
                     chi_std=None, k_std=None,
                     clamp_lo=None, clamp_hi=None)
```
MBACK：`μ_meas≈a·μ_th+b0+b1(E−E0)+…`，`a·Δμ_th` 做标尺；`a≤0`
拒绝。Fe K 短窗口实测 edge_step 相对误差 5.8e-16、平台 1.1353；
传统多项式 2.29% 且平台 NaN。chi_std 把低 R 泄漏 0.316→0.0618；
`ChiK.chi_std_used/clamp_lo/clamp_hi` 记录实际生效状态。

### B10-2 PCA / SPOIL / MCR-ALS
```python
xafs_lcf.run_pca(spectra, xmin, xmax)
# → {"x","components","scores","loadings","explained_variance_ratio",
#    "n_components","n_components_ind"}
xafs_lcf.target_transform(pca, candidates, n_components=None,
                          spoil_good=3.0, spoil_poor=6.0)
xafs_lcf.mcr_als(spectra, n_components, max_iter=200, tol=1e-6,
                 nonneg=True, closure=False)
# → {"components","concentrations","r_factor","iterations","converged"}
```
SPOIL=RET/REP；真实端元≈0.017、随机谱≈17。MCR-ALS 9 条混合谱
（0.5% 噪声、closure=True、1500 上限）621 次收敛，组分相关
>0.9997、浓度 RMSE 0.0069、r_factor=0.01953。

### B11-1 相移 / 异核共占位
```python
artemis.path_total_phase(datfile, k)     # 总相移 φ_tot(k)
athena.ft_to_r_space(..., phase_correct=φ)
larch_wt.continuous_wavelet_transform(..., phase_correct=φ)
artemis.cif_hetero_cocupancy_branches(cif, absorber, rmax=6.0,
                                      delta_z=3, occ_min=0.10)
```
相移校正把表观峰拉回 Reff：Fe–O FT 1.534→1.9942 Å（Reff
1.9808，误差 0.0134）、WT 1.48→2.0 Å（误差 0.019）。
异核分支按 Reff 聚壳；简并度按**不同晶位**数计。Ti(0.7)/Ru(0.3)
共占 1b 位 → n_geom=8、Ti degen=5.6、Ru degen=2.4、occ 和=1，
`single_element_approx=False`。不要用展开条目数（会得 16）。

### B11-2 关联 Debye / BVS
```python
artemis.correlated_debye_sigma2(T, theta_d, reff, m_absorber,
                                m_scatterer, k_debye=None,
                                sigma2_static=0.0)
artemis.bond_valence_sum(paths, r0_table=None, b=0.37,
                         absorber=None)
artemis.fit_r_space(..., bvs_restraint={"target": 3.0,
                                        "sigma": 0.05})
```
`σ²=σ²_static+[1−sin(2k_D R)/(2k_D R)]·σ²_un`；默认
`k_D=0.35 Å⁻¹`（`last_details.k_debye_source` 可核对）。300 K
近壳 0.003406、R=4 远壳 0.010127；1000 K 远壳/Einstein=1.309；
T→0 平台 0.0010828。BVS Fe–O 6×2.0=3.1280；约束
target=3/σ=0.05 → N=4.8567、R=1.9371、bvs=3.00097、
bvs_chi2=3.79e-4、Rf=0.0109，且与手工重算逐位一致。

### B11-3 / B12-1 三空间 / AIC-BIC / share_e0
```python
artemis.fit_r_space(..., fit_space="r"|"q"|"wt")
# FitResult 新增 aic/bic/fit_space
artemis.fit_multi_datasets(datasets, shared_bonds=None,
                           einstein_models=None,
                           share_e0=True|False|"auto")
artemis.profile_chi2_contour(fit_func, x_name, y_name, x_values,
                             y_values, fixed=None,
                             sigma_levels=(2.30, 6.18, 11.83))
```
三空间合成单壳 N/R/σ² 最大差 1.2e-9/1.1e-9/1e-17；默认 `"r"`
与显式 `"r"` 逐位一致。`share_e0="auto"` 按 `edge_group` 分组；
双吸收边真值 +2.0/−6.0 eV → 精确回收、global Rf≈2.1e-29；
`True` 单组 ΔE₀、Rf≈0.0203。BIC 对多参数惩罚更重
（aic_gain > bic_gain；每个数据集 BIC>AIC）。

### B12-2 SAC/DAC / 元素扩充
```python
artemis.sac_cluster_sites(bond_m_ligand=2.0, ligand="N",
                          n_coord=4, axial=None, axial_bond=2.10,
                          carbon_shell=True, ring_radius=2.90)
artemis.diatomic_m1m2_sites("Fe", "Co", bond_length=2.5,
                            ligand="N", n_coord=4, ligand_bond=2.0)
```
SAC：中心吸收原子不在返回列表；M–N₄+C 共 12 位点
（4×N@2.0、8×C@2.9），FEFF 生成 9 路径、第一壳 Reff=2.0。
DAC：M1–M2 2.5 Å + 8×N。`absorber.py` 新增 12 元素全部
`s02=None, calibrated=False`；Fe/Co/Ni/Cu/Pt 标定值逐位未改。

### CLI 补齐（run_xafs / xafs_report 同源）
`--fit-kweights 1,2,3`（B1-1）、`--include-ms`（B1-3）、
`--lcf-standards`（B2-2）、`--deglitch` / `--deadtime`（B2-3）、
`--vary-c3 first|all`（B4-2）、`--wt-sweep`（B3-1）、
`--reject-outlier-scans`（B6-1）、`--isotropic-deltar`（B7-1）、
`--vary-c4 first|all`（B9-1）、`--mback` / `--chi-std` /
`--clamp-lo` / `--clamp-hi`（B10-1）、`--fit-space r|q|wt`
（B11-3）、`--bvs-target` / `--bvs-sigma` / `--bvs-b`（B11-2）。
配置优先级：CLI > 配置文件 > 代码默认值。

### CI（run_ci.py + .github/workflows/xafs-ci.yml）
`python run_ci.py` 执行 py_compile + 60 项契约回归，产物
`ci_report.json` / `ci_report.md`，失败非零退出；`--stress` 追加
Co/Ni/Cu/NiO 压力测试。GitHub Actions 在 push/PR 自动运行。

### 批处理输出目录
`--batch-dir` 优先级：显式 `-o` > 配置 `batch.outdir` > `xafs_out_batch`。
（旧版忽略 `-o`，已修复并实测。）
