# Figure Spec — 五联图与 CN 表逐项规格

复刻已验收的图面（参考版式：三扫描 μ(E) + XANES + k²χ(k) + |χ(R)| + WT）。

## 全局 rcParams（出图脚本开头必设）

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "axes.linewidth": 1.0, "font.size": 11,
})
```

## 五联图布局

```python
fig = plt.figure(figsize=(14.5, 6.4))
gs = GridSpec(2, 3, figure=fig, hspace=0.46, wspace=0.33)
ax_a, ax_b, ax_c = gs[0, 0], gs[0, 1], gs[0, 2]
ax_d, ax_e = gs[1, 0], gs[1, 1:]   # e 跨两列
```

## 逐面板规格

### (a) μ(E) 多扫描
- **xlim 以 E0 为基准自动缩放**：`max(E_min, E0−150)` 到
  `min(E_max, E0+250)`——不再硬编码 Fe 区，Co/Ni/Cu/Pt 同样可用。
  绝不画全数据范围（全范围图与参考图"趋势走向"看起来不同）。
- 配色按 `SCAN_COLORS` 循环（深→浅：`#a50026`, `#f46d43`, `#fdae61`）；
  lw=0.9；label "scan 1..N"，legend frameon=False, loc="lower right",
  fontsize=9。
- E₀ 竖直蓝虚线（ls="--", color="b", lw=1.0）+ 两行标注
  `f"{edge_label}\n$E_0$ = {e0:.1f} eV"`（edge_label="Fe K-edge"，
  蓝色，放线右侧、y 底部上方 5% 处；**元素+边必须写明，防读图误判元素**）。
- xlabel "Energy (eV)"；ylabel r"$\mu(E)$"。

### (b) 归一化 XANES
- 合并谱蓝色 lw=1.1，label "merge"；legend loc="upper right"。
- **白线标注**：在 E0±(10,30) 窗口内找 μ_norm 最大值 (wl_e, wl_h)，
  画空心圆标记 `"o", ms=5, mfc="none", mec="b", mew=1.2"`，
  annotate `f"white line\n({wl_h:.2f})"`（蓝色, fontsize 8.5, 文字在峰右上）。
- xlim `max(E_min, E0−40)` 到 `min(E_max, E0+300)`；
  ylabel "Normalized $\\mu(E)$"。

### (c) k²χ(k)
- 黑色 lw=0.8；xlim `max(0, kmin−1)` 到 `min(k_data_max, kmax+1.5)`。
- xlabel r"$k$ ($\AA^{-1}$)"；ylabel r"$k^2\chi(k)$ ($\AA^{-3}$)"。

### (d) |χ(R)|
- 黑色 lw=1.0；xlim `0` 到 `max(3.5, 最远壳层 hi)+0.3`；xlabel r"$R$ ($\AA$)"；
  ylabel r"$|\chi(R)|$ ($\AA^{-3}$)"；ylim 取 `1.35×峰高` 给标签留空间。
- 标题 "Fourier transform (not phase-corrected)"，fontsize=10
  （**必须保留 not phase-corrected**——未修正 R 比真实键长短 ~0.4 Å）。
- **壳层归属标注**：对每个壳层窗口找 FT 峰 (rp, hp)，标签放在坐标轴
  顶部左/右两侧（axes fraction 0.04/0.96，交替上下），用引线指向峰顶；
  不再用“峰右侧 +0.18 Å”，否则相邻壳层标签会重叠。多个拟合壳层落到
  同一个未相位修正峰（Δrp<0.10 Å）时合并标签（如
  “Fe-Fe(1) / Fe-O(2)”），WT 面板同样只画一条线。

### (e) 小波 |W(R,k)|
- `contourf(k_centers, r, magnitude.T, levels=21, cmap="jet", vmin=0, vmax=1)`。
- colorbar：`fig.colorbar(im, ax=ax_e, pad=0.02)`，
  label r"$|W(R,k)|$ (normalized)"。
- xlim 与 (c) 相同并钳制到 `k_centers` 范围；ylim 到
  `max(4, 最远壳层 hi)+0.6` 并钳制到 `r` 范围。
- **标题**：`f"{sample}, {edge_label}"`（如 "γ-FeOOH (lepidocrocite), Fe
  K-edge"），fontsize 11——样品+元素写在标题，杜绝旧版角标被误读
  （"Fe_lepidocrocite"曾被读成 "P not included"）。
- **壳层键长标注（硬性要求）**：对每个壳层窗口
  `k_p, r_p, _ = larch_wt.locate_shell(wt, lo, hi)`，然后
  ① `axhline(r_p, color="k", ls="--", lw=0.8)`；
  ② 峰位空心圆 `"o", ms=6, mfc="none", mec="k", mew=1.3"`；
  ③ 文字 `f"{label}  R≈{r_p:.2f} Å"`（transform=get_yaxis_transform()，
  x=0.02、y=r_p+0.07，白底半透明 bbox 消读图歧义）。
  WT 峰位与 FT 峰位差 < 0.33 Å（0.5/σ）属正常，两者并列出现不是矛盾。
- 右下角注 "R not phase-corrected"（fontsize 8，白底半透明 bbox）。

### 面板标签
```python
for ax, lab in ((ax_a,"a"),(ax_b,"b"),(ax_c,"c"),(ax_d,"d"),(ax_e,"e")):
    ax.text(-0.13, 1.06, f"({lab})", transform=ax.transAxes,
            fontsize=13, fontweight="bold")
```

保存：`fig.savefig(out/"..._xas.png", dpi=300, bbox_inches="tight")`。

## CN 拟合参数表（Table S 风格）

```python
fig2, ax = plt.subplots(figsize=(10.5, 3.4))
ax.axis("off")
ax.set_title("Table S1. Fe K-edge EXAFS curve-fitting parameters for ...",
             fontsize=12, pad=18)
```

- 列：Sample / Shell / CN$^a$ / R$^b$ ($\AA$) / $\sigma^2$$^c$ ($\AA^2$) /
  $\Delta E_0$$^d$ (eV) / R-factor。
- colWidths `[0.16, 0.11, 0.13, 0.16, 0.14, 0.15, 0.11]`；
  `tbl.set_fontsize(10); tbl.scale(1, 1.9)`。
- 数值格式：CN `6.5±0.4`（固定 N 写 `6$^f$`）；R `1.994±0.004`
  （R=eff_r，含 ΔR）；σ² 4 位小数；ΔE₀/R-factor 仅组首行。
- 边框：默认虚线格线 `:`；表头行与每组首行实线 `#444444`；
  表头灰底 `#f0f0f0` 加粗。
- 同一样品多壳层：Sample 列仅首行写名，其余空串（matplotlib table
  无合并单元格，靠视觉分组）。
- 脚注（`fig2.text(0.02, 0.02, foot, fontsize=8.5, va="bottom")`）三行：
  1. 上标定义：CN/R/σ²/ΔE₀ + R-factor 是 goodness of fit。
  2. S₀² 值与来源、FEFF 结构来源、有效壳层合并说明
     （如 "Fe–O is an effective single shell (1.98 Å ×2 + 2.07 Å ×4)"）。
  3. 拟合窗口：k 范围、kweight、各组 R 范围。

保存：`dpi=300, bbox_inches="tight"`。

## 数据导出规格（与图一一对应）

全部编码 utf-8-sig（Excel 直开不乱码）；数值 `f"{v:.6g}"`。

| 文件 | 格式 | 列/结构 |
|------|------|---------|
| `<stem>_mu.csv` | 逐点 | energy_eV, mu_scan1..N, mu_merged（各扫描插值到合并网格） |
| `<stem>_xanes.csv` | 逐点 | energy_eV, mu, pre_edge_line, mu_norm |
| `<stem>_chik.csv` | 逐点 | k_inv_A, chi, k2chi |
| `<stem>_chir.csv` | 逐点 | R_A, chi_R_re, chi_R_im, magnitude |
| `<stem>_wavelet.csv` | 矩阵 | 首行 `k\R` + R 值；每行 k 值 + 该 k 下各 R 的 W（4 位小数，与 `WaveletMap.magnitude` 同向） |
| `<stem>_fit_table.csv` | 表格 | 与 PNG 同构 7 列（粘贴友好） |
| `<stem>_fit_table.txt` | TSV | 同上，制表符分隔——粘贴 Word/Excel 自动成表 |
| `<stem>_indicators.csv` | 汇总 | indicator, value, unit, 说明 |

**indicators.csv 必含行**（防读图猜测的最低集）：
absorbing_edge（如 "Fe K-edge"）、e0、edge_step、white_line、
post_edge_plateau、rbkg、k_window、ft_first_peak、
phase_correction（"not applied"）、`wt_<壳层>_peak` 每壳一行、
s02、delta_e0、r_factor、fit_points_vs_vars、fit_<壳层> 每壳一行、
structure_status/coverage_ok（S6.5 判定与覆盖审计）、
`cn_systematic_budget`（CN 绝对系统误差预算，不能只给形式误差）。

## 输出位置与命名

`xafs_out\<stem>\`（默认；`-o/--outdir` 可改），一次共 11 个文件
+ 可选诊断图 6 张（`diagnostics\`）：
`_xas.png`、`_fit_table.png/.csv/.txt`、`_indicators.csv`、
`_mu.csv`、`_xanes.csv`、`_chik.csv`、`_chir.csv`、`_wavelet.csv`、
`_structure_id.json`、`_fit.json`。
入口：`python xafs_report.py <文件...> --absorber Fe`（通用）；
无文件参数或 `--demo` 运行内置 lepidocrocite + Fe 箔示例。
FEFF 结果缓存在 `_feff_cache\`，重复运行不重算。

## 出图前数值检查

- (a) E₀ 标注 = determine_e0 实测值（不是 7112 名义值）。
- (b) 白线峰高 ≈1.5、边后平台 ≈1.0（归一化正确的外证）。
- (d) 第一峰 1.47 Å 量级（γ-FeOOH）；若 ~1.9 Å 说明误用了相位修正网格。
- (e) 两带结构：低 k(5–9) 的 Fe–O 带 + 高 k(8–12) 的 Fe–Fe 带。
