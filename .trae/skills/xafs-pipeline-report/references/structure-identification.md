# 结构识别与 motif 注册表（S6.5）

S6.5 在 S5/S6 之后、CN 拟合之前回答：**"实验 χ(R) 像库里的谁？"**
它只做候选排序与覆盖审计，**不做 CN 定量**——定量仍走逐候选全拟合。

## 1. 完备性哲学（为什么这样设计）

完备性不能靠 agent 的记性，靠三个确定性机制：

1. **生成式枚举**：motif 按核数分类（m1 孤立单点 / m2 双核 / m3-6 寡聚体 /
   np 纳米颗粒 / crystal 晶相）组织在 `_structure_lib/<元素>/`，
   新结构类型出现 = 往库里加一条 JSON，之后每次运行自动包含。
2. **覆盖矩阵**：实验 |χ(R)| 的每个显著峰（≥15% 窗口峰高、
   prominence ≥2%、半高宽 ≥0.08 Å、离窗口边缘 ≥0.1 Å）必须被最优
   候选的理论谱解释（|ΔR|≤0.15 Å），否则列入 `uncovered_peaks`——
   **未解释特征 = 假设空间不完备，应扩库，不得静默略过**。
   prominence/宽度/边缘门槛用于排除低 R 背景鼓包、单点噪声尖刺和
   FT 截断边界效应；真实壳层受 k 窗口限制，FWHM 不会窄到 0.05 Å。
3. **不可分辨并列**：与最优候选 Δscore < 0.05 的候选强制并列报出
   （`indistinguishable`），该分辨窗口内分不开时不装确定。

## 2. motif 注册表（structure_lib.py）

双源统一——晶体类靠"空间群骨架"，非晶类靠"直角坐标"：

```jsonc
// kind=crystal：空间群 + 晶胞 + 非对称单元（Wyckoff 特征点）
{
  "id": "hematite", "absorber": "Fe", "kind": "crystal",
  "nuclearity": "crystal",
  "cell": {"a": 5.0356, "c": 13.7489, "gamma": 120},
  "spacegroup": "R-3c",
  "asym_sites": [{"el": "Fe", "x": 0.0, "y": 0.0, "z": 0.3553},
                 {"el": "O", "x": 0.3059, "y": 0.0, "z": 0.25}],
  "source": {"type": "cod", "ref": "COD-xxxxxxx"}
}
// kind=cluster：吸收原子在原点（不入列），sites 只放近邻
{"id": "fen4_planar", "kind": "cluster", "nuclearity": "m1",
 "sites": [[1.93,0,0,"N"], ...]}
```

注册流程（`structure_lib.register`，一次性）：
骨架 → gemmi 展开（与 `artemis.read_cif_sites` 同一引擎，
`_structure_lib/<el>/_cif/<id>.cif` 留档）→ 近邻坐标表 →
几何壳层指纹 `shell_expect`（ΔR<0.08 Å 并层，自动计算）→
查重（指纹容差 0.05 Å）→ JSON 落盘。

入库通道：`structure_lib.import_cif(cif_path, ...)` 从任意来源的小 CIF
抽取骨架。**在线检索（COD REST API `crystallography.net/cod/result?...&format=cif`）
只是添条目的手段，失败不影响库内既有条目**——库的主体是人工精管的
种子条目，这是规避网页抓取脆弱性的架构决策。

## 3. S6.5 识别方法（structure_id.py）

```
每个 motif → FEFF 理论 χ(k)（自然简并度、s02、nleg≤4 含 MS）
→ σ² 网格（0.003/0.006/0.009 Å²）× ΔE0 网格（-8..+8 eV，
   k 轴重采样：k'² = k² + 0.2625·ΔE0）
→ 与实验同窗口 xftf
→ score = min_α ‖χe(R)−α·χt(R)‖²/‖χe(R)‖² (R∈[1,4] Å，复数)
```

- **复数谱比较，不是只看 |χ(R)|**：同时约束幅度与相位；允许一个
  复数尺度 α 吸收整体振幅/全局相位，但不允许逐点调相。
- **判定四态**：`confident`（best≤0.16、margin≥0.03，且覆盖清单为空）/
  `partial`（分数/margin 通过但有未解释峰，模型不完整）/
  `ambiguous`（与次优过近）/ `no_match`（库内无相似候选）。
  非 confident 默认只做 L1 第一壳，`--allow-ambiguous` 才强制 Top-1。
  实测：Hedenbergite/Vivianite 等非库矿物会被 `partial`/`ambiguous`
  拦下；磁铁矿在只有 16d 八面体位时被标 `partial`（正确物相、不完整
  模型）——**已由 `magnetite_tet`（8a 四面体位）+ `magnetite_mixed`
  联合条目解决，实测 confident（score 0.083）**。
- **larch 教训**：裸 feffpath 的 `e0` 属性实测不被 `_calc_chi()` 消费
  （χ 完全不变）——所以 ΔE0 用在求和后 χ(k) 上重采样实现，
  不依赖 larch 内部行为。改库前若有疑问先跑小诊断脚本验证。
- FEFF 缓存键 = motif id + 近邻坐标哈希前 8 位；非 K 边再加边名
  （`_feff_cache/lib_<id>_<hash>_l3/`），坐标/吸收边变缓存自动失效。

## 4. 已知边界（诚实声明）

- 指纹是**候选排序**不是结构确证；弱结晶/非晶（水铁矿）无可靠 motif，
  预期落空 → 只报第一壳（这是科学诚实的选择）
- `uncovered_peaks` 是**保守**机制：窗口干涉侧瓣仍可能触发误报；
  已用半高宽 ≥0.08 Å 过滤单点噪声尖刺，但侧瓣是有限 k 窗口的物理
  产物，触发后人工看一眼，不是必然错误
- 多位点吸收已支持：`Motif(kind="mixture", components=[{id,weight},...])`
  ——χ = Σ 权重·χ_位点（识别与拟合都支持；先注册各位点再注册联合条目，
  权重和=1、不嵌套）。覆盖判定对 mixture 按"复合谱 ∪ 各分量分别解释"：
  复合谱中分量峰受干涉位移/压制（实测磁铁矿 tet 位 Fe-Fe@3.48 峰在
  复合谱中相对峰高 14.7% 被过滤），不分分判定会把正确模型错标 partial。
  拟合侧 `_mixture_rows` 按权重缩放简并度、并把合并容差收紧到 ≤0.10 Å
  （磁铁矿 Fe-O 1.886/2.060 Å 只差 0.17 Å，0.20 会错并）。
- mixture 相分数（`--fit-weights`，2026-11 落地）：识别用固定晶体学
  权重，定量时保留 `N_crystal`、用 softmax 放开组分权重（和=1），
  输出 `mixture_weights`。合成磁铁矿 tet/oct 第一壳 0.25/0.75 可回收
  至 ≤0.03；真实磁铁矿默认窗 Rf=0.097 → `weights_reliable=False`，
  因此不能把两态权重当作独立物相定量，且 N 列是 N_crystal·w 不是配位数。
- cluster motif 若含共线构型（平面 FeN₄），MS 显著——已放开 legs≤4；
  但构造几何的 motif 只是假设，建议以 DFT 精化坐标替换后重注册
- N≈O（ΔZ=1）等背散射体不可分对，指纹识别不了，靠 XANES/CWT 联判
- **简并簇不可辨识**（针铁矿 goethite 实测，2026-11 五轮诊断）：S6.5 能
  正确 confident 命中（score 0.156, margin 0.085），但 S7 的 Rf≈0.17 与
  `cn_reliable=False` 是**数据素质 + 模型歧义的正确反映**，不是拟合 bug：
  - 3.0–3.44 Å 的 Fe/O 简并簇（Fe@3.02×2、O@3.21×1、Fe@3.31×2、
    Fe@3.44×4）互相干涉，选壳敏感：3 壳 O321→Rf 0.19；换 Fe331→Rf 0.11
    但 N 全撞上界 6.0（简并吞噬=失去定量意义，**拒绝**）；4 壳 Rf 0.117
    但 N_err ±6.1/±35/±16.5 爆炸（过参数化）。
  - 中距离壳 σ² 被压到 0.0001–0.0007 Å²（阈值 ≤5e-4 触发
    `sigma2_at_bound`）= 模型用**压缩 σ²** 补偿形状失配的强信号
    （不是"振动极小"）。该旗作为 `cn_reliable=False` 的独立原因写进
    S7 warning；goethite 实测触发壳层 Fe-O@3.21（σ²=0.0001）。
  - ΔE₀ 自由/固定对照（−13 vs 固定后 Rf 0.232）证明实验-FEFF 存在真实
    ~13 eV 相位偏移（数据边跳仅 0.162，低信噪放大归一化/E₀ 偏差）。
  - 第一壳可靠：合并单 O 壳 N=5.89±0.51（晶体学 6，偏差 −2%），R=1.952。
    但把它拆成 O1/O2 双峰虽把短窗 Rf 压到 0.077，N_err 却涨到 ±6.35
    （≈N 本身，同元素简并吞噬）——**短窗 Rf 优不等于 CN 可信**。
  - 窗口下沿同样不可自由调：rmin 1.0→1.15→1.25 时 Rf 0.117→0.079→0.067，
    但 N 6.0→6.9→7.8 膨胀——"压 Rf"与"CN 可信"在此数据上数学冲突，
    没有客观选窗准则能同时满足，故**不采纳**任何压 Rf 的窗口/换壳方案。
  结论：报告只给第一壳 CN 锚点，中远壳声明不可定量；不得为降低 Rf
  改选散射体或窗口。

## 5. 扩展流程

**新 motif**（3 步）：
```python
from structure_lib import Motif, register, import_cif
m = import_cif("xxx.cif", "new_phase", "Fe")   # 或 Motif(kind="cluster", ...)
register(m)                                     # 展开+指纹+查重+落盘
```

**新元素**（3 步，代码零改动）：
1. `absorber.ABSORBERS` 加条目（s02=None 起步）
2. 测标准样标定 S₀² 填入（未标定时 structure_id 拒绝运行）
3. `_structure_lib/<元素>/` 挂 motif

**种子库**（当前 11 条，`python _structure_lib/build_seed.py` 重建）：
lepidocrocite / hematite / magnetite_oct / **magnetite_tet** /
**magnetite_mixed**（mixture, ⅓χ_8a + ⅔χ_16d）/ wustite / goethite / fe_bcc
（晶体）+ fen4_planar / fe2_mu_oxo / fe2n6_dual（团簇，构造几何）。
真实 BMM Fe 标准验证：hematite / goethite / wustite / magnetite（mixed 模型）
为 confident，ferrihydrite 为 no_match（弱结晶，无可靠体相 motif）。
另有压力测试顺产的参考相种子（不在 build_seed 内）：
`Co/co_hcp`（stress_test_co.py）、`Ni/ni_fcc`（stress_test_ni.py）、
`Cu/cu_fcc`（stress_test_cu.py）、`Ni/nio_rocksalt`（stress_test_nio.py，
多元素/多壳层/强 MS 的复杂固体代表）——幂等可重跑，随对应压力测试维护。
**Pt/pt_fcc**（Fm-3m a=3.9231）由 `build_seed.py` 注册，是 Pt L3
S₀² 标定源；目前 Pt 库只有这一条，分子配合物会正确 no_match → L1。
单候选时 `confident` 只表示"没有竞争假设"（margin=n/a），不是强证据；
Pt 库扩到多条结构后判别力才有意义。

## 6. 特殊体系工具（structure_lib）

**掺杂替位**（ZSM-5 类场景，候选名册机械生成）：
```python
from structure_lib import substitute_series, register_series
cands = substitute_series(host, "Si", "Fe", coord="tet")   # 每个 T 位一条
made, skipped = register_series(cands)   # 等价位点指纹重复自动跳过
```
弛豫规则：第一壳整体径向缩放到 `BOND_REF` 参考键长（**保留层内劈裂**），
远程保持宿主几何；`BOND_REF` 查不到的离子对必须显式给 `relax_to`，不猜。
深度弛豫属 DFT 范畴——算完拿坐标重注册即可。

**纳米颗粒**（表面-核心平均，N 先验）：
```python
from structure_lib import nanoparticle_motif
m = nanoparticle_motif(fe_bcc, 2.0)   # 2 nm 球形：size_model 含各壳层 scale
```
只改 N 先验不改散射物理（FEFF 路径仍用体相）；硬球截断、
不含表面重构；<1 nm 统计噪声大。

**部分占据**：`read_cif_sites(..., with_occupancy=True)` →
指纹 `shell_expect` 带 `n_occ`（占据率加权计数，与几何计数 n 并存）。

**分布偏斜**：`fit_r_space(..., vary_c3=True)` 放开三阶累积量 C₃
（larch 的 third 参数）。仅用于纳米颗粒表面张力/强非谐；
C₃ 与 ΔR 部分相关，常规体系不要开。

## 7. L1 降级（识别非 confident 时）

`identify_result` 的 `status` 不是 `confident`（`partial`、`ambiguous`
或 `no_match`）时，自动追加 `l1_first_shell`：合成八面体/四面体团簇的
SS 路径拟 N/R/σ²/ΔE0，并标注 **"仅第一壳可信，长程结构未知"**。
这是弱结晶/非晶/库外体系唯一诚实的定量路径——长程拟合是强行给
数据安结构。

```python
from structure_id import l1_first_shell_fit
l1 = l1_first_shell_fit(res, ligand="O", geometry="oct")  # 或 tet/S/Cl
# 自测基准：lepidocrocite → N=5.84±0.39, R=1.948 Å, Rf=0.073
```

`ligand`/`geometry` 是显式结构假设，默认 O/oct 只适用于多数氧化物；
硫化物/氯化物/四面体体系必须由调用方覆盖（CLI: `--l1-ligand`,
`--l1-geometry`）。第一壳若被 FEFF 拆成多条 SS 路径，N 取壳层总和。

**选峰与窗口（2026-11 实测修正）**：r0 = 配体感知选出的第一壳峰 +0.45 Å
——用共价半径和（`artemis.COVALENT_RADII`）估预期键长，在显著峰里选
最近者。不能取 argmax（Fe 氧化物会选到更强的 Fe–Fe 第二壳 2.55），也
不能取"最低显著峰"（Pt L3 白线在 1.2–1.4 Å 的伪峰会被选走）。R 窗上限
`min(3.0, max(2.3, peak+0.15))`，否则 Pt–Pt 2.77 / Fe–Fe 2.48 的模型峰
会被截在窗外（曾出现 Rf=1.0 的伪失败）。实测：K₂PtCl₄ `--l1-ligand Cl`
→ N=3.98±0.14、R=2.302 Å、Rf=0.023；四氨合铂 `--l1-ligand N` → Rf=0.024。
默认 O 对 Pt–Cl 会 Rf=0.32 并正确判 `cn_reliable=False`——配体假设错了
必须显式换，管线不会替你猜。
