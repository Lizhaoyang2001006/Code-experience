# 催化复杂体系非理想结构 XAFS 建模与 CN 拟合执行方案

> **适用范围**：单原子催化剂（SACs）、界面负载催化剂（如 Pt/CeO₂）、晶格掺杂/置换固溶体、非高斯强静态无序结构。  
> **核心逻辑**：架构统一抽象（4 大原型引擎） + 参数按需装配（具体体系输入） + 门禁绝对统一（Stern 独立点数与全局约束）。

---

## 一、 快速执行命令速查表（Copy & Paste）

| 催化体系 | 核心难点 | 终端一键运行命令（Bash） |
|:---|:---|:---|
| **单原子 SAC**<br>($Fe-N_4-O_{\text{axial}}$) | 轻原子相移接近、轴向分子破坏平面称 | ```bash
python run_xafs.py data/Fe_SAC.xdi --absorber Fe --fit \
  --catalytic-model "sac:N4@1.98,axial=O@2.12,C@2.88" \
  --fit-rmin 1.0 --fit-rmax 3.0 -o xafs_out/Fe_SAC
``` |
| **界面负载**<br>($Pt/CeO_2(111)$ 界面) | 表面 Pt–O 与次表层 Pt–Ce 界面散射重叠 | ```bash
python run_xafs.py data/Pt_CeO2.xdi --absorber Pt --fit \
  --phase-cif models/Pt_CeO2_interface.cif \
  --fit-kweights 1,2,3 --fit-rmin 1.2 --fit-rmax 3.6 \
  --constraints "dPtCe = 1.62 * dPtO" -o xafs_out/Pt_CeO2
``` |
| **晶格掺杂**<br>($Fe$ 掺杂 $TiO_2$) | 离子半径失配引起局域畸变、变量易发散 | ```bash
python run_xafs.py data/Fe_TiO2.xdi --absorber Fe --fit \
  --phase-cif models/rutile_TiO2.cif \
  --isotropic-deltar --fit-rmax 3.2 -o xafs_out/Fe_TiO2
``` |
| **强静态无序**<br>(非晶水合氧化物 / 高温原位) | 键长非高斯分布，常规拟合使 N 虚低 20%~40% | ```bash
python run_xafs.py data/amorphous_Co.xdi --absorber Co --fit \
  --fit-motif co_oxide_nano --vary-c3 \
  --fit-rmin 1.0 --fit-rmax 2.2 -o xafs_out/amorphous_Co
``` |

---

## 二、 4 大典型场景深度建模机制与物理约束

### 1. 单原子催化剂（SAC Archetype）

#### 1.1 局域团簇生成逻辑
* **中心吸收原子**：固定在原点 $(0, 0, 0)$。
* **平面第一配位层**：$n$ 配位多边形（如 $N_4$），以半径 $d_{M-L} \approx 1.95 \sim 2.05\text{ \AA}$ 分布在 $xy$ 平面。
* **轴向反应配体**：沿 $z$ 轴正负方向添加（如吸附态 $O$、$OH$ 或外来配体 $Cl$），$z = \pm d_{\text{axial}}$。
* **次近邻载体碳环**：在半径 $R \approx 2.85 \sim 3.00\text{ \AA}$ 处布置 8 个载体基底 $C$ 原子。

#### 1.2 拟合约束卡
* **$\Delta E_0$ 绑定**：平面 $M-N$ 与轴向 $M-O$ 必须绑定全局同一 $\Delta E_0$（1 自由度）。
* **轻元素解耦红线**：严禁无约束自由放开同壳层 $M-C$、$M-N$、$M-O$ 多条独立路径。
* **配位数物理上限**：严格施加边界限制 $N(M-N) \le 4.2$，$N(M-O) \le 2.1$。

---

### 2. 界面负载催化剂（Interface Archetype，以 Pt/CeO₂ 为例）

#### 2.1 局域切片生成逻辑
1. 提取载体晶面超胞（如 $CeO_2(111)$）；
2. 在表面氧空位或三配位台阶位放置吸附原子 $Pt$；
3. 以 $Pt$ 为原点，切取 $R_{\text{cut}} = 5.0\text{ \AA}$ 范围内的所有近邻原子，直接生成非周期性 FEFF 输入文件。

#### 2.2 拟合解耦规则
* **路径解耦判据**：
  * 第一壳 $Pt-O_{\text{surf}}$：键长 $R \approx 2.00\text{ \AA}$；
  * 第二壳 $Pt-Ce_{\text{sub}}$：键长 $R \approx 3.25\text{ \AA}$（重原子散射，高 $k$ 强烈背散射包）；
  * 金属核 $Pt-Pt$：键长 $R \approx 2.75\text{ \AA}$（仅团簇或金属相出现）。
* **几何刚性约束**：利用载体骨架的刚性几何关系，将次表层距离与第一壳联动：
  $$\Delta R(Pt-Ce) = 1.62 \times \Delta R(Pt-O)$$
  不额外消耗独立拟合点数，彻底规避过拟合。
* **小波变换（WT）质检**：确认 $R \approx 3.2\text{ \AA}$ 处的信号极大值落在 $k \approx 8.5\text{ \AA}^{-1}$（$Ce$ 特征峰），而非 $k \approx 10.5\text{ \AA}^{-1}$（$Pt$ 金属峰）。

---

### 3. 晶格掺杂体系（Doping Archetype）

#### 3.1 母体置换与局域应变
* 将母体晶格（如锐钛矿 $TiO_2$）中心阳离子强制置换为掺杂元素（$Fe$）；
* 激活管线 `--isotropic-deltar`，使用单一各向同性应变参数 $\alpha$：
  $$\Delta R_j = \alpha \cdot R_{0,j}$$
  将八面体 6 条不同键长的伸缩收敛为 1 个自由度，配合全局 $\Delta E_0$ 与统一 $\sigma^2$，总拟合变量仅为 3 个，严格满足 $N_{\text{var}} \le \frac{2}{3} N_{\text{idp}}$。

---

### 4. 强静态无序与非晶体系（Disorder Archetype）

#### 4.1 非高斯展开（$C_3$ 偏斜度校正）
* 当体系存在键长不对称分布（软晶格、非晶、高温原位）时，激活三阶累积量：
  $$\chi(k) \propto \sin\left(2kR - \frac{4}{3}C_3 k^3 + \text{phase}\right) \exp(-2\sigma^2 k^2)$$
* **参数界限**：仅放开第一配位壳层的 $C_3$，且必须落入合理物理区间 $|C_3| \le 1.0 \times 10^{-3}\text{ \AA}^3$，挽回被常规高斯模型低估的有效配位数 $N$。

---

## 三、 论文填报模板与合规声明规范

### 1. 表格规范填报示例（Pt/CeO₂ 界面催化剂）

| 壳层 | 散射路径 | 配位数 N | 键长 R (Å) | $\sigma^2$ ($10^{-3}$ Å²) | $\Delta E_0$ (eV) | $R\text{-factor}$ |
|:---|:---|:---|:---|:---|:---|:---|
| 1 | Pt–O (界面表面氧) | 3.8 ± 0.4 | 2.01 ± 0.02 | 4.5 ± 1.1 | +4.8 ± 0.9 | 0.048 |
| 2 | Pt–Ce (载体次表层) | 2.9 ± 0.6 | 3.26 ± 0.03 | 7.8 ± 1.8 | +4.8 ± 0.9 | — |

> **论文脚注标准文案**：  
> 拟合采用基于 $CeO_2(111)$ 界面切片模型生成的局域非周期性团簇；拟合区间为 $k = 3.0 \sim 12.0\text{ \AA}^{-1}$，$R = 1.2 \sim 3.6\text{ \AA}$；独立数据点数 $N_{\text{idp}} = 15.3$，总自由变量数 $N_{\text{var}} = 6$（严格满足 $N_{\text{var}} \le \frac{2}{3}N_{\text{idp}}$）；多壳层共用单一全局 $\Delta E_0$；$S_0^2 = 0.8075$（由 Pt 箔同步标定）；报告的键长 $R$ 已经过 FEFF 理论相移函数修正，未引用偏短的 FT 表观峰位。
