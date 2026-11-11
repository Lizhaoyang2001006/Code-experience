# XAFS 光谱处理与定量分析管线 (XAFS Processing & Fitting Pipeline)

> **核心承诺**：基于严谨物理门禁的 XAFS/EXAFS 自动化分析管线。**不放宽门禁、不产生虚假拟合优度、遇到异常如实降级**，全面支持常规体相材料及前沿催化非理想结构（单原子 SACs、界面负载、纳米合金、掺杂固溶体），确保学术成果具备最高等级的可复现性。

---

## 🚀 极速体验（1 分钟上手）

```bash
# 安装依赖
pip install -r requirements.txt

# 一行命令跑完完整管线（数据预处理 + 结构识别 + 定量 CN 拟合）
python run_xafs.py data/sample_Fe.xdi --absorber Fe --fit -o xafs_out/sample_Fe
```

所有产物统一落在 `xafs_out/sample_Fe/`（不覆盖已有结果），核心重点查看 3 个文件：
* `*_xas.png`：五联图（$\mu(E)$、XANES 边形、$k^2\chi(k)$、未修正表观 FT $|\chi(R)|$、小波变换 WT）；
* `*_fit_table.csv` / `.png`：定量拟合结果表（经 FEFF 相位修正的真实键长 $R$、配位数 $N$、$\sigma^2$、误差棒）；
* `*_indicators.csv`：关键物理门禁核验表（$R\text{-factor}$、$\Delta E_0$、$S_0^2$ 来源、自吸收状态）。

---

## 📚 文档全景导航（按角色查阅）

为避免信息冗余与口径脱节，本仓库建立了分工明确的文档体系：

| 读者角色 | 推荐必读文档 | 包含核心内容 |
|:---|:---|:---|
| 🔬 **实验与科研人员** | [**`user_guid.md`**](./user_guid.md) | **用户操作总指南**：3 分钟跑通、结果红黄绿灯三步判读法、避坑决策树、标准样标定及论文表格规范填报。 |
| 🧪 **复杂催化体系进阶** | [**`catalytic_motif_pipeline.md`**](./artifacts/file_generation/ttl=63d/output/catalytic_motif_pipeline.md) | **非理想结构建模执行方案**：单原子（SACs）、界面负载（如 $Pt/CeO_2$）、晶格掺杂及非高斯静态无序的局域点阵建模与物理约束。 |
| 🤖 **AI 模型与核心开发者** | [**`AGENTS.md`**](./AGENTS.md) | **代码架构宪法**：12 锁定元素不可动、门禁常量不可放宽、Stern 独立点数硬截断、CI 零退化回归红线。 |

---

## 🌟 核心特性与设计原则

1. **绝对守恒的物理门禁（Anti-Overfitting）**：
   * 独立点数预算门禁：严格执行 Stern 规则 $N_{\text{idp}} \approx \frac{2\Delta k\Delta R}{\pi} + 2$，保证自由变量 $N_{\text{var}} \le \frac{2}{3}N_{\text{idp}}$，超出壳层强制截断，绝不欠约束过拟合；
   * 全局 $\Delta E_0$ 绑定：多壳层拟合严格共用单一 $\Delta E_0$ 自由度，杜绝壳层独立浮动导致的假键长；
   * 荧光自吸收硬红线：浓样 $\Delta\mu \ge 1.0$ 未校正自吸收时，强制标记 `cn_reliable=False` 并拒绝定量宣称。
2. **催化前沿非周期性体系原生支持**：
   * 内建动态拓扑生成器，支持 `--catalytic-model "sac:N4@2.00,axial=O@2.10,C@2.90"` 直接构建单原子局域点阵；
   * 原生支持金属纳米颗粒粒径从第一壳配位数降幅自动反推（`nanoparticle_size`）；
   * 混合相多位点采用 Softmax 权重参数化拟合（`--fit-weights`），保证相分数恒正且和为 1。
3. **真实键长与表观峰位防伪分离**：
   * 严格执行 FEFF 相位修正，论文表格严格采用修正后的真实 $R$，明确隔离未校正相移的 FT 表观距离（系统偏短 0.3~0.5 Å）。

---

## 🛠️ 项目结构与模块划分

```
xafs_pipeline/
├── run_xafs.py               # CLI 统一主入口（参数解析、管线调度、落盘）
├── xafs_pipeline.py          # S1-S6 预处理主干（读数、定边、归一化、AUTOBK扣背景、FT、WT）
├── structure_id.py           # S6.5 结构指纹识别与置信分级（confident/partial/ambiguous/no_match）
├── catalytic_builder.py      # v2.0 晶体库驱动与算术畸变催化引擎（免 DFT：切片锚定/掺杂/高斯无序）
├── structure_lib.py          # 结构库管理与动态催化 Motif 生成（原生支持 algo:/sac:/dac:/alloy:）
├── xafs_fit.py               # S7 Artemis 空间拟合编排与 CN 求解（多位点、约束、截断）
├── artemis.py                # 底层空间拟合数学引擎（FEFF 路径求和、雅可比矩阵、协方差、相关性）
├── athena.py                 # Athena 算法底层（能量对齐、去毛刺、样条寻优、自吸收校正）
├── larch_wt.py               # Morlet 高分辨率连续小波变换（轻重散射原子定性辨析）
├── xafs_report.py            # 多面板图像渲染、报表生成与指标导出
├── user_guid.md              # 面向科研人员的完整操作指南
└── AGENTS.md                 # 面向开发与 AI 模型的硬约束规范
```

---

## 🧪 自动化测试与质量保障

在修改任何代码后，必须在终端执行以下检查确保零退化：

```bash
# 运行完整 CI 测试套件
python run_ci.py

# 运行基线回归测试（确保标准参比数据逐位一致）
python regression_test.py
```
