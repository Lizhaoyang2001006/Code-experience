# XAFS 管线项目贡献指南 (Contributing Guide)

感谢你关注并愿意为本 XAFS 光谱处理与定量分析管线贡献代码！

本管线是一套面向高水平学术研究与发表级可复现性的严谨工具。为了维护项目的物理真实性与零退化基线，请在提交代码前仔细阅读本指南。

---

## 一、 核心科学哲学与不可动红线

在编写任何代码前，请牢记本项目的**三项最高原则**：

1. **物理真实性高于拟合残差**：严禁为了让拟合指标（如 $R\text{-factor}$）数值好看而人为放宽物理门禁，或盲目引入过多自由变量引发过拟合。
2. **遇到异常如实降级**：当拟合参数不满足物理合理性或超出置信区间时，代码逻辑必须如实降级（退守 L1 第一配位层或标记为不可信），严禁掩盖拟合失败。
3. **严格遵守 `AGENTS.md` 架构宪法**：
   * **12 锁定元素不可动**：`absorber.py` 中的 `Cr, Mn, Zn, Mo, Ru, Rh, Pd, Ag, Ir, Au, Ce, U` 必须严格保持 `s02=None, calibrated=False`，未经过严格金属箔实验标定前严禁擅填默认值；Fe 氧化物生产值锁定为 `0.803`；
   * **独立点数硬门禁（Stern 规则）**：严格保证自由变量数 $N_{\text{var}} \le \frac{2}{3} N_{\text{idp}}$，超出预算的高阶壳层必须自动截断；
   * **多壳层共用 $\Delta E_0$**：所有路径必须绑定全局单一 $\Delta E_0$ 自由度，严禁各壳层独立浮动；
   * **荧光自吸收拒收红线**：浓样 $\Delta\mu \ge 1.0$ 未校正自吸收时，强制触发 `cn_reliable=False`。

---

## 二、 本地开发环境搭建

### 1. 基础依赖安装
推荐在 Python 3.10 或更高版本的虚拟环境中开发：

```bash
# 克隆仓库
git clone https://github.com/your-org/xafs_pipeline.git
cd xafs_pipeline

# 创建并激活虚拟环境
python3 -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

# 安装开发与运行依赖
pip install -r requirements.txt
```

### 2. 核心科学计算依赖
管线核心数值计算依赖 `numpy` 和 `scipy`。若需运行全流程预处理与 FEFF 路径求解，推荐配置：
* `larch`（XAFS 背景扣除与傅里叶变换算法栈）；
* `FEFF 8.4 / 9.0`（用于散射路径理论计算）。

---

## 三、 代码扩展规范与开发指引

### 1. 添加新的催化 Motif 或扩展算法建模（`catalytic_builder.py`）
* **原点对齐约束**：所有局域非周期性点阵构建，吸收原子必须严格放置在笛卡尔坐标原点 $(0, 0, 0)$；
* **数值可复现性**：任何涉及随机高斯无序（`--disorder-sigma`）的算法，必须固定随机数种子（默认 `np.random.seed(42)`），确保科学计算与 CI 测试结果严格可复现；
* **缓存隔离**：动态生成的催化构型与 FEFF 计算必须隔离在 `_dynamic_cache/<hash>/`，严禁污染 `_structure_lib/` 主库静态指纹索引。

### 2. 报表与绘图输出（`xafs_report.py`）
* **键长防伪分离**：报表（CSV、PNG、HTML）输出的键长 $R$ 必须且只能是经过 FEFF 理论相移函数修正后的拟合值；
* **表观峰位标识**：五联图 $|\chi(R)|$ 上的峰位未经相移修正（系统偏短 0.3~0.5 Å），必须明确标注为“表观距离（Apparent R）”，严禁误导用户直接填入学术论文表格。

---

## 四、 提交流程与 CI 回归核验

在发起 Pull Request (PR) 之前，**必须在本地终端全绿通过以下自检套件**：

```bash
# 1. 运行核心单元测试与 CI 校验
python run_ci.py

# 2. 运行基线回归测试（保证标样 Lepidocrocite 与 Fe foil 数值逐位一致）
python regression_test.py
```

### Pull Request 提交清单 (Checklist)
提交 PR 时，请确保在描述中确认以下事项：
- [ ] 代码修改未放宽任何 $R\text{-factor}$、$\Delta E_0$、结构识别评分或自吸收门禁阈值；
- [ ] 未擅自修改 `absorber.py` 中的 12 个锁定元素定义；
- [ ] 本地 `python run_ci.py` 与 `python regression_test.py` 100% 通过；
- [ ] 新增功能已补充对应的单元测试或实例验证脚本；
- [ ] 对应更新了相关的文档（普通用户使用更新 `user_guid.md`，架构规范更新 `AGENTS.md`）。

---

## 五、 问题反馈与交流 (Issues)

如果你在处理实际催化数据时发现异常、算法缺陷或有新的非理想构型原型需求，欢迎在 GitHub Issues 中提交反馈。提交 Issue 时请提供：
1. 吸收元素与测试模式（透射 / 荧光）；
2. 运行的完整 CLI 命令与终端报错日志；
3. 输出的 `*_indicators.csv` 与 `*_meta.json` 关键门禁记录。
