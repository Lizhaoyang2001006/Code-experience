# AGENTS.md — XAFS 项目 agent 入口

本文件面向新加入的 AI agent。开始任何 XAFS 处理前，先确认工作目录是项目根目录，并阅读：

1. `README.md`（能力、CLI、数值基准）
2. `.trae/skills/xafs-pipeline-report/SKILL.md`（执行流程与交互规范）
3. `.trae/skills/xafs-pipeline-report/references/pipeline-reference.md`（参数与门禁）

## 工作方式

- 默认先做 S1–S6：读谱、E0、归一化、背景扣除、FT、结构识别。
- 先把 E0、边步、Rbkg、第一壳位置、结构候选给用户确认，再做 S7 定量拟合。
- 吸收元素或边不确定时先问用户，不要猜。
- 用户给原始谱文件时，先确认文件存在、格式（XDI / Athena `.prj` / 普通 dat）和吸收元素。
- 输出统一写到 `xafs_out/<样品名>/`，不要覆盖已有结果。
- 修改代码、提交、推送前先征得用户同意；不要擅自 `git clean`。

## CLI 入口

单文件报告：

```powershell
python xafs_report.py <文件> --absorber <元素> -o xafs_out/<样品名>
```

批量 / 结构识别 / 拟合：

```powershell
python run_xafs.py <文件> --absorber <元素> -o xafs_out/<样品名>
python run_xafs.py <文件> --absorber <元素> --fit -o xafs_out/<样品名>
python run_xafs.py <文件> --absorber <元素> --fit --fit-motif <motif> --allow-ambiguous -o xafs_out/<样品名>
```

多元素批处理：

```powershell
python run_xafs.py a.xdi b.xdi --multi-element -o xafs_out/batch
```

其他常用参数：`--fit-weights`（混合相权重）、`--fit-rmax`、`--e0-method half_height`、`--xanes-cross-check`、`--diagnostics`、`--figure-pack`、`--operando`。参数以 `--help` 为准。

## S₀² 规则（重要）

优先级：显式 `--s02` / `--multi-element-s02` > `--chemical-state` > `--s02-from-t5` > `absorber.py` 静态表。

- Fe 氧化物 / 羟基氧化物：默认 `Fe=0.803`，不要替换。
- 金属 Fe：`--chemical-state metal`（0.6566；等价旧写法 `--s02 0.6566`），并配 `--fit-rmax 2.60` 做单壳 bcc。
- T5 已标定元素：`--s02-from-t5`（只接受 `status=calibrated` 且 `static_absorber_py_unchanged=true` 的产物）。
- 其他未标定元素：显式 `--s02`，或先跑 `calibrate_t5.py` / `calibrate_s02.py`。

## 硬约束

- 不修改 `absorber.py` 的静态 S₀² 表；12 个锁定元素保持 `s02=None`、`calibrated=False`。
- 不放宽 `calibrate_s02.py` 的冻结门禁常量。
- 不把 `_standards_lib/`、`data/standards_raw/`、`logs/`、`patches/` 加入 Git（IIT-Lytle 再分发许可 / 本地产物）。
- `data/s02_calibration/` 是已通过产物；重新标定用 `calibrate_t5.py --only ...`，不要手改 JSON。
- 拟合失败或门禁不通过时，如实报告并降级，不要通过放宽门禁“修好”。

## 常用自检

```powershell
python env_check.py
python run_ci.py            # 编译 + 135 项回归
python cif_crawler.py --t5  # 重建本地 T5 标准谱（产物不入库）
```
