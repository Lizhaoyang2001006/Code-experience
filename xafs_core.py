# -*- coding: utf-8 -*-
"""
xafs_core.py — 公共基础层
================================================================
只放三样东西：物理常数、异常类型、阶段之间传递的数据结构。
本模块不包含任何计算，也不 import 其他任何业务模块，
因此 athena / artemis / larch_wt 三者互不依赖，可各自单独运行。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

# ================================================================
# 1. 物理常数
# ================================================================

# 自由电子波矢与能量的换算系数（Å⁻²·eV⁻¹）
#   k = sqrt(2 m_e (E - E0) / ħ²) = sqrt(K2 * (E - E0))
#   2 m_e / ħ² = 0.2624685 Å⁻² eV⁻¹
K2 = 0.2624685

# k = 0.5123 * sqrt(E - E0) 中的前置系数
K_COEFF = float(np.sqrt(K2))

# 常用吸收边标准能量（eV），用于参考通道对齐与自检
EDGE_ENERGY = {
    ("Ti", "K"): 4966.0, ("V", "K"): 5465.0, ("Cr", "K"): 5989.0,
    ("Mn", "K"): 6539.0, ("Fe", "K"): 7112.0, ("Co", "K"): 7709.0,
    ("Ni", "K"): 8333.0, ("Cu", "K"): 8979.0, ("Zn", "K"): 9659.0,
    ("Mo", "K"): 20000.0, ("Pt", "L3"): 11564.0,
}


# ================================================================
# 2. 异常
# ================================================================

class XAFSError(Exception):
    """所有可预期的处理错误。携带出错阶段名，便于定位。"""

    def __init__(self, stage: str, message: str):
        self.stage = stage
        self.message = message
        super().__init__(f"[{stage}] {message}")


class XAFSDataError(XAFSError):
    """数据本身不合格（列角色错、边步过小、归一化失败等）。"""


class XAFSConfigError(XAFSError):
    """调用参数不合法（未知 method、手工值越界等）。"""


# ================================================================
# 3. 阶段数据结构
# ================================================================
#
# 设计原则：每阶段只往下游传「它自己定义清楚」的字段，
# 下游不回改上游对象。所有数组统一为 np.ndarray，单位写进字段注释。

@dataclass
class RawScan:
    """S1 解析产物：一次扫描的原始光强。"""
    energy: np.ndarray                 # eV
    i0: Optional[np.ndarray] = None    # 入射光强
    it: Optional[np.ndarray] = None    # 透射光强
    ifluor: Optional[np.ndarray] = None  # 荧光光强
    iref: Optional[np.ndarray] = None  # 参考通道（Fe 箔等）
    dtime: Optional[np.ndarray] = None  # 每点积分时间
    mu: Optional[np.ndarray] = None    # 已算好的 μ(E)，或文件自带 xmu
    sigma_mu: Optional[np.ndarray] = None  # 多扫描合并的逐点标准差 σ_μ(E)
    meta: dict = field(default_factory=dict)

    def __len__(self):
        return len(self.energy)


@dataclass
class NormalizedScan:
    """S3 归一化产物。"""
    energy: np.ndarray
    mu: np.ndarray
    mu_norm: np.ndarray
    pre_edge_line: np.ndarray          # 边前线性外推 y_pre(E)
    edge_step: float                   # Δμ，边跳变高度
    plateau: float                     # 边后平台均值，判据用
    e0: float
    flattened: Optional[np.ndarray] = None  # 边后展平谱（Athena flattened 约定：
                                            # E<E0 与 mu_norm 相同；E≥E0 扣掉
                                            # 边后二次背景，平台中心≈1.0）
    engine: str = "larch"              # B10-1: larch / native / mback
    mback: Optional[dict] = None       # MBACK 拟合诊断（a/b/理论边跳）


@dataclass
class ChiK:
    """S4 背景扣除产物：k 空间振荡信号。

    长度契约（容易踩坑）：
        k / chi                 k 空间数组，二者等长；
        background / mu_norm    仍在**能量网格**上，长度等于 S3 的
                                NormalizedScan.energy，而不是 len(k)。
    下游若需要 k 轴上的 μ(E)，必须先按 E(k) 重采样。
    """
    k: np.ndarray                      # Å⁻¹
    chi: np.ndarray                    # 裸 χ(k)
    background: np.ndarray             # 被扣掉的 μ0(E)，能量网格
    mu_norm: np.ndarray                # 归一化谱，能量网格
    rbkg: float
    chi_std_used: bool = False         # B10-1: 是否用理论先验谱引导 AUTOBK
    clamp_lo: Optional[float] = None   # B10-1: 低 k 端钳制权重
    clamp_hi: Optional[float] = None   # B10-1: 高 k 端钳制权重


@dataclass
class FourierTransform:
    """S5 傅里叶变换产物。"""
    r: np.ndarray                      # Å
    chi_r: np.ndarray                  # 复数 χ(R)
    kmin: float
    kmax: float
    kweight: int
    dk: float

    @property
    def magnitude(self) -> np.ndarray:
        return np.abs(self.chi_r)

    @property
    def real(self) -> np.ndarray:
        return self.chi_r.real

    @property
    def imag(self) -> np.ndarray:
        return self.chi_r.imag


@dataclass
class BackFourierTransform:
    """S5b 反傅里叶变换产物：R 窗滤波后回到 q（波数）空间。

    字段:
        q         : Å⁻¹，`larch.xafs.xftr` 的 q 网格
        chi_q     : 复数 χ(q)（chiq_re + 1j·chiq_im）
        rmin/rmax : 滤波窗口（Å）
        dr        : R 窗 Hanning 过渡宽度（Å）
        r_window  : 与正变换 R 网格等长的窗函数值（供诊断图叠加）
    """
    q: np.ndarray
    chi_q: np.ndarray
    rmin: float
    rmax: float
    dr: float
    r_window: np.ndarray

    @property
    def magnitude(self) -> np.ndarray:
        return np.abs(self.chi_q)

    @property
    def real(self) -> np.ndarray:
        return self.chi_q.real

    @property
    def imag(self) -> np.ndarray:
        return self.chi_q.imag

    @property
    def phase(self) -> np.ndarray:
        return np.angle(self.chi_q)


@dataclass
class WaveletMap:
    """S6 小波变换产物：二维幅值图 + 可选复数矩阵。

    complex_map 与 magnitude 形状完全一致；normalize=True 时二者同步
    除以同一峰值，保证 ``np.abs(complex_map)`` 与 ``magnitude`` 逐点相等。
    仅保存模值会漏掉"幅值相等但相位差 π"的模型失配（见 B5-1）。
    """
    k_centers: np.ndarray              # Å⁻¹，滑窗中心
    r: np.ndarray                      # Å
    magnitude: np.ndarray              # |W(R,k)|，形状 (len(k_centers), len(r))
    sigma: float
    complex_map: np.ndarray | None = None  # 复数 W(R,k)，同 magnitude 形状

    @property
    def phase(self) -> np.ndarray | None:
        """arg W(R,k)，单位 rad；无复数矩阵时返回 None。"""
        return None if self.complex_map is None else np.angle(self.complex_map)


@dataclass
class PathFit:
    """单条散射路径的拟合结果。"""
    label: str
    reff: float                        # 晶体学半程长度，Å
    degeneracy: float                  # 简并度：参与拟合时为拟合值 N，否则为晶体学 N
    delta_r: float
    delta_r_err: float
    sigma2: float                      # Å²
    sigma2_err: float
    n_err: Optional[float] = None      # N 的拟合误差；N 未作为变量时为 None
    c3: Optional[float] = None         # 三阶累积量 C₃ (Å³)；未作为变量时为 None
    c3_err: Optional[float] = None
    c3_at_bound: bool = False          # |C₃| ≥ 0.8×0.01 Å³：非谐项在吞噬
                                       # 相位失配，CN 不可靠
    c4: Optional[float] = None         # 四阶累积量 C₄ (Å⁴)；未作为变量时为 None
    c4_err: Optional[float] = None
    c4_at_bound: bool = False          # |C₄| ≥ 0.8×2e-4 Å⁴：非高斯展宽补偿
                                       # 到边界，CN 不可靠
    sigma2_at_bound: bool = False      # σ² 压到下界附近（≤2×下界）——
                                       # 形状失配的强证据：失配不是被 σ²
                                       # 吸收而是反向压缩，cn_reliable 需看
    nleg: int = 2                      # 路径腿数；>2 为受约束多散射路径
    ms_parent: str | None = None       # MS 路径绑定的母单散射壳层标签

    @property
    def r(self) -> float:
        return self.reff + self.delta_r


@dataclass
class FitResult:
    """S7 拟合汇总。"""
    paths: list
    s02: float
    delta_e0: float
    delta_e0_err: float
    r_factor: float
    n_independent: float               # N_idp
    n_variables: int
    kmin: float
    kmax: float
    rmin: float
    rmax: float
    s02_err: float | None = None       # S₀² 误差；仅标定模式（vary_s02）非 None
    mixture_weights: list | None = None  # mixture 权重拟合结果：
                                         # [{"group", "weight", "weight_err"}, ...]
                                         # 仅 vary_weights=True 时非 None
    parameter_names: list | None = None  # 自由变量名（与 correlation 行列对齐）
    correlation: np.ndarray | None = None  # 参数相关矩阵（对角线 1）；用于识别
                                           # N–σ²–ΔE0 强相关导致的虚假小误差
    sigma2_groups: list | None = None    # σ² 共享组：[[0, 1], [2]] 表示路径 0、1
                                         # 共享一个 σ² 变量；None 时每条路径独立
    sigma2_group_names: list | None = None  # 共享组名称，如 ["shell1", "shell2"]
    kweights: tuple | None = None        # 实际参与联合最小二乘的 k 权重元组；
                                         # 单权时也为 (2,)（向后兼容只读字段）
    r_factors_by_kweight: dict | None = None  # {kw: R-factor}，逐权重分支
    epsilon_r: float | None = None       # 高 R 区（15–25 Å）噪声估计 (Å⁻²)
    chi_square: float | None = None      # IFEFFIT 定义 χ²
    reduced_chi_square: float | None = None  # χ²_ν = χ² / max(N_idp − N_var, 1)
    lattice_strain: float | None = None  # 各向同性晶格膨胀系数 α：ΔR_i = α·Reff_i
                                         # 仅 isotropic_deltar=True 时非 None
    lattice_strain_err: float | None = None  # α 的标准误（线性传播到 ΔR_i）
    model_r: np.ndarray | None = None    # 拟合模型在数据 R 网格上的 R 轴（Å）
    model_chi_r: np.ndarray | None = None  # 最优模型的复数 χ(R)（供 back-FT/WT 残差）
    model_k: np.ndarray | None = None    # 最优模型 χ(k)（与 model_k_grid 等长）
    model_k_grid: np.ndarray | None = None  # model_k 对应的 k 网格（Å⁻¹）
    aic: float | None = None             # B12-1: Akaike 信息准则
    bic: float | None = None             # B12-1: Bayesian 信息准则
    bvs: float | None = None             # B11-2: 拟合结构键价和
    bvs_target: float | None = None      # B11-2: 目标价态（软约束目标）
    bvs_chi2: float | None = None        # B11-2: ((BVS-target)/sigma)²
    fit_space: str = "r"                 # B11-3: r / q / wt

    @property
    def overfit(self) -> bool:
        """变量数超过独立点数即为过拟合。"""
        return self.n_variables > self.n_independent


@dataclass
class PipelineResult:
    """run_xafs.process_sample 的全流程产物：S1–S6 各阶段对象 + 关键结果。
    取代原先的裸 dict，字段访问拼写错误在导入期即暴露。"""
    scans: list                       # List[RawScan]，对齐后的扫描列表
    merged: RawScan                   # 合并谱
    e0: dict                          # determine_e0 的输出
    norm: NormalizedScan
    rbkg: dict                        # optimize_rbkg 选中的最优行
    chi_k: ChiK
    ft: FourierTransform
    wt: WaveletMap
    r1: float                         # |χ(R)| 第一峰位 (Å)，无峰时为 NaN
    r_wt: float                       # 小波第一峰位 (Å)
    k_wt: float                       # 小波第一峰处的 k (Å⁻¹)
    absorber_symbol: Optional[str] = None  # 吸收元素符号（用户指定；S6.5/S7 消费）
    outdir: Optional[Path] = None     # 落盘目录（由入口层赋值）
    stem: str = ""                    # 输出文件名前缀（由入口层赋值）
    rbkg_table: list = field(default_factory=list)  # optimize_rbkg 全网格（诊断用）
    self_absorption: Optional[dict] = None  # S4.5 荧光自吸收校正元数据（Booth）
    structure: Optional[dict] = None  # S6.5 结构识别结果（由 xafs_fit.run 填充——拟合链唯一写点）
    fit: Optional[dict] = None        # S7 拟合结果（由 xafs_fit.run 填充）
    xanes_features: dict = field(default_factory=dict)  # B2-1 XANES 定量特征/颈部指纹
    snr: Optional[dict] = None        # B9-2 高 k 信噪比诊断（high_k_snr 输出）
