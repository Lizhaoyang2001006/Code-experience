# -*- coding: utf-8 -*-
"""
absorber.py — 吸收元素配置（元素泛化的唯一入口）
================================================================
物理计算层（athena/artemis/larch_wt）本身元素无关；元素特异的全在这里：

    symbol        元素符号
    edge          吸收边（"K" / "L3"）
    e0_edge       边的标称能量 (eV)，E0 初值与对齐基准
    s02           S₀² 标定值。**None = 未标定**——此时自动拟合必须拒绝
                  启动（S₀² 与 N 全相关，没有标定就没有 CN 定量）
    kmax_ceiling  k 上限 (Å⁻¹)。L 边被下一条边硬截断
                  （如 Pt L₃ 受 L₂ 边限制），K 边取实用上限
    next_edge     下一条边的 (名称, 能量)，供诊断展示

新元素接入流程：
    1. 查边能量填表（s02=None）
    2. 测该元素标准样（金属箔或已知结构），标定 s02 填入
    3. 向 _structure_lib/<symbol>/ 挂 motif
    代码本身零改动——这是泛化设计是否成立的验收标准。
"""
from __future__ import annotations

from dataclasses import dataclass

from xafs_core import XAFSConfigError


@dataclass(frozen=True)
class AbsorberConfig:
    symbol: str
    z: int
    edge: str
    e0_edge: float
    s02: float | None
    kmax_ceiling: float
    next_edge: tuple | None = None

    @property
    def calibrated(self) -> bool:
        return self.s02 is not None


# K 边 kmax_ceiling 取 15（数据质量通常先于物理极限恶化）；
# L 边按下一吸收边估算后取保守值。
ABSORBERS = {
    # 生产值 0.803：氧化物价态验收锚定的跨材料折中值（原始流程不可复原）。
    # 金属 Fe 交叉标定（2026-11-10）: BMM Fe-K-IronFoil.xdi, fe_bcc 单壳
    # N=8, E0=half_height → 0.6566（3 窗口 0.6471/0.6566/0.6676，极差
    # 0.0205，Rf≤0.0422；IIT-Lytle 5 µm 箔交叉 0.6335）。该值只适用于
    # 金属 Fe，不替换 0.803——否则 Fe 氧化物 CN 会高约 22%。
    "Fe": AbsorberConfig("Fe", 26, "K", 7112.0, 0.803, 15.0),
    "Co": AbsorberConfig("Co", 27, "K", 7709.0, 0.7308, 15.0),
    # 标定: Co 箔 co_hcp N=12, Rf=0.0001, 2026-11-05
    "Ni": AbsorberConfig("Ni", 28, "K", 8333.0, 0.7795, 15.0),
    # 标定: Ni 箔 ni_fcc N=12, Rf=0.0002, 2026-11-05
    "Cu": AbsorberConfig("Cu", 29, "K", 8979.0, 0.8670, 15.0),
    # 标定: Cu 箔 cu_fcc N=12, Rf=0.0006, 2026-11-05
    "Pt": AbsorberConfig("Pt", 78, "L3", 11564.0, 0.8075, 16.0,
    # 标定: Pt 箔 pt_fcc N=12, Rf=0.0004, 2026-11-05
                         next_edge=("L2", 13273.0)),
    # ---- Phase III：常用催化/环境元素。全部 s02=None（未标定），
    # 自动 CN 定量必须拒绝；标定后逐元素填入实测值。
    "Cr": AbsorberConfig("Cr", 24, "K", 5989.0, None, 15.0),
    "Mn": AbsorberConfig("Mn", 25, "K", 6539.0, None, 15.0),
    "Zn": AbsorberConfig("Zn", 30, "K", 9659.0, None, 15.0),
    "Mo": AbsorberConfig("Mo", 42, "K", 20000.0, None, 16.0),
    "Ru": AbsorberConfig("Ru", 44, "K", 22117.0, None, 16.0),
    "Rh": AbsorberConfig("Rh", 45, "K", 23220.0, None, 16.0),
    "Pd": AbsorberConfig("Pd", 46, "K", 24350.0, None, 16.0),
    "Ag": AbsorberConfig("Ag", 47, "K", 25514.0, None, 16.0),
    "Ir": AbsorberConfig("Ir", 77, "L3", 11215.0, None, 16.0,
                         next_edge=("L2", 12824.0)),
    "Au": AbsorberConfig("Au", 79, "L3", 11919.0, None, 16.0,
                         next_edge=("L2", 13734.0)),
    "Ce": AbsorberConfig("Ce", 58, "L3", 5723.0, None, 12.0,
                         next_edge=("L2", 6164.0)),
    "U": AbsorberConfig("U", 92, "L3", 17166.0, None, 16.0,
                        next_edge=("L2", 20948.0)),
}


# 荧光自吸收校正默认荧光线：K→Kα，L3→Lα，L2→Lβ1，L1→Lβ3。
# 用户可用 --self-absorption-line 覆盖（例如探测器只收 Kβ）。
_FLUORESCENCE_LINES = {"K": "Ka", "L3": "La", "L2": "Lb1", "L1": "Lb3"}


def default_fluorescence_line(edge: str) -> str:
    try:
        return _FLUORESCENCE_LINES[edge.upper()]
    except KeyError:
        raise XAFSConfigError(
            "S4.5", f"边 {edge} 没有默认荧光线；用 --self-absorption-line "
                    f"显式指定（K→Ka, L3→La, L2→Lb1, L1→Lb3）") from None


def get_absorber(symbol: str) -> AbsorberConfig:
    """按元素符号取配置；未收录元素显式报错而不是静默猜参数。"""
    key = symbol.strip().capitalize()
    if key not in ABSORBERS:
        raise XAFSConfigError(
            "absorber", f"元素 {symbol} 未收录。已收录: {sorted(ABSORBERS)}。"
                        f"新元素需按 absorber.py 头注的三步流程接入")
    return ABSORBERS[key]
