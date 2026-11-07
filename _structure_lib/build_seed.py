# -*- coding: utf-8 -*-
"""
build_seed.py — 生成 Fe 种子 motif 库（一次性，可重复执行）
================================================================
覆盖核数分类的代表条目：

    crystal  lepidocrocite   γ-FeOOH（从仓库已有 CIF 导入，走完整入库通道）
             hematite        α-Fe2O₃（R-3c，标准文献值）
             magnetite_oct   Fe3O4（Fd-3m，多位点吸收：选 16d 八面体位，
                              占 Fe 总量 2/3）
             magnetite_tet   Fe3O4（同骨架，8a 四面体位，占 1/3）
             magnetite_mixed 近黑石多位点联合：χ=(1/3)χ_8a+(2/3)χ_16d
             wustite         FeO（Fm-3m 理想岩盐）
             goethite        α-FeOOH（Pnma，COD 2211652 单晶）
             fe_bcc          金属 Fe 箔（Im-3m，S₀² 标定源/颗粒代表）
             以及 Pt/pt_fcc   金属 Pt 箔（Fm-3m，Pt L3 S₀² 标定源）
    cluster  fen4_planar     FeN₄ 平面（Fe-N-C 单点假设，构造几何）
             fe2_mu_oxo      μ-oxo 双核 [Fe2O]（ZSM-5 活性位候选，构造几何）
             fe2n6_dual      Fe2N6 双原子位（构造几何，建议以 DFT 精化替换）

用法: python _structure_lib/build_seed.py
依赖: 在仓库根目录的模块（artemis/structure_lib），脚本自行加路径。
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from structure_lib import Motif, import_cif, register  # noqa: E402


def main():
    made = []

    # ---- 晶体：lepidocrocite（从已有 CIF 走导入器，验证入库通道）----
    m = import_cif(ROOT / "lepidocrocite.cif", "lepidocrocite", "Fe",
                   nuclearity="crystal",
                   source={"type": "cif", "file": "lepidocrocite.cif",
                           "ref": "COD-9002763 等"})
    m.notes = "γ-FeOOH 标样。管线基准结构。"
    made.append(register(m, overwrite=True))

    # ---- 晶体：bcc Fe（金属箔/颗粒代表）----
    m = Motif(id="fe_bcc", absorber="Fe", kind="crystal", nuclearity="crystal",
              source={"type": "manual", "ref": "a=2.8664 Å, Im-3m (#229)"},
              cell={"a": 2.8664}, spacegroup="Im-3m",
              asym_sites=[{"el": "Fe", "x": 0.0, "y": 0.0, "z": 0.0,
                           "wyckoff": "2a"}],
              notes="金属 Fe。S₀² 标定源；也作 α-Fe 颗粒的体相极限。")
    made.append(register(m, overwrite=True))

    # ---- 晶体：hematite α-Fe2O3 ----
    m = Motif(id="hematite", absorber="Fe", kind="crystal", nuclearity="crystal",
              source={"type": "manual", "ref": "R-3c (#167), a=5.0356 c=13.7489"},
              cell={"a": 5.0356, "c": 13.7489, "gamma": 120},
              spacegroup="R-3c",
              asym_sites=[
                  {"el": "Fe", "x": 0.0, "y": 0.0, "z": 0.3553,
                   "wyckoff": "12c"},
                  {"el": "O", "x": 0.3059, "y": 0.0, "z": 0.25,
                   "wyckoff": "18e"}],
              notes="α-Fe2O3。Fe-O 1.94x3/2.11x3，Fe-Fe 2.89/2.97/3.36。")
    made.append(register(m, overwrite=True))

    # ---- 晶体：magnetite Fe3O4（多位点吸收，选 16d 八面体位）----
    m = Motif(id="magnetite_oct", absorber="Fe", kind="crystal",
              nuclearity="crystal",
              source={"type": "manual",
                      "ref": "Fd-3m (#227) origin-2, a=8.396"},
              cell={"a": 8.396}, spacegroup="F d -3 m:2",
              asym_sites=[
                  {"el": "Fe", "x": 0.125, "y": 0.125, "z": 0.125,
                   "wyckoff": "8a"},
                  {"el": "Fe", "x": 0.5, "y": 0.5, "z": 0.5,
                   "wyckoff": "16d"},
                  {"el": "O", "x": 0.2548, "y": 0.2548, "z": 0.2548,
                   "wyckoff": "32e"}],
              absorber_fract=[0.5, 0.5, 0.5],
              notes="Fe3O4 反尖晶石。吸收位=16d 八面体（占 2/3）；"
                    "混合价 Fe2+/Fe3+，XANES 边位低于 Fe3O4 名义+2.67。"
                    "完整多位点模型见 magnetite_tet + magnetite_mixed。")
    made.append(register(m, overwrite=True))

    # ---- 晶体：magnetite 8a 四面体位（骨架与 magnetite_oct 相同）----
    m = Motif(id="magnetite_tet", absorber="Fe", kind="crystal",
              nuclearity="crystal",
              source={"type": "manual",
                      "ref": "Fd-3m (#227) origin-2, a=8.396"},
              cell={"a": 8.396}, spacegroup="F d -3 m:2",
              asym_sites=[
                  {"el": "Fe", "x": 0.125, "y": 0.125, "z": 0.125,
                   "wyckoff": "8a"},
                  {"el": "Fe", "x": 0.5, "y": 0.5, "z": 0.5,
                   "wyckoff": "16d"},
                  {"el": "O", "x": 0.2548, "y": 0.2548, "z": 0.2548,
                   "wyckoff": "32e"}],
              absorber_fract=[0.125, 0.125, 0.125],
              notes="Fe3O4 反尖晶石。吸收位=8a 四面体（占 1/3）；"
                    "第一壳 Fe-O 1.886×4。与 magnetite_oct 共用同一"
                    "晶体骨架，只换吸收位点（asym_sites 未变）。")
    made.append(register(m, overwrite=True))

    # ---- 混合位点：磁铁矿 χ = (1/3)·χ_8a + (2/3)·χ_16d ----
    m = Motif(id="magnetite_mixed", absorber="Fe", kind="mixture",
              nuclearity="mixed-site",
              source={"type": "composite",
                      "ref": "8a(1/3) + 16d(2/3) 位点联合；多位点吸收"
                             "的关键是权重=位点占据，与化学式一致"},
              components=[{"id": "magnetite_tet", "weight": 1.0 / 3.0},
                          {"id": "magnetite_oct", "weight": 2.0 / 3.0}],
              notes="多位点吸收：实验 χ(R) = 两个吸收位点 χ 的线性组合。"
                    "各位点共用 S₀²/σ²/ΔE₀（指纹级近似）；定量拟合时"
                    "Fe-O 子壳层 1.886/2.060 分开给 N。")
    made.append(register(m, overwrite=True))

    # ---- 晶体：wustite FeO（岩盐，理想化学计量）----
    m = Motif(id="wustite", absorber="Fe", kind="crystal",
              nuclearity="crystal",
              source={"type": "manual",
                      "ref": "FeO rock salt, Fm-3m (#225), a=4.310 Å"},
              cell={"a": 4.310}, spacegroup="Fm-3m",
              asym_sites=[
                  {"el": "Fe", "x": 0.0, "y": 0.0, "z": 0.0,
                   "wyckoff": "4a"},
                  {"el": "O", "x": 0.5, "y": 0.5, "z": 0.5,
                   "wyckoff": "4b"}],
              absorber_fract=[0.0, 0.0, 0.0],
              notes="FeO 理想岩盐结构（实际常为非整比 Fe(1-x)O）；"
                    "第一壳 Fe-O 2.155×6，Fe-Fe 3.048×12。")
    made.append(register(m, overwrite=True))

    # ---- 晶体：goethite α-FeOOH（COD 2211652 单晶）----
    m = Motif(id="goethite", absorber="Fe", kind="crystal",
              nuclearity="crystal",
              source={"type": "cod", "id": "2211652",
                      "ref": "Yang et al. Acta Cryst E 62 (2006) i250; "
                             "Pnma, a=4.5979 b=9.9510 c=3.0178"},
              cell={"a": 4.5979, "b": 9.9510, "c": 3.0178},
              spacegroup="P b n m",
              asym_sites=[
                  {"el": "Fe", "x": 0.04893, "y": 0.853657, "z": 0.25,
                   "wyckoff": "4c"},
                  {"el": "O", "x": 0.7057, "y": 0.19914, "z": 0.25,
                   "wyckoff": "4c"},
                  {"el": "O", "x": 0.1987, "y": 0.05298, "z": 0.25,
                   "wyckoff": "4c"}],
              absorber_fract=[0.04893, 0.853657, 0.25],
              notes="α-FeOOH（针铁矿）。第一壳 Fe-O 1.95×3 + 2.10×3，"
                    "Fe-Fe 3.02×2 / 3.31×2 / 3.44×4。H 坐标未入模"
                    "（背散射可忽略，且粉晶 H 常不可靠）。")
    made.append(register(m, overwrite=True))

    # ---- 团簇：FeN4 平面（Fe-N-C 单点假设）----
    # 约定：吸收 Fe 在原点（不入列），sites 只放近邻。
    n = 1.93
    sites = [(n, 0.0, 0.0, "N"), (-n, 0.0, 0.0, "N"),
             (0.0, n, 0.0, "N"), (0.0, -n, 0.0, "N")]
    for sx in (1, -1):
        for sy in (1, -1):
            sites.append((sx * 2.09, sy * 2.09, 0.0, "C"))  # 第二壳 ~2.95
    m = Motif(id="fen4_planar", absorber="Fe", kind="cluster", nuclearity="m1",
              source={"type": "constructed",
                      "ref": "FePc/DFT 典型 Fe-N≈1.93；建议以 DFT 精化替换"},
              sites=[list(s) for s in sites],
              notes="平面 FeN4 + 4 个第二壳 C。N-Fe-N 共线 MS 路径显著，"
                    "用于拟合时须放开 legs==2 限制（见 structure_id 头注）。")
    made.append(register(m, overwrite=True))

    # ---- 团簇：μ-oxo 双核 [Fe2O]（ZSM-5 活性位候选）----
    # Fe1 在原点；Fe2 在 +x 3.32 处；桥氧使 Fe-O=1.77、∠Fe-O-Fe≈139°。
    fe2 = 3.32
    sites = [(fe2, 0.0, 0.0, "Fe"),
             (fe2 / 2, 0.62, 0.0, "O")]                     # 桥氧 @1.77
    for dx, dy, dz in ((-1.0, 1.732, 0.0), (-1.0, -1.732, 0.0),
                       (-1.0, 0.0, 1.732)):                 # Fe1 端基 O @2.0
        sites.append((dx, dy, dz, "O"))
    for dx, dy, dz in ((1.0, 1.732, 0.0), (1.0, -1.732, 0.0),
                       (1.0, 0.0, 1.732)):                  # Fe2 端基 O @2.0
        sites.append((fe2 + dx, dy, dz, "O"))
    m = Motif(id="fe2_mu_oxo", absorber="Fe", kind="cluster", nuclearity="m2",
              source={"type": "constructed",
                      "ref": "Fe-O 1.77/∠139° 取自 μ-oxo 双核文献典型值"},
              sites=[list(s) for s in sites],
              notes="判定要点：Fe-Fe≈3.32 且 Fe-Fe CN=1；与 Fe2O3 颗粒的"
                    "边共享峰位近，必须靠 CN 定量区分（见 skill 文档）。")
    made.append(register(m, overwrite=True))

    # ---- 团簇：Fe2N6 双原子位 ----
    # Fe1 在原点，Fe2 在 +x 2.45（Fe-Fe 直接键）；每 Fe 三个端基 N @1.9。
    fe2 = 2.45
    import math
    sites = [(fe2, 0.0, 0.0, "Fe")]
    for deg in (90.0, 210.0, 330.0):                        # Fe1 端基 N
        th = math.radians(deg)
        sites.append((1.9 * math.cos(th), 1.9 * math.sin(th), 0.0, "N"))
    for deg in (30.0, 90.0, 330.0):                         # Fe2 端基 N（朝外）
        th = math.radians(deg)
        sites.append((fe2 + 1.9 * math.cos(th), 1.9 * math.sin(th), 0.0, "N"))
    m = Motif(id="fe2n6_dual", absorber="Fe", kind="cluster", nuclearity="m2",
              source={"type": "constructed",
                      "ref": "Fe-Fe≈2.45 直接键；N 排布近似，须 DFT 精化"},
              sites=[list(s) for s in sites],
              notes="Fe-N-C 双原子位候选。Fe-Fe≈2.45 与 α-Fe(2.48)/"
                    "Fe5C2(~2.6) 同区间——误判高危，需第二壳+XANES 联判。")
    made.append(register(m, overwrite=True))

    # ---- B16-2：Co / Ni / Cu / Pt 催化基准 Motif ----
    # 第一壳全部落在 1.85–2.31 Å，与金属箔第一壳（2.49–2.77 Å）及
    # NiO 主导峰（2.95 Å）错开，保证压力测试 top1 不被抢走。

    # 晶体：Cu2O 赤铜矿。Pn-3m 必须显式 origin choice 1：gemmi 默认
    # choice 2 会把 O(2a) 误展开成 4 个 O，给出 Cu-O×4 的错误指纹。
    m = Motif(id="cu2o_cuprite", absorber="Cu", kind="crystal",
              nuclearity="crystal",
              source={"type": "manual",
                      "ref": "Pn-3m (#224, origin choice 1), a=4.2696 Å"},
              cell={"a": 4.2696}, spacegroup="P n -3 m:1",
              asym_sites=[
                  {"el": "Cu", "x": 0.25, "y": 0.25, "z": 0.25,
                   "wyckoff": "4b"},
                  {"el": "O", "x": 0.0, "y": 0.0, "z": 0.0,
                   "wyckoff": "2a"}],
              absorber_fract=[0.25, 0.25, 0.25],
              notes="Cu2O 赤铜矿。第一壳 Cu-O 2×1.849（线性二配位），"
                    "第二壳 Cu-Cu 12×3.019。")
    made.append(register(m, overwrite=True))

    # 团簇：平面 M-N4 单原子位 + 碳壳（M = Co / Ni / Cu）。
    # 吸收原子在原点，sites 只放近邻；C 放在 (±c, ±c, 0)，
    # 使 C 壳距离 = sqrt(2)·c = 目标值。
    for el, motif_id, r_m, r_c in (
            ("Co", "con4_planar", 1.96, 2.88),
            ("Ni", "nin4_planar", 1.92, 2.85),
            ("Cu", "cun4_planar", 1.95, 2.88)):
        c = r_c / math.sqrt(2.0)
        sites = [(r_m, 0.0, 0.0, "N"), (-r_m, 0.0, 0.0, "N"),
                 (0.0, r_m, 0.0, "N"), (0.0, -r_m, 0.0, "N")]
        for sx in (1, -1):
            for sy in (1, -1):
                sites.append((sx * c, sy * c, 0.0, "C"))
        m = Motif(id=motif_id, absorber=el, kind="cluster", nuclearity="m1",
                  source={"type": "constructed",
                          "ref": f"{el}-N4 平面，M-N={r_m} Å；"
                                 f"碳壳 {r_c} Å（构造几何，DFT 可精化）"},
                  sites=[list(s) for s in sites],
                  notes=f"M-N-C 单原子催化剂候选：平面 {el}N4 + 4 个"
                        f"第二壳 C。第一壳 {el}-N={r_m} Å 与金属箔"
                        f"（≥2.49 Å）错开。")
        made.append(register(m, overwrite=True))

    # 团簇：平面四方 [PtCl4]²⁻ + α-PtCl2 链间第二壳 Cl。
    # 真实标样 Pt-L3-PlatinumChloride.xdi 的 Sample.name 是 PtCl2：
    # 第一壳 Pt-Cl 4×2.31，|χ(R)| 在 3.13 Å（表观）还有显著峰；
    # 第二壳 4×Cl@3.45 与自由拟合的 3.44 Å 一致（Rf<0.01），
    # 第一壳 CN 仍是主结论。
    r2 = 3.45
    s2 = r2 / math.sqrt(3.0)
    ptcl_sites = [[2.31, 0.0, 0.0, "Cl"], [-2.31, 0.0, 0.0, "Cl"],
                  [0.0, 2.31, 0.0, "Cl"], [0.0, -2.31, 0.0, "Cl"]]
    for d in ((s2, s2, s2), (s2, -s2, -s2), (-s2, s2, -s2),
              (-s2, -s2, s2)):
        ptcl_sites.append([d[0], d[1], d[2], "Cl"])
    m = Motif(id="ptcl4_square", absorber="Pt", kind="cluster",
              nuclearity="m1",
              source={"type": "constructed",
                      "ref": "平面四方 [PtCl4]2-，Pt-Cl≈2.31 Å；第二壳 "
                             "Cl≈3.45 Å 对应 α-PtCl2 链间近邻；"
                             "匹配 Pt-L3-PlatinumChloride.xdi 标样"},
              sites=ptcl_sites,
              notes="Pt(II) 平面四方氯配合物。第一壳 Pt-Cl 4×2.31 Å；"
                    "第二壳 4×Cl@3.45 Å 解释实测 3.13 Å（表观）峰。"
                    "用于区分 Pt 箔（Pt-Pt 2.774 Å）与 PtCl2 标样。")
    made.append(register(m, overwrite=True))

    # ---- 金属 Pt fcc（Pt L3 标定源 / 颗粒代表）----
    # 与 Fe/Co/Ni/Cu 种子同地位：真实室温 Pt 箔晶格常数，第一壳
    # Pt-Pt = a/√2 = 2.774 Å × 12（fcc 配位数）。
    m = Motif(id="pt_fcc", absorber="Pt", kind="crystal", nuclearity="crystal",
              source={"type": "manual",
                      "ref": "a=3.9231 Å, Fm-3m (#225), 室温 Pt 箔"},
              cell={"a": 3.9231}, spacegroup="Fm-3m",
              asym_sites=[{"el": "Pt", "x": 0.0, "y": 0.0, "z": 0.0,
                           "wyckoff": "4a"}],
              notes="金属 Pt。Pt L3（11564 eV）S₀² 标定源；"
                    "第一壳 Pt-Pt 2.774×12，第二壳 3.924×6。")
    made.append(register(m, overwrite=True))

    # ---- B18-2：Au / Mo / Cr / Zn 金属箔种子（运行时 --s02 标定用）----
    # 静态 absorber.py 保持 s02=None（#60 契约）；S₀² 只经 --s02 运行时
    # 覆盖，标定逻辑与 Fe/Co/Ni/Cu/Pt 完全同构。
    m = Motif(id="au_fcc", absorber="Au", kind="crystal", nuclearity="crystal",
              source={"type": "manual",
                      "ref": "a=4.0782 Å, Fm-3m (#225), 室温 Au 箔"},
              cell={"a": 4.0782}, spacegroup="Fm-3m",
              asym_sites=[{"el": "Au", "x": 0.0, "y": 0.0, "z": 0.0,
                           "wyckoff": "4a"}],
              notes="金属 Au。第一壳 Au-Au 2.8837×12，第二壳 4.0782×6。")
    made.append(register(m, overwrite=True))

    m = Motif(id="mo_bcc", absorber="Mo", kind="crystal", nuclearity="crystal",
              source={"type": "manual",
                      "ref": "a=3.1470 Å, Im-3m (#229), 室温 Mo 箔"},
              cell={"a": 3.1470}, spacegroup="Im-3m",
              asym_sites=[{"el": "Mo", "x": 0.0, "y": 0.0, "z": 0.0,
                           "wyckoff": "2a"}],
              notes="金属 Mo（bcc，N₁=8）。第一壳 Mo-Mo 2.7254×8，"
                    "第二壳 3.1470×6。")
    made.append(register(m, overwrite=True))

    m = Motif(id="cr_bcc", absorber="Cr", kind="crystal", nuclearity="crystal",
              source={"type": "manual",
                      "ref": "a=2.8846 Å, Im-3m (#229), 室温 Cr 箔"},
              cell={"a": 2.8846}, spacegroup="Im-3m",
              asym_sites=[{"el": "Cr", "x": 0.0, "y": 0.0, "z": 0.0,
                           "wyckoff": "2a"}],
              notes="金属 Cr（bcc）。第一壳 Cr-Cr 2.4981×8。")
    made.append(register(m, overwrite=True))

    m = Motif(id="zn_hcp", absorber="Zn", kind="crystal", nuclearity="crystal",
              source={"type": "manual",
                      "ref": "a=2.6649 c=4.9468 Å, P63/mmc (#194), 室温 Zn"},
              cell={"a": 2.6649, "b": 2.6649, "c": 4.9468,
                    "alpha": 90, "beta": 90, "gamma": 120},
              spacegroup="P 63/m m c",
              asym_sites=[{"el": "Zn", "x": 1.0 / 3.0, "y": 2.0 / 3.0,
                           "z": 0.25, "wyckoff": "2c"}],
              notes="金属 Zn（hcp）。面内 Zn-Zn 2.6649×6，"
                    "面外 2.9131×6（c/a=1.856 非理想）。")
    made.append(register(m, overwrite=True))

    # ---- B20-1：电催化核心氧化物 / 羟基氧化物 ----
    # 全部用仓库真实 BMM 标准谱预验证过的晶体学参数；第一壳落在
    # 1.92–2.13 Å，与金属箔（≥2.49 Å）错开，不抢金属 top1。

    # CoO 岩盐（Fm-3m）：Co-O 6×2.1306，Co-Co 12×3.0131。
    m = Motif(id="coo_rocksalt", absorber="Co", kind="crystal",
              nuclearity="crystal",
              source={"type": "manual",
                      "ref": "Fm-3m (#225), a=4.2612 Å"},
              cell={"a": 4.2612}, spacegroup="Fm-3m",
              asym_sites=[
                  {"el": "Co", "x": 0.0, "y": 0.0, "z": 0.0,
                   "wyckoff": "4a"},
                  {"el": "O", "x": 0.5, "y": 0.5, "z": 0.5,
                   "wyckoff": "4b"}],
              notes="CoO 岩盐氧化物。第一壳 Co-O 6×2.131 Å、"
                    "第二壳 Co-Co 12×3.013 Å；Co-K-CoO 标样基准相。")
    made.append(register(m, overwrite=True))

    # Co3O4 尖晶石双位点：8a 四面体（1/3）+ 16d 八面体（2/3）。
    co3o4_asym = [
        {"el": "Co", "x": 0.125, "y": 0.125, "z": 0.125,
         "wyckoff": "8a"},
        {"el": "Co", "x": 0.5, "y": 0.5, "z": 0.5,
         "wyckoff": "16d"},
        {"el": "O", "x": 0.263, "y": 0.263, "z": 0.263,
         "wyckoff": "32e"}]
    m = Motif(id="co3o4_tet", absorber="Co", kind="crystal",
              nuclearity="crystal",
              source={"type": "manual",
                      "ref": "Fd-3m:2 (#227 origin 2), a=8.084 Å"},
              cell={"a": 8.084}, spacegroup="F d -3 m:2",
              asym_sites=co3o4_asym, absorber_fract=[0.125, 0.125, 0.125],
              notes="Co3O4 四面体 8a 位（占 Co 的 1/3）。"
                    "第一壳 Co-O 4×1.932 Å；与八面体位共同构成混合谱。")
    made.append(register(m, overwrite=True))

    m = Motif(id="co3o4_oct", absorber="Co", kind="crystal",
              nuclearity="crystal",
              source={"type": "manual",
                      "ref": "Fd-3m:2 (#227 origin 2), a=8.084 Å"},
              cell={"a": 8.084}, spacegroup="F d -3 m:2",
              asym_sites=co3o4_asym, absorber_fract=[0.5, 0.5, 0.5],
              notes="Co3O4 八面体 16d 位（占 Co 的 2/3）。"
                    "第一壳 Co-O 6×1.922 Å、第二壳 Co-Co 6×2.858 Å。")
    made.append(register(m, overwrite=True))

    m = Motif(id="co3o4_mixed", absorber="Co", kind="mixture",
              nuclearity="mixed-site",
              source={"type": "manual",
                      "ref": "Co3O4 尖晶石 8a:16d = 1:2 位点加权"},
              components=[{"id": "co3o4_tet", "weight": 1.0 / 3.0},
                          {"id": "co3o4_oct", "weight": 2.0 / 3.0}],
              notes="Co3O4 多位点联合：χ=(1/3)χ_tet+(2/3)χ_oct。"
                    "真实 Co-K-Co3O4 标样盲测 top1（score≈0.025）。")
    made.append(register(m, overwrite=True))

    # β-Ni(OH)2 水镁石（P-3m1）：Ni-O 6×2.070，Ni-Ni 6×3.126。
    m = Motif(id="ni_oh2_brucite", absorber="Ni", kind="crystal",
              nuclearity="crystal",
              source={"type": "manual",
                      "ref": "P-3m1 (#164), a=3.126 c=4.605 Å"},
              cell={"a": 3.126, "b": 3.126, "c": 4.605,
                    "alpha": 90, "beta": 90, "gamma": 120},
              spacegroup="P -3 m 1",
              asym_sites=[
                  {"el": "Ni", "x": 0.0, "y": 0.0, "z": 0.0,
                   "wyckoff": "1a"},
                  {"el": "O", "x": 1.0 / 3.0, "y": 2.0 / 3.0, "z": 0.22,
                   "wyckoff": "2d"}],
              notes="β-Ni(OH)2 水镁石层状氢氧化物。第一壳 Ni-O 6×2.070、"
                    "第二壳 Ni-Ni 6×3.126；H 未入模（背散射可忽略）。")
    made.append(register(m, overwrite=True))

    # CuO 单斜黑铜矿（C2/c）：Cu-O 4×1.956。
    m = Motif(id="cuo_tenorite", absorber="Cu", kind="crystal",
              nuclearity="crystal",
              source={"type": "manual",
                      "ref": "C2/c (#15), a=4.6837 b=3.4226 c=5.1288 Å, "
                             "beta=99.54°"},
              cell={"a": 4.6837, "b": 3.4226, "c": 5.1288,
                    "alpha": 90, "beta": 99.54, "gamma": 90},
              spacegroup="C 1 2/c 1",
              asym_sites=[
                  {"el": "Cu", "x": 0.25, "y": 0.25, "z": 0.0,
                   "wyckoff": "4c"},
                  {"el": "O", "x": 0.0, "y": 0.4184, "z": 0.25,
                   "wyckoff": "4e"}],
              absorber_fract=[0.25, 0.25, 0.0],
              notes="CuO 黑铜矿。第一壳 Cu-O 4×1.956 Å（畸变平面四方），"
                    "Cu-Cu 2.900/3.083 Å；Cu-CuO 标样基准相。")
    made.append(register(m, overwrite=True))

    # β-PtO2 层状六方参照相（P-3m1）：Pt-O 6×2.047，Pt-Pt 6×3.100。
    m = Motif(id="pto2_beta", absorber="Pt", kind="crystal",
              nuclearity="crystal",
              source={"type": "manual",
                      "ref": "P-3m1 (#164), a=3.100 c=4.161 Å"},
              cell={"a": 3.100, "b": 3.100, "c": 4.161,
                    "alpha": 90, "beta": 90, "gamma": 120},
              spacegroup="P -3 m 1",
              asym_sites=[
                  {"el": "Pt", "x": 0.0, "y": 0.0, "z": 0.0,
                   "wyckoff": "1a"},
                  {"el": "O", "x": 1.0 / 3.0, "y": 2.0 / 3.0, "z": 0.24,
                   "wyckoff": "2d"}],
              notes="β-PtO2 层状六方参照相。第一壳 Pt-O 6×2.047、"
                    "第二壳 Pt-Pt 6×3.100；用于 Pt 氧化物/氢氧化物鉴别。")
    made.append(register(m, overwrite=True))

    print("已注册 motif:")
    for p in made:
        try:
            label = p.relative_to(ROOT)
        except ValueError:
            label = p
        print(f"  {label}")


if __name__ == "__main__":
    main()
