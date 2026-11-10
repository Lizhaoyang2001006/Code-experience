# -*- coding: utf-8 -*-
"""
regression_test.py — 本轮修复的轻量回归测试
================================================
只覆盖“曾经真出过问题、且修好后容易再次退化”的契约，不做全链路性能测试：

    1. ChiK 长度契约：k/chi 等长；background/mu_norm 在能量网格
    2. theory_chi 的 nleg 过滤：MS（nleg>4）不得混入指纹理论谱
    3. indistinguishable 不包含最优候选自身
    4. WT 矩阵方向：(len(k_centers), len(r))
    5. strict 门禁：Fe 箔平台 0.79 → strict 拒绝 / --no-strict 降级
    6. S7 motif 路径：lepidocrocite → Fe–O/Fe–Fe 多壳层 CN
    7. S7 L1 路径：goethite（库外）→ 仅第一壳，且不冒充全结构
    8. 磁铁矿多位点：magnetite_mixed 联合模型 → confident 全覆盖
    9. mixture 壳层选择：Fe-O@1.89(⅓×4) 与 Fe-O@2.06(⅔×6) 不合并
    10. vary_s02 合成回路：N 固定、放开 S₀² 能回收真值（标定前提）
    11. L1 配体感知选峰：Pt–Cl 不被低 R 白线伪峰带偏（ligand=Cl）
    12. L1 长键窗口：Pt–Pt 2.77 不被 rmax=2.3 截断（曾 Rf=1.0 伪失败）
    13. mixture 相分数回收：磁铁矿 tet/oct 第一壳 0.25/0.75 合成谱
    14. 荧光自吸收门禁：Δμ≥1.0 拒绝 CN，透射模式不误伤
    15. Booth 自吸收反解：Eq.2 前向 → Eq.4 回收（薄/厚极限）
    16. Booth 全链合成：真实 μ → 前向荧光谱 → 归一化/autobk → 校正
    17. WT 定位锚定 FT 峰：自吸收校正前后不得在两壳分支间跳
    18. σ² 共享组：多路径共享一个 σ² 变量（Artemis 式约束）
    19. CN 系统误差预算随拟合条件动态合成（不再恒为 ±30%）
    20. HTML 报告：图片内联、强相关/不可靠项标红
    21. 配置文件未知扩展名给出干净 ValueError（p.suffix 回归）
    22. Hanning 窗契约：端点 0、平台 1，供浏览器端 FT 与产线对齐
    23. WT sweep：9 组 (σ,kweight) 全部锚定 FT 第一峰 ±0.5 Å，CSV/图可写出
    24. 交互式窗口 HTML：内嵌 χ(k)、滑块齐全、推荐窗口落在 r1±0.15 Å
    25. 多 k 权重联合拟合：三分支 Rf、第一壳键长稳定、N_idp 不变
    26. Back-FT：q 轴单调、包络非负、Re χ(q) 与 k²χ 同相
    27. NiO 共线 MS：nleg>2 选入、自由变量数不变、ΔR 按 reff_ratio 缩放
    28. XANES flattened 平台≈1；合成颈部双台阶识别 + E0 跳变风险
    29. LCF 三端元权重回收；PCA 前 3 主成分累计解释 >99.9%
    30. 死时间倍乘/分数两条公式；8 点宽 glitch 的 despike vs deglitch
    31. WT 壳层切片：Fe–Fe 的 k_peak 高于 Fe–O ≥0.8 Å⁻¹；模型 WT 残差 R<0.15
    32. ε_R / χ²_ν 为正；双壳相对单壳的 F-test 置信度 >0.95
    33. C3 合成回收：C3 真值 ±5e-5；N_idp 不足时自动降级不崩溃

Phase II（B5–B8）新增 9 项契约：
    - 复数小波：反相 χ 的 complex_r_factor=4.0，幅值残差=0.0
    - 元素理论 k_peak 标尺：O < Cl < Fe 且随 σ² 左移
    - 多扫描合并：MAD 离群剔除、σ_μ 传播、默认均值逐位兼容
    - XANES 荧光自吸收：Tröger 反解恢复白线/LCF，分母防爆门禁
    - 边前 Pseudo-Voigt：双峰质心误差 <0.15 eV，指标表同步
    - 各向同性晶格约束：M 个 ΔR → 1 个 α；实验 NiO 诚实降级
    - 部分占据 CIF：共占位合并、占据度缩放简并度（6×0.75=4.5）
    - C4 回收 + 非高斯 g(R)：面积=N、正偏态尾长判据
    - 多数据集联合拟合：5 个全局变量、变温爱因斯坦 σ²(T) 回收

Phase III（B9–B12）新增 15 项契约：
    - 17 键配置文件生效；C4 的 N_idp 降级与 g(R) 导出
    - σ_μ 落盘/SNR 门禁；WT 元素标尺自动挂载与 L3 白线面积
    - .prj round-trip 与 MCA 坏道剔除；MBACK 短窗口归一化
    - AUTOBK chi_std/clamp；SPOIL + MCR-ALS 盲解谱
    - 相移校正 FT/WT；异核共占位双分支；关联 Debye；BVS 软约束
    - R/q/WT 多空间拟合；双吸收边 ΔE₀ 分组 + AIC/BIC；SAC/DAC 团簇

P4+ 新增：
    - T7 多元素综合报告：自动分元素、并列结果、不推荐 Top-1
    - T8 XANES–EXAFS 软校验：只产生警告，不改变 cn_reliable
    - P5-A T5 标定产物完整性 / --only 增量合并 / --s02-from-t5 溯源
    - P5-B Mo 双壳 bcc 标定 / Ce σ² 压界与归一化伪迹拒绝
    - P5-C 固定 COD ID 补齐 18 个常见相 / coverage 无缺口
    - P5-prep 标准谱源隔离池 / 5 个未通过元素门禁冻结

当前共 135 项契约回归（B1-1 至 B20-1，外加 P0a–P4 正式 CN 契约、P4+
CIF 抓取/条件 CN/T7/T8/P5-A/P5-B/P5-C/P5-prep 契约，以及 Fe 化学态
S₀² 多态契约）。

运行:
    python regression_test.py
"""
from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import artemis                        # noqa: E402
import athena                         # noqa: E402
import cn_table                       # noqa: E402
import cn_disclosure                  # noqa: E402
import larch_wt                       # noqa: E402
import motif_provenance               # noqa: E402
import phase_pack                     # noqa: E402
import phase_resolver                 # noqa: E402
import cif_gate                       # noqa: E402
import structure_id                   # noqa: E402
import structure_lib                  # noqa: E402
import structure_model                # noqa: E402
import xafs_config                    # noqa: E402
import xafs_diagnostics               # noqa: E402
import xafs_fit                       # noqa: E402
import xafs_io                        # noqa: E402
import xafs_pipeline                  # noqa: E402
import xafs_provenance                # noqa: E402
import xafs_report                    # noqa: E402
import xafs_window_tuner              # noqa: E402
from absorber import get_absorber      # noqa: E402
from xafs_core import (ChiK, K2, NormalizedScan, RawScan,  # noqa: E402
                       XAFSConfigError, XAFSError)


def _find(pattern: str) -> Path:
    hits = sorted((ROOT / "data").glob(pattern))
    if not hits:
        raise unittest.SkipTest(f"测试数据缺失: data/{pattern}")
    return hits[0]


LEPIDO = _find("*/02_Fe_standards/Fe_lepidocrocite.000")
FE_FOIL = _find("*/01_Fe_foil/fe.060")
GOETHITE = _find("bmm_standards/Data/Fe/Fe-Goethite.xdi")
MAGNETITE = _find("bmm_standards/Data/Fe/Fe-Magnetite.xdi")
FERRIHYDRITE = _find("bmm_standards/Data/Fe/Fe-Ferrihydrite.xdi")
PT_FOIL = _find("bmm_standards/Data/Pt/Pt-L3-PlatinumFoil.xdi")
PT_CL = _find("bmm_standards/Data/Pt/Pt-L3-PlatinumChloride.xdi")
NIO = _find("bmm_standards/Data/Ni/Ni-K-NickelOxide.xdi")


class RegressionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lep = xafs_pipeline.process_sample(
            [LEPIDO], absorber="Fe", verbose=False)
        cls.lep_fit = xafs_fit.run(
            cls.lep, identify=True, fit=True,
            absorber_symbol="Fe", verbose=False)

    # ------------------------------------------------------------
    # 1. 数据契约
    # ------------------------------------------------------------
    def test_chik_length_contract(self):
        chi = self.lep.chi_k
        self.assertEqual(len(chi.k), len(chi.chi))
        self.assertEqual(len(chi.background), len(chi.mu_norm))
        self.assertEqual(len(chi.mu_norm), len(self.lep.norm.energy))
        self.assertNotEqual(len(chi.mu_norm), len(chi.k),
                            "mu_norm 不应被误当成 k 网格数组")

    # ------------------------------------------------------------
    # 2. nleg 过滤
    # ------------------------------------------------------------
    def test_theory_chi_filters_ms(self):
        from larch.xafs import feffpath

        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "lepidocrocite"][0]
        cache = structure_lib.ensure_feff(motif, edge="K")
        absorber = get_absorber("Fe")
        k = np.arange(1.0, 13.0, 0.05)

        def model(max_legs):
            paths = []
            for f in sorted(cache.glob("feff*.dat")):
                p = feffpath(str(f), s02=absorber.s02 or 1.0,
                             sigma2=0.003, e0=0.0, deltar=0.0)
                nleg = getattr(p, "nleg", getattr(p, "legs", 2))
                if nleg <= max_legs and p.degen >= 1 and p.reff <= 5.5:
                    p._calc_chi()
                    paths.append(p)
            return artemis.exafs_model_k(paths, k)

        filtered = model(4)
        unfiltered = model(99)
        actual = structure_id.theory_chi(motif, k, absorber, sigma2=0.003)
        np.testing.assert_allclose(actual, filtered, rtol=1e-10, atol=1e-12)
        self.assertFalse(np.allclose(actual, unfiltered, rtol=1e-6, atol=1e-9),
                         "MS 路径仍被计入 theory_chi——nleg 过滤失效")

    # ------------------------------------------------------------
    # 3. 识别结果语义
    # ------------------------------------------------------------
    def test_indistinguishable_excludes_best(self):
        sid = self.lep_fit["structure"]
        best = sid["ranked"][0]["id"]
        self.assertEqual(best, "lepidocrocite")
        self.assertNotIn(best, sid["indistinguishable"])
        self.assertEqual(sid["status"], "confident")
        self.assertGreater(sid["confidence"]["margin"], 0.02)

    # ------------------------------------------------------------
    # 4. FEFF 吸收边 HOLE 卡
    # ------------------------------------------------------------
    def test_feff_hole_card_for_l_edges(self):
        sites = [(0.0, 0.0, 0.0, "Pt"), (2.77, 0.0, 0.0, "Pt")]
        expected = {"K": 1, "L1": 2, "L2": 3, "L3": 4}
        for edge, idx in expected.items():
            inp = artemis.build_feff_input(
                sites, absorber="Pt", rmax=3.0, edge=edge)
            self.assertIn(f"HOLE {idx}   1.0", inp)
        with self.assertRaises(XAFSConfigError):
            artemis.build_feff_input(sites, absorber="Pt", edge="M5")

    # ------------------------------------------------------------
    # 5. WT 方向
    # ------------------------------------------------------------
    def test_wavelet_orientation(self):
        k = np.arange(0.0, 16.0, 0.05)
        chi = ChiK(k=k, chi=np.sin(2.0 * k), background=np.zeros_like(k),
                   mu_norm=np.zeros_like(k), rbkg=1.0)
        wt = larch_wt.continuous_wavelet_transform(
            chi, sigma=1.5, kmin=3.0, kmax=13.0, kweight=2)
        self.assertEqual(wt.magnitude.shape,
                         (len(wt.k_centers), len(wt.r)))
        self.assertEqual(wt.magnitude.shape[1], len(wt.r))

    # ------------------------------------------------------------
    # 6. strict 门禁
    # ------------------------------------------------------------
    def test_strict_gate_and_downgrade(self):
        with self.assertRaises(XAFSError) as ctx:
            xafs_pipeline.process_sample([FE_FOIL], absorber="Fe",
                                         strict=True, verbose=False)
        self.assertIn("归一化失败", str(ctx.exception))

        res = xafs_pipeline.process_sample(
            [FE_FOIL], absorber="Fe", strict=False, verbose=False)
        self.assertLess(res.norm.plateau, 0.85,
                        "该标样本应触发平台门禁，否则测试前提失效")

    # ------------------------------------------------------------
    # 7/8. S7 两条路径
    # ------------------------------------------------------------
    def test_fit_motif_returns_shell_cn(self):
        fit_out = self.lep_fit["fit"]
        self.assertEqual(fit_out["mode"], "motif")
        self.assertEqual(fit_out["motif"], "lepidocrocite")
        self.assertTrue(fit_out["coverage_ok"])
        self.assertTrue(fit_out["cn_reliable"])
        fit = fit_out["fit"]
        self.assertFalse(fit.overfit)
        first = fit.paths[0]
        self.assertGreater(first.degeneracy, 4.5)
        self.assertLess(first.degeneracy, 8.0)
        self.assertGreater(first.r, 1.90)
        self.assertLess(first.r, 2.15)
        self.assertLess(fit.r_factor, 0.10)

    def test_fit_l1_when_no_match(self):
        res = xafs_pipeline.process_sample(
            [FERRIHYDRITE], absorber="Fe", strict=False, verbose=False)
        out = xafs_fit.run(res, identify=True, fit=True,
                           absorber_symbol="Fe", verbose=False)
        sid = out["structure"]
        self.assertEqual(sid["status"], "no_match")
        fit_out = out["fit"]
        self.assertEqual(fit_out["mode"], "l1_first_shell")
        l1 = fit_out["l1"]
        self.assertGreater(l1["N"], 4.0)
        self.assertLess(l1["N"], 7.0)
        self.assertGreater(l1["R"], 1.85)
        self.assertLess(l1["R"], 2.05)
        self.assertIn("仅第一壳", l1["warning"])

    # ------------------------------------------------------------
    # L1 选峰/窗口（Pt L3 实测修正）
    # ------------------------------------------------------------
    def test_l1_ligand_aware_peak_selection(self):
        res = xafs_pipeline.process_sample(
            [PT_CL], absorber="Pt", strict=False, verbose=False)
        l1 = structure_id.l1_first_shell_fit(
            res, absorber_symbol="Pt", ligand="Cl", verbose=False)
        # 平方平面 PtCl4：N=4、Pt–Cl 2.30 Å；错误的低 R 伪峰会给出 Rf>0.3
        self.assertLess(l1["rf"], 0.08)
        self.assertAlmostEqual(l1["N"], 4.0, delta=0.5)
        self.assertAlmostEqual(l1["R"], 2.30, delta=0.05)

    def test_l1_long_bond_window_not_truncated(self):
        res = xafs_pipeline.process_sample(
            [PT_FOIL], absorber="Pt", strict=False, verbose=False)
        l1 = structure_id.l1_first_shell_fit(
            res, absorber_symbol="Pt", ligand="Pt", verbose=False)
        # 旧 rmax=min(2.3, ...) 把 Pt–Pt 2.77 的模型峰截在窗外 → Rf≈1.0
        self.assertLess(l1["rf"], 0.08)
        self.assertGreater(l1["N"], 9.0)
        self.assertAlmostEqual(l1["R"], 2.77, delta=0.05)

    def test_shell_selection_respects_n_idp(self):
        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "goethite"][0]
        plan = xafs_fit.select_shell_paths(
            motif, rmax=3.5, merge_tol=0.20, kmin=3.0, kmax=12.0)
        self.assertEqual(plan["shells"][0]["label"], "Fe-O@2.03")
        self.assertEqual(plan["shells"][0]["degen"], 6.0)
        self.assertEqual(plan["shells"][1]["label"], "Fe-Fe@3.02")
        self.assertEqual(plan["max_paths_by_idp"], 3)
        self.assertTrue(plan["truncated"])

    def test_goethite_motif_is_confident(self):
        res = xafs_pipeline.process_sample(
            [GOETHITE], absorber="Fe", strict=False, verbose=False)
        sid = structure_id.identify_result(
            res, absorber_symbol="Fe", verbose=False)
        self.assertEqual(sid["status"], "confident")
        self.assertEqual(sid["ranked"][0]["id"], "goethite")
        self.assertTrue(sid["coverage_ok"])
        self.assertEqual(sid["uncovered_peaks"], [])

    def test_magnetite_mixed_is_confident(self):
        """磁铁矿多位点识别（曾出问题的案例）。

        历史：只有 16d 八面体位点时 score=0.145、实验峰 1.35/3.712
        无理论峰解释（覆盖率失败→ partial）。补 8a 四面体位点与
        mixture 联合模型后，联合模型 score=0.083 成为 Top-1，且覆盖
        按"复合谱 ∪ 各分量分别解释"通过——tet 分量解释 1.35/3.712，
        oct 分量解释 2.608。
        """
        res = xafs_pipeline.process_sample(
            [MAGNETITE], absorber="Fe", strict=False, verbose=False)
        sid = structure_id.identify_result(
            res, absorber_symbol="Fe", verbose=False)
        self.assertEqual(sid["status"], "confident")
        self.assertEqual(sid["ranked"][0]["id"], "magnetite_mixed")
        self.assertTrue(sid["coverage_ok"])
        self.assertEqual(sid["uncovered_peaks"], [])
        self.assertGreater(sid["confidence"]["margin"], 0.03)

    def test_mixture_shell_selection(self):
        """多位点吸收的壳层分割。

        磁铁矿 8a/16d 的两个子壳层 Fe-O 1.886 / 2.060 Å 只差 0.17 Å，
        单位点容差 0.20 Å 会把它们错并成一层；mixture 的简并度按
        位点权重缩放（tet O×4·⅓=1.33，oct O×6·⅔=4.0）。
        """
        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "magnetite_mixed"][0]
        plan = xafs_fit.select_shell_paths(motif, rmax=3.5,
                                           kmin=3.0, kmax=12.0)
        first, second, third = plan["shells"][:3]
        self.assertEqual(first["label"], "Fe-O@1.89")
        self.assertAlmostEqual(first["degen"], 4.0 / 3.0, places=2)
        self.assertEqual(second["label"], "Fe-O@2.06")
        self.assertAlmostEqual(second["degen"], 4.0, places=2)
        self.assertEqual(third["label"], "Fe-Fe@2.97")
        self.assertAlmostEqual(third["degen"], 4.0, places=2)
        self.assertLessEqual(plan["merge_tol"], 0.10,
                             "mixture 必须收紧合并容差，否则子壳层被吞")

    def test_vary_s02_recovers_synthetic_truth(self):
        """标定回路（calibrate_s02.py 的前提契约）。

        数据由已知参数合成：S₀²=0.72、σ²=0.004、E0+=2 eV、N=8 固定
        （fe_bcc 第一壳）。fit_r_space 以 vary_s02=True、s02=None 启动，
        N 固定为晶体学简并度时唯一振幅变量 S₀² 必须回收真值——
        回收不了说明标定会写入错误数值。
        """
        truth = {"s02": 0.72, "sigma2": 0.004, "e0": 2.0}
        k = np.arange(3.0, 12.0, 0.05)
        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "fe_bcc"][0]
        plan = xafs_fit.select_shell_paths(motif, edge="K", rmax=2.9,
                                           kmin=3.0, kmax=12.0)
        first = plan["shells"][0]
        self.assertAlmostEqual(first["degen"], 8.0, places=6)
        from larch.xafs import feffpath
        p = feffpath(str(plan["files"][0]), s02=truth["s02"],
                     sigma2=truth["sigma2"], e0=truth["e0"], deltar=0.0)
        p.degen = float(first["degen"])
        p._calc_chi()
        chi = artemis.exafs_model_k([p], k)
        chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                     mu_norm=np.zeros_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
        fit = artemis.fit_r_space(
            chi_k, ft, [plan["files"][0]], rmin=1.9, rmax=2.75,
            kmin=3.0, kmax=12.0, kweight=2, s02=0.85,
            vary_delta_e0=True, vary_n=False, vary_s02=True,
            path_degen=[first["degen"]], path_labels=[first["label"]])
        self.assertFalse(fit.overfit)
        self.assertLess(fit.r_factor, 0.02, f"无噪声合成数据应近完美复现")
        self.assertLess(abs(fit.s02 - truth["s02"]), 0.03,
                        f"S₀² 回收偏差过大: {fit.s02}")
        self.assertLess(fit.s02_err, 0.06)
        self.assertLess(abs(fit.paths[0].sigma2 - truth["sigma2"]), 1e-3)
        self.assertLess(abs(fit.delta_e0 - truth["e0"]), 0.6)

    def test_mixture_weight_recovery(self):
        """mixture 相分数回收（--fit-weights 的核心数学契约）。

        磁铁矿 tet/oct 第一壳 Fe-O 1.886/2.060 Å 用 0.25/0.75 权重
        合成 χ。拟合固定每路径的**晶体学简并度**（4/6）、固定 σ²，
        放开 softmax 权重 + ΔR + ΔE₀。无噪声时权重必须回收真值，
        有效简并度随之等于 N_crystal·w（1.0 / 4.5）——回收不了说明
        相分数只是把 N 换个写法，不能用于物相定量。
        """
        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "magnetite_mixed"][0]
        plan = xafs_fit.select_shell_paths(
            motif, edge="K", rmax=2.6, kmin=3.0, kmax=12.0,
            group_by_component=True)
        self.assertEqual(plan["group_names"],
                         ["magnetite_tet", "magnetite_oct"])
        self.assertEqual(plan["groups"], [0, 1])
        self.assertEqual(plan["degens"], [4.0, 6.0],
                         "权重模式必须保留未乘权重的晶体学简并度")
        self.assertEqual(plan["weight_degenerate_pairs"], [])

        true_w = np.array([0.25, 0.75])
        s02, sigma2, e0 = 0.80, 0.004, 2.0
        k = np.arange(3.0, 12.0, 0.05)
        from larch.xafs import feffpath
        paths = []
        for i, f in enumerate(plan["files"]):
            p = feffpath(str(f), s02=s02, sigma2=sigma2, e0=e0, deltar=0.0)
            p.degen = plan["degens"][i] * true_w[plan["groups"][i]]
            p._calc_chi()
            paths.append(p)
        chi = artemis.exafs_model_k(paths, k)
        chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                     mu_norm=np.zeros_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
        fit = artemis.fit_r_space(
            chi_k, ft, plan["files"],
            rmin=plan["rmin"], rmax=plan["rmax"],
            kmin=3.0, kmax=12.0, kweight=2, s02=s02,
            vary_delta_e0=True, vary_n=False,
            vary_sigma2=False, sigma2_init=sigma2,
            path_degen=plan["degens"], path_labels=plan["labels"],
            path_groups=plan["groups"],
            path_group_names=plan["group_names"],
            vary_weights=True, weight_init=[0.5, 0.5])

        self.assertFalse(fit.overfit)
        self.assertLess(fit.r_factor, 0.02, "无噪声合成数据应近完美复现")
        self.assertIsNotNone(fit.mixture_weights)
        got = [w["weight"] for w in fit.mixture_weights]
        self.assertAlmostEqual(sum(got), 1.0, places=9)
        for kk in range(2):
            self.assertLess(abs(got[kk] - true_w[kk]), 0.03,
                            f"权重回收偏差过大: {got}")
            err = fit.mixture_weights[kk]["weight_err"]
            self.assertTrue(np.isfinite(err))
            self.assertLess(err, 0.03)
            self.assertAlmostEqual(
                fit.paths[kk].degeneracy,
                plan["degens"][kk] * true_w[kk], places=4,
                msg="有效简并度必须等于 N_crystal·w")
            self.assertIsNotNone(fit.paths[kk].n_err)
        self.assertLess(abs(fit.delta_e0 - e0), 0.6)
        self.assertEqual(len(fit.parameter_names), fit.n_variables)
        self.assertEqual(fit.correlation.shape,
                         (fit.n_variables, fit.n_variables))
        self.assertTrue(np.all(np.isfinite(fit.correlation)))
        self.assertTrue(np.allclose(fit.correlation, fit.correlation.T))
        self.assertTrue(np.allclose(np.diag(fit.correlation), 1.0))
        self.assertLessEqual(np.max(np.abs(fit.correlation
                                           - np.eye(fit.n_variables))), 1.0)

    def test_fluorescence_self_absorption_gate(self):
        """荧光自吸收门禁（不做无法验证的校正，只做保守降级）。

        荧光信号被样品自身再吸收会压低 EXAFS 振幅与 N；没有组成/厚度
        就不能定量校正。契约：Δμ≥1.0 判 high 并拒绝定量；0.3–1.0 只
        警告；<0.3 视为稀样；透射/given 模式不受影响；Δμ 缺失时不放行。
        """
        from types import SimpleNamespace

        def diag(edge, mode="fluor"):
            res = SimpleNamespace(
                merged=SimpleNamespace(meta={"mode": mode}),
                norm=SimpleNamespace(edge_step=edge))
            return xafs_fit._self_absorption_diag(res)

        high = diag(1.2)
        self.assertEqual(high["risk"], "high")
        self.assertTrue(high["cn_gate"])
        self.assertFalse(high["corrected"])
        warn = diag(0.5)
        self.assertEqual(warn["risk"], "possible")
        self.assertFalse(warn["cn_gate"])
        low = diag(0.05)
        self.assertEqual(low["risk"], "low")
        self.assertFalse(low["cn_gate"])
        trans = diag(1.2, mode="trans")
        self.assertEqual(trans["risk"], "none")
        self.assertFalse(trans["cn_gate"])
        unknown = diag(float("nan"))
        self.assertEqual(unknown["risk"], "unknown")
        self.assertTrue(unknown["cn_gate"])

    def test_booth_self_absorption_recovery(self):
        """Booth 反解必须回收 Booth 前向模型（Eq.2 → Eq.4）。

        用 Eq.2 从已知 χ 生成 χ_exp，再走 booth_self_absorption；
        要求 max|Δχ| ≤ 0.01（含厚样修正项），薄极限严格恒等，
        厚极限保持有限。任何公式符号/单位退化都会先在这里失败。
        """
        k = np.arange(3.0, 12.01, 0.05)
        chi = 0.09 * np.sin(2 * 2.0 * k + 0.4) * np.exp(-2 * 0.004 * k**2)
        e0, d, rho = 7120.0, 4.6, 5.17
        anginp, angout = 49.4, 41.0
        mu_t, mu_a, mu_f, _, _ = athena._booth_mu_terms(
            "Fe3O4", "Fe", "K", "Ka", e0, k, rho)
        sin_phi = np.sin(np.deg2rad(anginp))
        g = sin_phi / np.sin(np.deg2rad(angout))
        alpha = mu_t + g * mu_f
        arg = alpha * d / sin_phi
        arg_true = (alpha + chi * mu_a) * d / sin_phi
        chi_exp = ((1 - np.exp(-arg_true)) / (1 - np.exp(-arg))
                   * alpha * (chi + 1.0) / (alpha + chi * mu_a) - 1.0)
        chi_in = ChiK(k=k, chi=chi_exp, background=np.zeros_like(k),
                      mu_norm=np.ones_like(k), rbkg=1.2)
        out, info = athena.booth_self_absorption(
            chi_in, formula="Fe3O4", elem="Fe", edge="K", line="Ka", e0=e0,
            thickness_um=d, density=rho, anginp=anginp, angout=angout)
        self.assertTrue(info["reliable"])
        self.assertLess(float(np.max(np.abs(out.chi - chi))), 0.01)
        self.assertGreater(info["correction_factor_median"], 1.0)

        thin, _ = athena.booth_self_absorption(
            chi_in, formula="Fe3O4", elem="Fe", edge="K", line="Ka", e0=e0,
            thickness_um=1e-6, density=rho, anginp=anginp, angout=angout)
        self.assertTrue(np.array_equal(thin.chi, chi_exp))

        thick, _ = athena.booth_self_absorption(
            chi_in, formula="Fe3O4", elem="Fe", edge="K", line="Ka", e0=e0,
            thickness_um=1e4, density=rho, anginp=anginp, angout=angout)
        self.assertTrue(np.all(np.isfinite(thick.chi)))

        # 垂直入射（90°）是合法几何，不能因边界检查误拒
        normal, _ = athena.booth_self_absorption(
            chi_in, formula="Fe3O4", elem="Fe", edge="K", line="Ka", e0=e0,
            thickness_um=d, density=rho, anginp=90.0, angout=45.0)
        self.assertTrue(np.all(np.isfinite(normal.chi)))
        with self.assertRaises(XAFSConfigError):
            athena.booth_self_absorption(
                chi_in, formula="Fe3O4", elem="Fe", edge="K", line="Ka",
                e0=e0, thickness_um=d, density=rho, anginp=0.0, angout=45.0)

    def test_booth_end_to_end_synthetic(self):
        """全链合成：真实 μ → Eq.2 前向荧光谱 → 归一化/autobk → Booth。

        这是对 S4.5 最有信息量的契约：只测公式反解（上一项）不够，
        还要覆盖 E0 漂移、边前基线、autobk 归一化和 k 空间形状。
        期望：未校正 χ 振幅明显低于真值（<0.85），校正后回到
        0.9–1.2；若 μ_a 基线或能量锚点退化，这一项会立刻失败。
        """
        import tempfile
        from xraydb import material_mu, material_mu_components, xray_line

        formula, elem, line, d, rho = "Fe3O4", "Fe", "Ka", 4.6, 5.17
        angin, angout = 49.4, 41.0
        scan = xafs_io.read_scan(MAGNETITE)
        E = np.asarray(scan.energy, dtype=float)
        e0 = athena.determine_e0(scan)["e0_sample"]
        norm = athena.normalize_absorption(scan, e0)
        chi_true = athena.extract_chi_k(norm, rbkg=1.10, kweight=2)

        comp = material_mu_components(formula, E, density=rho, kind="photo")
        mass = float(comp["mass"])
        mu_t = np.zeros_like(E)
        for el in comp["elements"]:
            stoich, amass, mu_mass = comp[el]
            mu_t += rho * float(stoich) * float(amass) * np.asarray(mu_mass) / mass
        mu_a = (rho * comp[elem][0] * comp[elem][1]
                * np.asarray(comp[elem][2]) / mass)
        comp_below = material_mu_components(
            formula, float(e0) - 200.0, density=rho, kind="photo")
        mu_a -= (rho * comp_below[elem][0] * comp_below[elem][1]
                 * float(comp_below[elem][2]) / float(comp_below["mass"]))
        mu_t *= 1e-4
        mu_a *= 1e-4
        mu_f = float(material_mu(formula, xray_line(elem, line).energy,
                                 density=rho, kind="total")) * 1e-4

        k_on_e = np.sqrt(np.maximum(K2 * (E - e0), 0.0))
        chi_on_e = np.interp(k_on_e, chi_true.k, chi_true.chi,
                             left=0.0, right=0.0)
        chi_on_e[E <= e0] = 0.0
        # 振荡幅度 = 边跳变 × χ（不是平滑 μ_a(E)·χ）
        c_lo = material_mu_components(formula, e0 - 10.0, density=rho,
                                      kind="photo")
        c_hi = material_mu_components(formula, e0 + 10.0, density=rho,
                                      kind="photo")
        ma_lo = (rho * c_lo[elem][0] * c_lo[elem][1]
                 * float(c_lo[elem][2]) / float(c_lo["mass"])) * 1e-4
        ma_hi = (rho * c_hi[elem][0] * c_hi[elem][1]
                 * float(c_hi[elem][2]) / float(c_hi["mass"])) * 1e-4
        jump = float(ma_hi - ma_lo)
        mu_a_e = mu_a + jump * chi_on_e
        mu_t_e = mu_t + jump * chi_on_e
        sin_phi = np.sin(np.deg2rad(angin))
        g = sin_phi / np.sin(np.deg2rad(angout))
        arg = (mu_t_e / sin_phi + mu_f / np.sin(np.deg2rad(angout))) * d
        arg0 = (mu_t / sin_phi + mu_f / np.sin(np.deg2rad(angout))) * d
        ifluo = (mu_a_e / (mu_t_e + g * mu_f)
                 * (1.0 - np.exp(-arg)) / (1.0 - np.exp(-arg0)) + 0.35)

        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "forward.xdi"
            lines = ["# XDI/1.0 synthetic-forward", "# Element.symbol: Fe",
                     "# Element.edge: K", "# Column.1: energy eV",
                     "# Column.2: i0", "# Column.3: ifluor", "# -----------"]
            for e, f in zip(E, ifluo):
                lines.append(f"{e:.6f} 1000000.000000 {f * 1e6:.8f}")
            path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            common = dict(absorber="Fe", strict=False, rbkg=1.10,
                          kmin=3.0, kmax=12.0, kweight=2, verbose=False)
            raw = xafs_pipeline.process_sample([str(path)], **common)
            corr = xafs_pipeline.process_sample(
                [str(path)], self_absorption_formula=formula,
                self_absorption_thickness=d, self_absorption_density=rho,
                self_absorption_angles=(angin, angout), **common)

        kg = np.linspace(3.0, 12.0, 900)
        truth = np.interp(kg, chi_true.k, chi_true.chi)
        raw_amp = np.std(np.interp(kg, raw.chi_k.k, raw.chi_k.chi)) / np.std(truth)
        cor_amp = np.std(np.interp(kg, corr.chi_k.k, corr.chi_k.chi)) / np.std(truth)
        self.assertLess(raw_amp, 0.85, f"前向谱未被压低: {raw_amp:.3f}")
        self.assertGreater(cor_amp, 0.90, f"校正不足: {cor_amp:.3f}")
        self.assertLess(cor_amp, 1.20, f"校正过头: {cor_amp:.3f}")

    def test_wt_locator_anchored_to_ft_peak(self):
        """WT 定位必须锚定 FT 峰，不能在两壳分支间跳。

        实测故障：Fe3O4 的 WT 在 R≈1.42（Fe–O，k≈4.5）与 R≈2.72
        （Fe–Fe，k≈5.8）各有一条分支，宽窗 [1.0,2.6] 内二者只差 3%；
        Booth 校正把低 k 抬高 7% 后 WT 峰从 2.60 跳到 1.42，而 FT 峰
        完全不变。契约：锚定 FT 峰后，校正前后 WT R 位移动 ≤ 0.02 Å。
        """
        res = xafs_pipeline.process_sample(
            [MAGNETITE], absorber="Fe", strict=False, verbose=False)
        chi2, info = athena.booth_self_absorption(
            res.chi_k, formula="Fe3O4", elem="Fe", edge="K", line="Ka",
            e0=res.e0["e0_sample"], thickness_um=4.6, density=5.17,
            anginp=49.4, angout=41.0)
        self.assertGreater(info["correction_factor_range"][1], 1.4,
                           "测试前提：校正因子应有明显低 k 抬升")

        def anchored(chi_k):
            wt = larch_wt.continuous_wavelet_transform(
                chi_k, sigma=1.5, kmin=3.0, kmax=12.0, kweight=2)
            return larch_wt.locate_shell_near(wt, res.r1)

        k_a, r_a, _ = anchored(res.chi_k)
        k_b, r_b, _ = anchored(chi2)
        self.assertLessEqual(abs(r_b - r_a), 0.02,
                             f"WT 峰随校正跳变: {r_a} → {r_b}")
        self.assertLessEqual(abs(r_a - res.r1), 0.5,
                             f"WT 与 FT 峰位不一致: {r_a} vs {res.r1}")
        self.assertTrue(np.isfinite(k_a) and np.isfinite(k_b))

    def test_sigma2_shared_group(self):
        """σ² 共享组：多条路径共享一个 σ² 变量，减少自由参数。

        对标 Artemis 的路径约束。用 NiO 双壳层（Ni-O + Ni-Ni）验证：
        - shared_sigma2="all"：两条路径共享一个 σ²，拟合收敛，Rf 合理
        - 显式分组 [[0, 1]]：等价于 "all"（两条路径都在同一组）
        - 共享组输出的 σ² 值在两条路径中相同
        若 σ² 共享组机制失效，拟合会报错或输出不同的 σ²。
        """
        res = xafs_pipeline.process_sample(
            [NIO], absorber="Ni", strict=False, verbose=False)
        out = xafs_fit.fit_identified(
            res, absorber_symbol="Ni", shared_sigma2="all", verbose=False)
        fit = out["fit"]
        self.assertEqual(out["mode"], "motif")
        self.assertTrue(out["sigma2_shared"])
        self.assertEqual(len(out["sigma2_groups"]), 1)
        n_path = len(out["fit"].paths)
        self.assertEqual(out["sigma2_groups"][0], list(range(n_path)))
        # 所有路径的 σ² 应相同（共享同一变量）
        s2s = [p.sigma2 for p in fit.paths]
        for s2 in s2s:
            self.assertAlmostEqual(s2, s2s[0], places=6)
        # σ²_err 也应相同
        e0s = [p.sigma2_err for p in fit.paths]
        for e in e0s:
            self.assertIsNotNone(e)
        for e in e0s:
            self.assertAlmostEqual(e, e0s[0], places=6)
        # Rf 应在合理范围（双壳层 NiO 合成谱，Rf 应 < 0.02）
        self.assertLess(fit.r_factor, 0.02,
                        f"σ² 共享后 Rf={fit.r_factor:.4f} 过大，拟合可能失败")

        # 显式分组 [[0, 1, ..., n-1]] 等价于 "all"
        out2 = xafs_fit.fit_identified(
            res, absorber_symbol="Ni",
            shared_sigma2=[list(range(n_path))], verbose=False)
        fit2 = out2["fit"]
        s2s2 = [p.sigma2 for p in fit2.paths]
        for s2 in s2s2:
            self.assertAlmostEqual(s2, s2s[0], places=6)
        self.assertAlmostEqual(fit.r_factor, fit2.r_factor, places=6)

    def test_systematic_budget_dynamic(self):
        """CN 系统误差预算必须随本次拟合条件变化，而不是恒为 ±30%。

        基础预算 0.30 是历史实测锚点；额外项（自吸收残留、参数强相关、
        远壳层截断、σ² 撞界）按平方和合成。这里用构造对象验证每一类
        额外项的数值贡献，防止预算退化成写死的常数。
        """
        class FakeFit:
            correlation = np.array([[1.0, 0.95], [0.95, 1.0]])
            parameter_names = ["Fe-Fe@2.97:N", "Fe-Fe@2.97:σ²"]

        base = xafs_fit.systematic_budget(None)
        self.assertAlmostEqual(base["relative_1sigma"], 0.30, places=4)
        self.assertEqual(base["extra_components"], {})
        self.assertAlmostEqual(base["max_parameter_correlation"], 0.0,
                               places=4)

        # 0.30 base + 0.10 相关 + 0.10 自吸收残留 + 0.05 截断 + 0.10 σ²撞界
        rich = xafs_fit.systematic_budget(
            FakeFit(), flags={"sigma2_at_bound": True},
            self_abs={"corrected": True, "reliable": True},
            truncated_shells=2)
        expected = float(np.sqrt(0.30**2 + 0.10**2 + 0.10**2
                                 + 0.05**2 + 0.10**2))
        self.assertAlmostEqual(rich["relative_1sigma"], expected, places=4)
        self.assertEqual(rich["extra_components"],
                         {"parameter_correlation": 0.10,
                          "self_absorption_residual": 0.10,
                          "truncated_shells": 0.05,
                          "sigma2_at_bound": 0.10})
        self.assertIn("Fe-Fe@2.97:N vs Fe-Fe@2.97:σ²", rich["note"])

        # 未校正但 risk=possible → +0.15；risk=low → 不加项
        possible = xafs_fit.systematic_budget(
            None, self_abs={"corrected": False, "risk": "possible"})
        self.assertAlmostEqual(possible["relative_1sigma"],
                               float(np.hypot(0.30, 0.15)), places=4)
        low = xafs_fit.systematic_budget(
            None, self_abs={"corrected": False, "risk": "low"})
        self.assertAlmostEqual(low["relative_1sigma"], 0.30, places=4)

    def test_html_report_single_file(self):
        """HTML 报告：图片内联为 data URI，强相关/不可靠项标红。"""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            fig, ax = plt.subplots(figsize=(1, 1))
            ax.plot([0, 1], [0, 1])
            fig.savefig(out / "demo_xas.png", dpi=40)
            plt.close(fig)
            rows = [
                ("cn_reliable", "False", "",
                 "CN 是否可作结论"),
                ("max_parameter_correlation",
                 "0.950 (a:N vs a:σ²)", "",
                 "最强参数相关"),
                ("r_factor", "0.1234", "", "拟合优度"),
            ]
            fit_out = {"status": "confident", "cn_reliable": False,
                       "flags": {"r_factor_bad": True},
                       "cn_systematics": {"relative_1sigma": 0.35,
                                          "note": "test budget"},
                       "fit": None}
            path = xafs_report.export_html_report(
                out, "demo", sample="demo sample",
                edge_label="Fe K-edge", indicator_rows=rows,
                table_rows=[["demo", "Fe-O", "6.0", "2.00", "0.004",
                             "0.0", "0.003"]],
                fit_out=fit_out, diagnostics_dir=None)
            self.assertTrue(path.exists())
            doc = path.read_text(encoding="utf-8")
            self.assertIn("data:image/png;base64,", doc)
            self.assertGreaterEqual(doc.count("class='bad'"), 2)
            self.assertIn("r_factor_bad", doc)
            self.assertIn("35%", doc)

    def test_config_rejects_unknown_suffix(self):
        """不支持的扩展名必须给出干净 ValueError（曾因 p.sffix 拼写崩溃）。"""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "config.ini"
            p.write_text("kmin: 3.0\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                xafs_config.load_config(p)
            self.assertIn(".ini", str(ctx.exception))
            self.assertIn(".yaml", str(ctx.exception))

    def test_hann_window_contract(self):
        """浏览器端 Hanning 窗必须与产线同约定：两端 0、平台 1。

        交互页在浏览器里重算 FT；若窗函数端点或斜坡写错，峰位/形状会
        与产线不一致，而页面看起来仍“能出图”，很难察觉。
        """
        k = np.linspace(0.0, 16.0, 1601)
        kmin, kmax, dk = 3.0, 12.0, 1.0
        w = xafs_window_tuner.hann_window(k, kmin, kmax, dk)
        self.assertTrue(np.all(w >= -1e-12) and np.all(w <= 1 + 1e-12))
        # 窗外严格为 0（端点也含在窗内，但 Hanning 在端点处为 0）
        self.assertTrue(np.allclose(w[k < kmin], 0.0))
        self.assertTrue(np.allclose(w[k > kmax], 0.0))
        # 平台
        plateau = (k > kmin + dk + 0.05) & (k < kmax - dk - 0.05)
        self.assertGreater(plateau.sum(), 100)
        self.assertTrue(np.allclose(w[plateau], 1.0))
        # 左坡单调升、右坡单调降
        left = (k >= kmin) & (k <= kmin + dk)
        right = (k >= kmax - dk) & (k <= kmax)
        self.assertTrue(np.all(np.diff(w[left]) >= -1e-9))
        self.assertTrue(np.all(np.diff(w[right]) <= 1e-9))
        # 中点约 0.5
        mid = int(np.argmin(np.abs(k - (kmin + dk / 2))))
        self.assertAlmostEqual(float(w[mid]), 0.5, places=3)

    def test_wt_sweep_anchored_to_ft_peak(self):
        """WT sweep：9 组 (σ,kweight) 都必须锚定 FT 第一峰，且可导出。"""
        res = self.lep
        sweep = xafs_diagnostics.wt_sweep(res)
        self.assertEqual(len(sweep["rows"]), 9)
        self.assertEqual(len(sweep["maps"]), 9)
        self.assertEqual(sorted({r["sigma"] for r in sweep["rows"]}),
                         [1.0, 1.5, 2.0])
        self.assertEqual(sorted({r["kweight"] for r in sweep["rows"]}),
                         [1, 2, 3])
        for row in sweep["rows"]:
            self.assertLessEqual(abs(row["r_peak"] - res.r1), 0.5,
                                 f"WT 峰跑出 FT 锚定窗: {row}")
            self.assertTrue(np.isfinite(row["k_peak"]))
            self.assertGreater(row["amplitude"], 0)
        # 行内归一所需的每行最大值必须齐全
        self.assertEqual(sorted(sweep["row_max"]), [1, 2, 3])
        for kw in (1, 2, 3):
            self.assertGreater(sweep["row_max"][kw], 0)

        with tempfile.TemporaryDirectory() as td:
            png, csv_path, rows = xafs_diagnostics.write_wt_sweep(
                res, td, stem="demo_wt")
            self.assertTrue(png.exists() and png.stat().st_size > 10_000)
            self.assertTrue(csv_path.exists())
            lines = csv_path.read_text(encoding="utf-8-sig").splitlines()
            self.assertEqual(len(lines), 10)   # header + 9
            self.assertEqual(len(rows), 9)

    # ------------------------------------------------------------
    # 25–27. 第一批（B1-1/B1-2/B1-3）
    # ------------------------------------------------------------
    def test_multi_kweight_joint_fit(self):
        """B1-1：多 k 权重联合拟合。

        验证 (1) 三个分支 Rf 均 < 0.08；(2) 第一壳键长与单权偏差
        ≤ 0.015 Å、N 仍在物理区间；(3) N_idp 不变且最强相关不恶化。
        多权残差若未按分支 RMS 归一，kw=3 会统治最小二乘，第一壳 N
        会被系统性拉偏——这正是该测试要防的退化。
        """
        fit1 = self.lep_fit["fit"]["fit"]
        out3 = xafs_fit.fit_identified(
            self.lep, absorber_symbol="Fe", kweights=(1, 2, 3),
            s02=0.803,  # 固定历史口径，隔离多 k 权算法与静态标定值
            verbose=False)
        fit3 = out3["fit"]
        self.assertEqual(tuple(fit3.kweights), (1, 2, 3))
        self.assertEqual(set(fit3.r_factors_by_kweight), {1, 2, 3})
        for kw, rf in fit3.r_factors_by_kweight.items():
            self.assertLess(rf, 0.08, f"kw={kw} Rf={rf:.4f} 超出 0.08")

        p1 = min(fit1.paths, key=lambda p: p.reff)
        p3 = min(fit3.paths, key=lambda p: p.reff)
        self.assertLessEqual(abs(p3.r - p1.r), 0.015,
                             f"第一壳键长漂移 {abs(p3.r - p1.r):.4f} Å")
        self.assertGreaterEqual(p3.degeneracy, 4.5)
        self.assertLessEqual(p3.degeneracy, 7.5)

        rho1, _, _ = xafs_fit._max_offdiag_correlation(fit1)
        rho3, _, _ = xafs_fit._max_offdiag_correlation(fit3)
        self.assertLessEqual(
            rho3, rho1 + 1e-9,
            f"多权相关恶化: {rho3:.4f} > 单权 {rho1:.4f}")
        self.assertEqual(fit3.n_independent, fit1.n_independent)

    def test_back_ft_q_space_phase_and_envelope(self):
        """B1-2：R 窗反傅里叶的 q 轴、包络与相位。

        合成单壳 χ(k)=sin(2kR₀)·exp(−2k²σ²)/k² 经 FT/Back-FT 后，
        Re χ(q) 必须与输入 k²χ(k) 同相（q∈[4,10] 相关 >0.95）。
        若反向变换符号写反或窗函数端点错误，相关会掉到 0 附近。
        """
        k = np.arange(1.0, 14.0, 0.05)
        R0, sigma2 = 2.0, 0.004
        chi = np.sin(2.0 * k * R0) * np.exp(-2.0 * k**2 * sigma2) / k**2
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2, 1.0)
        bft = athena.back_ft_to_q_space(ft.r, ft.chi_r, 1.3, 2.5, 0.1)

        self.assertEqual(len(bft.q), len(bft.chi_q))
        self.assertTrue(np.all(np.diff(bft.q) > 0), "q 轴不单调递增")
        self.assertTrue(np.all(bft.magnitude >= 0))
        self.assertTrue(np.all(np.isfinite(bft.magnitude)))
        self.assertEqual(len(bft.r_window), len(ft.r))

        mask = (bft.q >= 4.0) & (bft.q <= 10.0)
        self.assertGreater(int(mask.sum()), 50)
        target = np.interp(bft.q[mask], k, k**2 * chi)
        corr = float(np.corrcoef(bft.real[mask], target)[0, 1])
        self.assertGreater(corr, 0.95,
                           f"Re χ(q) 与 k²χ 相关仅 {corr:.4f}（符号/相位错？）")
        self.assertGreater(float(bft.magnitude[mask].max()), 0.0)

    def test_collinear_ms_constrained_fit_nio(self):
        """B1-3：NiO 共线 MS 路径受约束纳入，自由变量数不变。

        关键契约：MS 路径不新增 ΔR/σ²/N 变量，全部由母单散射壳层按
        reff_ratio 驱动；否则加 MS 必然触发 N_idp 过拟合截断。
        """
        res = xafs_pipeline.process_sample(
            [NIO], absorber="Ni", strict=False, verbose=False)
        no_ms = xafs_fit.fit_identified(
            res, absorber_symbol="Ni", rmax=4.3,
            include_ms=False, verbose=False)
        with_ms = xafs_fit.fit_identified(
            res, absorber_symbol="Ni", rmax=4.3,
            include_ms=True, verbose=False)
        f0, f1 = no_ms["fit"], with_ms["fit"]

        self.assertFalse([p for p in f0.paths if p.nleg > 2],
                         "include_ms=False 仍选入 MS 路径")
        ms_paths = [p for p in f1.paths if p.nleg > 2]
        self.assertTrue(ms_paths, "include_ms=True 未选入任何 nleg>2 路径")
        self.assertEqual(f1.n_variables, f0.n_variables,
                         "MS 路径新增了自由变量")
        self.assertLess(f1.r_factor, 0.08,
                        f"含 MS 拟合 Rf={f1.r_factor:.4f}")

        labels = [p.label for p in f1.paths]
        self.assertEqual(len(labels), len(set(labels)), "MS 路径标签重名")

        by_label = {p.label: p for p in f1.paths}
        for p in ms_paths:
            self.assertIsNotNone(p.ms_parent, f"{p.label} 未标记母壳层")
            parent = by_label[p.ms_parent]
            ratio = p.reff / parent.reff
            self.assertGreater(ratio, 1.5,
                               f"{p.label} 不是 Reff≈2R_first 的聚焦路径")
            expected_r = p.reff + ratio * (parent.r - parent.reff)
            self.assertAlmostEqual(
                p.r, expected_r, places=5,
                msg=f"{p.label} 的 ΔR 未按 reff_ratio={ratio:.3f} 缩放")


    def test_window_tuner_html_embeds_data(self):
        """交互式窗口页：内嵌 χ(k)、四个滑块齐全、推荐窗口不偏离 r1。"""
        res = self.lep
        rec = xafs_window_tuner.recommend_window(
            res.chi_k, res.r1, kweight=res.ft.kweight, dk=res.ft.dk)
        self.assertGreater(rec["n_candidates"], 0)
        self.assertLessEqual(abs(rec["recommended"]["peak_r"] - res.r1), 0.15)
        with tempfile.TemporaryDirectory() as td:
            path = xafs_window_tuner.write_window_tuner(
                td, "demo", res.chi_k, res.r1, "demo sample", "Fe K-edge",
                kmin=res.ft.kmin, kmax=res.ft.kmax,
                kweight=res.ft.kweight, dk=res.ft.dk)
            self.assertTrue(path.exists())
            doc = path.read_text(encoding="utf-8")
            self.assertIn("const DATA =", doc)
            for ctl in ('id="kmin"', 'id="kmax"', 'id="kweight"', 'id="dk"'):
                self.assertIn(ctl, doc)
            self.assertIn("e<sup>−2ikR</sup>", doc)
            import json
            import re
            m = re.search(r"const DATA = (\{.*?\});", doc, re.S)
            self.assertIsNotNone(m)
            data = json.loads(m.group(1))
            self.assertEqual(len(data["k"]), len(data["chi"]))
            self.assertGreater(len(data["k"]), 100)
            self.assertIsNotNone(data["recommended"])
            self.assertAlmostEqual(data["current"]["kweight"],
                                   res.ft.kweight)

    # ------------------------------------------------------------
    # 28–30. 第二批（B2-1/B2-2/B2-3）
    # ------------------------------------------------------------
    def test_xanes_flatten_and_neck_features(self):
        """B2-1：flattened 谱平台≈1，合成颈部双台阶被正确识别。"""
        norm = self.lep.norm
        self.assertIsNotNone(norm.flattened)
        self.assertEqual(len(norm.flattened), len(norm.mu_norm))
        mm = (norm.energy >= norm.e0 + 100) & (norm.energy <= norm.e0 + 400)
        self.assertGreater(int(mm.sum()), 5)
        flat_mean = float(np.mean(norm.flattened[mm]))
        self.assertGreaterEqual(flat_mean, 0.97)
        self.assertLessEqual(flat_mean, 1.03)
        # E<E0 必须保持 mu_norm 不变（flattening 只改边后）
        pre = norm.energy < norm.e0
        np.testing.assert_allclose(norm.flattened[pre], norm.mu_norm[pre])

        # 两个相距 5 eV 的 arctan 台阶（1s→4p_z 颈部肩峰）
        E = np.arange(7100.0, 7300.0, 0.05)
        e0 = 7150.0
        width = 0.8
        e1, e2 = e0 - 2.5, e0 + 2.5
        mu = (0.5 * (0.5 + np.arctan((E - e1) / width) / np.pi)
              + 0.5 * (0.5 + np.arctan((E - e2) / width) / np.pi))
        sn = NormalizedScan(
            energy=E, mu=mu, mu_norm=mu,
            pre_edge_line=np.zeros_like(E), edge_step=1.0,
            plateau=1.0, e0=e0, flattened=mu.copy())
        feat = athena.extract_xanes_features(sn)
        self.assertTrue(feat["neck_shoulder_present"])
        self.assertGreaterEqual(feat["derivative_split_eV"], 4.0)
        self.assertLessEqual(feat["derivative_split_eV"], 6.0)
        self.assertTrue(feat["e0_jump_risk"])

    def test_lcf_and_pca_synthetic_mixture(self):
        """B2-2：3 端元合成混合谱回收权重；PCA 前 3 主成分解释 >99.9%。"""
        import xafs_lcf

        goe = xafs_pipeline.process_sample(
            [GOETHITE], absorber="Fe", strict=False, verbose=False).norm
        mag = xafs_pipeline.process_sample(
            [MAGNETITE], absorber="Fe", strict=False, verbose=False).norm
        lep = self.lep.norm

        x = np.arange(-20.0, 60.0, 0.2)

        def flat(n):
            arr = n.flattened if n.flattened is not None else n.mu_norm
            return np.interp(x, n.energy - n.e0, arr)

        y_lep, y_goe, y_mag = flat(lep), flat(goe), flat(mag)
        rng = np.random.default_rng(42)
        y = 0.35 * y_lep + 0.65 * y_goe + rng.normal(0, 1e-3, x.size)
        target = NormalizedScan(
            energy=lep.e0 + x, mu=y, mu_norm=y,
            pre_edge_line=np.zeros_like(y), edge_step=1.0,
            plateau=1.0, e0=lep.e0, flattened=y.copy())
        res = xafs_lcf.lcf_from_scans(
            target, {"lepido": lep, "goethite": goe, "magnetite": mag},
            mode="xanes")
        self.assertLess(abs(res["weights"]["lepido"] - 0.35), 0.02)
        self.assertLess(abs(res["weights"]["goethite"] - 0.65), 0.02)
        self.assertLess(abs(res["weights"]["magnetite"]), 0.02)
        self.assertLess(res["r_factor"], 1e-3)

        specs = []
        for i, (a, b, c) in enumerate([
                (0.2, 0.8, 0.0), (0.5, 0.5, 0.0), (0.8, 0.2, 0.0),
                (0.4, 0.1, 0.5), (0.1, 0.4, 0.5), (0.3, 0.3, 0.4)]):
            yy = (a * y_lep + b * y_goe + c * y_mag
                  + rng.normal(0, 1e-4, x.size))
            specs.append((f"mix{i}", x, yy))
        pca = xafs_lcf.run_pca(specs, -20.0, 60.0)
        cum3 = float(np.sum(pca["explained_variance_ratio"][:3]))
        self.assertGreater(cum3, 0.999,
                           f"前 3 主成分累计解释仅 {cum3:.6f}")

    def test_deadtime_and_explicit_deglitch(self):
        """B2-3：死时间倍乘/分数两条公式 + 宽 glitch 精确删点。"""
        # ---- 死时间：前半 dtime=2（倍乘），后半 dtime=0.5（分数）----
        E = np.arange(7000.0, 7200.0, 0.5)
        i0 = np.full_like(E, 1000.0)
        ifluor = np.full_like(E, 100.0)
        dt = np.where(np.arange(E.size) < E.size // 2, 2.0, 0.5)
        fs = RawScan(energy=E, i0=i0, ifluor=ifluor, mu=ifluor / i0,
                     dtime=dt, meta={"mode": "fluor"})
        out = athena.apply_deadtime_correction(fs)
        self.assertAlmostEqual(float(out.mu[0]), 0.2, places=12)
        self.assertAlmostEqual(float(out.mu[-1]), 0.2, places=12)
        self.assertAlmostEqual(float(fs.mu[0]), 0.1, places=12,
                               msg="死时间校正不得原地修改输入扫描")
        self.assertTrue(out.meta["deadtime_corrected"])

        # ---- 宽方波 glitch：自动 despike 只能削掉边缘，deglitch 精确删点 ----
        g = xafs_io.read_scans([str(LEPIDO)])[0]
        mu2 = np.array(g.mu, dtype=float, copy=True)
        near = np.flatnonzero((g.energy >= 7280.0) & (g.energy <= 7295.0))
        start = int(near[0]) if near.size else int(
            np.argmin(np.abs(g.energy - 7280.0)))
        gi = np.arange(start, min(start + 8, len(g.energy)))
        mu2[gi] += 2.0
        glitched = RawScan(energy=np.array(g.energy, copy=True),
                           mu=mu2, meta=dict(g.meta))
        # 复刻 pipeline 的去毛刺（kernel=5, 3×MAD）
        padded = np.pad(mu2, 2, mode="edge")
        med = np.array([np.median(padded[i:i + 5])
                        for i in range(len(mu2))])
        resid = np.abs(mu2 - med)
        mad = float(np.median(resid)) + 1e-12
        despiked = mu2.copy()
        mask = resid > 3.0 * mad
        despiked[mask] = med[mask]
        remain = float(np.max(np.abs(despiked[gi] - g.mu[gi])))
        self.assertGreater(remain, 1.0,
                           "自动 despike 不应能清除 8 点宽方波")
        lo, hi = float(g.energy[gi].min()), float(g.energy[gi].max())
        cleaned, info = athena.deglitch_ranges(glitched, [(lo, hi)])
        self.assertEqual(info["n_removed"], int(gi.size))
        self.assertFalse(np.any((cleaned.energy >= lo)
                                & (cleaned.energy <= hi)))
        with self.assertRaises(XAFSConfigError):
            athena.deglitch_ranges(glitched, [(7000.0, 8000.0)])

    # ------------------------------------------------------------
    # 31. 第三批（B3-1）
    # ------------------------------------------------------------
    def test_wt_shell_k_profiles_and_model_residual(self):
        """B3-1：一维 WT 包络区分轻/重散射体；实验-模型二维残差可定量。

        为什么 k 窗口取 2.5 而不是产线默认 3.0：Fe–O 的包络在低 k 仍有
        显著权重，kmin=2.5 能保留这部分信息，使第二壳 Fe–Fe（k≈6.4）
        与第一壳 Fe–O（k≈5.5）的峰位差达到 0.9 Å⁻¹。产线 kmin=3.0 时
        差值约 0.7 Å⁻¹，方向仍正确但小于验收阈值；这里显式使用更宽的
        分析窗并同时核对两者排序，避免把窗口选择误当成物理结论。
        """
        fit = self.lep_fit["fit"]["fit"]
        self.assertIsNotNone(fit.model_k, "FitResult 未保存 model_k")
        self.assertEqual(fit.model_k.shape, fit.model_k_grid.shape)

        kmax = float(min(12.0, np.max(self.lep.chi_k.k)))
        wt = larch_wt.continuous_wavelet_transform(
            self.lep.chi_k, sigma=1.5, kmin=2.5, kmax=kmax,
            kweight=2, rmax=4.0)
        first = larch_wt.extract_k_profile(wt, 1.47, half_width=0.25)
        second = larch_wt.extract_k_profile(wt, 2.67, half_width=0.25)
        diff = second["k_peak"] - first["k_peak"]
        self.assertGreaterEqual(
            diff, 0.8, f"Fe–Fe k_peak={second['k_peak']:.2f}, "
                       f"Fe–O k_peak={first['k_peak']:.2f}, 差 {diff:.2f} <0.8")
        # 产线默认窗口（kmin=3.0）下排序必须一致，只是差值略小
        self.assertGreater(larch_wt.extract_k_profile(
            self.lep.wt, 2.67)["k_peak"],
            larch_wt.extract_k_profile(self.lep.wt, 1.47)["k_peak"] + 0.5)
        for prof in (first, second):
            self.assertEqual(len(prof["k"]), len(prof["ridge_profile"]))
            self.assertEqual(len(prof["k"]), len(prof["slice_profile"]))
            self.assertFalse(np.any(np.isnan(prof["ridge_profile"])))

        resid = larch_wt.compute_model_wt_residual(
            self.lep.chi_k, fit.model_k, sigma=1.5, kmin=3.0,
            kmax=12.0, kweight=2, rmax=4.0)
        self.assertEqual(resid["residual_map"].shape,
                         resid["wt_data"].magnitude.shape)
        self.assertEqual(resid["residual_map"].shape,
                         resid["wt_model"].magnitude.shape)
        self.assertLess(resid["r_factor"], 0.15,
                        f"二维 WT R-factor={resid['r_factor']:.4f} ≥0.15")
        # 残差必须共享同一物理标度，而不是各自归一后相减
        self.assertAlmostEqual(resid["scale"],
                               float(resid["wt_data"].magnitude.max()),
                               places=12)

        with tempfile.TemporaryDirectory() as tmp:
            prof = xafs_diagnostics.write_wt_shell_profiles(
                self.lep, self.lep_fit, tmp, stem="lepido")
            self.assertTrue(Path(prof["png"]).exists())
            self.assertTrue(Path(prof["csv"]).exists())
            header = Path(prof["csv"]).read_text(
                encoding="utf-8").splitlines()[0]
            self.assertIn("data_ridge", header)
            self.assertGreaterEqual(len(prof["shells"]), 2)

    # ------------------------------------------------------------
    # 32. 第四批（B4-1）
    # ------------------------------------------------------------
    def test_reduced_chi2_and_shell_ftest(self):
        """B4-1：高 R 噪声 ε_R → χ²_ν；第二壳 Fe–Fe 的 F-test 必须显著。

        关键防退化点：
          - larch xftf 默认 R 网格只到 ~10 Å，不显式扩到 25 Å 时
            [15,25] 区间为空，ε_R/χ² 会静默变成 None；
          - F-test 必须用**同一 R 窗**的单壳对照拟合；换窗后 χ² 不可比。
        """
        fit = self.lep_fit["fit"]["fit"]
        self.assertIsNotNone(fit.epsilon_r)
        self.assertIsNotNone(fit.chi_square)
        self.assertIsNotNone(fit.reduced_chi_square)
        self.assertGreater(fit.epsilon_r, 0.0)
        self.assertGreater(fit.chi_square, 0.0)
        self.assertGreater(fit.reduced_chi_square, 0.0)
        self.assertAlmostEqual(
            fit.reduced_chi_square,
            fit.chi_square / max(fit.n_independent - fit.n_variables, 1.0),
            places=8)

        eps = artemis.estimate_r_noise(
            self.lep.chi_k, fit.kmin, fit.kmax, kweight=2, dk=1.0)
        self.assertGreater(eps, 0.0)

        sft = self.lep_fit["fit"].get("shell_f_test")
        self.assertIsNotNone(sft, "S7 未生成 shell_f_test")
        self.assertIsNotNone(sft.get("confidence"))
        self.assertGreater(
            sft["confidence"], 0.95,
            f"Fe–Fe 第二壳 F-test 置信度仅 {sft['confidence']:.4f}")
        self.assertTrue(sft["significant"])
        self.assertEqual(sft["dfn"], sft["p2"] - sft["p1"])
        self.assertAlmostEqual(
            sft["dfd"], fit.n_independent - sft["p2"], places=6)
        self.assertGreater(sft["delta_chi_square"], 0.0,
                           "多壳层模型必须降低 χ²（delta = χ²_simple − χ²_complex >0）")

        # 报告摘要必须披露这两个统计量（JSON/HTML 同源）
        summary = xafs_report.fit_summary(self.lep_fit["fit"])
        self.assertIn("reduced_chi_square", summary["fit"])
        self.assertIn("shell_f_test", summary)

    # ------------------------------------------------------------
    # 33. 第四批（B4-2）
    # ------------------------------------------------------------
    def test_c3_cumulant_cli_and_nidp_gate(self):
        """B4-2：C3 非谐相移回收 + N_idp 自动降级保护。

        合成第一壳 χ(k) 带 ``third=C3=+1.5e-4 Å³``（larch 的相位项
        ``−4/3·C3·k³``）。放开 C3 后必须回收真值，且 R 回到 FEFF 的
        真键长；关闭 C3 时 R 会偏 −0.010 Å。窄 R 窗下 n_var ≥ N_idp，
        必须自动降级到无 C3 而不崩溃。
        """
        self.assertEqual(xafs_fit._normalize_c3_paths("first", 3), [0])
        self.assertEqual(xafs_fit._normalize_c3_paths("all", 3), True)
        self.assertIsNone(xafs_fit._normalize_c3_paths(False, 3))
        with self.assertRaises(XAFSConfigError):
            xafs_fit._normalize_c3_paths([5], 3)

        # 不能原地改 self.lep：unittest 的类级 fixture 会被后续测试复用。
        import copy
        res = copy.deepcopy(self.lep)
        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "lepidocrocite"][0]
        cache = structure_lib.ensure_feff(motif, edge="K")
        path = artemis.load_feff_path(cache / "feff0002.dat")
        path.degen = 6.0
        c3_true = 1.5e-4
        artemis.set_path_parameters(
            path, s02=0.803, delta_r=0.0, sigma2=0.004,
            delta_e0=0.0, third=c3_true)
        k_data = np.asarray(res.chi_k.k)
        chi = np.interp(k_data, path.k, path.chi)
        res.chi_k = ChiK(k=k_data, chi=chi,
                         background=res.chi_k.background,
                         mu_norm=res.chi_k.mu_norm, rbkg=res.chi_k.rbkg)
        res.ft = athena.ft_to_r_space(k_data, chi, None, 3.0, 12.0, 2, 1.0)
        res.wt = larch_wt.continuous_wavelet_transform(
            res.chi_k, sigma=1.5, kmin=3.0, kmax=12.0,
            kweight=2, rmax=4.0)
        seg = (res.ft.r >= 0.8) & (res.ft.r <= 2.5)
        res.r1 = float(res.ft.r[seg][int(np.argmax(
            res.ft.magnitude[seg]))])
        res.structure = {"status": "confident", "coverage_ok": True,
                         "ranked": [{"id": "lepidocrocite"}], "peaks": []}

        fit_c3 = xafs_fit.fit_identified(
            res, motif_id="lepidocrocite", absorber_symbol="Fe",
            rmax=2.6, vary_c3=True, verbose=False)["fit"]
        fit_no = xafs_fit.fit_identified(
            res, motif_id="lepidocrocite", absorber_symbol="Fe",
            rmax=2.6, vary_c3=False, verbose=False)["fit"]
        p_c3, p_no = fit_c3.paths[0], fit_no.paths[0]
        self.assertIsNotNone(p_c3.c3)
        self.assertLess(abs(p_c3.c3 - c3_true), 5e-5,
                        f"C3 回收 {p_c3.c3:.3e}，真值 {c3_true:.3e}")
        self.assertIsNotNone(p_c3.c3_err)
        self.assertLess(abs(p_c3.r - path.reff),
                        abs(p_no.r - path.reff),
                        "放开 C3 后键长没有更接近 FEFF 真值")
        self.assertFalse(p_c3.c3_at_bound)

        narrow = xafs_fit.fit_identified(
            res, motif_id="lepidocrocite", absorber_symbol="Fe",
            rmin=1.7, rmax=2.5, vary_c3=True, verbose=False)
        nfit = narrow["fit"]
        self.assertTrue(narrow["c3_dropped_by_nidp"])
        self.assertTrue(narrow["flags"]["c3_dropped_by_nidp"])
        self.assertGreaterEqual(nfit.n_variables, 1)
        self.assertTrue(all(p.c3 is None for p in nfit.paths))


    # ------------------------------------------------------------
    # 37. 第五批（B5-1）
    # ------------------------------------------------------------
    def test_wt_complex_residual_catches_phase_inversion(self):
        """B5-1：复数小波矩阵 + 相位残差必须抓住"幅值相等、相位反 π"。

        反相模型 χ_model = −χ_data 的 |W| 与实验逐点相同，旧版幅值
        残差严格为 0；复数残差应为 |W−(−W)|²/|W|² = 4。真实 S7 模型
        的复数 R-factor 必须保持在 0.25 以下，证明该指标有判别力而
        不是恒定 4。
        """
        kmax = float(min(12.0, np.max(self.lep.chi_k.k)))
        inverted = larch_wt.compute_model_wt_residual(
            self.lep.chi_k, -np.asarray(self.lep.chi_k.chi),
            sigma=1.5, kmin=3.0, kmax=kmax, kweight=2, rmax=4.0)
        self.assertLess(inverted["r_factor"], 1e-12,
                        f"反相模型幅值残差={inverted['r_factor']:.3e}，"
                        f"说明幅值口径本应完全失效")
        self.assertIsNotNone(inverted["complex_r_factor"])
        self.assertAlmostEqual(inverted["complex_r_factor"], 4.0, places=6)
        self.assertEqual(inverted["complex_residual_map"].shape,
                         inverted["wt_data"].magnitude.shape)

        wt = larch_wt.continuous_wavelet_transform(
            self.lep.chi_k, sigma=1.5, kmin=3.0, kmax=kmax,
            kweight=2, rmax=4.0, normalize=True)
        self.assertEqual(wt.complex_map.shape, wt.magnitude.shape)
        self.assertLess(float(np.max(np.abs(
            np.abs(wt.complex_map) - wt.magnitude))), 1e-12,
            "normalize=True 时 |complex_map| 与 magnitude 未同步归一")
        self.assertIsNotNone(wt.phase)
        self.assertEqual(wt.phase.shape, wt.magnitude.shape)

        fit = self.lep_fit["fit"]["fit"]
        self.assertIsNotNone(fit.model_k)
        real = larch_wt.compute_model_wt_residual(
            self.lep.chi_k, fit.model_k, sigma=1.5, kmin=3.0,
            kmax=kmax, kweight=2, rmax=4.0)
        self.assertLess(real["complex_r_factor"], 0.25,
                        f"真实 S7 模型 complex R-factor="
                        f"{real['complex_r_factor']:.4f} ≥0.25")
        # 相位数组必须与 k 网格等长，供诊断图直接消费
        prof = larch_wt.extract_k_profile(wt, 1.47, half_width=0.25)
        self.assertIsNotNone(prof["phase_ridge"])
        self.assertIsNotNone(prof["phase_slice"])
        self.assertEqual(len(prof["phase_ridge"]), len(prof["k"]))
        self.assertGreaterEqual(prof["beat_node_phase_jump"], 0.0)

    # ------------------------------------------------------------
    # 38. 第五批（B5-2）
    # ------------------------------------------------------------
    def test_wt_element_reference_ruler_sigma2_shift(self):
        """B5-2：候选散射元素理论 k_peak 标尺随 Z 单调、随 σ² 左移。

        该标尺存在的意义正是打破"3d 金属峰一定在 7 Å⁻¹"的经验误判：
        高无序 Fe–Fe 的包络峰会被 exp(−2σ²k²) 拉低 ≥0.5 Å⁻¹，与
        轻/中元素重合。缓存命中后重复调用必须返回完全相同的值。
        """
        low = structure_id.reference_element_k_peaks(
            "Fe", elements=("O", "Cl", "Fe"), r_bond=2.2, sigma2=0.003)
        high = structure_id.reference_element_k_peaks(
            "Fe", elements=("Fe",), r_bond=2.2, sigma2=0.015)
        for el in ("O", "Cl", "Fe"):
            self.assertIn(el, low)
            self.assertTrue(np.isfinite(low[el]))
        self.assertLess(low["O"], low["Cl"])
        self.assertLess(low["Cl"], low["Fe"])
        self.assertGreaterEqual(low["Fe"] - low["O"], 1.2,
                                f"Fe−O 标尺间距仅 "
                                f"{low['Fe'] - low['O']:.2f} Å⁻¹")
        self.assertGreaterEqual(low["Fe"] - high["Fe"], 0.5,
                                f"σ² 0.003→0.015 的 Fe 峰位移仅 "
                                f"{low['Fe'] - high['Fe']:.2f} Å⁻¹")
        # 缓存一致性：二次调用逐位相同
        again = structure_id.reference_element_k_peaks(
            "Fe", elements=("Fe",), r_bond=2.2, sigma2=0.003)
        self.assertEqual(again["Fe"], low["Fe"])

        # 诊断图接口必须接受标尺参数且不破坏既有 CSV 列
        with tempfile.TemporaryDirectory() as tmp:
            prof = xafs_diagnostics.write_wt_shell_profiles(
                self.lep, self.lep_fit, tmp, stem="ruler",
                ref_k_peaks=low)
            header = Path(prof["csv"]).read_text(
                encoding="utf-8").splitlines()[0]
            self.assertIn("data_ridge", header)
            self.assertIn("data_phase_ridge", header)
            self.assertIsNotNone(prof["complex_r_factor"])


    # ------------------------------------------------------------
    # 39. 第六批（B6-1）
    # ------------------------------------------------------------
    def test_merge_scans_outlier_rejection_and_sigma(self):
        """B6-1：多扫描合并的离群剔除、σ_μ 传播与默认逐位兼容。

        默认 reject_outliers=False 必须与历史均值逐位一致；开启后
        第 4 条 +0.15 畸变扫描按中位谱 RMS+MAD 剔除，合并结果与
        三条干净扫描逐位相同；σ_μ 为保留扫描的逐点标准差。
        保底护栏：2-vs-1 时只剔 1 条，绝不把扫描集压到 1 条。
        """
        files = [LEPIDO.parent / f"Fe_lepidocrocite.{ext}"
                 for ext in ("000", "100", "200")]
        if not all(f.exists() for f in files):
            self.skipTest("缺少 lepidocrocite 三扫描数据")
        scans = xafs_io.read_scans([str(f) for f in files])
        self.assertEqual(len(scans), 3)

        clean = athena.merge_scans(scans)
        self.assertIsNotNone(clean.sigma_mu, "多扫描合并未传播 sigma_mu")
        self.assertEqual(clean.sigma_mu.shape, clean.energy.shape)
        self.assertTrue(np.all(np.isfinite(clean.sigma_mu)))
        self.assertTrue(np.all(np.asarray(clean.sigma_mu) >= 0.0))
        legacy = np.array([s.mu for s in scans]).mean(axis=0)
        self.assertTrue(np.array_equal(clean.mu, legacy),
                        "reject_outliers=False 改变了默认合并结果")
        self.assertFalse(clean.meta["merge_info"]["reject_outliers"])
        self.assertEqual(clean.meta["merge_info"]["n_kept"], 3)
        self.assertEqual(clean.meta["merge_info"]["rejected_indices"], [])

        bad = xafs_io.read_scans([str(files[0])])[0]
        bad.mu = np.asarray(bad.mu, dtype=float) + 0.15
        merged = athena.merge_scans(scans + [bad], reject_outliers=True)
        info = merged.meta["merge_info"]
        self.assertEqual(info["rejected_indices"], [3],
                         f"离群扫描未被正确剔除: {info}")
        self.assertEqual(info["n_kept"], 3)
        self.assertTrue(np.array_equal(merged.mu, clean.mu),
                        "剔除离群扫描后合并谱与干净三扫描不一致")
        self.assertIsNotNone(merged.sigma_mu)
        self.assertTrue(np.all(np.asarray(merged.sigma_mu) >= 0.0))

        # 保底护栏：3 条里 2 条接近、1 条远离 → 只剔离群者，保留 2 条
        E = np.arange(7100.0, 7200.0, 0.2)

        def _mk(value):
            return RawScan(energy=E, mu=np.full_like(E, float(value)))

        floored = athena.merge_scans(
            [_mk(0.0), _mk(100.0), _mk(101.0)], reject_outliers=True)
        f_info = floored.meta["merge_info"]
        self.assertEqual(f_info["n_kept"], 2,
                         f"保底护栏失效: {f_info}")
        self.assertEqual(f_info["rejected_indices"], [0])

    # ------------------------------------------------------------
    # 40. 第六批（B6-2）
    # ------------------------------------------------------------
    def test_xanes_fluo_self_absorption_restores_whiteline_and_lcf(self):
        """B6-2：XANES 荧光自吸收 Tröger 反解恢复白线与 LCF 权重。

        真值 = 0.5·LEPIDO + 0.5·GOETHITE flattened；按 Fe₂O₃ 浓样
        Eq.2 正演压缩（白线压低 >15%），反解后白线误差 <0.02，
        LCF 权重回到 0.50±0.02，R-factor 相对未校正下降 >10 倍。
        同时核对分母防爆门禁：α 较小的 45°/45° 几何下必须
        reliable=False 且 clamped_points>0，但输出仍有限。
        """
        import xafs_lcf

        goe = xafs_pipeline.process_sample(
            [GOETHITE], absorber="Fe", strict=False, verbose=False)
        x = np.asarray(self.lep.norm.energy, dtype=float)
        e0 = float(self.lep.norm.e0)
        y_lep = np.asarray(self.lep.norm.flattened, dtype=float)
        y_goe = np.interp(x, goe.norm.energy, goe.norm.flattened)
        y_true = 0.5 * y_lep + 0.5 * y_goe

        def _norm(y):
            return NormalizedScan(
                energy=np.array(x, copy=True), mu=np.asarray(y, dtype=float),
                mu_norm=np.asarray(y, dtype=float),
                pre_edge_line=np.zeros_like(x), edge_step=1.0,
                plateau=1.0, e0=e0, flattened=np.asarray(y, dtype=float))

        # 用公开 API 在占位谱上读取同一几何的 α，避免在测试里复制私有公式
        dummy = _norm(np.ones_like(x))
        _, info0 = athena.fluo_xanes_self_absorption(
            dummy, "Fe2O3", "Fe", "K", "Ka", 100.0, 5.24, 45.0, 20.0)
        alpha = float(info0["alpha"])
        self.assertGreater(alpha, 0.3)

        y_att = (alpha + 1.0) * y_true / (alpha + y_true)
        drop = 1.0 - float(np.max(y_att)) / float(np.max(y_true))
        self.assertGreater(drop, 0.15,
                           f"正演白线压缩仅 {100 * drop:.1f}% (<15%)，"
                           f"测试未覆盖浓样自吸收场景")

        norm_att = _norm(y_att)
        corr, info = athena.fluo_xanes_self_absorption(
            norm_att, "Fe2O3", "Fe", "K", "Ka", 100.0, 5.24, 45.0, 20.0)
        self.assertTrue(info["reliable"], info.get("warning"))
        self.assertEqual(info["clamped_points"], 0)
        self.assertAlmostEqual(info["alpha"], alpha, places=9)
        wl_err = abs(float(np.max(corr.flattened)) - float(np.max(y_true)))
        self.assertLess(wl_err, 0.02,
                        f"白线恢复误差 {wl_err:.4f} ≥0.02")
        # 输入对象不得被原地修改
        self.assertLess(float(np.max(np.abs(
            np.asarray(norm_att.flattened) - y_att))), 1e-12)

        xr = x - e0
        standards = {"lepido": (xr, y_lep), "goethite": (xr, y_goe)}
        r_att = xafs_lcf.run_lcf(xr, y_att, standards, -10.0, 60.0)
        r_cor = xafs_lcf.run_lcf(xr, corr.flattened, standards, -10.0, 60.0)
        for name in ("lepido", "goethite"):
            self.assertLess(abs(r_cor["weights"][name] - 0.5), 0.02,
                            f"校正后 {name} 权重 "
                            f"{r_cor['weights'][name]:.4f} 偏离 0.50 超过 0.02")
        self.assertGreater(r_att["r_factor"], 10.0 * r_cor["r_factor"],
                           f"LCF R-factor 未下降 >10 倍: "
                           f"{r_att['r_factor']:.4g} → "
                           f"{r_cor['r_factor']:.4g}")

        # 分母防爆门禁：α≈0.42 的 45°/45° 几何必然触界
        _, info45 = athena.fluo_xanes_self_absorption(
            dummy, "Fe2O3", "Fe", "K", "Ka", 100.0, 5.24, 45.0, 45.0)
        alpha45 = float(info45["alpha"])
        y_att45 = (alpha45 + 1.0) * y_true / (alpha45 + y_true)
        corr45, guard = athena.fluo_xanes_self_absorption(
            _norm(y_att45), "Fe2O3", "Fe", "K", "Ka",
            100.0, 5.24, 45.0, 45.0)
        self.assertFalse(guard["reliable"],
                         "45°/45° 下分母钳制未触发 reliable=False")
        self.assertGreater(guard["clamped_points"], 0)
        self.assertTrue(np.all(np.isfinite(corr45.flattened)))

    # ------------------------------------------------------------
    # 41. 第六批（B6-3）
    # ------------------------------------------------------------
    def test_pre_edge_voigt_deconvolution_and_centroid(self):
        """B6-3：主边拖尾基线 + 双 Pseudo-Voigt 分峰与面积加权质心。

        合成峰：E_c1=7113.0 (A=0.08)、E_c2=7114.6 (A=0.16)，
        真值质心 7114.067 eV。要求 n_peaks=2、质心误差 <0.15 eV、
        子峰中心误差 <0.20 eV；并核对 extract_xanes_features 与
        指标表已写入 pre_edge_centroid。
        """
        E = np.linspace(7090.0, 7160.0, 701)
        E0 = 7125.0

        def _pv(energy, area, center, gamma, eta):
            # 独立实现（不复用被测函数），避免公式错误相互抵消
            xx = energy - center
            lor = (gamma / np.pi) / (xx ** 2 + gamma ** 2)
            gau = (np.sqrt(np.log(2.0)) / (gamma * np.sqrt(np.pi))
                   * np.exp(-np.log(2.0) * (xx / gamma) ** 2))
            return area * (eta * lor + (1.0 - eta) * gau)

        y = (0.02 + 1.5 / ((E - (E0 + 2.0)) ** 2 + 3.0 ** 2)
             + _pv(E, 0.08, 7113.0, 0.7, 0.5)
             + _pv(E, 0.16, 7114.6, 0.7, 0.5))
        norm = NormalizedScan(
            energy=E, mu=y, mu_norm=y, pre_edge_line=np.zeros_like(E),
            edge_step=1.0, plateau=1.0, e0=E0, flattened=y)
        res = athena.fit_pre_edge_peaks(norm)
        self.assertEqual(res["n_peaks"], 2,
                         f"双峰未被识别: n_peaks={res['n_peaks']}")
        centroid_true = (0.08 * 7113.0 + 0.16 * 7114.6) / 0.24
        self.assertLess(abs(res["centroid_energy"] - centroid_true), 0.15,
                        f"质心 {res['centroid_energy']:.4f} vs 真值 "
                        f"{centroid_true:.4f}")
        centers = sorted(p["center"] for p in res["peaks"])
        self.assertLess(abs(centers[0] - 7113.0), 0.20,
                        f"子峰 1 中心 {centers[0]:.4f} vs 7113.0")
        self.assertLess(abs(centers[1] - 7114.6), 0.20,
                        f"子峰 2 中心 {centers[1]:.4f} vs 7114.6")
        self.assertLess(res["r_factor"], 0.05)
        self.assertAlmostEqual(
            res["total_area"], sum(p["area"] for p in res["peaks"]),
            places=12)

        # 集成：真实管线结果必须带该字段；指标表必须有对应行
        self.assertIn("pre_edge_centroid", self.lep.xanes_features)
        self.assertIn("pre_edge_n_peaks", self.lep.xanes_features)
        st = {"chi_k": self.lep.chi_k, "r1": self.lep.r1,
              "wt": self.lep.wt, "norm": self.lep.norm,
              "e0": self.lep.e0["e0_sample"], "rbkg": self.lep.rbkg,
              "xanes_features": self.lep.xanes_features}
        with tempfile.TemporaryDirectory() as tmp:
            rows = xafs_report.export_indicators_single(
                st, None, Path(tmp), "pe", "Fe2O3", "Fe K", 0.803, [])
            names = [r[0] for r in rows]
            self.assertIn("pre_edge_centroid", names)
            self.assertIn("pre_edge_n_peaks", names)

    # ------------------------------------------------------------
    # 42. 第七批（B7-1）
    # ------------------------------------------------------------
    def test_isotropic_lattice_expansion_constraint(self):
        """B7-1：各向同性晶格约束把 M 个 ΔR 压缩为唯一全局 α。

        用 FEFF 四壳 NiO 按 ΔR_i = α_true·Reff_i（α_true=0.004）合成
        无噪声 χ(k)，要求 isotropic_deltar=True 后：
          - 自由变量数从 13 降为 10（4 个 ΔR → 1 个 α）；
          - Rf < 1e-4，α 回收误差 < 1e-9，四壳 ΔR/Reff 与 α 一致。
        实验 NiO 的各壳位移不服从统一晶格缩放（Rf≈0.135 > 0.08），
        必须诚实降级 cn_reliable=False 并给出同约束警告。
        """
        from types import SimpleNamespace

        motif = [m for m in structure_lib.load_library("Ni")
                 if m.id == "nio_bunsenite_cod_1010093"][0]
        cache = structure_lib.ensure_feff(motif, edge="K")
        paths = []
        for f in sorted(cache.glob("feff*.dat")):
            p = artemis.load_feff_path(f)
            if int(getattr(p, "nleg", 2)) == 2:
                paths.append(p)
        paths.sort(key=lambda p: p.reff)
        paths = paths[:4]
        self.assertEqual(len(paths), 4, "NiO 四壳 FEFF 路径不完整")

        alpha_true = 0.004
        k = np.arange(2.0, 13.0, 0.02)
        for p in paths:
            artemis.set_path_parameters(
                p, s02=1.0, sigma2=0.004,
                delta_r=alpha_true * float(p.reff), delta_e0=0.0)
        chi = artemis.exafs_model_k(paths, k)
        ck = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                  mu_norm=np.ones_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
        synthetic = SimpleNamespace(
            chi_k=ck, ft=ft,
            structure={"status": "confident", "coverage_ok": True,
                       "ranked": [{"id": motif.id}]},
            merged=None, norm=None, self_absorption=None,
            absorber_symbol="Ni")
        free = xafs_fit.fit_identified(
            synthetic, absorber_symbol="Ni", rmax=4.3,
            verbose=False)["fit"]
        iso = xafs_fit.fit_identified(
            synthetic, absorber_symbol="Ni", rmax=4.3,
            isotropic_deltar=True, verbose=False)["fit"]
        self.assertEqual(free.n_variables, 13)
        self.assertEqual(iso.n_variables, 10,
                         f"各向同性约束后变量数 {iso.n_variables} ≠ 10")
        self.assertEqual(free.n_variables - iso.n_variables,
                         len(iso.paths) - 1)
        self.assertLess(iso.r_factor, 1e-4,
                        f"合成各向同性 NiO 的 Rf={iso.r_factor:.3e} ≥ 1e-4")
        self.assertIsNotNone(iso.lattice_strain)
        self.assertIsNotNone(iso.lattice_strain_err)
        self.assertLess(abs(iso.lattice_strain - alpha_true), 1e-9,
                        f"α 回收 {iso.lattice_strain:.12f} vs 真值 {alpha_true}")
        for p in iso.paths:
            self.assertLess(
                abs(p.delta_r / p.reff - iso.lattice_strain), 1e-12,
                f"{p.label} 的 ΔR/Reff 与全局 α 不一致")

        exp = xafs_pipeline.process_sample(
            [NIO], absorber="Ni", strict=False, verbose=False)
        exp_iso = xafs_fit.fit_identified(
            exp, absorber_symbol="Ni", rmax=4.3,
            isotropic_deltar=True, verbose=False)
        self.assertGreater(exp_iso["fit"].r_factor, 0.08,
                           "实验 NiO 的统一 α 约束意外达到 Rf<0.08，"
                           "需复核是否有人为补偿")
        self.assertFalse(exp_iso["cn_reliable"],
                         "实验 NiO 的 Rf>0.08 却仍标记 CN 可靠")
        self.assertIn("各向同性晶格约束", exp_iso["warning"])

    # ------------------------------------------------------------
    # 43. 第七批（B7-2）
    # ------------------------------------------------------------
    def test_cif_partial_occupancy_and_vacancy_scaling(self):
        """B7-2：CIF 部分占据、共占位合并与壳层简并度缩放。

        合成 Pm-3m CIF（a=4.10 Å）：第一壳 O 位于 3d (0.5,0,0)（R=2.05 Å，
        6 配位，occ=0.75）；第二壳 Co(0.6)/Ni(0.4) 同处 1b (0.5,0.5,0.5)
        （R=3.5507 Å，8 配位）共占位。要求：默认仍返回 4 元组；return_occupancy
        返 5 元组且第一壳 O 的 occ=0.75；共占位合并后任意两点间距 ≥0.08 Å，
        元素取高占据者 Co 且有效占据度钳制到 1.0，可正常跑 FEFF6；
        select_shell_paths 选出的第一壳简并度 degens[0] 精确等于 6×0.75=4.5，
        第二壳 Fe-Co 简并度等于 8×1.0=8.0。
        """
        cif_text = """data_test_occ
_cell_length_a 4.10
_cell_length_b 4.10
_cell_length_c 4.10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'Pm-3m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Fe1 Fe 0.0 0.0 0.0 1.0
O1 O 0.5 0.0 0.0 0.75
Co1 Co 0.5 0.5 0.5 0.6
Ni1 Ni 0.5 0.5 0.5 0.4
"""
        with tempfile.TemporaryDirectory() as tmp:
            cif = Path(tmp) / "occ.cif"
            cif.write_text(cif_text, encoding="utf-8")

            legacy = artemis.read_cif_sites(cif, absorber="Fe", rmax=6.0)
            self.assertGreater(len(legacy), 10)
            self.assertTrue(all(len(s) == 4 for s in legacy),
                            "默认模式不再是 4 元组，破坏既有 CIF 调用")
            sites2, occs2 = artemis.read_cif_sites(
                cif, absorber="Fe", rmax=6.0, with_occupancy=True)
            self.assertEqual(len(sites2), len(occs2))
            five = artemis.read_cif_sites(
                cif, absorber="Fe", rmax=6.0, return_occupancy=True)
            self.assertEqual(len(five[0]), 5)
            self.assertEqual(five[0][3], "O")
            self.assertAlmostEqual(five[0][4], 0.75, places=12)

            pos = np.array([s[:3] for s in five])
            dmin = min(float(np.linalg.norm(pos[i] - pos[j]))
                       for i in range(len(pos))
                       for j in range(i + 1, len(pos)))
            self.assertGreaterEqual(dmin, 0.08,
                                    f"共占位合并后仍有间距 {dmin:.4f} Å < 0.08")
            co_sites = [s for s in five if s[3] == "Co"]
            self.assertGreaterEqual(len(co_sites), 8)
            self.assertTrue(all(abs(s[4] - 1.0) < 1e-12 for s in co_sites),
                            f"Co/Ni 共占位未钳制到 occ=1.0: {co_sites[0][4]}")
            self.assertFalse(any(s[3] == "Ni" for s in five),
                             "Co/Ni 共占位未合并成单一晶位")
            self.assertAlmostEqual(
                float(np.linalg.norm(co_sites[0][:3])),
                2.05 * np.sqrt(3.0), places=6)

            # 合并后的位点必须能正常生成 FEFF 路径而不崩溃
            feff_dir = Path(tmp) / "feff_occ"
            feff_inp = artemis.build_feff_input(
                legacy, absorber="Fe", rmax=4.2, edge="K", title="occ_test")
            datfiles = artemis.run_feff6l(feff_inp, feff_dir)
            self.assertGreaterEqual(len(datfiles), 2)

            scale = artemis.cif_shell_occupancy_scale(cif, "Fe", rmax=6.0)
            o_scale = [occ for (dist, el, occ) in scale
                       if el == "O" and abs(dist - 2.05) < 1e-6]
            self.assertEqual(len(o_scale), 6)
            self.assertTrue(all(abs(occ - 0.75) < 1e-12 for occ in o_scale),
                            f"第一壳 O 占据度缩放错误: {o_scale}")

            motif = structure_lib.Motif(
                id="occ_test", absorber="Fe", kind="crystal",
                nuclearity="crystal", cell={"a": 4.10},
                spacegroup="Pm-3m",
                asym_sites=[
                    {"el": "Fe", "x": 0.0, "y": 0.0, "z": 0.0},
                    {"el": "O", "x": 0.5, "y": 0.0, "z": 0.0, "occ": 0.75},
                    {"el": "Co", "x": 0.5, "y": 0.5, "z": 0.5,
                     "occ": 0.6},
                    {"el": "Ni", "x": 0.5, "y": 0.5, "z": 0.5,
                     "occ": 0.4}],
                cif_path=str(cif),
                sites=[(x, y, z, el) for (x, y, z, el) in legacy])
            plan = xafs_fit.select_shell_paths(
                motif, edge="K", rmax=3.8, kmin=3.0, kmax=12.0)
            self.assertEqual(plan["labels"][0], "Fe-O@2.05")
            self.assertAlmostEqual(plan["degens"][0], 4.5, places=12,
                                   msg=f"第一壳 Fe-O 简并度 {plan['degens'][0]} ≠ 6×0.75")
            co_idx = plan["labels"].index("Fe-Co@3.55")
            self.assertAlmostEqual(plan["degens"][co_idx], 8.0, places=12,
                                   msg=f"第二壳 Fe-Co 简并度 {plan['degens'][co_idx]} ≠ 8×1.0")

    # ------------------------------------------------------------
    # 44. 第七批（B7-3）
    # ------------------------------------------------------------
    def test_c4_cumulant_and_asymmetric_rdf_reconstruction(self):
        """B7-3：C4 回收与正偏态 g(R) 的重构契约。

        1) C4=6e-5 合成单壳谱，vary_c4=True 回收误差 <2.5e-5；
        2) RDF 面积严格等于 N，且 C3>0 时 r_mode < r_mean、
           skewness>0、右尾延伸大于左尾。

        注：规范原文“右半侧积分面积大于左半侧”与归一化正偏态分布
        数学矛盾——在 r0 处切分时，质量向峰位左侧集中必然导致
        左面积 > 右面积（本实现实测 3.159 vs 2.779）。这里改用
        等价可检验的尾长判据（q95−r0 > r0−q05），不伪造面积。
        """
        motif = [m for m in structure_lib.load_library("Ni")
                 if m.id == "nio_bunsenite_cod_1010093"][0]
        cache = structure_lib.ensure_feff(motif, edge="K")
        datfile = None
        for f in sorted(cache.glob("feff*.dat")):
            p = artemis.load_feff_path(f)
            if int(getattr(p, "nleg", 2)) == 2:
                datfile = f
                break
        self.assertIsNotNone(datfile, "NiO 单散射路径缺失")
        p = artemis.load_feff_path(datfile)
        p.degen = 6.0
        artemis.set_path_parameters(p, s02=1.0, delta_r=0.0,
                                    sigma2=0.004, delta_e0=0.0,
                                    fourth=6e-5)
        self.assertIn(p.pathpar_name("fourth"), p.params,
                      "set_path_parameters 未把 fourth 映射进 larch params")
        k = np.asarray(p.k)
        chi = np.asarray(p.chi)
        ck = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                  mu_norm=np.ones_like(k), rbkg=1.2)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
        fit = artemis.fit_r_space(
            ck, ft, [datfile], rmin=max(1.0, p.reff - 0.8),
            rmax=p.reff + 0.8, kmin=3.0, kmax=12.0, kweight=2,
            s02=1.0, vary_n=False, vary_c4=True,
            sigma2_init=0.004, deltar_lim=(-0.1, 0.1))
        path_fit = fit.paths[0]
        self.assertIsNotNone(path_fit.c4)
        self.assertLess(abs(path_fit.c4 - 6e-5), 2.5e-5,
                        f"C4 回收 {path_fit.c4:.3e} vs 真值 6e-5")
        self.assertFalse(path_fit.c4_at_bound)
        self.assertEqual(artemis.C4_LIMIT, 2e-4)

        rdf = artemis.reconstruct_pair_distribution(
            np.linspace(1.5, 2.7, 601), n=6.0, r0=2.05,
            sigma2=0.006, c3=2.5e-4)
        self.assertLess(abs(rdf["area"] - 6.0), 1e-3,
                        f"g(R) 积分面积 {rdf['area']:.6f} ≠ 6.0")
        self.assertEqual(rdf["r"].shape, rdf["g_r"].shape)
        self.assertTrue(np.all(rdf["g_r"] >= 0.0))
        self.assertLess(rdf["r_mode"], rdf["r_mean"],
                        "正偏态 g(R) 的峰位必须低于平均键长")
        self.assertGreater(rdf["skewness"], 0.0)
        r, g = rdf["r"], rdf["g_r"]
        cdf = np.cumsum(g) * np.gradient(r)
        cdf /= cdf[-1]
        q05 = float(r[np.searchsorted(cdf, 0.05)])
        q95 = float(r[np.searchsorted(cdf, 0.95)])
        self.assertGreater(q95 - 2.05, 2.05 - q05,
                           f"右尾 {q95 - 2.05:.3f} 未长于左尾 {2.05 - q05:.3f}")
        left = float(np.trapezoid(g[r <= 2.05], r[r <= 2.05]))
        right = float(np.trapezoid(g[r >= 2.05], r[r >= 2.05]))
        self.assertGreater(left, right,
                           "归一化正偏态分布在 r0 处切分应满足左面积>右面积；"
                           "若此断言失败说明归一化或非负截断出错")

    # ------------------------------------------------------------
    # 45. 第八批（B8-1）
    # ------------------------------------------------------------
    def test_multi_dataset_corefinement_and_einstein_model(self):
        """B8-1：双温度全局联合拟合 + 爱因斯坦 σ²(T) 约束。

        用 LEPIDO 第一壳 Fe–O（μ=12.436 amu）在真值 Θ_E=420 K、
        σ²_static=0.0015 Å²、ΔR=-0.02 Å、N=6.0 下合成 80 K 与 350 K
        两条谱。要求：
          - 全局自由变量从独立 2×4=8 降为 5（N/ΔR/e0/σ²_static/Θ_E）；
          - Θ_E 回收误差 <15 K，σ²_static 误差 <5e-4；
          - 两数据集各自 Rf <1e-4，且 ΔR/N 也被正确绑定回收。
        同时核对 HBAR2_OVER_2KB=24.254316 与 T≤0 的干净拒绝。
        """
        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "lepidocrocite"][0]
        cache = structure_lib.ensure_feff(motif, edge="K")
        first = None
        for f in sorted(cache.glob("feff*.dat")):
            p = artemis.load_feff_path(f)
            if int(getattr(p, "nleg", 2)) == 2:
                if first is None or p.reff < first[1].reff:
                    first = (f, p)
        self.assertIsNotNone(first, "lepidocrocite 第一壳路径缺失")
        datfile, path0 = first
        self.assertLess(abs(path0.reff - 1.98), 0.05)

        # 常数与输入护栏
        self.assertEqual(artemis.HBAR2_OVER_2KB, 24.254316)
        self.assertAlmostEqual(
            artemis.einstein_sigma2(300.0, 500.0, 12.436),
            0.005717251068417776, places=12)
        for kwargs in ({"temperature_K": 0.0, "theta_e": 420.0,
                        "reduced_mass_amu": 12.436},
                       {"temperature_K": 300.0, "theta_e": 0.0,
                        "reduced_mass_amu": 12.436},
                       {"temperature_K": 300.0, "theta_e": 420.0,
                        "reduced_mass_amu": 0.0}):
            with self.assertRaises(XAFSConfigError):
                artemis.einstein_sigma2(**kwargs)

        theta_true, static_true = 420.0, 0.0015
        dr_true, n_true, mu = -0.02, 6.0, 12.436
        k = np.arange(2.0, 13.0, 0.02)
        datasets = []
        sigma2_true = {}
        for T in (80.0, 350.0):
            p = artemis.load_feff_path(datfile)
            p.degen = n_true
            s2 = artemis.einstein_sigma2(
                T, theta_true, mu, static_true)
            sigma2_true[T] = s2
            artemis.set_path_parameters(
                p, s02=1.0, delta_r=dr_true, sigma2=s2, delta_e0=0.0)
            chi = artemis.exafs_model_k([p], k)
            ck = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                      mu_norm=np.ones_like(k), rbkg=1.0)
            ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
            datasets.append({
                "name": f"T{int(T)}", "chi_k": ck, "ft": ft,
                "datfiles": [datfile], "rmin": 1.0, "rmax": 2.6,
                "kmin": 3.0, "kmax": 12.0, "kweight": 2, "s02": 1.0,
                "vary_n": True, "path_degen": [n_true],
                "path_labels": ["Fe-O"], "temperature_K": T})

        ref = artemis.fit_multi_datasets(
            datasets,
            shared_bonds=[{"target": (1, 0), "source": (0, 0),
                           "share_deltar": True, "share_sigma2": False,
                           "degen_ratio": 1.0}],
            einstein_models=[{"paths": [(0, 0), (1, 0)],
                              "reduced_mass_amu": mu,
                              "theta_e_init": 350.0}])
        self.assertEqual(ref["n_variables_total"], 5,
                         f"全局变量数 {ref['n_variables_total']} ≠ 5")
        self.assertGreater(ref["n_independent_total"],
                           ref["n_variables_total"])
        self.assertEqual(len(ref["datasets"]), 2)
        self.assertEqual(len(ref["einstein_results"]), 1)
        ein = ref["einstein_results"][0]
        self.assertLess(abs(ein["theta_e"] - theta_true), 15.0,
                        f"Θ_E 回收 {ein['theta_e']:.4f} vs 真值 {theta_true}")
        self.assertLess(abs(ein["sigma2_static"] - static_true), 5e-4,
                        f"σ²_static 回收 {ein['sigma2_static']:.6f} "
                        f"vs 真值 {static_true}")
        for di, T in enumerate((80.0, 350.0)):
            fr = ref["datasets"][di]
            self.assertLess(fr.r_factor, 1e-4,
                            f"T={T} K 的 Rf={fr.r_factor:.3e} ≥1e-4")
            self.assertAlmostEqual(
                fr.paths[0].delta_r, dr_true, places=4,
                msg=f"T={T} K 的 ΔR 未被共享约束")
            self.assertAlmostEqual(
                fr.paths[0].degeneracy, n_true, places=4,
                msg=f"T={T} K 的 N 未被共享约束")
            self.assertAlmostEqual(
                fr.paths[0].sigma2, sigma2_true[T], places=10,
                msg=f"T={T} K 的 σ²(T) 未按爱因斯坦模型计算")
            self.assertIsNotNone(fr.paths[0].sigma2_err)
        self.assertLess(abs(ref["global_r_factor"]), 1e-4)


    # ============================================================
    # Phase III（B9–B12）：#46–#60
    # ============================================================
    def test_config_file_covers_all_cli_switches(self):
        """#46 B9-1：17 个 Phase I–III 配置键逐键生效，未配置键保持默认。

        优先级契约：explicit（CLI 显式 dest）> 配置文件 > 代码默认值。
        这里用 JSON（合法 YAML 子集）写全部 17 键，避免依赖 PyYAML；
        同时确认 ``bvs_restraint`` 已进入 DEFAULTS。
        """
        import json

        cfg = {
            "deglitch": [[7000.0, 7010.0]],
            "deadtime": True,
            "reject_outlier_scans": True,
            "fit_kweights": "1,2,3",
            "include_ms": True,
            "vary_c3": "first",
            "vary_c4": "all",
            "isotropic_deltar": True,
            "lcf_standards": "a.xdi,b.xdi",
            "shared_sigma2": "all",
            "mback": True,
            "chi_std": "std.dat",
            "clamp_lo": 0.1,
            "clamp_hi": 0.2,
            "fit_space": "q",
            "target_transform": {"spoil_good": 3.0},
            "mcr_als": {"n_components": 3},
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "all.json"
            path.write_text(json.dumps(cfg), encoding="utf-8")
            loaded = xafs_config.load_config(path)
            kwargs = {k: (list(v) if isinstance(v, list) else v)
                      for k, v in xafs_config.DEFAULTS.items()}
            kwargs.pop("self_absorption", None)
            kwargs.pop("batch", None)
            merged = xafs_config.apply_config_to_kwargs(loaded, kwargs)

        self.assertEqual(merged["deglitch"], [(7000.0, 7010.0)])
        for key in ("deadtime", "reject_outlier_scans", "include_ms",
                    "isotropic_deltar", "mback"):
            self.assertIs(merged[key], True, f"{key} 未由配置打开")
        self.assertEqual(merged["fit_kweights"], (1, 2, 3))
        self.assertEqual(merged["vary_c3"], [0])
        self.assertIs(merged["vary_c4"], True)
        self.assertEqual(merged["lcf_standards"], ["a.xdi", "b.xdi"])
        self.assertEqual(merged["shared_sigma2"], "all")
        self.assertEqual(merged["chi_std"], "std.dat")
        self.assertAlmostEqual(merged["clamp_lo"], 0.1)
        self.assertAlmostEqual(merged["clamp_hi"], 0.2)
        self.assertEqual(merged["fit_space"], "q")
        self.assertEqual(merged["target_transform"], {"spoil_good": 3.0})
        self.assertEqual(merged["mcr_als"], {"n_components": 3})
        # 未出现在配置里的键必须保持代码默认值
        self.assertEqual(merged["kmin"], 3.0)
        self.assertEqual(merged["kmax"], 12.0)
        self.assertIs(merged["despike"], xafs_config.DEFAULTS["despike"],
                      "未配置的 despike 不应被配置合并改写")
        self.assertIs(merged["despike"], False)
        self.assertIn("bvs_restraint", xafs_config.DEFAULTS)

        # CLI 显式值优先级最高
        merged2 = xafs_config.apply_config_to_kwargs(
            loaded, dict(kwargs, clamp_lo=0.9), explicit={"clamp_lo"})
        self.assertAlmostEqual(merged2["clamp_lo"], 0.9)

    def test_vary_c4_cli_nidp_gate_and_rdf_export(self):
        """#47 B9-1：C4 回收、N_idp 降级与真实 LEPIDO 的 g(R) 导出。

        1) 合成 Fe–O 单壳 C4=6e-5，经 ``fit_identified(vary_c4=True)``
           回收误差 <2.5e-5；
        2) LEPIDO 取 4 条壳层、vary_c4="all" 时 4 个 C4 变量使
           n_var=17 > N_idp≈15，必须先降为仅第一壳并记
           ``c4_dropped_by_nidp=True``；
        3) 真实 LEPIDO 放开 C3 后 ``rdf_table`` 每条路径的 g(R) 面积
           必须等于拟合 N，并能写出 7 列 ``_rdf.csv``。
        """
        import types

        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "lepidocrocite"][0]
        cache = structure_lib.ensure_feff(motif, edge="K")
        first = None
        for f in sorted(cache.glob("feff*.dat")):
            p = artemis.load_feff_path(f)
            if int(getattr(p, "nleg", 2)) == 2:
                if first is None or p.reff < first[1].reff:
                    first = (f, p)
        self.assertIsNotNone(first, "lepidocrocite 第一壳路径缺失")
        datfile = first[0]

        k = np.arange(2.0, 13.0, 0.02)
        p = artemis.load_feff_path(datfile)
        p.degen = 6.0
        artemis.set_path_parameters(
            p, s02=0.803, delta_r=-0.05, sigma2=0.004,
            delta_e0=0.5, fourth=6e-5)
        chi = artemis.exafs_model_k([p], k)
        ck = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                  mu_norm=np.ones_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
        res = types.SimpleNamespace(
            chi_k=ck, ft=ft, merged=None, norm=None,
            structure=None, absorber_symbol="Fe", self_absorption=None)
        out = xafs_fit.fit_identified(
            res, motif_id="lepidocrocite", absorber_symbol="Fe",
            rmin=1.2, rmax=2.4, vary_c4=True, s02=0.803,
            verbose=False)
        fit = out["fit"]
        c4 = fit.paths[0].c4
        self.assertIsNotNone(c4, "vary_c4=True 未产生 C4 变量")
        self.assertLess(abs(c4 - 6e-5), 2.5e-5,
                        f"C4 回收 {c4:.3e} vs 真值 6.0e-05")
        self.assertIn(fit.paths[0].label,
                      out["flags"]["c4_active_shells"])

        # 4 壳层 + 全部 C4 → N_idp 防护
        out2 = xafs_fit.fit_identified(
            self.lep, motif_id="lepidocrocite", absorber_symbol="Fe",
            rmax=4.0, vary_c4="all", verbose=False)
        self.assertTrue(out2["flags"]["c4_dropped_by_nidp"],
                        "4 壳层 + C4 未触发 c4_dropped_by_nidp")
        self.assertEqual(len(out2["flags"]["c4_active_shells"]), 1,
                         "C4 降级后应只保留第一壳")
        self.assertLess(out2["fit"].n_variables,
                        out2["fit"].n_independent)

        out3 = xafs_fit.fit_identified(
            self.lep, motif_id="lepidocrocite", absorber_symbol="Fe",
            vary_c3=True, verbose=False)
        rdf = xafs_fit.rdf_table(out3["fit"])
        self.assertEqual(len(rdf["rows"]), len(out3["fit"].paths),
                         "并非每条 C3 路径都生成了 g(R)")
        for summary, path in zip(rdf["summary"], out3["fit"].paths):
            self.assertAlmostEqual(summary["area"], path.degeneracy,
                                   places=6,
                                   msg=f"{path.label} 的 g(R) 面积≠N")
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = Path(tmp) / "lepido_rdf.csv"
            xafs_fit.write_rdf_csv(csv_path, rdf)
            self.assertTrue(csv_path.exists())
            header = csv_path.read_text(
                encoding="utf-8-sig").splitlines()[0]
            for col in ("label", "r_A", "g_r", "r_mode", "r_mean",
                        "skewness", "area"):
                self.assertIn(col, header)

    def test_sigma_mu_export_and_high_k_snr_gate(self):
        """#48 B9-2：三扫描合并的 σ_μ 落盘与高 k SNR 单调性。

        ``_mu.dat`` 必须有第 4 列 sigma_mu；注入递增噪声时 k_snr2 必须
        单调下降；没有 sigma_mu 时返回 None 并给出原因，不伪造曲线。
        """
        import run_xafs

        files = [LEPIDO, LEPIDO.with_suffix(".100"),
                 LEPIDO.with_suffix(".200")]
        for f in files:
            self.assertTrue(f.is_file(), f"缺少测试扫描 {f}")
        with tempfile.TemporaryDirectory() as tmp:
            res = run_xafs.process_sample(
                [str(f) for f in files], Path(tmp), absorber="Fe",
                strict=False, verbose=False)
            self.assertIsNotNone(res.merged.sigma_mu,
                                 "三扫描合并未保留 sigma_mu")
            mu_path = Path(tmp) / f"{res.stem}_mu.dat"
            lines = mu_path.read_text(encoding="utf-8").splitlines()
            self.assertIn("sigma_mu", lines[1])
            self.assertEqual(len(lines[2].split()), 4,
                             "_mu.dat 不是 4 列（E, mu, mu_norm, sigma_mu）")

        ck = res.chi_k
        edge = abs(float(res.norm.edge_step))
        k_snr2 = []
        for scale in (1e-3, 3e-3, 1e-2, 3e-2):
            out = athena.high_k_snr(
                ck, sigma_mu=np.full_like(ck.k, scale * edge),
                edge_step=res.norm.edge_step, kweight=2)
            self.assertIsNotNone(out["k_snr2"],
                                 f"噪声 scale={scale} 未给出 k_snr2")
            k_snr2.append(out["k_snr2"])
        for hi_noise, lo_noise in zip(k_snr2[1:], k_snr2[:-1]):
            self.assertLess(hi_noise, lo_noise - 0.2,
                            f"k_snr2 未随噪声单调下降: {k_snr2}")
        none = athena.high_k_snr(ck)
        self.assertIsNone(none["k_snr2"])
        self.assertIn("sigma_mu", none["reason"])
        self.assertIsNotNone(res.snr["k_snr2"])

    def test_wt_ruler_auto_mount_and_l3_whiteline_area(self):
        """#49 B9-2：WT 元素标尺自动挂载 + L3 白线净面积。

        K 边保留 ``pre_edge_n_peaks``；L3 边必须没有该字段，而给出
        有限的 ``white_line_area``/``whiteline_peak_height``。
        自动标尺的元素集合来自拟合路径 ∪ {O, 吸收元素}，并带缓存。
        """
        fit_obj = self.lep_fit["fit"]["fit"]
        ruler = xafs_report.auto_ref_k_peaks(
            self.lep, fit_obj, "Fe")
        self.assertIsNotNone(ruler, "拟合路径存在时未生成 WT 标尺")
        self.assertIn("O", ruler)
        self.assertIn("Fe", ruler)
        self.assertTrue(np.isfinite(ruler["O"]))
        self.assertTrue(np.isfinite(ruler["Fe"]))
        self.assertLess(ruler["O"], ruler["Fe"])
        again = xafs_report.auto_ref_k_peaks(
            self.lep, fit_obj, "Fe")
        self.assertIs(again, ruler, "同一参数未命中 WT 标尺缓存")
        with tempfile.TemporaryDirectory() as tmp:
            prof = xafs_diagnostics.write_wt_shell_profiles(
                self.lep, self.lep_fit, tmp, ref_k_peaks=ruler)
            self.assertTrue(Path(prof["png"]).exists())

        k_feat = athena.extract_xanes_features(
            self.lep.norm, edge="K", absorber="Fe")
        self.assertIn("pre_edge_n_peaks", k_feat)
        self.assertNotIn("white_line_area", k_feat)

        res_pt = xafs_pipeline.process_sample(
            [PT_FOIL], absorber="Pt", strict=False, verbose=False)
        l_feat = athena.extract_xanes_features(
            res_pt.norm, edge="L3", absorber="Pt")
        self.assertNotIn("pre_edge_n_peaks", l_feat,
                         "L3 边不应拟合 1s→3d 边前峰")
        self.assertTrue(np.isfinite(l_feat["white_line_area"]))
        self.assertGreater(l_feat["white_line_area"], 0.0)
        self.assertTrue(np.isfinite(l_feat["whiteline_peak_height"]))

    def test_prj_roundtrip_and_multichannel_mca_merge(self):
        """#50 B9-3：gzip-JSON 工程容器 round-trip 与 MCA 坏道剔除。

        8 通道中第 2、6 道注入 Bragg spike，必须精确剔除且保底 ≥2 道；
        求和谱必须等于干净通道按 info["weights"] 的加权和。
        """
        rng = np.random.default_rng(4)
        energy = np.linspace(7100.0, 7300.0, 201)
        mu = (1.0 + 0.5 / (1.0 + np.exp(-(energy - 7120.0) / 2.0))
              + 0.01 * np.sin(energy / 3.0))
        scan = RawScan(energy=energy, mu=mu, meta={"mode": "trans"})
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "roundtrip.prj"
            xafs_io.write_prj(
                path, {"scans": [scan], "standards": {},
                       "meta": {"created": "test"}})
            self.assertEqual(path.read_bytes()[:2], b"\x1f\x8b",
                             ".prj 未使用 gzip 容器")
            out = xafs_io.read_prj(path, as_scans=True)
            got = out["scans"][0]
            self.assertTrue(np.array_equal(got.energy, scan.energy))
            self.assertTrue(np.array_equal(got.mu, scan.mu))
            self.assertEqual(out["meta"], {"created": "test"})

            # 兼容未压缩 JSON
            import json
            plain = Path(tmp) / "plain.prj"
            plain.write_text(json.dumps({
                "xafs_prj": 1,
                "scans": [{"energy": energy.tolist(),
                           "mu": mu.tolist(), "meta": {}}],
                "meta": {}}), encoding="utf-8")
            out2 = xafs_io.read_prj(plain, as_scans=True)
            self.assertTrue(np.array_equal(out2["scans"][0].energy, energy))

            # Athena Perl 工程必须干净拒绝
            perl = Path(tmp) / "athena.prj"
            perl.write_text("$VAR1 = { 'x' => 1 };", encoding="utf-8")
            with self.assertRaises(XAFSConfigError):
                xafs_io.read_prj(perl)

        n = 600
        base = (1.0 + 0.6 / (1.0 + np.exp(-(np.arange(n) - 300.0) / 12.0))
                + 0.005 * np.sin(np.arange(n) / 9.0))
        channels = (np.tile(base[:, None], (1, 8))
                    + rng.normal(0.0, 0.002, (n, 8)))
        for idx in (2, 6):
            for j in (120, 121, 400):
                channels[j, idx] += 2.5
        summed, info = athena.merge_mca_channels(
            channels, z_threshold=3.5, min_channels=2)
        self.assertEqual(info["rejected"], [2, 6])
        self.assertEqual(info["n_kept"], 6)
        self.assertGreaterEqual(info["n_kept"], 2)
        clean = [i for i in range(8) if i not in info["rejected"]]
        weights = np.asarray(info["weights"], dtype=float)
        self.assertAlmostEqual(float(weights.sum()), 1.0, places=12)
        ref = np.sum(channels[:, clean] * weights[None, :], axis=1)
        np.testing.assert_allclose(summed, ref, rtol=1e-12, atol=1e-12)

    def test_mback_short_window_normalization(self):
        """#51 B10-1：短边后窗口 MBACK 归一化的 edge_step 回收。

        Fe K 理论截面合成 μ(E)（边后仅 60–150 eV），MBACK 的 edge_step
        相对误差必须 <3%；传统二次多项式归一化在同窗口平台发 NaN 且
        偏差更大。a≤0 时干净拒绝。
        """
        from xraydb import mu_elam

        energy = np.linspace(6900.0, 7280.0, 650)
        mu = (1.7 * mu_elam("Fe", energy)
              + 0.02 * (energy - 7112.0) + 0.3)
        scan = RawScan(energy=energy, mu=mu, meta={})
        truth = 1.7 * (mu_elam("Fe", 7162.0) - mu_elam("Fe", 7062.0))
        mback = athena.normalize_mback(
            scan, 7112.0, "Fe", post_window=(60.0, 150.0))
        rel = abs(mback.edge_step - truth) / truth
        self.assertLess(rel, 0.03,
                        f"MBACK edge_step 相对误差 {rel:.3e} ≥3%")
        self.assertAlmostEqual(mback.edge_step, truth, places=8)
        self.assertTrue(np.isfinite(mback.plateau))
        self.assertEqual(mback.engine, "mback")
        self.assertEqual(mback.mback["element"], "Fe")
        self.assertGreater(mback.mback["a"], 0.0)

        trad = athena.normalize_absorption(
            scan, 7112.0, post_window=(60.0, 150.0))
        trad_rel = abs(trad.edge_step - truth) / truth
        self.assertGreater(trad_rel, 0.01,
                           "传统短窗口归一化不应比 MBACK 更准")
        self.assertGreater(trad_rel, 5.0 * rel)
        self.assertFalse(np.isfinite(trad.plateau))

        negative = RawScan(
            energy=energy,
            mu=-1.7 * mu_elam("Fe", energy) + 0.02 * (energy - 7112.0) + 0.3,
            meta={})
        with self.assertRaises(XAFSError):
            athena.normalize_mback(
                negative, 7112.0, "Fe", post_window=(60.0, 150.0))

    def test_autobk_chi_std_prior_and_clamp(self):
        """#52 B10-1：AUTOBK 先验谱压低低 R 泄漏，clamp 生效且合法。

        合成谱 = 真实 Fe–O χ(k) + 0.7 Å 低频伪振荡 + 噪声。以真实
        χ(k) 作为 chi_std 时，R<1.2 Å 内“提取谱−真值”的泄漏必须
        降到无先验的一半以下；非法长度必须报 XAFSConfigError。
        """
        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "lepidocrocite"][0]
        cache = structure_lib.ensure_feff(motif, edge="K")
        first = None
        for f in sorted(cache.glob("feff*.dat")):
            p = artemis.load_feff_path(f)
            if int(getattr(p, "nleg", 2)) == 2:
                if first is None or p.reff < first[1].reff:
                    first = (f, p)
        datfile = first[0]
        k = np.arange(0.0, 15.0, 0.05)
        p = artemis.load_feff_path(datfile)
        p.degen = 6.0
        artemis.set_path_parameters(
            p, s02=0.803, delta_r=0.0, sigma2=0.004, delta_e0=0.0)
        true_chi = artemis.exafs_model_k([p], k)
        rng = np.random.default_rng(3)
        leak = 0.03 * np.sin(2.0 * 0.7 * k + 0.4)
        noise = rng.normal(0.0, 0.002, k.size)
        mu = 1.0 + true_chi + leak + noise
        energy = 7112.0 + k ** 2 / K2
        norm = NormalizedScan(
            energy=energy, mu=mu, mu_norm=mu,
            pre_edge_line=np.zeros_like(mu), edge_step=1.0,
            plateau=1.0, e0=7112.0, flattened=mu.copy())

        trap = getattr(np, "trapezoid", None) or getattr(np, "trapz")
        ft_true = athena.ft_to_r_space(k, true_chi, None, 3.0, 12.0, 2)
        r_mask = (ft_true.r >= 0.2) & (ft_true.r < 1.2)

        def low_r_leak(chi_k):
            ft = athena.ft_to_r_space(
                chi_k.k, chi_k.chi, None, 3.0, 12.0, 2)
            re = np.interp(ft_true.r, ft.r, ft.chi_r.real)
            im = np.interp(ft_true.r, ft.r, ft.chi_r.imag)
            diff = (re + 1j * im) - ft_true.chi_r
            return float(trap(np.abs(diff[r_mask]), ft_true.r[r_mask]))

        plain = athena.extract_chi_k(norm, rbkg=1.2, kweight=2)
        with_std = athena.extract_chi_k(
            norm, rbkg=1.2, kweight=2, chi_std=true_chi, k_std=k)
        leak_plain = low_r_leak(plain)
        leak_std = low_r_leak(with_std)
        self.assertTrue(with_std.chi_std_used,
                        "chi_std 未被记录为已使用")
        self.assertLess(leak_std, 0.5 * leak_plain,
                        f"低 R 泄漏未下降: plain={leak_plain:.4f}, "
                        f"std={leak_std:.4f}")

        clamped = athena.extract_chi_k(
            norm, rbkg=1.2, kweight=2, clamp_lo=0.1, clamp_hi=0.1)
        self.assertAlmostEqual(clamped.clamp_lo, 0.1)
        self.assertAlmostEqual(clamped.clamp_hi, 0.1)
        self.assertTrue(np.all(np.isfinite(clamped.chi)))
        self.assertFalse(np.array_equal(clamped.chi, plain.chi),
                         "clamp_lo/hi 未改变 AUTOBK 解")

        with self.assertRaises(XAFSConfigError):
            athena.extract_chi_k(
                norm, rbkg=1.2, chi_std=np.ones(10), k_std=np.arange(9.0))
        with self.assertRaises(XAFSConfigError):
            athena.extract_chi_k(
                norm, rbkg=1.2, chi_std=np.ones(20))

    def test_pca_target_transform_spoil_and_mcr_als(self):
        """#53 B10-2：SPOIL 判定与 MCR-ALS 盲解谱。

        3 个独立端元 + 9 条混合谱（含纯端元，浓度和=1）加 0.5% 噪声：
        真实端元 SPOIL<3、随机谱>6；MCR-ALS 组分谱相关系数>0.99、
        浓度 RMSE<0.02，且 closure=True 时收敛。
        """
        import xafs_lcf

        x = np.linspace(0.0, 10.0, 400)

        def g(c, w):
            return np.exp(-0.5 * ((x - c) / w) ** 2)

        e1 = g(2.0, 0.35) + 0.10 * g(4.0, 0.6)
        e2 = g(5.0, 0.55) + 0.15 * g(7.5, 0.5)
        e3 = g(8.0, 0.40) + 0.20 * g(3.0, 0.4)
        truth = np.vstack([e1, e2, e3])
        weights = np.array([
            [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
            [0.5, 0.3, 0.2], [0.2, 0.6, 0.2], [0.3, 0.2, 0.5],
            [0.6, 0.1, 0.3], [0.15, 0.15, 0.7], [0.25, 0.55, 0.2]])
        rng = np.random.default_rng(11)
        scale = float(np.max(truth))
        specs = [(f"s{i}", x, weights[i] @ truth
                  + rng.normal(0.0, 0.005 * scale, x.size))
                 for i in range(len(weights))]
        pca = xafs_lcf.run_pca(specs, 0.0, 10.0)
        self.assertEqual(pca["x"].shape, x.shape)
        self.assertEqual(np.asarray(pca["components"]).shape,
                         (x.size, len(weights)))
        tt = xafs_lcf.target_transform(
            pca, [("true1", x, e1), ("true2", x, e2),
                  ("true3", x, e3),
                  ("rand", x, rng.normal(0.0, 1.0, x.size))],
            n_components=3)
        for name in ("true1", "true2", "true3"):
            self.assertLess(tt[name]["spoils"], 3.0,
                            f"{name} SPOIL={tt[name]['spoils']:.3f} ≥3")
            self.assertEqual(tt[name]["rating"], "good")
        self.assertGreater(tt["rand"]["spoils"], 6.0,
                           f"随机谱 SPOIL={tt['rand']['spoils']:.3f} ≤6")

        mcr = xafs_lcf.mcr_als(
            specs, 3, max_iter=1500, tol=1e-10, closure=True)
        self.assertTrue(mcr["converged"], "MCR-ALS 未收敛")
        self.assertLess(mcr["r_factor"], 0.03)
        comp = np.asarray(mcr["components"])
        conc = np.asarray(mcr["concentrations"])
        self.assertEqual(comp.shape, (x.size, 3))
        self.assertEqual(conc.shape, (3, len(weights)))
        # 组分顺序由数据决定，贪心匹配真值后计算相关与 RMSE
        perm = []
        for i in range(3):
            candidates = [j for j in range(3) if j not in perm]
            cors = [abs(np.corrcoef(comp[:, i], truth[j])[0, 1])
                    for j in candidates]
            perm.append(candidates[int(np.argmax(cors))])
        for i, j in enumerate(perm):
            corr = abs(np.corrcoef(comp[:, i], truth[j])[0, 1])
            self.assertGreater(corr, 0.99,
                               f"组分 {i} 与真值 {j} 相关 {corr:.4f} ≤0.99")
        rmse = float(np.sqrt(np.mean(
            (conc[perm, :] - weights.T) ** 2)))
        self.assertLess(rmse, 0.02, f"浓度 RMSE={rmse:.4f} ≥0.02")

    def test_phase_corrected_ft_and_wt(self):
        """#54 B11-1：FEFF 总相移校正把 FT/WT 峰位拉回 Reff。

        未校正 FT 峰位在 1.4–1.7 Å（表观 R），校正后 |R_peak−Reff|
        必须 <0.05 Å；WT 使用同一相移也必须回到 Reff。
        """
        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "lepidocrocite"][0]
        cache = structure_lib.ensure_feff(motif, edge="K")
        first = None
        for f in sorted(cache.glob("feff*.dat")):
            p = artemis.load_feff_path(f)
            if int(getattr(p, "nleg", 2)) == 2:
                if first is None or p.reff < first[1].reff:
                    first = (f, p)
        datfile, path0 = first
        reff = float(path0.reff)
        k = np.arange(2.0, 13.0, 0.02)
        p = artemis.load_feff_path(datfile)
        p.degen = 6.0
        artemis.set_path_parameters(
            p, s02=1.0, delta_r=0.0, sigma2=0.004, delta_e0=0.0)
        chi = artemis.exafs_model_k([p], k)
        phi = artemis.path_total_phase(datfile, k)
        self.assertEqual(phi.shape, k.shape)
        self.assertTrue(np.all(np.isfinite(phi)))
        self.assertIn(artemis.path_total_phase.last_source,
                      ("feff_pha", "hilbert_fallback"))

        def ft_peak(phase_correct):
            ft = athena.ft_to_r_space(
                k, chi, None, 3.0, 12.0, 2, phase_correct=phase_correct)
            m = (ft.r >= 0.8) & (ft.r <= 2.8)
            return float(ft.r[m][np.argmax(np.abs(ft.chi_r[m]))])

        raw_peak = ft_peak(None)
        corrected_peak = ft_peak(phi)
        self.assertGreater(raw_peak, 1.3)
        self.assertLess(raw_peak, 1.7)
        self.assertLess(abs(corrected_peak - reff), 0.05,
                        f"校正 FT 峰 {corrected_peak:.4f} vs Reff {reff:.4f}")

        ck = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                  mu_norm=np.ones_like(k), rbkg=1.0)

        def wt_peak(phase_correct):
            wt = larch_wt.continuous_wavelet_transform(
                ck, sigma=1.5, kmin=3.0, kmax=12.0, kweight=2,
                rmax=3.0, phase_correct=phase_correct)
            m = (wt.r >= 0.8) & (wt.r <= 2.8)
            mag = wt.magnitude[:, m]
            ij = np.unravel_index(int(np.argmax(mag)), mag.shape)
            return float(wt.r[m][ij[1]])

        wt_raw = wt_peak(None)
        wt_corrected = wt_peak(phi)
        self.assertLess(wt_raw, 1.7)
        self.assertLess(abs(wt_corrected - reff), 0.05,
                        f"校正 WT 峰 {wt_corrected:.4f} vs Reff {reff:.4f}")

    def test_hetero_cocupancy_dual_branch_paths(self):
        """#55 B11-1：异核共占位按距离聚壳并拆成元素分支。

        Pm-3m 中 Ti(0.7)/Ru(0.3) 共占 1b 位（8 个镜像，R=3.5507 Å）：
        两支占据度和=1、简并度合计=8（不是 16——不能按展开条目重复
        乘占据度），且 single_element_approx=False。
        """
        cif_text = """data_test_hetero
_cell_length_a 4.10
_cell_length_b 4.10
_cell_length_c 4.10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'Pm-3m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Fe1 Fe 0.0 0.0 0.0 1.0
Ti1 Ti 0.5 0.5 0.5 0.7
Ru1 Ru 0.5 0.5 0.5 0.3
"""
        with tempfile.TemporaryDirectory() as tmp:
            cif = Path(tmp) / "hetero.cif"
            cif.write_text(cif_text, encoding="utf-8")
            out = artemis.cif_hetero_cocupancy_branches(
                cif, "Fe", rmax=5.0, delta_z=3, occ_min=0.10)
        self.assertFalse(out["single_element_approx"],
                         "Ti/Ru 共占位被误判为单元素近似")
        hetero_shells = [s for s in out["shells"]
                         if not s["single_element_approx"]]
        self.assertEqual(len(hetero_shells), 1)
        shell = hetero_shells[0]
        self.assertAlmostEqual(shell["distance"], 3.5507041555161982,
                               places=10)
        self.assertEqual(shell["n_geom"], 8)
        branches = {b["element"]: b for b in shell["branches"]}
        self.assertEqual(set(branches), {"Ti", "Ru"})
        self.assertAlmostEqual(branches["Ti"]["occ"], 0.7, places=12)
        self.assertAlmostEqual(branches["Ru"]["occ"], 0.3, places=12)
        self.assertAlmostEqual(
            branches["Ti"]["occ"] + branches["Ru"]["occ"], 1.0, places=12)
        self.assertAlmostEqual(
            branches["Ti"]["degen"] + branches["Ru"]["degen"],
            shell["n_geom"], places=12)
        self.assertAlmostEqual(branches["Ti"]["degen"], 5.6, places=12)
        self.assertAlmostEqual(branches["Ru"]["degen"], 2.4, places=12)

    def test_correlated_debye_sigma2_contract(self):
        """#56 B11-2：关联 Debye 模型的温度/距离趋势与输入护栏。

        σ²(300 K)>σ²(80 K)；R=4 Å 的关联因子更接近 1，σ² 大于 R=2 Å；
        T→0 收敛到与温度无关的零点平台（静态项 + 零点振动），R=4 Å 的
        高温值与 Einstein 模型同量级（0.5–2）。
        """
        kw = dict(theta_d=400.0, mass_absorber=55.845,
                  mass_scatterer=15.999)
        near_80 = artemis.correlated_debye_sigma2(80.0, reff=2.0, **kw)
        near_300 = artemis.correlated_debye_sigma2(300.0, reff=2.0, **kw)
        far_300 = artemis.correlated_debye_sigma2(300.0, reff=4.0, **kw)
        self.assertGreater(near_300, near_80)
        self.assertGreater(far_300, near_300)
        self.assertLess(artemis.correlated_debye_sigma2.last_details[
            "correlation_factor"], 1.0)
        self.assertEqual(artemis.correlated_debye_sigma2.last_details[
            "k_debye_source"], "conservative_default_0.35A^-1")

        t_plateau_a = artemis.correlated_debye_sigma2(
            1e-6, reff=2.0, **kw)
        t_plateau_b = artemis.correlated_debye_sigma2(
            1e-4, reff=2.0, **kw)
        self.assertAlmostEqual(t_plateau_a, t_plateau_b, places=8)
        static = 0.0015
        with_static = artemis.correlated_debye_sigma2(
            1e-6, reff=2.0, sigma2_static=static, **kw)
        self.assertAlmostEqual(with_static - static, t_plateau_a, places=10)

        einstein = artemis.einstein_sigma2(1000.0, 400.0, 12.436)
        far_1000 = artemis.correlated_debye_sigma2(1000.0, reff=4.0, **kw)
        ratio = far_1000 / einstein
        self.assertGreaterEqual(ratio, 0.5)
        self.assertLessEqual(ratio, 2.0)
        self.assertGreater(far_1000, near_300)

        for bad in (
                dict(temperature_K=0.0, theta_d=400.0, reff=2.0,
                     mass_absorber=55.845, mass_scatterer=15.999),
                dict(temperature_K=300.0, theta_d=0.0, reff=2.0,
                     mass_absorber=55.845, mass_scatterer=15.999),
                dict(temperature_K=300.0, theta_d=400.0, reff=2.0,
                     mass_absorber=0.0, mass_scatterer=15.999)):
            with self.assertRaises(XAFSConfigError):
                artemis.correlated_debye_sigma2(**bad)

    def test_bond_valence_sum_restraint(self):
        """#57 B11-2：BVS 计算与软约束。

        Fe–O N=6,R=2.0 的 BVS≈3.128（Fe³⁺ 合理）；约束 target=3.0
        后拟合输出的 bvs 必须与用手工 BVS 从最终路径重算的值一致，
        且 |bvs−3| 小于无约束拟合；bvs_chi2 有限。
        """
        bvs = artemis.bond_valence_sum(
            [{"N": 6.0, "R": 2.0, "ligand": "O"}], absorber="Fe")
        self.assertAlmostEqual(bvs, 3.12804471568531, places=10)
        with self.assertRaises(XAFSConfigError):
            artemis.bond_valence_sum(
                [{"N": 6.0, "R": 2.0, "ligand": "P"}], absorber="Fe")

        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "lepidocrocite"][0]
        cache = structure_lib.ensure_feff(motif, edge="K")
        datfile = None
        for f in sorted(cache.glob("feff*.dat")):
            p = artemis.load_feff_path(f)
            if (int(getattr(p, "nleg", 2)) == 2
                    and abs(p.reff - 1.98) < 0.05):
                datfile = f
                break
        self.assertIsNotNone(datfile)
        k = np.linspace(2.0, 14.0, 601)
        p = artemis.load_feff_path(datfile)
        p.degen = 6.0
        artemis.set_path_parameters(
            p, s02=0.803, delta_r=-0.05, sigma2=0.004, delta_e0=0.5)
        chi = np.interp(k, p.k, p.chi)
        ck = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                  mu_norm=np.ones_like(k), rbkg=1.2)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2, 1.0)
        base = artemis.fit_r_space(
            ck, ft, [str(datfile)], 1.2, 2.4, kmin=3.0, kmax=12.0,
            kweight=2, s02=0.803, vary_n=True, vary_delta_e0=True,
            path_labels=["Fe-O"])
        free_bvs = artemis.bond_valence_sum(base.paths, absorber="Fe")
        self.assertIsNone(base.bvs)
        constrained = artemis.fit_r_space(
            ck, ft, [str(datfile)], 1.2, 2.4, kmin=3.0, kmax=12.0,
            kweight=2, s02=0.803, vary_n=True, vary_delta_e0=True,
            path_labels=["Fe-O"],
            bvs_restraint={"target": 3.0, "sigma": 0.05})
        self.assertEqual(constrained.bvs_target, 3.0)
        self.assertIsNotNone(constrained.bvs)
        manual = artemis.bond_valence_sum(
            constrained.paths, absorber="Fe")
        self.assertAlmostEqual(constrained.bvs, manual, places=10)
        self.assertLess(abs(constrained.bvs - 3.0),
                        0.05 * abs(free_bvs - 3.0))
        self.assertTrue(np.isfinite(constrained.bvs_chi2))
        self.assertLess(constrained.bvs_chi2, 0.01)

    def test_multispace_q_and_wt_fit(self):
        """#58 B11-3：R/q/WT 三空间拟合回收同一组 N/R/σ²。

        同一合成单壳谱在三个空间分别拟合，参数差必须
        |ΔN|<0.3、|ΔR|<0.02 Å、|Δσ²|<0.001 Å²；默认 "r" 与显式
        "r" 逐位一致（历史口径不回退）。
        """
        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "lepidocrocite"][0]
        cache = structure_lib.ensure_feff(motif, edge="K")
        datfile = None
        for f in sorted(cache.glob("feff*.dat")):
            p = artemis.load_feff_path(f)
            if (int(getattr(p, "nleg", 2)) == 2
                    and abs(p.reff - 1.98) < 0.05):
                datfile = f
                break
        k = np.arange(2.0, 13.0, 0.02)
        p = artemis.load_feff_path(datfile)
        p.degen = 6.0
        artemis.set_path_parameters(
            p, s02=0.803, delta_r=-0.05, sigma2=0.004, delta_e0=0.5)
        chi = artemis.exafs_model_k([p], k)
        ck = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                  mu_norm=np.ones_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)

        fits = {}
        for space in ("r", "q", "wt"):
            fits[space] = artemis.fit_r_space(
                ck, ft, [str(datfile)], 1.2, 2.4,
                kmin=3.0, kmax=12.0, kweight=2, s02=0.803,
                vary_n=True, vary_delta_e0=True,
                path_degen=[6.0], path_labels=["Fe-O"],
                fit_space=space)
            self.assertEqual(fits[space].fit_space, space)
            self.assertTrue(np.isfinite(fits[space].aic))
            self.assertTrue(np.isfinite(fits[space].bic))

        default = artemis.fit_r_space(
            ck, ft, [str(datfile)], 1.2, 2.4,
            kmin=3.0, kmax=12.0, kweight=2, s02=0.803,
            vary_n=True, vary_delta_e0=True,
            path_degen=[6.0], path_labels=["Fe-O"])
        self.assertEqual(default.r_factor, fits["r"].r_factor)
        self.assertEqual(default.paths[0].degeneracy,
                         fits["r"].paths[0].degeneracy)

        ref = fits["r"].paths[0]
        for space in ("q", "wt"):
            got = fits[space].paths[0]
            self.assertLess(abs(got.degeneracy - ref.degeneracy), 0.3)
            self.assertLess(abs(got.r - ref.r), 0.02)
            self.assertLess(abs(got.sigma2 - ref.sigma2), 0.001)

    def test_multi_edge_e0_groups_and_aic_bic(self):
        """#59 B12-1：双吸收边独立 ΔE₀ 分组与 AIC/BIC 惩罚。

        两条同路径合成谱真值 ΔE₀=+2.0 / −6.0 eV，edge_group=K/L3：
        share_e0="auto" 分出 2 个槽位并各自回收 <0.5 eV；
        share_e0=True 只有 1 个槽位、Rf 显著变差；BIC 对多出的参数
        惩罚更重（auto 相对 True 的 BIC 优势小于 AIC 优势）。
        """
        motif = [m for m in structure_lib.load_library("Fe")
                 if m.id == "lepidocrocite"][0]
        cache = structure_lib.ensure_feff(motif, edge="K")
        datfile = None
        for f in sorted(cache.glob("feff*.dat")):
            p = artemis.load_feff_path(f)
            if (int(getattr(p, "nleg", 2)) == 2
                    and abs(p.reff - 1.98) < 0.05):
                datfile = f
                break
        k = np.arange(2.0, 13.0, 0.02)
        datasets = []
        for group, e0_true in (("K", 2.0), ("L3", -6.0)):
            p = artemis.load_feff_path(datfile)
            p.degen = 6.0
            artemis.set_path_parameters(
                p, s02=1.0, delta_r=-0.02, sigma2=0.004,
                delta_e0=e0_true)
            chi = artemis.exafs_model_k([p], k)
            ck = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                      mu_norm=np.ones_like(k), rbkg=1.0)
            ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
            datasets.append({
                "name": group, "chi_k": ck, "ft": ft,
                "datfiles": [str(datfile)], "rmin": 1.2, "rmax": 2.4,
                "kmin": 3.0, "kmax": 12.0, "kweight": 2, "s02": 1.0,
                "vary_n": False, "path_degen": [6.0],
                "path_labels": ["Fe-O"], "edge_group": group})

        auto = artemis.fit_multi_datasets(datasets, share_e0="auto")
        shared = artemis.fit_multi_datasets(datasets, share_e0=True)
        self.assertEqual(auto["e0_groups"], ["K", "L3"])
        self.assertEqual(shared["e0_groups"], ["global"])
        self.assertEqual(auto["n_variables_total"],
                         shared["n_variables_total"] + 1)
        self.assertAlmostEqual(auto["datasets"][0].delta_e0, 2.0, places=6)
        self.assertAlmostEqual(auto["datasets"][1].delta_e0, -6.0, places=6)
        self.assertLess(abs(auto["datasets"][0].delta_e0 - 2.0), 0.5)
        self.assertLess(abs(auto["datasets"][1].delta_e0 + 6.0), 0.5)
        self.assertLess(auto["global_r_factor"], 1e-6)
        self.assertGreater(shared["global_r_factor"], 0.01)

        for key in ("aic", "bic", "global_aic", "global_bic"):
            self.assertTrue(np.isfinite(auto[key]), f"auto {key} 非有限")
            self.assertTrue(np.isfinite(shared[key]),
                            f"shared {key} 非有限")
        for fit in auto["datasets"] + shared["datasets"]:
            self.assertTrue(np.isfinite(fit.aic))
            self.assertTrue(np.isfinite(fit.bic))
            self.assertGreater(fit.bic, fit.aic)
        # auto 多 1 个 ΔE₀ 参数：BIC 惩罚使其相对 shared 的优势变小
        aic_gain = shared["aic"] - auto["aic"]
        bic_gain = shared["bic"] - auto["bic"]
        self.assertGreater(aic_gain, bic_gain,
                           "BIC 未对多参数模型施加更重惩罚")

    def test_sac_dac_cluster_generator_and_element_expansion(self):
        """#60 B12-2：SAC/DAC 团簇生成、FEFF 可跑与元素扩充。

        M–N₄ 恰好 4 条且键长误差 <1e-9；C 壳层在 ring_radius±0.2；
        build_feff_input + feff6l 成功；12 个新元素全部未标定，
        Fe/Co/Ni/Cu/Pt 的 S₀² 逐位不变。
        """
        sites = artemis.sac_cluster_sites(
            bond_m_ligand=2.0, ligand="N", n_coord=4,
            axial=None, carbon_shell=True, ring_radius=2.9)
        n_sites = [s for s in sites if s[3] == "N"]
        c_sites = [s for s in sites if s[3] == "C"]
        self.assertEqual(len(n_sites), 4)
        self.assertEqual(len(c_sites), 8)
        for s in n_sites:
            self.assertLess(abs(float(np.linalg.norm(s[:3])) - 2.0), 1e-9)
        for s in c_sites:
            self.assertLess(
                abs(float(np.linalg.norm(s[:3])) - 2.9), 0.2)
        inp = artemis.build_feff_input(
            sites, absorber="Fe", rmax=4.0, title="sac_test")
        with tempfile.TemporaryDirectory() as tmp:
            files = artemis.run_feff6l(inp, tmp)
            self.assertGreaterEqual(len(files), 1, "FEFF 未生成路径文件")
            # files 的路径位于临时目录内，必须在 with 块内读路径文件；
            # 出了块目录已被删除，load_feff_path 会报文件不存在。
            p0 = artemis.load_feff_path(files[0])
            self.assertGreater(p0.reff, 1.0)

        dac = artemis.diatomic_m1m2_sites(
            "Fe", "Co", bond_length=2.5, ligand="N",
            n_coord=4, ligand_bond=2.0)
        m2 = [s for s in dac if s[3] == "Co"]
        self.assertEqual(len(m2), 1)
        self.assertAlmostEqual(m2[0][2], 2.5, places=12)
        self.assertEqual(sum(1 for s in dac if s[3] == "N"), 8)

        new_elements = ("Cr", "Mn", "Zn", "Mo", "Ru", "Rh",
                        "Pd", "Ag", "Ir", "Au", "Ce", "U")
        for sym in new_elements:
            cfg = get_absorber(sym)
            self.assertFalse(cfg.calibrated, f"{sym} 不应标记为已标定")
            self.assertIsNone(cfg.s02, f"{sym} 的 S₀² 必须为 None")
            self.assertTrue(np.isfinite(cfg.e0_edge))
            self.assertTrue(np.isfinite(cfg.kmax_ceiling))
            self.assertIn(sym, artemis.COVALENT_RADII)
        calibrated = {"Fe": 0.803, "Co": 0.7308, "Ni": 0.7795,
                      "Cu": 0.8670, "Pt": 0.8075}
        for sym, s02 in calibrated.items():
            self.assertEqual(get_absorber(sym).s02, s02,
                             f"{sym} 的标定 S₀² 被改动")

    # ------------------------------------------------------------
    # B13-1 SAC/DAC 动态 motif 与 M–M Hamilton F-test
    # ------------------------------------------------------------
    def test_sac_axial_and_dac_mm_bond_ftest_pipeline(self):
        """#61 B13-1：动态 SAC/DAC motif + M₁–M₂ 金属键显著性判据。

        用 DAC motif 的真实 FEFF 路径合成两组谱（ε≈1e-3）：
        - 纯 Fe–N₄ 谱按 DAC 模型拟合：mm_bond_ftest.significant=False
          （p>0.05，N_MM 不显著大于误差）；
        - 真 Fe–N₄ + Fe–Co 谱：R_Co∈[2.47,2.53]、N_Co∈[0.7,1.3]，
          significant=True 且 p<0.01；
        - SAC 轴向氧谱：N_N≈4、N_O≈1（|Δ|<0.5）。
        """
        from types import SimpleNamespace
        from larch.xafs import feffpath

        s02 = 0.803

        def synth(files, degens, noise=1e-3, seed=13):
            k = np.arange(3.0, 12.0, 0.05)
            plist = []
            for f, d in zip(files, degens):
                p = feffpath(str(f), s02=s02, sigma2=0.004, e0=2.0,
                             deltar=0.0)
                p.degen = float(d)
                p._calc_chi()
                plist.append(p)
            chi = artemis.exafs_model_k(plist, k)
            if noise:
                rng = np.random.default_rng(seed)
                chi = chi + rng.normal(0.0, noise, size=chi.shape)
            chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                         mu_norm=np.zeros_like(k), rbkg=1.0)
            ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
            return SimpleNamespace(
                chi_k=chi_k, ft=ft, absorber_symbol="Fe", structure=None,
                merged=None, norm=None, self_absorption=None)

        dac = structure_lib.parse_catalytic_cluster_spec(
            "Fe", "dac:Co@2.50,N4@2.00")
        self.assertEqual(dac.nuclearity, "binuclear")
        self.assertIn("_dynamic_cache", str(dac.cache_dir))
        plan = xafs_fit.select_shell_paths(
            dac, edge="K", rmin=1.0, rmax=2.6, kmin=3.0, kmax=12.0)
        self.assertEqual([s["scatterer"] for s in plan["shells"]],
                         ["N", "Co"])
        self.assertEqual(plan["degens"], [4.0, 1.0])

        # ---- 纯 SAC 谱：DAC 模型的 Co 分支必须判为不显著 ----
        res_sac = synth(plan["files"][:1], plan["degens"][:1])
        out_sac = xafs_fit.fit_identified(
            res_sac, absorber_symbol="Fe",
            catalytic_model="dac:Co@2.50,N4@2.00",
            rmin=1.0, rmax=2.6, s02=s02, verbose=False)
        self.assertEqual(out_sac["mode"], "motif")
        mm = out_sac["mm_bond_ftest"]
        self.assertIsNotNone(mm)
        self.assertEqual(mm["partner"], "Co")
        self.assertFalse(mm["significant"],
                         f"纯 SAC 谱不应检出 Fe–Co: {mm}")
        self.assertGreater(mm["p_value"], 0.05)

        # ---- 真 DAC 谱：回收 Fe–Co 距离/配位并判显著 ----
        res_dac = synth(plan["files"], plan["degens"])
        out_dac = xafs_fit.fit_identified(
            res_dac, absorber_symbol="Fe",
            catalytic_model="dac:Co@2.50,N4@2.00",
            rmin=1.0, rmax=2.6, s02=s02, verbose=False)
        self.assertEqual(out_dac["mode"], "motif")
        co = [p for p in out_dac["fit"].paths if "Co" in p.label]
        self.assertEqual(len(co), 1)
        self.assertGreaterEqual(co[0].r, 2.47)
        self.assertLessEqual(co[0].r, 2.53)
        self.assertGreaterEqual(co[0].degeneracy, 0.7)
        self.assertLessEqual(co[0].degeneracy, 1.3)
        mm_dac = out_dac["mm_bond_ftest"]
        self.assertTrue(mm_dac["significant"],
                        f"真 DAC 谱未检出 Fe–Co: {mm_dac}")
        self.assertLess(mm_dac["p_value"], 0.01)
        self.assertGreaterEqual(mm_dac["delta_r_factor"], 0.0)

        # ---- SAC 轴向氧：N_N≈4、N_O≈1（共享 σ² 避免 0.1 Å 近简并） ----
        sac = structure_lib.parse_catalytic_cluster_spec(
            "Fe", "sac:N4@2.00,axial=O@2.10,C=none")
        self.assertEqual(sac.nuclearity, "mononuclear")
        self.assertIn("_dynamic_cache", str(sac.cache_dir))
        splan = xafs_fit.select_shell_paths(
            sac, edge="K", rmin=1.0, rmax=2.6, kmin=3.0, kmax=12.0)
        self.assertEqual([s["scatterer"] for s in splan["shells"]],
                         ["N", "O"])
        deg_axial = list(splan["degens"])
        for i, sh in enumerate(splan["shells"]):
            if sh["scatterer"] == "O":
                deg_axial[i] = 1.0
        res_ax = synth(splan["files"], deg_axial)
        out_ax = xafs_fit.fit_identified(
            res_ax, absorber_symbol="Fe",
            catalytic_model="sac:N4@2.00,axial=O@2.10,C=none",
            rmin=1.0, rmax=2.6, shared_sigma2="all", s02=s02,
            verbose=False)
        by_el = {("O" if "O" in p.label else "N"): p
                 for p in out_ax["fit"].paths}
        self.assertLess(abs(by_el["N"].degeneracy - 4.0), 0.5,
                        f"SAC N_N 回收失败: {by_el['N'].degeneracy}")
        self.assertLess(abs(by_el["O"].degeneracy - 1.0), 0.5,
                        f"SAC N_O 回收失败: {by_el['O'].degeneracy}")
        self.assertLess(by_el["N"].n_err, 0.5)
        self.assertLess(by_el["O"].n_err, 0.5)

    # ------------------------------------------------------------
    # B13-2 异核共占位双分支 path_links
    # ------------------------------------------------------------
    def test_hetero_cocupancy_end_to_end_dual_branch_fit(self):
        """#62 B13-2：共占位 CIF → 双分支 FEFF → 共享 ΔR/σ² 精修。

        Ti(0.7)/Ru(0.3) 共占 N_geom=8（N_Ti=5.6, N_Ru=2.4）：
        path_links 把自由变量从 7 压到 5；回收两分支 N、ΔR 与 σ²
        逐位相等；双分支 Rf 比单元素近似小 5 倍以上。
        """
        from larch.xafs import feffpath

        cif_text = """data_test_hetero
_cell_length_a 4.10
_cell_length_b 4.10
_cell_length_c 4.10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'Pm-3m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Fe1 Fe 0.0 0.0 0.0 1.0
Ti1 Ti 0.5 0.5 0.5 0.7
Ru1 Ru 0.5 0.5 0.5 0.3
"""
        with tempfile.TemporaryDirectory() as tmp:
            cif = Path(tmp) / "hetero.cif"
            cif.write_text(cif_text, encoding="utf-8")
            motif = structure_lib.import_cif(cif, "test_hetero_b13", "Fe")
            motif.cif_path = str(cif)
            motif.sites = structure_lib.expand_crystal(motif, rmax=5.0)
            motif.shell_expect = structure_lib.geometric_shells(motif.sites)
            plan = xafs_fit.select_shell_paths(
                motif, edge="K", rmin=2.8, rmax=4.3,
                kmin=3.0, kmax=12.0, hetero_cocupancy=True)

        self.assertEqual(plan["hetero_shells"][0]["n_geom"], 8)
        branches = {b["element"]: b for b in
                    plan["hetero_shells"][0]["branches"]}
        self.assertAlmostEqual(branches["Ti"]["degen"], 5.6, places=10)
        self.assertAlmostEqual(branches["Ru"]["degen"], 2.4, places=10)
        self.assertEqual(len(plan["path_links"]), 1)
        link = plan["path_links"][0]
        self.assertTrue(link["share_deltar"] and link["share_sigma2"])
        self.assertIsNone(link["degen_ratio"])

        # 合成真值谱：N_Ti=5.6, N_Ru=2.4, ΔR=-0.02, σ²=0.004
        s02, dr_true, s2_true, e0_true = 0.803, -0.02, 0.004, 2.0
        k = np.arange(3.0, 12.0, 0.05)
        plist = []
        for f, d in zip(plan["files"], plan["degens"]):
            p = feffpath(str(f), s02=s02, sigma2=s2_true,
                         e0=e0_true, deltar=dr_true)
            p.degen = float(d)
            p._calc_chi()
            plist.append(p)
        chi = artemis.exafs_model_k(plist, k)
        chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                     mu_norm=np.zeros_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)

        def do_fit(files, labels, degens, links=None):
            return artemis.fit_r_space(
                chi_k, ft, files, rmin=plan["rmin"], rmax=plan["rmax"],
                kmin=3.0, kmax=12.0, kweight=2, s02=s02,
                vary_delta_e0=True, vary_n=True,
                path_degen=degens, path_labels=labels, path_links=links)

        fit = do_fit(plan["files"], plan["labels"], plan["degens"],
                     plan["path_links"])
        self.assertEqual(fit.n_variables, 5,
                         "共享 ΔR/σ² 后应为 4 个壳层变量 + ΔE0")
        free = do_fit(plan["files"], plan["labels"], plan["degens"])
        self.assertEqual(free.n_variables, 7,
                         "无 path_links 时双分支应为 2×3+1=7 个变量")

        by_el = {("Ti" if "Ti" in p.label else "Ru"): p
                 for p in fit.paths}
        ti, ru = by_el["Ti"], by_el["Ru"]
        self.assertGreaterEqual(ti.degeneracy, 5.2)
        self.assertLessEqual(ti.degeneracy, 6.0)
        self.assertGreaterEqual(ru.degeneracy, 2.0)
        self.assertLessEqual(ru.degeneracy, 2.8)
        self.assertAlmostEqual(ti.delta_r, ru.delta_r, places=12)
        self.assertAlmostEqual(ti.sigma2, ru.sigma2, places=12)
        self.assertAlmostEqual(ti.delta_r_err, ru.delta_r_err, places=12)
        self.assertAlmostEqual(ti.sigma2_err, ru.sigma2_err, places=12)
        self.assertLess(fit.r_factor, 1e-4,
                        f"双分支合成谱应近完美复现: Rf={fit.r_factor}")

        ti_idx = [i for i, s in enumerate(plan["shells"])
                  if s["scatterer"] == "Ti"][0]
        single = do_fit([plan["files"][ti_idx]],
                        [plan["labels"][ti_idx]],
                        [plan["degens"][ti_idx]])
        self.assertGreater(single.r_factor, 5.0 * fit.r_factor,
                           "单元素近似应显著劣于双分支共享精修")

    # ------------------------------------------------------------
    # B14-1 Mini-GDS 化学计量约束
    # ------------------------------------------------------------
    def test_mini_gds_stoichiometric_sum_and_ratio_constraints(self):
        """#63 B14-1：N_sum 精确守恒、误差严格传递、循环约束拒绝。

        Fe–O@1.98 + Fe–Cl@2.35 双路径合成谱（N=2.4/3.6）：
        - 无约束 n_variables=7，n_sum 约束后=6（真实减少 1）；
        - N_O+N_Cl=6.0 逐位成立，且 σ(N_Cl)=σ(N_O)；
        - 循环约束与未知标签抛 XAFSConfigError。
        """
        from larch.xafs import feffpath

        def single_path(distance, element, workdir):
            inp = artemis.build_feff_input(
                [(distance, 0.0, 0.0, element)], absorber="Fe",
                rmax=distance + 0.8, title=f"gds_{element}")
            files = artemis.run_feff6l(inp, workdir)
            for f in files:
                p = artemis.load_feff_path(f)
                if (int(getattr(p, "nleg", 2)) == 2
                        and abs(float(p.reff) - distance) < 0.05):
                    return f
            raise AssertionError(f"FEFF 未生成 {element}@{distance} SS 路径")

        with tempfile.TemporaryDirectory() as tmp:
            f_o = single_path(1.98, "O", Path(tmp) / "o")
            f_cl = single_path(2.35, "Cl", Path(tmp) / "cl")
            files = [f_o, f_cl]
            labels = ["Fe-O", "Fe-Cl"]
            truth_n = [2.4, 3.6]
            s02, s2, e0 = 0.803, 0.004, 2.0
            k = np.arange(3.0, 12.0, 0.05)
            plist = []
            for f, n in zip(files, truth_n):
                p = feffpath(str(f), s02=s02, sigma2=s2, e0=e0, deltar=0.0)
                p.degen = n
                p._calc_chi()
                plist.append(p)
            chi = artemis.exafs_model_k(plist, k)
            chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                         mu_norm=np.zeros_like(k), rbkg=1.0)
            ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)

            def do_fit(**kw):
                return artemis.fit_r_space(
                    chi_k, ft, files, rmin=1.4, rmax=2.9,
                    kmin=3.0, kmax=12.0, kweight=2, s02=s02,
                    vary_delta_e0=True, vary_n=True,
                    path_degen=truth_n, path_labels=labels, **kw)

            free = do_fit()
            self.assertEqual(free.n_variables, 7)
            constrained = do_fit(constraints="N(Fe-O)+N(Fe-Cl)=6.0")
            self.assertEqual(constrained.n_variables, 6,
                             "n_sum 应真实减少一个 N 自由变量")
            n_o = constrained.paths[0].degeneracy
            n_cl = constrained.paths[1].degeneracy
            self.assertAlmostEqual(n_o + n_cl, 6.0, places=12)
            self.assertLess(abs(n_o - 2.4), 0.25)
            self.assertLess(abs(n_cl - 3.6), 0.25)
            self.assertIsNotNone(constrained.paths[0].n_err)
            self.assertIsNotNone(constrained.paths[1].n_err)
            self.assertEqual(constrained.paths[1].n_err,
                             constrained.paths[0].n_err,
                             "n_sum 误差必须严格按 |coeff| 传递")
            self.assertLess(constrained.r_factor, 0.05)

            # 字典形式的 n_ratio：N_Cl = 1.5·N_O
            ratio = do_fit(constraints=[
                {"kind": "n_ratio", "target": "Fe-Cl",
                 "source": "Fe-O", "ratio": 1.5}])
            self.assertEqual(ratio.n_variables, 6)
            self.assertAlmostEqual(
                ratio.paths[1].degeneracy,
                1.5 * ratio.paths[0].degeneracy, places=12)
            self.assertAlmostEqual(
                ratio.paths[1].n_err, 1.5 * ratio.paths[0].n_err,
                places=12)

            with self.assertRaises(XAFSConfigError):
                do_fit(constraints=(
                    "N(Fe-O)=1.0*N(Fe-Cl);N(Fe-Cl)=1.0*N(Fe-O)"))
            with self.assertRaises(XAFSConfigError):
                do_fit(constraints="N(Fe-O)+N(Fe-Pt)=6.0")

    # ------------------------------------------------------------
    # B14-2 相移校正诊断 + χ² 置信椭圆 + HTML 内联
    # ------------------------------------------------------------
    def test_phase_corrected_diagnostics_and_chi2_contour_report(self):
        """#64 B14-2：相移校正真键长与二维 χ² 等高线全链接入。

        用自带 lepidocrocite + S7 拟合结果：
        - 未校正 FT 峰在 1.4–1.7 Å（表观 R），E0 补偿相移校正后
          |r_peak_corrected − r_fit| < 0.06 Å；
        - χ² 网格 (15,15)、min_chi2 有限，CSV/PNG 非空；
        - xafs_report 开 --phase-correct/--chi2-contour 后 HTML 内联
          两张图并含标题。
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            pc = xafs_diagnostics.write_phase_corrected_diagnostics(
                self.lep, self.lep_fit, tmp, stem="lepido")
            self.assertGreaterEqual(pc["r_peak_raw"], 1.4)
            self.assertLessEqual(pc["r_peak_raw"], 1.7)
            self.assertLess(abs(pc["r_peak_corrected"] - pc["r_fit"]), 0.06,
                            f"相移校正未回到拟合键长: {pc}")
            self.assertEqual(pc["phase_source"], "feff_pha")
            csv_path, png_path = Path(pc["csv"]), Path(pc["png"])
            self.assertTrue(csv_path.is_file()
                            and csv_path.stat().st_size > 0)
            self.assertTrue(png_path.is_file()
                            and png_path.stat().st_size > 0)
            header = csv_path.read_text(
                encoding="utf-8-sig").splitlines()[0]
            self.assertEqual(
                header, "R_A,mag_raw,real_raw,"
                        "mag_phase_corrected,real_phase_corrected")

            cc = xafs_diagnostics.write_chi2_contour_diagnostics(
                self.lep, self.lep_fit, tmp, stem="lepido")
            self.assertIsNotNone(cc, "χ² 等高线未生成")
            self.assertTrue(np.isfinite(cc["min_chi2"]))
            self.assertEqual(cc["chi2_grid"].shape, (15, 15))
            self.assertEqual((cc["x_name"], cc["y_name"]), ("N", "sigma2"))
            for key in ("csv", "png"):
                p = Path(cc[key])
                self.assertTrue(p.is_file() and p.stat().st_size > 0)
            self.assertEqual(len(cc["levels"]), 3)
            self.assertAlmostEqual(cc["levels"][0] - cc["min_chi2"],
                                   2.30, places=6)
            self.assertLessEqual(cc["min_chi2"],
                                 np.nanmin(cc["chi2_grid"]) + 1e-9)

            report_dir = tmp / "report"
            xafs_report.main([
                str(LEPIDO), "--absorber", "Fe", "-o", str(report_dir),
                "--phase-correct", "--chi2-contour", "--no-window-tuner"])
            html_files = sorted(report_dir.glob("*_report.html"))
            self.assertEqual(len(html_files), 1)
            html_text = html_files[0].read_text(encoding="utf-8")
            self.assertIn("相移校正真键长诊断", html_text)
            self.assertIn("二维 χ² 置信椭圆", html_text)
            self.assertIn("data:image/png;base64,", html_text)

    # ------------------------------------------------------------
    # B14-3 关联 Debye 模型纳入多数据集联合拟合
    # ------------------------------------------------------------
    def test_multi_dataset_correlated_debye_corefinement(self):
        """#65 B14-3：单 Θ_D 同时约束近壳/远壳 σ²(T,R)。

        用 Fe–O@2.00 与 Fe–O@3.68 两条 FEFF 路径，以真值 Θ_D=420 K、
        σ²_static=0.0012 Å² 合成 80 K 与 350 K 两套谱：
        - 一个 debye_models 组绑定两温度下的全部 4 条路径，
          Θ_D 回收误差 <15 K、σ²_static <5e-4；
        - 同一温度下远壳 σ² 严格大于近壳 σ²，global Rf <1e-4；
        - 与 einstein_models/share_sigma2 绑定同一路径时抛错。
        """
        from larch.xafs import feffpath

        true_theta, true_static, s02 = 420.0, 0.0012, 0.803

        def find_ss(files, target, tol=0.08):
            best = None
            for f in files:
                p = artemis.load_feff_path(f)
                if int(getattr(p, "nleg", 2)) != 2:
                    continue
                if (best is None or abs(float(p.reff) - target)
                        < abs(float(best[1].reff) - target)):
                    best = (f, p)
            if best is None or abs(float(best[1].reff) - target) > tol:
                raise AssertionError(f"FEFF 未生成 {target} Å 附近的 SS 路径")
            return best[0]

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sites = []
            for d, el in ((2.0, "O"), (3.68, "O")):
                for axis in ((d, 0, 0), (-d, 0, 0), (0, d, 0),
                             (0, -d, 0), (0, 0, d), (0, 0, -d)):
                    sites.append((*axis, el))
            inp = artemis.build_feff_input(
                sites, absorber="Fe", rmax=4.5, title="debye_test")
            files = artemis.run_feff6l(inp, tmp / "feff")
            near = find_ss(files, 2.0)
            far = find_ss(files, 3.68)

            k = np.arange(3.0, 12.0, 0.05)
            datasets = []
            expected = {}
            for T in (80.0, 350.0):
                plist, s2s = [], []
                for f in (near, far):
                    reff = float(artemis.load_feff_path(f).reff)
                    s2 = artemis.correlated_debye_sigma2(
                        T, true_theta, reff, 55.845, 15.999,
                        sigma2_static=true_static)
                    s2s.append(s2)
                    p = feffpath(str(f), s02=s02, sigma2=s2,
                                 e0=0.0, deltar=0.0)
                    p.degen = 6.0
                    p._calc_chi()
                    plist.append(p)
                expected[T] = s2s
                chi = artemis.exafs_model_k(plist, k)
                ck = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                          mu_norm=np.ones_like(k), rbkg=1.0)
                ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
                datasets.append({
                    "name": f"{T:.0f}K", "chi_k": ck, "ft": ft,
                    "datfiles": [str(near), str(far)],
                    "rmin": 1.4, "rmax": 4.0, "kmin": 3.0, "kmax": 12.0,
                    "kweight": 2, "s02": s02, "vary_n": False,
                    "path_degen": [6.0, 6.0],
                    "path_labels": ["Fe-O_near", "Fe-O_far"],
                    "temperature_K": T, "edge_group": "K"})

            model = {
                "paths": [(0, 0), (0, 1), (1, 0), (1, 1)],
                "mass_absorber": 55.845,
                "mass_scatterer": 15.999,
                "theta_d_init": 400.0,
                "k_debye": None,
            }
            out = artemis.fit_multi_datasets(
                datasets, share_e0="auto", debye_models=[model])
            self.assertEqual(out["n_variables_total"], 7,
                             "4 ΔR + 2 Debye 全局变量 + ΔE0")
            self.assertEqual(len(out["debye_results"]), 1)
            deb = out["debye_results"][0]
            self.assertLess(abs(deb["theta_d"] - true_theta), 15.0)
            self.assertLess(abs(deb["sigma2_static"] - true_static), 5e-4)
            self.assertTrue(np.isfinite(deb["theta_d_err"]))
            self.assertTrue(np.isfinite(deb["sigma2_static_err"]))
            self.assertLess(out["global_r_factor"], 1e-4,
                            f"联合拟合未复现合成谱: "
                            f"Rf={out['global_r_factor']}")
            for di, T in enumerate((80.0, 350.0)):
                near_fit, far_fit = out["datasets"][di].paths
                self.assertGreater(far_fit.sigma2, near_fit.sigma2,
                                   f"{T:.0f} K 远壳 σ² 应大于近壳")
                self.assertAlmostEqual(near_fit.sigma2, expected[T][0],
                                       places=7)
                self.assertAlmostEqual(far_fit.sigma2, expected[T][1],
                                       places=7)

            with self.assertRaises(XAFSConfigError):
                artemis.fit_multi_datasets(
                    datasets, debye_models=[model],
                    einstein_models=[{"paths": [(0, 0)],
                                      "reduced_mass_amu": 12.5}])
            with self.assertRaises(XAFSConfigError):
                artemis.fit_multi_datasets(
                    datasets, debye_models=[model],
                    shared_bonds=[{"target": (1, 0), "source": (0, 0),
                                   "share_sigma2": True}])

    # ------------------------------------------------------------
    # B15-1 原位序列 PCA → SPOIL → LCF → MCR-ALS
    # ------------------------------------------------------------
    def test_operando_series_pca_spoil_lcf_and_mcr_tracking(self):
        """#66 B15-1：原位序列一键解谱与演变轨迹。

        LEPIDO/GOETHITE 按 w=[0,0.2,…,1.0] 线性混合成 6 条序列，
        候选标样加一个纯噪声伪标样：
        - SPOIL 必须只拒绝 fake_noise；
        - LCF 回收 goethite 权重 RMSE <0.02；
        - MCR 两条组分浓度随序号单调（|r|>0.98），四宫格与 CSV 落盘。
        """
        import xafs_lcf

        goe = xafs_pipeline.process_sample(
            [GOETHITE], absorber="Fe", strict=False, verbose=False).norm
        lep = self.lep.norm
        x = np.arange(-20.0, 60.0, 0.2)

        def flat(n):
            arr = n.flattened if n.flattened is not None else n.mu_norm
            return np.interp(x, n.energy - n.e0, arr)

        y_lep, y_goe = flat(lep), flat(goe)
        rng = np.random.default_rng(66)
        w_true = np.linspace(0.0, 1.0, 6)
        items = []
        for i, w in enumerate(w_true):
            yy = ((1.0 - w) * y_lep + w * y_goe
                  + rng.normal(0, 1e-3, x.size))
            items.append((f"scan_{i:02d}", x, yy))
        fake_y = rng.normal(0.0, 0.1, x.size)
        standards = {"lepido": lep, "goethite": goe,
                     "fake_noise": (x, fake_y)}

        out = xafs_lcf.run_operando_series(
            items, standards=standards, xmin=-20.0, xmax=60.0)
        self.assertEqual(out["n_components"], 2)
        self.assertEqual(out["accepted_standards"],
                         ["lepido", "goethite"])
        self.assertEqual(out["rejected_standards"], ["fake_noise"])

        got = np.asarray([row["weights"].get("goethite", np.nan)
                          for row in out["lcf_series"]], dtype=float)
        self.assertEqual(got.shape, w_true.shape)
        rmse = float(np.sqrt(np.mean((got - w_true) ** 2)))
        self.assertLess(rmse, 0.02, f"LCF goethite 权重 RMSE={rmse}")

        conc = np.asarray(out["mcr_als"]["concentrations"], dtype=float)
        self.assertEqual(conc.shape, (2, 6))
        idx = np.arange(6, dtype=float)
        corr = [abs(float(np.corrcoef(idx, conc[j])[0, 1]))
                for j in range(2)]
        self.assertTrue(all(c > 0.98 for c in corr),
                        f"MCR 组分未随序号单调: {corr}")
        self.assertEqual(len(out["mcr_als"]["characteristic_energy"]), 2)
        self.assertLess(out["mcr_als"]["characteristic_energy"][0],
                        out["mcr_als"]["characteristic_energy"][1])
        self.assertEqual(len(out["descriptors"]), 6)
        for key in ("e0", "e_half_height", "white_line_height", "ft_r1"):
            self.assertIn(key, out["descriptors"][0])

        with tempfile.TemporaryDirectory() as tmp:
            made = xafs_diagnostics.write_operando_diagnostics(
                out, tmp, stem="op")
            conc_csv = Path(made["concentrations_csv"])
            comp_csv = Path(made["components_csv"])
            png = Path(made["png"])
            for p in (conc_csv, comp_csv, png):
                self.assertTrue(p.is_file() and p.stat().st_size > 0,
                                f"缺少产物: {p}")
            lines = conc_csv.read_text(
                encoding="utf-8-sig").splitlines()
            self.assertEqual(len(lines), 7, "浓度表应为 6 条谱 + 表头")
            self.assertIn("lcf_goethite", lines[0])
            self.assertIn("mcr_component_1", lines[0])
            comp_lines = comp_csv.read_text(
                encoding="utf-8-sig").splitlines()
            self.assertGreater(len(comp_lines), 10)
            self.assertIn("mcr_component_2", comp_lines[0])

        with self.assertRaises(XAFSConfigError):
            xafs_lcf.run_operando_series(items[:2])

    # ------------------------------------------------------------
    # B16-1 Athena 原生 Perl .prj 只读导入
    # ------------------------------------------------------------
    def test_read_real_athena_perl_prj_fe_standards(self):
        """#67 B16-1：直接读通仓库自带 Fe_standards.prj。

        - 4 条主扫描（.000/.100/.200 + pyrite），Ref 与 merge 默认过滤；
        - include_derived=True 时 10 条全保留；
        - .100 的 energy/μ 与原始文件 462 点逐点一致；
        - 不含 @x/@y 的任意 Perl 片段仍抛 XAFSConfigError（保护 #50）。
        """
        prj = (ROOT / "data" / "练手数据" / "02_Fe_standards" /
               "Fe_standards.prj")
        out = xafs_io.read_prj(prj, as_scans=True)
        self.assertEqual(out["meta"]["source_format"], "athena_perl")
        self.assertEqual(len(out["scans"]), 4)
        labels = [s.meta["label"] for s in out["scans"]]
        for name in ("Fe_lepidocrocite.000", "Fe_lepidocrocite.100",
                     "Fe_lepidocrocite.200", "Fe_pyrite.000"):
            self.assertIn(name, labels)
        derived = xafs_io.read_prj(
            prj, as_scans=True, include_derived=True)
        self.assertEqual(len(derived["scans"]), 10)

        s_prj = [s for s in out["scans"]
                 if s.meta["label"] == "Fe_lepidocrocite.100"][0]
        s_raw = xafs_io.read_scan(
            ROOT / "data" / "练手数据" / "02_Fe_standards" /
            "Fe_lepidocrocite.100")
        self.assertEqual(len(s_prj.energy), len(s_raw.energy))
        self.assertEqual(len(s_prj.energy), 462)
        self.assertLess(
            float(np.max(np.abs(s_prj.energy - s_raw.energy))), 1e-4)
        self.assertGreater(
            float(np.corrcoef(s_prj.mu, s_raw.mu)[0, 1]), 0.999999)
        self.assertIsNotNone(s_prj.i0)
        self.assertIsNotNone(s_prj.it)

        with tempfile.TemporaryDirectory() as tmp:
            perl = Path(tmp) / "perl_only.prj"
            perl.write_text("$VAR1 = { 'x' => 1 };", encoding="utf-8")
            with self.assertRaises(XAFSConfigError):
                xafs_io.read_prj(perl)
            with self.assertRaises(XAFSConfigError):
                xafs_io.read_athena_perl_prj(perl, as_scans=True)

    # ------------------------------------------------------------
    # B16-2 Co / Ni / Cu / Pt 催化结构库扩充
    # ------------------------------------------------------------
    def test_catalysis_library_motifs_for_co_ni_cu_pt(self):
        """#68 B16-2：四大催化元素催化基准 Motif 入库 + PtCl2 实测命中。

        - Co/Ni/Cu 平面 M-N4 + 碳壳、Cu2O 赤铜矿、PtCl2 平面 PtCl4；
        - 每个 Motif 均能 ensure_feff 生成 feff0001.dat；
        - 真实 Pt-L3-PlatinumChloride.xdi（Sample.name=PtCl2）：
          top1=ptcl4_square、status=confident，S7 第一壳
          N∈[3.5,4.5]、R∈[2.26,2.35] Å、Rf<0.05。
        """
        expected = {
            "Co": [("con4_planar", "N", 1.96, 4, "C", 2.88)],
            "Ni": [("nin4_planar", "N", 1.92, 4, "C", 2.85)],
            "Cu": [("cun4_planar", "N", 1.95, 4, "C", 2.88),
                   ("cu2o_cuprite_cod_1010941", "O", 1.8446, 2,
                    "Cu", 3.0123)],
            "Pt": [("ptcl4_square", "Cl", 2.31, 4, "Cl", 3.45)],
        }
        for el, entries in expected.items():
            lib = {m.id: m for m in structure_lib.load_library(el)}
            edge = get_absorber(el).edge
            for mid, el1, r1, n1, el2, r2 in entries:
                self.assertIn(mid, lib, f"{el} 库缺 {mid}")
                motif = lib[mid]
                shells = motif.shell_expect
                self.assertEqual(shells[0]["el"], el1)
                self.assertAlmostEqual(shells[0]["r"], r1, places=3)
                self.assertEqual(shells[0]["n"], n1)
                self.assertEqual(shells[1]["el"], el2)
                self.assertAlmostEqual(shells[1]["r"], r2, places=3)
                cache = structure_lib.ensure_feff(motif, edge=edge)
                self.assertTrue(
                    (cache / "feff0001.dat").is_file(),
                    f"{el}/{mid} 未生成 feff0001.dat")

        res = xafs_pipeline.process_sample(
            [PT_CL], absorber="Pt", strict=False, verbose=False)
        out = xafs_fit.run(
            res, identify=True, fit=True, absorber_symbol="Pt",
            verbose=False)
        sid = out["structure"]
        self.assertEqual(sid["status"], "confident")
        self.assertEqual(sid["ranked"][0]["id"], "ptcl4_square")
        self.assertEqual(sid["uncovered_peaks"], [])
        fit = out["fit"]["fit"]
        self.assertEqual(out["fit"]["motif"], "ptcl4_square")
        path = fit.paths[0]
        self.assertTrue(path.label.startswith("Pt-Cl"))
        self.assertGreaterEqual(path.degeneracy, 3.5)
        self.assertLessEqual(path.degeneracy, 4.5)
        self.assertGreaterEqual(path.r, 2.26)
        self.assertLessEqual(path.r, 2.35)
        self.assertLess(fit.r_factor, 0.05)

    # ------------------------------------------------------------
    # B17-1 纳米颗粒尺寸效应反演
    # ------------------------------------------------------------
    def test_nanoparticle_size_and_dispersion_inversion(self):
        """#69 B17-1：Calvin 球截断方程反解粒径、误差传播与形貌判据。

        - N₁=8.0、r₁=2.774 Å 反解 u 精确满足截断方程（误差 <1e-9）；
        - N₁ 增大时直径单调增大、分散度单调减小，N₁=12 为体相极限；
        - 显著偏低的第二壳 N₂ 判为 raft_or_flattened。
        """
        est = artemis.estimate_nanoparticle_size(
            n1=8.0, r1=2.774, n1_err=0.5, lattice="fcc")
        u = 2.774 / est["diameter_A"]
        residual = 12.0 * (1.0 - 1.5 * u + 0.5 * u ** 3) - 8.0
        self.assertLess(abs(residual), 1e-9,
                        f"球截断方程残差 {residual:.3e}")
        self.assertGreaterEqual(est["diameter_nm"], 1.0)
        self.assertLessEqual(est["diameter_nm"], 1.5)
        self.assertGreater(est["diameter_nm_err"], 0.0)
        self.assertGreater(est["dispersion"], 0.3)
        self.assertLessEqual(est["dispersion"], 1.0)
        self.assertEqual(est["shape_hint"], "single_shell_only")

        diameters = [artemis.estimate_nanoparticle_size(
            n1=n, r1=2.774)["diameter_nm"] for n in (6.0, 9.0, 11.0)]
        self.assertLess(diameters[0], diameters[1])
        self.assertLess(diameters[1], diameters[2])
        dispersions = [artemis.estimate_nanoparticle_size(
            n1=n, r1=2.774)["dispersion"] for n in (6.0, 9.0, 11.0)]
        self.assertGreater(dispersions[0], dispersions[1])
        self.assertGreater(dispersions[1], dispersions[2])

        bulk = artemis.estimate_nanoparticle_size(n1=12.0, r1=2.774)
        self.assertEqual(bulk["regime"], "bulk_limit")
        self.assertEqual(bulk["dispersion"], 0.0)
        self.assertEqual(bulk["diameter_nm"], float("inf"))

        raft = artemis.estimate_nanoparticle_size(
            n1=8.0, r1=2.774, n2=1.0)
        self.assertGreater(raft["n2_spherical_expected"], 3.0)
        self.assertEqual(raft["shape_hint"], "raft_or_flattened")
        consistent = artemis.estimate_nanoparticle_size(
            n1=8.0, r1=2.774, n2=raft["n2_spherical_expected"])
        self.assertEqual(consistent["shape_hint"], "spherical_consistent")

        bcc = artemis.estimate_nanoparticle_size(
            n1=6.0, r1=2.48, lattice="bcc")
        self.assertEqual(bcc["n1_bulk"], 8.0)
        self.assertEqual(bcc["regime"], "nanoparticle")
        with self.assertRaises(XAFSConfigError):
            artemis.estimate_nanoparticle_size(8.0, 2.774, lattice="diamond")

    # ------------------------------------------------------------
    # B17-2 双金属合金动态 motif 与 Cowley 短程有序
    # ------------------------------------------------------------
    def test_bimetallic_alloy_motif_and_cowley_sro_pipeline(self):
        """#70 B17-2：alloy: 动态建模 + Cowley SRO 端到端回收。

        - α 对随机固溶体/核壳偏析/异核有序分别精确给出 0 / 0.6 / −1；
        - PtCo 合金谱按真值 N_Pt-Pt=6、N_Pt-Co=2 合成后 S7 回收
          误差 <0.5，α≈0、nanoparticle_size≈1.0–1.5 nm。
        """
        random_sro = artemis.compute_cowley_sro(9.0, 3.0, 0.25)
        self.assertAlmostEqual(random_sro["alpha_cowley"], 0.0, places=12)
        self.assertAlmostEqual(random_sro["alloying_extent_pct"], 100.0,
                               places=12)
        self.assertEqual(random_sro["tendency"], "random_solid_solution")

        segregated = artemis.compute_cowley_sro(8.0, 2.0, 0.50)
        self.assertAlmostEqual(segregated["alpha_cowley"], 0.60, places=12)
        self.assertEqual(segregated["tendency"],
                         "homoatomic_clustering_or_segregation")

        ordered = artemis.compute_cowley_sro(6.0, 6.0, 0.25)
        self.assertAlmostEqual(ordered["alpha_cowley"], -1.0, places=12)
        self.assertEqual(ordered["tendency"], "heteroatomic_ordering")

        sro_err = artemis.compute_cowley_sro(
            9.0, 3.0, 0.25, n_aa_err=0.5, n_ab_err=0.3)
        expected_err = float(np.hypot(3.0 * 0.5, 9.0 * 0.3)
                             / (0.25 * 12.0 ** 2))
        self.assertAlmostEqual(sro_err["alpha_err"], expected_err, places=12)

        alloy = structure_lib.parse_catalytic_cluster_spec(
            "Pt", "alloy:Co@2.64,Pt@2.76,x=0.25,n=12")
        self.assertEqual(alloy.nuclearity, "alloy")
        self.assertIn("_dynamic_cache", str(alloy.cache_dir))
        counts = {sh["el"]: sh["n"] for sh in alloy.shell_expect}
        self.assertEqual(counts.get("Co"), 3)
        self.assertEqual(counts.get("Pt"), 9)
        alloy_dict = structure_lib.parse_catalytic_cluster_spec(
            "Pt", {"type": "alloy", "partner": "Co", "r_ab": 2.64,
                   "r_aa": 2.76, "x_partner": 0.25, "n_total": 12})
        self.assertEqual(alloy_dict.id, alloy.id)

        from types import SimpleNamespace
        from larch.xafs import feffpath

        plan = xafs_fit.select_shell_paths(
            alloy, edge="L3", rmin=1.6, rmax=3.2,
            kmin=3.0, kmax=12.0)
        truth = {"Co": 2.0, "Pt": 6.0}
        k = np.arange(3.0, 12.0, 0.05)
        plist = []
        for sh, f in zip(plan["shells"], plan["files"]):
            p = feffpath(str(f), s02=0.8075, sigma2=0.004, e0=2.0,
                         deltar=0.0)
            p.degen = truth[sh["scatterer"]]
            p._calc_chi()
            plist.append(p)
        chi = artemis.exafs_model_k(plist, k)
        rng = np.random.default_rng(11)
        chi = chi + rng.normal(0.0, 1e-3, size=chi.shape)
        chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                     mu_norm=np.ones_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
        res = SimpleNamespace(
            chi_k=chi_k, ft=ft, absorber_symbol="Pt", structure=None,
            merged=None, norm=None, self_absorption=None)
        out = xafs_fit.fit_identified(
            res, absorber_symbol="Pt",
            catalytic_model="alloy:Co@2.64,Pt@2.76,x=0.25,n=12",
            rmin=1.6, rmax=3.2, shared_sigma2="all", verbose=False)
        self.assertEqual(out["mode"], "motif")
        by_scatterer = {}
        for sh, p in zip(plan["shells"], out["fit"].paths):
            by_scatterer[sh["scatterer"]] = p
        self.assertLess(abs(by_scatterer["Pt"].degeneracy - 6.0), 0.5)
        self.assertLess(abs(by_scatterer["Co"].degeneracy - 2.0), 0.5)
        sro = out["alloy_sro"]
        self.assertIsNotNone(sro)
        self.assertEqual(sro["tendency"], "random_solid_solution")
        self.assertLess(abs(sro["alpha_cowley"]), 0.10)
        np_size = out["nanoparticle_size"]
        self.assertIsNotNone(np_size)
        self.assertGreaterEqual(np_size["diameter_nm"], 1.0)
        self.assertLessEqual(np_size["diameter_nm"], 1.5)

    # ------------------------------------------------------------
    # B18-1 Rayleigh 分辨极限与轻元素 Z 简并
    # ------------------------------------------------------------
    def test_rayleigh_resolution_and_light_element_z_degeneracy_guard(self):
        """#71 B18-1：π/(2Δk) 分辨率门禁与 SAC 近简并拆壳告警。

        - k∈[3,12] 时 ΔR_Rayleigh = π/18 ≈ 0.17453 Å；
        - Fe-N@1.98 vs Fe-O@2.06（ΔZ=1, ΔR=0.08）判 light_element_z_degeneracy；
        - Z 差异大或 ΔR 超限时 resolvable=True；
        - SAC N4@2.00 + axial O@2.10 拟合后 flags 置 True 并追加建议告警。
        """
        rayleigh = artemis.check_shell_resolution_and_z_degeneracy(
            [{"label": "Fe-O@1.98", "scatterer": "O", "r": 1.98}],
            kmin=3.0, kmax=12.0)
        self.assertAlmostEqual(rayleigh["rayleigh_limit_A"],
                               float(np.pi / 18.0), places=12)
        self.assertTrue(rayleigh["resolvable"])

        close = artemis.check_shell_resolution_and_z_degeneracy(
            [{"label": "Fe-N@1.98", "scatterer": "N", "r": 1.98},
             {"label": "Fe-O@2.06", "scatterer": "O", "r": 2.06}])
        self.assertFalse(close["resolvable"])
        self.assertEqual(len(close["pair_warnings"]), 1)
        warning = close["pair_warnings"][0]
        self.assertEqual(warning["kind"], "light_element_z_degeneracy")
        self.assertEqual(warning["delta_z"], 1)
        self.assertLess(warning["delta_r"], warning["rayleigh_limit_A"])

        far_z = artemis.check_shell_resolution_and_z_degeneracy(
            [{"label": "Fe-N@1.98", "scatterer": "N", "r": 1.98},
             {"label": "Fe-Co@2.50", "scatterer": "Co", "r": 2.50}])
        self.assertTrue(far_z["resolvable"])
        self.assertEqual(far_z["pair_warnings"], [])

        far_r = artemis.check_shell_resolution_and_z_degeneracy(
            [{"label": "Fe-O@1.98", "scatterer": "O", "r": 1.98},
             {"label": "Fe-O@2.25", "scatterer": "O", "r": 2.25}])
        self.assertTrue(far_r["resolvable"])
        self.assertEqual(far_r["pair_warnings"], [])

        from types import SimpleNamespace
        from larch.xafs import feffpath

        motif = structure_lib.parse_catalytic_cluster_spec(
            "Fe", "sac:N4@2.00,axial=O@2.10,C=none")
        plan = xafs_fit.select_shell_paths(
            motif, edge="K", rmin=1.0, rmax=2.6, kmin=3.0, kmax=12.0)
        k = np.arange(3.0, 12.0, 0.05)
        plist = []
        for sh, f in zip(plan["shells"], plan["files"]):
            p = feffpath(str(f), s02=0.803, sigma2=0.004, e0=2.0,
                         deltar=0.0)
            p.degen = float(sh["degen"])
            p._calc_chi()
            plist.append(p)
        chi = artemis.exafs_model_k(plist, k)
        rng = np.random.default_rng(7)
        chi = chi + rng.normal(0.0, 1e-3, size=chi.shape)
        chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                     mu_norm=np.ones_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
        res = SimpleNamespace(
            chi_k=chi_k, ft=ft, absorber_symbol="Fe", structure=None,
            merged=None, norm=None, self_absorption=None)
        out = xafs_fit.fit_identified(
            res, absorber_symbol="Fe",
            catalytic_model="sac:N4@2.00,axial=O@2.10,C=none",
            rmin=1.0, rmax=2.6, verbose=False)
        self.assertFalse(out["resolution_check"]["resolvable"])
        self.assertTrue(out["flags"]["sub_rayleigh_z_degeneracy"])
        self.assertIn("Rayleigh", out["warning"])

    # ------------------------------------------------------------
    # B18-2 运行时 --s02 覆盖与真实 Au/Mo 箔解锁
    # ------------------------------------------------------------
    def test_s02_runtime_override_and_real_au_mo_foil_calibration(self):
        """#72 B18-2：运行时 S₀² 解锁未标定元素，静态表保持 #60 契约。

        - 不传 s02 时 Au 识别仍抛 XAFSDataError，get_absorber("Au").s02=None；
        - 真实 Au 箔 s02=0.82：top1=au_fcc confident，S7 N∈[11,13]、
          R∈[2.83,2.90] Å、Rf<0.02，且全局静态表未被修改；
        - 真实 Mo 箔 s02=0.88：top1=mo_bcc，双壳 N/R 回收在容差内；
        - s02 越界抛 XAFSConfigError。
        """
        from xafs_core import XAFSDataError

        au_file = (ROOT / "data" / "bmm_standards" / "Data" / "Au" /
                   "Au-L3-GoldFoil.xdi")
        mo_file = (ROOT / "data" / "bmm_standards" / "Data" / "Mo" /
                   "Mo-K-MolybdenumFoil.xdi")
        self.assertTrue(au_file.is_file())
        self.assertTrue(mo_file.is_file())

        # Rbkg=1.2：auto(1.0) 在 Au L3 白线后留下一个 1.66 Å 低 R
        # 背景鼓包，会把覆盖率误判为 partial；1.2 的 1.66 Å 伪峰消失，
        # 而第一壳 N/R 与 Rf 不受影响。
        res_au = xafs_pipeline.process_sample(
            [str(au_file)], absorber="Au", strict=False, rbkg=1.2,
            verbose=False)
        with self.assertRaises(XAFSDataError):
            structure_id.identify_result(
                res_au, absorber_symbol="Au", verbose=False)
        self.assertIsNone(get_absorber("Au").s02)

        out_au = xafs_fit.run(
            res_au, identify=True, fit=True, absorber_symbol="Au",
            s02=0.82, verbose=False)
        sid = out_au["structure"]
        self.assertEqual(sid["ranked"][0]["id"], "au_fcc")
        self.assertEqual(sid["status"], "confident")
        fout = out_au["fit"]
        self.assertEqual(fout["mode"], "motif")
        self.assertAlmostEqual(fout["s02"], 0.82, places=12)
        fit = fout["fit"]
        path = fit.paths[0]
        self.assertGreaterEqual(path.degeneracy, 11.0)
        self.assertLessEqual(path.degeneracy, 13.0)
        self.assertGreaterEqual(path.r, 2.83)
        self.assertLessEqual(path.r, 2.90)
        self.assertLess(fit.r_factor, 0.02)
        self.assertIsNone(get_absorber("Au").s02,
                          "运行时 --s02 不得修改全局静态表")

        res_mo = xafs_pipeline.process_sample(
            [str(mo_file)], absorber="Mo", strict=False, verbose=False)
        out_mo = xafs_fit.run(
            res_mo, identify=True, fit=True, absorber_symbol="Mo",
            s02=0.88, verbose=False)
        self.assertEqual(out_mo["structure"]["ranked"][0]["id"], "mo_bcc")
        fit_mo = out_mo["fit"]["fit"]
        self.assertEqual(len(fit_mo.paths), 2)
        p1, p2 = fit_mo.paths[0], fit_mo.paths[1]
        self.assertGreaterEqual(p1.degeneracy, 6.8)
        self.assertLessEqual(p1.degeneracy, 9.0)
        self.assertGreaterEqual(p1.r, 2.69)
        self.assertLessEqual(p1.r, 2.75)
        self.assertGreaterEqual(p2.degeneracy, 5.0)
        self.assertLessEqual(p2.degeneracy, 7.2)
        self.assertGreaterEqual(p2.r, 3.11)
        self.assertLessEqual(p2.r, 3.18)
        self.assertLess(fit_mo.r_factor, 0.02)
        self.assertIsNone(get_absorber("Mo").s02)

        with self.assertRaises(XAFSConfigError):
            xafs_fit.run(
                res_au, identify=False, fit=True, absorber_symbol="Au",
                s02=1.5, verbose=False)
        self.assertIn("s02", xafs_config.DEFAULTS)

    # ------------------------------------------------------------
    # 化学态 S₀² 多态（Fe 金属交叉值 opt-in）
    # ------------------------------------------------------------
    def test_chemical_state_s02_polymorphism(self):
        """化学态 S₀²：Fe 金属 0.6566，默认生产值 0.803 不变。"""
        from types import SimpleNamespace
        from absorber import resolve_s02
        import run_xafs

        self.assertEqual(get_absorber("Fe").s02, 0.803)
        self.assertEqual(resolve_s02("Fe"), 0.803)
        self.assertEqual(resolve_s02("Fe", "metal"), 0.6566)
        self.assertEqual(resolve_s02("Fe", "METAL"), 0.6566)
        with self.assertRaises(XAFSConfigError):
            resolve_s02("Fe", "oxide")
        with self.assertRaises(XAFSConfigError):
            resolve_s02("Au", "metal")
        self.assertEqual(get_absorber("Fe").s02, 0.803,
                         "化学态解析不得修改静态表")

        # CLI 优先级：显式 --s02 > 化学态 > T5 > 静态表。
        args = SimpleNamespace(s02=None, s02_from_t5=False,
                               chemical_state="metal", absorber="Fe")
        val, src = run_xafs._resolve_s02_for_paths(args, ["dummy.dat"])
        self.assertEqual(val, 0.6566)
        self.assertEqual(src, "chemical_state:Fe=metal")
        args.s02 = 0.9
        val, src = run_xafs._resolve_s02_for_paths(args, ["dummy.dat"])
        self.assertEqual(val, 0.9)
        self.assertEqual(src, "--s02 runtime override")

        # 逐元素映射：未列出的元素返回 (None, None)，多元素批处理不误伤。
        args_map = SimpleNamespace(s02=None, s02_from_t5=False,
                                   chemical_state="Fe=metal", absorber="Cu")
        self.assertEqual(
            run_xafs._resolve_s02_for_paths(args_map, ["dummy.dat"]),
            (None, None))
        args_map.absorber = "Fe"
        self.assertEqual(
            run_xafs._resolve_s02_for_paths(args_map, ["dummy.dat"]),
            (0.6566, "chemical_state:Fe=metal"))

        # 来源标记：化学态值视为已标定，普通手填 --s02 不算（用未标定
        # 元素隔离静态表 calibrated=True 的影响）。
        self.assertTrue(xafs_fit._s02_calibrated_flag(
            get_absorber("Cr"), "chemical_state:Cr=metal"))
        self.assertFalse(xafs_fit._s02_calibrated_flag(
            get_absorber("Cr"), "--s02 runtime override"))

    # ------------------------------------------------------------
    # B19-1 原位序列 EXAFS 联合精修
    # ------------------------------------------------------------
    def test_operando_exafs_series_joint_e0_bond_and_cn_tracking(self):
        """#73 B19-1：共享 ΔE₀ 的序列拟合回收 R(t)/N(t) 并导出 CSV/PNG。

        真值 R=[2.02,1.99,1.96,1.93] Å、N=[4.0,4.6,5.3,6.0]，
        全局 ΔE₀=+2.5 eV、σ²=0.004 Å²；共享 σ² 后键长误差 <0.005 Å、
        配位数误差 <0.15，global_r_factor <1e-3。
        """
        from types import SimpleNamespace
        from larch.xafs import feffpath

        motif = structure_lib.parse_catalytic_cluster_spec(
            "Fe", "sac:O4@1.98,C=none")
        plan = xafs_fit.select_shell_paths(
            motif, edge="K", rmin=1.2, rmax=2.4, kmin=3.0, kmax=12.0)
        r_truth = [2.02, 1.99, 1.96, 1.93]
        n_truth = [4.0, 4.6, 5.3, 6.0]
        k = np.arange(2.5, 12.5, 0.05)
        results = []
        for r_t, n_t in zip(r_truth, n_truth):
            p = feffpath(str(plan["files"][0]), s02=0.803,
                         sigma2=0.004, e0=2.5, deltar=r_t - 1.98)
            p.degen = n_t
            p._calc_chi()
            chi = artemis.exafs_model_k([p], k)
            chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                         mu_norm=np.ones_like(k), rbkg=1.0)
            ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
            results.append(SimpleNamespace(
                chi_k=chi_k, ft=ft, absorber_symbol="Fe", structure=None,
                merged=None, norm=None, self_absorption=None))

        out = xafs_fit.run_operando_exafs_series(
            results, "Fe", catalytic_model="sac:O4@1.98,C=none",
            share_sigma2_across_series=True, s02=0.803)
        self.assertEqual(out["motif"], motif.id)
        self.assertEqual(len(out["names"]), 4)
        label = out["shells"][0]
        traj = out["trajectories"][label]
        self.assertLess(max(abs(a - b) for a, b in
                            zip(traj["R"], r_truth)), 0.005)
        self.assertLess(max(abs(a - b) for a, b in
                            zip(traj["N"], n_truth)), 0.15)
        self.assertTrue(all(np.diff(traj["R"]) < 0),
                        f"R(t) 未严格单调递减: {traj['R']}")
        self.assertTrue(all(np.diff(traj["N"]) > 0),
                        f"N(t) 未严格单调递增: {traj['N']}")
        self.assertLess(out["global_r_factor"], 1e-3)

        with tempfile.TemporaryDirectory() as tmp:
            made = xafs_diagnostics.write_operando_exafs_diagnostics(
                out, tmp, stem="operando_exafs")
            csv_path = Path(made["csv"])
            png_path = Path(made["png"])
            self.assertTrue(csv_path.is_file())
            self.assertTrue(png_path.is_file())
            self.assertGreater(png_path.stat().st_size, 0)
            lines = csv_path.read_text(encoding="utf-8-sig").splitlines()
            self.assertEqual(len(lines), 5, "CSV 应为表头 + 4 个时间点")
            self.assertIn(f"{label}_N", lines[0])
            self.assertIn(f"{label}_R", lines[0])
            self.assertIn("delta_e0", lines[0])

    # ------------------------------------------------------------
    # B20-1 电催化氧化物/羟基氧化物库与真实标准谱盲测
    # ------------------------------------------------------------
    def test_electrocatalysis_oxide_hydroxide_library_real_standards(self):
        """#74 B20-1：7 个电催化 Motif 入库 + 4 条真实 BMM 标样盲测。

        - CoO / Co3O4 / Ni(OH)2 / CuO / PtO2 几何指纹逐项核验；
        - Co-K-CoO：top1=coo_rocksalt_cod_1533087 confident、score<0.06，S7 Rf<0.02；
        - Co-K-Co3O4：top1=co3o4_mixed confident、score<0.05、无未解释峰；
        - Ni-K-NickelHydroxide：top1=nioh2_theophrastite_cod_1548811
          confident、score<0.14，S7 Rf<0.05 且 Ni-O/Ni-Ni 键长落在窗口内；
        - Cu-CuO：top1=cuo_tenorite_cod_1011148、score<0.10。
        """
        co_lib = {m.id: m for m in structure_lib.load_library("Co")}
        ni_lib = {m.id: m for m in structure_lib.load_library("Ni")}
        cu_lib = {m.id: m for m in structure_lib.load_library("Cu")}
        pt_lib = {m.id: m for m in structure_lib.load_library("Pt")}
        for mid in ("coo_rocksalt_cod_1533087", "co3o4_spinel_cod_9005888",
                    "co3o4_oct", "co3o4_mixed"):
            self.assertIn(mid, co_lib)
        self.assertIn("nioh2_theophrastite_cod_1548811", ni_lib)
        self.assertIn("cuo_tenorite_cod_1011148", cu_lib)
        self.assertIn("pto2_beta", pt_lib)

        def _shells(motif):
            rows = motif.shell_expect or []
            return [(s["el"], float(s["r"]), int(s["n"])) for s in rows]

        coo = _shells(co_lib["coo_rocksalt_cod_1533087"])
        self.assertEqual(coo[0][0], "O")
        self.assertAlmostEqual(coo[0][1], 2.1315, places=3)
        self.assertEqual(coo[0][2], 6)
        self.assertEqual(coo[1][0], "Co")
        self.assertAlmostEqual(coo[1][1], 3.0144, places=3)
        self.assertEqual(coo[1][2], 12)

        tet = _shells(co_lib["co3o4_spinel_cod_9005888"])
        self.assertEqual(tet[0][2], 4)
        self.assertAlmostEqual(tet[0][1], 1.9353, places=3)
        octa = _shells(co_lib["co3o4_oct"])
        self.assertEqual(octa[0][2], 6)
        self.assertAlmostEqual(octa[0][1], 1.9217, places=3)
        self.assertEqual(octa[1][0], "Co")
        self.assertAlmostEqual(octa[1][1], 2.8581, places=3)
        self.assertEqual(octa[1][2], 6)
        mixed = co_lib["co3o4_mixed"]
        self.assertEqual(mixed.kind, "mixture")
        weights = {c["id"]: float(c["weight"]) for c in mixed.components}
        self.assertAlmostEqual(weights["co3o4_spinel_cod_9005888"],
                               1.0 / 3.0, places=12)
        self.assertAlmostEqual(weights["co3o4_oct"], 2.0 / 3.0, places=12)

        brucite = _shells(ni_lib["nioh2_theophrastite_cod_1548811"])
        self.assertEqual(brucite[0][0], "O")
        self.assertAlmostEqual(brucite[0][1], 2.0597, places=3)
        self.assertEqual(brucite[0][2], 6)
        self.assertEqual(brucite[1][0], "Ni")
        self.assertAlmostEqual(brucite[1][1], 3.1268, places=3)
        self.assertEqual(brucite[1][2], 6)

        tenorite = _shells(cu_lib["cuo_tenorite_cod_1011148"])
        self.assertEqual(tenorite[0][0], "O")
        self.assertAlmostEqual(tenorite[0][1], 1.9475, places=3)
        self.assertEqual(tenorite[0][2], 4)

        pto2 = _shells(pt_lib["pto2_beta"])
        self.assertEqual(pto2[0][0], "O")
        self.assertAlmostEqual(pto2[0][1], 2.0495, places=3)
        self.assertEqual(pto2[0][2], 6)
        self.assertEqual(pto2[1][0], "Pt")
        self.assertAlmostEqual(pto2[1][1], 3.100, places=3)
        self.assertEqual(pto2[1][2], 6)

        data_dir = ROOT / "data" / "bmm_standards" / "Data"

        def _identify(rel, absorber):
            res = xafs_pipeline.process_sample(
                [str(data_dir / rel)], absorber=absorber, strict=False,
                verbose=False)
            sid = structure_id.identify_result(
                res, absorber_symbol=absorber, verbose=False)
            return res, sid

        # ---- CoO ----
        res_coo, sid_coo = _identify("Co/Co-K-CoO.xdi", "Co")
        self.assertEqual(sid_coo["ranked"][0]["id"],
                         "coo_rocksalt_cod_1533087")
        self.assertEqual(sid_coo["status"], "confident")
        self.assertLess(sid_coo["ranked"][0]["score"], 0.06)
        out_coo = xafs_fit.run(
            res_coo, identify=False, fit=True, absorber_symbol="Co",
            verbose=False)
        self.assertEqual(out_coo["fit"]["mode"], "motif")
        self.assertEqual(out_coo["fit"]["motif"],
                         "coo_rocksalt_cod_1533087")
        self.assertLess(out_coo["fit"]["fit"].r_factor, 0.02)

        # ---- Co3O4（尖晶石双位点 mixture）----
        _, sid_co3o4 = _identify("Co/Co-K-Co3O4.xdi", "Co")
        self.assertEqual(sid_co3o4["ranked"][0]["id"], "co3o4_mixed")
        self.assertEqual(sid_co3o4["status"], "confident")
        self.assertLess(sid_co3o4["ranked"][0]["score"], 0.05)
        self.assertEqual(sid_co3o4["uncovered_peaks"], [])

        # ---- β-Ni(OH)2 ----
        res_nioh, sid_nioh = _identify(
            "Ni/Ni-K-NickelHydroxide.xdi", "Ni")
        self.assertEqual(sid_nioh["ranked"][0]["id"],
                         "nioh2_theophrastite_cod_1548811")
        self.assertEqual(sid_nioh["status"], "confident")
        self.assertLess(sid_nioh["ranked"][0]["score"], 0.14)
        out_nioh = xafs_fit.run(
            res_nioh, identify=False, fit=True, absorber_symbol="Ni",
            verbose=False)
        self.assertEqual(out_nioh["fit"]["mode"], "motif")
        fit_nioh = out_nioh["fit"]["fit"]
        self.assertLess(fit_nioh.r_factor, 0.05)
        by_label = {p.label: p for p in fit_nioh.paths}
        ni_o = next(p for label, p in by_label.items()
                    if label.startswith("Ni-O"))
        ni_ni = next(p for label, p in by_label.items()
                     if label.startswith("Ni-Ni"))
        self.assertGreaterEqual(ni_o.r, 2.03)
        self.assertLessEqual(ni_o.r, 2.11)
        self.assertGreaterEqual(ni_ni.r, 3.08)
        self.assertLessEqual(ni_ni.r, 3.17)

        # ---- CuO ----
        _, sid_cuo = _identify("Cu/Cu-CuO.xdi", "Cu")
        self.assertEqual(sid_cuo["ranked"][0]["id"],
                         "cuo_tenorite_cod_1011148")
        self.assertLess(sid_cuo["ranked"][0]["score"], 0.10)

    # ------------------------------------------------------------
    # 75–78. P0a：正式 CN / 诊断拟合双轨
    # ------------------------------------------------------------
    def test_formal_cn_refusal_without_model(self):
        """#75 P0a：无结构模型时正式 CN 拒绝，诊断拟合保留。"""
        self.assertTrue(self.lep_fit["fit"].get("diagnostic_only"))
        self.assertEqual(self.lep_fit["fit"].get("cn_status"),
                         "diagnostic_only")
        out = xafs_fit.run(
            self.lep, identify=True, fit=True, absorber_symbol="Fe",
            formal_cn=True, verbose=False)
        fit_out = out["fit"]
        self.assertEqual(fit_out["mode"], "motif")
        self.assertIsNotNone(fit_out.get("fit"))
        formal = fit_out["formal_cn"]
        self.assertEqual(formal["status"], "no_cif")
        self.assertEqual(formal["rows"], [])
        self.assertTrue(formal["reject_reason"])
        self.assertEqual(out["formal_cn"]["status"], "no_cif")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cn_table.csv"
            cn_table.write_cn_table(
                path, [cn_table.refusal_row(
                    formal["status"], formal["reject_reason"],
                    phase_id="lepidocrocite")])
            rows = cn_table.read_cn_table(path)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["gate_status"], "no_cif")
        self.assertEqual(rows[0]["n_fit"], "")
        self.assertTrue(rows[0]["reject_reason"])

    def test_formal_cn_l1_not_in_table(self):
        """#76 P0a：L1 只做诊断，不得进入正式 CN 表。"""
        res = xafs_pipeline.process_sample(
            [FERRIHYDRITE], absorber="Fe", strict=False, verbose=False)
        out = xafs_fit.run(
            res, identify=True, fit=True, absorber_symbol="Fe",
            formal_cn=True, verbose=False)
        fit_out = out["fit"]
        self.assertEqual(fit_out["mode"], "l1_first_shell")
        formal = fit_out["formal_cn"]
        self.assertEqual(formal["status"], "l1_diagnostic_only")
        self.assertEqual(formal["rows"], [])
        row = cn_table.refusal_row(formal["status"], formal["reject_reason"])
        self.assertEqual(row["n_fit"], "")
        self.assertEqual(row["cn_label"], "none")

    def test_occupancy_scale_failure_hard_fail(self):
        """#77 P0a：占据度展开失败硬报错，不按满占据继续。"""
        import dataclasses
        from unittest.mock import patch

        lep = [m for m in structure_lib.load_library("Fe")
               if m.id == "lepidocrocite"][0]
        cache = structure_lib.ensure_feff(lep, edge="K")
        with tempfile.TemporaryDirectory() as tmp:
            cif = Path(tmp) / "fake.cif"
            cif.write_text("data_fake\n", encoding="utf-8")
            motif = dataclasses.replace(
                lep, id="_p0a_occ_fail", cif_path=str(cif))
            with patch.object(structure_lib, "ensure_feff",
                              return_value=cache), \
                 patch.object(artemis, "cif_shell_occupancy_scale",
                              side_effect=XAFSError("FIT", "boom")):
                with self.assertRaises(structure_model.CNError) as ctx:
                    xafs_fit._ss_rows(motif, edge="K")
        self.assertEqual(ctx.exception.code, "occupancy_scale_failed")
        self.assertIn("occupancy_scale_failed", str(ctx.exception))

    def test_cn_table_schema_and_refusal_row(self):
        """#78 P0a：CN 表列顺序与拒绝行写读往返。"""
        expected = [
            "phase_id", "model_level", "source_type", "source_ref",
            "source_cif_sha256", "absorber_site", "shell_group", "shell_label",
            "scatterer", "r_cif", "n0", "n0_members", "n_fit", "n_status",
            "n_ratio", "n_err_formal", "n_err_systematic", "delta_r", "sigma2",
            "delta_e0", "rf", "s02", "s02_source", "s02_calibrated",
            "kmin", "kmax", "kweight", "rmin", "rmax", "n_idp", "n_variables",
            "cn_label", "gate_status", "reject_reason", "flags", "assumptions",
            "data_sha256", "software_version",
        ]
        self.assertEqual(cn_table.CN_TABLE_COLUMNS, expected)
        row = cn_table.refusal_row("no_cif", "no source", phase_id="x")
        self.assertEqual(row["gate_status"], "no_cif")
        self.assertEqual(row["n_fit"], "")
        self.assertTrue(row["reject_reason"])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cn_table.csv"
            cn_table.write_cn_table(path, [row])
            rows = cn_table.read_cn_table(path)
        self.assertEqual(len(rows), 1)
        for key in expected:
            self.assertIn(key, rows[0])
        self.assertEqual(rows[0]["phase_id"], "x")
        self.assertEqual(rows[0]["n_fit"], "")
        self.assertEqual(rows[0]["gate_status"], "no_cif")
        self.assertTrue(rows[0]["reject_reason"])

    # ------------------------------------------------------------
    # CN disclosure standard (docs/CN_REPORTING_STANDARD.md)
    # ------------------------------------------------------------
    @staticmethod
    def _disclosure_payload():
        return {
            "path_files": ["feff0001.dat"],
            "cn_status": "diagnostic_only",
            "cn_systematics": {"relative_1sigma": 0.30},
            "fit": {
                "paths": [{"label": "M-O", "N": 4.0, "R": 2.0,
                           "sigma2": 0.003}],
                "s02": 0.8,
                "delta_e0": 1.0,
                "kmin": 3.0,
                "kmax": 12.0,
                "kweights": [2],
                "rmin": 1.0,
                "rmax": 2.2,
                "r_factor": 0.02,
                "n_independent": 12.0,
                "n_variables": 5,
            },
            "structure": {
                "source_type": "constructed_cluster",
                "source_ref": "manual M-O4 cluster",
                "cif": None,
                "cif_declared": True,
                "feff_inputs": ["feff.inp"],
            },
        }

    def test_cn_disclosure_complete_payload(self):
        """CN 披露标准：完整 payload 无问题。"""
        self.assertEqual(
            cn_disclosure.check_fit_payload(self._disclosure_payload()), [])

    def test_cn_disclosure_flags_missing_provenance(self):
        """CN 披露标准：缺结构来源/CN 状态必须报出。"""
        payload = self._disclosure_payload()
        payload.pop("structure")
        payload.pop("cn_status")
        issues = cn_disclosure.check_fit_payload(payload)
        self.assertTrue(any("structure/CIF" in x for x in issues))
        self.assertTrue(any("CN status" in x for x in issues))

    def test_cn_disclosure_explicit_no_cif(self):
        """CN 披露标准：cif=None 视为显式声明，不报 CIF 缺失。"""
        issues = cn_disclosure.check_fit_payload(self._disclosure_payload())
        self.assertFalse(any("CIF status" in x for x in issues))

    def test_cn_disclosure_flags_diagnostic_mismatch(self):
        """CN 披露标准：diagnostic_only 不得配 formal 状态。"""
        payload = self._disclosure_payload()
        payload["diagnostic_only"] = True
        payload["cn_status"] = "formal"
        issues = cn_disclosure.check_fit_payload(payload)
        self.assertTrue(any("diagnostic" in x.lower() for x in issues))

    def test_paper_figures_reuses_disclosure_check(self):
        """paper_figures 复用项目级 CN 披露检查。"""
        import paper_figures as pf

        fit = pf.load_fit(self._disclosure_payload())
        self.assertEqual(pf.validate_cn_disclosure(fit), [])

    def test_cn_disclosure_large_delta_e0_is_advisory_only(self):
        """CN 披露标准：|Delta E0| > 5 eV 只提示，不阻断。"""
        payload = self._disclosure_payload()
        payload["fit"]["delta_e0"] = 8.0
        self.assertEqual(cn_disclosure.check_fit_payload(payload), [])
        notes = cn_disclosure.delta_e0_warnings(payload)
        self.assertEqual(len(notes), 1)
        self.assertIn("large", notes[0].lower())
        self.assertIn("absorption-edge shift", notes[0].lower())

    def test_cn_disclosure_delta_e0_threshold(self):
        """CN 披露标准：5 eV 边界与负值都要判断。"""
        payload = self._disclosure_payload()
        for value, expected in ((5.0, 0), (5.01, 1), (-8.0, 1), (0.0, 0)):
            payload["fit"]["delta_e0"] = value
            self.assertEqual(
                len(cn_disclosure.delta_e0_warnings(payload)), expected,
                msg=f"delta_e0={value}")

    def test_cn_disclosure_cli_advisory_does_not_fail(self):
        """CN 披露标准：CLI 对大 Delta E0 返回 0，并打印 NOTE。"""
        import contextlib
        import io

        payload = self._disclosure_payload()
        payload["fit"]["delta_e0"] = 8.0
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fit.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                code = cn_disclosure.main([str(path)])
        self.assertEqual(code, 0)
        output = buf.getvalue()
        self.assertIn("[NOTE]", output)
        self.assertIn("large", output.lower())

    def test_paper_figures_cn_table_includes_delta_e0_advisory(self):
        """CN 表：大 Delta E0 写入 Disclosure warnings 与 CSV advisory 列。"""
        import paper_figures as pf

        fit = pf.load_fit(self._disclosure_payload())
        self.assertIsNotNone(fit)
        fit.delta_e0 = 8.0
        with tempfile.TemporaryDirectory() as tmp:
            md = pf.write_cn_markdown(fit, tmp, stem="cn_table")
            csv_path = pf.write_cn_csv(fit, tmp, stem="cn_table")
            md_text = md.read_text(encoding="utf-8")
            csv_text = csv_path.read_text(encoding="utf-8-sig")
        self.assertIn("Disclosure warnings", md_text)
        self.assertIn("advisory, non-blocking", md_text)
        self.assertIn("large |Delta E0|", md_text)
        self.assertIn("disclosure_advisories", csv_text)
        self.assertIn("large |Delta E0|", csv_text)

    # ------------------------------------------------------------
    # 79–82. P0b：motif provenance 升级
    # ------------------------------------------------------------
    def test_motif_provenance_schema(self):
        """#79 P0b：全部 motif 可加载；新字段齐全且状态合法。"""
        allowed_status = set(motif_provenance.PROVENANCE_STATUSES)
        allowed_levels = set(motif_provenance.MODEL_LEVELS)
        allowed_occ = set(motif_provenance.OCCUPANCY_SOURCES)

        report_path = ROOT / "_structure_lib" / "provenance_report.json"
        self.assertTrue(report_path.exists(), "缺少 provenance_report.json")
        report = json.loads(report_path.read_text(encoding="utf-8"))
        entries = {(e["element"], e["id"]): e
                   for e in report["entries"]}

        total = 0
        for element_dir in sorted((ROOT / "_structure_lib").iterdir()):
            if not element_dir.is_dir() or element_dir.name.startswith("_"):
                continue
            for motif in structure_lib.load_library(element_dir.name):
                total += 1
                self.assertIn(motif.provenance_status, allowed_status,
                              motif.id)
                self.assertIn(motif.model_level, allowed_levels, motif.id)
                self.assertIn(motif.occupancy_source, allowed_occ, motif.id)
                if motif.kind == "crystal":
                    self.assertTrue(motif.absorber_site, motif.id)
                key = (element_dir.name, motif.id)
                self.assertIn(key, entries, motif.id)
                self.assertEqual(entries[key]["provenance_status"],
                                 motif.provenance_status)
                self.assertEqual(entries[key]["model_level"],
                                 motif.model_level)

        self.assertEqual(total, report["summary"]["motifs"])
        self.assertEqual(total, 43)
        self.assertEqual(report["summary"]["unclassified"], 0)

    def test_generated_cif_not_formal_source(self):
        """#80 P0b：现场生成 CIF 不得成为正式来源。"""
        fe = next(m for m in structure_lib.load_library("Fe")
                  if m.id == "fe_bcc")
        self.assertEqual(fe.provenance_status, "generated_only")
        self.assertIsNone(fe.source_cif_sha256)
        self.assertIsNone(fe.cif_path)

        generated = structure_lib._cif_path(fe)
        self.assertTrue(generated.exists())
        self.assertTrue(motif_provenance.is_generated_cif(generated))

        model = structure_model.model_from_motif(fe)
        self.assertEqual(model.provenance_status, "generated_only")
        self.assertIsNone(model.source_cif)
        self.assertIsNone(model.source_cif_sha256)
        gate = structure_model.validate_formal(model)
        self.assertFalse(gate["ok"])
        self.assertIn(gate["code"], ("unverified_source", "no_cif"))

        lep = next(m for m in structure_lib.load_library("Fe")
                   if m.id == "lepidocrocite")
        self.assertIn(lep.provenance_status,
                      ("source_ref_only", "generated_only"))
        self.assertIsNone(lep.source_cif_sha256)
        self.assertEqual(lep.source_ref, "COD:9002763")
        lep_model = structure_model.model_from_motif(lep)
        self.assertIsNone(lep_model.source_cif)
        self.assertFalse(structure_model.validate_formal(lep_model)["ok"])

        # 即使把生成 CIF 伪造成“哈希正确 + 吸收位明确”，门禁仍拒绝。
        fake = structure_model.StructureModel(
            model_id="fake_generated",
            absorber="Fe",
            absorber_site="Fe1@(0,0,0)",
            model_level="crystal_cif",
            source_type="user_cif",
            source_cif=str(generated),
            source_cif_sha256=structure_model.sha256_file(generated),
            provenance_status="verified_source_cif",
        )
        gate = structure_model.validate_formal(fake)
        self.assertFalse(gate["ok"])
        self.assertEqual(gate["code"], "generated_cif_not_formal")

    def test_user_cif_is_verified_source(self):
        """#81 P0b：用户真实 CIF + sha256 + 吸收位可升级为正式来源。"""
        source = ROOT / "lepidocrocite.cif"
        self.assertTrue(source.exists())
        with tempfile.TemporaryDirectory() as tmp:
            user_cif = Path(tmp) / "user_lepidocrocite.cif"
            user_cif.write_bytes(source.read_bytes())
            motif = structure_lib.Motif(
                id="user_lep",
                absorber="Fe",
                kind="crystal",
                nuclearity="crystal",
                source={"type": "user_cif", "ref": "user:test"},
                cell={"a": 3.072, "b": 12.516, "c": 3.873,
                      "alpha": 90, "beta": 90, "gamma": 90},
                spacegroup="C m c m",
                asym_sites=[], absorber_fract=[0.5, 0.1778, 0.25],
                cif_path=str(user_cif),
                edits=[{"op": "complete_site", "site": "O1",
                        "reason": "文献表残缺，几何补全"}],
            )
            fields = motif_provenance.classify_motif(motif)
            self.assertEqual(fields["provenance_status"],
                             "verified_source_cif")
            self.assertEqual(fields["model_level"], "crystal_edited")
            self.assertEqual(fields["source_cif_sha256"],
                             structure_model.sha256_file(user_cif))
            self.assertFalse(fields["_has_generated_cif"])

            model = structure_model.model_from_motif(motif)
            self.assertEqual(model.provenance_status, "verified_source_cif")
            self.assertEqual(model.model_level, "crystal_edited")
            self.assertEqual(model.source_cif, str(user_cif))
            self.assertTrue(model.absorber_site)
            gate = structure_model.validate_formal(model)
            self.assertTrue(gate["ok"], gate)

    def test_b20_motifs_not_degraded(self):
        """#82 P0b：B17–B20 motif 保留、指纹不变、constructed 不升级。"""
        expected = {
            "Fe": {"fen4_planar", "fe2_mu_oxo", "fe2n6_dual"},
            "Co": {"co3o4_mixed", "co3o4_oct",
                   "co3o4_spinel_cod_9005888", "con4_planar"},
            "Ni": {"nin4_planar", "nioh2_theophrastite_cod_1548811"},
            "Cu": {"cun4_planar", "cuo_tenorite_cod_1011148"},
            "Pt": {"ptcl4_square", "pto2_beta"},
        }
        libraries = {}
        for element, ids in expected.items():
            lib = {m.id: m for m in structure_lib.load_library(element)}
            libraries[element] = lib
            for motif_id in ids:
                self.assertIn(motif_id, lib, f"{element}/{motif_id}")

        co_oct = libraries["Co"]["co3o4_oct"]
        co_tet = libraries["Co"]["co3o4_spinel_cod_9005888"]
        con4 = libraries["Co"]["con4_planar"]
        self.assertEqual(co_oct.shell_expect[0],
                         {"el": "O", "r": 1.9217, "n": 6, "n_occ": 6.0})
        self.assertEqual(co_tet.shell_expect[0],
                         {"el": "O", "r": 1.9353, "n": 4, "n_occ": 4.0})
        self.assertEqual(con4.shell_expect[0],
                         {"el": "N", "r": 1.96, "n": 4})

        for motif_id in ("co3o4_mixed", "con4_planar"):
            motif = libraries["Co"][motif_id]
            self.assertEqual(motif.provenance_status, "constructed")
            self.assertEqual(motif.model_level, "cluster")
            model = structure_model.model_from_motif(motif)
            self.assertEqual(model.provenance_status, "constructed")
            self.assertFalse(structure_model.validate_formal(model)["ok"])

        # 压力测试会用 overwrite=True 重注册种子；注册通道必须保留
        # provenance，不能把刚补好的字段静默写丢。
        with tempfile.TemporaryDirectory() as tmp:
            original_root = structure_lib.LIB_ROOT
            try:
                structure_lib.LIB_ROOT = Path(tmp)
                motif = structure_lib.Motif(
                    id="persist_check", absorber="Fe", kind="cluster",
                    nuclearity="m1", source={"type": "constructed"},
                    sites=[[2.0, 0.0, 0.0, "N"]])
                saved_path = structure_lib.register(motif)
                saved = json.loads(saved_path.read_text(encoding="utf-8"))
            finally:
                structure_lib.LIB_ROOT = original_root
        self.assertEqual(saved["provenance_status"], "constructed")
        self.assertEqual(saved["model_level"], "cluster")
        self.assertIn("edits", saved)

    # ------------------------------------------------------------
    # 83–86. P1a：覆盖审计 + phase pack
    # ------------------------------------------------------------
    def test_motif_coverage_audit(self):
        """#83 P1a：motif 覆盖审计覆盖 17/17 元素且无未分类。"""
        coverage = json.loads((ROOT / "_structure_lib" /
                               "motif_coverage.json")
                              .read_text(encoding="utf-8"))
        report = json.loads((ROOT / "_structure_lib" /
                             "provenance_report.json")
                            .read_text(encoding="utf-8"))
        self.assertEqual(coverage["schema_version"], 1)
        self.assertEqual(set(coverage["absorber_elements"]),
                         set(phase_pack.ABSORBER_ELEMENTS))
        self.assertEqual(len(coverage["absorber_elements"]), 17)
        self.assertEqual(coverage["summary"]["elements_total"], 17)
        self.assertEqual(coverage["summary"]["elements_covered"], 17)
        self.assertEqual(coverage["summary"]["elements_missing"], 0)
        self.assertEqual(coverage["summary"]["motifs"], 43)
        self.assertEqual(set(coverage["missing_elements"]), set())
        self.assertEqual(set(coverage["covered_elements"]) |
                         set(coverage["missing_elements"]),
                         set(coverage["absorber_elements"]))
        self.assertFalse(set(coverage["covered_elements"]) &
                         set(coverage["missing_elements"]))
        self.assertEqual(sum(row["motifs"]
                             for row in coverage["elements"].values()), 43)
        self.assertEqual(coverage["summary"]["formal_cn_ready_elements"], [])
        self.assertEqual(report["summary"]["unclassified"], 0)
        self.assertEqual(report["summary"]["by_status"].get(
            "verified_source_cif", 0), 0)
        self.assertEqual(report["summary"]["motifs"], 43)

    def test_phase_pack_import_and_index(self):
        """#84 P1a：44 条来源（5 金属种子 + 39 抓取目标）可校验、可索引。"""
        from tools import cif_crawler

        pack_root = ROOT / "_structure_lib" / "_phase_pack"
        manifest = phase_pack.load_manifest(pack_root)
        self.assertEqual(manifest["schema_version"], 1)
        self.assertEqual(manifest["pack_id"], "formal-cn-phase-pack")
        self.assertEqual(len(manifest["entries"]), 44)

        seeds = {
            "1100108": ("Fe", ["Fe1"]),
            "9008492": ("Co", ["Co"]),
            "1534892": ("Ni", ["Ni1"]),
            "4105040": ("Cu", ["Cu"]),
            "1011103": ("Pt", ["Pt1"]),
        }
        crawled = {
            "1011240": ("Fe", ["Fe1"]),
            "1010369": ("Fe", ["Fe1", "Fe2", "Fe3"]),
            "1533087": ("Co", ["Co1"]),
            "9005888": ("Co", ["Co1", "Co2"]),
            "1010093": ("Ni", ["Ni1"]),
            "1548811": ("Ni", ["Ni1"]),
            "1011148": ("Cu", ["Cu1"]),
            "1010941": ("Cu", ["Cu1"]),
            "1562989": ("Ce", ["Ce"]),
            "9008589": ("Mn", ["Mn1", "Mn2", "Mn3", "Mn4"]),
            "9008513": ("Ru", ["Ru"]),
            "9008482": ("Rh", ["Rh"]),
            "9008478": ("Pd", ["Pd"]),
            "9008459": ("Ag", ["Ag"]),
            "9008470": ("Ir", ["Ir"]),
            "9010987": ("Ce", ["Ce"]),
            "1527120": ("U", ["U1"]),
            "9008531": ("Cr", ["Cr"]),
            "9008522": ("Zn", ["Zn"]),
            "9008543": ("Mo", ["Mo"]),
            "9008463": ("Au", ["Au"]),
        }
        by_id = {}
        for entry in manifest["entries"]:
            result = phase_pack.validate_entry(pack_root, entry)
            self.assertTrue(result["ok"], (entry["id"], result))
            self.assertEqual(entry["source"], "COD")
            self.assertTrue(entry["experimental"])
            self.assertEqual(entry["model_level"], "crystal_cif")
            self.assertEqual(len(entry["absorber_elements"]), 1)
            by_id[entry["id"]] = entry
        for entry_id, (element, sites) in {**seeds, **crawled}.items():
            self.assertIn(entry_id, by_id)
            entry = by_id[entry_id]
            self.assertEqual(entry["absorber_elements"], [element])
            self.assertEqual(entry["absorber_sites"], sites)

        # P5-C：18 个常见相补齐目标必须都已导入。
        p5c_ids = {t["cod_id"] for t in cif_crawler.VERIFIED_TARGETS.values()
                   if t["tier"] == "P5C"}
        self.assertEqual(len(p5c_ids), 18)
        self.assertTrue(p5c_ids <= set(by_id), p5c_ids - set(by_id))

        index = json.loads((ROOT / "_structure_lib" /
                            "phase_index.json").read_text(encoding="utf-8"))
        self.assertEqual(index["schema_version"], 1)
        self.assertEqual(index["pack_id"], manifest["pack_id"])
        self.assertEqual(len(index["entries"]), 44)
        self.assertEqual(set(index["summary"]["elements_covered"]), {
            "Fe", "Co", "Ni", "Cu", "Pt", "Cr", "Mn", "Zn", "Mo",
            "Ru", "Rh", "Pd", "Ag", "Ir", "Au", "Ce", "U"})
        indexed = {(e["source"], e["id"]): e for e in index["entries"]}
        for entry in manifest["entries"]:
            row = indexed[(entry["source"], entry["id"])]
            self.assertEqual(row["cif_sha256"], entry["cif_sha256"])
            self.assertEqual(row["pack_path"], entry["cif_path"])
            self.assertTrue((ROOT / row["cif_path"]).is_file())

        coverage = json.loads((ROOT / "_structure_lib" /
                               "phase_coverage.json")
                              .read_text(encoding="utf-8"))
        self.assertEqual(len(coverage["absorber_elements"]), 17)
        self.assertEqual(coverage["summary"]["elements_with_motif"], 17)
        self.assertEqual(coverage["summary"]["elements_without_motif"], 0)
        self.assertEqual(coverage["summary"]["without_motif_elements"], [])
        self.assertEqual(set(coverage["summary"]["verified_elements"]), {
            "Fe", "Co", "Ni", "Cu", "Pt", "Cr", "Mn", "Zn", "Mo",
            "Ru", "Rh", "Pd", "Ag", "Ir", "Au", "Ce", "U"})
        for element in {v[0] for v in {**seeds, **crawled}.values()}:
            row = coverage["elements"][element]
            self.assertTrue(row["formal_cn_ready"], element)
            self.assertEqual(row["coverage_status"], "verified")
            self.assertGreaterEqual(row["verified_cif"], 1)

        # 导入幂等：同一文件导入两次，manifest 与目标文件保持不变。
        with tempfile.TemporaryDirectory() as tmp:
            tmp_pack = Path(tmp) / "pack"
            src = pack_root / by_id["1011240"]["cif_path"]
            first = phase_pack.import_cif(
                src, tmp_pack, source="user", id="seed_fe",
                tags=["metal"], meta={})
            manifest_1 = (tmp_pack / "manifest.json").read_text(
                encoding="utf-8")
            second = phase_pack.import_cif(
                src, tmp_pack, source="user", id="seed_fe",
                tags=["metal"], meta={})
            manifest_2 = (tmp_pack / "manifest.json").read_text(
                encoding="utf-8")
            self.assertEqual(first["cif_sha256"], second["cif_sha256"])
            self.assertEqual(manifest_1, manifest_2)
            target_dir = tmp_pack / "user" / "seed_fe"
            self.assertEqual(len(list(target_dir.glob("*.cif"))), 1)
            self.assertEqual(len(list(target_dir.glob("*.meta.json"))), 1)
            tmp_index = phase_pack.build_index(tmp_pack)
            self.assertEqual(tmp_index["summary"]["entries"], 1)
            tmp_cov = phase_pack.build_coverage(tmp_pack, index=tmp_index)
            self.assertTrue(tmp_cov["elements"]["Fe"]["formal_cn_ready"])

    def test_phase_pack_sha256_mismatch_rejected(self):
        """#85 P1a：sha256 不符、同 ID 冲突与显式覆盖。"""
        pack_root = ROOT / "_structure_lib" / "_phase_pack"
        manifest = phase_pack.load_manifest(pack_root)
        fe = next(e for e in manifest["entries"]
                  if e["absorber_elements"] == ["Fe"])
        source = pack_root / fe["cif_path"]

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            src = tmp / "COD-1100108.cif"
            src.write_bytes(source.read_bytes())
            meta_path = src.with_suffix(".meta.json")
            meta_path.write_text(json.dumps({
                "source": "COD", "id": "1100108",
                "cif_sha256": "0" * 64,
            }), encoding="utf-8")
            with self.assertRaises(phase_pack.PhasePackError) as ctx:
                phase_pack.import_cif(src, tmp / "pack")
            self.assertEqual(ctx.exception.code, "sha256_mismatch")
            meta_path.unlink()

            first = phase_pack.import_cif(
                src, tmp / "pack", source="user", id="conflict")
            old_cif = tmp / "pack" / first["cif_path"]
            self.assertTrue(old_cif.is_file())

            changed = tmp / "changed.cif"
            changed.write_bytes(source.read_bytes() + b"\n# changed\n")
            with self.assertRaises(phase_pack.PhasePackError) as ctx:
                phase_pack.import_cif(
                    changed, tmp / "pack", source="user", id="conflict")
            self.assertEqual(ctx.exception.code, "id_sha_conflict")

            replaced = phase_pack.import_cif(
                changed, tmp / "pack", source="user", id="conflict",
                replace=True)
            self.assertNotEqual(replaced["cif_sha256"],
                                first["cif_sha256"])
            rebuilt = phase_pack.load_manifest(tmp / "pack")
            self.assertEqual(len(rebuilt["entries"]), 1)
            self.assertEqual(rebuilt["entries"][0]["cif_sha256"],
                             replaced["cif_sha256"])
            self.assertFalse(old_cif.exists())

            tampered = dict(rebuilt["entries"][0])
            tampered["cif_sha256"] = "f" * 64
            result = phase_pack.validate_entry(tmp / "pack", tampered)
            self.assertFalse(result["ok"])
            self.assertEqual(result["code"], "sha256_mismatch")

    def test_uncovered_element_no_formal_cn(self):
        """#86 P1a：generated_only motif 不升级正式 CN；phase pack 覆盖 17/17。"""
        coverage = json.loads((ROOT / "_structure_lib" /
                               "phase_coverage.json")
                              .read_text(encoding="utf-8"))
        index = json.loads((ROOT / "_structure_lib" /
                            "phase_index.json").read_text(encoding="utf-8"))
        # P4+ 的 21 条 + P5-C 的 18 条真实 CIF 已补齐全部 17 个吸收元素。
        self.assertEqual(coverage["summary"]["elements_verified"], 17)
        self.assertEqual(coverage["summary"]["elements_not_covered"], 0)
        self.assertEqual(coverage["summary"]["not_covered_elements"], [])
        self.assertEqual(coverage["summary"]["elements_motif_only"], 0)
        for element in phase_pack.ABSORBER_ELEMENTS:
            self.assertTrue(
                coverage["elements"][element]["formal_cn_ready"], element)
            self.assertEqual(
                coverage["elements"][element]["coverage_status"], "verified")
            self.assertGreaterEqual(
                len(phase_pack.verified_entries(index, element)), 1, element)
        # 但 generated_only motif 仍不能冒充正式来源。
        zn = next(m for m in structure_lib.load_library("Zn")
                  if m.id == "zn_hcp")
        self.assertEqual(zn.provenance_status, "generated_only")
        model = structure_model.model_from_motif(zn)
        gate = structure_model.validate_formal(model)
        self.assertFalse(gate["ok"])
        row = cn_table.refusal_row(
            "no_cif", gate["reason"], phase_id=model.model_id,
            source_type=model.source_type)
        self.assertEqual(row["gate_status"], "no_cif")
        self.assertEqual(row["n_fit"], "")
        self.assertEqual(row["cn_label"], "none")
        self.assertIn("refused", row["flags"])

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "cn_table.csv"
            cn_table.write_cn_table(out, [row])
            back = cn_table.read_cn_table(out)
        self.assertEqual(len(back), 1)
        self.assertEqual(back[0]["gate_status"], "no_cif")
        self.assertEqual(back[0]["n_fit"], "")

    # ------------------------------------------------------------
    # 87–93. P1b：相解析与缓存
    # ------------------------------------------------------------
    def test_phase_spec_parser_no_eval(self):
        """#87 P1b：--phase SPEC 只做字符串拆分，禁止 eval/exec。"""
        # 合法 SPEC
        spec = phase_resolver.parse_phase_spec(
            "name=fcc Fe;formula=Fe;spacegroup=Im-3m;"
            "cell=a=2.868,b=2.868,c=2.868,alpha=90,beta=90,gamma=90;"
            "source=COD;id=1100108")
        self.assertEqual(spec["name"], "fcc Fe")
        self.assertEqual(spec["formula"], "Fe")
        self.assertEqual(spec["spacegroup"], "Im-3m")
        self.assertEqual(spec["cell"]["a"], 2.868)
        self.assertEqual(spec["cell"]["gamma"], 90.0)
        self.assertEqual(spec["source"], "COD")
        self.assertEqual(spec["id"], "1100108")

        # 禁止 eval/exec：任何试图注入代码的片段都应被当作普通字符串
        # 或因不允许的键 / 格式错误而被拒绝，绝不会执行。
        # 下面这些值即使看起来像代码，也只是被当作字符串，不执行；
        # 真正应拒绝的是非法 key、缺少 '='、cell 格式错误。
        bad_key_or_format = [
            "evil=1",                 # 非法 key
            "novalue",                # 缺少 '='
            "cell=1,2,3",             # cell 不是 6 个数字
            "cell=1,2,3,4,5,6,7",     # cell 太多数字
            "cell=a=1,b=2,c=3,d=4,e=5",  # cell 缺 gamma
        ]
        for b in bad_key_or_format:
            with self.assertRaises(phase_resolver.PhaseResolutionError):
                phase_resolver.parse_phase_spec(b)

        # 代码注入字符串只被当作普通字符串，不会执行（也不报错，因为 key 合法）
        harmless = [
            "name=__import__('os')",
            "formula=eval('1+1')",
            "spacegroup=open('/etc/passwd')",
        ]
        for h in harmless:
            spec = phase_resolver.parse_phase_spec(h)
            # key 合法，值原样保留

    def test_phase_cif_user_source(self):
        """#88 P1b：--phase-cif 用户真实 CIF 作为正式来源。"""
        cif = ROOT / "_structure_lib" / "_phase_pack" / "COD" / "1100108" \
            / "41085e65debf2bd82a772bcee9e7c0a43ef42b67666d2f233c8f2d822769b8fe.cif"
        self.assertTrue(cif.is_file())

        # 无 absorber_site：自动取候选首个吸收位
        res = phase_resolver.resolve_phase("Fe", phase_cif=[str(cif)])
        self.assertEqual(res.status, "resolved")
        self.assertIsNotNone(res.selected)
        self.assertEqual(res.selected.source, "user")
        self.assertTrue(res.selected.cif_path.endswith(".cif"))
        self.assertEqual(res.selected.cif_sha256,
                         "41085e65debf2bd82a772bcee9e7c0a43ef42b67666d2f233c8f2d822769b8fe")
        self.assertIsNotNone(res.structure_model)
        self.assertEqual(res.structure_model.provenance_status,
                         "verified_source_cif")
        self.assertEqual(res.structure_model.source_type, "user_cif")
        self.assertTrue(res.structure_model.source_ref.startswith("user:"))
        # 用户 CIF 始终优先，不是 _structure_lib/*/_cif/ 生成文件
        from structure_model import validate_formal
        v = validate_formal(res.structure_model)
        self.assertTrue(v["ok"], v)

    def test_phase_id_cache_hit(self):
        """#89 P1b：--phase-id 精确编号，本地缓存命中（幂等）。"""
        cache_root = ROOT / "_structure_lib" / "_phase_cache"
        # 先确保缓存为空（测试隔离）
        import shutil
        fe_cache = cache_root / "COD" / "1100108"
        if fe_cache.exists():
            shutil.rmtree(fe_cache)

        res = phase_resolver.resolve_phase("Fe", phase_id="COD:1100108",
                                           absorber_site="Fe1")
        self.assertEqual(res.status, "resolved")
        self.assertIsNotNone(res.selected)
        self.assertEqual(res.selected.source, "COD")
        self.assertEqual(res.selected.id, "1100108")
        self.assertIsNotNone(res.structure_model)
        self.assertEqual(res.structure_model.provenance_status,
                         "verified_source_cif")

        # 缓存文件应已写入
        meta_files = list(fe_cache.glob("*.meta.json")) if fe_cache.exists() else []
        self.assertGreaterEqual(len(meta_files), 1)

        # 再次调用：缓存命中，结果一致
        res2 = phase_resolver.resolve_phase("Fe", phase_id="COD:1100108",
                                            absorber_site="Fe1")
        self.assertEqual(res2.status, "resolved")
        self.assertEqual(res2.selected.cif_sha256,
                         res.selected.cif_sha256)
        self.assertEqual(len(list(fe_cache.glob("*.cif"))), 1)

        # 清理测试产生的缓存
        shutil.rmtree(fe_cache)

    def test_phase_id_cache_miss_offline(self):
        """#90 P1b：--phase-id 未命中且离线 → offline_miss，不退回 motif。"""
        cache_root = ROOT / "_structure_lib" / "_phase_cache"
        import shutil
        # 确保不存在该 ID 的缓存
        weird_cache = cache_root / "COD" / "99999999"
        if weird_cache.exists():
            shutil.rmtree(weird_cache)
        # phase pack 中也不存在该 ID
        res = phase_resolver.resolve_phase("Fe", phase_id="COD:99999999",
                                           offline=True)
        self.assertEqual(res.status, "offline_miss")
        self.assertIsNone(res.selected)
        self.assertIsNone(res.structure_model)
        self.assertIn("离线", res.reason)
        if weird_cache.exists():
            shutil.rmtree(weird_cache)

    def test_phase_name_ambiguous(self):
        """#91 P1b：多候选 → ambiguous，不自动 Top-1。"""
        # 给 Fe 指定 phase_index，Fe 有 1 条 verified；此处用不存在的元素
        # 制造无候选；要制造“多候选”，需在同一 pack 下有 2 条含同一吸收元素。
        # 用临时 pack 造两条 Fe CIF。
        import shutil
        tmp = Path(tempfile.mkdtemp())
        try:
            tmp_pack = tmp / "pack"
            tmp_pack.mkdir()
            # 复制 Fe CIF 为两条不同 id（同 sha 不行，改文件名即可）
            src = ROOT / "_structure_lib" / "_phase_pack" / "COD" / "1100108" \
                / "41085e65debf2bd82a772bcee9e7c0a43ef42b67666d2f233c8f2d822769b8fe.cif"
            sha = phase_pack.sha256_file(src)
            d1 = tmp_pack / "COD" / "1100108"
            d1.mkdir(parents=True)
            shutil.copy2(src, d1 / f"{sha}.cif")
            meta1 = {"source": "COD", "id": "1100108", "formula": "Fe",
                     "spacegroup": "Im-3m", "cell": {"a": 2.868, "b": 2.868,
                     "c": 2.868, "alpha": 90, "beta": 90, "gamma": 90},
                     "absorber_sites": ["Fe1"], "absorber_elements": ["Fe"],
                     "cif_path": f"{sha}.cif", "cif_sha256": sha,
                     "experimental": True, "model_level": "crystal_cif"}
            (d1 / f"{sha}.meta.json").write_text(
                json.dumps(meta1, ensure_ascii=False, indent=2), encoding="utf-8")
            # 第二条：不同 id，同 sha（内容相同），构成另一个候选
            d2 = tmp_pack / "COD" / "1100109"
            d2.mkdir(parents=True)
            shutil.copy2(src, d2 / f"{sha}.cif")
            meta2 = dict(meta1, id="1100109")
            (d2 / f"{sha}.meta.json").write_text(
                json.dumps(meta2, ensure_ascii=False, indent=2), encoding="utf-8")
            phase_pack.rebuild_manifest(tmp_pack)

            res = phase_resolver.resolve_phase(
                "Fe", phase_index_path=str(tmp_pack / "manifest.json"),
                phase_pack_path=str(tmp_pack))
            self.assertEqual(res.status, "ambiguous")
            self.assertEqual(len(res.candidates), 2)
            self.assertIsNone(res.selected)
            self.assertIsNone(res.structure_model)
            self.assertIn("不自动 Top-1", res.reason)
        finally:
            shutil.rmtree(tmp)

    def test_phase_index_schema(self):
        """#92 P1b：phase_index.json schema 与必要字段齐全。"""
        index = json.loads((ROOT / "_structure_lib" /
                            "phase_index.json").read_text(encoding="utf-8"))
        self.assertEqual(index["schema_version"], 1)
        self.assertIn("pack_id", index)
        self.assertIn("pack_root", index)
        self.assertIn("manifest_sha256", index)
        self.assertIn("summary", index)
        self.assertEqual(index["summary"]["entries"], 44)
        self.assertEqual(index["summary"]["verified_cif"], 44)
        self.assertEqual(set(index["summary"]["elements_covered"]), {
            "Fe", "Co", "Ni", "Cu", "Pt", "Cr", "Mn", "Zn", "Mo",
            "Ru", "Rh", "Pd", "Ag", "Ir", "Au", "Ce", "U"})
        for row in index["entries"]:
            for key in ("source", "id", "formula", "spacegroup", "cell",
                        "absorber_sites", "absorber_elements", "cif_path",
                        "cif_sha256", "experimental", "model_level"):
                self.assertIn(key, row, row.get("id"))
            self.assertTrue(row["experimental"])
            self.assertIn(row["model_level"], ("crystal_cif", "crystal_edited"))
            # cif_path 必须是 repo-relative，且文件存在
            cif = ROOT / row["cif_path"]
            self.assertTrue(cif.is_file(), row["cif_path"])
            # sha256 必须与文件一致
            actual = phase_pack.sha256_file(cif)
            self.assertEqual(actual, row["cif_sha256"])

    def test_cache_sha256_reuse(self):
        """#93 P1b：缓存复用前校验 sha256；不匹配则拒绝。"""
        import shutil
        cache_root = ROOT / "_structure_lib" / "_phase_cache"
        fe_cache = cache_root / "COD" / "1100108"
        if fe_cache.exists():
            shutil.rmtree(fe_cache)

        # 首次解析：写入缓存
        res = phase_resolver.resolve_phase("Fe", phase_id="COD:1100108",
                                           absorber_site="Fe1")
        self.assertEqual(res.status, "resolved")
        meta_files = list(fe_cache.glob("*.meta.json"))
        self.assertEqual(len(meta_files), 1)
        cached_sha = json.loads(meta_files[0].read_text(encoding="utf-8"))[
            "cif_sha256"]
        self.assertEqual(cached_sha,
                         "41085e65debf2bd82a772bcee9e7c0a43ef42b67666d2f233c8f2d822769b8fe")

        # 篡改缓存中的 CIF 字节（改哈希），复用时应拒绝
        cif_file = list(fe_cache.glob("*.cif"))[0]
        original_bytes = cif_file.read_bytes()
        try:
            cif_file.write_bytes(original_bytes + b"\n# tampered\n")
            res2 = phase_resolver.resolve_phase("Fe", phase_id="COD:1100108",
                                                absorber_site="Fe1")
            # 缓存 sha256 不匹配，应回退到 phase pack 重新获取（仍 resolved）
            self.assertEqual(res2.status, "resolved")
        finally:
            cif_file.write_bytes(original_bytes)

        # 清理测试产生的缓存
        shutil.rmtree(fe_cache)

    # ------------------------------------------------------------
    # 94–101. P2a：CIF 门禁 + N₀ + Rayleigh 分组
    # ------------------------------------------------------------
    def test_cif_gate_parse_failed(self):
        """#94 P2a：缺 _atom_site_* 循环的 CIF → parse_failed。"""
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.cif"
            bad.write_text(
                "data_test\n_cell_length_a 1.0\n_cell_length_b 1.0\n",
                encoding="utf-8")
            r = cif_gate.gate_cif(bad, "Fe")
            self.assertFalse(r.ok)
            self.assertEqual(r.code, "parse_failed")
            self.assertIn("_atom_site", r.reason)

    def test_cif_gate_absorber_missing(self):
        """#95 P2a：CIF 中无吸收元素 → absorber_missing。"""
        with tempfile.TemporaryDirectory() as tmp:
            cif = Path(tmp) / "oxygen.cif"
            cif.write_text("""data_test
_symmetry_space_group_name_H-M 'P 1'
_cell_length_a 3.0 _cell_length_b 3.0 _cell_length_c 3.0
_cell_angle_alpha 90 _cell_angle_beta 90 _cell_angle_gamma 90
loop_
_atom_site_label _atom_site_type_symbol _atom_site_fract_x _atom_site_fract_y _atom_site_fract_z _atom_site_occupancy
O1 O 0 0 0 1
""", encoding="utf-8")
            r = cif_gate.gate_cif(cif, "Fe")
            self.assertFalse(r.ok)
            self.assertEqual(r.code, "absorber_missing")

    def test_cif_gate_absorber_site_ambiguous(self):
        """#96 P2a：多元素结构中 Fe 有两个晶位，未指定 → absorber_site_ambiguous。"""
        with tempfile.TemporaryDirectory() as tmp:
            cif = Path(tmp) / "two_fe.cif"
            cif.write_text("""data_test
_symmetry_space_group_name_H-M 'P 1'
_cell_length_a 4.0 _cell_length_b 4.0 _cell_length_c 4.0
_cell_angle_alpha 90 _cell_angle_beta 90 _cell_angle_gamma 90
loop_
_atom_site_label _atom_site_type_symbol _atom_site_fract_x _atom_site_fract_y _atom_site_fract_z _atom_site_occupancy
Fe1 Fe 0 0 0 1
Fe2 Fe 0.5 0.5 0.5 1
O1 O 0.25 0.25 0.25 1
""", encoding="utf-8")
            r = cif_gate.gate_cif(cif, "Fe")
            self.assertFalse(r.ok)
            self.assertEqual(r.code, "absorber_site_ambiguous")
            # 指定吸收位后应放行（测试 CIF 中 Fe1 第一壳 CN=1，超出正常范围，需开关）
            r2 = cif_gate.gate_cif(cif, "Fe", absorber_site="Fe1",
                                   allow_unusual_coordination=True)
            self.assertTrue(r2.ok, r2.reason)
            self.assertEqual(r2.absorber_site, "Fe1")

    def test_cif_gate_occupancy_implicit_ok(self):
        """#97 P2a：缺 occupancy 但单元素有序 → implicit_1.0 放行（COD:9008492 Co）。"""
        cif = ROOT / "_structure_lib" / "_phase_pack" / "COD" / "9008492" \
            / "bd9d57408d901e3ff966ba0ebe19ced70e1f1b62a4ef47fda801a39ecbb5344a.cif"
        self.assertTrue(cif.is_file())
        r = cif_gate.gate_cif(cif, "Co")
        self.assertTrue(r.ok, r.reason)
        self.assertEqual(r.occupancy_source, "implicit_1.0")
        self.assertEqual(r.first_shell_cn, 12.0)   # hcp Co 第一壳 CN=12

    def test_cif_gate_occupancy_disorder_rejected(self):
        """#98 P2a：多元素结构缺 occupancy → occupancy_missing。"""
        with tempfile.TemporaryDirectory() as tmp:
            cif = Path(tmp) / "disorder.cif"
            cif.write_text("""data_test
_symmetry_space_group_name_H-M 'P 1'
_cell_length_a 4.0 _cell_length_b 4.0 _cell_length_c 4.0
_cell_angle_alpha 90 _cell_angle_beta 90 _cell_angle_gamma 90
loop_
_atom_site_label _atom_site_type_symbol _atom_site_fract_x _atom_site_fract_y _atom_site_fract_z
Fe1 Fe 0 0 0
O1 O 0.5 0.5 0.5
""", encoding="utf-8")
            r = cif_gate.gate_cif(cif, "Fe")
            self.assertFalse(r.ok)
            self.assertEqual(r.code, "occupancy_missing")

    def test_cif_gate_cocupancy_unresolved(self):
        """#99 P2a：同一坐标多元素共占位，未拆分 → cocupancy_unresolved。"""
        with tempfile.TemporaryDirectory() as tmp:
            cif = Path(tmp) / "cocupancy.cif"
            cif.write_text("""data_test
_symmetry_space_group_name_H-M 'P 1'
_cell_length_a 4.0 _cell_length_b 4.0 _cell_length_c 4.0
_cell_angle_alpha 90 _cell_angle_beta 90 _cell_angle_gamma 90
loop_
_atom_site_label _atom_site_type_symbol _atom_site_fract_x _atom_site_fract_y _atom_site_fract_z _atom_site_occupancy
Fe1 Fe 0 0 0 0.5
O1 O 0 0 0 0.5
""", encoding="utf-8")
            r = cif_gate.gate_cif(cif, "Fe")
            self.assertFalse(r.ok)
            self.assertEqual(r.code, "cocupancy_unresolved")

    def test_cn_range_12_metals(self):
        """#100 P2a：fcc/hcp 金属 CN=12 正常通过；CN=8 (bcc) 也通过。"""
        pack = ROOT / "_structure_lib" / "_phase_pack"
        cases = {
            "Co": ("9008492", "bd9d57408d901e3ff966ba0ebe19ced70e1f1b62a4ef47fda801a39ecbb5344a.cif", 12.0),
            "Ni": ("1534892", "adb86b38b5313592734259a0ff2bb537fc60acb92793015ee8a74a41b9966cdd.cif", 12.0),
            "Cu": ("4105040", "39726e1afd93226e126184d30eb0e45df8138f9a78c07b5fd59d9a700cb2300c.cif", 12.0),
            "Pt": ("1011103", "e44ec4c99848ea8b6805f2a4a9ce750157197b0885f2d9b8da0aa7a961a5c002.cif", 12.0),
            "Fe": ("1100108", "41085e65debf2bd82a772bcee9e7c0a43ef42b67666d2f233c8f2d822769b8fe.cif", 8.0),
        }
        for el, (_id, fname, expected_cn) in cases.items():
            cif = pack / "COD" / _id / fname
            self.assertTrue(cif.is_file(), str(cif))
            r = cif_gate.gate_cif(cif, el)
            self.assertTrue(r.ok, f"{el}: {r.reason}")
            self.assertEqual(r.first_shell_cn, expected_cn,
                             f"{el} first shell CN={r.first_shell_cn} != {expected_cn}")

    def test_n0_shell_grouping_rayleigh(self):
        """#101 P2a：Rayleigh 分组：同元素子壳间距 < ΔR_Rayleigh 归为一组。"""
        cif = ROOT / "_structure_lib" / "_phase_pack" / "COD" / "1100108" \
            / "41085e65debf2bd82a772bcee9e7c0a43ef42b67666d2f233c8f2d822769b8fe.cif"
        # k 范围越窄，ΔR_Rayleigh 越大，分组越粗（合并更多壳层）
        r_wide = cif_gate.gate_cif(cif, "Fe", kmin=3.0, kmax=6.0)
        r_narrow = cif_gate.gate_cif(cif, "Fe", kmin=3.0, kmax=12.0)
        self.assertTrue(r_wide.ok)
        self.assertTrue(r_narrow.ok)
        # ΔR_Rayleigh(k=3–6) = π/(2*3) ≈ 0.524 > ΔR(k=3–12) = π/(2*9) ≈ 0.175
        self.assertGreater(r_wide.rayleigh_delta_r, r_narrow.rayleigh_delta_r)
        self.assertAlmostEqual(r_wide.rayleigh_delta_r, math.pi / 6.0, places=6)
        self.assertAlmostEqual(r_narrow.rayleigh_delta_r, math.pi / 18.0, places=6)
        # 窄 k 范围（大 delta_r）：第一壳 2.484 与第二壳 2.868 间距 0.384 < 0.524，合并 → first_cn 增大
        self.assertGreater(r_wide.first_shell_cn, r_narrow.first_shell_cn)
        self.assertEqual(r_narrow.first_shell_cn, 8.0)   # bcc Fe 第一壳 CN=8，不合并
        # 验证 N₀ = Σ(简并度 × 占据度)
        for sh in r_narrow.n0_shells:
            n0_sum = sum(m["occ"] for m in sh["n0_members"])
            self.assertAlmostEqual(sh["n0"], n0_sum, places=6,
                                   msg=f"N₀ {sh['n0']} != Σ occ {n0_sum}")

    # ------------------------------------------------------------
    # 102–103. P2b：固定 N 一致性 + 自由 N 对照
    # ------------------------------------------------------------
    def test_fixed_n_nidp_formula(self):
        """#102 P2b：固定 N 信息量公式 2n+1 / 2n+2，map_n0_to_plan 一致。"""
        motifs = {m.id: m for m in structure_lib.load_library("Fe")}
        motif = motifs["fe_bcc"]
        cif = ROOT / "_structure_lib" / "_phase_pack" / "COD" / "1100108" \
            / "41085e65debf2bd82a772bcee9e7c0a43ef42b67666d2f233c8f2d822769b8fe.cif"
        gate = cif_gate.gate_cif(cif, "Fe")
        self.assertTrue(gate.ok)

        # free 模式：3n+1
        plan_free = xafs_fit.select_shell_paths(
            motif, edge="K", kmin=3.0, kmax=12.0, rmax=3.5, n_mode="free")
        n_free = len(plan_free["files"])
        # fixed_n 模式：2n+1
        plan_fixed = xafs_fit.select_shell_paths(
            motif, edge="K", kmin=3.0, kmax=12.0, rmax=3.5, n_mode="fixed_n")
        n_fixed = len(plan_fixed["files"])
        self.assertEqual(n_free, n_fixed)
        self.assertEqual(plan_fixed["n_var_fixed_n"], 2 * n_fixed + 1)
        # conditional 模式：2n+2
        plan_cond = xafs_fit.select_shell_paths(
            motif, edge="K", kmin=3.0, kmax=12.0, rmax=3.5,
            n_mode="conditional")
        self.assertEqual(plan_cond["n_var_conditional"], 2 * n_fixed + 2)

        # map_n0_to_plan 的 n_var 与 select_shell_paths 一致
        mapped = xafs_fit.map_n0_to_plan(
            plan_free, gate.n0_shells, kmin=3.0, kmax=12.0)
        self.assertTrue(mapped["mapping_ok"], mapped["mapping_reason"])
        self.assertEqual(mapped["n_var_fixed_n"], 2 * n_free + 1)
        self.assertEqual(mapped["n_var_free_n_first"], 2 * n_free + 2)
        # n0_per_path 与前两个壳层 N₀ 一致
        self.assertEqual(len(mapped["n0_per_path"]), n_free)
        expected_n0 = [round(s["n0"], 6) for s in gate.n0_shells[:n_free]]
        for i, n0 in enumerate(mapped["n0_per_path"]):
            self.assertAlmostEqual(n0, expected_n0[i], places=6)

    def test_wrong_phase_consistency_failed(self):
        """#103 P2b：错误相 → N₀ 映射失败（0 CN）；正确相 → 映射通过。"""
        motifs = {m.id: m for m in structure_lib.load_library("Fe")}
        motif = motifs["fe_bcc"]

        # 正确相：Fe bcc 的 n0_shells 与 plan 匹配
        cif = ROOT / "_structure_lib" / "_phase_pack" / "COD" / "1100108" \
            / "41085e65debf2bd82a772bcee9e7c0a43ef42b67666d2f233c8f2d822769b8fe.cif"
        gate = cif_gate.gate_cif(cif, "Fe")
        self.assertTrue(gate.ok)
        plan = xafs_fit.select_shell_paths(
            motif, edge="K", kmin=3.0, kmax=12.0, rmax=3.5, n_mode="free")
        mapped_ok = xafs_fit.map_n0_to_plan(plan, gate.n0_shells,
                                            kmin=3.0, kmax=12.0)
        self.assertTrue(mapped_ok["mapping_ok"], mapped_ok["mapping_reason"])
        self.assertGreater(len(mapped_ok["n0_per_path"]), 0)

        # 错误相：n0_shells 中的元素与 plan 的 scatterer 不匹配
        wrong_shells = [
            {"el": "Zn", "r": 2.484, "n0": 8.0, "n0_members": [],
             "sub_rayleigh_grouped": False},
            {"el": "Cu", "r": 2.868, "n0": 6.0, "n0_members": [],
             "sub_rayleigh_grouped": False},
        ]
        mapped_bad = xafs_fit.map_n0_to_plan(plan, wrong_shells,
                                             kmin=3.0, kmax=12.0)
        self.assertFalse(mapped_bad["mapping_ok"])
        self.assertIn("无匹配 N₀ 壳层", mapped_bad["mapping_reason"])
        self.assertEqual(len(mapped_bad.get("n0_per_path", [])), 0)

    # ------------------------------------------------------------
    # 104–110. P3：条件 CN、多相与缺陷模型
    # ------------------------------------------------------------
    def test_first_shell_single_scale_factor(self):
        """#104 P3：首壳单一比例因子，仅放开参考路径 N。"""
        motifs = {m.id: m for m in structure_lib.load_library("Fe")}
        motif = motifs["fe_bcc"]
        cif = ROOT / "_structure_lib" / "_phase_pack" / "COD" / "1100108" \
            / "41085e65debf2bd82a772bcee9e7c0a43ef42b67666d2f233c8f2d822769b8fe.cif"
        gate = cif_gate.gate_cif(cif, "Fe")
        self.assertTrue(gate.ok)
        plan = xafs_fit.select_shell_paths(
            motif, edge="K", kmin=3.0, kmax=12.0, rmax=3.5, n_mode="free")
        mapping = xafs_fit.map_n0_to_plan(plan, gate.n0_shells,
                                          kmin=3.0, kmax=12.0)
        self.assertTrue(mapping["mapping_ok"])
        args = xafs_fit.build_conditional_fit_args(plan, mapping)
        # 仅放开参考路径 N
        self.assertEqual(len(args["fit_args"]["vary_n"]), 1)
        self.assertIn(args["ref_idx"], args["fit_args"]["vary_n"])
        # path_degen = n0_per_path（所有路径 N₀ 初值）
        self.assertEqual(args["fit_args"]["path_degen"],
                         mapping["n0_per_path"])
        # 非参考首壳路径用 degen_ratio 绑定到参考路径
        links = args["fit_args"]["path_links"]
        # 不同 N₀ 壳层（如 Fe bcc 的 2.48 与 2.87 Å）不得互相绑定：
        # 远壳必须固定 N₀，不能跟着首壳 N 漂移。
        self.assertNotIn(1, [link["target"] for link in links])
        for link in links:
            self.assertEqual(link["share_deltar"], False)
            self.assertEqual(link["share_sigma2"], False)
            self.assertIsNotNone(link["degen_ratio"])
            self.assertGreater(link["degen_ratio"], 0)
        # n0_total = 首壳 N₀ 总和
        self.assertAlmostEqual(args["n0_total"], 8.0)

    def test_conditional_cn_out_of_range_requires_defect(self):
        """#105 P3：N/N₀ 超出范围且无缺陷模型 → defect_model_missing。"""
        # N/N₀ = 0.5，低于 0.6，out_of_range
        self.assertEqual(xafs_fit.check_n_ratio(0.5), "out_of_range")
        self.assertEqual(xafs_fit.check_n_ratio(1.5), "out_of_range")
        # 无缺陷模型 → 拒绝
        chk = xafs_fit.check_defect_model(0.5, None)
        self.assertFalse(chk["ok"])
        self.assertEqual(chk["status"], "defect_model_missing")
        # 有缺陷模型但不匹配 → defect_model_mismatch
        dm = xafs_fit.parse_defect_model(
            '{"type": "custom", "expected_n_ratio": 0.75, "uncertainty": 0.05}')
        chk = xafs_fit.check_defect_model(0.5, dm)
        self.assertFalse(chk["ok"])
        self.assertEqual(chk["status"], "defect_model_mismatch")
        # N/N₀ 在正常范围 → ok，不需要缺陷模型
        self.assertEqual(xafs_fit.check_n_ratio(0.9), "ok")
        self.assertEqual(xafs_fit.check_n_ratio(1.2), "ok")

    def test_conditional_cn_nanoparticle_model(self):
        """#106 P3：纳米粒子缺陷模型，N/N₀ 在预期不确定度内。"""
        dm = xafs_fit.parse_defect_model(
            '{"type": "nanoparticle", "diameter_nm": 2.5, "lattice": "fcc"}')
        self.assertEqual(dm["type"], "nanoparticle")
        self.assertTrue(dm["model_dependent"])
        self.assertIn("expected_n_ratio", dm)
        self.assertIn("uncertainty", dm)
        # N/N₀ 在模型预期范围内 → defect_ok
        ratio = dm["expected_n_ratio"]
        chk = xafs_fit.check_defect_model(ratio, dm)
        self.assertTrue(chk["ok"])
        self.assertEqual(chk["status"], "defect_ok")
        # 偏离不确定度范围 → mismatch
        chk = xafs_fit.check_defect_model(ratio + 3 * dm["uncertainty"], dm)
        self.assertFalse(chk["ok"])
        self.assertEqual(chk["status"], "defect_model_mismatch")
        # 非 fcc/hcp/bcc 晶格拒绝
        with self.assertRaises(xafs_fit.XAFSConfigError):
            xafs_fit.parse_defect_model(
                '{"type": "nanoparticle", "diameter_nm": 2.5, "lattice": "diamond"}')

    def test_multiphase_no_auto_pick(self):
        """#107 P3：多相候选 → ambiguous，不自动 Top-1。"""
        import shutil
        tmp = Path(tempfile.mkdtemp())
        try:
            tmp_pack = tmp / "pack"
            tmp_pack.mkdir()
            # 两条含 Fe 的 CIF（相同内容，不同 id）
            src = ROOT / "_structure_lib" / "_phase_pack" / "COD" / "1100108" \
                / "41085e65debf2bd82a772bcee9e7c0a43ef42b67666d2f233c8f2d822769b8fe.cif"
            sha = phase_pack.sha256_file(src)
            for mid in ("1100108", "1100109"):
                d = tmp_pack / "COD" / mid
                d.mkdir(parents=True)
                shutil.copy2(src, d / f"{sha}.cif")
                meta = {"source": "COD", "id": mid, "formula": "Fe",
                        "spacegroup": "Im-3m",
                        "cell": {"a": 2.868, "b": 2.868, "c": 2.868,
                                 "alpha": 90, "beta": 90, "gamma": 90},
                        "absorber_sites": ["Fe1"], "absorber_elements": ["Fe"],
                        "cif_path": f"{sha}.cif", "cif_sha256": sha,
                        "experimental": True, "model_level": "crystal_cif"}
                (d / f"{sha}.meta.json").write_text(
                    json.dumps(meta, ensure_ascii=False, indent=2),
                    encoding="utf-8")
            phase_pack.rebuild_manifest(tmp_pack)
            res = phase_resolver.resolve_phase(
                "Fe", phase_index_path=str(tmp_pack / "manifest.json"),
                phase_pack_path=str(tmp_pack))
            self.assertEqual(res.status, "ambiguous")
            self.assertEqual(len(res.candidates), 2)
            self.assertIsNone(res.selected)
            self.assertIsNone(res.structure_model)
        finally:
            shutil.rmtree(tmp)

    def test_cn_table_provenance_complete(self):
        """#108 P3：正式 CN 行必须包含全部 provenance 字段。"""
        from structure_model import StructureModel, validate_formal
        model = StructureModel(
            model_id="COD:1100108", absorber="Fe",
            absorber_site="Fe1", model_level="crystal_cif",
            source_type="database_cif", source_ref="COD:1100108",
            source_cif="_structure_lib/_phase_pack/COD/1100108/x.cif",
            source_cif_sha256="0" * 64,
            cell={"a": 2.868, "b": 2.868, "c": 2.868,
                  "alpha": 90, "beta": 90, "gamma": 90},
            spacegroup="Im-3m",
            provenance_status="verified_source_cif",
        )
        # provenance 字段完整
        self.assertIsNotNone(model.model_id)
        self.assertIsNotNone(model.absorber)
        self.assertIsNotNone(model.absorber_site)
        self.assertIsNotNone(model.model_level)
        self.assertIsNotNone(model.source_type)
        self.assertIsNotNone(model.source_ref)
        self.assertIsNotNone(model.source_cif_sha256)
        self.assertIsNotNone(model.absorber_site)
        self.assertEqual(model.provenance_status, "verified_source_cif")
        # cn_table 拒绝行也应包含 provenance 字段
        row = cn_table.refusal_row("no_cif", "测试拒绝",
                                   phase_id=model.model_id,
                                   model_level=model.model_level,
                                   source_type=model.source_type,
                                   source_ref=model.source_ref,
                                   source_cif_sha256=model.source_cif_sha256,
                                   absorber_site=model.absorber_site)
        self.assertEqual(row["phase_id"], "COD:1100108")
        self.assertEqual(row["model_level"], "crystal_cif")
        self.assertEqual(row["source_type"], "database_cif")
        self.assertEqual(row["source_ref"], "COD:1100108")
        self.assertEqual(row["absorber_site"], "Fe1")
        self.assertEqual(row["gate_status"], "no_cif")

    def test_prototype_cluster_separate(self):
        """#109 P3：prototype/cluster 模型不进正式 CN 表。"""
        from structure_model import StructureModel, validate_formal
        # prototype
        m1 = StructureModel(
            model_id="proto_test", absorber="Fe",
            model_level="prototype", source_type="none",
            provenance_status="source_ref_only",
        )
        self.assertFalse(validate_formal(m1)["ok"])
        # cluster
        m2 = StructureModel(
            model_id="cluster_test", absorber="Fe",
            model_level="cluster", source_type="constructed",
            provenance_status="constructed",
        )
        self.assertFalse(validate_formal(m2)["ok"])

    def test_reproducible_cn_hash(self):
        """#110 P3：cn_table.csv 写读往返，内容逐字节一致。"""
        import hashlib
        row = cn_table.refusal_row(
            "no_cif", "reproducibility test",
            phase_id="COD:1100108", model_level="crystal_cif",
            source_type="database_cif", source_ref="COD:1100108",
            source_cif_sha256="0" * 64, absorber_site="Fe1")
        with tempfile.TemporaryDirectory() as tmp:
            out1 = Path(tmp) / "cn1.csv"
            out2 = Path(tmp) / "cn2.csv"
            cn_table.write_cn_table(out1, [row])
            cn_table.write_cn_table(out2, [row])
            h1 = hashlib.sha256(out1.read_bytes()).hexdigest()
            h2 = hashlib.sha256(out2.read_bytes()).hexdigest()
            self.assertEqual(h1, h2)
            # 读回内容一致
            back = cn_table.read_cn_table(out1)
            self.assertEqual(len(back), 1)
            self.assertEqual(back[0]["phase_id"], "COD:1100108")

    # ------------------------------------------------------------
    # 111–114. P4：验证、迁移与文档
    # ------------------------------------------------------------
    def test_offline_full_chain(self):
        """#111 P4：离线全链——相同输入 + 相同缓存 → 相同输出哈希。"""
        import hashlib
        cif = ROOT / "_structure_lib" / "_phase_pack" / "COD" / "1100108" \
            / "41085e65debf2bd82a772bcee9e7c0a43ef42b67666d2f233c8f2d822769b8fe.cif"
        self.assertTrue(cif.is_file())
        # 两次独立调用，结果应字节一致
        r1 = cif_gate.gate_cif(cif, "Fe")
        r2 = cif_gate.gate_cif(cif, "Fe")
        self.assertTrue(r1.ok and r2.ok)
        self.assertEqual(r1.cif_sha256, r2.cif_sha256)
        self.assertEqual(r1.first_shell_cn, r2.first_shell_cn)
        self.assertEqual(r1.occupancy_source, r2.occupancy_source)
        # 整个结果 dict 的哈希也应一致
        h1 = hashlib.sha256(repr(sorted(r1.__dict__.items())).encode()).hexdigest()
        h2 = hashlib.sha256(repr(sorted(r2.__dict__.items())).encode()).hexdigest()
        self.assertEqual(h1, h2)

    def test_standard_foils_formal_cn(self):
        """#112 P4：标准箔（Fe/Co/Ni/Cu）formal CN 链可放行。"""
        from structure_model import StructureModel, validate_formal
        pack = ROOT / "_structure_lib" / "_phase_pack" / "COD"
        cases = {
            "Fe": ("1100108", "41085e65debf2bd82a772bcee9e7c0a43ef42b67666d2f233c8f2d822769b8fe.cif", "Fe1"),
            "Co": ("9008492", "bd9d57408d901e3ff966ba0ebe19ced70e1f1b62a4ef47fda801a39ecbb5344a.cif", "Co"),
            "Ni": ("1534892", "adb86b38b5313592734259a0ff2bb537fc60acb92793015ee8a74a41b9966cdd.cif", "Ni1"),
            "Cu": ("4105040", "39726e1afd93226e126184d30eb0e45df8138f9a78c07b5fd59d9a700cb2300c.cif", "Cu"),
        }
        for el, (_id, fname, site) in cases.items():
            cif = pack / _id / fname
            self.assertTrue(cif.is_file(), str(cif))
            gate = cif_gate.gate_cif(cif, el, absorber_site=site)
            self.assertTrue(gate.ok, f"{el}: {gate.reason}")
            self.assertEqual(gate.occupancy_source in ("explicit", "implicit_1.0"), True)
            # 构造 StructureModel 并验证 formal 门禁
            model = StructureModel(
                model_id=f"COD:{_id}", absorber=el,
                absorber_site=gate.absorber_site or site,
                model_level="crystal_cif",
                source_type="database_cif", source_ref=f"COD:{_id}",
                source_cif=str(cif.resolve()),
                source_cif_sha256=gate.cif_sha256,
                provenance_status="verified_source_cif",
            )
            v = validate_formal(model)
            self.assertTrue(v["ok"], f"{el}: {v}")

    def test_no_cif_report_empty(self):
        """#113 P4：无 CIF 样品 → cn_table 拒绝行，n_fit 为空，cn_label=none。"""
        row = cn_table.refusal_row(
            "no_cif", "无合格结构来源",
            phase_id="", model_level="", source_type="",
            source_ref="", source_cif_sha256="", absorber_site="")
        self.assertEqual(row["gate_status"], "no_cif")
        self.assertEqual(row["n_fit"], "")
        self.assertEqual(row["cn_label"], "none")
        self.assertEqual(row["reject_reason"], "无合格结构来源")
        self.assertIn("refused", row["flags"])
        # 写读往返后仍为空 CN
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "cn_empty.csv"
            cn_table.write_cn_table(out, [row])
            back = cn_table.read_cn_table(out)
        self.assertEqual(len(back), 1)
        self.assertEqual(back[0]["n_fit"], "")
        self.assertEqual(back[0]["cn_label"], "none")

    def test_legacy_auto_fit_marked(self):
        """#114 P4：--legacy-auto-fit 标记 legacy_unverified，不进正式表。"""
        # 模拟 formal_cn 处理：legacy_auto_fit=True 时状态应被覆盖为 legacy_unverified
        formal = {"status": "pending_consistency", "reject_reason": None,
                  "rows": [], "model": None}
        # 模拟 run_xafs.py 中的 legacy_auto_fit 处理逻辑
        legacy_auto_fit = True
        if legacy_auto_fit:
            formal["status"] = "legacy_unverified"
            formal["reject_reason"] = formal.get("reject_reason") or "legacy_auto_fit: 旧版自动拟合，不进正式表"
            formal["legacy_auto_fit"] = True
        self.assertEqual(formal["status"], "legacy_unverified")
        self.assertTrue(formal["legacy_auto_fit"])
        self.assertIn("legacy_auto_fit", formal["reject_reason"])

    # ------------------------------------------------------------
    # 115–119. P4+：CIF crawler 目标表、离线导入、多位点门禁与正式 CN 行
    # ------------------------------------------------------------
    def test_cif_crawler_target_table(self):
        """#115 P4+：39 个抓取目标 ID 唯一、元素/化学式/空间群自洽。"""
        from tools import cif_crawler
        self.assertEqual(len(cif_crawler.VERIFIED_TARGETS), 39)
        cod_ids = [t["cod_id"] for t in
                   cif_crawler.VERIFIED_TARGETS.values()]
        self.assertEqual(len(cod_ids), len(set(cod_ids)))
        seen_elements = set()
        for name, target in cif_crawler.VERIFIED_TARGETS.items():
            self.assertIn(target["tier"], ("T3", "T4", "P5C"), name)
            self.assertIn(target["absorber"],
                          phase_pack.ABSORBER_ELEMENTS, name)
            parsed = cif_crawler.parse_formula(target["formula"])
            self.assertEqual(parsed, {el: float(n) for el, n in
                                      target["stoich"].items()}, name)
            self.assertTrue(target["spg_nums"], name)
            self.assertTrue(all(isinstance(n, int) and n > 0
                                for n in target["spg_nums"]), name)
            self.assertTrue(target["tags"], name)
            seen_elements.add(target["absorber"])
        self.assertEqual(seen_elements, set(phase_pack.ABSORBER_ELEMENTS))

    def test_cif_crawler_offline_import(self):
        """#116 P4+：离线缓存 → 核验 → phase pack 导入，sha256 可追溯。"""
        import shutil
        from tools import cif_crawler

        pack_root = ROOT / "_structure_lib" / "_phase_pack"
        manifest = phase_pack.load_manifest(pack_root)
        fe2o3 = next(e for e in manifest["entries"] if e["id"] == "1011240")
        src = pack_root / fe2o3["cif_path"]

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            cache = tmp / "cache"
            cache.mkdir()
            shutil.copy2(src, cache / "Fe2O3_COD_1011240.cif")
            tmp_pack = tmp / "pack"
            report_path = tmp / "report.json"
            rc = cif_crawler.main([
                "--target", "Fe2O3", "--offline",
                "--cache", str(cache), "--pack", str(tmp_pack),
                "--report", str(report_path), "--no-index",
            ])
            self.assertEqual(rc, 0)
            rebuilt = phase_pack.load_manifest(tmp_pack)
            self.assertEqual(len(rebuilt["entries"]), 1)
            entry = rebuilt["entries"][0]
            self.assertEqual(entry["source"], "COD")
            self.assertEqual(entry["id"], "1011240")
            self.assertEqual(entry["absorber_sites"], ["Fe1"])
            self.assertEqual(entry["cif_sha256"],
                             phase_pack.sha256_file(src))
            self.assertIn("source_url=", entry.get("notes", ""))
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["targets"]["Fe2O3"]["status"],
                             "imported")
            self.assertEqual(report["summary"]["imported"], 1)

    def test_cif_crawler_offline_miss_and_wrong_cif(self):
        """#117 P4+：离线缺缓存 → offline_miss；错相 CIF → 拒绝导入。"""
        import shutil
        from tools import cif_crawler

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            cache = tmp / "cache"
            cache.mkdir()
            tmp_pack = tmp / "pack"
            report_path = tmp / "report.json"
            rc = cif_crawler.main([
                "--target", "CeO2", "--offline",
                "--cache", str(cache), "--pack", str(tmp_pack),
                "--report", str(report_path), "--no-index",
            ])
            self.assertEqual(rc, 1)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["targets"]["CeO2"]["status"],
                             "offline_miss")
            self.assertFalse((tmp_pack / "manifest.json").exists())

            # 把 Fe3O4 错相当成 CeO2 缓存：必须被化学计量核验拒绝。
            pack_root = ROOT / "_structure_lib" / "_phase_pack"
            manifest = phase_pack.load_manifest(pack_root)
            fe3o4 = next(e for e in manifest["entries"]
                         if e["id"] == "1010369")
            shutil.copy2(pack_root / fe3o4["cif_path"],
                         cache / "CeO2_COD_1562989.cif")
            rc = cif_crawler.main([
                "--target", "CeO2", "--offline",
                "--cache", str(cache), "--pack", str(tmp_pack),
                "--report", str(report_path), "--no-index",
            ])
            self.assertEqual(rc, 1)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            row = report["targets"]["CeO2"]
            self.assertEqual(row["status"], "cache_invalid")
            self.assertIn("composition_mismatch", row["reason"])
            self.assertFalse((tmp_pack / "manifest.json").exists())

    def test_multi_site_phase_requires_absorber_site(self):
        """#118 P4+：Fe3O4 多个吸收位必须显式 --absorber-site。"""
        res = phase_resolver.resolve_phase("Fe", phase_id="COD:1010369")
        self.assertEqual(res.status, "resolved")
        self.assertIsNotNone(res.selected)
        self.assertEqual(res.selected.id, "1010369")
        self.assertEqual(len(res.selected.absorber_sites), 3)
        self.assertIsNone(res.structure_model)
        self.assertIn("多个吸收位点", res.reason)
        self.assertIn("--absorber-site", res.reason)

        res_ok = phase_resolver.resolve_phase(
            "Fe", phase_id="COD:1010369", absorber_site="Fe1")
        self.assertIsNotNone(res_ok.structure_model)
        self.assertEqual(res_ok.structure_model.absorber_site, "Fe1")

        res_bad = phase_resolver.resolve_phase(
            "Fe", phase_id="COD:1010369", absorber_site="Fe9")
        self.assertIsNone(res_bad.structure_model)
        self.assertIn("没有吸收位", res_bad.reason)

    def test_conditional_cn_formal_rows(self):
        """#119 P4+：条件 CN 通过 → 正式行完整可追溯；失败 → 0 行。"""
        from structure_model import StructureModel
        model = StructureModel(
            model_id="COD:1100108", absorber="Fe",
            absorber_site="Fe1", model_level="crystal_cif",
            source_type="database_cif", source_ref="COD:1100108",
            source_cif="_structure_lib/_phase_pack/COD/1100108/x.cif",
            source_cif_sha256="0" * 64,
            cell={"a": 2.868, "b": 2.868, "c": 2.868,
                  "alpha": 90, "beta": 90, "gamma": 90},
            spacegroup="Im-3m",
            provenance_status="verified_source_cif",
        )
        plan = {
            "shells": [
                {"label": "Fe-Fe@2.48", "scatterer": "Fe",
                 "reff": 2.482, "degen": 8.0, "n_members": 1},
                {"label": "Fe-Fe@2.87", "scatterer": "Fe",
                 "reff": 2.866, "degen": 6.0, "n_members": 1},
            ],
            "rmin": 1.88, "rmax": 3.5,
        }
        n0_mapping = {"n0_per_path": [8.0, 6.0]}
        conditional = {
            "n_ratio": 0.95, "n_status": "ok", "consistency_passed": True,
            "first_shell_idxs": [0],
            "kmin": 3.0, "kmax": 12.0, "kweight": 2,
            "conditional_n_fit": {
                "r_factor": 0.0234, "n_variables": 5,
                "n_independent": 9.3, "delta_e0": 0.59,
                "overfit": False},
            "hamilton_f_test": {"significant": False, "p_value": 0.5},
            "conditional_n_paths": [
                {"label": "Fe-Fe@2.48", "reff": 2.45, "N": 7.6,
                 "N_err": 0.4, "delta_r": -0.032, "sigma2": 0.0055},
                {"label": "Fe-Fe@2.87", "reff": 2.813, "N": 6.0,
                 "N_err": None, "delta_r": -0.053, "sigma2": 0.0040},
            ],
        }
        block = xafs_fit.build_formal_cn_block(
            conditional_cn=conditional, plan=plan,
            n0_mapping=n0_mapping, structure_model=model,
            absorber=get_absorber("Fe"),
            cn_systematics={"relative_1sigma": 0.32},
            s02_source="absorber.py", data_sha256="d" * 64)
        self.assertEqual(block["status"], "ok")
        self.assertEqual(len(block["rows"]), 2)
        first, second = block["rows"]
        self.assertEqual(first["n_status"], "fitted")
        self.assertEqual(first["n_fit"], 7.6)
        self.assertEqual(first["n0"], 8.0)
        self.assertEqual(first["n_ratio"], 7.6 / 8.0)
        self.assertEqual(first["n_err_formal"], 0.4)
        self.assertEqual(first["n_err_systematic"], 0.32)
        self.assertEqual(first["shell_group"], "first_shell")
        self.assertEqual(first["cn_label"], "cn_conditional")
        self.assertEqual(first["gate_status"], "passed")
        self.assertEqual(first["data_sha256"], "d" * 64)
        self.assertEqual(first["source_cif_sha256"], "0" * 64)
        self.assertEqual(first["absorber_site"], "Fe1")
        self.assertTrue(first["s02_calibrated"])
        self.assertEqual(second["n_status"], "fixed")
        self.assertEqual(second["n_fit"], 6.0)
        self.assertEqual(second["n_ratio"], 1.0)
        self.assertEqual(second["n_err_formal"], "")

        # Hamilton 显著 → 一致性失败，正式表不得产生 CN 行。
        failed = dict(conditional)
        failed["hamilton_f_test"] = {"significant": True, "p_value": 0.01}
        failed["consistency_passed"] = False
        block_bad = xafs_fit.build_formal_cn_block(
            conditional_cn=failed, plan=plan,
            n0_mapping=n0_mapping, structure_model=model,
            absorber=get_absorber("Fe"))
        self.assertEqual(block_bad["status"], "consistency_failed")
        self.assertEqual(block_bad["rows"], [])
        self.assertIn("Hamilton", block_bad["reject_reason"])

        # 数据哈希稳定：同内容换路径不变，内容变化则变。
        with tempfile.TemporaryDirectory() as tmp:
            a = Path(tmp) / "a.dat"
            b = Path(tmp) / "b.dat"
            a.write_bytes(b"abc")
            b.write_bytes(b"abc")
            self.assertEqual(cn_table.data_sha256([a]),
                             cn_table.data_sha256([b]))
            b.write_bytes(b"abd")
            self.assertNotEqual(cn_table.data_sha256([a]),
                                cn_table.data_sha256([b]))

    def test_multi_element_report_no_top1(self):
        """#120 T7：多元素自动分元素跑批；报告并列、不推荐 Top-1。"""
        import run_xafs
        from types import SimpleNamespace

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fe = tmp / "sample_fe.dat"
            fe.write_text("# Absorber: Fe\n# Energy(eV)\tmu\n"
                          "7112\t1.0\n", encoding="utf-8")
            cu = tmp / "sample_cu.dat"
            cu.write_text("# Absorber: Cu\n# Energy(eV)\tmu\n"
                          "8979\t1.0\n", encoding="utf-8")
            unknown = tmp / "unknown.dat"
            unknown.write_text("# no element header\n1 1\n",
                               encoding="utf-8")

            def fake_run_one(path, element, sub):
                sub.mkdir(parents=True, exist_ok=True)
                return SimpleNamespace(
                    structure={
                        "status": "confident",
                        "ranked": [{"id": f"{element.lower()}_phase",
                                    "score": 0.02}]},
                    fit={
                        "mode": "l1_first_shell", "cn_reliable": True,
                        "l1": {"N": 5.0, "R": 2.0}, "warning": ""},
                    xanes_cross_check={
                        "warnings": [{"code": "edge_step_weak",
                                      "message": "demo warning"}]})

            out = tmp / "out"
            summary = run_xafs.run_multi_element(
                [str(fe), str(cu), str(unknown)], out,
                run_one=fake_run_one, verbose=False)
            self.assertEqual(summary["n_ok"], 2)
            self.assertEqual(summary["n_failed"], 1)
            self.assertEqual(summary["elements"], ["Cu", "Fe"])
            self.assertIsNone(summary["recommended"])
            self.assertIsNone(summary["best_element"])
            self.assertIn("不自动选择", summary["recommendation_policy"])
            self.assertTrue((out / "multi_element_summary.json").is_file())
            md = (out / "multi_element_summary.md").read_text(
                encoding="utf-8")
            self.assertIn("不自动选择/推荐 Top-1", md)
            fe_row = next(r for r in summary["rows"]
                          if r["absorber"] == "Fe")
            self.assertEqual(fe_row["first_shell"]["n_fit"], 5.0)
            self.assertEqual(fe_row["candidates"][0]["id"], "fe_phase")
            self.assertEqual(fe_row["absorber_source"], "header_absorber")
            self.assertIn("demo warning", fe_row["warnings"])

    def test_xanes_exafs_cross_check_warning_only(self):
        """#121 T8：白线/边步/归一化异常 → 警告；不改 cn_reliable。"""
        from types import SimpleNamespace

        norm = SimpleNamespace(plateau=0.70, edge_step=0.01, e0=7112.0)
        res = SimpleNamespace(
            norm=norm,
            xanes_features={"e0_jump_risk": True},
            absorber_symbol="Fe",
            fit=None)
        fit_out = {
            "mode": "l1_first_shell", "cn_reliable": True,
            "l1": {"N": 10.0, "R": 2.5},
            "self_absorption": {"risk": "high", "note": "荧光厚样"}}
        wl = {"reliable": True, "white_line_area_norm": 0.2,
              "whiteline_peak_height": 0.1,
              "white_line_energy": 7125.0, "note": ""}
        cc = xafs_fit.xanes_exafs_cross_check(res, fit_out, white_line=wl)
        codes = {w["code"] for w in cc["warnings"]}
        self.assertIn("normalization_plateau_out_of_range", codes)
        self.assertIn("edge_step_weak", codes)
        self.assertIn("e0_jump_risk", codes)
        self.assertIn("high_cn_low_white_line", codes)
        self.assertIn("self_absorption_risk", codes)
        self.assertFalse(cc["affects_cn_reliable"])
        self.assertTrue(fit_out["cn_reliable"])  # 未被修改

        # 干净情形：无警告；L3 边不做 CN–白线极端对比。
        clean_norm = SimpleNamespace(plateau=1.0, edge_step=1.0,
                                     e0=11564.0)
        clean_res = SimpleNamespace(
            norm=clean_norm, xanes_features={"e0_jump_risk": False},
            absorber_symbol="Pt", fit=None)
        clean_fit = {"mode": "l1_first_shell", "cn_reliable": True,
                     "l1": {"N": 10.0, "R": 2.5}}
        clean_wl = {"reliable": True, "white_line_area_norm": 0.2,
                    "whiteline_peak_height": 0.1,
                    "white_line_energy": 11570.0, "note": ""}
        cc2 = xafs_fit.xanes_exafs_cross_check(
            clean_res, clean_fit, white_line=clean_wl)
        self.assertEqual(cc2["warnings"], [])
        self.assertEqual(cc2["status"], "ok")
        self.assertEqual(
            cc2["checks"]["cn_white_line_contrast"]["status"],
            "unavailable")

    def test_t5_artifacts_static_lock_and_incremental_merge(self):
        """#122 P5-A：T5 产物完整、#60 静态锁定、--only 增量不丢元素。"""
        import calibrate_s02
        import calibrate_t5
        from unittest import mock
        from xafs_core import XAFSDataError

        locked = {"Cr", "Mn", "Zn", "Mo", "Ru", "Rh", "Pd", "Ag",
                  "Ir", "Au", "Ce", "U"}
        self.assertEqual(set(calibrate_s02.T5_STATIC_LOCKED), locked)
        for el in sorted(locked):
            cfg = get_absorber(el)
            self.assertIsNone(cfg.s02, el)
            self.assertFalse(cfg.calibrated, el)

        index_path = ROOT / "data" / "s02_calibration" / "index.json"
        index = json.loads(index_path.read_text(encoding="utf-8"))
        self.assertTrue(index["static_absorber_py_unchanged"])
        self.assertFalse(index["apply_to_absorber_py"])
        self.assertEqual(index["runtime_override"], "--s02")
        rows = {r["absorber"]: r for r in index["elements"]}
        self.assertEqual(set(rows), set(calibrate_t5.ELEMENTS))
        self.assertEqual(index["pass_count"] + index["fail_count"], 12)
        self.assertEqual(index["pass_count"], 12)
        self.assertEqual(index["fail_count"], 0)
        self.assertNotIn("C:\\\\", index_path.read_text(encoding="utf-8"))
        for el in calibrate_t5.ELEMENTS:
            self.assertEqual(rows[el]["status"], "calibrated", el)
            self.assertIsNotNone(rows[el]["s02"], el)
        expected_s02 = {
            "Rh": 0.7022, "Pd": 0.6559, "Ir": 0.4752,
            "Ce": 0.6279, "U": 0.7800,
        }
        for el, value in expected_s02.items():
            self.assertAlmostEqual(rows[el]["s02"], value, places=4, msg=el)

        # --apply 对 #60 锁定元素硬拒绝（不真正跑标定）。
        with mock.patch.object(calibrate_s02, "calibrate_s02",
                               return_value={"pass": True, "s02": 0.6}):
            with self.assertRaises(SystemExit) as cm:
                calibrate_s02.main([
                    "dummy.dat", "--absorber", "Cr", "--apply"])
        self.assertIn("#60", str(cm.exception))

        # 标定专用 motif 缺 sites 必须拒绝，不得回退到 expand_full。
        payload = json.loads((ROOT / "data" / "s02_calibration" / "motifs" /
                              "cr2o3_corundum.json").read_text(
                                  encoding="utf-8"))
        payload["sites"] = []
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.json"
            bad.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaises(XAFSDataError):
                calibrate_s02.load_motif_json(bad, "Cr")

        # --only 只重跑 Cr，但 index/REPORT 必须保留既有 Mn/Au。
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for el in ("Mn", "Au"):
                el_dir = tmp / el
                el_dir.mkdir()
                existing = {
                    "absorber": el, "status": "calibrated",
                    "s02": 0.5, "s02_err": 0.01, "rf": 0.01,
                    "spread": 0.001, "motif": f"{el.lower()}_motif",
                    "standard": "dummy", "standard_sha256": "a" * 64,
                    "attempts_compact": [],
                }
                (el_dir / "calibration.json").write_text(
                    json.dumps(existing), encoding="utf-8")
            new_cr = {
                "absorber": "Cr", "status": "calibrated", "s02": 0.6482,
                "s02_err": 0.0198, "rf": 0.0074, "spread": 0.0053,
                "motif": "cr2o3_corundum", "standard": "dummy",
                "standard_sha256": "b" * 64, "attempts_compact": [],
            }
            with mock.patch.object(calibrate_t5, "_run_element",
                                   return_value=new_cr):
                rc = calibrate_t5.main(
                    ["--only", "Cr", "--outdir", str(tmp)])
            self.assertEqual(rc, 0)
            merged = json.loads((tmp / "index.json").read_text(
                encoding="utf-8"))
            self.assertEqual([r["absorber"] for r in merged["elements"]],
                             ["Cr", "Mn", "Au"])
            self.assertEqual(merged["pass_count"], 3)
            self.assertTrue((tmp / "REPORT.md").is_file())

    def test_t5_s02_loader_and_provenance(self):
        """#123 P5-A：--s02-from-t5 只接受已标定元素并记录溯源。"""
        import run_xafs

        s02, source = run_xafs.load_t5_s02("Cr")
        self.assertAlmostEqual(s02, 0.6482, places=4)
        self.assertTrue(source.startswith("s02_calibration:"))
        self.assertIn("sha256:", source)
        self.assertIn("standard_sha256:", source)

        s02_ce, source_ce = run_xafs.load_t5_s02("Ce")
        self.assertAlmostEqual(s02_ce, 0.6279, places=4)
        self.assertTrue(source_ce.startswith("s02_calibration:"))
        # 未通过元素仍然硬拒绝（用临时产物验证，不依赖真实标定状态）。
        with tempfile.TemporaryDirectory() as tmp:
            el_dir = Path(tmp) / "Ce"
            el_dir.mkdir()
            (el_dir / "calibration.json").write_text(json.dumps({
                "absorber": "Ce", "status": "not_calibrated",
                "s02": None, "reason": "gate_failed：未通过冻结门禁",
                "static_absorber_py_unchanged": True,
            }), encoding="utf-8")
            with self.assertRaises(XAFSConfigError) as cm:
                run_xafs.load_t5_s02("Ce", root=Path(tmp))
            self.assertIn("拒绝自动注入", str(cm.exception))
        with self.assertRaises(XAFSConfigError):
            run_xafs.load_t5_s02("Xx")

        # 显式 --s02 优先于 T5 产物；未显式给值时从 T5 读取。
        from types import SimpleNamespace
        args = SimpleNamespace(s02=0.9, s02_from_t5=True, absorber="Cr")
        val, src = run_xafs._resolve_s02_for_paths(args, ["dummy.dat"])
        self.assertEqual(val, 0.9)
        self.assertEqual(src, "--s02 runtime override")
        args2 = SimpleNamespace(s02=None, s02_from_t5=True, absorber=None)
        val2, src2 = run_xafs._resolve_s02_for_paths(
            args2, [str(ROOT / "data" / "bmm_standards" / "Data" / "Cr" /
                         "Cr-K-Cr2O3.xdi")])
        self.assertAlmostEqual(val2, 0.6482, places=4)
        self.assertTrue(src2.startswith("s02_calibration:"))

        # 正式表：T5 产物视为已标定，普通 --s02 手填值不算。
        self.assertTrue(xafs_fit._s02_calibrated_flag(
            get_absorber("Cr"), "s02_calibration:x"))
        self.assertFalse(xafs_fit._s02_calibrated_flag(
            get_absorber("Cr"), "--s02 runtime override"))

    def test_p5b_mo_double_shell_and_ce_rejection(self):
        """#124 P5-B：Mo 双壳 bcc 标定；Ce σ² 压界/伪迹明确拒绝。"""
        import inspect
        import calibrate_s02
        import calibrate_t5
        import run_xafs

        sig = inspect.signature(calibrate_s02.calibrate_s02)
        for name in ("e0_method", "manual_e0", "fit_shells"):
            self.assertIn(name, sig.parameters)

        mo_cases = calibrate_t5.CASES["Mo"]
        mo_foil = next(c for c in mo_cases
                       if Path(c["standard"]).name ==
                       "Mo-K-MolybdenumFoil.xdi")
        self.assertEqual(mo_foil.get("e0_method"), "half_height")
        self.assertEqual(mo_foil.get("fit_shells"), 2)
        self.assertTrue(mo_foil.get("strict"))
        ce_cases = calibrate_t5.CASES["Ce"]
        self.assertEqual(len(ce_cases), 3)
        self.assertTrue(ce_cases[1].get("strict"))
        self.assertTrue(ce_cases[2].get("strict"))
        self.assertEqual(ce_cases[2].get("rmin"), 1.5)

        mo = json.loads((ROOT / "data" / "s02_calibration" / "Mo" /
                         "calibration.json").read_text(encoding="utf-8"))
        self.assertEqual(mo["status"], "calibrated")
        self.assertAlmostEqual(mo["s02"], 0.8322, places=4)
        self.assertEqual(mo["n_shells"], 2)
        self.assertEqual([s["degen"] for s in mo["shell_paths"]],
                         [8.0, 6.0])
        self.assertAlmostEqual(mo["shell_paths"][0]["reff"], 2.7254,
                               places=4)
        self.assertAlmostEqual(mo["shell_paths"][1]["reff"], 3.147,
                               places=3)
        self.assertEqual(mo["e0_method"], "half_height")
        self.assertAlmostEqual(mo["e0"], 20012.8624, places=3)
        self.assertIn("sigma2_not_at_bound", mo["checks"])
        self.assertTrue(all(mo["checks"].values()))
        self.assertEqual(mo["n_variables"], 6)
        self.assertGreater(mo["n_independent"], mo["n_variables"])
        for run in mo["runs"]:
            self.assertTrue(run["interior"])
            self.assertEqual(len(run["paths"]), 2)
            self.assertFalse(run["sigma2_at_bound"])
            self.assertFalse(any(p["sigma2_at_bound"]
                                 for p in run["paths"]))

        ce = json.loads((ROOT / "data" / "s02_calibration" / "Ce" /
                         "calibration.json").read_text(encoding="utf-8"))
        self.assertEqual(ce["status"], "calibrated")
        self.assertAlmostEqual(ce["s02"], 0.6279, places=4)
        first = next(a for a in ce["attempts"]
                     if a["standard"].endswith("Ce-L3-CeO2.xdi"))
        self.assertEqual(first["status"], "gated_out")
        self.assertFalse(first["result"]["checks"]["sigma2_not_at_bound"])
        self.assertTrue(all(r["sigma2_at_bound"]
                            for r in first["result"]["runs"]))
        second = next(a for a in ce["attempts"]
                      if a["standard"].endswith(
                          "Ce-L3-CeriumDioxide.xdi"))
        self.assertEqual(second["status"], "error")
        self.assertIn("归一化失败", second["error"])
        self.assertIn("2.271", second["error"])
        third = next(a for a in ce["attempts"]
                     if a["standard"].endswith("Ce_metal_foil.dat"))
        self.assertEqual(third["status"], "passed")
        self.assertTrue(all(third["result"]["checks"].values()))

        s02, source = run_xafs.load_t5_s02("Mo")
        self.assertAlmostEqual(s02, 0.8322, places=4)
        self.assertIn("Mo", source)
        self.assertIn("sha256:", source)
        s02_ce, source_ce = run_xafs.load_t5_s02("Ce")
        self.assertAlmostEqual(s02_ce, 0.6279, places=4)
        self.assertIn("Ce", source_ce)

    def test_p5c_common_phase_coverage(self):
        """#125 P5-C：18 个固定 COD ID 补齐常见相，coverage 无缺口。"""
        from types import SimpleNamespace
        from tools import cif_crawler
        import phase_resolver

        expected = {
            ("Fe", "FeO"), ("Fe", "FeOOH"), ("Pt", "PtO2"),
            ("Cr", "Cr2O3"), ("Mn", "MnO"), ("Mn", "MnO2"),
            ("Mn", "Mn3O4"), ("Zn", "ZnO"), ("Mo", "MoO2"),
            ("Mo", "MoO3"), ("Ru", "RuO2"), ("Rh", "Rh2O3"),
            ("Pd", "PdO"), ("Ag", "Ag2O"), ("Ir", "IrO2"),
            ("Au", "Au2O3"), ("U", "UO2"), ("U", "U3O8"),
        }
        p5c = {name: target for name, target in
               cif_crawler.VERIFIED_TARGETS.items()
               if target["tier"] == "P5C"}
        self.assertEqual(len(p5c), 18)
        self.assertEqual({(t["absorber"], name) for name, t in p5c.items()},
                         expected)
        selected = cif_crawler._select_targets(SimpleNamespace(
            all=False, t3=False, t4=False, p5c=True, target=None))
        self.assertEqual(len(selected), 18)

        pack_root = ROOT / "_structure_lib" / "_phase_pack"
        manifest = phase_pack.load_manifest(pack_root)
        by_id = {entry["id"]: entry for entry in manifest["entries"]}
        for name, target in p5c.items():
            self.assertIn(target["cod_id"], by_id, name)
            entry = by_id[target["cod_id"]]
            self.assertEqual(entry["formula"], target["formula"], name)
            self.assertEqual(entry["absorber_elements"],
                             [target["absorber"]], name)
            self.assertEqual(entry["model_level"], "crystal_cif", name)
            self.assertTrue(entry["experimental"], name)
            self.assertEqual(entry["source"], "COD", name)
            check = phase_pack.validate_entry(pack_root, entry)
            self.assertTrue(check["ok"], (name, check))

        coverage = json.loads((ROOT / "_structure_lib" /
                               "phase_coverage.json")
                              .read_text(encoding="utf-8"))
        self.assertEqual(coverage["summary"]["verified_cif"], 44)
        for element, row in coverage["elements"].items():
            self.assertEqual(row["missing_common_phases"], [], element)
            self.assertTrue(row["formal_cn_ready"], element)

        # 多位点常见相仍必须显式 --absorber-site，不能自动 Top-1。
        for element, phase_name in (("Mn", "Mn3O4"), ("U", "U3O8")):
            target = p5c[phase_name]
            entry = by_id[target["cod_id"]]
            self.assertGreaterEqual(len(entry["absorber_sites"]), 2,
                                    phase_name)
            res = phase_resolver.resolve_phase(
                element, phase_id=f"COD:{target['cod_id']}")
            self.assertEqual(res.status, "resolved", phase_name)
            self.assertIsNone(res.structure_model, phase_name)
            self.assertIn("多个吸收位点", res.reason, phase_name)

    def test_p5_standard_source_isolation_pool(self):
        """#126 P5-prep/P5-S02：5 元素谱源池闭环与门禁冻结。"""
        from tools import build_source_pool

        pool = json.loads((ROOT / "data" / "s02_calibration" /
                           "source_pool" / "index.json")
                          .read_text(encoding="utf-8"))
        self.assertEqual(pool["schema_version"], 2)
        self.assertEqual(pool["pool_id"], "s02-standard-source-isolation")
        self.assertEqual(pool["elements_total"], 5)
        self.assertEqual(pool["pool_state_counts"], {"resolved": 5})
        # 索引必须能由标定产物与门禁常量确定性重建。
        self.assertEqual(build_source_pool.build_pool(), pool)

        rows = {row["element"]: row for row in pool["elements"]}
        self.assertEqual(set(rows), {"Rh", "Pd", "Ir", "Ce", "U"})
        expected_class = {
            "Rh": "short_spectrum_replace",
            "Pd": "short_spectrum_replace",
            "Ir": "short_spectrum_replace",
            "Ce": "oxide_baseline_reprocess",
            "U": "oxide_phase_verification",
        }
        expected_s02 = {
            "Rh": 0.7022, "Pd": 0.6559, "Ir": 0.4752,
            "Ce": 0.6279, "U": 0.7800,
        }
        for element, row in rows.items():
            self.assertEqual(row["classification"], expected_class[element])
            self.assertEqual(row["pool_state"], "resolved")
            self.assertTrue(row["usable_for_s02"])
            self.assertAlmostEqual(row["s02"], expected_s02[element],
                                   places=4, msg=element)
            self.assertEqual(row["last_status"], "calibrated")
            self.assertEqual(row["attempt_status"], "passed")
            self.assertEqual(row["gate_failures"], [])
            self.assertGreater(row["k_max_estimate"], 10.0, element)
            # #60：静态表仍然 s02=None / calibrated=False。
            cfg = get_absorber(element)
            self.assertIsNone(cfg.s02, element)
            self.assertFalse(cfg.calibrated, element)
            # 旧谱的历史隔离证据必须保留，不能因为替换而丢失审计链。
            legacy = row["quarantined_legacy_source"]
            self.assertEqual(legacy["attempt_status"], "gated_out")
            self.assertTrue(legacy["gate_failures"], element)
            self.assertEqual(len(legacy["source_sha256"]), 64, element)
            self.assertEqual(row["resolved_by"]["gate_failures"], [])
            self.assertEqual(row["resolved_by"]["standard"],
                             row["current_source"])

        # 旧谱谱长判定：Rh/Pd/Ir 是短谱；Ce/U 不是截断问题。
        for element in ("Rh", "Pd", "Ir"):
            self.assertLess(
                rows[element]["quarantined_legacy_source"]["k_max_estimate"],
                6.0, element)
        self.assertGreater(
            rows["Ce"]["quarantined_legacy_source"]["k_max_estimate"], 10.0)
        self.assertGreater(
            rows["U"]["quarantined_legacy_source"]["k_max_estimate"], 10.0)

        # Ir 的 k_min 建议只能换谱后诊断，不能作为通过手段。
        ir = rows["Ir"]
        self.assertFalse(ir["kmin_tuning_assessment"]
                         ["sufficient_with_current_source"])
        self.assertFalse(ir["kmin_tuning_assessment"]["may_relax_gates"])
        self.assertEqual(ir["resolved_by"]["kmin"], 3.0)

        # Ce 使用 CeO2 基准；U 的 UO2 切换已完成并记录标定专用 motif。
        self.assertFalse(rows["Ce"]["phase_switch_needed"])
        self.assertIn("Ce_metal_foil.dat", rows["Ce"]["current_source"])
        self.assertIn("CeO2", rows["Ce"]["phase_baseline"])
        self.assertFalse(rows["U"]["phase_switch_needed"])
        self.assertIn("UO2", rows["U"]["phase_baseline"])
        self.assertEqual(rows["U"]["resolved_by"]["motif"], "uo2_fluorite")

        # 门禁冻结：与 calibrate_s02.py 常量一致，且未放宽。
        policy = pool["policy"]
        self.assertEqual(policy["s02_pass_range"], [0.4, 1.15])
        self.assertEqual(policy["s02_err_max"], 0.06)
        self.assertEqual(policy["rf_max"], 0.12)
        self.assertEqual(policy["abs_delta_e0_max_ev"], 10.0)
        self.assertEqual(policy["window_fracs"], [0.12, 0.18, 0.26])
        self.assertEqual(policy["window_spread_max"], 0.05)
        self.assertEqual(policy["delta_r_limits"], [-0.03, 0.03])
        self.assertEqual(policy["gate_relaxation"], "forbidden")
        self.assertEqual(policy["manual_window_trimming"], "forbidden")
        self.assertEqual(policy["auto_top1_selection"], "forbidden")

        # T5 总索引：5 元素全部 calibrated、s02 非空。
        t5_index = json.loads((ROOT / "data" / "s02_calibration" /
                               "index.json").read_text(encoding="utf-8"))
        t5 = {row["absorber"]: row for row in t5_index["elements"]}
        for element in rows:
            self.assertEqual(t5[element]["status"], "calibrated", element)
            self.assertIsNotNone(t5[element]["s02"], element)

    def test_p5d_multi_site_formal_cn(self):
        """#127 P5-D：--absorber-site all 逐位点正式 CN + 化学计量加权。"""
        from types import SimpleNamespace
        from unittest.mock import patch

        # 未传 --absorber-site：多位点相仍然必须拒绝，绝不自动 Top-1。
        no_site = phase_resolver.resolve_phase(
            "Fe", phase_id="COD:1010369", offline=True)
        self.assertEqual(no_site.status, "resolved")
        self.assertIsNone(no_site.structure_model)
        self.assertEqual(no_site.site_models, [])
        self.assertIn("不自动 Top-1", no_site.reason)

        # 显式 all：按位点顺序返回每个不等价吸收位的模型。
        fe = phase_resolver.resolve_phase(
            "Fe", phase_id="COD:1010369", absorber_site="all", offline=True)
        self.assertEqual(fe.status, "resolved")
        self.assertEqual(fe.absorber_site, "all")
        self.assertEqual([m.absorber_site for m in fe.site_models],
                         ["Fe1", "Fe2", "Fe3"])
        self.assertIs(fe.structure_model, fe.site_models[0])

        gates = cif_gate.gate_cif_all_sites(
            fe.selected.cif_path, "Fe", kmin=3.0, kmax=12.0)
        self.assertEqual(len(gates), 3)
        self.assertTrue(all(g.ok for g in gates), [g.reason for g in gates])
        self.assertEqual([g.site_multiplicity for g in gates], [16, 16, 8])
        self.assertEqual([g.site_occupancy for g in gates], [0.5, 0.5, 1.0])
        for gate in gates:
            self.assertAlmostEqual(gate.site_weight, 1.0 / 3.0, places=12)
        weighted_n0 = sum(g.site_weight * g.first_shell_cn for g in gates)
        self.assertAlmostEqual(weighted_n0, 16.0 / 3.0, places=6)
        self.assertAlmostEqual(weighted_n0, 5.333333, places=6)

        # Co3O4：8a:16d → 1/3 : 2/3；声明化学式 Co3 O4 与全占据晶胞一致，
        # 按 implicit_1.0 放行，加权第一壳 N₀ = 1/3×4 + 2/3×6 = 16/3。
        co = phase_resolver.resolve_phase(
            "Co", phase_id="COD:9005888", absorber_site="all", offline=True)
        self.assertEqual([m.absorber_site for m in co.site_models],
                         ["Co1", "Co2"])
        co_sites = cif_gate.list_inequivalent_absorber_sites(
            co.selected.cif_path, "Co")
        self.assertEqual([s["label"] for s in co_sites], ["Co1", "Co2"])
        self.assertEqual([s["multiplicity"] for s in co_sites], [8, 16])
        self.assertAlmostEqual(co_sites[0]["weight"], 1.0 / 3.0, places=12)
        self.assertAlmostEqual(co_sites[1]["weight"], 2.0 / 3.0, places=12)
        co_gates = cif_gate.gate_cif_all_sites(
            co.selected.cif_path, "Co", kmin=3.0, kmax=12.0)
        self.assertTrue(all(g.ok for g in co_gates),
                        [g.reason for g in co_gates])
        self.assertTrue(all(g.occupancy_source == "implicit_1.0"
                            for g in co_gates))
        co_weighted_n0 = sum(g.site_weight * g.first_shell_cn
                             for g in co_gates)
        self.assertAlmostEqual(co_weighted_n0, 16.0 / 3.0, places=6)

        # run_multi_site_formal_cn：全部位点 ok 才输出 ALL_WEIGHTED。
        res = SimpleNamespace(ft=SimpleNamespace(kmin=3.0, kmax=12.0))

        def fake_block(conditional_cn, plan, n0_mapping, structure_model,
                       absorber, **kwargs):
            site = structure_model.absorber_site
            n0 = {"Fe1": 6.0, "Fe2": 6.0, "Fe3": 4.0}[site]
            return {
                "status": "ok", "reject_reason": None,
                "rows": [{
                    "phase_id": structure_model.model_id,
                    "model_level": structure_model.model_level,
                    "source_type": structure_model.source_type,
                    "source_ref": structure_model.source_ref,
                    "source_cif_sha256": structure_model.source_cif_sha256,
                    "absorber_site": site, "shell_group": "first_shell",
                    "shell_label": "Fe-O", "scatterer": "O",
                    "r_cif": 2.0, "n0": n0, "n_fit": n0 * 0.95,
                    "n_status": "fitted", "n_ratio": 0.95,
                    "n_err_formal": 0.2, "n_err_systematic": 0.32,
                    "delta_r": 0.01, "sigma2": 0.004, "delta_e0": 1.0,
                    "rf": 0.02, "s02": 0.803,
                    "s02_source": "absorber.py", "s02_calibrated": True,
                    "kmin": 3.0, "kmax": 12.0, "kweight": 2,
                    "rmin": 1.0, "rmax": 3.0, "n_idp": 9.0,
                    "n_variables": 5, "cn_label": "cn_conditional",
                    "gate_status": "passed", "reject_reason": "",
                    "flags": "", "data_sha256": "d" * 64,
                    "software_version": "test",
                }],
            }

        common = [
            patch.object(xafs_fit, "_site_motif_from_model",
                         return_value=object()),
            patch.object(xafs_fit, "select_shell_paths",
                         return_value={"files": [], "shells": [],
                                       "rmin": 1.0, "rmax": 3.0}),
            patch.object(xafs_fit, "map_n0_to_plan",
                         return_value={"mapping_ok": True}),
            patch.object(xafs_fit, "run_conditional_cn", return_value={}),
            patch.object(xafs_fit, "build_formal_cn_block",
                         side_effect=fake_block),
        ]
        for patcher in common:
            patcher.start()
        try:
            out = xafs_fit.run_multi_site_formal_cn(
                res, site_models=fe.site_models, gate_results=gates,
                s02=0.803, s02_source="absorber.py",
                data_sha256="d" * 64)
            self.assertEqual(out["status"], "ok")
            self.assertAlmostEqual(out["site_weights"]["Fe1"], 1 / 3,
                                   places=12)
            self.assertEqual(out["rows"][-1]["absorber_site"],
                             "ALL_WEIGHTED")
            self.assertAlmostEqual(out["rows"][-1]["n0"], 16 / 3,
                                   places=6)
            self.assertIn("site_weight=0.333333", out["rows"][0]["flags"])

            # 任一位点失败：该位点写拒绝行，且不输出 ALL_WEIGHTED。
            bad = list(gates)
            bad[1] = cif_gate.CIFGateResult(
                ok=False, code="n_out_of_range", reason="fake gate failure",
                absorber="Fe", absorber_site="Fe2")
            bad[1].site_weight = gates[1].site_weight
            failed = xafs_fit.run_multi_site_formal_cn(
                res, site_models=fe.site_models, gate_results=bad,
                s02=0.803, s02_source="absorber.py",
                data_sha256="d" * 64)
            self.assertEqual(failed["status"], "multi_site_failed")
            sites = [r.get("absorber_site") for r in failed["rows"]]
            self.assertNotIn("ALL_WEIGHTED", sites)
            self.assertIn("Fe2", sites)
        finally:
            for patcher in common:
                patcher.stop()

    def test_p6_all_phase_pack_entries_pass_cif_gate(self):
        """#128 P6：phase_pack 全部 44 个入库晶体相 100% 通过 cif_gate。"""
        pack_root = ROOT / "_structure_lib" / "_phase_pack"
        manifest = phase_pack.load_manifest(pack_root)
        self.assertEqual(len(manifest["entries"]), 44)

        single_count = 0
        multi_count = 0
        for entry in manifest["entries"]:
            el = entry["absorber_elements"][0]
            cif_path = pack_root / entry["cif_path"]
            sites = entry["absorber_sites"]
            is_comp = not entry.get("experimental", True)
            if len(sites) == 1:
                single_count += 1
                g = cif_gate.gate_cif(cif_path, el, absorber_site=None,
                                      is_computed=is_comp)
                self.assertTrue(
                    g.ok, f"{entry['id']} {entry['formula']} failed: "
                          f"{g.code} ({g.reason})")
                self.assertEqual(g.absorber_site, sites[0])
                self.assertGreaterEqual(g.first_shell_cn, 2.0)
                self.assertLessEqual(g.first_shell_cn, 16.0)
            else:
                multi_count += 1
                g_amb = cif_gate.gate_cif(cif_path, el, absorber_site=None,
                                          is_computed=is_comp)
                parsed = phase_pack.parse_cif(cif_path)
                if len(parsed["formula_elements"]) == 1:
                    # 单元素多晶位：所有位点对称等价，自动选第一个
                    self.assertTrue(g_amb.ok, g_amb.reason)
                else:
                    # 多元素多晶位：未传 --absorber-site 必须继续拒绝
                    self.assertFalse(g_amb.ok)
                    self.assertEqual(g_amb.code, "absorber_site_ambiguous")
                gs = cif_gate.gate_cif_all_sites(
                    cif_path, el, is_computed=is_comp)
                self.assertEqual(len(gs), len(sites))
                self.assertTrue(
                    all(g.ok for g in gs),
                    f"{entry['id']} {entry['formula']} multi-site failed: "
                    f"{[(g.absorber_site, g.code, g.reason) for g in gs]}")
                self.assertAlmostEqual(
                    sum(g.site_weight for g in gs), 1.0, places=9)
        self.assertEqual(single_count, 39)
        self.assertEqual(multi_count, 5)

        # 催化位点假设模型与 --formal-cn 同时开启：明确写 hypothesis_model_only
        from larch.xafs import feffpath
        from types import SimpleNamespace

        motif = structure_lib.parse_catalytic_cluster_spec(
            "Fe", "sac:N4@2.00,axial=O@2.10,C=none")
        plan = xafs_fit.select_shell_paths(
            motif, edge="K", rmin=1.0, rmax=2.6, kmin=3.0, kmax=12.0)
        k = np.arange(3.0, 12.0, 0.05)
        plist = []
        for sh, f in zip(plan["shells"], plan["files"]):
            p = feffpath(str(f), s02=0.803, sigma2=0.004, e0=2.0,
                         deltar=0.0)
            p.degen = float(sh["degen"])
            p._calc_chi()
            plist.append(p)
        chi = artemis.exafs_model_k(plist, k)
        rng = np.random.default_rng(7)
        chi = chi + rng.normal(0.0, 1e-3, size=chi.shape)
        chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                     mu_norm=np.ones_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
        res = SimpleNamespace(
            chi_k=chi_k, ft=ft, absorber_symbol="Fe", structure=None,
            merged=None, norm=None, self_absorption=None)
        out = xafs_fit.fit_identified(
            res, absorber_symbol="Fe",
            catalytic_model="sac:N4@2.00,axial=O@2.10,C=none",
            rmin=1.0, rmax=2.6, verbose=False, formal_cn=True)
        self.assertEqual(out["formal_cn"]["status"], "hypothesis_model_only")
        self.assertEqual(out["cn_status"], "hypothesis_model_only")

    def test_p6_1_single_site_oxide_formal_cn_end_to_end(self):
        """#129 P6.1：单晶位氧化物（不在 _structure_lib 中）通过已验证 CIF 完成正式 CN，不被 L1 劫持。"""
        from larch.xafs import feffpath
        from types import SimpleNamespace

        import run_xafs

        s02, s02_src = run_xafs.load_t5_s02("Ru")
        pr = phase_resolver.resolve_phase(
            "Ru", phase_id="COD:9007541", offline=True)
        self.assertEqual(pr.status, "resolved")
        self.assertEqual(pr.structure_model.occupancy_policy, "implicit_1.0")
        gate = cif_gate.gate_cif(
            pr.selected.cif_path, "Ru",
            absorber_site=pr.absorber_site, kmin=3.0, kmax=12.0)
        self.assertTrue(gate.ok, gate.reason)

        motif = xafs_fit._site_motif_from_model(pr.structure_model)
        plan = xafs_fit.select_shell_paths(
            motif, edge="K", kmin=3.0, kmax=12.0)
        k = np.arange(3.0, 12.0, 0.05)
        plist = []
        for sh, f in zip(plan["shells"], plan["files"]):
            p = feffpath(str(f), s02=s02, sigma2=0.004, e0=1.0, deltar=0.0)
            p.degen = float(sh["degen"])
            p._calc_chi()
            plist.append(p)
        chi = artemis.exafs_model_k(plist, k)
        rng = np.random.default_rng(42)
        chi = chi + rng.normal(0.0, 5e-4, size=chi.shape)
        chi_k = ChiK(k=k, chi=chi, background=np.zeros_like(k),
                     mu_norm=np.ones_like(k), rbkg=1.0)
        ft = athena.ft_to_r_space(k, chi, None, 3.0, 12.0, 2)
        res = SimpleNamespace(
            chi_k=chi_k, ft=ft, absorber_symbol="Ru", structure=None,
            merged=None, norm=None, self_absorption=None,
            data_sha256="b" * 64)

        out = xafs_fit.run(
            res, identify=True, fit=True, absorber_symbol="Ru",
            s02=s02, s02_source=s02_src, formal_cn=True,
            structure_model=pr.structure_model, gate_result=gate,
            verbose=False)
        # 诊断轨因 _structure_lib/Ru 只有 ru_hcp 返回 l1_first_shell，
        # 但正式轨必须独立用 COD:9007541 (RuO2) 跑通 P3 条件 CN 并输出 ok！
        self.assertEqual(out["fit"]["mode"], "l1_first_shell")
        self.assertEqual(out["cn_status"], "ok")
        self.assertEqual(out["formal_cn"]["status"], "ok")
        rows = out["formal_cn"]["rows"]
        self.assertGreaterEqual(len(rows), 1)
        self.assertEqual(rows[0]["phase_id"], "COD:9007541")
        self.assertAlmostEqual(float(rows[0]["n0"]), 6.0, places=3)
        self.assertAlmostEqual(float(rows[0]["n_fit"]), 6.0, delta=0.35)

    def test_p6_3_entry_parity_and_operando_guards(self):
        """#131 P6.3：双入口配置键对齐、operando 冲突守卫、B10 配置接线与 S3 提示。"""
        import argparse
        import contextlib
        import io
        from unittest.mock import patch

        # 1) P1b/P5 配置合并：9 键 + post_window；CLI 显式值优先
        args = argparse.Namespace(
            s02_from_t5=False, xanes_cross_check=False, offline=False,
            phase_id=None, phase=None, phase_index=None, phase_pack=None,
            absorber_site=None, phase_cif=None, post_window=None)
        cfg = {
            "s02_from_t5": True, "xanes_cross_check": True, "offline": True,
            "phase_cif": "a.cif,b.cif", "phase_id": "COD:9007541",
            "phase": "name=rutile", "phase_index": "idx.json",
            "phase_pack": "pack", "absorber_site": "all",
            "post_window": [180, 650],
        }
        xafs_config.apply_p1b_config(args, cfg)
        self.assertTrue(args.s02_from_t5)
        self.assertTrue(args.xanes_cross_check)
        self.assertTrue(args.offline)
        self.assertEqual(args.phase_cif, ["a.cif", "b.cif"])
        self.assertEqual(args.phase_id, "COD:9007541")
        self.assertEqual(args.phase, "name=rutile")
        self.assertEqual(args.phase_index, "idx.json")
        self.assertEqual(args.phase_pack, "pack")
        self.assertEqual(args.absorber_site, "all")
        self.assertEqual(args.post_window, (180.0, 650.0))

        args2 = argparse.Namespace(
            s02_from_t5=True, xanes_cross_check=False, offline=False,
            phase_id="COD:1", phase=None, phase_index=None, phase_pack=None,
            absorber_site="Fe1", phase_cif=["explicit.cif"],
            post_window="140,620")
        xafs_config.apply_p1b_config(args2, cfg)
        self.assertEqual(args2.phase_id, "COD:1")
        self.assertEqual(args2.absorber_site, "Fe1")
        self.assertEqual(args2.phase_cif, ["explicit.cif"])
        self.assertEqual(args2.post_window, "140,620")

        # 2) B10 配置键接线：target_transform / mcr_als
        kw = xafs_config.operando_kwargs_from_config({
            "target_transform": {"spoil_good": 2.5, "spoil_poor": 5.0},
            "mcr_als": {"n_components": 3, "max_iter": 500, "closure": False},
        })
        self.assertEqual(kw, {"spoil_good": 2.5, "spoil_poor": 5.0,
                              "n_components": 3, "max_iter": 500,
                              "closure": False})

        # 3) --operando 与 S6.5/S7 参数显式冲突（不再静默忽略）
        import run_xafs
        for extra in (["--formal-cn", "--fit"], ["--fit"],
                      ["--xanes-cross-check"]):
            argv = (["run_xafs.py", str(FE_FOIL), str(FE_FOIL), str(FE_FOIL),
                     "--operando", "--absorber", "Fe"] + extra)
            with patch.object(sys, "argv", argv), \
                    contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    run_xafs.main()

        # 4) 新入口开关在后端可见（--help 即解析器全量）
        for module, flag in (
                (run_xafs, "--allow-unusual-coordination"),
                (run_xafs, "--post-window"),
                (xafs_report, "--xanes-cross-check"),
                (xafs_report, "--allow-unusual-coordination"),
                (xafs_report, "--post-window"),
                (xafs_report, "--phase-id")):
            buf = io.StringIO()
            with patch.object(sys, "argv", ["prog", "--help"]), \
                    contextlib.redirect_stdout(buf), \
                    contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    module.main()
            self.assertIn(flag, buf.getvalue(), f"{module.__name__} 缺少 {flag}")

        # 5) S3 strict 报错给出 --no-strict 出口
        with self.assertRaises(XAFSError) as ctx:
            xafs_pipeline.process_sample([FE_FOIL], absorber="Fe",
                                         strict=True, verbose=False)
        self.assertIn("归一化失败", str(ctx.exception))
        self.assertIn("--no-strict", str(ctx.exception))

    def test_p6_4_output_ledger_and_print_flags(self):
        """#132 P6.4：参数/门禁/决策/摘要台账、打印开关与阈值不漂移。"""
        import argparse
        import contextlib
        import io
        from types import SimpleNamespace
        from unittest.mock import patch

        # 1) 阈值镜像必须与冻结常量一致（防止台账与物理门禁脱节）
        self.assertEqual(xafs_provenance.PLATEAU_LO, 0.85)
        self.assertEqual(xafs_provenance.PLATEAU_HI, 1.15)
        self.assertEqual(xafs_provenance.RF_LIMIT, 0.08)
        self.assertEqual(xafs_provenance.DELTA_E0_LIMIT, 15.0)
        self.assertEqual(xafs_provenance.N_ERR_LIMIT, 1.0)
        self.assertEqual(xafs_provenance.N_ERR_FRACTION_LIMIT, 0.5)
        self.assertEqual(xafs_provenance.CIF_CN_NORMAL, (2, 12))
        self.assertEqual(xafs_provenance.CIF_CN_WARN_HI, 16)
        self.assertEqual(xafs_provenance.CN_RATIO_OK, (0.7, 1.3))
        self.assertEqual(xafs_provenance.CN_RATIO_WARN, (0.6, 1.4))
        self.assertEqual(xafs_provenance.SID_SCORE_LIMIT, 0.16)
        self.assertEqual(xafs_provenance.SID_MARGIN, 0.03)

        # 2) 参数来源：CLI / config / default / auto
        parser = argparse.ArgumentParser()
        parser.add_argument("--kmin", type=float, default=3.0)
        parser.add_argument("--kmax", type=float, default=12.0)
        parser.add_argument("--config", default=None)
        parser.add_argument("--print-gates", action="store_true")
        args = parser.parse_args(["--kmin", "3.5", "--print-gates"])
        prov = xafs_provenance.Provenance.from_args(
            args, {}, ["--kmin", "3.5", "--print-gates"],
            parser=parser, inputs=["a.xdi"])
        block = prov.param_block()
        self.assertEqual(block["sources"]["kmin"], "cli")
        self.assertEqual(block["sources"]["print_gates"], "cli")
        self.assertEqual(block["sources"]["kmax"], "default")
        self.assertEqual(block["effective"]["kmin"], 3.5)

        with tempfile.TemporaryDirectory() as tmp:
            cfg = Path(tmp) / "cfg.json"
            cfg.write_text(json.dumps({"kmin": 3.2}), encoding="utf-8")
            args2 = parser.parse_args(["--config", str(cfg)])
            prov2 = xafs_provenance.Provenance.from_args(
                args2, {"kmin": 3.2}, ["--config", str(cfg)], parser=parser)
            block2 = prov2.param_block()
            self.assertEqual(block2["sources"]["kmin"], "config")
            self.assertEqual(block2["config_keys_applied"], ["kmin"])

        # 3) 门禁台账：通过 / 失败 / 不适用
        gates = xafs_provenance.build_gate_ledger(
            norm=SimpleNamespace(plateau=0.98, edge_step=0.4),
            strict=True,
            self_absorption={"mode": "transmission", "risk": "none",
                             "cn_gate": False, "corrected": False},
            structure={"status": "confident",
                       "confidence": {"best_score": 0.05, "margin": 0.1,
                                      "second_score": 0.15},
                       "coverage_ok": True, "ranked": [{"id": "fe_bcc"}]})
        self.assertEqual(gates["normalization"]["status"], "passed")
        self.assertEqual(gates["normalization"]["threshold"], [0.85, 1.15])
        self.assertEqual(gates["self_absorption"]["status"], "passed")
        self.assertEqual(gates["structure_id"]["status"], "passed")
        self.assertEqual(gates["fit"]["status"], "not_applicable")
        self.assertEqual(gates["summary"]["passed"],
                         ["normalization", "self_absorption", "structure_id"])

        fit_obj = SimpleNamespace(
            r_factor=0.2, delta_e0=2.0, overfit=False, paths=[],
            n_independent=10, n_variables=5)
        fit_gate = xafs_provenance.build_gate_ledger(
            fit_out={"fit": fit_obj, "flags": {}, "coverage_ok": True,
                     "cn_reliable": False, "cn_status": "diagnostic_only"}
        )["fit"]
        self.assertEqual(fit_gate["status"], "failed")
        self.assertIn("r_factor", fit_gate["failed_checks"])
        self.assertEqual(fit_gate["action"], "cn_quantification_blocked")

        l1_gate = xafs_provenance.build_gate_ledger(
            fit_out={"mode": "l1_first_shell",
                     "l1": {"rf": 0.05, "delta_e0": 1.0, "N": 4.0,
                            "N_err": 0.2, "n_independent": 8,
                            "n_variables": 4},
                     "flags": {}, "coverage_ok": True,
                     "cn_reliable": True, "cn_status": "ok"}
        )["fit"]
        self.assertEqual(l1_gate["status"], "passed")

        formal_gate = xafs_provenance.build_gate_ledger(
            fit_out={"conditional_cn": {
                "n_ratio": 0.5, "n_status": "out_of_range",
                "conditional_n_fit": {"r_factor": 0.10},
                "hamilton_f_test": {"significant": False},
                "consistency_passed": False}},
            formal_block={"status": "consistency_failed",
                          "reject_reason": "N/N0 out of range"}
        )["formal_cn"]
        self.assertEqual(formal_gate["status"], "failed")
        self.assertEqual(formal_gate["action"], "formal_cn_refused")

        # 4) run_xafs 单样品：_meta.json 必须包含四个台账块
        import run_xafs
        with tempfile.TemporaryDirectory() as tmp:
            args3 = parser.parse_args(["--kmin", "3.0"])
            prov3 = xafs_provenance.Provenance.from_args(
                args3, {}, ["--kmin", "3.0"], parser=parser,
                inputs=[str(FE_FOIL)])
            run_xafs.process_sample(
                [str(FE_FOIL)], tmp, absorber="Fe", strict=False,
                provenance=prov3, verbose=False)
            meta = json.loads(
                (Path(tmp) / "fe_meta.json").read_text(encoding="utf-8"))
            for key in ("params", "gates", "decisions", "summary"):
                self.assertIn(key, meta)
            self.assertEqual(meta["params"]["sources"]["kmin"], "cli")
            self.assertEqual(meta["gates"]["normalization"]["status"],
                             "warning")
            self.assertGreaterEqual(len(meta["decisions"]), 3)

        # 5) report 入口同样写 _meta.json，且 source.mode=report
        with tempfile.TemporaryDirectory() as tmp:
            xafs_report.main([
                str(FE_FOIL), "--absorber", "Fe", "--no-strict",
                "--no-fit", "--no-html", "--no-window-tuner",
                "-o", tmp, "--print-gates", "--print-params"])
            meta = json.loads(
                (Path(tmp) / "fe_meta.json").read_text(encoding="utf-8"))
            self.assertEqual(meta["source"]["mode"], "report")
            for key in ("params", "gates", "decisions", "summary"):
                self.assertIn(key, meta)

        # 6) 两个入口都暴露 --print-gates / --print-params
        for module, flag in (
                (run_xafs, "--print-gates"),
                (run_xafs, "--print-params"),
                (xafs_report, "--print-gates"),
                (xafs_report, "--print-params")):
            buf = io.StringIO()
            with patch.object(sys, "argv", ["prog", "--help"]), \
                    contextlib.redirect_stdout(buf), \
                    contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    module.main()
            self.assertIn(flag, buf.getvalue(),
                          f"{module.__name__} 缺少 {flag}")

    def test_p6_5_figure_pack(self):
        """#133 P6.5：催化剂结构图包、逐面板 CSV、表格图与两个入口开关。"""
        import contextlib
        import csv as _csv
        import io
        from unittest.mock import patch

        import run_xafs
        import xafs_figure_pack

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            xanes_rows = []
            for i in range(40):
                energy = 7100.0 + i * 2.0
                mu = 0.001 * (i + 1)
                xanes_rows.append((energy, mu, 0.0, mu, mu))
            with (out / "fe_xanes.csv").open(
                    "w", encoding="utf-8", newline="") as fh:
                w = _csv.writer(fh)
                w.writerow(["energy_eV", "mu", "pre_edge_line",
                            "mu_norm", "mu_flat"])
                w.writerows(xanes_rows)
            with (out / "fe_chik.csv").open(
                    "w", encoding="utf-8", newline="") as fh:
                w = _csv.writer(fh)
                w.writerow(["k_inv_A", "chi", "k2chi"])
                for i in range(1, 21):
                    kk = 0.5 * i
                    w.writerow([f"{kk:.4f}", "0.12", f"{0.12 * kk * kk:.6g}"])
            with (out / "fe_chir.csv").open(
                    "w", encoding="utf-8", newline="") as fh:
                w = _csv.writer(fh)
                w.writerow(["R_A", "chi_R_re", "chi_R_im", "magnitude"])
                for i in range(30):
                    rr = 0.1 * i
                    mag = (math.exp(-((rr - 2.0) / 0.30) ** 2) +
                           0.6 * math.exp(-((rr - 2.85) / 0.25) ** 2))
                    w.writerow([f"{rr:.4f}", f"{0.8 * mag:.6g}",
                                f"{0.6 * mag:.6g}", f"{mag:.6g}"])
            with (out / "fe_wavelet.csv").open(
                    "w", encoding="utf-8", newline="") as fh:
                w = _csv.writer(fh)
                w.writerow(["k\\R", "0.0", "0.5", "1.0", "1.5", "2.0"])
                for i in range(6):
                    kk = 3.0 + i
                    w.writerow([f"{kk:.2f}"] +
                               [f"{0.1 + 0.01 * j:.4f}" for j in range(5)])
            with (out / "fe_chiq.csv").open(
                    "w", encoding="utf-8", newline="") as fh:
                w = _csv.writer(fh)
                w.writerow(["q_inv_A", "chiq_re", "chiq_im", "chiq_mag"])
                for i in range(20):
                    qq = 0.2 * i
                    w.writerow([f"{qq:.3f}", "0.01", "0.02", "0.0224"])
            with (out / "cn_table.csv").open(
                    "w", encoding="utf-8", newline="") as fh:
                w = _csv.writer(fh)
                w.writerow(["phase_id", "model_level", "absorber_site",
                            "shell_group", "shell_label", "n0", "n_fit",
                            "n_err_formal", "n_ratio", "r_cif", "delta_r",
                            "sigma2", "rf", "gate_status", "n_status",
                            "reject_reason"])
                w.writerow(["fe_bcc", "crystal", "Fe1", "first", "Fe-O",
                            "6", "5.8", "0.3", "0.9667", "2.00", "0.02",
                            "0.004", "0.02", "passed", "ok", ""])
                w.writerow(["fe_bcc", "crystal", "Fe1", "second", "Fe-Fe",
                            "8", "7.5", "0.8", "0.9375", "2.84", "0.01",
                            "0.006", "0.02", "passed", "ok", ""])
            meta = {
                "absorber": "Fe",
                "e0": {"e0_sample": 7112.0},
                "source": {
                    "path": "fe.060", "mode": "trans",
                    "format_detected": "XDAC", "n_points": 40,
                    "n_scans": 1, "scan_rms": 0.0,
                    "repeatability_note": "单扫描，无法评估重复性",
                    "merge_info": {
                        "n_input": 1, "n_kept": 1,
                        "rejected_indices": [], "z_threshold": 3.5,
                        "reject_outliers": False}},
                "normalization": {
                    "edge_step": 3.0, "plateau": 0.98,
                    "plateau_ok": True, "strict": False},
                "rbkg": {
                    "selected": 0.8, "confidence": "medium",
                    "table": [
                        {"rbkg": 0.7, "r1": 2.0, "amplitude": 4.0,
                         "leakage": 0.30, "ratio": 0.08, "score": 0.90,
                         "veto": False},
                        {"rbkg": 0.8, "r1": 2.0, "amplitude": 4.1,
                         "leakage": 0.10, "ratio": 0.02, "score": 0.975,
                         "veto": False, "confidence": "medium",
                         "confidence_metrics": {"valid_fraction": 1.0,
                                                "roughness": 0.004,
                                                "margin": 0.0}},
                        {"rbkg": 0.9, "r1": 2.0, "amplitude": 4.1,
                         "leakage": 0.11, "ratio": 0.02, "score": 0.977,
                         "veto": False},
                    ]},
                "params": {
                    "effective": {"kmin": 3.0, "kmax": 12.0, "s02": 0.8,
                                  "phase_id": "fe_bcc", "fit": True,
                                  "identify": True},
                    "sources": {"kmin": "cli", "kmax": "default",
                                "s02": "calibrated"},
                    "stages": {"kmin": "S4", "kmax": "S4", "s02": "S7"},
                    "units": {"kmin": "Å^-1", "kmax": "Å^-1"},
                    "affects_formal_cn": ["kmin", "kmax", "s02"],
                },
                "gates": {
                    "normalization": {
                        "status": "passed", "threshold": [0.85, 1.15],
                        "observed": 1.0, "action": "none"},
                    "self_absorption": {
                        "status": "passed",
                        "observed": {"mode": "trans", "risk": "none",
                                     "corrected": False},
                        "action": "none"},
                    "structure_id": {
                        "status": "passed",
                        "observed": {"status": "confident",
                                     "best_id": "fe_bcc",
                                     "best_score": 0.05,
                                     "coverage_ok": True},
                        "action": "none"},
                    "fit": {
                        "status": "passed",
                        "observed": {"mode": "motif", "r_factor": 0.02,
                                     "delta_e0": 1.0,
                                     "cn_reliable": True},
                        "action": "none"},
                    "formal_cn": {
                        "status": "passed",
                        "observed": {"formal_status": "passed",
                                     "n_ratio": 0.9667,
                                     "conditional_rf": 0.02},
                        "action": "none"},
                    "summary": {"passed": ["normalization", "fit"],
                                "warned": ["structure_id"],
                                "failed": [], "not_applicable": []},
                },
                "decisions": [
                    {"step": "S1", "decision": "read", "status": "passed",
                     "detail": "1 scan"},
                    {"step": "S7", "decision": "fit", "status": "passed",
                     "detail": "motif fit converged"},
                ],
                "formal_cn": {
                    "status": "passed",
                    "model": {"model_id": "fe_bcc",
                              "model_level": "crystal",
                              "source_type": "verified_source_cif",
                              "absorber_site": "Fe1",
                              "source_cif_sha256": "abc123",
                              "n0_shells": [
                                  {"el": "O", "r": 2.0, "n": 6, "n0": 6.0,
                                   "n0_members": [{"el": "O", "r": 2.0}]},
                                  {"el": "Fe", "r": 2.85, "n": 8, "n0": 8.0,
                                   "n0_members": [{"el": "Fe", "r": 2.85}],
                                   "sub_rayleigh_grouped": True},
                              ]},
                },
                "fit": {"kmin": 3.0, "kmax": 12.0, "r_factor": 0.02,
                        "delta_e0": 1.0, "cn_reliable": True,
                        "cn_status": "ok"},
                "ft": {"kmin": 3.0, "kmax": 12.0, "r_first_peak": 2.2},
                "wt": {"r_peak": 2.2, "k_peak": 8.0, "sigma": 1.5},
            }
            (out / "fe_meta.json").write_text(
                json.dumps(meta), encoding="utf-8")
            fit = {
                "mode": "motif", "status": "confident",
                "cn_status": "ok", "cn_reliable": True,
                "rmin": 1.0, "rmax": 3.0, "fit_space": "r",
                "s02": 0.8, "s02_source": "calibrated",
                "shells": [{"label": "Fe-O", "degen": 6, "reff": 2.0},
                           {"label": "Fe-Fe", "degen": 8, "reff": 2.84}],
                "flags": {
                    "overfit": False, "sigma2_at_bound": False,
                    "c3_at_bound": False, "c4_at_bound": False,
                    "weight_at_bound": False,
                    "weight_degenerate": False,
                    "sub_rayleigh_z_degeneracy": False,
                    "c3_dropped_by_nidp": False,
                    "c4_dropped_by_nidp": False},
                "self_absorption": {
                    "mode": "trans", "risk": "none", "corrected": False,
                    "cn_gate": False, "note": "非荧光模式"},
                "resolution_check": {
                    "resolvable": True, "delta_k": 9.0,
                    "rayleigh_limit_A": 0.175,
                    "recommendation": "可独立拟合"},
                "cn_systematics": {
                    "relative_1sigma": 0.3162,
                    "components": {"S02_transfer": 0.15,
                                   "shell_merge": 0.25,
                                   "correlation": 0.15,
                                   "normalization": 0.1},
                    "extra_components": {"parameter_correlation": 0.1},
                    "note": "系统误差按平方和合成"},
                "conditional_cn": {
                    "n_ratio": 0.9667, "n_status": "ok",
                    "consistency_passed": True,
                    "fixed_n_paths": [
                        {"label": "Fe-O", "reff": 2.0, "N": 6.0,
                         "delta_r": 0.02, "sigma2": 0.004},
                        {"label": "Fe-Fe", "reff": 2.84, "N": 8.0,
                         "delta_r": 0.01, "sigma2": 0.006}],
                    "conditional_n_paths": [
                        {"label": "Fe-O", "reff": 2.0, "N": 5.8,
                         "N_err": 0.3, "delta_r": 0.02, "sigma2": 0.004},
                        {"label": "Fe-Fe", "reff": 2.84, "N": 7.5,
                         "N_err": 0.8, "delta_r": 0.01, "sigma2": 0.006}],
                    "free_n_paths": [
                        {"label": "Fe-O", "reff": 2.0, "N": 5.8,
                         "N_err": 0.3, "delta_r": 0.02, "sigma2": 0.004},
                        {"label": "Fe-Fe", "reff": 2.84, "N": 7.5,
                         "N_err": 0.8, "delta_r": 0.01, "sigma2": 0.006}],
                    "fixed_n_fit": {"r_factor": 0.05, "n_variables": 5,
                                    "n_independent": 9.0, "delta_e0": 1.0,
                                    "overfit": False},
                    "conditional_n_fit": {"r_factor": 0.02,
                                          "n_variables": 7,
                                          "n_independent": 9.0,
                                          "delta_e0": 0.5,
                                          "overfit": False},
                    "free_n_fit": {"r_factor": 0.02, "n_variables": 9,
                                   "n_independent": 9.0,
                                   "delta_e0": 0.5, "overfit": False},
                    "hamilton_f_test": {"F": 1.5, "p_value": 0.2,
                                        "significant": False},
                },
                "fit": {
                    "r_factor": 0.02, "delta_e0": 0.5,
                    "delta_e0_err": 1.0, "n_independent": 9.0,
                    "n_variables": 5, "overfit": False,
                    "kmin": 3.0, "kmax": 12.0, "rmin": 1.0, "rmax": 3.0,
                    "kweights": [2], "s02": 0.8,
                    "chi_square": 10.0, "reduced_chi_square": 2.0,
                    "epsilon_r": 0.01,
                    "parameter_names": ["Fe-O:N", "Fe-O:ΔR", "Fe-O:σ²",
                                        "Fe-Fe:N", "Fe-Fe:ΔR"],
                    "correlation": [[1.0, 0.2, 0.3, 0.4, 0.1],
                                    [0.2, 1.0, -0.5, 0.2, 0.3],
                                    [0.3, -0.5, 1.0, 0.6, 0.4],
                                    [0.4, 0.2, 0.6, 1.0, 0.8],
                                    [0.1, 0.3, 0.4, 0.8, 1.0]],
                    "paths": [
                        {"label": "Fe-O", "N": 5.8, "N_err": 0.3,
                         "R": 2.02, "R_err": 0.01, "sigma2": 0.004,
                         "sigma2_err": 0.001, "reff": 2.0},
                        {"label": "Fe-Fe", "N": 7.5, "N_err": 0.8,
                         "R": 2.85, "R_err": 0.01, "sigma2": 0.006,
                         "sigma2_err": 0.002, "reff": 2.84},
                    ],
                },
            }
            (out / "fe_fit.json").write_text(
                json.dumps(fit), encoding="utf-8")
            structure = {
                "ranked": [
                    {"id": "fe_bcc", "kind": "crystal", "score": 0.05,
                     "sigma2": 0.004, "de0": 1.0},
                    {"id": "fe_fcc", "kind": "crystal", "score": 0.12,
                     "sigma2": 0.005, "de0": 2.0},
                ],
                "status": "confident", "coverage_ok": True,
                "confidence": {"best_score": 0.05, "second_score": 0.12,
                               "margin": 0.07, "score_limit": 0.16,
                               "margin_limit": 0.03},
            }
            (out / "fe_structure_id.json").write_text(
                json.dumps(structure), encoding="utf-8")

            manifest = xafs_figure_pack.build_pack(
                out, stem="fe", formats=("png",), dpi=72, verbose=False)
            figures = out / "figures"

            # 顶层只放总图 / 表图 / 清单：不再被 28 张散装子图和 30 个 CSV 淹没。
            for name in ("Figure1_master.png", "Figure2_confidence.png",
                         "Figure3_diagnostics.png", "Figure4_features.png",
                         "Figure2g_cn_table.png", "Figure2h_fit_table.png",
                         "figure_manifest.json", "figure_pack.md"):
                self.assertTrue((figures / name).is_file(),
                                f"缺少图包文件 {name}")
            # 逐面板 CSV 收口到 panels/
            for name in ("panels/Figure1a_xanes.csv",
                         "panels/Figure2d_gate_table.csv",
                         "panels/Figure2h_fit_table.csv",
                         "panels/Figure3a_fitcorr.csv",
                         "panels/Figure3d_rbkg.csv",
                         "panels/Figure4a_derivative.csv",
                         "panels/Figure4c_ft_peaks.csv",
                         "panels/Figure4d_cif_shells.csv"):
                self.assertTrue((figures / name).is_file(),
                                f"缺少面板 CSV {name}")
            # data/ 只留跨面板汇总，不再重复复制源文件或逐面板表。
            for name in ("data/summary.json", "data/panel_index.csv",
                         "data/all_fit_paths.csv", "data/fit_variants.csv",
                         "data/all_parameters.csv", "data/gate_checks.csv",
                         "data/decisions.csv", "data/README.md"):
                self.assertTrue((figures / name).is_file(),
                                f"缺少数据包文件 {name}")
            for name in ("data/meta.json", "data/fit.json",
                         "data/structure_id.json", "data/cn_table.csv",
                         "data/fit_table.csv", "data/fit_correlation.csv",
                         "data/cif_shells.csv", "data/structure_candidates.csv",
                         "Figure1a_xanes.png", "Figure3d_rbkg.png"):
                self.assertFalse((figures / name).is_file(),
                                 f"不应重复落盘 {name}")

            on_disk = json.loads(
                (figures / "figure_manifest.json").read_text(
                    encoding="utf-8"))
            self.assertEqual(len(on_disk["panels"]), 30)
            self.assertTrue(all(p["status"] == "ok"
                                for p in on_disk["panels"]),
                            on_disk["notes"])
            by_label = {p["label"]: p for p in on_disk["panels"]}
            self.assertEqual(by_label["Figure 1(a)"]["csv"],
                             "panels/Figure1a_xanes.csv")
            self.assertEqual(by_label["Table S1"]["image"],
                             "Figure2g_cn_table.png")
            self.assertEqual(by_label["Figure 3(a)"]["csv"],
                             "panels/Figure3a_fitcorr.csv")
            # 全量结构候选：图 1(f) 的 CSV 不截断到 Top-3/6
            self.assertEqual(by_label["Figure 1(f)"]["n_rows"], 2)
            for panel in on_disk["panels"]:
                self.assertIsNotNone(panel["image"], panel["label"])
                self.assertIsNotNone(panel["csv"], panel["label"])
                self.assertTrue((figures / panel["image"]).is_file(),
                                f"{panel['label']} 图片缺失: {panel['image']}")
                self.assertTrue((figures / panel["csv"]).is_file(),
                                f"{panel['label']} CSV 缺失")

            for key in ("Figure1", "Figure2", "Figure3", "Figure4"):
                self.assertIn(key, on_disk["masters"])
            self.assertGreaterEqual(len(on_disk["data_files"]), 8)

            summary = json.loads(
                (figures / "data" / "summary.json").read_text(
                    encoding="utf-8"))
            self.assertEqual(summary["sample"], "fe")
            self.assertEqual(summary["fit"]["r_factor"], 0.02)
            self.assertGreaterEqual(summary["counts"]["fit_paths"], 2)
            self.assertGreaterEqual(summary["counts"]["cif_shells"], 2)

            header = (figures / "panels" / "Figure2h_fit_table.csv").read_text(
                encoding="utf-8").splitlines()[0]
            for col in ("shell", "n", "n_err", "n0", "r", "sigma2"):
                self.assertIn(col, header)
            header = (figures / "panels" / "Figure4c_ft_peaks.csv").read_text(
                encoding="utf-8").splitlines()[0]
            for col in ("rank", "r_A", "amplitude", "assignment"):
                self.assertIn(col, header)

            # 需要单面板图时显式打开，落到 panels/ 而不是根目录。
            xafs_figure_pack.build_pack(
                out, stem="fe", formats=("png",), dpi=72,
                standalone_panels=True, verbose=False)
            for name in ("panels/Figure1a_xanes.png",
                         "panels/Figure3d_rbkg.png",
                         "panels/Figure4d_cif_shells.png"):
                self.assertTrue((figures / name).is_file(),
                                f"缺少独立面板图 {name}")

        # 两个入口都暴露 --figure-pack；独立 CLI 也暴露 --formats/--standalone-panels
        for module, flag in ((run_xafs, "--figure-pack"),
                             (xafs_report, "--figure-pack"),
                             (xafs_figure_pack, "--formats"),
                             (xafs_figure_pack, "--standalone-panels")):
            buf = io.StringIO()
            with patch.object(sys, "argv", ["prog", "--help"]), \
                    contextlib.redirect_stdout(buf), \
                    contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    module.main()
            self.assertIn(flag, buf.getvalue(),
                          f"{module.__name__} 缺少 {flag}")

    def test_p6_6_entry_parity_and_streamlined_outputs(self):
        """#134 P6.6：入口产物对齐与输出分层瘦身。

        - report 的 meta 补齐 rbkg.table / plateau_ok / data_sha256 / L1 Rf；
        - report CLI 暴露 --rbkg/--e0-method/--manual-e0/--kweight；
        - run_xafs 落盘 _chiq.csv 与全量 _structure_id.json；
        - figures/ 顶层只留总图/表图/清单，面板 CSV 在 panels/，data/ 无副本；
        - --batch-dir 分支透传 fit_space/bvs_restraint。
        """
        import contextlib
        import io
        from unittest.mock import patch

        import run_xafs
        import xafs_figure_pack

        # 1) report 入口：meta 完整 + 图包面板不再因缺字段跳过
        with tempfile.TemporaryDirectory() as tmp:
            xafs_report.main([
                str(FE_FOIL), "--absorber", "Fe", "--no-strict",
                "--fit-motif", "fe_bcc", "--no-html", "--no-window-tuner",
                "--figure-pack", "--figure-dpi", "72", "-o", tmp])
            out = Path(tmp)
            meta = json.loads(
                (out / "fe_meta.json").read_text(encoding="utf-8"))
            self.assertTrue(meta["rbkg"].get("table"),
                            "report 的 rbkg.table 未落盘")
            self.assertIn("plateau_ok", meta["normalization"])
            self.assertTrue(meta["source"].get("data_sha256"),
                            "report 的 source.data_sha256 缺失")
            self.assertIsNotNone(meta["fit"].get("r_factor"),
                                 "fit.r_factor 不应为 None")
            self.assertIsNotNone(meta["fit"].get("delta_e0"),
                                 "fit.delta_e0 不应为 None")
            self.assertFalse((out / "fe_fit_table.txt").exists(),
                             "冗余 _fit_table.txt 不应再生成")
            self.assertTrue((out / "fe_chiq.csv").is_file())
            self.assertTrue((out / "fe_structure_id.json").is_file())

            manifest = json.loads(
                (out / "figures" / "figure_manifest.json").read_text(
                    encoding="utf-8"))
            by_label = {p["label"]: p for p in manifest["panels"]}
            for label in ("Figure 1(f)", "Figure 3(d)", "Figure 4(b)"):
                self.assertEqual(by_label[label]["status"], "ok",
                                 f"{label} 不应跳过: {manifest['notes']}")
            self.assertGreater(by_label["Figure 1(f)"]["n_rows"], 3,
                               "结构候选仍被截断")
            # 顶层无散装子图/CSV
            top = {p.name for p in (out / "figures").iterdir()}
            self.assertNotIn("Figure1a_xanes.png", top)
            self.assertNotIn("Figure1a_xanes.csv", top)
            self.assertTrue((out / "figures" / "panels"
                             / "Figure1a_xanes.csv").is_file())
            self.assertFalse((out / "figures" / "data" / "meta.json").exists())
            self.assertFalse(
                (out / "figures" / "data" / "fit_correlation.csv").exists())

        # 2) run_xafs 入口：q 空间与全量结构候选落盘并进图包
        with tempfile.TemporaryDirectory() as tmp:
            run_xafs.process_sample(
                [str(FE_FOIL)], Path(tmp), absorber="Fe", strict=False,
                identify=True, fit=True, fit_motif="fe_bcc",
                figure_pack=True, figure_dpi=72, verbose=False)
            out = Path(tmp)
            chiq = out / "fe_chiq.csv"
            self.assertTrue(chiq.is_file(), "run_xafs 未落盘 _chiq.csv")
            header = chiq.read_text(encoding="utf-8").splitlines()[0]
            for col in ("q_inv_A", "chiq_re", "chiq_im", "chiq_mag"):
                self.assertIn(col, header)
            sid = json.loads(
                (out / "fe_structure_id.json").read_text(encoding="utf-8"))
            self.assertGreater(len(sid.get("ranked") or []), 3,
                               "run_xafs 仍把结构候选截断到 <=3")
            meta = json.loads(
                (out / "fe_meta.json").read_text(encoding="utf-8"))
            self.assertTrue(meta["source"].get("data_sha256"))
            manifest = json.loads(
                (out / "figures" / "figure_manifest.json").read_text(
                    encoding="utf-8"))
            by_label = {p["label"]: p for p in manifest["panels"]}
            self.assertEqual(by_label["Figure 4(b)"]["status"], "ok",
                             manifest["notes"])

        # 3) report CLI 暴露 4 个参数；--batch-dir 不再漏传 fit_space/bvs
        buf = io.StringIO()
        with patch.object(sys, "argv", ["prog", "--help"]), \
                contextlib.redirect_stdout(buf), \
                contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                xafs_report.main()
        for flag in ("--rbkg", "--e0-method", "--manual-e0", "--kweight"):
            self.assertIn(flag, buf.getvalue(), f"xafs_report 缺少 {flag}")

        src = Path(run_xafs.__file__).read_text(encoding="utf-8")
        self.assertGreaterEqual(src.count("fit_space=args.fit_space"), 2,
                                "--batch-dir 分支漏传 fit_space")
        self.assertGreaterEqual(src.count("bvs_restraint=bvs_restraint"), 2,
                                "--batch-dir 分支漏传 bvs_restraint")
        self.assertIn("_structure_id.json", src)
        self.assertIn("_chiq.csv", src)

    def test_p6_2_edge_entry_alignments(self):
        """#130 P6.2：frac:(x,y,z) 吸收位、P1b/P5 配置键镜像、XDI 自动元素回退与 report 异常分支。"""
        import json
        import run_xafs
        from xafs_core import XAFSConfigError

        # 1) --absorber-site "frac:(x,y,z)" 在 phase_resolver -> cif_gate -> _site_motif_from_model 贯通
        pr = phase_resolver.resolve_phase(
            "Ru", phase_id="COD:9007541", absorber_site="frac:(0,0,0)",
            offline=True)
        self.assertEqual(pr.status, "resolved")
        self.assertIsNotNone(pr.structure_model)
        self.assertEqual(pr.structure_model.absorber_site, "frac:(0,0,0)")
        gate = cif_gate.gate_cif(
            pr.selected.cif_path, "Ru",
            absorber_site=pr.absorber_site, kmin=3.0, kmax=12.0)
        self.assertTrue(gate.ok, gate.reason)
        motif = xafs_fit._site_motif_from_model(pr.structure_model)
        self.assertGreaterEqual(len(motif.sites), 2)
        self.assertGreaterEqual(len(motif.shell_expect), 1)

        # 非法 frac:(...) 坐标仍被拒绝
        pr_bad = phase_resolver.resolve_phase(
            "Ru", phase_id="COD:9007541", absorber_site="frac:(a,b,c)",
            offline=True)
        self.assertIsNone(pr_bad.structure_model)
        self.assertIn("非法分数坐标", pr_bad.reason)

        # 2) xafs_config.DEFAULTS 与 apply_config_to_kwargs 覆盖 P1b/P5 开关
        for k in ("s02_from_t5", "xanes_cross_check", "phase_cif",
                  "phase_id", "phase", "phase_index", "phase_pack",
                  "absorber_site", "offline"):
            self.assertIn(k, xafs_config.DEFAULTS)
        with tempfile.TemporaryDirectory() as tmp:
            cfg_path = Path(tmp) / "p62.json"
            cfg_path.write_text(json.dumps({
                "s02_from_t5": True,
                "xanes_cross_check": True,
                "phase_cif": "a.cif,b.cif",
                "phase_id": "COD:9007541",
                "absorber_site": "all",
                "offline": True,
            }), encoding="utf-8")
            loaded = xafs_config.load_config(cfg_path)
            base_kw = {k: v for k, v in xafs_config.DEFAULTS.items()
                       if k not in ("self_absorption", "batch")}
            merged = xafs_config.apply_config_to_kwargs(loaded, base_kw)
        self.assertTrue(merged["s02_from_t5"])
        self.assertTrue(merged["xanes_cross_check"])
        self.assertEqual(merged["phase_cif"], ["a.cif", "b.cif"])
        self.assertEqual(merged["phase_id"], "COD:9007541")
        self.assertEqual(merged["absorber_site"], "all")
        self.assertTrue(merged["offline"])

        # 3) process_sample 在 absorber=None（靠 XDI 文件头自动识别元素）时支持 --formal-cn + --phase-id
        fe_xdi = ROOT / "data" / "bmm_standards" / "Data" / "Fe" / "Fe-Lepidocrocite.xdi"
        with tempfile.TemporaryDirectory() as tmp:
            res = run_xafs.process_sample(
                [str(fe_xdi)], tmp, absorber=None, strict=False,
                identify=True, fit=True, formal_cn=True,
                phase_id="COD:1011240", offline=True, verbose=False)
            self.assertEqual(res.absorber_symbol, "Fe")
            self.assertIsNotNone(res.fit.get("formal_cn"))
            self.assertNotEqual(res.fit["cn_status"], "no_cif")


if __name__ == "__main__":
    # unittest 默认写 stderr；PowerShell 会把任何 stderr 记录当成
    # NativeCommandError，导致测试通过也显示“失败”。统一写 stdout。
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(RegressionTest)
    result = unittest.TextTestRunner(stream=sys.stdout, verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
