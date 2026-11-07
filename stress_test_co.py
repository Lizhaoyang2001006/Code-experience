# -*- coding: utf-8 -*-
"""
stress_test_co.py — Co K 边元素普适性压力测试
================================================================
验证对象（2026-11 元素泛化重构的声明）：Co K 数据喂入管线后，
要么正确运行、要么在明确护栏处报错——绝无一处静默用 Fe 值。

测试矩阵（5 阶段，全部断言收尾，任一失败即非零退出并指明阶段）：

    P1  合成已知真值的 Co K 边数据：
        注册 hcp Co motif（gemmi 展开非 Fe 元素；真实参考相，幂等可重跑）
        → FEFF6 Co K 路径 → 以 s02=0.85、σ²=0.005 合成 μ(E)
        （两列、无文件头——顺路覆盖 SIMPLE 格式与 e0_nominal 缺失路径）
    P2  S1–S6 归约（absorber=Co）：E0≈7709（绝不在 Fe 区 7112/7128）、
        edge_step≈0.5、Rbkg 落在寻优网格、第一峰 ≈2.2 Å（Co–Co 2.50 表观位）
    P3  护栏必须触发（错误信息逐条匹配）：
        a) Co 未标定     → S6.5 拒绝
        b) 无元素+无文件头 → S2 拒绝（不得静默搜 Fe 窗口）
        c) 无 s02        → S7 拒绝
    P4  内存注入 s02=0.85（不落盘 absorber.py）→ identify 走通：top1=co_hcp
    P5  CN 回路：第一壳拟合回收 P1 埋入的真值（N≈12、R≈2.50、σ²≈0.005）

运行:
    & "C:\\Users\\lizha\\.workbuddy\\binaries\\python\\versions\\3.13.12\\python.exe" stress_test_co.py

副产物（有意为之，非污染）：
    _structure_lib/Co/co_hcp.json —— 真实 hcp Co 参考相（与 Fe 种子同地位）
    xafs_out/_co_stress/*.dat     —— 合成数据留档可查
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import absorber as absorber_mod       # noqa: E402
import artemis                        # noqa: E402
import structure_lib                 # noqa: E402
import xafs_fit                      # noqa: E402
import xafs_pipeline                 # noqa: E402
from absorber import AbsorberConfig   # noqa: E402

OUT = ROOT / "xafs_out" / "_co_stress"
E0_CO = 7709.0        # Co K 标称边能（与 absorber.py 一致）
S02_GEN = 0.85        # 合成时埋入的振幅因子真值
SIGMA2_GEN = 0.005    # 合成时埋入的无序度真值 (Å²)
C_K2E = 0.2625        # k²[Å⁻²] = C·(E−E0)[eV]，与 structure_id 一致


def _msg(exc) -> str:
    return getattr(exc, "message", str(exc))


def phase1() -> tuple:
    print("== P1 注册 hcp Co motif 并合成 Co K 边数据 ==")
    motif = structure_lib.Motif(
        id="co_hcp", absorber="Co", kind="crystal", nuclearity="crystal",
        source={"type": "manual",
                "ref": "α-Co (hcp) a=2.507 c=4.070 Å, P6_3/mmc, Co@2c"},
        cell={"a": 2.507, "c": 4.070, "gamma": 120},
        spacegroup="P6_3/mmc",
        asym_sites=[{"el": "Co", "x": 1 / 3, "y": 2 / 3, "z": 0.25,
                     "wyckoff": "2c"}],
        absorber_fract=[1 / 3, 2 / 3, 0.25],
        notes="压力测试种子：hcp 金属 Co，第一壳 12×Co–Co≈2.50 Å。")
    p = structure_lib.register(motif, overwrite=True)
    m = [x for x in structure_lib.load_library("Co") if x.id == "co_hcp"][0]
    s0 = m.shell_expect[0]
    assert s0["el"] == "Co" and s0["n"] == 12 and 2.48 < s0["r"] < 2.52, s0
    print(f"  motif: {p.relative_to(ROOT)}")
    print(f"  几何指纹第一壳: Co@{s0['r']:.3f}×{s0['n']}（hcp 期望 12×2.50）")

    cache = structure_lib.ensure_feff(m)
    from larch.xafs import feffpath
    k = np.arange(0.02, 13.8, 0.035)
    paths = []
    for f in sorted(cache.glob("feff*.dat")):
        pth = feffpath(str(f), s02=S02_GEN, sigma2=SIGMA2_GEN,
                       e0=0.0, deltar=0.0)
        nleg = getattr(pth, "nleg", getattr(pth, "legs", 2))
        if nleg <= 4 and 0.5 < getattr(pth, "reff", 0) <= 4.5:
            pth._calc_chi()
            paths.append(pth)
    assert len(paths) >= 3, f"FEFF Co 路径过少: {len(paths)}"
    chi = artemis.exafs_model_k(paths, k)
    print(f"  FEFF 路径 {len(paths)} 条 → χ(k) 合成完成")

    OUT.mkdir(parents=True, exist_ok=True)
    e_post = E0_CO + k ** 2 / C_K2E
    e = np.unique(np.concatenate([np.arange(E0_CO - 160, E0_CO, 0.7), e_post]))
    made = []
    for seed in (42, 43):        # 两个"扫描"，噪声不同（S1.5 多文件路径）
        mu = (0.40 + 2.0e-5 * (e - E0_CO)                        # 边前基线
              + 0.5 * 0.5 * (1 + np.tanh((e - E0_CO) / 1.2))     # 台阶 0.5
              - 0.02 * np.log1p(np.maximum(e - E0_CO, 0) / 400)  # 慢衰减
              + 0.5 * np.interp(e, e_post, chi, left=0.0))       # EXAFS
        mu += np.random.default_rng(seed).normal(0, 8e-4, e.size)
        f = OUT / f"Co_foil_synthetic.{seed:03d}"
        np.savetxt(f, np.column_stack([e, mu]), fmt="%.6f")
        made.append(f)
    print(f"  合成数据: {made[0].relative_to(ROOT)} + .043（{e.size} 点/文件，"
          f"{e.min():.0f}–{e.max():.0f} eV，两列无文件头）")
    return made


def phase2(files: list):
    print("== P2 S1–S6 归约（absorber=Co，元素无关验证）==")
    # align=True → S1.5 需查 Co 参考边能 7709（若元素表失效此处直接炸）
    res = xafs_pipeline.process_sample([str(f) for f in files], align=True,
                                       absorber="Co", verbose=True)
    assert res.absorber_symbol == "Co", "absorber 未写入 PipelineResult 契约"
    e0 = res.e0["e0_sample"]
    assert 7700 < e0 < 7720, f"E0={e0:.2f} 偏离 Co K 7709（疑似 Fe 值泄漏）"
    assert abs(res.norm.edge_step - 0.5) < 0.10, f"edge_step={res.norm.edge_step}"
    assert 0.85 <= res.norm.plateau <= 1.15, f"plateau={res.norm.plateau}"
    assert 0.5 <= res.rbkg["rbkg"] <= 1.6, f"Rbkg={res.rbkg['rbkg']} 出网格"
    assert 1.8 < res.r1 < 2.4, f"第一峰 {res.r1}（Co–Co 2.50 表观位应 ≈2.2）"
    assert abs(res.r_wt - res.r1) < 0.1, f"WT/FT 峰位打架: {res.r_wt} vs {res.r1}"
    print(f"  ✓ E0={e0:.2f} eV（Co K ✓，非 Fe 区）  "
          f"edge_step={res.norm.edge_step:.3f}  Rbkg={res.rbkg['rbkg']:.2f}  "
          f"r1={res.r1:.3f} Å")
    return res


def phase3(data: Path, res):
    print("== P3 护栏触发（错误信息必须逐条匹配）==")
    # a) 未标定 Co：S6.5 必须拒绝。
    #    absorber.py 已写入真实箔标定值（2026-11），护栏不能依赖
    #    "生产值恰好缺"这一可变前提——内存中显式清掉 s02 再验证。
    orig = absorber_mod.ABSORBERS["Co"]
    absorber_mod.ABSORBERS["Co"] = AbsorberConfig(
        "Co", 27, "K", E0_CO, None, 15.0)
    try:
        xafs_fit.run(res, absorber_symbol="Co")
        raise AssertionError("护栏失效：未标定 Co 竟通过 identify")
    except AssertionError:
        raise
    except Exception as exc:
        assert "未标定" in _msg(exc), _msg(exc)
        print(f"  ✓ [S6.5] 未标定拒绝: {_msg(exc)}")
    finally:
        absorber_mod.ABSORBERS["Co"] = orig

    # b) 无元素 + 无文件头：S2 必须拒绝（不得静默搜 Fe 窗口）
    try:
        xafs_pipeline.process_sample([str(data)])
        raise AssertionError("护栏失效：无元素+无文件头竟跑过 S2")
    except AssertionError:
        raise
    except Exception as exc:
        assert "吸收边搜索窗口" in _msg(exc), _msg(exc)
        print(f"  ✓ [S2] 无元素+无文件头拒绝: {_msg(exc)}")

    # c) 无 s02：S7 必须拒绝
    m = [x for x in structure_lib.load_library("Co") if x.id == "co_hcp"][0]
    cache = structure_lib.ensure_feff(m)
    rows = artemis.describe_paths(sorted(cache.glob("feff*.dat")), 60)
    cands = [cache / r["file"] for r in rows
             if r["legs"] == 2 and 2.45 < r["reff"] < 2.55]
    assert cands, f"第一壳 SS 路径未找到: {rows[:6]}"
    # 不可分辨壳层合并（hcp 两子路径 ΔR=0.01 ≪ 0.26 Å）：只拟一条代表路径，
    # N 语义 = 壳层总配位数——与 Fe 基准同一条铁律，两路径各自放开必简并漂移
    first = cands[:1]
    try:
        artemis.fit_r_space(res.chi_k, res.ft, first, 1.5, 2.8)
        raise AssertionError("护栏失效：无 s02 竟通过拟合")
    except AssertionError:
        raise
    except Exception as exc:
        assert "s02" in _msg(exc), _msg(exc)
        print(f"  ✓ [S7] s02 未提供拒绝: {_msg(exc)}")
    return first


def phase4(res):
    print("== P4 内存完成接入（注入 s02=0.85 合成真值，不落盘 absorber.py）==")
    orig = absorber_mod.ABSORBERS["Co"]
    absorber_mod.ABSORBERS["Co"] = AbsorberConfig(
        "Co", 27, "K", E0_CO, S02_GEN, 15.0)
    try:
        sid = xafs_fit.run(res, absorber_symbol="Co")["structure"]
    finally:
        absorber_mod.ABSORBERS["Co"] = orig
    top = sid["ranked"][0]
    assert top["id"] == "co_hcp", f"top1 应为 co_hcp: {top}"
    assert top["score"] < 0.20, f"score 过大（理论自洽数据应很小）: {top}"
    assert not sid["uncovered_peaks"], \
        f"单相合成数据存在未解释峰（覆盖矩阵失效）: {sid['uncovered_peaks']}"
    print(f"  ✓ identify: top1=co_hcp score={top['score']:.3f} "
          f"ΔE0={top['de0']:+.0f} eV  "
          f"并列={sid['indistinguishable']}  未解释峰={sid['uncovered_peaks']}")
    return sid


def phase5(res, first: list):
    print("== P5 CN 回路（回收 P1 埋入的真值 N=12, R≈2.50, σ²≈0.005）==")
    # 单代表路径拟合（合并壳层铁律）：N 即壳层总配位数
    fit = artemis.fit_r_space(res.chi_k, res.ft, first, rmin=1.5, rmax=2.8,
                              kmin=3.0, kmax=12.0, s02=S02_GEN, vary_n=True)
    p = fit.paths[0]
    assert 10.0 < p.degeneracy < 14.0, \
        f"N={p.degeneracy:.2f}（真值 12；合并壳层语义）"
    assert 2.45 < p.r < 2.55, f"R={p.r:.4f}（真值 2.50）"
    assert 0.002 < p.sigma2 < 0.010, f"σ²={p.sigma2:.4f}（真值 0.005）"
    assert fit.r_factor < 0.10, f"Rf={fit.r_factor}（合成数据回路应 <0.10）"
    assert not fit.overfit, \
        f"变量数 {fit.n_variables} ≥ N_idp {fit.n_independent:.1f}"
    n_err = f"±{p.n_err:.2f}" if p.n_err else ""
    print(f"  ✓ N={p.degeneracy:.2f}{n_err}（真值 12）  "
          f"R={p.r:.4f} Å（真值 2.50）  σ²={p.sigma2:.4f}（真值 0.005）  "
          f"Rf={fit.r_factor:.4f}  ΔE0={fit.delta_e0:+.2f} eV")


def main():
    files = phase1()
    res = phase2(files)
    first = phase3(files[0], res)
    phase4(res)
    phase5(res, first)
    print("\n压力测试全部通过：管线对 Co K 边元素无关，护栏如预期触发，"
          "合成真值回收在容差内。")


if __name__ == "__main__":
    main()
