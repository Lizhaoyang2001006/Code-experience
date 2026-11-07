# -*- coding: utf-8 -*-
"""
stress_test_ni.py — Ni K 边元素普适性压力测试（模板：stress_test_co.py）
================================================================
验证管线对 Ni K 边元素无关：正确运行或明确护栏报错，绝无静默 Fe/Co 值。

与 Co 模板的差异点（其余逐条对应）：
    - fcc 金属 Ni（Fm-3m, a=3.524 Å）——gemmi 立方 F 心展开是新增路径
    - 第一壳为单一距离 12×Ni–Ni 2.492 Å（无 hcp 亚壳层劈裂，
      合并壳层铁律在此自然退化为单路径）
    - 真值三件套：s02=0.88、σ²=0.006、E0(Ni K)=8333 eV

测试矩阵（5 阶段，全部断言收尾，任一失败即非零退出并指明阶段）：
    P1 注册 fcc Ni motif → FEFF6 Ni K 路径 → 合成 μ(E)（两列无文件头 ×2 扫描）
    P2 S1–S6 归约（absorber=Ni）：E0≈8333（绝不在 Fe 7112 / Co 7709 区）、
       edge_step≈0.5、r1≈2.1（Ni–Ni 2.492 表观位）
    P3 三护栏逐条触发：未标定拒绝 / 无元素+无文件头 S2 拒绝 / 无 s02 S7 拒绝
    P4 内存注入 s02=0.88（不落盘 absorber.py）→ identify：top1=ni_fcc
    P5 CN 回路：回收真值 N≈12、R≈2.492、σ²≈0.006

运行:
    & "C:\\Users\\lizha\\.workbuddy\\binaries\\python\\versions\\3.13.12\\python.exe" stress_test_ni.py

副产物（有意为之，非污染）：
    _structure_lib/Ni/ni_fcc.json —— 真实 fcc Ni 参考相种子（与 Fe/Co 种子同地位）
    xafs_out/_ni_stress/*.dat     —— 合成数据留档可查
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

OUT = ROOT / "xafs_out" / "_ni_stress"
E0_NI = 8333.0        # Ni K 标称边能（与 absorber.py 一致）
S02_GEN = 0.88        # 合成时埋入的振幅因子真值
SIGMA2_GEN = 0.006    # 合成时埋入的无序度真值 (Å²)
C_K2E = 0.2625        # k²[Å⁻²] = C·(E−E0)[eV]，与 structure_id 一致


def _msg(exc) -> str:
    return getattr(exc, "message", str(exc))


def phase1() -> list:
    print("== P1 注册 fcc Ni motif 并合成 Ni K 边数据 ==")
    motif = structure_lib.Motif(
        id="ni_fcc", absorber="Ni", kind="crystal", nuclearity="crystal",
        source={"type": "manual", "ref": "fcc Ni, a=3.524 Å, Fm-3m, Ni@4a"},
        cell={"a": 3.524}, spacegroup="Fm-3m",
        asym_sites=[{"el": "Ni", "x": 0.0, "y": 0.0, "z": 0.0,
                     "wyckoff": "4a"}],
        notes="压力测试种子：fcc 金属 Ni，第一壳 12×Ni–Ni=2.492 Å 单一距离。")
    p = structure_lib.register(motif, overwrite=True)
    m = [x for x in structure_lib.load_library("Ni") if x.id == "ni_fcc"][0]
    s0 = m.shell_expect[0]
    assert s0["el"] == "Ni" and s0["n"] == 12 and 2.48 < s0["r"] < 2.50, s0
    print(f"  motif: {p.relative_to(ROOT)}")
    print(f"  几何指纹第一壳: Ni@{s0['r']:.3f}×{s0['n']}（fcc 期望 12×2.492）")

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
    assert len(paths) >= 3, f"FEFF Ni 路径过少: {len(paths)}"
    chi = artemis.exafs_model_k(paths, k)
    print(f"  FEFF 路径 {len(paths)} 条 → χ(k) 合成完成")

    OUT.mkdir(parents=True, exist_ok=True)
    e_post = E0_NI + k ** 2 / C_K2E
    e = np.unique(np.concatenate([np.arange(E0_NI - 160, E0_NI, 0.7), e_post]))
    made = []
    for seed in (42, 43):        # 两个"扫描"，噪声不同（S1.5 多文件路径）
        mu = (0.40 + 2.0e-5 * (e - E0_NI)                        # 边前基线
              + 0.5 * 0.5 * (1 + np.tanh((e - E0_NI) / 1.2))     # 台阶 0.5
              - 0.02 * np.log1p(np.maximum(e - E0_NI, 0) / 400)  # 慢衰减
              + 0.5 * np.interp(e, e_post, chi, left=0.0))       # EXAFS
        mu += np.random.default_rng(seed).normal(0, 8e-4, e.size)
        f = OUT / f"Ni_foil_synthetic.{seed:03d}"
        np.savetxt(f, np.column_stack([e, mu]), fmt="%.6f")
        made.append(f)
    print(f"  合成数据: {made[0].relative_to(ROOT)} + .043（{e.size} 点/文件，"
          f"{e.min():.0f}–{e.max():.0f} eV，两列无文件头）")
    return made


def phase2(files: list):
    print("== P2 S1–S6 归约（absorber=Ni，元素无关验证）==")
    # align=True → S1.5 需查 Ni 参考边能 8333（若元素表失效此处直接炸）
    res = xafs_pipeline.process_sample([str(f) for f in files], align=True,
                                       absorber="Ni", verbose=True)
    assert res.absorber_symbol == "Ni", "absorber 未写入 PipelineResult 契约"
    e0 = res.e0["e0_sample"]
    assert 8325 < e0 < 8345, f"E0={e0:.2f} 偏离 Ni K 8333（疑似 Fe/Co 值泄漏）"
    assert abs(res.norm.edge_step - 0.5) < 0.10, f"edge_step={res.norm.edge_step}"
    assert 0.85 <= res.norm.plateau <= 1.15, f"plateau={res.norm.plateau}"
    assert 0.5 <= res.rbkg["rbkg"] <= 1.6, f"Rbkg={res.rbkg['rbkg']} 出网格"
    assert 1.8 < res.r1 < 2.4, f"第一峰 {res.r1}（Ni–Ni 2.492 表观位应 ≈2.1）"
    assert abs(res.r_wt - res.r1) < 0.1, f"WT/FT 峰位打架: {res.r_wt} vs {res.r1}"
    print(f"  ✓ E0={e0:.2f} eV（Ni K ✓）  edge_step={res.norm.edge_step:.3f}  "
          f"Rbkg={res.rbkg['rbkg']:.2f}  r1={res.r1:.3f} Å")
    return res


def phase3(data: Path, res):
    print("== P3 护栏触发（错误信息必须逐条匹配）==")
    # a) 未标定 Ni：S6.5 必须拒绝。
    #    absorber.py 已写入真实箔标定值（2026-11），护栏不能依赖
    #    "生产值恰好缺"这一可变前提——内存中显式清掉 s02 再验证。
    orig = absorber_mod.ABSORBERS["Ni"]
    absorber_mod.ABSORBERS["Ni"] = AbsorberConfig(
        "Ni", 28, "K", E0_NI, None, 15.0)
    try:
        xafs_fit.run(res, absorber_symbol="Ni")
        raise AssertionError("护栏失效：未标定 Ni 竟通过 identify")
    except AssertionError:
        raise
    except Exception as exc:
        assert "未标定" in _msg(exc), _msg(exc)
        print(f"  ✓ [S6.5] 未标定拒绝: {_msg(exc)}")
    finally:
        absorber_mod.ABSORBERS["Ni"] = orig

    # b) 无元素 + 无文件头：S2 必须拒绝（不得静默搜 Fe/Co 窗口）
    try:
        xafs_pipeline.process_sample([str(data)])
        raise AssertionError("护栏失效：无元素+无文件头竟跑过 S2")
    except AssertionError:
        raise
    except Exception as exc:
        assert "吸收边搜索窗口" in _msg(exc), _msg(exc)
        print(f"  ✓ [S2] 无元素+无文件头拒绝: {_msg(exc)}")

    # c) 无 s02：S7 必须拒绝
    m = [x for x in structure_lib.load_library("Ni") if x.id == "ni_fcc"][0]
    cache = structure_lib.ensure_feff(m)
    rows = artemis.describe_paths(sorted(cache.glob("feff*.dat")), 60)
    cands = [cache / r["file"] for r in rows
             if r["legs"] == 2 and 2.44 < r["reff"] < 2.55]
    assert cands, f"第一壳 SS 路径未找到: {rows[:6]}"
    # fcc 第一壳本就是单一距离（12×2.492 一条路径）；沿用合并壳层约定取代表路径
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
    print("== P4 内存完成接入（注入 s02=0.88 合成真值，不落盘 absorber.py）==")
    orig = absorber_mod.ABSORBERS["Ni"]
    absorber_mod.ABSORBERS["Ni"] = AbsorberConfig(
        "Ni", 28, "K", E0_NI, S02_GEN, 15.0)
    try:
        sid = xafs_fit.run(res, absorber_symbol="Ni")["structure"]
    finally:
        absorber_mod.ABSORBERS["Ni"] = orig
    top = sid["ranked"][0]
    assert top["id"] == "ni_fcc", f"top1 应为 ni_fcc: {top}"
    assert top["score"] < 0.20, f"score 过大（理论自洽数据应很小）: {top}"
    assert not sid["uncovered_peaks"], \
        f"单相合成数据存在未解释峰（覆盖矩阵失效）: {sid['uncovered_peaks']}"
    print(f"  ✓ identify: top1=ni_fcc score={top['score']:.3f} "
          f"ΔE0={top['de0']:+.0f} eV  "
          f"并列={sid['indistinguishable']}  未解释峰={sid['uncovered_peaks']}")
    return sid


def phase5(res, first: list):
    print("== P5 CN 回路（回收 P1 埋入的真值 N=12, R≈2.492, σ²≈0.006）==")
    fit = artemis.fit_r_space(res.chi_k, res.ft, first, rmin=1.5, rmax=2.8,
                              kmin=3.0, kmax=12.0, s02=S02_GEN, vary_n=True)
    p = fit.paths[0]
    assert 10.0 < p.degeneracy < 14.0, \
        f"N={p.degeneracy:.2f}（真值 12）"
    assert 2.44 < p.r < 2.54, f"R={p.r:.4f}（真值 2.492）"
    assert 0.002 < p.sigma2 < 0.011, f"σ²={p.sigma2:.4f}（真值 0.006）"
    assert fit.r_factor < 0.10, f"Rf={fit.r_factor}（合成数据回路应 <0.10）"
    assert not fit.overfit, \
        f"变量数 {fit.n_variables} ≥ N_idp {fit.n_independent:.1f}"
    n_err = f"±{p.n_err:.2f}" if p.n_err else ""
    print(f"  ✓ N={p.degeneracy:.2f}{n_err}（真值 12）  "
          f"R={p.r:.4f} Å（真值 2.492）  σ²={p.sigma2:.4f}（真值 0.006）  "
          f"Rf={fit.r_factor:.4f}  ΔE0={fit.delta_e0:+.2f} eV")


def main():
    files = phase1()
    res = phase2(files)
    first = phase3(files[0], res)
    phase4(res)
    phase5(res, first)
    print("\n压力测试全部通过：管线对 Ni K 边元素无关，护栏如预期触发，"
          "合成真值回收在容差内。")


if __name__ == "__main__":
    main()
