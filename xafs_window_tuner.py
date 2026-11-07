# -*- coding: utf-8 -*-
"""
xafs_window_tuner.py — 交互式 FT 窗口调节器（单文件 HTML）
============================================================
对标 Athena 的「拖动 k 窗口看 |χ(R)| 怎么变」。不依赖 GUI/服务器：
把 χ(k) 与 FT 内核嵌进 HTML，浏览器里实时拖动 kmin/kmax/kweight/dk。

关键设计:
    - 浏览器端 FT 与产线同约定：Hanning 窗 + e^{-2ikR} + k^n 加权，
      只看峰位/形状，幅度按各自最大值归一（页面明确标注）。
    - 后端先用 larch.xftf 扫描一组候选窗口，按
      「峰位稳定 + 振幅大 + 不偏离产线第一峰」给出**参考**推荐窗口；
      推荐不会自动改产线，页面给出对应 CLI 命令。
    - 不改变 Rbkg/归一化：本工具只调 FT/拟合的 k 窗口。

用法（由 xafs_report 调用）:
    write_window_tuner(out, stem, chi_k, r1, ...)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import athena
from xafs_core import ChiK, XAFSConfigError

R_DISPLAY_MAX = 4.0
R_SEARCH_MIN = 0.8
R_SEARCH_MAX = 3.4
DEFAULT_GRID_KMIN = (2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0)
DEFAULT_GRID_KMAX = (9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0)


def hann_window(k, kmin: float, kmax: float, dk: float) -> np.ndarray:
    """与 larch.xftf ``window='hanning'`` 相同的分段 Hanning 窗。

    左端 kmin→kmin+dk 从 0 升到 1；中间为 1；右端 kmax-dk→kmax 降回 0。
    测试用它验证浏览器 FT 与产线窗函数一致。
    """
    k = np.asarray(k, dtype=float)
    w = np.zeros_like(k)
    dk = max(float(dk), 1e-9)
    left = (k >= kmin) & (k < kmin + dk)
    w[left] = 0.5 * (1.0 - np.cos(np.pi * (k[left] - kmin) / dk))
    mid = (k >= kmin + dk) & (k <= kmax - dk)
    w[mid] = 1.0
    right = (k > kmax - dk) & (k <= kmax)
    w[right] = 0.5 * (1.0 + np.cos(np.pi * (k[right] - kmax + dk) / dk))
    return w


def _peak_near(r, mag, r_center: float, half: float = 0.6):
    """在 r_center ± half 内找 |χ(R)| 最大值，避免跳到相邻壳。"""
    r = np.asarray(r, dtype=float)
    mag = np.asarray(mag, dtype=float)
    lo = max(R_SEARCH_MIN, float(r_center) - half)
    hi = min(R_SEARCH_MAX, float(r_center) + half)
    mask = (r >= lo) & (r <= hi)
    if mask.sum() < 3:
        mask = (r >= R_SEARCH_MIN) & (r <= R_SEARCH_MAX)
    if not mask.any():
        return float("nan"), 0.0
    i = int(np.argmax(mag[mask]))
    return float(r[mask][i]), float(mag[mask][i])


def recommend_window(chi_k: ChiK, r1: float, kweight: int = 2,
                     dk: float = 1.0,
                     kmin_grid=DEFAULT_GRID_KMIN,
                     kmax_grid=DEFAULT_GRID_KMAX) -> dict:
    """扫描候选 k 窗口，给一个**参考**推荐（不自动写入产线）。

    评分:
        只保留峰位落在 r1 ± 0.15 Å 的候选（防止把窗口推到伪峰上）；
        stability = 与相邻网格点的峰位最大差（Å）；
        score = 振幅 / (1 + 20·stability + 5·|R−r1|)。
    这样选出来的是「振幅大、且对窗口微扰不敏感」的区域，而不是
    恰好把某个噪声峰顶到最高的窗口。
    """
    r1 = float(r1) if np.isfinite(r1) else 1.8
    grid = []
    for kmin in kmin_grid:
        for kmax in kmax_grid:
            if kmax - kmin < 2.0:
                continue
            try:
                ft = athena.ft_to_r_space(
                    chi_k.k, chi_k.chi, kmin=float(kmin),
                    kmax=float(kmax), kweight=int(kweight), dk=float(dk))
            except XAFSConfigError:
                continue
            rp, amp = _peak_near(ft.r, ft.magnitude, r1)
            if not np.isfinite(rp):
                continue
            grid.append({"kmin": float(kmin), "kmax": float(kmax),
                         "peak_r": rp, "amplitude": amp})
    if not grid:
        raise XAFSConfigError("TUNER", "候选窗口全部无有效 FT 数据")

    # 邻居稳定性：同一 kmin 行/列上 ±1 格的峰位差
    kmins = sorted({g["kmin"] for g in grid})
    kmaxs = sorted({g["kmax"] for g in grid})
    lookup = {(g["kmin"], g["kmax"]): g for g in grid}
    for g in grid:
        diffs = []
        i = kmins.index(g["kmin"])
        j = kmaxs.index(g["kmax"])
        for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            if 0 <= i + di < len(kmins) and 0 <= j + dj < len(kmaxs):
                nb = lookup.get((kmins[i + di], kmaxs[j + dj]))
                if nb is not None:
                    diffs.append(abs(nb["peak_r"] - g["peak_r"]))
        g["stability"] = max(diffs) if diffs else 0.0
        g["offset"] = abs(g["peak_r"] - r1)
        g["feasible"] = g["offset"] <= 0.15
        g["score"] = (g["amplitude"]
                      / (1.0 + 20.0 * g["stability"] + 5.0 * g["offset"]))

    feasible = [g for g in grid if g["feasible"]]
    pool = feasible or grid
    ranked = sorted(pool, key=lambda g: g["score"], reverse=True)
    best = dict(ranked[0])
    top = []
    for g in ranked[:5]:
        top.append({k: (round(v, 4) if isinstance(v, float) else v)
                    for k, v in g.items()})
    return {"recommended": {k: (round(v, 4) if isinstance(v, float) else v)
                            for k, v in best.items()},
            "top": top, "n_feasible": len(feasible),
            "n_candidates": len(grid), "anchor_r": r1,
            "kweight": int(kweight), "dk": float(dk)}


def _json_safe(obj):
    """NaN/Inf → None，保证 JSON.parse 与 Python json 都能读。"""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        return v if np.isfinite(v) else None
    if isinstance(obj, np.integer):
        return int(obj)
    return obj


def build_tuner_html(chi_k: ChiK, r1: float, sample: str, edge_label: str,
                     kmin: float, kmax: float, kweight: int,
                     recommendation: dict) -> str:
    payload = _json_safe({
        "k": np.asarray(chi_k.k, dtype=float).tolist(),
        "chi": np.asarray(chi_k.chi, dtype=float).tolist(),
        "r1": float(r1) if np.isfinite(r1) else None,
        "current": {"kmin": float(kmin), "kmax": float(kmax),
                    "kweight": int(kweight), "dk": 1.0},
        "recommended": recommendation.get("recommended"),
        "top": recommendation.get("top", []),
        "n_feasible": recommendation.get("n_feasible"),
        "n_candidates": recommendation.get("n_candidates"),
        "sample": sample, "edge_label": edge_label,
    })
    data_json = json.dumps(payload, ensure_ascii=False,
                           separators=(",", ":"), allow_nan=False)
    rec = recommendation.get("recommended") or {}
    rec_text = (f"kmin={rec.get('kmin')}, kmax={rec.get('kmax')} → "
                f"R={rec.get('peak_r'):.3f} Å（稳定性 "
                f"{rec.get('stability'):.3f} Å）"
                if rec else "无推荐（候选不足）")
    top_rows = "".join(
        "<tr><td>{kmin:.2f}</td><td>{kmax:.2f}</td><td>{peak_r:.3f}</td>"
        "<td>{stability:.3f}</td><td>{amplitude:.3g}</td></tr>".format(**g)
        for g in recommendation.get("top", []))
    return _HTML_TEMPLATE.replace("__DATA__", data_json) \
        .replace("__SAMPLE__", sample) \
        .replace("__EDGE__", edge_label) \
        .replace("__REC_TEXT__", rec_text) \
        .replace("__TOP_ROWS__", top_rows) \
        .replace("__NFEAS__", str(recommendation.get("n_feasible"))) \
        .replace("__NCAND__", str(recommendation.get("n_candidates")))


def write_window_tuner(out, stem: str, chi_k: ChiK, r1: float, sample: str,
                       edge_label: str, kmin: float, kmax: float,
                       kweight: int, dk: float = 1.0) -> Path:
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rec = recommend_window(chi_k, r1, kweight=kweight, dk=dk)
    doc = build_tuner_html(chi_k, r1, sample, edge_label,
                           kmin, kmax, kweight, rec)
    path = out / f"{stem}_window_tuner.html"
    path.write_text(doc, encoding="utf-8")
    return path


_HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__SAMPLE__ — 交互式 k 窗口</title>
<style>
:root { color-scheme: light; }
* { box-sizing: border-box; }
body { margin:0; padding:22px; background:#f5f6f8; color:#1c2430;
       font:14px/1.55 -apple-system,"Segoe UI","Microsoft YaHei",sans-serif; }
main { max-width:1080px; margin:0 auto; background:#fff; padding:26px 32px;
       border-radius:10px; box-shadow:0 1px 4px rgba(0,0,0,.08); }
h1 { margin:0 0 4px; font-size:22px; }
h2 { margin:26px 0 8px; font-size:16px; border-bottom:1px solid #e3e7ec;
     padding-bottom:5px; }
.sub { color:#5a6675; margin:0 0 16px; }
.controls { display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr));
            gap:12px 22px; background:#fafbfc; border:1px solid #e3e7ec;
            border-radius:8px; padding:14px 18px; }
.ctl label { display:flex; justify-content:space-between; font-size:13px;
             color:#34404e; margin-bottom:3px; }
.ctl output { font-variant-numeric:tabular-nums; font-weight:600; }
input[type=range] { width:100%; }
.readout { display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
           gap:10px; margin:14px 0 4px; }
.card { border:1px solid #e3e7ec; border-radius:8px; padding:8px 12px;
        background:#fafbfc; }
.card .label { color:#68737f; font-size:12px; }
.card .value { font-size:17px; font-weight:600; margin-top:1px; }
canvas { width:100%; height:auto; display:block; border:1px solid #e3e7ec;
         border-radius:6px; background:#fff; }
.buttons { margin:12px 0 2px; display:flex; gap:10px; flex-wrap:wrap; }
button { border:1px solid #c9d3dd; background:#fff; color:#1c2430;
         border-radius:6px; padding:7px 14px; font-size:13px; cursor:pointer; }
button.primary { background:#1b6ca8; border-color:#1b6ca8; color:#fff; }
button:hover { filter:brightness(.96); }
pre { background:#0f1720; color:#d7e2ec; padding:10px 13px; border-radius:6px;
      overflow-x:auto; font-size:12.5px; }
.note { background:#f7f9fb; border-left:4px solid #9fb3c8; padding:9px 13px;
        border-radius:0 6px 6px 0; color:#34404e; margin:10px 0; }
table { border-collapse:collapse; width:100%; font-size:12.5px; margin-top:6px; }
th,td { border:1px solid #e1e5ea; padding:5px 8px; text-align:left; }
th { background:#f0f3f6; }
.rec { border:1px solid #b9d7c6; background:#eef8f2; color:#175c3b;
       padding:9px 13px; border-radius:6px; margin:8px 0; }
</style>
</head>
<body>
<main>
<h1>__SAMPLE__</h1>
<p class="sub">__EDGE__ · 交互式 k 窗口（FT 峰位/形状实时预览）</p>

<div class="controls">
  <div class="ctl"><label>k<sub>min</sub>（Å<sup>−1</sup>）<output id="kminOut"></output></label>
    <input type="range" id="kmin" min="1.5" max="6" step="0.05"></div>
  <div class="ctl"><label>k<sub>max</sub>（Å<sup>−1</sup>）<output id="kmaxOut"></output></label>
    <input type="range" id="kmax" min="6" max="16" step="0.1"></div>
  <div class="ctl"><label>kweight<output id="kwOut"></output></label>
    <input type="range" id="kweight" min="0" max="3" step="1"></div>
  <div class="ctl"><label>dk（Å<sup>−1</sup>）<output id="dkOut"></output></label>
    <input type="range" id="dk" min="0.3" max="3" step="0.1"></div>
</div>

<div class="readout">
  <div class="card"><div class="label">FT 第一峰 R</div><div class="value" id="peakR">—</div></div>
  <div class="card"><div class="label">峰位 vs 产线 r1</div><div class="value" id="peakShift">—</div></div>
  <div class="card"><div class="label">峰高（归一）</div><div class="value" id="peakAmp">—</div></div>
  <div class="card"><div class="label">窗口内数据点</div><div class="value" id="npoints">—</div></div>
</div>

<h2>k 空间：加权 χ(k) 与窗函数</h2>
<canvas id="kplot" width="1000" height="240"></canvas>

<h2>R 空间：|χ(R)|</h2>
<canvas id="rplot" width="1000" height="320"></canvas>

<div class="buttons">
  <button class="primary" id="applyRec">应用推荐窗口</button>
  <button id="reset">恢复产线窗口</button>
  <button id="copyCmd">复制对应 CLI 命令</button>
</div>

<div class="rec"><b>后端参考推荐</b>：__REC_TEXT__<br>
<span style="font-size:12.5px">候选 __NCAND__ 个，峰位落在产线 r1±0.15 Å 的可行窗口 __NFEAS__ 个。
推荐只考虑峰位稳定与振幅，不替代你的物理判断。</span></div>

<details>
<summary>推荐候选 Top-5（后端 larch.xftf 计算）</summary>
<table><thead><tr><th>kmin</th><th>kmax</th><th>R 峰位 (Å)</th>
<th>邻域稳定性 (Å)</th><th>相对振幅</th></tr></thead>
<tbody>__TOP_ROWS__</tbody></table>
</details>

<h2>对应产线命令</h2>
<pre id="cmd"></pre>

<div class="note">
本页只重算 FT：χ(k) 已由产线的归一化 + AUTOBK 得到。改背景请用
<code>--rbkg</code> 重跑 S1–S4；本页的 k 窗口对应
<code>--kmin/--kmax/--kweight</code>。浏览器 FT 与产线同用 Hanning 窗和
e<sup>−2ikR</sup>，幅度按各自最大值归一，故只比较峰位与形状。
</div>
</main>
<script>
const DATA = __DATA__;
const K = DATA.k, CHI = DATA.chi;
const R_GRID = [];
for (let r = 0; r <= 4.0001; r += 0.02) R_GRID.push(r);

const $ = id => document.getElementById(id);
const fmt = (v, n=2) => Number(v).toFixed(n);

function hann(k, kmin, kmax, dk) {
  if (k < kmin || k > kmax) return 0;
  if (k < kmin + dk) return 0.5 * (1 - Math.cos(Math.PI * (k - kmin) / dk));
  if (k > kmax - dk) return 0.5 * (1 + Math.cos(Math.PI * (k - kmax + dk) / dk));
  return 1;
}
function frame() {
  return {
    kmin: +$('kmin').value, kmax: +$('kmax').value,
    kw: +$('kweight').value, dk: +$('dk').value
  };
}
function transform(f) {
  const re = new Array(R_GRID.length).fill(0);
  const im = new Array(R_GRID.length).fill(0);
  let n = 0;
  for (let i = 0; i < K.length; i++) {
    const k = K[i];
    if (k < f.kmin || k > f.kmax) continue;
    const w = hann(k, f.kmin, f.kmax, f.dk);
    if (w <= 0) continue;
    n++;
    const a = Math.pow(k, f.kw) * CHI[i] * w;
    for (let j = 0; j < R_GRID.length; j++) {
      const ph = -2 * k * R_GRID[j];
      re[j] += a * Math.cos(ph);
      im[j] += a * Math.sin(ph);
    }
  }
  const mag = re.map((v, j) => Math.hypot(v, im[j]));
  const peak = Math.max(...mag) || 1;
  return {re, im, mag, peak, n};
}
function findPeak(mag, center) {
  const lo = Math.max(0.8, center - 0.6), hi = Math.min(3.4, center + 0.6);
  let best = -1, bi = 0;
  for (let j = 0; j < R_GRID.length; j++) {
    const r = R_GRID[j];
    if (r < lo || r > hi) continue;
    if (mag[j] > best) { best = mag[j]; bi = j; }
  }
  return {r: R_GRID[bi], mag: best};
}
function axes(ctx, W, H, pad, xlab, ylab) {
  ctx.clearRect(0, 0, W, H);
  ctx.strokeStyle = '#8a97a6'; ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, H - pad.b);
  ctx.lineTo(W - pad.r, H - pad.b); ctx.stroke();
  ctx.fillStyle = '#5a6675'; ctx.font = '12px sans-serif';
  ctx.fillText(xlab, W - pad.r - 40, H - 6);
  ctx.save(); ctx.translate(14, pad.t + 90); ctx.rotate(-Math.PI/2);
  ctx.fillText(ylab, 0, 0); ctx.restore();
}
function drawK(f) {
  const c = $('kplot'), ctx = c.getContext('2d');
  const W = c.width, H = c.height, pad = {l: 52, r: 14, t: 12, b: 26};
  axes(ctx, W, H, pad, 'k (Å⁻¹)', 'k^w χ(k)');
  const k0 = 0, k1 = Math.max(...K);
  let maxA = 0;
  const vals = K.map((k, i) => {
    const a = Math.pow(k, f.kw) * CHI[i];
    if (k >= f.kmin && k <= f.kmax) maxA = Math.max(maxA, Math.abs(a));
    return a;
  });
  maxA = maxA || 1;
  const x = k => pad.l + (k - k0) / (k1 - k0) * (W - pad.l - pad.r);
  const y = a => H - pad.b - (a / maxA) * 0.5 * (H - pad.t - pad.b) - 0.25*(H-pad.t-pad.b);
  // window shading
  ctx.fillStyle = 'rgba(240,160,40,0.13)';
  ctx.fillRect(x(f.kmin), pad.t, x(f.kmax) - x(f.kmin), H - pad.t - pad.b);
  ctx.strokeStyle = '#c77c1a'; ctx.lineWidth = 1.2;
  ctx.beginPath();
  const step = Math.max(1, Math.floor(K.length / 900));
  for (let i = 0; i < K.length; i += step) {
    const px = x(K[i]), py = y(vals[i]);
    i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
  }
  ctx.stroke();
  ctx.strokeStyle = '#b3261e'; ctx.setLineDash([4,3]);
  ctx.beginPath();
  for (let i = 0; i < K.length; i += step) {
    const w = hann(K[i], f.kmin, f.kmax, f.dk);
    const py = H - pad.b - w * 0.18 * (H - pad.t - pad.b);
    i === 0 ? ctx.moveTo(x(K[i]), py) : ctx.lineTo(x(K[i]), py);
  }
  ctx.stroke(); ctx.setLineDash([]);
}
function drawR(res, f, peak) {
  const c = $('rplot'), ctx = c.getContext('2d');
  const W = c.width, H = c.height, pad = {l: 52, r: 14, t: 12, b: 26};
  axes(ctx, W, H, pad, 'R (Å)', '|χ(R)| (norm.)');
  const x = r => pad.l + r / 4 * (W - pad.l - pad.r);
  const y = v => H - pad.b - v / (res.peak || 1) * (H - pad.t - pad.b);
  ctx.strokeStyle = '#1b6ca8'; ctx.lineWidth = 1.4;
  ctx.beginPath();
  for (let j = 0; j < R_GRID.length; j++) {
    const px = x(R_GRID[j]), py = y(res.mag[j]);
    j === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
  }
  ctx.stroke();
  if (DATA.r1) {
    ctx.strokeStyle = '#68737f'; ctx.setLineDash([5,4]);
    ctx.beginPath(); ctx.moveTo(x(DATA.r1), pad.t);
    ctx.lineTo(x(DATA.r1), H - pad.b); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = '#68737f'; ctx.font = '12px sans-serif';
    ctx.fillText('产线 r1=' + fmt(DATA.r1,3), x(DATA.r1)+4, pad.t+12);
  }
  ctx.strokeStyle = '#b3261e'; ctx.setLineDash([4,3]);
  ctx.beginPath(); ctx.moveTo(x(peak.r), pad.t);
  ctx.lineTo(x(peak.r), H - pad.b); ctx.stroke(); ctx.setLineDash([]);
  ctx.fillStyle = '#b3261e'; ctx.font = '13px sans-serif';
  ctx.fillText('R=' + fmt(peak.r,3) + ' Å', x(peak.r)+5, pad.t+28);
}
function update() {
  const f = frame();
  $('kminOut').value = fmt(f.kmin); $('kmaxOut').value = fmt(f.kmax);
  $('kwOut').value = f.kw; $('dkOut').value = fmt(f.dk);
  if (f.kmax - f.kmin < 1.0) {
    $('cmd').textContent = 'kmax 必须比 kmin 至少大 1 Å⁻¹';
    return;
  }
  const res = transform(f);
  const center = DATA.r1 || 1.8;
  const peak = findPeak(res.mag, center);
  $('peakR').textContent = fmt(peak.r, 3) + ' Å';
  $('peakShift').textContent = (DATA.r1 ? (peak.r - DATA.r1 >= 0 ? '+' : '') +
      fmt(peak.r - DATA.r1, 3) + ' Å' : '—');
  $('peakAmp').textContent = fmt(peak.mag / (res.peak || 1), 3);
  $('npoints').textContent = res.n;
  drawK(f); drawR(res, f, peak);
  $('cmd').textContent =
    'python xafs_report.py <files> --absorber <El> --kmin ' +
    fmt(f.kmin,2) + ' --kmax ' + fmt(f.kmax,2) +
    ' --kweight ' + f.kw;
}
function setFrame(kmin, kmax, kw, dk) {
  $('kmin').value = kmin; $('kmax').value = kmax;
  $('kweight').value = kw; $('dk').value = dk;
  update();
}
['kmin','kmax','kweight','dk'].forEach(id =>
  $(id).addEventListener('input', update));
$('applyRec').addEventListener('click', () => {
  const r = DATA.recommended;
  if (r) setFrame(r.kmin, r.kmax, DATA.current.kweight, DATA.current.dk);
});
$('reset').addEventListener('click', () =>
  setFrame(DATA.current.kmin, DATA.current.kmax,
           DATA.current.kweight, DATA.current.dk));
$('copyCmd').addEventListener('click', async () => {
  try { await navigator.clipboard.writeText($('cmd').textContent);
        $('copyCmd').textContent = '已复制'; }
  catch (e) { $('copyCmd').textContent = '复制失败，请手动选择'; }
  setTimeout(() => { $('copyCmd').textContent = '复制对应 CLI 命令'; }, 1500);
});
setFrame(DATA.current.kmin, DATA.current.kmax,
         DATA.current.kweight, DATA.current.dk);
</script>
</body>
</html>
"""
