# -*- coding: utf-8 -*-
"""
xafs_io.py — 数据进出层
================================================================
职责边界：**只做格式转换，不做任何物理计算**。
  读：XDAC / XDI / 纯文本  →  RawScan
  写：ndarray → 文本表 / JSON

本模块不 import athena / artemis / larch_wt，因此可被任意阶段复用。
μ(E) 的换算也放在这里，因为它是「通道角色」问题，不是「谱学处理」问题：
透射取对数、荧光直接比值，选哪条由文件里有哪些通道决定。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from xafs_core import RawScan, XAFSConfigError, XAFSDataError

# ================================================================
# 1. 格式识别
# ================================================================

def sniff_format(path) -> str:
    """
    按首行特征判断格式，不做推测性解析。

    返回: "XDAC" | "XDI" | "SIMPLE"
    """
    path = Path(path)
    with path.open("rb") as fh:
        head = fh.read(2048).decode("latin-1", errors="replace")
    first = head.splitlines()[0] if head else ""
    if first.upper().startswith("XDAC"):
        return "XDAC"
    if first.lstrip().upper().startswith("# XDI"):
        return "XDI"
    if "# Column." in head:
        return "XDI"
    return "SIMPLE"


def _safe_float(value):
    """容忍 '7112.0 eV' 这类带单位的字段；解析失败返回 None 而不是崩。"""
    if value is None:
        return None
    try:
        return float(str(value).split()[0]) or None
    except (ValueError, IndexError):
        return None


def _to_float(token: str):
    """'0.05k' / '15k' / '2.80' → float；无法转换返回 None。"""
    try:
        return float(token.strip().rstrip("kK"))
    except ValueError:
        return None


# ================================================================
# 2. 通道角色校验
# ================================================================

def assert_channel_roles(i0, it, energy, name=""):
    """
    物理合理性检查：挡住「把 I0 当成 μ」这类静默错误。

    判据:
        I0/It 中位数 ∈ (1.01, 1e5)   —— 入射必须大于透射
        μ 中位数 ∈ (0.05, 20)        —— 吸收系数不可能跨几个数量级

    这两条来自实测：把 I0 当 μ 时，E0 只偏 1.25 eV 看不出来，
    但 edge_step 会从 0.47 变成 979.9。
    """
    if np.any(i0 <= 0) or np.any(it <= 0):
        raise XAFSDataError("S1", f"{name}: I0/It 含非正值，列角色可能判错")
    ratio = float(np.median(i0 / it))
    if not (1.01 < ratio < 1e5):
        raise XAFSDataError(
            "S1", f"{name}: I0/It 中位数 = {ratio:.4g}，不像入射/透射比。"
                  f"检查列顺序（XDAC 为 Energy,I0,It,Ir,IntTime）")
    mu = np.log(i0 / it)
    if not (0.05 < float(np.median(mu)) < 20):
        raise XAFSDataError(
            "S1", f"{name}: μ 中位数 = {np.median(mu):.4g}，超出合理范围 0.05–20")


# ================================================================
# 3. 各格式解析器
# ================================================================

def read_xdac(path) -> RawScan:
    """
    XDAC V1.x（SSRL / NSLS 老格式）。

    结构: 14 行头部 + 分隔线 + 数值表
    坑位: ① CRLF 行尾（splitlines 可解）
          ② SRB/SRSS 里有 '15k'、'0.05k'，float() 会抛异常
          ③ 列表固定为 Energy, I0, It, Ir, IntTime
    """
    path = Path(path)
    lines = path.read_text(encoding="latin-1", errors="replace").splitlines()
    header, rows = {}, []
    for line in lines:
        s = line.strip()
        if not s or s.startswith("-") or s.lower().startswith("energy"):
            continue
        tokens = s.split()
        values = [_to_float(t) for t in tokens]
        # ≥3 列即可（E, I0, It 为最小组）；无参考通道的 XDAC（如 01_Fe_foil/fe.*）
        # 只有 3 列，旧门槛 >=4 会把数据行全部丢掉。
        if all(v is not None for v in values) and len(values) >= 3:
            rows.append(values)
        elif "=" in s:
            key, val = s.split("=", 1)
            header[key.strip()] = val.strip()

    arr = np.asarray(rows, dtype=float)
    if arr.size == 0:
        raise XAFSDataError("S1", f"{path.name}: 未解析到数值行")
    if arr.shape[1] < 3:
        raise XAFSDataError("S1", f"{path.name}: 仅 {arr.shape[1]} 列，无法区分 I0/It")

    energy = arr[:, 0]
    i0, it = arr[:, 1], arr[:, 2]
    iref = arr[:, 3] if arr.shape[1] >= 4 else None
    dtime = arr[:, 4] if arr.shape[1] >= 5 else None
    assert_channel_roles(i0, it, energy, path.name)

    line2 = lines[1] if len(lines) > 1 else ""
    meta = {
        "format": "XDAC",
        "mode": "trans",
        "path": str(path),
        "n_points": len(energy),
        "beamline": line2.rsplit(" on ", 1)[-1].strip() if " on " in line2 else None,
        "e0_nominal": _to_float(header.get("E0", "")),
        "header": {k: v for k, v in header.items()
                   if k in ("E0", "NUM_REGIONS", "SRB", "SRSS", "Offsets", "Gains")},
        "channels": ["Energy", "I0", "It"] + (["Ir"] if iref is not None else []),
    }
    return RawScan(energy=energy, i0=i0, it=it, iref=iref, dtime=dtime,
                   mu=np.log(i0 / it), meta=meta)


def read_xdi(path) -> RawScan:
    """
    XDI 1.0（NSLS-II / APS / ESRF 现代输出）。

    结构: '#' 注释 + '# Column.N: name unit' 声明 + 数值表
    坑位: **列名行本身也被 '#' 注释**，按「跳过所有 # 行」实现就读不到列名。

    通道角色按优先级判定：
        xmu 列存在        → mode="given"
        I0 与 It 都在     → mode="trans"
        I0 与 If/FF 都在  → mode="fluor"
    """
    path = Path(path)
    kv, colmap, colnames, rows = {}, {}, [], []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            body = s.lstrip("# ").strip()
            if body.startswith("Column.") and ":" in body:
                key, val = body.split(":", 1)
                parts = val.strip().split()
                if parts:
                    colmap[parts[0].lower()] = int(key.split(".")[1]) - 1
            elif ":" in body:
                key, val = body.split(":", 1)
                kv[key.strip()] = val.strip()
            elif " " in body:
                cand = body.split()
                if any(c.lower() in ("energy", "xmu") for c in cand):
                    colnames = [c.lower() for c in cand]
            continue
        try:
            rows.append([float(t) for t in s.split()])
        except ValueError:
            continue

    arr = np.asarray(rows, dtype=float)
    if arr.size == 0:
        raise XAFSDataError("S1", f"{path.name}: 未解析到数值行")
    names = colnames or sorted(colmap, key=colmap.get)

    def pick(*candidates):
        for name in candidates:
            if name in names:
                return arr[:, names.index(name)]
        return None

    energy = pick("energy")
    if energy is None:
        raise XAFSDataError("S1", f"{path.name}: 找不到 energy 列")

    # 列名兼容 XDI 规范的两套写法：
    #   Demeter 系 it/if/ir/xmu ；larch/BMM 系 itrans/ifluor/irefer/mutrans。
    #   BMM 标准库(mds2-4032)中透射谱为 energy,...,mutrans,i0,itrans,irefer，
    #   缺这些别名时 70/359 个谱全部误判为「无 xmu，无法配对通道」。
    i0 = pick("i0")
    it = pick("it", "itrans")
    iref = pick("ir", "irefer")
    ifluor = pick("if", "ff", "i1", "ifluor", "fluor", "fluorescence")
    mu_given = pick("xmu", "mu", "mutrans", "mufluor")

    # ---- B9-3 多通道荧光探测器（MCA / 多 If 通道）----
    def _numbered(prefix):
        pairs = []
        for nm in names:
            tail = nm[len(prefix):] if nm.startswith(prefix) else ""
            if tail.isdigit():
                pairs.append((int(tail), nm))
        return [nm for _, nm in sorted(pairs)]

    mca_names = _numbered("mca") or _numbered("if")
    mca_channels = dtime_channels = None
    if len(mca_names) >= 2:
        mca_channels = np.column_stack([arr[:, names.index(nm)]
                                        for nm in mca_names])
        dt_names = _numbered("dt")
        if len(dt_names) == len(mca_names):
            dtime_channels = np.column_stack([arr[:, names.index(nm)]
                                              for nm in dt_names])
        ifluor = np.sum(mca_channels, axis=1)
        # 通道求和只是"物理总荧光"的第一步；坏道剔除/死时间校正由
        # athena.merge_mca_channels 在 S1.2 完成。此处不得提前平均，
        # 否则坏道的 Bragg spike 会被稀释到无法识别。
        if i0 is not None and mu_given is None:
            mode = "fluor"
            mu = ifluor / i0
        else:
            mode = "given" if mu_given is not None else "fluor"
            mu = mu_given if mu_given is not None else ifluor

    if mca_channels is None:
        if mu_given is not None:
            mode, mu = "given", mu_given
        elif i0 is not None and it is not None:
            mode = "trans"
            assert_channel_roles(i0, it, energy, path.name)
            mu = np.log(i0 / it)
        elif i0 is not None and ifluor is not None:
            mode = "fluor"
            if np.any(ifluor <= 0):
                raise XAFSDataError("S1", f"{path.name}: 荧光通道 If 含非正值")
            mu = ifluor / i0            # 荧光模式绝不取对数
        else:
            raise XAFSDataError("S1", f"{path.name}: 无 xmu，也无法配对 I0/It 或 I0/If")

    meta = {
        "format": "XDI",
        "mode": mode,
        "path": str(path),
        "n_points": len(energy),
        "beamline": kv.get("Beamline.name"),
        "sample": kv.get("Sample.name"),
        "element": kv.get("Element.symbol"),
        "e0_nominal": _safe_float(kv.get("Scan.edge_energy")),
        "channels": names,
    }
    if mca_channels is not None:
        meta["mca"] = {
            "channels": mca_channels,
            "dtime": dtime_channels,
            "channel_names": mca_names,
        }
        meta["mode"] = "fluor_multichannel"
        meta["channel_sum_note"] = ("mu 为各通道求和（坏道剔除由 "
                                    "athena.merge_mca_channels 执行）")
    if mu_given is not None and i0 is not None and it is not None:
        meta["xmu_selfcheck"] = float(np.max(np.abs(np.log(i0 / it) - mu_given)))
    return RawScan(energy=energy, i0=i0, it=it, ifluor=ifluor, iref=iref,
                   mu=mu, meta=meta)


def read_simple(path) -> RawScan:
    """
    无格式声明的纯文本。
        2 列 → 视作 (E, μ) 已算好（查中位数量级，挡住把 I0 当 μ）
        3~4 列 → 视作 (E, I0, It[, Ir])，按透射处理
    """
    path = Path(path)
    rows, header_names = [], []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            body = s.lstrip("# ").strip()
            if body and not any(ch.isdigit() for ch in body[:3]):
                cand = body.replace(",", " ").split()
                if len(cand) >= 2:
                    header_names = [c.lower() for c in cand]
            continue
        try:
            rows.append([float(t) for t in s.replace(",", " ").split()])
        except ValueError:
            continue

    arr = np.asarray(rows, dtype=float)
    ncol = arr.shape[1] if arr.ndim == 2 else 0
    base = {"path": str(path), "n_points": len(arr), "channels": header_names}

    if ncol == 2:
        mu = arr[:, 1]
        mu_median = float(np.median(mu))
        # 量级检查：μ(E)（无论透射对数还是荧光比值）中位数不可能到几百以上；
        # 若给的是 I0/It 等原始光强（典型 1e4–1e6），在此挡住而不是静默当 μ 用。
        if not (1e-3 < mu_median < 50.0):
            raise XAFSDataError(
                "S1", f"{path.name}: 第 2 列中位数 = {mu_median:.4g}，不像 μ(E)"
                      f"（合理范围约 0.001–50）。若该列是 I0/It 等原始光强，"
                      f"请改用 3 列格式 (E, I0, It)，由程序取对数算 μ")
        return RawScan(energy=arr[:, 0], mu=mu,
                       meta={**base, "format": "SIMPLE_2COL", "mode": "given"})
    if ncol >= 3:
        assert_channel_roles(arr[:, 1], arr[:, 2], arr[:, 0], path.name)
        return RawScan(energy=arr[:, 0], i0=arr[:, 1], it=arr[:, 2],
                       iref=arr[:, 3] if ncol > 3 else None,
                       mu=np.log(arr[:, 1] / arr[:, 2]),
                       meta={**base, "format": "SIMPLE_3COL", "mode": "trans"})
    raise XAFSDataError("S1", f"{path.name}: 无法识别的列数 {ncol}")


# ================================================================
# 4. 统一入口
# ================================================================

_PARSERS = {"XDAC": read_xdac, "XDI": read_xdi, "SIMPLE": read_simple}


def read_scan(path) -> RawScan:
    """按格式分发解析，返回 RawScan。"""
    path = Path(path)
    fmt = sniff_format(path)
    scan = _PARSERS[fmt](path)
    scan.meta.setdefault("format_detected", fmt)
    return scan


def read_scans(paths) -> list:
    """批量读取；单文件也返回长度为 1 的列表。"""
    return [read_scan(p) for p in paths]


# ================================================================
# 5. 写出
# ================================================================

def write_columns(path, header: str, columns, column_names) -> Path:
    """
    写多列文本表。

    columns     : 可迭代的 ndarray，长度必须一致
    column_names: 列名列表，写进表头
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = [np.asarray(c, dtype=float) for c in columns]
    n = len(arrays[0])
    if any(len(a) != n for a in arrays):
        raise XAFSConfigError("write", f"{path.name}: 各列长度不一致")

    with open(path, "w", encoding="utf-8") as f:
        f.write(header if header.endswith("\n") else header + "\n")
        f.write("# " + "  ".join(column_names) + "\n")
        for i in range(n):
            f.write("  ".join(f"{a[i]:.6e}" for a in arrays) + "\n")
    return path


def write_wide_csv(path, index, index_name: str, series: dict) -> Path:
    """
    宽表：一列自变量 + 每条谱一列，可直接导入 Origin。
    series 形如 {"scan1": arr1, "merge": arr2}
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(series)
    with open(path, "w", encoding="utf-8-sig") as f:
        f.write(",".join([index_name] + names) + "\n")
        for i in range(len(index)):
            f.write(",".join([f"{index[i]:.6f}"] + [f"{series[n][i]:.6e}" for n in names]) + "\n")
    return path


def write_json(path, payload) -> Path:
    """把阶段状态落盘，便于断点续跑与复现。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    def _default(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()
        if isinstance(obj, Path):
            return str(obj)
        raise TypeError(f"无法序列化 {type(obj)}")

    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False,
                               default=_default), encoding="utf-8")
    return path


# ================================================================
# 5.5 `.prj` 工程容器（B9-3）
# ================================================================

_PRJ_SCAN_FIELDS = ("energy", "i0", "it", "ifluor", "iref", "dtime",
                    "mu", "sigma_mu")


def _prj_default(obj):
    """gzip-JSON 容器的序列化器：numpy / Path / dataclass 都落成纯 JSON。"""
    import dataclasses

    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, Path):
        return str(obj)
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return {k: v for k, v in dataclasses.asdict(obj).items()}
    raise TypeError(f"无法序列化 {type(obj)}")


def write_prj(path, payload) -> Path:
    """写 gzip 压缩的 JSON 工程容器（`.prj`）。

    容器结构::

        {"xafs_prj": 1, "scans": [...], "standards": {...}, "meta": {...}}

    ``payload`` 可为任意 JSON 可序列化字典；``scans`` 里的 RawScan
    dataclass 会被自动展开。使用 gzip 是为了避免大数组把工程文件
    膨胀到几百 MB（实测 3 条 1000 点扫描 JSON≈70 KB，gzip≈18 KB）。
    注意：这是本项目的容器，**不是** Athena/Demeter 的原生 Perl `.prj`；
    后者请先导出为文本/JSON。
    """
    import gzip

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not isinstance(payload, dict):
        raise XAFSConfigError("prj", "write_prj 的 payload 必须是字典")
    doc = {"xafs_prj": 1}
    doc.update(payload)
    raw = json.dumps(doc, ensure_ascii=False, default=_prj_default,
                     separators=(",", ":")).encode("utf-8")
    with gzip.open(path, "wb") as fh:
        fh.write(raw)
    return path


def _athena_split_top(text: str) -> list[str]:
    """按顶层逗号切分 Perl 列表，跳过引号与嵌套 ``{}`` / ``[]`` / ``()``。"""
    out: list[str] = []
    buf: list[str] = []
    depth = 0
    quote = None
    esc = False
    for ch in text:
        if quote:
            buf.append(ch)
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == quote:
                quote = None
            continue
        if ch in "'\"":
            quote = ch
            buf.append(ch)
        elif ch in "([{":
            depth += 1
            buf.append(ch)
        elif ch in ")]}":
            depth -= 1
            buf.append(ch)
        elif ch == "," and depth == 0:
            out.append("".join(buf).strip())
            buf = []
        else:
            buf.append(ch)
    if buf:
        out.append("".join(buf).strip())
    return out


def _athena_unquote(token: str) -> str:
    """去掉 Perl 单/双引号并还原最小转义；非字符串 token 原样返回。"""
    token = token.strip()
    if len(token) >= 2 and token[0] == token[-1] and token[0] in "'\"":
        quote = token[0]
        body = token[1:-1]
        return body.replace("\\" + quote, quote).replace("\\\\", "\\")
    return token


def _athena_balanced(text: str, start: int) -> str | None:
    """从 ``@args = (`` / ``@x = (`` 的左括号起取平衡括号内容。"""
    depth = 0
    quote = None
    esc = False
    for i in range(start, len(text)):
        ch = text[i]
        if quote:
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == quote:
                quote = None
            continue
        if ch in "'\"":
            quote = ch
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return text[start + 1:i]
    return None


def _athena_args(block: str) -> dict:
    """解析 ``@args = (key, value, ...)``；嵌套结构整体保留为字符串值。"""
    import re

    match = re.search(r"@args\s*=\s*\(", block)
    if match is None:
        return {}
    body = _athena_balanced(block, match.end() - 1)
    if body is None:
        return {}
    tokens = _athena_split_top(body)
    args: dict = {}
    for i in range(0, len(tokens) - 1, 2):
        args[_athena_unquote(tokens[i])] = _athena_unquote(tokens[i + 1])
    return args


def _athena_array(block: str, name: str):
    """解析 ``@x = (...)`` 浮点数组；缺失返回 None，坏 token 记 NaN。"""
    import re

    match = re.search(r"@" + re.escape(name) + r"\s*=\s*\(", block)
    if match is None:
        return None
    body = _athena_balanced(block, match.end() - 1)
    if body is None:
        return None
    values = []
    for token in _athena_split_top(body):
        try:
            values.append(float(_athena_unquote(token)))
        except ValueError:
            values.append(np.nan)
    if not values:
        return None
    return np.asarray(values, dtype=float)


def read_athena_perl_prj(path_or_text, as_scans: bool = False,
                         include_derived: bool = False) -> dict:
    """B16-1：只读解析 Athena 原生 Perl ``.prj`` 工程文件。

    判定为合法 Athena 工程的充要条件是文本同时含 ``@args =`` 与至少一组
    非空的 ``@x =`` / ``@y =`` 数组；否则抛 ``XAFSConfigError``（保护
    既有 #50：任意 ``$VAR1`` Perl 片段不得被猜成扫描）。

    记录以行首 ``[record]`` 分隔（后面常带 ``# create`` 注释）。每条记录：
      - ``@args`` 偶数下标为键、奇数下标为值（嵌套 ``{}``/``[]`` 跳过）；
      - ``label`` 缺省时退到 ``tag`` / ``old_group``；
      - 默认过滤参比通道派生记录（``label`` 去空格以 ``Ref `` 开头或
        ``tag`` 以 ``d_Ref_`` 开头）与合并派生记录（``is_merge∉{"", "0"}``）；
      - ``ln=1`` 且有等长非零 ``@i0``/``@signal`` 时按 μ=ln(I₀/I_signal)
        还原（实测 Ref 记录 ``@y`` 为负值）；否则直接取 ``@y``。
    """
    import gzip
    import re

    if isinstance(path_or_text, bytes):
        blob = path_or_text
        source = "<bytes>"
    elif isinstance(path_or_text, Path):
        if not path_or_text.is_file():
            raise XAFSConfigError("prj", f"工程文件不存在: {path_or_text}")
        blob = path_or_text.read_bytes()
        source = str(path_or_text)
    elif isinstance(path_or_text, str):
        if "\n" in path_or_text or "\r" in path_or_text:
            blob = path_or_text.encode("utf-8")
            source = "<text>"
        else:
            p = Path(path_or_text)
            if p.is_file():
                blob = p.read_bytes()
                source = str(p)
            else:
                blob = path_or_text.encode("utf-8")
                source = "<text>"
    else:
        raise XAFSConfigError(
            "prj", f"read_athena_perl_prj 不支持 {type(path_or_text)}")
    if blob[:2] == b"\x1f\x8b":
        try:
            blob = gzip.decompress(blob)
        except OSError as exc:
            raise XAFSConfigError("prj", f"{source}: gzip 解压失败: {exc}")
    try:
        text = blob.decode("utf-8")
    except UnicodeDecodeError:
        text = blob.decode("latin-1")

    if not (re.search(r"@args\s*=", text)
            and re.search(r"@x\s*=", text)
            and re.search(r"@y\s*=", text)):
        raise XAFSConfigError(
            "prj", f"{source}: 检测到 Perl 文本但缺少 @args/@x/@y 扫描记录，"
                   f"不是可识别的 Athena 原生工程（拒绝猜测格式）")

    blocks = re.split(r"(?m)^\[record\].*$", text)
    records = []
    for idx, block in enumerate(blocks):
        if "@args" not in block:
            continue
        args = _athena_args(block)
        x = _athena_array(block, "x")
        y = _athena_array(block, "y")
        if x is None or y is None or x.size == 0 or y.size == 0:
            continue
        label = str(args.get("label") or args.get("tag")
                    or args.get("old_group") or f"record_{idx}")
        tag = str(args.get("tag", ""))
        is_ref = (label.strip().startswith("Ref ")
                  or tag.startswith("d_Ref_"))
        is_merge = str(args.get("is_merge", "")) not in ("", "0")
        records.append({"index": idx, "block": block, "args": args,
                        "x": x, "y": y, "label": label,
                        "is_ref": is_ref, "is_merge": is_merge})

    if not records:
        raise XAFSConfigError(
            "prj", f"{source}: @args/@x/@y 标记存在但没有非空扫描记录")

    scans = []
    for rec in records:
        if not include_derived and (rec["is_ref"] or rec["is_merge"]):
            continue
        args = rec["args"]
        x, y = rec["x"], rec["y"]
        if x.size != y.size:
            raise XAFSConfigError(
                "prj", f"{source}: 记录 {rec['label']!r} 的 @x/@y 长度"
                       f"不一致（{x.size} vs {y.size}）")
        i0 = _athena_array(rec["block"], "i0")
        signal = _athena_array(rec["block"], "signal")
        if i0 is not None and i0.size != x.size:
            i0 = None
        if signal is not None and signal.size != x.size:
            signal = None
        mu = y
        if (i0 is not None and signal is not None
                and np.any(i0 != 0) and np.any(signal != 0)
                and str(args.get("ln", "1")) == "1"
                and str(args.get("is_nor", "0")) != "1"):
            ratio = np.divide(i0, signal,
                              out=np.full_like(i0, np.nan),
                              where=signal != 0)
            if np.all(np.isfinite(ratio)) and np.all(ratio > 0):
                mu = np.log(ratio)
        e0_value = args.get("bkg_e0")
        try:
            e0_float = float(e0_value) if e0_value is not None else None
        except (TypeError, ValueError):
            e0_float = None
        meta = {
            "format_detected": "athena_perl",
            "source_format": "athena_perl",
            "label": rec["label"],
            "tag": str(args.get("tag", "")),
            "absorber": str(args.get("bkg_z", "")),
            "e0": e0_float,
            "is_merge": str(args.get("is_merge", "")),
            "is_nor": str(args.get("is_nor", "")),
            "ln": str(args.get("ln", "")),
            "source": str(args.get("source", "")),
            "record_index": rec["index"],
        }
        scans.append(RawScan(energy=x, i0=i0, it=signal, mu=mu, meta=meta))

    doc = {
        "xafs_prj": "athena_perl_v1",
        "scans": ([_prj_default(s) for s in scans] if not as_scans else scans),
        "standards": {},
        "meta": {"source_format": "athena_perl",
                 "n_records": len(records),
                 "n_scans": len(scans),
                 "source": source},
    }
    return doc


def read_prj(path, as_scans: bool = False,
             include_derived: bool = False) -> dict:
    """读 `.prj` 工程容器；兼容未压缩 JSON 与 Athena 原生 Perl 文本。

    JSON 通道保持 B9-3 契约；非 JSON 文本只有在同时含 ``@args =``、
    ``@x =``、``@y =`` 且至少一条非空扫描时才转发给
    ``read_athena_perl_prj``，否则抛 ``XAFSConfigError``（保护 #50）。
    """
    import gzip

    path = Path(path)
    if not path.is_file():
        raise XAFSConfigError("prj", f"工程文件不存在: {path}")
    blob = path.read_bytes()
    if blob[:2] == b"\x1f\x8b":
        try:
            blob = gzip.decompress(blob)
        except OSError as exc:
            raise XAFSConfigError("prj", f"{path.name}: gzip 解压失败: {exc}")
    try:
        text = blob.decode("utf-8")
    except UnicodeDecodeError:
        text = blob.decode("latin-1")
    try:
        doc = json.loads(text)
    except json.JSONDecodeError as exc:
        import re
        if (re.search(r"@args\s*=", text)
                and re.search(r"@x\s*=", text)
                and re.search(r"@y\s*=", text)):
            return read_athena_perl_prj(
                text, as_scans=as_scans, include_derived=include_derived)
        raise XAFSConfigError(
            "prj", f"{path.name}: 不是有效 JSON（{exc.msg}），也不是含 "
                   f"@args/@x/@y 扫描记录的 Athena 原生工程。"
                   f"若为其他 Perl 序列化，请先在 Athena 导出为文本/JSON")
    if not isinstance(doc, dict) or doc.get("xafs_prj") is None:
        raise XAFSConfigError(
            "prj", f"{path.name}: 缺少 xafs_prj 版本标记，拒绝猜测格式")
    if as_scans and isinstance(doc.get("scans"), list):
        scans = []
        for item in doc["scans"]:
            if not isinstance(item, dict) or "energy" not in item:
                continue
            kwargs = {}
            for key in _PRJ_SCAN_FIELDS:
                if item.get(key) is not None:
                    kwargs[key] = np.asarray(item[key], dtype=float)
            kwargs["meta"] = dict(item.get("meta") or {})
            scans.append(RawScan(**kwargs))
        doc["scans"] = scans
    return doc


# ================================================================
# 6. 自测
# ================================================================

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        target = sys.argv[1]
    else:
        target = (Path(__file__).resolve().parent / "data" / "练手数据" /
                  "02_Fe_standards" / "Fe_lepidocrocite.000")
    scan = read_scan(target)
    print(f"文件      : {Path(target).name}")
    print(f"格式      : {scan.meta['format_detected']}  mode={scan.meta['mode']}")
    print(f"点数      : {len(scan)}")
    print(f"能量范围  : {scan.energy.min():.1f} – {scan.energy.max():.1f} eV")
    print(f"μ 范围    : {scan.mu.min():.4f} – {scan.mu.max():.4f}")
