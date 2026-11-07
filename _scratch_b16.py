# -*- coding: utf-8 -*-
import gzip
import re
import sys
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

p = Path(r"data/练手数据/02_Fe_standards/Fe_standards.prj")
text = gzip.decompress(p.read_bytes()).decode("utf-8", errors="replace")
print("total len", len(text))
print("record markers", text.count("[record]"))
print("args markers", text.count("@args ="))
print("x markers", text.count("@x ="), "y", text.count("@y ="))

blocks = re.split(r"(?m)^\[record\].*$", text)
print("blocks", len(blocks), [len(b) for b in blocks][:12])


def arg_string(block):
    m = re.search(r"@args\s*=\s*\(", block)
    if not m:
        return None
    i = m.end() - 1
    depth = 0
    quote = None
    esc = False
    j = i
    while j < len(block):
        ch = block[j]
        if quote:
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == quote:
                quote = None
        else:
            if ch in "'\"":
                quote = ch
            elif ch in "([{":
                depth += 1
            elif ch in ")]}":
                depth -= 1
                if depth == 0:
                    return block[i + 1:j]
        j += 1
    return None


def split_top(s):
    out, buf, depth, quote, esc = [], [], 0, None, False
    for ch in s:
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


def unquote(tok):
    tok = tok.strip()
    if len(tok) >= 2 and tok[0] == tok[-1] and tok[0] in "'\"":
        body = tok[1:-1]
        body = body.replace("\\" + tok[0], tok[0]).replace("\\\\", "\\")
        return body
    return tok


def get_array(block, name):
    m = re.search(r"@" + re.escape(name) + r"\s*=\s*\(", block)
    if not m:
        return None
    i = m.end() - 1
    depth = 0
    j = i
    while j < len(block):
        ch = block[j]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                body = block[i + 1:j]
                vals = []
                for tok in split_top(body):
                    try:
                        vals.append(float(tok.strip("'\"")))
                    except ValueError:
                        pass
                return vals
        j += 1
    return None


for bi, block in enumerate(blocks):
    if "@args" not in block:
        continue
    body = arg_string(block)
    toks = split_top(body) if body else []
    args = {}
    for k in range(0, len(toks) - 1, 2):
        key = unquote(toks[k])
        val = unquote(toks[k + 1])
        args[key] = val
    label = args.get("label") or args.get("tag") or args.get("old_group")
    x = get_array(block, "x")
    y = get_array(block, "y")
    i0 = get_array(block, "i0")
    sig = get_array(block, "signal")
    yneg = sum(1 for v in (y or []) if v < 0)
    print(bi, repr(label), "tag=", repr(args.get("tag")),
          "merge=", repr(args.get("is_merge")),
          "nor=", repr(args.get("is_nor")), "ln=", repr(args.get("ln")),
          "z=", repr(args.get("bkg_z")), "e0=", repr(args.get("bkg_e0")),
          "len x/y/i0/sig=", [len(v) if v is not None else None
                               for v in (x, y, i0, sig)],
          "yneg=", yneg)

# ---- mu 还原公式对比 ----
import numpy as np
sys.path.insert(0, str(Path.cwd()))
import xafs_io

for bi in (0, 2, 4, 8):
    block = blocks[bi]
    body = arg_string(block)
    toks = split_top(body)
    args = {}
    for k in range(0, len(toks) - 1, 2):
        args[unquote(toks[k])] = unquote(toks[k + 1])
    x = np.asarray(get_array(block, "x"))
    y = np.asarray(get_array(block, "y"))
    i0 = np.asarray(get_array(block, "i0"))
    sig = np.asarray(get_array(block, "signal"))
    mu_log = np.log(i0 / sig)
    label = args.get("label")
    raw = xafs_io.read_scan(
        Path("data/练手数据/02_Fe_standards") / label.strip())
    print(label, "y vs log maxdiff", float(np.max(np.abs(y - mu_log))),
          "corr y/raw", float(np.corrcoef(y, raw.mu)[0, 1]),
          "corr log/raw", float(np.corrcoef(mu_log, raw.mu)[0, 1]),
          "maxdiff log/raw", float(np.max(np.abs(mu_log - raw.mu))),
          "energy diff", float(np.max(np.abs(x - raw.energy))))

