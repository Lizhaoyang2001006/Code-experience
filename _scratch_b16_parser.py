# -*- coding: utf-8 -*-
"""Prototype of the Athena Perl .prj reader (B16-1) for offline validation."""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def _read_text(path_or_text):
    if isinstance(path_or_text, bytes):
        blob = path_or_text
    elif isinstance(path_or_text, Path):
        blob = path_or_text.read_bytes()
    elif isinstance(path_or_text, str):
        if "\n" in path_or_text or "\r" in path_or_text:
            return path_or_text
        p = Path(path_or_text)
        if not p.is_file():
            raise FileNotFoundError(path_or_text)
        blob = p.read_bytes()
    else:
        raise TypeError(type(path_or_text))
    if blob[:2] == b"\x1f\x8b":
        blob = gzip.decompress(blob)
    try:
        return blob.decode("utf-8")
    except UnicodeDecodeError:
        return blob.decode("latin-1")


def _split_perl_top(s):
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


def _unquote(token):
    token = token.strip()
    if len(token) >= 2 and token[0] == token[-1] and token[0] in "'\"":
        q = token[0]
        body = token[1:-1]
        body = body.replace("\\" + q, q).replace("\\\\", "\\")
        return body
    return token


def _balanced_body(text, match):
    i = match.end() - 1
    depth = 0
    quote = None
    esc = False
    for j in range(i, len(text)):
        ch = text[j]
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
                return text[i + 1:j]
    return None


def _args(block):
    m = re.search(r"@args\s*=\s*\(", block)
    if not m:
        return {}
    body = _balanced_body(block, m)
    toks = _split_perl_top(body or "")
    out = {}
    for k in range(0, len(toks) - 1, 2):
        out[_unquote(toks[k])] = _unquote(toks[k + 1])
    return out


def _floats(block, name):
    m = re.search(r"@" + re.escape(name) + r"\s*=\s*\(", block)
    if not m:
        return None
    body = _balanced_body(block, m)
    if body is None:
        return None
    vals = []
    for tok in _split_perl_top(body):
        try:
            vals.append(float(_unquote(tok)))
        except ValueError:
            vals.append(np.nan)
    return np.asarray(vals, dtype=float) if vals else None


def read_athena_perl_prj(path_or_text, as_scans=False, include_derived=False):
    text = _read_text(path_or_text)
    if not (re.search(r"@args\s*=", text)
            and re.search(r"@x\s*=", text)
            and re.search(r"@y\s*=", text)):
        raise ValueError("not an Athena prj")
    blocks = re.split(r"(?m)^\[record\].*$", text)
    records = []
    for idx, block in enumerate(blocks):
        if "@args" not in block:
            continue
        args = _args(block)
        x = _floats(block, "x")
        y = _floats(block, "y")
        if x is None or y is None or x.size == 0 or y.size == 0:
            continue
        label = (args.get("label") or args.get("tag")
                 or args.get("old_group") or f"record_{idx}")
        label = str(label)
        tag = str(args.get("tag", ""))
        is_ref = label.strip().startswith("Ref ") or tag.startswith("d_Ref_")
        is_merge = str(args.get("is_merge", "")) not in ("", "0")
        records.append({"block": block, "args": args, "x": x, "y": y,
                        "label": label, "is_ref": is_ref,
                        "is_merge": is_merge})
    kept = [r for r in records
            if include_derived or not (r["is_ref"] or r["is_merge"])]
    out = []
    for rec in kept:
        args = rec["args"]
        x, y = rec["x"], rec["y"]
        i0 = _floats(rec["block"], "i0")
        signal = _floats(rec["block"], "signal")
        mu = y
        it = signal
        if (i0 is not None and signal is not None
                and i0.size == x.size and signal.size == x.size
                and np.any(i0 != 0) and np.any(signal != 0)
                and str(args.get("ln", "1")) == "1"
                and str(args.get("is_nor", "0")) != "1"):
            ratio = np.divide(i0, signal, out=np.full_like(i0, np.nan),
                             where=signal != 0)
            if np.all(np.isfinite(ratio)) and np.all(ratio > 0):
                mu = np.log(ratio)
            else:
                it = None
        else:
            if i0 is None or i0.size != x.size:
                i0 = None
            if signal is None or signal.size != x.size:
                it = None
        if mu.size != x.size:
            raise ValueError(f"label={rec['label']} y/x 长度不一致")
        out.append({"label": rec["label"], "args": args, "energy": x,
                    "mu": np.asarray(mu, float), "i0": i0, "it": it})
    return {"xafs_prj": "athena_perl_v1", "scans": out,
            "standards": {},
            "meta": {"source_format": "athena_perl",
                     "n_records": len(records), "n_scans": len(out)}}


if __name__ == "__main__":
    p = Path(r"data/练手数据/02_Fe_standards/Fe_standards.prj")
    out = read_athena_perl_prj(p)
    print("n_records/n_scans:", out["meta"])
    print([s["label"] for s in out["scans"]])
    allout = read_athena_perl_prj(p, include_derived=True)
    print("include_derived:", len(allout["scans"]))
    sys.path.insert(0, str(Path.cwd()))
    import xafs_io
    for s in out["scans"]:
        if s["label"] in ("Fe_lepidocrocite.000", "Fe_lepidocrocite.100",
                          "Fe_lepidocrocite.200", "Fe_pyrite.000"):
            raw = xafs_io.read_scan(
                Path("data/练手数据/02_Fe_standards") / s["label"])
            print(s["label"], len(s["energy"]), len(raw.energy),
                  float(np.max(np.abs(s["energy"] - raw.energy))),
                  float(np.corrcoef(s["mu"], raw.mu)[0, 1]),
                  float(np.max(np.abs(s["mu"] - raw.mu))))
    try:
        read_athena_perl_prj("$VAR1 = { 'x' => 1 };")
        print("BAD: $VAR1 accepted")
    except ValueError:
        print("$VAR1 rejected OK")
