# -*- coding: utf-8 -*-
"""
run_ci.py — 跨平台 CI 入口（编译 + 契约回归 + 可选压力测试）
============================================================
为什么需要它：
    回归测试直接跑 ``unittest`` 时，PowerShell 会把 stderr 当错误；
    GitHub Actions / 本地 / 服务器又都用不同的 shell。本脚本用
    ``subprocess`` 统一调用，按步骤收集输出与耗时，失败即非零退出。

用法:
    python run_ci.py               # 编译 + regression_test.py
    python run_ci.py --stress      # 额外跑 Co/Ni/Cu/NiO 压力测试
    python run_ci.py --json-only   # 不打印子进程输出（只写报告）

产物:
    ci_report.json   机器可读（步骤/耗时/退出码/输出尾部）
    ci_report.md     人读摘要（可直接贴 PR）

退出码: 0 全部通过；1 任一步骤失败。
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import platform
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PY = sys.executable
TAIL_CHARS = 6000


def _run(name: str, cmd: list[str], timeout: int = 3600,
         echo: bool = True) -> dict:
    start = _dt.datetime.now()
    print(f"\n===== {name} =====")
    print("$ " + " ".join(cmd))
    try:
        proc = subprocess.run(
            cmd, cwd=ROOT, capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=timeout)
        out = (proc.stdout or "") + (proc.stderr or "")
        code = proc.returncode
    except subprocess.TimeoutExpired as exc:
        out = ((exc.stdout or "") if isinstance(exc.stdout, str)
               else "") + f"\n[timeout] 超过 {timeout}s"
        code = 124
    except OSError as exc:
        out = f"[launch error] {exc}"
        code = 127
    elapsed = (_dt.datetime.now() - start).total_seconds()
    if echo and out.strip():
        print(out.strip())
    print(f"--> exit={code}  {elapsed:.1f}s")
    return {"name": name, "cmd": cmd, "exit_code": code,
            "seconds": round(elapsed, 2), "output_tail": out[-TAIL_CHARS:]}


def _compile_step() -> dict:
    """逐个 py_compile 根目录的 .py（避免 compileall 扫进无关目录）。"""
    files = sorted(p.name for p in ROOT.glob("*.py") if p.is_file())
    start = _dt.datetime.now()
    print(f"\n===== py_compile（{len(files)} 个文件）=====")
    bad = []
    messages = []
    for name in files:
        proc = subprocess.run(
            [PY, "-m", "py_compile", name], cwd=ROOT,
            capture_output=True, text=True, encoding="utf-8",
            errors="replace")
        if proc.returncode != 0:
            bad.append(name)
            messages.append(f"{name}:\n{proc.stderr or proc.stdout}")
    elapsed = (_dt.datetime.now() - start).total_seconds()
    code = 0 if not bad else 1
    out = "全部通过" if not bad else "失败: " + ", ".join(bad) + "\n" + \
        "\n".join(messages)
    print(out)
    return {"name": "py_compile", "cmd": [PY, "-m", "py_compile", "*.py"],
            "exit_code": code, "seconds": round(elapsed, 2),
            "output_tail": out[-TAIL_CHARS:]}


def _write_reports(steps: list[dict], seconds: float) -> tuple[Path, Path]:
    env = {
        "python": sys.version.split()[0],
        "executable": PY,
        "platform": platform.platform(),
        "machine": platform.machine(),
    }
    reg = next((s for s in steps if s["name"] == "regression_test.py"), None)
    match = (re.search(r"Ran\s+(\d+)\s+tests?", reg["output_tail"])
             if reg is not None else None)
    n_tests = int(match.group(1)) if match else None
    n_stress = sum(1 for s in steps if s["name"].startswith("stress_test_"))
    ok = all(s["exit_code"] == 0 for s in steps)
    report = {
        "generated_at": _dt.datetime.now().isoformat(timespec="seconds"),
        "ok": ok,
        "total_seconds": round(seconds, 2),
        "contract_tests": n_tests,
        "stress_tests": n_stress,
        "environment": env,
        "steps": steps,
    }
    jpath = ROOT / "ci_report.json"
    jpath.write_text(json.dumps(report, ensure_ascii=False, indent=2),
                     encoding="utf-8")

    lines = [
        "# XAFS CI 报告",
        "",
        f"- 结果: **{'PASS' if ok else 'FAIL'}**",
        f"- 时间: {report['generated_at']}",
        f"- 总耗时: {seconds:.1f}s",
        f"- 契约回归: **{n_tests if n_tests is not None else '未知'} 项**",
        f"- 压力测试: **{n_stress} 项**",
        f"- 环境: Python {env['python']} / {env['platform']}",
        "",
        "| 步骤 | 退出码 | 耗时 (s) |",
        "|------|--------|----------|",
    ]
    for s in steps:
        lines.append(f"| {s['name']} | {s['exit_code']} | {s['seconds']} |")
    if not ok:
        lines += ["", "## 失败输出", ""]
        for s in steps:
            if s["exit_code"] != 0:
                lines += [f"### {s['name']}", "", "```",
                          s["output_tail"][-3000:], "```", ""]
    mpath = ROOT / "ci_report.md"
    mpath.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return jpath, mpath


def main(argv=None) -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, OSError):
        pass
    ap = argparse.ArgumentParser(description="XAFS 管线 CI")
    ap.add_argument("--stress", action="store_true",
                    help="额外运行 Co/Ni/Cu/NiO 合成压力测试")
    ap.add_argument("--json-only", action="store_true",
                    help="不打印子进程完整输出（只写报告）")
    args = ap.parse_args(argv)
    echo = not args.json_only

    start = _dt.datetime.now()
    steps = [_compile_step()]
    steps.append(_run("regression_test.py",
                      [PY, "regression_test.py"], timeout=3600, echo=echo))
    if args.stress:
        for name in ("stress_test_co.py", "stress_test_ni.py",
                     "stress_test_cu.py", "stress_test_nio.py"):
            steps.append(_run(name, [PY, name], timeout=3600, echo=echo))

    seconds = (_dt.datetime.now() - start).total_seconds()
    jpath, mpath = _write_reports(steps, seconds)
    ok = all(s["exit_code"] == 0 for s in steps)
    print(f"\n{'=' * 50}")
    print(f"CI {'PASS' if ok else 'FAIL'}  ({seconds:.1f}s)")
    print(f"  → {jpath}")
    print(f"  → {mpath}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
