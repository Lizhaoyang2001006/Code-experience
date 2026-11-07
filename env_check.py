# -*- coding: utf-8 -*-
"""
env_check.py — 环境自检（SKILL.md §0 配套，agent 接手任何任务前必跑）
================================================================
用法（必须用 WorkBuddy 解释器；PATH 上的 python 是系统/Store stub，无依赖库）:

    & "C:\\Users\\lizha\\.workbuddy\\binaries\\python\\versions\\3.13.12\\python.exe" env_check.py

打印 Python 与各依赖库版本；缺任何一个库 → exit 1，禁止跑管线。
"""
import sys

REQUIRED = [
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("matplotlib", "matplotlib"),
    ("gemmi", "gemmi"),
    ("larch", "xraylarch"),
]


def _version(mod) -> str:
    v = getattr(mod, "__version__", None)
    if isinstance(v, str):
        return v
    try:
        from larch.version import __version__ as lv   # larch 特例
        return lv
    except Exception:
        return "?"


def main():
    print(f"python      {sys.version.split()[0]}   {sys.executable}")
    missing = []
    for name, pkg in REQUIRED:
        try:
            print(f"{pkg:<11} {_version(__import__(name))}")
        except Exception as exc:      # ImportError 及其连带失败都算缺库
            missing.append(pkg)
            print(f"{pkg:<11} 缺失 —— {type(exc).__name__}: {exc}")
    if missing:
        sys.exit(f"\n环境不完整，缺: {sorted(missing)}"
                 f"（larch/gemmi 缺失 = S1–S7 全瘫）")
    print("\n环境自检通过")


if __name__ == "__main__":
    main()
