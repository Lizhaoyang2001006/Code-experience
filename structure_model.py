# -*- coding: utf-8 -*-
"""structure_model.py — 正式 CN 的结构模型与来源校验（P0a）

设计边界
========
P0a 只建立“结构模型 -> 是否允许进入正式 CN”的最小契约：

- 诊断拟合继续使用现有 motif/L1 路径，不受影响；
- 正式 CN 必须显式携带 ``StructureModel``；
- P0a 不计算 N₀、不做固定 N 一致性、不输出任何正式 N 值；
- 只有 ``provenance_status == "verified_source_cif"`` 且 CIF 哈希匹配的
  模型才可能进入后续阶段（P0b/P2 再补齐完整门禁）。

这里刻意保持保守：没有真实源 CIF、没有 sha256、没有明确吸收位，一律
返回拒绝码，而不是猜测。
"""
from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass, field
from pathlib import Path

from xafs_core import XAFSError


class CNError(XAFSError):
    """正式 CN 链的可预期错误；``code`` 用于机器可读的拒绝原因。"""

    def __init__(self, code: str, message: str):
        self.code = code
        super().__init__("CN", f"{code}: {message}")


@dataclass
class ShellN0:
    """一个 CN shell 的结构参考配位数。

    ``n0_members`` 保存参与求和的子壳层，便于 P2 复核 Rayleigh 分组。
    """

    group_id: str
    shell_label: str
    scatterer: str
    r_cif: float
    n0: float
    n0_members: list[dict] = field(default_factory=list)
    sub_rayleigh_grouped: bool = False
    merge_reason: str | None = None


@dataclass
class StructureEdit:
    """对源 CIF 的显式编辑记录；P0a 只定义结构，不做编辑。"""

    op: str
    site: str
    old_value: object
    new_value: object
    reason: str
    approved_by: str
    timestamp: str


@dataclass
class StructureModel:
    """正式 CN 使用的结构模型。

    ``source_cif`` 是真正用于计算 N₀ 的 CIF 路径；motif 自带的
    ``cif_path`` 只有升级为 ``verified_source_cif`` 后才能放这里。
    """

    model_id: str
    absorber: str
    absorber_site: str | None = None
    model_level: str = "prototype"
    source_type: str = "none"
    source_ref: str | None = None
    source_cif: str | None = None
    source_cif_sha256: str | None = None
    cell: dict | None = None
    spacegroup: str | None = None
    asym_sites: list | None = None
    cif_path: str | None = None
    occupancy_policy: str = "rejected"
    edits: list[StructureEdit] = field(default_factory=list)
    assumptions: dict = field(default_factory=dict)
    n0_shells: list[ShellN0] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    provenance_status: str = "unverified"

    def to_json(self) -> dict:
        """转成可写入 JSON 的纯 dict；不包含 Path/ndarray 等对象。"""
        return asdict(self)


def sha256_file(path: str | Path) -> str:
    """流式计算文件 SHA256，避免把大 CIF 一次性读入内存。"""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def model_from_motif(motif, *, source_type: str | None = None,
                     source_ref: str | None = None,
                     source_cif_sha256: str | None = None,
                     absorber_site: str | None = None,
                     occupancy_policy: str | None = None,
                     model_level: str | None = None) -> StructureModel:
    """从 motif 构造结构模型；P0a 只做保守分类，不自动升级来源。

    关键约束：
    - ``source_type`` 必须由调用方显式声明为 ``user_cif`` / ``database_cif``；
      仅凭 ``motif.source.ref="COD-..."`` 不能算真实源 CIF；
    - 现场生成的 CIF 默认 ``source_ref_only`` 或 ``unverified``；
    - 没有 sha256 的模型不会通过 ``validate_formal``。
    """
    src = dict(getattr(motif, "source", None) or {})
    src_type = (source_type or src.get("type") or "none")
    ref = (source_ref if source_ref is not None
           else src.get("ref") or src.get("id"))
    cif_path = getattr(motif, "cif_path", None)
    sha = source_cif_sha256
    if sha is None and src_type in ("user_cif", "database_cif") and cif_path:
        path = Path(cif_path)
        if path.exists():
            sha = sha256_file(path)

    if model_level is None:
        if (src_type in ("user_cif", "database_cif") and sha
                and getattr(motif, "kind", None) == "crystal"):
            model_level = "crystal_cif"
        elif getattr(motif, "kind", None) == "cluster":
            model_level = "cluster"
        else:
            model_level = "prototype"

    if occupancy_policy is None:
        occs_lt_one = any(
            float(site.get("occ", 1.0)) < 0.999
            for site in (getattr(motif, "asym_sites", None) or []))
        occupancy_policy = ("explicit" if (cif_path or occs_lt_one)
                            else "implicit_1.0")

    if src_type in ("user_cif", "database_cif") and sha and cif_path:
        provenance = "verified_source_cif"
    elif ref:
        provenance = "source_ref_only"
    elif getattr(motif, "kind", None) == "cluster":
        provenance = "constructed"
    else:
        provenance = "unverified"

    shells: list[ShellN0] = []
    for i, shell in enumerate(getattr(motif, "shell_expect", None) or []):
        element = str(shell.get("el", "?"))
        r_cif = float(shell.get("r", 0.0))
        shells.append(ShellN0(
            group_id=f"g{i}",
            shell_label=f"{element}@{r_cif:.3f}",
            scatterer=element,
            r_cif=r_cif,
            n0=float(shell.get("n", 0.0)),
            n0_members=[dict(shell)],
        ))

    warnings: list[str] = []
    if provenance != "verified_source_cif":
        warnings.append(
            "P0a 保守分类：该 motif 未提供可校验的源 CIF + sha256，"
            "不得进入正式 CN 表")
    if cif_path and provenance != "verified_source_cif":
        warnings.append(
            "motif.cif_path 可能是现场生成文件；只有显式声明为 "
            "user_cif/database_cif 并提供 sha256 才能升级")

    return StructureModel(
        model_id=str(getattr(motif, "id", "unknown")),
        absorber=str(getattr(motif, "absorber", "")),
        absorber_site=absorber_site,
        model_level=model_level,
        source_type=src_type,
        source_ref=str(ref) if ref is not None else None,
        source_cif=str(cif_path) if cif_path else None,
        source_cif_sha256=sha,
        cell=dict(getattr(motif, "cell", None) or {}) or None,
        spacegroup=getattr(motif, "spacegroup", None),
        asym_sites=(list(getattr(motif, "asym_sites", None) or [])
                    or None),
        cif_path=str(cif_path) if cif_path else None,
        occupancy_policy=occupancy_policy,
        edits=[],
        assumptions={"p0a": "provenance scaffold only; no N0 gate yet"},
        n0_shells=shells,
        warnings=warnings,
        provenance_status=provenance,
    )


def validate_formal(model: StructureModel | None) -> dict:
    """P0a 的正式来源校验；返回 ``{ok, code, reason}``，不抛异常。

    P0a 只校验“来源是否可能合格”，不校验 CIF 门禁、占据度、键长、
    一致性拟合；这些在 P0b/P2 实现。
    """
    if model is None:
        return {"ok": False, "code": "no_cif",
                "reason": "未提供 StructureModel（无合格结构来源）"}
    if model.provenance_status != "verified_source_cif":
        return {
            "ok": False,
            "code": "unverified_source",
            "reason": (f"provenance_status={model.provenance_status!r}；"
                       "只有 verified_source_cif 可进入正式 CN"),
        }
    if model.model_level not in ("crystal_cif", "crystal_edited"):
        return {
            "ok": False,
            "code": "unverified_source",
            "reason": (f"model_level={model.model_level!r}；"
                       "P0a 只允许 crystal_cif/crystal_edited"),
        }
    if not model.source_cif:
        return {"ok": False, "code": "no_cif",
                "reason": "模型没有 source_cif 路径"}
    path = Path(model.source_cif)
    if not path.exists():
        return {"ok": False, "code": "source_cif_missing",
                "reason": f"源 CIF 不存在: {path}"}
    if not model.source_cif_sha256:
        return {"ok": False, "code": "source_hash_missing",
                "reason": "缺少 source_cif_sha256，无法复核来源"}
    try:
        actual = sha256_file(path)
    except OSError as exc:
        return {"ok": False, "code": "source_cif_unreadable",
                "reason": f"源 CIF 读取失败: {exc}"}
    if actual.lower() != str(model.source_cif_sha256).lower():
        return {
            "ok": False,
            "code": "source_hash_mismatch",
            "reason": ("源 CIF 当前哈希与记录不一致："
                       f"actual={actual} recorded={model.source_cif_sha256}"),
        }
    if not model.absorber_site:
        return {"ok": False, "code": "absorber_site_missing",
                "reason": "缺少显式 absorber_site；多晶位/正式链不得自动选位"}
    return {"ok": True, "code": None, "reason": None}
