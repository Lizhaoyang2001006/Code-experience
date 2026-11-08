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

from dataclasses import asdict, dataclass, field
from pathlib import Path

from motif_provenance import (
    CRYSTAL_SOURCE_TYPES,
    classify_motif,
    is_generated_cif,
    sha256_file,
)
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


def model_from_motif(motif, *, source_type: str | None = None,
                     source_ref: str | None = None,
                     source_cif_sha256: str | None = None,
                     absorber_site: str | None = None,
                     occupancy_policy: str | None = None,
                     model_level: str | None = None) -> StructureModel:
    """从 motif 构造结构模型；分类规则见 ``motif_provenance``。

    P0b 起，只有 ``Motif.cif_path`` 指向**外部真实 CIF**（不在
    ``_structure_lib/<元素>/_cif/`` 下）时，才有可能升级为
    ``verified_source_cif``。现场生成的 CIF 永远不进入 ``source_cif``。
    """
    src = dict(getattr(motif, "source", None) or {})
    fields = classify_motif(motif)
    real_cif = fields["_source_cif"]

    src_type = (source_type or src.get("type") or "none")
    ref = (source_ref if source_ref is not None
           else fields["source_ref"])
    site = (absorber_site if absorber_site is not None
            else fields["absorber_site"])
    sha = (source_cif_sha256 if source_cif_sha256 is not None
           else fields["source_cif_sha256"])
    status = fields["provenance_status"]
    level = model_level if model_level is not None else fields["model_level"]

    # 调用方显式补充吸收位/哈希时，允许把真实 CIF 升级；但绝不从
    # 现场生成 CIF 升级。
    if real_cif is not None and sha and site:
        status = "verified_source_cif"
        if level not in ("crystal_cif", "crystal_edited"):
            level = "crystal_edited" if fields.get("edits") else "crystal_cif"
    elif status == "verified_source_cif":
        status = "unverified"
        level = "prototype"

    if status != "verified_source_cif":
        sha = None
    if src_type not in CRYSTAL_SOURCE_TYPES and status == "verified_source_cif":
        src_type = "user_cif"

    if occupancy_policy is None:
        occupancy_policy = fields["occupancy_source"]

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
    if status != "verified_source_cif":
        warnings.append(
            f"P0b 来源分层：provenance_status={status!r}，"
            "不得进入正式 CN 表")
    if fields["_has_generated_cif"]:
        warnings.append(
            "motif 的 _cif 文件由 structure_lib 现场生成，"
            "不能作为 source_cif")

    return StructureModel(
        model_id=str(getattr(motif, "id", "unknown")),
        absorber=str(getattr(motif, "absorber", "")),
        absorber_site=site,
        model_level=level,
        source_type=src_type,
        source_ref=str(ref) if ref is not None else None,
        source_cif=real_cif if status == "verified_source_cif" else None,
        source_cif_sha256=sha,
        cell=dict(getattr(motif, "cell", None) or {}) or None,
        spacegroup=getattr(motif, "spacegroup", None),
        asym_sites=(list(getattr(motif, "asym_sites", None) or [])
                    or None),
        cif_path=(str(getattr(motif, "cif_path", None))
                  if getattr(motif, "cif_path", None) else None),
        occupancy_policy=occupancy_policy,
        edits=list(fields.get("edits") or []),
        assumptions={"p0b": "provenance classified; "
                             "formal CN requires verified source CIF"},
        n0_shells=shells,
        warnings=warnings,
        provenance_status=status,
    )


def validate_formal(model: StructureModel | None) -> dict:
    """P0b 的正式来源校验；返回 ``{ok, code, reason}``，不抛异常。

    P0b 只校验“来源是否可能合格”，不校验 CIF 门禁、占据度、键长、
    一致性拟合；这些在 P2 实现。
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
    if is_generated_cif(path):
        return {
            "ok": False,
            "code": "generated_cif_not_formal",
            "reason": ("source_cif 位于 _structure_lib/*/_cif/ 或带生成标记；"
                       "现场生成 CIF 不得作为正式来源"),
        }
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
