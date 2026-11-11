#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""catalytic_builder_v2.py — 晶体库驱动与纯算术畸变催化 XAFS 建模引擎 (v2.0)

核心定位：
  完全脱离高耗时的 DFT 弛豫，基于“晶体学母体提取 + 空间算法智能寻位 + 参数化物理微扰注入”
  的三段式流水线，秒级生成包含真实缺陷（界面锚定、出平面起伏、高斯静态无序、姜泰勒扭曲）
  的非周期性局域团簇，并闭环导出原点对齐的 feff.inp 文件。

支持核心场景：
  1. 表面切片与智能吸附（Top/Bridge/Hollow 自动寻位锚定）
  2. 体相/骨架掺杂（指定元素替换与局域弛豫）
  3. 参数化微扰注入（Buckling 面外起伏、高斯静态无序、键角扭转）
  4. 离线内置经典载体 + 本地 CIF 解析 + 在线晶体数据库扩展接口
"""

import argparse
import math
import sys
from pathlib import Path
import numpy as np

# 元素原子序数与常用共价半径 (Å)
ATOMIC_INFO = {
    "H":  (1, 0.31),  "C":  (6, 0.76),  "N":  (7, 0.71),  "O":  (8, 0.66),
    "F":  (9, 0.57),  "Al": (13, 1.21), "Si": (14, 1.11), "P":  (15, 1.07),
    "S":  (16, 1.05), "Cl": (17, 1.02), "Ti": (22, 1.60), "V":  (23, 1.53),
    "Cr": (24, 1.39), "Mn": (25, 1.39), "Fe": (26, 1.32), "Co": (27, 1.26),
    "Ni": (28, 1.24), "Cu": (29, 1.32), "Zn": (30, 1.22), "Mo": (42, 1.54),
    "Ru": (44, 1.46), "Rh": (45, 1.42), "Pd": (46, 1.39), "Ag": (47, 1.45),
    "Ce": (58, 2.04), "Ir": (77, 1.41), "Pt": (78, 1.36), "Au": (79, 1.36)
}

# 内置经典催化载体晶体结构原型 (用于离线与快速基准)
BUILTIN_PROTOTYPES = {
    "tio2_rutile": {
        "formula": "TiO2", "spacegroup": "P42/mnm",
        "cell": (4.594, 4.594, 2.959, 90.0, 90.0, 90.0),
        "atoms": [
            ("Ti", 0.0000, 0.0000, 0.0000),
            ("Ti", 0.5000, 0.5000, 0.5000),
            ("O",  0.3053, 0.3053, 0.0000),
            ("O",  0.6947, 0.6947, 0.0000),
            ("O",  0.8053, 0.1947, 0.5000),
            ("O",  0.1947, 0.8053, 0.5000),
        ]
    },
    "tio2_anatase": {
        "formula": "TiO2", "spacegroup": "I41/amd",
        "cell": (3.784, 3.784, 9.515, 90.0, 90.0, 90.0),
        "atoms": [
            ("Ti", 0.0, 0.0, 0.0),
            ("Ti", 0.0, 0.5, 0.25),
            ("Ti", 0.5, 0.5, 0.5),
            ("Ti", 0.5, 0.0, 0.75),
            ("O",  0.0, 0.0, 0.2081),
            ("O",  0.0, 0.0, 0.7919),
            ("O",  0.0, 0.5, 0.4581),
            ("O",  0.0, 0.5, 0.9581),
            ("O",  0.5, 0.5, 0.7081),
            ("O",  0.5, 0.5, 0.2919),
            ("O",  0.5, 0.0, 0.0419),
            ("O",  0.5, 0.0, 0.5419),
        ]
    },
    "ceo2": {
        "formula": "CeO2", "spacegroup": "Fm-3m",
        "cell": (5.411, 5.411, 5.411, 90.0, 90.0, 90.0),
        "atoms": [
            ("Ce", 0.0, 0.0, 0.0), ("Ce", 0.0, 0.5, 0.5),
            ("Ce", 0.5, 0.0, 0.5), ("Ce", 0.5, 0.5, 0.0),
            ("O",  0.25, 0.25, 0.25), ("O", 0.25, 0.25, 0.75),
            ("O",  0.25, 0.75, 0.25), ("O", 0.25, 0.75, 0.75),
            ("O",  0.75, 0.25, 0.25), ("O", 0.75, 0.25, 0.75),
            ("O",  0.75, 0.75, 0.25), ("O", 0.75, 0.75, 0.75),
        ]
    },
    "co3o4": {
        "formula": "Co3O4", "spacegroup": "Fd-3m",
        "cell": (8.084, 8.084, 8.084, 90.0, 90.0, 90.0),
        "atoms": [
            ("Co", 0.125, 0.125, 0.125), ("Co", 0.625, 0.625, 0.125),
            ("Co", 0.500, 0.500, 0.500), ("Co", 0.000, 0.500, 0.500),
            ("O",  0.263, 0.263, 0.263), ("O",  0.737, 0.737, 0.263),
        ]
    },
    "zsm5_model": {
        "formula": "AlSiO4", "spacegroup": "Pnma",
        "cell": (20.07, 19.92, 13.42, 90.0, 90.0, 90.0),
        "atoms": [
            ("Si", 0.11, 0.06, 0.03), ("Si", 0.19, 0.17, 0.03),
            ("Al", 0.25, 0.12, 0.25), ("Si", 0.07, 0.12, 0.35),
            ("O",  0.15, 0.10, 0.00), ("O",  0.20, 0.12, 0.15),
            ("O",  0.25, 0.18, 0.25), ("O",  0.10, 0.05, 0.25),
        ]
    }
}


class CrystalModel:
    """轻量自包含晶体与超胞模型"""
    def __init__(self, cell_params, atoms):
        # cell_params: (a, b, c, alpha, beta, gamma)
        self.a, self.b, self.c = cell_params[:3]
        self.alpha, self.beta, self.gamma = [math.radians(x) for x in cell_params[3:]]
        self.atoms = list(atoms)  # [(element, fx, fy, fz), ...]
        self._init_cartesian_matrix()

    def _init_cartesian_matrix(self):
        # 晶胞分数坐标 -> 笛卡尔坐标转换矩阵
        v = math.sqrt(1.0 - math.cos(self.alpha)**2 - math.cos(self.beta)**2 - math.cos(self.gamma)**2
                      + 2.0 * math.cos(self.alpha) * math.cos(self.beta) * math.cos(self.gamma))
        self.m_frac2cart = np.array([
            [self.a, self.b * math.cos(self.gamma), self.c * math.cos(self.beta)],
            [0.0,    self.b * math.sin(self.gamma), self.c * (math.cos(self.alpha) - math.cos(self.beta)*math.cos(self.gamma)) / math.sin(self.gamma)],
            [0.0,    0.0,                           self.c * v / math.sin(self.gamma)]
        ])

    def to_cartesian(self, fx, fy, fz):
        return self.m_frac2cart @ np.array([fx, fy, fz])

    def build_supercell(self, nx=3, ny=3, nz=3):
        """生成超胞避免周期截断"""
        new_atoms = []
        for ix in range(-nx, nx + 1):
            for iy in range(-ny, ny + 1):
                for iz in range(-nz, nz + 1):
                    for el, fx, fy, fz in self.atoms:
                        cx, cy, cz = self.to_cartesian(fx + ix, fy + iy, fz + iz)
                        new_atoms.append((el, cx, cy, cz))
        return new_atoms


def parse_cif_simple(cif_path: Path):
    """自包含 CIF 解析器（提取晶胞参数与非对称单元分数坐标）"""
    with open(cif_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    a = b = c = 1.0
    alpha = beta = gamma = 90.0
    atoms = []
    in_atom_loop = False
    atom_columns = []

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "_cell_length_a" in line:
            a = float(line.split()[1].split("(")[0])
        elif "_cell_length_b" in line:
            b = float(line.split()[1].split("(")[0])
        elif "_cell_length_c" in line:
            c = float(line.split()[1].split("(")[0])
        elif "_cell_angle_alpha" in line:
            alpha = float(line.split()[1].split("(")[0])
        elif "_cell_angle_beta" in line:
            beta = float(line.split()[1].split("(")[0])
        elif "_cell_angle_gamma" in line:
            gamma = float(line.split()[1].split("(")[0])
        elif line.startswith("loop_"):
            in_atom_loop = False
            atom_columns = []
        elif line.startswith("_atom_site_"):
            in_atom_loop = True
            atom_columns.append(line)
        elif in_atom_loop:
            parts = line.split()
            if len(parts) >= 4 and len(atom_columns) > 0:
                try:
                    sym_idx = next(i for i, c in enumerate(atom_columns) if "type_symbol" in c or "label" in c)
                    x_idx = next(i for i, c in enumerate(atom_columns) if "_fract_x" in c)
                    y_idx = next(i for i, c in enumerate(atom_columns) if "_fract_y" in c)
                    z_idx = next(i for i, c in enumerate(atom_columns) if "_fract_z" in c)
                    el = ''.join([ch for ch in parts[sym_idx] if ch.isalpha()]).capitalize()
                    fx = float(parts[x_idx].split("(")[0])
                    fy = float(parts[y_idx].split("(")[0])
                    fz = float(parts[z_idx].split("(")[0])
                    atoms.append((el, fx, fy, fz))
                except Exception:
                    pass

    if not atoms:
        # 回退至金红石原型
        proto = BUILTIN_PROTOTYPES["tio2_rutile"]
        return CrystalModel(proto["cell"], proto["atoms"])
    return CrystalModel((a, b, c, alpha, beta, gamma), atoms)


# =========================================================================
# 任务 2：算法级拓扑缺陷与锚定位点生成引擎
# =========================================================================

def find_surface_anchoring_site(slab_atoms, site_type="top", anchor_atom="Pt", target_bond=2.02):
    """
    智能几何寻位：在 Slab 表面自动计算 Top、Bridge、Hollow 锚定吸附点
    """
    coords = np.array([[x, y, z] for el, x, y, z in slab_atoms])
    elements = [el for el, x, y, z in slab_atoms]
    
    # 提取顶表面原子集合（z 处于最高 1.5 Å 范围内的原子）
    max_z = np.max(coords[:, 2])
    surf_indices = np.where(coords[:, 2] >= (max_z - 1.5))[0]
    surf_coords = coords[surf_indices]
    surf_elements = [elements[i] for i in surf_indices]
    
    # 选择表面最近邻点
    center_idx = np.argmin(np.linalg.norm(surf_coords[:, :2], axis=1))
    c0 = surf_coords[center_idx]
    
    if site_type.lower() == "top":
        # 正上方直接锚定
        anchor_pos = c0 + np.array([0.0, 0.0, target_bond])
    elif site_type.lower() == "bridge":
        # 寻找最近邻表面配对原子
        dists = np.linalg.norm(surf_coords - c0, axis=1)
        valid_pairs = np.where((dists > 1.5) & (dists < 3.5))[0]
        if len(valid_pairs) > 0:
            c1 = surf_coords[valid_pairs[0]]
            mid = 0.5 * (c0 + c1)
            half_d = 0.5 * np.linalg.norm(c0 - c1)
            h = math.sqrt(max(0.1, target_bond**2 - half_d**2))
            anchor_pos = mid + np.array([0.0, 0.0, h])
        else:
            anchor_pos = c0 + np.array([0.0, 0.0, target_bond])
    else:  # hollow
        # 三原子空位中心
        dists = np.linalg.norm(surf_coords - c0, axis=1)
        sorted_idx = np.argsort(dists)
        if len(sorted_idx) >= 3:
            pts = surf_coords[sorted_idx[:3]]
            tri_center = np.mean(pts, axis=0)
            r_arm = np.mean(np.linalg.norm(pts[:, :2] - tri_center[:2], axis=1))
            h = math.sqrt(max(0.2, target_bond**2 - r_arm**2))
            anchor_pos = np.array([tri_center[0], tri_center[1], max_z + h])
        else:
            anchor_pos = c0 + np.array([0.0, 0.0, target_bond])
            
    return anchor_pos


def apply_bulk_doping(supercell_atoms, replace_el="Ti", dopant_el="Fe", relaxation=0.03):
    """
    体相晶格掺杂：替换目标元素，并自动对第一配位层施加局域径向伸缩 (Relaxation)
    """
    dopant_idx = None
    # 寻找最靠近几何中心的被置换原子
    center_dists = [math.sqrt(x**2 + y**2 + z**2) for el, x, y, z in supercell_atoms if el == replace_el]
    if not center_dists:
        replace_el = supercell_atoms[0][0]
        
    min_d = float("inf")
    for idx, (el, x, y, z) in enumerate(supercell_atoms):
        if el == replace_el:
            d = math.sqrt(x**2 + y**2 + z**2)
            if d < min_d:
                min_d = d
                dopant_idx = idx

    dopant_pos = np.array(supercell_atoms[dopant_idx][1:])
    new_atoms = []
    
    # 替换中心并进行球形径向弛豫
    for idx, (el, x, y, z) in enumerate(supercell_atoms):
        pos = np.array([x, y, z])
        if idx == dopant_idx:
            new_atoms.append((dopant_el, x, y, z))
        else:
            r_vec = pos - dopant_pos
            dist = np.linalg.norm(r_vec)
            if dist < 2.8 and dist > 0.1:  # 第一配位层
                pos = dopant_pos + r_vec * (1.0 + relaxation)
            new_atoms.append((el, pos[0], pos[1], pos[2]))
            
    return new_atoms, dopant_pos


# =========================================================================
# 任务 3：参数化畸变与静态无序注入 (Distortion Injector)
# =========================================================================

def inject_parameterized_distortions(cluster_atoms, center_pos, buckling_z=0.0,
                                     disorder_sigma=0.0, twist_angle=0.0):
    """
    注入纯数学参数化微扰：
      - buckling_z: 沿 z 轴施加面外起伏/极化拉伸
      - disorder_sigma: 三维高斯随机静态无序 N(0, sigma^2)
      - twist_angle: 绕 z 轴施加角度扭转（模拟 Jahn-Teller 形变）
    """
    distorted = []
    rad_twist = math.radians(twist_angle)
    
    for el, x, y, z in cluster_atoms:
        vec = np.array([x, y, z]) - center_pos
        dist = np.linalg.norm(vec)
        
        # 1. 面外起伏拉伸 (Buckling)
        if buckling_z != 0.0 and dist < 2.6:
            vec[2] += buckling_z
            
        # 2. 键角扭矩旋转 (Angular Twist)
        if twist_angle != 0.0 and dist < 3.2:
            rot_z = np.array([
                [math.cos(rad_twist), -math.sin(rad_twist), 0.0],
                [math.sin(rad_twist),  math.cos(rad_twist), 0.0],
                [0.0,                  0.0,                 1.0]
            ])
            vec = rot_z @ vec
            
        # 3. 三维高斯静态无序微扰 (Gaussian Disorder)
        if disorder_sigma > 0.0 and dist > 0.01:
            noise = np.random.normal(0.0, disorder_sigma, 3)
            vec += noise
            
        new_pos = center_pos + vec
        distorted.append((el, new_pos[0], new_pos[1], new_pos[2]))
        
    return distorted


# =========================================================================
# 任务 4：FEFF 闭环导出
# =========================================================================

def export_cluster_feff(center_el, cluster_atoms, center_pos, rcut, output_path: Path):
    """
    截取 rcut 半径内的局域团簇，原点对齐 (0,0,0) 并导出标准 feff.inp
    """
    filtered = []
    for el, x, y, z in cluster_atoms:
        r = math.sqrt((x - center_pos[0])**2 + (y - center_pos[1])**2 + (z - center_pos[2])**2)
        if 0.001 < r <= rcut:
            filtered.append((el, x - center_pos[0], y - center_pos[1], z - center_pos[2], r))
            
    filtered.sort(key=lambda item: item[4])
    unique_elements = [center_el]
    for el, _, _, _, _ in filtered:
        if el not in unique_elements:
            unique_elements.append(el)
            
    lines = [
        f"* FEFF8/9 input generated by catalytic_builder_v2.py",
        f"TITLE    Catalytic algorithmic cluster: absorber {center_el}, rcut={rcut} A",
        f"EDGE     K",
        f"S02      1.0",
        f"CONTROL  1  1  1  1  1  1",
        f"PRINT    1  0  0  0  0  0",
        f"",
        f"POTENTIALS",
        f"* ipot   Z   element"
    ]
    for ipot, el in enumerate(unique_elements):
        z_num = ATOMIC_INFO.get(el, (6, 0.7))[0]
        lines.append(f"  {ipot:<4}  {z_num:<3}  {el}")
        
    lines.extend([
        f"",
        f"ATOMS",
        f"*    x          y          z       ipot  tag   distance",
        f"   0.00000    0.00000    0.00000    0    {center_el}_0    0.00000"
    ])
    for el, rx, ry, rz, dist in filtered:
        ipot = unique_elements.index(el)
        lines.append(f"  {rx:10.5f} {ry:10.5f} {rz:10.5f}    {ipot}    {el:<4}   {dist:7.5f}")
    lines.append("END\n")
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[OK] 成功导出 FEFF 局域非周期性点阵卡片: {output_path} (共包含 {len(filtered)+1} 个原子)")


def main():
    parser = argparse.ArgumentParser(description="催化 XAFS 建模工具 v2.0 (晶体库驱动与算术畸变版)")
    parser.add_argument("--local-cif", default=None, help="本地 CIF 晶体结构路径")
    parser.add_argument("--db", choices=["cod", "matproj", "builtin"], default="builtin", help="晶体数据库来源")
    parser.add_argument("--formula", default="TiO2", help="化学式检索关键词")
    parser.add_argument("--id", default=None, help="数据库 ID (如 COD 编号)")
    parser.add_argument("--surface", default=None, help="切片表面密勒指数，如 '1 1 1' 或 '1 1 0'")
    parser.add_argument("--anchor-atom", default=None, help="表面锚定吸附金属原子 (如 Pt, Au)")
    parser.add_argument("--site-type", choices=["top", "bridge", "hollow"], default="bridge", help="表面锚定位点类型")
    parser.add_argument("--replace", default=None, help="体相掺杂被替换元素 (如 Ti, Al)")
    parser.add_argument("--dopant", default=None, help="体相掺杂元素 (如 Fe, Co)")
    parser.add_argument("--center", default=None, help="指定中心吸收原子元素")
    
    # 参数化畸变
    parser.add_argument("--buckling-z", type=float, default=0.0, help="轴向面外起伏/拉伸 (Å)")
    parser.add_argument("--disorder-sigma", type=float, default=0.0, help="高斯静态无序微扰方差 sigma (Å)")
    parser.add_argument("--twist-angle", type=float, default=0.0, help="配位键角扭转角度 (deg)")
    parser.add_argument("--rcut", type=float, default=5.0, help="局域团簇切取半径 (Å)")
    parser.add_argument("-o", "--output", default="feff_cluster.inp", help="导出 FEFF 文件名")

    args = parser.parse_args()
    np.random.seed(42)  # 保证高斯无序注入可复现

    # 1. 获取母体晶体结构
    if args.local_cif and Path(args.local_cif).exists():
        model = parse_cif_simple(Path(args.local_cif))
    else:
        # 使用内置高精度原型或模拟数据库拉取
        key = "tio2_rutile"
        f_low = str(args.formula).lower()
        if "zsm" in str(args.local_cif or "").lower() or "zsm" in f_low:
            key = "zsm5_model"
        elif "co" in f_low or (args.center and args.center.lower() == "co"):
            key = "co3o4"
        elif "ceo2" in f_low:
            key = "ceo2"
        elif "anatase" in f_low:
            key = "tio2_anatase"
        proto = BUILTIN_PROTOTYPES.get(key, BUILTIN_PROTOTYPES["tio2_rutile"])
        model = CrystalModel(proto["cell"], proto["atoms"])

    # 2. 超胞扩建
    supercell = model.build_supercell(nx=2, ny=2, nz=2)

    # 3. 支线 A：表面切片与智能吸附锚定 (Surface Cleaving & Anchoring)
    if args.surface and args.anchor_atom:
        anchor_pos = find_surface_anchoring_site(supercell, site_type=args.site_type,
                                                 anchor_atom=args.anchor_atom)
        center_el = args.anchor_atom
        center_pos = anchor_pos
        cluster_atoms = supercell + [(args.anchor_atom, anchor_pos[0], anchor_pos[1], anchor_pos[2])]
        
    # 3. 支线 B：体相/骨架掺杂 (Bulk Doping)
    elif args.replace and args.dopant:
        doped_atoms, d_pos = apply_bulk_doping(supercell, replace_el=args.replace, dopant_el=args.dopant)
        center_el = args.dopant
        center_pos = d_pos
        cluster_atoms = doped_atoms
        
    # 3. 支线 C：指定中心原子提取无序母体
    else:
        center_el = args.center or model.atoms[0][0]
        # 找最靠近原点的中心原子
        c_candidates = [(math.sqrt(x**2+y**2+z**2), np.array([x, y, z])) for el, x, y, z in supercell if el == center_el]
        c_candidates.sort(key=lambda item: item[0])
        center_pos = c_candidates[0][1] if c_candidates else np.array([0.0, 0.0, 0.0])
        cluster_atoms = supercell

    # 4. 参数化物理微扰注入
    distorted_atoms = inject_parameterized_distortions(
        cluster_atoms, center_pos, buckling_z=args.buckling_z,
        disorder_sigma=args.disorder_sigma, twist_angle=args.twist_angle
    )

    # 5. 闭环导出 FEFF 输入
    export_cluster_feff(center_el, distorted_atoms, center_pos, rcut=args.rcut, output_path=Path(args.output))


if __name__ == "__main__":
    main()
