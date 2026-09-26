#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
实时监看 JOREK HDF5 文件中选定物理场在 (R,Z) 截面的剖面形状。

特性：
- 周期性扫描目录中匹配 pattern 的 jorek*.h5 文件
- 只用 numpy + h5py + matplotlib，不依赖 torch
- 使用 JorekDataParser 直接读取非结构化网格与物理量
- 使用 matplotlib.tri 在原始非结构化网格上绘制等值面（无光栅化）
- 命令行参数选择目录、文件模式 (pattern)、物理场、刷新时间等

依赖：
- numpy
- h5py
- matplotlib

用法示例：
    cd /home/ac/work_dir/EHL-2/wallgrid_const
    python util/monitor_jorek_field.py \
        --dir . \
        --field psi \
        --interval 60.0 \
        --pattern "jorek*.h5"

默认只监看最新的 jorek*.h5，每 60 秒刷新一次；
也可以通过 --pattern 选择其他文件模式（例如 "jorek*.h5" 会始终选最新的文件）。
支持的场名：psi, phi, j, omega, rho, T, v_par, T_e, rho_n
plane 使用 1-based 索引，对应第 plane 个物理环向平面，环向角为
phi = 2*pi*(plane-1)/(n_period*n_plane)；该平面上的场由 HDF5 中
n_tor 个环向谐波按 Moivre 基实时重建（与 mode_moivre / jorek2vtk_3d 一致）。
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import h5py
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np


@dataclass
class GraphData:
    """简化版数据结构，仅包含可视化所需信息（全部为 numpy 数组）。"""

    x: np.ndarray          # [n_nodes, n_fields] 节点物理量矩阵
    field_names: List[str]  # x 的列名（规范化后的场名）
    pos: np.ndarray        # [n_nodes, 2] 节点坐标 (R, Z)
    triangles: Optional[np.ndarray]  # [n_tri, 3] 三角形拓扑 (0-based)，或 None
    time: float            # 当前物理时间 t_now
    plane: Optional[int] = None   # 当前环向平面 (1-based)
    phi: Optional[float] = None   # 当前环向角 (rad)


def normalize_field_name(field_name: str) -> str:
    return field_name.strip().lower().replace("-", "_").replace(" ", "_")


def decode_h5_strings(raw_names) -> List[str]:
    names: List[str] = []
    for value in np.asarray(raw_names).reshape(-1):
        if isinstance(value, bytes):
            names.append(value.decode("utf-8", errors="ignore"))
        else:
            names.append(str(value))
    return names


def find_string_list_in_h5(handle: h5py.File, candidates: List[str]) -> Optional[List[str]]:
    for key in candidates:
        if key in handle:
            return decode_h5_strings(handle[key][()])
        if key in handle.attrs:
            return decode_h5_strings(handle.attrs[key])
    return None


CANONICAL_FIELD_ORDER = ["psi", "phi", "j", "omega", "rho", "t", "v_par", "t_e", "rho_n"]

FIELD_ALIASES = {
    "psi": ["psi"],
    "phi": ["phi", "u", "potential"],
    "j": ["j", "zj", "current", "current_density"],
    "omega": ["omega", "w"],
    "rho": ["rho"],
    "t": ["t", "ti", "temperature"],
    "v_par": ["v_par", "vpar"],
    "t_e": ["t_e", "te"],
    "rho_n": ["rho_n", "rhon"],
}

LEGACY_FIELD_ORDER = ["psi", "phi", "j", "omega", "rho", "t", "v_par", "rho_n", "t_e"]


def moivre_weights(n_tor: int, n_period: int, phi: float) -> np.ndarray:
    """JOREK 环向 Moivre 基在角度 phi 处的取值（与 mode_moivre 一致）。

    HZ(1)    = 1                      (n=0 模)
    HZ(2i)   = cos(i*n_period*phi)    (第 i 个谐波的余弦分量)
    HZ(2i+1) = sin(i*n_period*phi)    (第 i 个谐波的正弦分量)
    """
    weights = np.zeros(max(n_tor, 0), dtype=np.float64)
    if n_tor > 0:
        weights[0] = 1.0
    for i in range(1, (n_tor - 1) // 2 + 1):
        phase = float(i) * float(n_period) * float(phi)
        # Fortran 1-based: HZ(2i)=cos, HZ(2i+1)=sin -> Python 0-based: 2i-1 / 2i
        weights[2 * i - 1] = np.cos(phase)
        weights[2 * i] = np.sin(phase)
    return weights


def infer_total_planes(h5_filepath: str, plane_axis: Optional[int] = None) -> int:
    """返回文件中可监看的环向平面总数。

    JOREK restart 文件里 n_plane 是环向平面的离散数（例如 marker_spi_source 的
    SPI 源使用 32 个环向平面），而 values 数组本身只存 n_tor 个环向谐波；
    任意 n_plane 个物理平面都可以由这些谐波重建，因此直接采用元数据 n_plane。
    """
    with h5py.File(h5_filepath, "r") as f:
        for dataset_name in ("n_plane", "n_tor"):
            if dataset_name in f:
                arr = np.asarray(f[dataset_name][...]).reshape(-1)
                if arr.size > 0:
                    value = int(arr[0])
                    if value > 0:
                        return value
    return 1


def resolve_field_layout(variable_names: Optional[List[str]]) -> Tuple[List[str], Dict[str, int]]:
    if variable_names is None:
        field_names = LEGACY_FIELD_ORDER.copy()
        name_to_index = {name: idx for idx, name in enumerate(field_names)}
        return field_names, name_to_index

    normalized_names = [normalize_field_name(name) for name in variable_names]
    name_to_index_raw = {name: idx for idx, name in enumerate(normalized_names)}

    field_names: List[str] = []
    name_to_index: Dict[str, int] = {}
    for canonical in CANONICAL_FIELD_ORDER:
        for alias in FIELD_ALIASES.get(canonical, [canonical]):
            alias_key = normalize_field_name(alias)
            if alias_key in name_to_index_raw:
                field_names.append(canonical)
                name_to_index[canonical] = name_to_index_raw[alias_key]
                break

    return field_names, name_to_index


def resolve_requested_field(field_name: str, available_field_names: List[str]) -> str:
    canonical = normalize_field_name(field_name)
    if canonical in available_field_names:
        return canonical
    raise ValueError(
        f"Unsupported field '{field_name}'. Available fields: {available_field_names}",
    )


class JorekDataParser:
    """解析 JOREK HDF5 输出文件，提取网格与物理量（纯 numpy 版本）。"""

    def __init__(self, h5_filepath: str) -> None:
        self.filepath = h5_filepath
        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f"Cannot find HDF5 file: {self.filepath}")

    def parse_to_graph(self, plane: int = 1, plane_axis: Optional[int] = None) -> GraphData:
        """解析 JOREK restart HDF5 并返回 GraphData。

        marker_spi_source 运行产生的 jorek*.h5 真实布局（Fortran 列主序写盘后
        在 HDF5 中维度反转，以真实文件为准）：

          x      : (n_dim, n_degrees, n_coord_tor, n_nodes)
          values : (n_var,  n_degrees, n_tor,       n_nodes)
          vertex : (n_vertex_max, n_elements)          ! 1-based 节点索引

        节点坐标取 degree 0、第 1 个坐标谐波（网格在坐标上是轴对称的，
        n_coord_tor=1）；物理量取 degree 0 后，按环向角
        phi = 2*pi*(plane-1)/(n_period*n_plane) 对 n_tor 个环向谐波做
        Moivre 展开求和，得到该物理平面上的场值（与 mode_moivre、
        jorek2vtk_3d 完全一致）。
        """
        with h5py.File(self.filepath, "r") as f:
            def _get_int(name: str, default: int) -> int:
                if name in f:
                    arr = np.asarray(f[name][...]).reshape(-1)
                    if arr.size > 0:
                        return int(arr[0])
                return default

            n_tor = _get_int("n_tor", 1)
            n_period = _get_int("n_period", 1)
            n_plane = _get_int("n_plane", max(n_tor, 1))
            n_var = _get_int("n_var", -1)

            if plane < 1:
                raise ValueError(f"Plane index must be >= 1, got {plane}.")
            phi = 2.0 * np.pi * float(plane - 1) / float(max(n_period, 1) * max(n_plane, 1))
            hz = moivre_weights(n_tor, n_period, phi)

            # 1. 节点坐标: x (n_dim, n_degrees, n_coord_tor, n_nodes) -> (n_nodes, n_dim)
            if "x" not in f:
                raise KeyError("Coordinate array 'x' not found in HDF5.")
            x_raw = np.asarray(f["x"][...])
            if x_raw.ndim == 0:
                raise ValueError(f"Unexpected scalar shape for 'x': {x_raw.shape}")
            if x_raw.ndim == 4:
                pos = x_raw[:, 0, 0, :].T.copy()      # degree 0, 坐标谐波 0
            elif x_raw.ndim == 3:
                pos = x_raw[:, 0, :].T.copy()         # rst v1: (n_dim, n_degrees, n_nodes)
            else:
                raise ValueError(f"Unexpected rank for 'x': {x_raw.shape}")
            if pos.shape[1] < 2:
                raise ValueError(f"Coordinate array 'x' has too few spatial components: {x_raw.shape}")
            pos = np.array(pos[:, :2], dtype=np.float32)
            n_nodes = pos.shape[0]

            # 2. 物理量: values (n_var, n_degrees, n_tor, n_nodes)
            #    取 degree 0，再对环向谐波做 Moivre 求和 -> (n_nodes, n_var)
            if "values" not in f:
                raise KeyError("Physical values array 'values' not found.")
            values_raw = np.asarray(f["values"][...])
            if values_raw.ndim != 4:
                raise ValueError(f"Unexpected rank for 'values': {values_raw.shape}")
            var_tor_node = values_raw[:, 0, :, :]
            if var_tor_node.shape[1] != n_tor:
                raise ValueError(
                    f"'values' toroidal axis size {var_tor_node.shape[1]} "
                    f"!= n_tor metadata {n_tor}.",
                )
            if var_tor_node.shape[2] != n_nodes:
                raise ValueError(
                    f"'values' node axis size {var_tor_node.shape[2]} "
                    f"!= x node count {n_nodes}.",
                )
            field_var_node = np.einsum("vtn,t->vn", var_tor_node, hz)

            variable_names = find_string_list_in_h5(f, ["variables", "variable_names", "var_names"])
            if variable_names is not None:
                field_names, name_to_index = resolve_field_layout(variable_names)
                n_chan = field_var_node.shape[0]
                bad = [
                    name
                    for name in field_names
                    if name_to_index.get(name, -1) < 0 or name_to_index[name] >= n_chan
                ]
                if bad:
                    raise ValueError(
                        "Failed to resolve one or more fields from HDF5 variable names: "
                        f"{bad}. Raw names: {variable_names}",
                    )
            else:
                # 文件不含变量名时使用 model 502 的 9 通道约定顺序
                # psi, u(phi), zj(j), w(omega), rho, Ti, Vpar, rhon, Te
                if n_var > 0 and field_var_node.shape[0] != n_var:
                    raise ValueError(
                        f"'values' variable axis size {field_var_node.shape[0]} "
                        f"!= n_var metadata {n_var}.",
                    )
                field_names = LEGACY_FIELD_ORDER[: field_var_node.shape[0]]
                if not field_names:
                    raise ValueError("No variable channels found in 'values'.")
                name_to_index = {name: idx for idx, name in enumerate(field_names)}

            x = np.stack(
                [
                    np.array(field_var_node[name_to_index[fn], :], dtype=np.float32)
                    for fn in field_names
                ],
                axis=1,
            )

            # 3. 时间戳
            if "t_now" in f:
                t_val = f["t_now"][()]
                if isinstance(t_val, np.ndarray):
                    t_val = float(t_val.reshape(-1)[0])
                else:
                    t_val = float(t_val)
            else:
                t_val = 0.0

            # 4. 三角形拓扑（由 vertex 拆分为两个三角形）
            triangles = None
            if "vertex" in f:
                vertex = np.array(f["vertex"][...])
                if vertex.ndim != 2 or vertex.shape[0] < 3:
                    raise ValueError(f"Unexpected shape for 'vertex': {vertex.shape}")
                valid_mask = (vertex > 0).all(axis=0)
                valid_vertex = vertex[:, valid_mask].astype(np.int64) - 1  # -> 0-based

                if valid_vertex.shape[0] == 3:
                    triangles = valid_vertex.T.copy()
                elif valid_vertex.shape[0] >= 4:
                    tri1 = valid_vertex[[0, 1, 2], :].T
                    tri2 = valid_vertex[[0, 2, 3], :].T
                    triangles = np.vstack([tri1, tri2])
                else:
                    triangles = None

            # 5. 简单检查
            if pos.shape[0] != x.shape[0]:
                raise ValueError(
                    f"pos node count {pos.shape[0]} != values node count {x.shape[0]}",
                )

            return GraphData(
                x=x,
                field_names=field_names,
                pos=pos,
                triangles=triangles,
                time=t_val,
                plane=plane,
                phi=phi,
            )


FIELD_CYCLE = ["psi", "phi", "j", "omega", "rho", "t", "v_par", "t_e", "rho_n"]


def find_latest_file(directory: str, pattern: str) -> Optional[str]:
    """在给定目录中查找匹配 pattern 的最新文件（按修改时间排序）。"""
    paths = glob.glob(os.path.join(directory, pattern))
    if not paths:
        return None
    # 选取修改时间最新的文件
    return max(paths, key=os.path.getmtime)


def list_files_by_sequence(directory: str, pattern: str) -> List[str]:
    """按文件名中的最后一个数字序列排序（升序），用于 playback。"""
    paths = glob.glob(os.path.join(directory, pattern))
    if not paths:
        return []

    def seq_key(path: str):
        basename = os.path.basename(path)
        nums = re.findall(r"\d+", basename)
        if nums:
            return (0, int(nums[-1]), basename)
        return (1, basename)

    return sorted(paths, key=seq_key)


def monitor_field(
    directory: str,
    field_name: str,
    interval: float = 60.0,
    cmap: str = "jet",
    pattern: str = "jorek_restart.h5",
    selected_file: Optional[str] = None,
    plane: int = 1,
    plane_axis: Optional[int] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    levels: int = 50,
    log_scale: bool = False,
    playback: bool = False,
    playback_interval: float = 0.3,
    playback_loop: bool = False,
) -> None:
    """循环监看目录中最新匹配 pattern 的 JOREK HDF5 文件的指定物理场。"""
    current_field_name = normalize_field_name(field_name)
    current_field_index: Optional[int] = None
    current_plane = int(plane)
    total_planes = infer_total_planes(find_latest_file(directory, pattern) or selected_file or "", plane_axis) if (selected_file or find_latest_file(directory, pattern)) else 1

    plt.ion()
    fig, ax = plt.subplots(figsize=(8, 10))
    cbar = None

    forced_file: Optional[str] = None
    if selected_file:
        forced_file = selected_file

    last_file: Optional[str] = None
    last_mtime: Optional[float] = None
    force_refresh = False
    available_field_names = FIELD_CYCLE.copy()

    def print_hotkey_help() -> None:
        mode = "playback" if playback else ("fixed-file" if forced_file else "auto-pattern")
        print("\n=== 操作指南 (Hotkeys) ===")
        print("  h       : 显示本帮助与当前状态")
        print("  ,  /  . : 上一个 / 下一个环向平面 (physical toroidal plane)")
        print("            plane 由 HDF5 中 n_tor 个环向谐波按 Moivre 基实时重建,")
        print("            对应环向角 phi = 2*pi*(plane-1)/(n_period*n_plane),")
        print("            标题栏会显示 plane 序号与 phi 角度")
        print("  ]  /  [ : 下一个 / 上一个物理场 (psi, phi, j, omega, rho, ...)")
        print("  z       : 切换到指定的 HDF5 文件")
        print("  a       : 返回自动 pattern 模式")
        print("  Ctrl+C  : 退出监看 (保留当前图像窗口)")
        print("\n=== 当前状态 (Current Status) ===")
        print(f"  mode             : {mode}")
        print(f"  directory        : {directory}")
        print(f"  pattern          : {pattern}")
        print(f"  fixed file       : {forced_file if forced_file else '(none)'}")
        print(f"  field            : {current_field_name}")
        print(f"  fields available : {', '.join(available_field_names)}")
        print(f"  plane            : {current_plane}/{total_planes}")
        print(f"  colormap         : {cmap}")
        print(f"  levels           : {levels}")
        print(f"  log scale        : {log_scale}")
        print(f"  vmin             : {vmin}")
        print(f"  vmax             : {vmax}")
        print(f"  monitor interval : {interval}s")
        print(f"  playback interval: {playback_interval}s")
        print(f"  playback loop    : {playback_loop}")

    def render_file(latest: str) -> bool:
        nonlocal cbar, available_field_names, current_field_index, current_plane, total_planes
        print(f"\nReading {os.path.basename(latest)}", flush=True)
        try:
            total_planes = infer_total_planes(latest, plane_axis)
            if current_plane < 1 or current_plane > total_planes:
                current_plane = 1
            parser = JorekDataParser(latest)
            graph = parser.parse_to_graph(plane=current_plane, plane_axis=plane_axis)
        except (OSError, KeyError, ValueError, AssertionError) as e:
            print(f"Failed to parse {latest}: {e}")
            return False

        available_field_names = graph.field_names
        if current_field_name not in available_field_names:
            raise ValueError(
                f"Field '{current_field_name}' is not available in this file. Available fields: {available_field_names}",
            )
        current_field_index = available_field_names.index(current_field_name)

        # 提取位置和选定物理场（全为 numpy 数组）
        pos = graph.pos
        values = graph.x[:, current_field_index]
        triangles = graph.triangles

        R = pos[:, 0]
        Z = pos[:, 1]
        if triangles is not None:
            triang = tri.Triangulation(R, Z, triangles)
        else:
            print("Warning: No explicit triangles provided. Using Delaunay triangulation.")
            triang = tri.Triangulation(R, Z)

        # 更新图像
        ax.clear()
        contour_kwargs = {"cmap": cmap}

        if log_scale:
            positive_values = values[values > 0.0]
            if positive_values.size == 0:
                raise ValueError("Log scale requested but selected field has no positive values.")
            data_min = float(vmin) if vmin is not None else float(np.min(positive_values))
            data_max = float(vmax) if vmax is not None else float(np.max(positive_values))
            if data_min <= 0.0:
                raise ValueError("For log scale, vmin must be > 0.")
            if data_max <= data_min:
                raise ValueError("For log scale, vmax must be greater than vmin.")
            level_values = np.logspace(np.log10(data_min), np.log10(data_max), levels)
            contour_kwargs["levels"] = level_values
            contour_kwargs["norm"] = mcolors.LogNorm(vmin=data_min, vmax=data_max)
        else:
            if vmin is not None:
                contour_kwargs["vmin"] = float(vmin)
            if vmax is not None:
                contour_kwargs["vmax"] = float(vmax)
            contour_kwargs["levels"] = levels

        contour = ax.tricontourf(triang, values, **contour_kwargs)

        if cbar is not None:
            cbar.remove()
        cbar = fig.colorbar(contour, ax=ax)
        cbar.set_label(current_field_name)

        t_now = float(graph.time)
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        title = (
            f"{current_field_name} at t={t_now:.4e} "
            f"(file: {os.path.basename(latest)}, plane {current_plane}/{total_planes}"
        )
        if graph.phi is not None:
            title += f", phi={np.degrees(graph.phi):.1f} deg"
        title += ")"
        ax.set_title(title)
        ax.axis("equal")

        fig.canvas.draw()
        fig.canvas.flush_events()
        return True

    def on_key_press(event) -> None:
        nonlocal forced_file, last_file, last_mtime, current_field_name, force_refresh, current_field_index, current_plane, total_planes
        if event.key == "h":
            print_hotkey_help()
        elif event.key == "z":
            user_input = input("\nEnter file path to monitor (absolute or relative): ").strip()
            if not user_input:
                print("Keeping current monitoring target.")
                return
            candidate = os.path.abspath(user_input)
            if not os.path.exists(candidate):
                print(f"File does not exist: {candidate}")
                return
            forced_file = candidate
            last_file = None
            last_mtime = None
            print(f"Switched to fixed file: {forced_file}")
        elif event.key == "a":
            forced_file = None
            last_file = None
            last_mtime = None
            print(f"Switched to auto mode with pattern {pattern!r}")
        elif event.key in (",", "."):
            if total_planes <= 1:
                print(f"Only one plane available in the current file ({total_planes}).")
                return
            step = -1 if event.key == "," else 1
            current_plane = ((current_plane - 1 + step) % total_planes) + 1
            force_refresh = True
            print(
                f"Switched plane to {current_plane}/{total_planes} "
                "(reconstructed from toroidal harmonics). "
                "Will redraw immediately; phi is shown in the plot title.",
            )
        elif event.key in ("]", "["):
            cycle_list = available_field_names if available_field_names else FIELD_CYCLE
            cycle_idx = cycle_list.index(current_field_name)
            if event.key == "]":
                cycle_idx = (cycle_idx + 1) % len(cycle_list)
            else:
                cycle_idx = (cycle_idx - 1) % len(cycle_list)
            current_field_name = cycle_list[cycle_idx]
            current_field_index = None
            force_refresh = True
            print(
                f"Switched field to {current_field_name!r}. "
                "Will redraw immediately.",
            )

    fig.canvas.mpl_connect("key_press_event", on_key_press)

    print(
        f"Monitoring '{current_field_name}' in directory {directory!r} "
        f"with pattern {pattern!r}, interval={interval}s",
    )
    print(
        "操作指南: 'h' 帮助, ','/'.' 切换环向平面 (由 n_tor 个谐波按 Moivre 基重建), "
        "']'/'[' 切换物理场, 'z' 指定文件, 'a' 返回自动模式, Ctrl+C 退出。",
    )
    if forced_file:
        print(f"Start in fixed-file mode: {forced_file}")

    if playback:
        print(
            f"Playback mode enabled: interval={playback_interval}s, "
            f"loop={'on' if playback_loop else 'off'}",
        )
        if forced_file:
            print("Note: --file is ignored in playback mode; playback always uses --pattern.")

    try:
        while True:
            if playback:
                files = list_files_by_sequence(directory, pattern)
                if not files:
                    print("No matching HDF5 files found for playback; waiting...", end="\r", flush=True)
                    plt.pause(playback_interval)
                    continue
                for latest in files:
                    try:
                        render_file(latest)
                    except (OSError, KeyError, ValueError, AssertionError) as e:
                        print(f"Failed to render {latest}: {e}")
                    plt.pause(playback_interval)
                if not playback_loop:
                    print("Playback completed once from minimum to maximum sequence.")
                    break
                continue

            if forced_file is not None:
                latest = forced_file if os.path.exists(forced_file) else None
            else:
                latest = find_latest_file(directory, pattern)
            if latest is None:
                print("No matching HDF5 files found; waiting...", end="\r", flush=True)
                plt.pause(interval)
                continue

            try:
                mtime = os.path.getmtime(latest)
            except OSError:
                plt.pause(interval)
                continue

            # 仅当文件有变化时才重新读取
            if latest != last_file or mtime != last_mtime or force_refresh:
                try:
                    render_file(latest)
                except (OSError, KeyError, ValueError, AssertionError) as e:
                    print(f"Failed to render {latest}: {e}")
                    plt.pause(interval)
                    last_file = latest
                    last_mtime = mtime
                    continue

                last_file = latest
                last_mtime = mtime
                force_refresh = False

            plt.pause(interval)
    except KeyboardInterrupt:
        print("\nMonitoring stopped by user.")
    finally:
        # 关闭交互模式但保留当前图像
        plt.ioff()
        plt.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Monitor JOREK HDF5 fields in real time",
    )
    parser.add_argument(
        "--dir",
        default=".",
        help=(
            "Directory containing jorek HDF5 files "
            "(default: current directory)"
        ),
    )
    parser.add_argument(
        "--pattern",
        default="jorek0*",
        help=(
            "Glob pattern for HDF5 files "
            "(default: jorek0*)"
        ),
    )
    parser.add_argument(
        "--file",
        default=None,
        help=(
            "Specific HDF5 file to monitor. If set, starts in fixed-file mode; "
            "press 'a' to return to auto pattern mode."
        ),
    )
    parser.add_argument(
        "--field",
        default="psi",
        help=(
            "Field to visualize: one of "
            + ", ".join(CANONICAL_FIELD_ORDER)
        ),
    )
    parser.add_argument(
        "--plane",
        type=int,
        default=1,
        help=(
            "1-based physical toroidal plane to visualize (default: 1); "
            "planes are reconstructed from the toroidal harmonics stored in the HDF5 file"
        ),
    )
    parser.add_argument(
        "--plane-axis",
        type=int,
        choices=(1, 2),
        default=None,
        help=(
            "Deprecated: kept only for backward compatibility. The parser now reads "
            "the canonical JOREK restart layout directly and ignores this option."
        ),
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=20.0,
        help="Refresh interval in seconds (default: 20.0)",
    )
    parser.add_argument(
        "--cmap",
        default="jet",
        help="Matplotlib colormap name (default: jet)",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Optional minimum value for color range",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Optional maximum value for color range",
    )
    parser.add_argument(
        "--levels",
        type=int,
        default=50,
        help="Number of contour levels (default: 50)",
    )
    parser.add_argument(
        "--log-scale",
        action="store_true",
        help="Use logarithmic color scaling with log-spaced contour levels",
    )
    parser.add_argument(
        "--playback",
        action="store_true",
        help=(
            "Playback mode: iterate files from minimum to maximum sequence once "
            "using --pattern (ignores --file)"
        ),
    )
    parser.add_argument(
        "--playback-interval",
        type=float,
        default=0.3,
        help="Frame interval in seconds for playback mode (default: 0.3)",
    )
    parser.add_argument(
        "--playback-loop",
        action="store_true",
        help="Loop playback continuously after reaching the maximum sequence",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    directory = os.path.abspath(args.dir)
    selected_file = os.path.abspath(args.file) if args.file else None
    monitor_field(
        directory,
        args.field,
        args.interval,
        args.cmap,
        args.pattern,
        selected_file,
        plane=args.plane,
        plane_axis=args.plane_axis,
        vmin=args.vmin,
        vmax=args.vmax,
        levels=args.levels,
        log_scale=args.log_scale,
        playback=args.playback,
        playback_interval=args.playback_interval,
        playback_loop=args.playback_loop,
    )


if __name__ == "__main__":
    main()
