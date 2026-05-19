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
        --pattern "jorek_restart.h5"

默认只监看 jorek_restart.h5，每 60 秒刷新一次；
也可以通过 --pattern 选择其他文件模式（例如 "jorek*.h5" 会始终选最新的文件）。
支持的场名：psi, phi, omega, rho, T, v_par, T_e
"""

import argparse
import glob
import os
import re
import time
from dataclasses import dataclass
from typing import List, Optional

import h5py
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np


@dataclass
class GraphData:
    """简化版数据结构，仅包含可视化所需信息（全部为 numpy 数组）。"""

    x: np.ndarray          # [n_nodes, 7] 节点物理量 (psi, Phi, omega, rho, T, v_par, T_e)
    pos: np.ndarray        # [n_nodes, 2] 节点坐标 (R, Z)
    triangles: Optional[np.ndarray]  # [n_tri, 3] 三角形拓扑 (0-based)，或 None
    time: float            # 当前物理时间 t_now


class JorekDataParser:
    """解析 JOREK HDF5 输出文件，提取网格与物理量（纯 numpy 版本）。"""

    def __init__(self, h5_filepath: str) -> None:
        self.filepath = h5_filepath
        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f"Cannot find HDF5 file: {self.filepath}")

    def parse_to_graph(self) -> GraphData:
        """解析 HDF5 并返回 GraphData。"""
        with h5py.File(self.filepath, "r") as f:
            # 1. 坐标 pos (R,Z)
            if "x" not in f:
                raise KeyError("Coordinate array 'x' not found in HDF5.")
            x_data = f["x"]
            if len(x_data.shape) == 4:
                coords_raw = x_data[:, 0, 0, :]
            elif len(x_data.shape) == 3:
                coords_raw = x_data[:, 0, :]
            elif len(x_data.shape) == 2:
                coords_raw = x_data[:, :]
            else:
                raise ValueError(f"Unexpected shape for 'x': {x_data.shape}")
            # coords_raw 形状 (2, n_nodes) -> 转为 (n_nodes, 2)
            pos = np.array(coords_raw, dtype=np.float32).T.copy()

            # 2. 物理量 values -> 节点特征 x
            if "values" not in f:
                raise KeyError("Physical values array 'values' not found.")
            values_data = f["values"]
            if len(values_data.shape) == 4:
                values_raw = values_data[:, 0, 0, :]
            elif len(values_data.shape) == 3:
                values_raw = values_data[:, 0, :]
            elif len(values_data.shape) == 2:
                values_raw = values_data[:, :]
            else:
                raise ValueError(f"Unexpected shape for 'values': {values_data.shape}")

            # 按照约定提取 7 个变量: psi, Phi, omega, rho, T (Ti), v_par, T_e
            if values_raw.shape[0] <= 8:
                raise ValueError(
                    "Channel index 8 (T_e) is not available in 'values'. "
                    f"Current channel count: {values_raw.shape[0]}",
                )
            psi = np.array(values_raw[2, :], dtype=np.float32)
            Phi = np.array(values_raw[1, :], dtype=np.float32)
            omega = np.array(values_raw[3, :], dtype=np.float32)
            rho = np.array(values_raw[4, :], dtype=np.float32)
            T = np.array(values_raw[5, :], dtype=np.float32)
            v_par = np.array(values_raw[6, :], dtype=np.float32)
            T_e = np.array(values_raw[8, :], dtype=np.float32)
            x = np.stack([psi, Phi, omega, rho, T, v_par, T_e], axis=1)  # [n_nodes, 7]

            # 3. 时间戳
            if "t_now" in f:
                t_val = f["t_now"][()]
                if isinstance(t_val, np.ndarray):
                    if t_val.size == 1:
                        t_val = float(t_val.reshape(-1)[0])
                    else:
                        t_val = float(t_val.reshape(-1)[0])
                else:
                    t_val = float(t_val)
            else:
                t_val = 0.0

            # 4. 三角形拓扑（由 vertex 拆分为两个三角形）
            triangles = None
            if "vertex" in f:
                vertex = np.array(f["vertex"][...])  # [4, n_elements]
                if vertex.ndim != 2 or vertex.shape[0] < 3:
                    raise ValueError(f"Unexpected shape for 'vertex': {vertex.shape}")
                valid_mask = (vertex > 0).all(axis=0)
                valid_vertex = vertex[:, valid_mask].astype(np.int64) - 1  # -> 0-based
                # 每个四边形拆分为两个三角形
                tri1 = valid_vertex[[0, 1, 2], :].T
                tri2 = valid_vertex[[0, 2, 3], :].T
                triangles = np.vstack([tri1, tri2])

            # 5. 简单检查
            if pos.shape[0] != x.shape[0]:
                raise ValueError(
                    f"pos node count {pos.shape[0]} != values node count {x.shape[0]}",
                )

            return GraphData(x=x, pos=pos, triangles=triangles, time=t_val)


FIELD_NAME_TO_INDEX = {
    "psi": 0,
    "phi": 1,  # 对应 Phi
    "omega": 2,
    "w": 2,
    "rho": 3,
    "t": 4,    # 温度 Ti
    "ti": 4,
    "temperature": 4,
    "v_par": 5,
    "vpar": 5,
    "t_e": 6,
    "te": 6,
    "electron_temperature": 6,
}

INDEX_TO_CANONICAL_FIELD = {
    0: "psi",
    1: "phi",
    2: "omega",
    3: "rho",
    4: "t",
    5: "v_par",
    6: "t_e",
}

FIELD_CYCLE = ["psi", "phi", "omega", "rho", "t", "v_par", "t_e"]


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
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    levels: int = 50,
    log_scale: bool = False,
    playback: bool = False,
    playback_interval: float = 0.3,
    playback_loop: bool = False,
) -> None:
    """循环监看目录中最新匹配 pattern 的 JOREK HDF5 文件的指定物理场。"""
    fname_key = field_name.lower()
    if fname_key not in FIELD_NAME_TO_INDEX:
        raise ValueError(
            f"Unsupported field '{field_name}'. Supported: {sorted(FIELD_NAME_TO_INDEX.keys())}",
        )
    field_index = FIELD_NAME_TO_INDEX[fname_key]
    current_field_name = INDEX_TO_CANONICAL_FIELD[field_index]

    plt.ion()
    fig, ax = plt.subplots(figsize=(8, 10))
    cbar = None

    forced_file: Optional[str] = None
    if selected_file:
        forced_file = selected_file

    last_file: Optional[str] = None
    last_mtime: Optional[float] = None
    force_refresh = False

    def print_hotkey_help() -> None:
        mode = "playback" if playback else ("fixed-file" if forced_file else "auto-pattern")
        print("\n=== Hotkeys ===")
        print("  h : show this help and current status")
        print("  z : switch to a specific file")
        print("  a : return to auto pattern mode")
        print("  ] : next field")
        print("  [ : previous field")
        print("\n=== Current Status ===")
        print(f"  mode             : {mode}")
        print(f"  directory        : {directory}")
        print(f"  pattern          : {pattern}")
        print(f"  fixed file       : {forced_file if forced_file else '(none)'}")
        print(f"  field            : {current_field_name}")
        print(f"  colormap         : {cmap}")
        print(f"  levels           : {levels}")
        print(f"  log scale        : {log_scale}")
        print(f"  vmin             : {vmin}")
        print(f"  vmax             : {vmax}")
        print(f"  monitor interval : {interval}s")
        print(f"  playback interval: {playback_interval}s")
        print(f"  playback loop    : {playback_loop}")

    def render_file(latest: str) -> bool:
        nonlocal cbar
        print(f"\nReading {os.path.basename(latest)}", flush=True)
        try:
            parser = JorekDataParser(latest)
            graph = parser.parse_to_graph()
        except (OSError, KeyError, ValueError, AssertionError) as e:
            print(f"Failed to parse {latest}: {e}")
            return False

        # 提取位置和选定物理场（全为 numpy 数组）
        pos = graph.pos
        values = graph.x[:, field_index]
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
        ax.set_title(
            f"{current_field_name} at t={t_now:.4e} "
            f"(file: {os.path.basename(latest)})",
        )
        ax.axis("equal")

        fig.canvas.draw()
        fig.canvas.flush_events()
        return True

    def on_key_press(event) -> None:
        nonlocal forced_file, last_file, last_mtime, field_index, current_field_name, force_refresh
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
        elif event.key in ("]", "["):
            cycle_idx = FIELD_CYCLE.index(current_field_name)
            if event.key == "]":
                cycle_idx = (cycle_idx + 1) % len(FIELD_CYCLE)
            else:
                cycle_idx = (cycle_idx - 1) % len(FIELD_CYCLE)
            current_field_name = FIELD_CYCLE[cycle_idx]
            field_index = FIELD_NAME_TO_INDEX[current_field_name]
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
        "Hotkeys: 'h' help, 'z' switch to a specific file, 'a' return to auto pattern mode, "
        "']' next field, '[' previous field.",
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
            + ", ".join(sorted(FIELD_NAME_TO_INDEX.keys()))
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
