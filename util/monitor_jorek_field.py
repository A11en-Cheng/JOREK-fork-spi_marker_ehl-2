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
支持的场名：psi, phi, rho, T, v_par
"""

import argparse
import glob
import os
import time
from dataclasses import dataclass
from typing import Optional

import h5py
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np


@dataclass
class GraphData:
    """简化版数据结构，仅包含可视化所需信息（全部为 numpy 数组）。"""

    x: np.ndarray          # [n_nodes, 5] 节点物理量 (psi, Phi, rho, T, v_par)
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

            # 按照约定提取 5 个变量: psi, Phi, rho, T (Ti), v_par
            psi = np.array(values_raw[2, :], dtype=np.float32)
            Phi = np.array(values_raw[1, :], dtype=np.float32)
            rho = np.array(values_raw[4, :], dtype=np.float32)
            T = np.array(values_raw[5, :], dtype=np.float32)
            v_par = np.array(values_raw[6, :], dtype=np.float32)
            x = np.stack([psi, Phi, rho, T, v_par], axis=1)  # [n_nodes, 5]

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
    "rho": 2,
    "t": 3,    # 温度 Ti
    "ti": 3,
    "temperature": 3,
    "v_par": 4,
    "vpar": 4,
}


def find_latest_file(directory: str, pattern: str) -> Optional[str]:
    """在给定目录中查找匹配 pattern 的最新文件（按修改时间排序）。"""
    paths = glob.glob(os.path.join(directory, pattern))
    if not paths:
        return None
    # 选取修改时间最新的文件
    return max(paths, key=os.path.getmtime)


def monitor_field(
    directory: str,
    field_name: str,
    interval: float = 60.0,
    cmap: str = "jet",
    pattern: str = "jorek_restart.h5",
) -> None:
    """循环监看目录中最新匹配 pattern 的 JOREK HDF5 文件的指定物理场。"""
    fname_key = field_name.lower()
    if fname_key not in FIELD_NAME_TO_INDEX:
        raise ValueError(
            f"Unsupported field '{field_name}'. Supported: {sorted(FIELD_NAME_TO_INDEX.keys())}",
        )
    field_index = FIELD_NAME_TO_INDEX[fname_key]

    plt.ion()
    fig, ax = plt.subplots(figsize=(8, 10))
    cbar: Optional[plt.colorbar] = None

    last_file: Optional[str] = None
    last_mtime: Optional[float] = None

    print(
        f"Monitoring '{field_name}' in directory {directory!r} "
        f"with pattern {pattern!r}, interval={interval}s",
    )
    try:
        while True:
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
            if latest != last_file or mtime != last_mtime:
                print(f"\nReading {os.path.basename(latest)}", flush=True)
                try:
                    parser = JorekDataParser(latest)
                    graph = parser.parse_to_graph()
                except (OSError, KeyError, ValueError, AssertionError) as e:
                    print(f"Failed to parse {latest}: {e}")
                    plt.pause(interval)
                    last_file = latest
                    last_mtime = mtime
                    continue

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
                contour = ax.tricontourf(triang, values, levels=50, cmap=cmap)

                if cbar is not None:
                    cbar.remove()
                cbar = fig.colorbar(contour, ax=ax)
                cbar.set_label(field_name)

                t_now = float(graph.time)
                ax.set_xlabel("R [m]")
                ax.set_ylabel("Z [m]")
                ax.set_title(
                    f"{field_name} at t={t_now:.4e} "
                    f"(file: {os.path.basename(latest)})",
                )
                ax.axis("equal")

                fig.canvas.draw()
                fig.canvas.flush_events()

                last_file = latest
                last_mtime = mtime

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
        default="jorek_restart.h5",
        help=(
            "Glob pattern for HDF5 files "
            "(default: jorek_restart.h5)"
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
        default=60.0,
        help="Refresh interval in seconds (default: 60.0)",
    )
    parser.add_argument(
        "--cmap",
        default="jet",
        help="Matplotlib colormap name (default: jet)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    directory = os.path.abspath(args.dir)
    monitor_field(directory, args.field, args.interval, args.cmap, args.pattern)


if __name__ == "__main__":
    main()
