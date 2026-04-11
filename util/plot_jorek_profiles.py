#!/usr/bin/env python3
"""Plot JOREK density/temperature derived profiles.

Reads two-column profile files (psi, value), computes:
- pressure profile: p(psi) = rho(psi) * T(psi)
- first derivative of density with respect to psi: d rho / d psi

Usage example:
  python3 tools/plot_jorek_profiles.py \
    --rho jorek_density \
    --temperature jorek_temperature \
    --output rho_T_p_profiles.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np


def read_two_column_profile(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Read a 2-column ASCII profile file.

    Supports comma-separated and whitespace-separated formats.
    Lines beginning with '#' or '!' are ignored.
    """
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    psi = []
    val = []

    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#") or line.startswith("!"):
                continue

            # Accept both CSV and whitespace separators.
            line = line.replace(",", " ")
            parts = line.split()
            if len(parts) < 2:
                continue

            try:
                x = float(parts[0])
                y = float(parts[1])
            except ValueError:
                continue

            psi.append(x)
            val.append(y)

    if not psi:
        raise ValueError(f"No valid profile data found in: {path}")

    psi_arr = np.asarray(psi, dtype=float)
    val_arr = np.asarray(val, dtype=float)

    # Sort by psi and remove duplicate psi by keeping last occurrence.
    order = np.argsort(psi_arr)
    psi_sorted = psi_arr[order]
    val_sorted = val_arr[order]

    unique_psi, unique_indices = np.unique(psi_sorted, return_index=True)
    # np.unique keeps first index; choose first is fine after sorting.
    unique_val = val_sorted[unique_indices]

    return unique_psi, unique_val


def common_grid(
    psi_rho: np.ndarray,
    rho: np.ndarray,
    psi_t: np.ndarray,
    t: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map rho and T onto a common psi grid using linear interpolation."""
    psi_min = max(float(np.min(psi_rho)), float(np.min(psi_t)))
    psi_max = min(float(np.max(psi_rho)), float(np.max(psi_t)))

    if psi_min >= psi_max:
        raise ValueError("Density and temperature psi ranges do not overlap.")

    # Use the denser source as base grid inside overlap region.
    base = psi_rho if psi_rho.size >= psi_t.size else psi_t
    mask = (base >= psi_min) & (base <= psi_max)
    psi = base[mask]

    if psi.size < 3:
        # Fallback to a uniform grid if overlap has too few points.
        psi = np.linspace(psi_min, psi_max, 256)

    rho_i = np.interp(psi, psi_rho, rho)
    t_i = np.interp(psi, psi_t, t)
    return psi, rho_i, t_i


def plot_profiles(
    psi: np.ndarray,
    rho: np.ndarray,
    temp: np.ndarray,
    pressure: np.ndarray,
    drho_dpsi: np.ndarray,
    output: Path,
    show: bool,
) -> None:
    """Create and save the plot figure."""
    fig, axes = plt.subplots(2, 1, figsize=(9, 9), sharex=True)

    axes[0].plot(psi, rho, label="rho", lw=2)
    axes[0].plot(psi, temp, label="T", lw=2)
    axes[0].plot(psi, pressure, label="p = rho*T", lw=2)
    axes[0].set_ylabel("Profile value")
    axes[0].set_title("JOREK profiles")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(psi, drho_dpsi, label="d(rho)/d(psi)", lw=2, color="tab:red")
    axes[1].axhline(0.0, color="k", lw=1, alpha=0.5)
    axes[1].set_xlabel("psi")
    axes[1].set_ylabel("Derivative")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(output, dpi=160)

    if show:
        plt.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read JOREK density/temperature profiles and plot rho, T, p=rho*T, and d(rho)/d(psi)."
    )
    parser.add_argument("--rho", required=True, help="Path to jorek_density file")
    parser.add_argument("--temperature", required=True, help="Path to jorek_temperature file")
    parser.add_argument(
        "--output",
        default="jorek_profiles.png",
        help="Output figure path (default: jorek_profiles.png)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display interactive plot window in addition to saving file",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    rho_path = Path(args.rho)
    t_path = Path(args.temperature)
    out_path = Path(args.output)

    psi_rho, rho = read_two_column_profile(rho_path)
    psi_t, temp = read_two_column_profile(t_path)

    psi, rho_i, temp_i = common_grid(psi_rho, rho, psi_t, temp)

    pressure = rho_i * temp_i
    drho_dpsi = np.gradient(rho_i, psi)

    plot_profiles(psi, rho_i, temp_i, pressure, drho_dpsi, out_path, args.show)

    print(f"Saved plot: {out_path}")
    print(f"psi range: [{psi.min():.6g}, {psi.max():.6g}], points: {psi.size}")


if __name__ == "__main__":
    main()
