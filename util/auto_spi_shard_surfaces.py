#!/usr/bin/env python3
"""
auto_spi_shard_surfaces.py

SPI 碎片"火箭效应"三维诊断：把每个碎片**看到的**温度与它**被推出来**的速度/消融
量画在同一张 (碎片 x 时间) 网格上。

默认（--layout merged）只出一张合并图：

    x = 碎片编号 i_spi
    y = 自消融开始的时间 [ms]
    z = |ΔV_R| = |V_R(t) - V_R(t_inj)|   （速度衰减幅度，单位 m/s）
    颜色 = 碎片在该位置"看到"的 T_e（--color te/tn/psin 可切）
    透明度 = r/r0（刚注入不透明，烧蚀完变淡）—— 半径信息不再单独占一张图

用 --layout split 可以退回两张分开的图（ΔV_R 与 Δr 各一张，颜色标度一致）。

关于"每个弹丸是一条线还是曲面"：在 (碎片编号, 时间) 网格上，每个碎片本来只是
一条沿时间方向的**线**；曲面是通过连接相邻碎片把线连成整片。所以：
  --render surface （默认）画网格曲面，相邻碎片被连成整片（视觉上像一排帘子）；
  --render lines   每个碎片画一条 3D 折线，语义上就是"一个弹丸一条线"；
  --render both    两者叠加。

x 轴排序 --x_axis 除默认的碎片编号 index 外，还可按初始半径 r0、初始 Psi_N psi0，
或按"半径置零"（烧蚀完成）的先后 death 排列：时间窗内越早烧蚀完越靠前；窗口末
仍未烧蚀完的碎片按末时刻半径升序排在其后。用于观察"烧蚀完成顺序 vs 火箭加速"的关联。

除 ΔV_R 合并图外，默认（--accel on）还会额外出"碎片加速度"合并图（同样的
x/y/颜色/透明度，只把 z 换成加速度）：
    z = |a_R| = |dV_R/dt|                       [m/s²]
它是 ΔV_R 图的"被积函数"（ΔV_R = ∫ a_R dt）：对与 ΔV_R 相同的速度序列按真实时间
做非等距中心差分。速度带 f12.3 的 1 mm/s 量化阶梯，微分会放大噪声，必要时用
--smooth 先平滑；--acc_sign {abs,minus,plus} 控制 z 的正负。

物理含义：ΔV_R = V_R(t) - V_R(t_inj) 严格等于模型施加的火箭加速度积分 ∫a_R dt；
Δr = r(t_inj) - r(t) 是消融掉掉的半径（m = 4/3 π r^3 ρ_atom 因此变小）。
两者合起来就是 ΔV_R = ∫ p_asymm(T_e,n_e) π r² / m dt 的被积函数分解。

数据来源（沿用 util/auto_spi_diagnostics.sh 的产出）
--------------------------------------------------
* fragments_position_<step>.csv      列: i_spi, R, Z, phi, radius, species
* fragments_plume_total_<step>.csv   列: R, -phi, Z, VR, -VRxZ, VZ, radius, X_D2
  两者都由 jorek2_diagno_spi 以 `do i_spi = 1, n_spi_tot` 顺序写出，因此**行号即碎片编号**，
  本脚本按行对齐索引（plume_total 文件本身没有 i_spi 列，i_spi 从 position 文件取）。
  注意 plume_total 存的是 -phi / -VRxZ（JOREK 手性取负），脚本内部已翻回。
* postproc/exprs_Rmin*_s<step>.h5 (jorek2_postproc rectangle) -> T_e / ne / Psi_N，用于颜色
* macroscopic_vars.dat + abl_history.dat -> 时间归一化与消融起点 t0

时间约定（与 util/track_largest_pellets.py 完全一致）
----------------------------------------------------
    t_si(step) = time_map[step] * time_scaling        # time_scaling 来自 @sqrt_mu0_rho0:
    t0         = abl_history.dat 中第一个 abl_rate>0 的前一行  (= 消融开始)
    y 轴       = (t_si - t0) * 1e3                    # [ms]，自消融开始计时

用法
----
  python3 util/auto_spi_shard_surfaces.py \\
      --csv_dir      auto_diagno_spi/csv_files \\
      --rectangle_dir auto_diagno_spi/postproc \\
      --abl_file     auto_diagno_spi/abl_history.dat \\
      --live_data    macroscopic_vars.dat \\
      --out_dir      auto_diagno_spi/surfaces

只看速度/半径（不插值 T_e）时省略 --rectangle_dir 并加 --color none。
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors

# --------------------------------------------------------------------------------------
# 常量 / 数据列定义（与 Fortran 的 write 语句一一对应）
# --------------------------------------------------------------------------------------
POSITION_COLS = ["i_spi", "spi_R", "spi_Z", "spi_phi", "spi_radius", "spi_species"]
# fragments_plume_total: write(20,'(6f12.3,e14.6,f12.3)') R, -phi, Z, VR, -VRxZ, VZ, radius, X_D2
PLUME_COLS = ["R", "neg_phi", "Z", "VR", "neg_VRxZ", "VZ", "radius", "X_D2"]

COLOR_CHOICES = ["te", "tn", "psin", "none"]
COLOR_EXPR = {"te": "T_e", "tn": "T_e*ne", "psin": "Psi_N"}
# 变量名容错：不同 postproc build / 表达式注册表里可能有细微差别
VAR_ALIASES = {
    "T_e": ["Te", "T_E", "T", "te", "T_e"],
    "ne": ["n_e", "ne_main", "nn_main"],
    "Psi_N": ["PsiN", "psi_N", "psi_n"],
}
RECT_BOUNDS_RE = re.compile(
    r"Rmin(?P<rmin>[-+0-9.eEdD]+)_Rmax(?P<rmax>[-+0-9.eEdD]+)"
    r"_Zmin(?P<zmin>[-+0-9.eEdD]+)_Zmax(?P<zmax>[-+0-9.eEdD]+)"
)
STEP_RE = re.compile(r"_s(\d{6})\.h5$")   # rectangle h5: ..._s000100.h5
CSV_STEP_RE = re.compile(r"_(\d{6})\.csv$")  # fragments_<...>_000100.csv


def log(msg: str) -> None:
    print(msg, flush=True)


# --------------------------------------------------------------------------------------
# 时间归一化（与 track_largest_pellets.py::load_time_mapping 同逻辑）
# --------------------------------------------------------------------------------------
def load_time_mapping(live_data_file: str):
    """解析 macroscopic_vars.dat -> (step -> t_jorek 映射, time_scaling)。

    time_scaling 来自 `@sqrt_mu0_rho0:`，即 1 个 JOREK 时间单位 = time_scaling 秒。
    """
    time_map: dict[int, float] = {}
    time_scaling = None
    with open(live_data_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("@sqrt_mu0_rho0:"):
                try:
                    time_scaling = float(line.strip().split()[1])
                except Exception:
                    pass
            elif line.startswith("@times:"):
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        time_map[int(parts[1])] = float(parts[2])
                    except ValueError:
                        continue
    if time_scaling is None:
        raise ValueError(
            f"FATAL: 在 {live_data_file} 中找不到 '@sqrt_mu0_rho0:'，无法做时间归一化。"
        )
    if not time_map:
        raise ValueError(f"FATAL: 在 {live_data_file} 中找不到任何 '@times:' 记录。")
    return time_map, time_scaling


def load_ablation_start(abl_file: str, time_scaling: float, start_time=None):
    """自 abl_history.dat 定位消融开始时刻 t0（SI 秒）。

    与 track_largest_pellets.py 一致：取第一个 abl_rate>0 记录的**前一行**作为 t0
    （即消融率还是 0、紧接着就起跳的那一步）。
    """
    if start_time is not None:
        log(f"[time] t0 由 --start_time 指定: {start_time:.6e} s")
        return float(start_time)

    if not os.path.exists(abl_file):
        log(f"[time] WARNING: 找不到 {abl_file}，无法自动定位消融起点，"
            f"将退化为使用碎片数据的第一个时间步。")
        return None

    # 注意：jorek2_diagno_spi 写的表头是
    #   write(20,'(A11)') 'time','abl_rate_inj','abl_number_inj',...
    # 格式 (A11) 会被复用 -> **每个表头名字各占一行**，所以不能直接交给 pandas
    # （它会拿第一行推断列数）。这里逐行解析，跳过所有非数值行。
    steps_list, times_list, rates_list = [], [], []
    n_skipped = 0
    with open(abl_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            toks = line.split()
            if len(toks) < 3:
                n_skipped += 1
                continue
            try:
                step_v = int(float(toks[0]))
                time_v = float(toks[1])
                rate_v = float(toks[2])       # 第 1 个注入器的消融率
            except ValueError:
                n_skipped += 1                # 表头行
                continue
            steps_list.append(step_v)
            times_list.append(time_v)
            rates_list.append(rate_v)

    if not steps_list:
        log(f"[time] WARNING: {abl_file} 中找不到数值行，退化到碎片数据首步。")
        return None

    step_arr = np.asarray(steps_list, dtype=int)
    time_si = np.asarray(times_list, dtype=float) * time_scaling
    rate_arr = np.asarray(rates_list, dtype=float)
    log(f"[time] abl_history.dat: 解析出 {len(steps_list)} 行数值记录"
        f"（跳过 {n_skipped} 行表头/空行）")

    active = np.nonzero(rate_arr > 0.0)[0]
    if active.size == 0:
        log(f"[time] WARNING: {abl_file} 中没有任何 abl_rate>0 的记录（注入器 1），"
            f"退化到碎片数据首步。")
        return None

    first = int(active.min())
    idx = max(0, first - 1)          # 往前推一步，保留消融率刚等于 0 的时刻
    t0 = float(time_si[idx])
    log(f"[time] 自动识别消融起点: step = {int(step_arr[idx])}, "
        f"t0 = {t0:.6e} s (abl_rate = {rate_arr[idx]:.3e} -> {rate_arr[first]:.3e})")
    return t0


# --------------------------------------------------------------------------------------
# 碎片 CSV 读取（按行对齐索引）
# --------------------------------------------------------------------------------------
def collect_steps(csv_dir: str):
    """返回同时存在 position 与 plume_total 两个文件的时间步列表（升序）。"""
    pos = {}
    for f in glob.glob(os.path.join(csv_dir, "fragments_position_*.csv")):
        m = CSV_STEP_RE.search(os.path.basename(f))
        if m and re.match(r"fragments_position_\d+\.csv$", os.path.basename(f)):
            pos[int(m.group(1))] = f
    plum = {}
    for f in glob.glob(os.path.join(csv_dir, "fragments_plume_total_*.csv")):
        m = CSV_STEP_RE.search(os.path.basename(f))
        if m and re.match(r"fragments_plume_total_\d+\.csv$", os.path.basename(f)):
            plum[int(m.group(1))] = f
    common = sorted(set(pos) & set(plum))
    missing_plume = sorted(set(pos) - set(plum))
    if missing_plume:
        log(f"[csv] WARNING: {len(missing_plume)} 个 step 缺少 fragments_plume_total_*.csv，"
            f"已跳过（含速度的文件才会被使用）。")
    return common, pos, plum


def _read_columns(path: str, ncol: int, ncol_expected: int, tag: str):
    """读取空白分隔的 CSV -> (nrow, ncol) 的 float 数组。

    两个层次的容错：
    1) 常规路径用 pandas 快速解析，并用 to_numeric(errors='coerce') 把无法解析的
       字段（例如 Fortran 固定宽度溢出写出的星号）置为 NaN；
    2) 若出现**字段粘连**（某个值超出固定宽度、把相邻两个数连在一起，导致某些行的
       token 数与其它行不一致），则退化为逐行解析：行数保持不变（保证"行号=碎片编号"
       的对齐关系），只把出问题的那一行整行置为 NaN，并给出明确警告。
    """
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    if df.empty or df.shape[1] < ncol_expected:
        log(f"[csv] WARNING: {tag} {os.path.basename(path)} 列数不足 "
            f"({df.shape[1]} < {ncol_expected})，跳过。")
        return None

    if df.shape[1] == ncol:
        arr = df.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        if np.isnan(arr).any():
            log(f"[csv] WARNING: {tag} {os.path.basename(path)} 含无法解析的字段，已置为 NaN。")
        return arr

    # --- 字段粘连的兜底路径 ---
    log(f"[csv] WARNING: {tag} {os.path.basename(path)} 每行列数不一致 "
        f"(多数为 {df.shape[1]}，期望 {ncol})——通常意味着某个数值超出了 Fortran 固定宽度。"
        f"将逐行解析，并把异常行整行置为 NaN。")
    raw_lines = [ln for ln in open(path, "r", encoding="utf-8", errors="ignore")
                 if ln.strip()]
    arr = np.full((len(raw_lines), ncol), np.nan)
    n_bad = 0
    for i, ln in enumerate(raw_lines):
        toks = ln.split()
        if len(toks) != ncol:
            n_bad += 1
            continue
        try:
            arr[i, :] = [float(t) for t in toks]
        except ValueError:
            n_bad += 1
    if n_bad:
        log(f"[csv] WARNING: {tag} {os.path.basename(path)}: {n_bad}/{len(raw_lines)} "
            f"行无法解析（已置 NaN）。若数量很大，请检查该时间步的 f12.3 字段是否溢出。")
    return arr


def load_shard_time_series(csv_dir: str, steps, debug: bool = False):
    """把逐 step 的碎片表堆成 (n_shard, n_step) 的时间序列。

    Returns
    -------
    dict with keys:
        steps      (n_step,)      时间步数组
        i_spi      (n_shard,)     碎片编号（取自 position 文件第 0 列）
        R, Z, phi, radius, species, VR, VRxZ, VZ  : (n_shard, n_step)
    """
    n_step = len(steps)
    # 第一遍：确定最大碎片数，并按行对齐
    per_step = []
    n_shard_max = 0
    for step in steps:
        pos_path = os.path.join(csv_dir, f"fragments_position_{step:06d}.csv")
        plum_path = os.path.join(csv_dir, f"fragments_plume_total_{step:06d}.csv")
        pos = _read_columns(pos_path, len(POSITION_COLS), 6, "position")
        plum = _read_columns(plum_path, len(PLUME_COLS), 8, "plume_total")
        if pos is None or plum is None:
            per_step.append(None)
            continue
        n_pos, n_plum = pos.shape[0], plum.shape[0]
        if n_pos != n_plum:
            log(f"[csv] WARNING: step {step} 两个文件行数不一致 "
                f"(position {n_pos} vs plume_total {n_plum})，按较小者对齐。")
            n = min(n_pos, n_plum)
            pos, plum = pos[:n], plum[:n]
        per_step.append((pos, plum))
        n_shard_max = max(n_shard_max, pos.shape[0])

    if n_shard_max == 0:
        raise RuntimeError(f"{csv_dir} 中没有可用的碎片数据。")

    # 第二遍：填充数组（不足的碎片补 NaN，便于后续 mask）
    keys_f = ["R", "Z", "phi", "radius", "species", "VR", "VRxZ", "VZ", "X_D2"]
    out = {k: np.full((n_shard_max, n_step), np.nan) for k in keys_f}
    i_spi = np.arange(1, n_shard_max + 1, dtype=int)

    for j, item in enumerate(per_step):
        if item is None:
            continue
        pos, plum = item
        n = pos.shape[0]
        i_spi[:n] = pos[:, 0].astype(int)
        out["R"][:n, j] = pos[:, 1]
        out["Z"][:n, j] = pos[:, 2]
        out["phi"][:n, j] = pos[:, 3]
        out["radius"][:n, j] = pos[:, 4]
        out["species"][:n, j] = pos[:, 5]
        # plume_total: R, -phi, Z, VR, -VRxZ, VZ, radius, X_D2
        out["VR"][:n, j] = plum[:, 3]
        out["VRxZ"][:n, j] = -plum[:, 4]      # 翻回 JOREK 手性
        out["VZ"][:n, j] = plum[:, 5]
        out["X_D2"][:n, j] = plum[:, 7]
        # radius 以 position 文件为准（两者都是 e14.6），交叉检查
        if debug:
            dmax = np.nanmax(np.abs(plum[:, 6] - pos[:, 4]))
            if dmax > 1e-9:
                log(f"[csv][debug] step {step}: radius 字段最大差 {dmax:.3e}")

    out["steps"] = np.asarray(steps, dtype=int)
    out["i_spi"] = i_spi
    return out


# --------------------------------------------------------------------------------------
# rectangle HDF5（jorek2_postproc rectangle）读取与插值
# --------------------------------------------------------------------------------------
class RectangleFields:
    """读取 postproc `rectangle` 导出的 HDF5，并按需在 (R,Z) 上插值。

    文件内部结构（见 diagnostics/new_diag/mod_diag_output.f90::write_hdf5_2d）:
        variables : 变量名列表,  n_var : 变量个数
        dim       : (nR, nZ)
        values    : 场数据（维度顺序取决于 HDF5 写入方式，本类自动识别）
        time/index_now : SI 时间与 step
    坐标轴不写入文件，但 R/Z 范围可从文件名 (Rmin..._Rmax..._Zmin..._Zmax...) 精确解析。
    """

    def __init__(self, rectangle_dir: str, bounds=None, debug: bool = False,
                 expr_order=None):
        self.dir = rectangle_dir
        self.bounds = bounds
        self.debug = debug
        # postproc 写出 h5 时把变量名数组声明成了 character(len=1)
        # (mod_diag_output.f90: `character(len=1), allocatable :: varnames(:)`),
        # 所以文件里的 'variables' 可能只剩每个名字的首字母（例如 P/T/n/n）。
        # 这时用 --expr_order（= 生成 h5 时 postproc 的 expressions 顺序）做确定性映射。
        self.expr_order = [s for s in re.split(r"[,\s]+", expr_order) if s] if expr_order else None
        self.step_map: dict[int, str] = {}
        for f in sorted(glob.glob(os.path.join(rectangle_dir, "exprs_Rmin*.h5"))):
            m = STEP_RE.search(os.path.basename(f))
            if m:
                self.step_map[int(m.group(1))] = f
        self._var_cache: dict[str, dict] = {}   # path -> {'names': [...], 'index': {...}}
        self._interp_cache: dict[tuple, object] = {}
        self._noted: set = set()                # 已打印过的提示（避免刷屏）
        # 插值器缓存上限：每步每变量一份 200x200 网格，不限制会吃掉几个 GB
        self.max_interp_cache = 6

    # -- helpers ------------------------------------------------------------------
    @staticmethod
    def _bounds_from_name(path: str):
        m = RECT_BOUNDS_RE.search(os.path.basename(path))
        if not m:
            return None
        vals = []
        for k in ("rmin", "rmax", "zmin", "zmax"):
            s = m.group(k).replace("D", "E").replace("d", "e")
            vals.append(float(s))
        return tuple(vals)

    @staticmethod
    def _decode_varnames(raw, n_var=None):
        """把 h5py 读到的 'variables' 解成名字列表。

        不同写盘方式会让 h5py 返回不同形状，这里全部兼容——否则会出现
        "按名字查不到变量"的假故障（原 track_largest_pellets.py 从不按名字查，
        所以一直没暴露；它用 values[0] 恰好就是 Psi_N）：
          * 1D 定长字符串数组      -> 每个元素一个名字（最常见）
          * 2D 字符数组 (n_var, w) -> 每行拼起来是一个名字
          * 1D 单字符数组 (n_var*w,) -> 按 n_var 分块后再拼
        """
        def _elem(x):
            if isinstance(x, (bytes, bytearray)):
                return x.decode("utf-8", "ignore")
            return str(x)

        arr = np.asarray(raw)
        if arr.size == 0:
            return []
        if arr.ndim >= 2:
            names = ["".join(_elem(x) for x in row) for row in arr]
        else:
            flat = [_elem(x) for x in arr.ravel()]
            if n_var and n_var > 0 and len(flat) % n_var == 0 and all(len(s) <= 1 for s in flat):
                w = len(flat) // n_var
                names = ["".join(flat[i * w:(i + 1) * w]) for i in range(n_var)]
            else:
                names = flat
        return [n.replace("\x00", "").strip() for n in names]

    def _var_index(self, path: str):
        """读取该文件的变量元信息（名字 + 维度），带缓存。"""
        if path in self._var_cache:
            return self._var_cache[path]
        import h5py
        with h5py.File(path, "r") as f:
            raw = np.asarray(f["variables"][()]) if "variables" in f else np.array([])
            n_var = int(f["n_var"][()]) if "n_var" in f else None
            shape = tuple(int(s) for s in f["values"].shape) if "values" in f else None
            dim = tuple(int(d) for d in np.atleast_1d(f["dim"][()])) if "dim" in f else None
        names = self._decode_varnames(raw, n_var=n_var)
        info = {
            "names": names,
            "index": {n: i for i, n in enumerate(names)},
            "n_var": n_var,
            "values_shape": shape,
            "dim": dim,
            "raw_dtype": str(raw.dtype),
            "raw_shape": tuple(raw.shape),
        }
        self._var_cache[path] = info
        return info

    def describe(self, path: str) -> str:
        """一行描述文件里的变量清单，排错时直接打印它。"""
        info = self._var_index(path)
        degenerate = (info["names"] and all(len(n) <= 1 for n in info["names"]))
        hint = ("  <-- 名字退化成首字母（postproc 的 character(len=1) 写盘所致），"
                "请用 --expr_order 指定表达式顺序" if degenerate else "")
        return (f"{os.path.basename(path)}:\n"
                f"           variables = {info['names']}{hint}\n"
                f"           n_var = {info['n_var']}, values.shape = {info['values_shape']}, "
                f"dim = {info['dim']}, raw variables dtype/shape = "
                f"{info['raw_dtype']}{info['raw_shape']}")

    def resolve(self, path: str, var_name: str):
        """把请求的名字解析成文件里实际存在的名字 -> (actual_name, index)。

        三条路径，按可信度递减：
          1) 名字（含别名）在 h5 的 'variables' 里直接命中；
          2) 名字数组退化成首字母（Fortran character(len=1) 的锅）时，
             用 --expr_order 给出的表达式顺序做位置映射；
          3) 都不行 -> 返回 (None, None)，调用方打印变量清单并提示用
             --color_var / --color_index。
        """
        info = self._var_index(path)
        for cand in [var_name] + VAR_ALIASES.get(var_name, []):
            if cand in info["index"]:
                return cand, info["index"][cand]

        if self.expr_order and var_name in self.expr_order:
            idx = self.expr_order.index(var_name)
            if idx < (info["n_var"] or len(info["names"])):
                self._note_once(
                    f"expr_order:{var_name}",
                    f"[h5] NOTE: 文件里的 variables={info['names']} 无法按名字匹配，"
                    f"改用 --expr_order 的位置映射: '{var_name}' -> index {idx}。")
                return var_name, idx
        return None, None

    def _note_once(self, key: str, msg: str):
        """同一条提示只打印一次（否则每个 step/变量都会刷屏）。"""
        if key not in self._noted:
            self._noted.add(key)
            log(msg)

    def nearest_step(self, step: int, max_gap: int = 10):
        """找最近的可用的 rectangle step。"""
        if not self.step_map:
            return None
        if step in self.step_map:
            return step
        avail = sorted(self.step_map.keys())
        idx = int(np.searchsorted(avail, step))
        cands = []
        if idx < len(avail):
            cands.append(avail[idx])
        if idx > 0:
            cands.append(avail[idx - 1])
        if not cands:
            return None
        best = min(cands, key=lambda s: abs(s - step))
        if abs(best - step) > max_gap:
            log(f"[h5] WARNING: step {step} 无 rectangle 文件，最近为 {best} "
                f"(Δ={abs(best - step)})，插值结果可能不准确。")
        return best

    def read_plane(self, path: str, var_name: str = None, index: int = None):
        """读取单个变量在整个矩形区域上的 2D 场，返回 (plane, (Rmin,Rmax,Zmin,Zmax))。

        plane 统一为 shape (nR, nZ)。变量可用名字（含别名）或直接给 0-based 序号。
        """
        import h5py
        info = self._var_index(path)
        if index is None:
            actual, index = self.resolve(path, var_name)
            if actual is None:
                log(f"[h5] WARNING: 变量 '{var_name}' 不存在。文件内容：\n"
                    f"           {self.describe(path)}\n"
                    f"           -> 可用 --color_var <实际名字> 或 --color_index <0-based 序号>；"
                    f"用 --probe_h5 可单独查看变量清单。")
                return None, None
            if actual != var_name:
                self._note_once(f"alias:{var_name}->{actual}",
                                f"[h5] NOTE: 变量 '{var_name}' 在文件中写作 '{actual}'，"
                                f"已自动匹配。")
        i_var = int(index)

        with h5py.File(path, "r") as f:
            values = f["values"][()]
            dim = tuple(int(d) for d in f["dim"][()]) if "dim" in f else None

        if values.ndim == 2:
            plane = values
        elif values.ndim == 3:
            # 找到长度等于 n_var 的那个轴作为变量轴（n_var 优先取数据集，其次用名字个数）
            n_var = info["n_var"] or len(info["names"])
            var_axes = [ax for ax in range(3) if values.shape[ax] == n_var]
            if not var_axes:
                log(f"[h5] WARNING: {os.path.basename(path)} values.shape={values.shape} "
                    f"与 n_var={n_var} 不匹配，假定第 0 轴为变量轴。")
                axis = 0
            else:
                axis = var_axes[0]
            if i_var >= values.shape[axis]:
                log(f"[h5] WARNING: 变量序号 {i_var} 超出范围 "
                    f"(该轴长度 {values.shape[axis]})。")
                return None, None
            plane = np.take(values, i_var, axis=axis)
        else:
            raise ValueError(f"不支持的 values 维度: {values.shape}")

        # 把 plane 统一成 (nR, nZ)
        if dim is not None and dim[0] != dim[1] and plane.shape == dim[::-1]:
            plane = plane.T

        bounds = self._bounds_from_name(path) or self.bounds
        return plane, bounds

    def interpolator(self, step: int, var_name: str = None, index: int = None):
        """返回 (R,Z) -> value 的插值器（线性，区域外为 NaN），失败返回 None。"""
        real_step = self.nearest_step(step)
        if real_step is None:
            return None
        key = (real_step, var_name, index)
        if key in self._interp_cache:
            return self._interp_cache[key]

        path = self.step_map[real_step]
        plane, bounds = self.read_plane(path, var_name, index)
        if plane is None:
            return None
        if bounds is None:
            log(f"[h5] WARNING: 无法从文件名/--bounds 确定 {os.path.basename(path)} 的 "
                f"R/Z 范围，跳过该变量插值。")
            return None

        if self.debug:
            log(f"[h5][debug] step {real_step} {var_name}: plane.shape={plane.shape}, "
                f"bounds={tuple(round(b, 4) for b in bounds)}")

        rmin, rmax, zmin, zmax = bounds
        n1, n2 = plane.shape
        r_ax = np.linspace(rmin, rmax, n1)
        z_ax = np.linspace(zmin, zmax, n2)

        from scipy.interpolate import RegularGridInterpolator
        interp = RegularGridInterpolator(
            (r_ax, z_ax), plane, method="linear",
            bounds_error=False, fill_value=np.nan,
        )
        if len(self._interp_cache) >= self.max_interp_cache:
            self._interp_cache.pop(next(iter(self._interp_cache)))   # FIFO 淘汰
        self._interp_cache[key] = interp
        return interp


# --------------------------------------------------------------------------------------
# 场量组装
# --------------------------------------------------------------------------------------
def _time_derivative(Y, t):
    """对时间做数值导数 dY/dt（非等距中心差分；NaN/断段自动跳过）。

    Y: (n_shard, n_t)；t: (n_t,) 秒。返回同形状数组。
    只在 "t 与 Y 都有限" 的连续段内求导，孤立点/断段保持 NaN；
    段边缘退化为单侧差分（np.gradient 的默认行为）。
    """
    Y = np.asarray(Y, dtype=float)
    t = np.asarray(t, dtype=float)
    D = np.full(Y.shape, np.nan)
    for i in range(Y.shape[0]):
        ok = np.isfinite(Y[i]) & np.isfinite(t)
        if ok.sum() < 2:
            continue
        idx = np.nonzero(ok)[0]
        for seg in np.split(idx, np.nonzero(np.diff(idx) != 1)[0] + 1):
            if seg.size < 2:
                continue
            tt = t[seg]
            if np.any(np.diff(tt) <= 0):
                continue          # 时间不严格递增（数据异常），该段跳过
            D[i, seg] = np.gradient(Y[i, seg], tt)
    return D


def build_fields(series, time_map, time_scaling, t0, fields: RectangleFields | None,
                 color_mode: str, smooth: int, mask_dead: str,
                 color_var: str = None, color_index: int = None, debug: bool = False):
    """组装绘图所需数组。

    Returns
    -------
    dict with:
        steps, i_spi, t_ms, dV_R, acc, dr_mm, color, color_label, n_alive
    """
    steps = series["steps"]
    n_shard, n_step = series["VR"].shape
    radius = series["radius"]
    R, Z = series["R"], series["Z"]

    # --- 时间轴：t_si = time_map[step] * time_scaling，再相对 t0 取 ms ---
    t_si = np.array([time_map.get(int(s), np.nan) * time_scaling for s in steps])
    if np.all(np.isnan(t_si)):
        log("[time] WARNING: 所有 step 都不在 macroscopic_vars.dat 的 @times 表里，"
            "改用 step 号作为时间轴（单位: step）。")
        t_axis = steps.astype(float)
        t_label = "Time step (index)"
        t_ms = t_axis
    else:
        if t0 is None:
            if np.isnan(t_si).all():
                t0 = 0.0
            else:
                alive_any = (radius > 0).any(axis=0)
                first_alive = np.nonzero(alive_any)[0]
                j0 = int(first_alive[0]) if first_alive.size else 0
                t0 = float(t_si[j0])
                log(f"[time] WARNING: 未指定 t0，退化为第一个有存活碎片的时间步 "
                    f"t0 = {t0:.6e} s")
        t_ms = (t_si - t0) * 1e3
        t_label = "Time since ablation start (ms)"

    # --- 参考量：以第一个有效时间步（= 注入时刻）为基准 ---
    finite_cols = np.isfinite(series["VR"]).any(axis=0)
    ref_idx = int(np.argmax(finite_cols)) if finite_cols.any() else 0
    VR0 = series["VR"][:, ref_idx][:, None]
    r0 = radius[:, ref_idx][:, None]
    log(f"[ref] 参考时间步 = {int(steps[ref_idx])} (step), t = {t_ms[ref_idx]:.4f} ms (relative)")
    log(f"[ref] 参考时的碎片速度 V_R0 范围: "
        f"[{np.nanmin(VR0):.3f}, {np.nanmax(VR0):.3f}] m/s")

    dead = np.zeros_like(radius, dtype=bool)
    if mask_dead in ("dr", "both"):
        dead = radius <= 0.0

    dV_R = series["VR"] - VR0          # [m/s]，= ∫ a_R dt
    dr = r0 - radius                   # [m]，消融掉掉半径
    dr_mm = dr * 1e3                   # [mm]

    if mask_dead in ("dr", "both"):
        dr_mm = np.where(dead, np.nan, dr_mm)
    if mask_dead == "both":
        dV_R = np.where(dead, np.nan, dV_R)

    # --- 可选的时间平滑（f12.3 的 1 mm/s 量化会带来阶梯） ---
    if smooth and smooth > 1:
        dV_R = pd.DataFrame(dV_R.T).rolling(smooth, min_periods=1, center=True).mean().to_numpy().T
        dr_mm = pd.DataFrame(dr_mm.T).rolling(smooth, min_periods=1, center=True).mean().to_numpy().T

    # --- 火箭加速度 a_R = d(ΔV_R)/dt [m/s²] ---
    # 与 ΔV_R 图共用同一套（可选平滑后的）速度序列；t_ms 已相对 t0，差分不受零点影响。
    # 时间轴若退化为 step 序号（未匹配到 @times），无法定标，置 NaN 并提示。
    if t_label.startswith("Time since ablation start"):
        acc = _time_derivative(dV_R, t_ms / 1e3)
    else:
        acc = np.full_like(dV_R, np.nan)
        log("[acc] WARNING: 时间轴不是真实时间（未匹配到 @times），"
            "无法计算加速度 a_R=dV_R/dt，加速度图将跳过。")

    # --- 颜色场：每个碎片在每个时刻"看到的"局域量 ---
    color = np.full((n_shard, n_step), np.nan)
    color_label = ""
    if color_mode != "none":
        if fields is None or not fields.step_map:
            log("[color] WARNING: 未提供有效的 --rectangle_dir，颜色字段留空。"
                "（可加 --color none 或补上 rectangle 目录）")
        else:
            # specs: (键, 请求的变量名或 None, 备用 0-based 序号)
            primary = {"te": "T_e", "tn": "T_e", "psin": "Psi_N"}.get(color_mode, "color")
            custom_key = None
            if color_var or color_index is not None:
                # 显式指定变量名或序号 -> 单变量"自定义"模式，直接覆盖 --color 的默认变量
                custom_key = color_var or primary
                specs = [(custom_key, color_var, color_index)]
                eff_mode = "custom"
            elif color_mode == "te":
                specs, eff_mode = [("T_e", "T_e", None)], "te"
            elif color_mode == "tn":
                specs, eff_mode = [("T_e", "T_e", None), ("ne", "ne", None)], "tn"
            else:
                specs, eff_mode = [("Psi_N", "Psi_N", None)], "psin"

            probe = {k: fields.interpolator(int(steps[0]), vn, ix)
                     for k, vn, ix in specs}
            if all(v is None for v in probe.values()):
                log(f"[color] WARNING: rectangle 文件中取不到 "
                    f"{[vn for _, vn, _ in specs]}，颜色字段留空。"
                    "（检查 auto 脚本 Step 3 的 expressions 行；"
                    "或先用 --probe_h5 看变量清单）")
            else:
                for j, step in enumerate(steps):
                    pts = np.column_stack([R[:, j], Z[:, j]])
                    ok = np.isfinite(pts).all(axis=1)
                    vals = {}
                    for k, vn, ix in specs:
                        # 场随时间是演化的，必须逐 step 取插值器（内部有 FIFO 缓存）
                        itp = fields.interpolator(int(step), vn, ix)
                        if itp is None:
                            continue
                        cur = np.full(n_shard, np.nan)
                        cur[ok] = itp(pts[ok])
                        vals[k] = cur
                    if eff_mode == "te" and "T_e" in vals:
                        color[:, j] = vals["T_e"]
                    elif eff_mode == "tn" and "T_e" in vals and "ne" in vals:
                        color[:, j] = vals["T_e"] * vals["ne"]   # [eV * m^-3]
                    elif eff_mode == "psin" and "Psi_N" in vals:
                        color[:, j] = vals["Psi_N"]
                    elif eff_mode == "custom" and custom_key in vals:
                        color[:, j] = vals[custom_key]

                if eff_mode == "te":
                    color_label = "T_e seen by shard [eV]"
                elif eff_mode == "tn":
                    color_label = r"$T_e \cdot n_e$ seen by shard [eV$\cdot$m$^{-3}$]"
                elif eff_mode == "psin":
                    color_label = r"$\Psi_N$ at shard position"
                else:
                    color_label = f"{custom_key} at shard position (raw file units)"

    n_alive = (radius > 0).sum(axis=0)
    return {
        "steps": steps,
        "i_spi": series["i_spi"],
        "t_ms": t_ms,
        "t_label": t_label,
        "dV_R": dV_R,
        "acc": acc,
        "dr_mm": dr_mm,
        "color": color,
        "color_label": color_label,
        "n_alive": n_alive,
        "radius": radius,
        "r0": r0,
        "dead": dead,
    }


# --------------------------------------------------------------------------------------
# 绘图
# --------------------------------------------------------------------------------------
def _make_norm(color, mode: str, vmin=None, vmax=None, log_mode: str = "auto"):
    """构造（两张图共享的）颜色 norm —— 保证 Fig.1 / Fig.2 配色完全一致。"""
    finite = color[np.isfinite(color)]
    if finite.size == 0:
        return None, None
    if vmin is None:
        vmin = float(np.nanmin(finite))
    if vmax is None:
        vmax = float(np.nanmax(finite))
    use_log = (log_mode == "log") or (log_mode == "auto" and mode in ("te", "tn"))
    if use_log and vmin > 0 and vmax > vmin:
        return mcolors.LogNorm(vmin=vmin, vmax=vmax), "log"
    return mcolors.Normalize(vmin=vmin, vmax=vmax), "lin"


def _mask_for_norm(C, norm):
    """对数/线性 norm 之前的清洗：对数标度下必须剔除 <=0 与 NaN。"""
    C = np.asarray(C, dtype=float)
    if isinstance(norm, mcolors.LogNorm):
        return np.where(np.isfinite(C) & (C > 0), C, np.nan)
    return np.where(np.isfinite(C), C, np.nan)


def _alpha_from_radius(radius, r0, a_min: float, a_max: float):
    """透明度映射：用 r/r0 线性映射到 [a_min, a_max]。

    刚注入 (r=r0) -> 不透明；完全烧蚀 (r=0) -> 最透明（但保留 a_min 以便仍可见）。
    这样"碎片还剩多少"这个信息就不必再占一张图。
    """
    r0 = np.asarray(r0, dtype=float)
    r0 = np.where(np.isfinite(r0) & (r0 > 0), r0, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        rel = np.clip(np.asarray(radius, dtype=float) / r0, 0.0, 1.0)
    rel = np.nan_to_num(rel, nan=0.0, posinf=1.0, neginf=0.0)
    return a_min + (a_max - a_min) * rel


def _lines_from_grid(xs, ts, Z, C, A):
    """把 (n_shard, n_t) 网格转成逐碎片的 3D 折线段 + 每段颜色/透明度。

    用折线渲染时，每个碎片就是**一条线**（x 固定、沿时间走），
    不会被 plot_surface 的网格四边形连成整片曲面。
    """
    segs, cols, alphas = [], [], []
    nt = len(ts)
    for i in range(len(xs)):
        xv = float(xs[i])
        z, c, a = Z[i], C[i], A[i]
        finite = np.isfinite(z)
        for j in range(nt - 1):
            if not (finite[j] and finite[j + 1]):
                continue
            segs.append([(xv, ts[j], z[j]), (xv, ts[j + 1], z[j + 1])])
            if np.isfinite(c[j]) or np.isfinite(c[j + 1]):
                cols.append(float(np.nanmean([c[j], c[j + 1]])))
            else:
                cols.append(np.nan)
            alphas.append(float(np.nanmean([a[j], a[j + 1]])))
    return segs, np.asarray(cols, dtype=float), np.asarray(alphas, dtype=float)


def _set_limits_3d(ax, xs, ts, Z):
    """显式设定 3D 轴范围。

    只用 Line3DCollection 时 mplot3d **不会**自动缩放（会画出一张空图），
    所以无论曲面还是折线，都在最后统一设一遍。
    """
    if np.isfinite(xs).any():
        ax.set_xlim(float(np.nanmin(xs)), float(np.nanmax(xs)))
    if np.isfinite(ts).any():
        ax.set_ylim(float(np.nanmin(ts)), float(np.nanmax(ts)))
    zf = np.isfinite(Z)
    if zf.any():
        z0, z1 = float(np.nanmin(Z)), float(np.nanmax(Z))
        pad = 0.05 * (z1 - z0) if z1 > z0 else abs(z0) * 0.05 + 1e-9
        ax.set_zlim(z0 - pad, z1 + pad)


def plot_merged_fig(xs, ts, Z, C, A, norm, cmap, title, z_label, x_label, t_label,
                    color_label, alpha_note, out_png, elev=28, azim=-120,
                    render="surface", line_width=1.0):
    """合并图：z 轴 = 速度衰减量，颜色 = 碎片见到的局域量，透明度 = 半径。

    render='surface' 画网格曲面（相邻碎片会被连成整片，视觉上像一排帘子）；
    render='lines' 每个碎片画一条 3D 折线，语义上就是"每个弹丸一条线"。
    """
    if not np.isfinite(Z).any():
        log(f"[out] WARNING: {os.path.basename(out_png)} 无有效数据，跳过。")
        return
    Xg, Yg = np.meshgrid(xs, ts, indexing="ij")
    fig = plt.figure(figsize=(11.5, 7.6))
    ax = fig.add_subplot(111, projection="3d")

    has_color = (norm is not None) and bool(np.isfinite(np.asarray(C, dtype=float)).any())

    if render in ("surface", "both"):
        surf_kw = dict(rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)
        if has_color:
            rgba = cmap(norm(_mask_for_norm(C, norm)))
            rgba[..., 3] = A                      # 透明度 = r/r0 映射
            surf_kw["facecolors"] = rgba
        ax.plot_surface(Xg, Yg, Z, **surf_kw)

    if render in ("lines", "both"):
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        segs, cols, alphas = _lines_from_grid(xs, ts, Z, C, A)
        if segs:
            lc = Line3DCollection(segs, linewidths=line_width)
            if has_color:
                base = cmap(norm(np.nan_to_num(cols, nan=norm.vmin)))
                base[:, 3] = np.nan_to_num(alphas, nan=0.15)
                lc.set_color(base)
            else:
                lc.set_color([(0.15, 0.35, 0.9, 0.55)])
            ax.add_collection3d(lc)
            log(f"[out] 折线渲染: {len(xs)} 条碎片轨迹 / {len(segs)} 段")

    ax.set_xlabel(x_label, labelpad=10)
    ax.set_ylabel(t_label, labelpad=10)
    ax.set_zlabel(z_label, labelpad=8)
    ax.set_title(title)
    ax.view_init(elev=elev, azim=azim)
    _set_limits_3d(ax, xs, ts, Z)

    if has_color:
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, shrink=0.62, pad=0.10, aspect=22)
        cb.set_label(color_label)

    if alpha_note:
        fig.text(0.015, 0.015, alpha_note, fontsize=9, va="bottom")

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    log(f"[out] {out_png}")


def plot_merged_2d(xs, ts, Z, C, A, norm, cmap, title, z_label, color_label,
                   t_label, out_png):
    """合并图的 2D 伴侣三联图：被画量 / 颜色场 / 透明度(=半径)，便于读数值。"""
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.0), constrained_layout=True)
    Xg, Yg = np.meshgrid(xs, ts, indexing="ij")

    p0 = axes[0].pcolormesh(Xg, Yg, Z, shading="nearest", cmap="viridis")
    fig.colorbar(p0, ax=axes[0], label=z_label)
    axes[0].set_title(title)

    if norm is not None and np.isfinite(np.asarray(C, dtype=float)).any():
        p1 = axes[1].pcolormesh(Xg, Yg, _mask_for_norm(C, norm), shading="nearest",
                                cmap=cmap, norm=norm)
        fig.colorbar(p1, ax=axes[1], label=color_label)
        axes[1].set_title("Color field")
    else:
        axes[1].set_title("Color field unavailable")

    p2 = axes[2].pcolormesh(Xg, Yg, A, shading="nearest", cmap="Greys", vmin=0.0, vmax=1.0)
    fig.colorbar(p2, ax=axes[2], label="opacity = r/r0")
    axes[2].set_title("Opacity map (shard radius)")

    for a in axes:
        a.set_xlabel("Shard index")
        a.set_ylabel(t_label)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    log(f"[out] {out_png}")


def _stride_indices(n: int, stride: int):
    """抽稀索引：等间隔取样，并保证末点被保留。

    多张图（z / 颜色 / 透明度 / 2D 伴侣图）共用同一组索引，避免网格错位。
    """
    stride = max(1, int(stride))
    idx = list(range(0, n, stride))
    if idx[-1] != n - 1:
        idx.append(n - 1)
    return np.asarray(idx)


def _downsample(x_vals, t_vals, Z, C, max_x=250, max_t=200,
                x_stride=None, t_stride=None):
    nx = len(x_vals)
    nt = len(t_vals)
    if x_stride is None:
        x_stride = max(1, int(np.ceil(nx / max_x)))
    if t_stride is None:
        t_stride = max(1, int(np.ceil(nt / max_t)))
    # 保证末尾点被保留
    xi = list(range(0, nx, x_stride))
    if xi[-1] != nx - 1:
        xi.append(nx - 1)
    ti = list(range(0, nt, t_stride))
    if ti[-1] != nt - 1:
        ti.append(nt - 1)
    xi = np.asarray(xi)
    ti = np.asarray(ti)
    return x_vals[xi], t_vals[ti], Z[np.ix_(xi, ti)], C[np.ix_(xi, ti)]


def plot_surface_fig(x_vals, t_vals, Z, C, norm, cmap, title, z_label, x_label,
                     t_label, color_label, out_png, elev=28, azim=-120):
    if not np.isfinite(Z).any():
        log(f"[out] WARNING: {os.path.basename(out_png)} 无有效数据，跳过。")
        return
    Xg, Yg = np.meshgrid(x_vals, t_vals, indexing="ij")
    fig = plt.figure(figsize=(11, 7.5))
    ax = fig.add_subplot(111, projection="3d")

    # 注意：显式传 facecolors=None 会让 matplotlib 抛
    # UnboundLocalError: cannot access local variable 'color'
    # 所以只有在真的有颜色数组时才把 facecolors 放进 kwargs。
    surf_kw = dict(rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)
    if norm is not None:
        # 对数标度下必须把 <=0 / NaN 置为 NaN，否则会被画成 bad 色
        facecolors = cmap(norm(_mask_for_norm(C, norm)))
        if np.isfinite(np.asarray(C, dtype=float)).any():
            surf_kw["facecolors"] = facecolors
        else:
            log("[out] NOTE: 颜色场全为 NaN，3D 曲面将用单色绘制。")
    ax.plot_surface(Xg, Yg, Z, **surf_kw)
    ax.set_xlabel(x_label, labelpad=10)
    ax.set_ylabel(t_label, labelpad=10)
    ax.set_zlabel(z_label, labelpad=8)
    ax.set_title(title)
    ax.view_init(elev=elev, azim=azim)
    _set_limits_3d(ax, x_vals, t_vals, Z)

    if norm is not None:
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, shrink=0.62, pad=0.10, aspect=22)
        cb.set_label(color_label)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    log(f"[out] {out_png}")


def plot_2d_companion(x_vals, t_vals, Z, C, norm, cmap, title, z_label,
                      color_label, t_label, out_png):
    """(时间, 碎片) 平面上的两张热力图：左=被画量, 右=颜色场（与 3D 图同色标）。"""
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.2), constrained_layout=True)
    Xg, Yg = np.meshgrid(x_vals, t_vals, indexing="ij")

    p0 = axes[0].pcolormesh(Xg, Yg, Z, shading="nearest", cmap="RdBu_r")
    fig.colorbar(p0, ax=axes[0], label=z_label)
    axes[0].set_title(title)
    axes[0].set_xlabel("Shard index")
    axes[0].set_ylabel(t_label)

    if norm is None or not np.isfinite(np.asarray(C, dtype=float)).any():
        axes[1].set_title("Color field unavailable (no rectangle data / variable)")
        axes[1].set_xlabel("Shard index")
        axes[1].set_ylabel(t_label)
    else:
        p1 = axes[1].pcolormesh(Xg, Yg, _mask_for_norm(C, norm), shading="nearest",
                               cmap=cmap, norm=norm)
        fig.colorbar(p1, ax=axes[1], label=color_label)
        axes[1].set_title("Color field (same scale as 3D figure)")
        axes[1].set_xlabel("Shard index")
        axes[1].set_ylabel(t_label)

    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    log(f"[out] {out_png}")


def plot_html(x_vals, t_vals, Z, C, norm, title, z_label, color_label,
              t_label, out_html):
    """可交互 HTML（需要 plotly；未安装则跳过）。"""
    try:
        import plotly.graph_objects as go
    except Exception:  # noqa: BLE001
        log("[out] plotly 未安装，跳过交互 HTML（pip install plotly 后可启用）。")
        return
    if not np.isfinite(Z).any():
        log(f"[out] WARNING: {os.path.basename(out_html)} 无有效数据，跳过。")
        return

    Xg, Yg = np.meshgrid(x_vals, t_vals, indexing="ij")
    if norm is not None:
        # norm() 对 log 标度可能返回 masked array，统一填成 NaN 交给 plotly
        norm_C = np.ma.filled(norm(_mask_for_norm(C, norm)), np.nan)
    else:
        norm_C = np.zeros_like(C)
    fig = go.Figure(data=[go.Surface(
        x=Xg, y=Yg, z=Z,
        surfacecolor=norm_C,
        cmin=0.0, cmax=1.0,
        colorscale="Turbo",
        colorbar=dict(title=color_label),
    )])
    fig.update_layout(
        title=title,
        scene=dict(xaxis_title="Shard index i_spi", yaxis_title=t_label, zaxis_title=z_label),
        margin=dict(l=0, r=0, t=50, b=0),
    )
    fig.write_html(out_html, include_plotlyjs="cdn")
    log(f"[out] {out_html}")


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        description="SPI 碎片火箭效应三维曲面诊断（ΔV_R / a_R / Δr vs 碎片编号 vs 时间，"
                    "颜色 = 碎片见到的 T_e）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--csv_dir", type=str, default=None,
                   help="含 fragments_position_*.csv 与 fragments_plume_total_*.csv 的目录"
                        "（除 --probe_h5 外必需）")
    p.add_argument("--rectangle_dir", type=str, default=None,
                   help="含 exprs_Rmin*.h5（postproc rectangle）的目录，用于颜色场")
    p.add_argument("--abl_file", type=str, default="abl_history.dat",
                   help="jorek2_diagno_spi 产出的 abl_history.dat，用于自动定位消融起点")
    p.add_argument("--live_data", type=str, default="macroscopic_vars.dat",
                   help="macroscopic_vars.dat，用于时间归一化")
    p.add_argument("--start_time", type=float, default=None,
                   help="显式指定 t0（SI 秒）；缺省则从 abl_history.dat 自动识别")
    # 与 track_largest_pellets.py 保持一致的命名：--min_time_ms / --max_time_ms
    p.add_argument("--t_min_ms", "--min_time_ms", dest="t_min_ms", type=float, default=0.0,
                   help="只绘制 t-t0 >= 该值（ms）的部分")
    p.add_argument("--t_max_ms", "--max_time_ms", dest="t_max_ms", type=float, default=0.0,
                   help="只绘制 t-t0 <= 该值（ms）的部分；0 = 全部")

    p.add_argument("--color", type=str, default="te", choices=COLOR_CHOICES,
                   help="颜色场：te=T_e[eV], tn=T_e*ne, psin=Psi_N, none=不画颜色")
    p.add_argument("--color_var", type=str, default=None,
                   help="直接指定 rectangle 文件里的变量名作为颜色场（覆盖 --color 的默认变量）")
    p.add_argument("--color_index", type=int, default=None,
                   help="按 0-based 序号指定颜色变量（变量名解析失败时的兜底）")
    p.add_argument("--expr_order", type=str, default=None,
                   help="生成 rectangle h5 时 postproc 'expressions' 的顺序，"
                        "例如 'Psi_N T_e nimp ne'。"
                        "postproc 把 h5 里的变量名数组写成了 character(len=1)，"
                        "变量名会退化成首字母，此时用本参数做位置映射")
    p.add_argument("--probe_h5", action="store_true",
                   help="只检查 rectangle HDF5 的变量清单与各变量取值范围后退出（排错用）")
    p.add_argument("--color_log", type=str, default="auto", choices=["auto", "log", "lin"],
                   help="颜色的对数/线性标度（te/tn 默认 log）")
    p.add_argument("--vmin", type=float, default=None, help="颜色下限（两张图共用）")
    p.add_argument("--vmax", type=float, default=None, help="颜色上限（两张图共用）")

    p.add_argument("--x_axis", type=str, default="index",
                   choices=["index", "r0", "psi0", "death"],
                   help="x 轴排序：index=碎片编号 i_spi（默认）, r0=按初始半径, psi0=按初始 Psi_N, "
                        "death=按半径置零(烧蚀完成)的先后（时间窗内越早置零越靠前；"
                        "窗口末仍未置零的碎片按末时刻半径升序排在其后）")
    p.add_argument("--mask_dead", type=str, default="dr", choices=["dr", "both", "none"],
                   help="把已烧蚀完(r=0)的碎片置 NaN：dr=只作用在 Δr 图, both=两图都置, none=都不置")
    p.add_argument("--dvr_sign", type=str, default="abs", choices=["abs", "minus", "plus"],
                   help="Fig.1 的 z 量：abs=|ΔV_R|=|V_R(t)-V_R(注入)|（默认，即速度衰减幅度）、"
                        "minus=-ΔV_R、plus=+ΔV_R")
    p.add_argument("--accel", type=str, default="on", choices=["on", "off"],
                   help="是否额外出碎片加速度合并图（z=|a_R|=|dV_R/dt|，m/s²）")
    p.add_argument("--acc_sign", type=str, default="abs", choices=["abs", "minus", "plus"],
                   help="加速度图的 z 量：abs=|a_R|（默认，加速度幅度）、minus=-a_R、plus=+a_R")
    p.add_argument("--layout", type=str, default="merged", choices=["merged", "split", "both"],
                   help="merged=只出合并图（z=|ΔV_R|，颜色=T_e，透明度=半径）；"
                        "split=分别出 ΔV_R 与 Δr 两张图；both=都出")
    p.add_argument("--render", type=str, default="surface", choices=["surface", "lines", "both"],
                   help="渲染方式：surface=网格曲面（相邻碎片会被连成整片）；"
                        "lines=每个碎片一条 3D 折线；both=两者叠加")
    p.add_argument("--alpha", type=str, default="radius", choices=["radius", "none"],
                   help="曲面/折线的透明度映射：radius=按 r/r0（新鲜=不透明，烧蚀完=透明）")
    p.add_argument("--alpha_min", type=float, default=0.06, help="透明度下限（已烧蚀完的碎片）")
    p.add_argument("--alpha_max", type=float, default=1.0, help="透明度上限（刚注入的碎片）")
    p.add_argument("--smooth", type=int, default=0,
                   help="沿时间做移动平均的窗口（0/1=不平滑）；速度字段是 f12.3 量化过的，"
                        "时间微分出的加速度图对阶梯噪声尤其敏感，可适当加大该窗口")
    p.add_argument("--max_plotted_shards", type=int, default=0,
                   help="最多绘制多少个碎片（0=全部，>0 时均匀抽稀）")
    p.add_argument("--time_stride", type=int, default=1, help="时间步抽稀间隔")
    p.add_argument("--max_x_points", type=int, default=250, help="x 方向最多网格点数")
    p.add_argument("--max_t_points", type=int, default=200, help="y 方向最多网格点数")
    p.add_argument("--max_cells", type=int, default=60000,
                   help="3D 曲面最大格子数（超出会自动加大抽稀步长，避免渲染过慢）")
    p.add_argument("--elev", type=float, default=28.0, help="3D 视角仰角")
    p.add_argument("--azim", type=float, default=-120.0, help="3D 视角方位角")
    p.add_argument("--bounds", type=float, nargs=4, metavar=("RMIN", "RMAX", "ZMIN", "ZMAX"),
                   default=None, help="rectangle 的 R/Z 范围（文件名解析失败时的兜底）")

    p.add_argument("--out_dir", type=str, default=".", help="输出目录")
    p.add_argument("--out_prefix", type=str, default="shard_surfaces", help="输出文件名前缀")
    p.add_argument("--png", dest="png", action="store_true", default=True, help="输出 PNG")
    p.add_argument("--no-png", dest="png", action="store_false", help="不输出 PNG")
    p.add_argument("--html", type=str, default="auto", choices=["auto", "on", "off"],
                   help="交互 HTML：auto=有 plotly 才出")
    p.add_argument("--also2d", dest="also2d", action="store_true", default=True,
                   help="额外输出 2D 伴侣热力图")
    p.add_argument("--no-also2d", dest="also2d", action="store_false", help="不输出 2D 伴侣图")
    p.add_argument("--cache_npz", type=str, default=None,
                   help="把组装好的数组存/读 npz（重画时秒级，避免重复读 h5）")
    p.add_argument("--debug", type=int, default=0, help="1=打印调试信息")
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    debug = int(args.debug) == 1
    os.makedirs(args.out_dir, exist_ok=True)

    # ---------- 0. 排错模式：只看 rectangle HDF5 的变量清单 ----------
    if args.probe_h5:
        if not args.rectangle_dir:
            raise SystemExit("--probe_h5 需要同时给出 --rectangle_dir")
        fields = RectangleFields(args.rectangle_dir, bounds=args.bounds, debug=True,
                                 expr_order=args.expr_order)
        if not fields.step_map:
            raise SystemExit(f"{args.rectangle_dir} 中没有 exprs_Rmin*.h5")
        step = min(fields.step_map)
        path = fields.step_map[step]
        log(f"[probe] 样本文件: {path}")
        log(f"[probe] {fields.describe(path)}")
        info = fields._var_index(path)
        for i in range(max(1, info["n_var"] or len(info["names"]))):
            try:
                plane, bounds = fields.read_plane(path, index=i)
            except Exception as exc:  # noqa: BLE001
                log(f"[probe]   [{i}] 读取失败: {exc}")
                continue
            if plane is None:
                continue
            nm = info["names"][i] if i < len(info["names"]) else f"<index {i}>"
            finite = plane[np.isfinite(plane)]
            if finite.size:
                log(f"[probe]   [{i}] {nm:<14s} shape={plane.shape} "
                    f"min={finite.min():.6g}  max={finite.max():.6g}")
            else:
                log(f"[probe]   [{i}] {nm:<14s} shape={plane.shape} 全为 NaN")
        log("[probe] 提示: 用 --color_var <名字> 或 --color_index <序号> 指定颜色变量；"
            "若上面的 variables 只剩首字母，请用 --expr_order \"Psi_N T_e nimp ne\" "
            "（与生成 h5 时 postproc 的 expressions 行一致）恢复顺序。")
        return 0

    if not args.csv_dir:
        parser.error("--csv_dir 是必需的（除非只做 --probe_h5）")

    # ---------- 1. 时间归一化 + 消融起点（与 track_largest_pellets.py 一致） ----------
    if not os.path.exists(args.live_data):
        raise FileNotFoundError(f"FATAL: 找不到 {args.live_data}（时间归一化必需）")
    time_map, time_scaling = load_time_mapping(args.live_data)
    log(f"[time] 载入 {len(time_map)} 个时间步，time_scaling = {time_scaling:.6e} s/step")
    t0 = load_ablation_start(args.abl_file, time_scaling, args.start_time)

    # ---------- 2. 缓存 / 读取碎片时间序列 ----------
    if args.cache_npz and os.path.exists(args.cache_npz):
        log(f"[cache] 从 {args.cache_npz} 读取已组装的数组")
        with np.load(args.cache_npz, allow_pickle=True) as data:
            fields_arrays = {k: data[k] for k in data.files}
        t_label = str(fields_arrays.get("t_label",
                                        np.array("Time since ablation start (ms)")))
        color_label = str(fields_arrays.get("color_label", np.array("")))
    else:
        steps, _, _ = collect_steps(args.csv_dir)
        if not steps:
            raise RuntimeError(f"{args.csv_dir} 中没有成对的 fragments CSV 文件。")
        log(f"[csv] 找到 {len(steps)} 个时间步：{steps[0]} ... {steps[-1]}")
        series = load_shard_time_series(args.csv_dir, steps, debug=debug)
        log(f"[csv] 碎片数 = {series['VR'].shape[0]}，时间步数 = {series['VR'].shape[1]}")

        fields = None
        if args.color != "none" and args.rectangle_dir:
            fields = RectangleFields(args.rectangle_dir, bounds=args.bounds, debug=debug,
                                     expr_order=args.expr_order)
            if not fields.step_map:
                log(f"[h5] WARNING: {args.rectangle_dir} 中没有 exprs_Rmin*.h5")
        arrays = build_fields(series, time_map, time_scaling, t0, fields,
                              args.color, args.smooth, args.mask_dead,
                              color_var=args.color_var, color_index=args.color_index,
                              debug=debug)
        fields_arrays = arrays
        t_label = arrays["t_label"]
        color_label = arrays["color_label"]

        if args.cache_npz:
            # fields_arrays 已含 t_label / color_label，直接整体落盘
            np.savez_compressed(args.cache_npz, **fields_arrays)
            log(f"[cache] 数组已保存到 {args.cache_npz}")

    steps = fields_arrays["steps"]
    i_spi = fields_arrays["i_spi"]
    t_ms = fields_arrays["t_ms"]
    dV_R = fields_arrays["dV_R"]
    dr_mm = fields_arrays["dr_mm"]
    color = fields_arrays["color"]
    n_alive = fields_arrays["n_alive"]
    if "acc" in fields_arrays:
        acc_full = np.asarray(fields_arrays["acc"], dtype=float)
    elif str(t_label).startswith("Time since ablation start"):
        # 兼容旧版缓存：缺少加速度数组时用 ΔV_R 与时间轴现算（差分只关心时间差）
        acc_full = _time_derivative(np.asarray(dV_R, dtype=float),
                                    np.asarray(t_ms, dtype=float) / 1e3)
        log("[cache] NOTE: 缓存中没有加速度数组，已用 ΔV_R 与时间轴重算 a_R。")
    else:
        acc_full = np.full_like(np.asarray(dV_R, dtype=float), np.nan)
        log("[acc] WARNING: 缓存缺少加速度数组且时间轴不是真实时间，跳过加速度图。")

    # ---------- 3. 时间窗裁剪 ----------
    tsel = np.ones_like(t_ms, dtype=bool)
    tsel &= t_ms >= args.t_min_ms
    if args.t_max_ms and args.t_max_ms > 0:
        tsel &= t_ms <= args.t_max_ms
    if tsel.sum() < 2:
        log(f"[time] WARNING: 时间窗筛选后只剩 {int(tsel.sum())} 步，忽略时间窗设置。")
        tsel = np.ones_like(t_ms, dtype=bool)
        tsel &= np.isfinite(t_ms)
    t_sel = t_ms[tsel]
    dV_sel = dV_R[:, tsel]
    dr_sel = dr_mm[:, tsel]
    C_sel = color[:, tsel]
    acc_sel = acc_full[:, tsel]
    # 半径（用于透明度映射）：r/r0，r0 = 注入时刻每个碎片自己的半径
    if "radius" in fields_arrays:
        r_sel = np.asarray(fields_arrays["radius"], dtype=float)[:, tsel]
    else:
        r_sel = np.full_like(C_sel, np.nan)
    r0_src = fields_arrays.get("r0") if hasattr(fields_arrays, "get") else None
    if r0_src is None:
        r0_src = np.asarray(fields_arrays["radius"], dtype=float)[:, [0]]
    r0_vec = np.asarray(r0_src, dtype=float).reshape(-1)
    log(f"[time] 绘图时间范围: {t_sel.min():.4f} ~ {t_sel.max():.4f} ms (相对消融起点)，"
        f"{int(tsel.sum())} 个时间步")

    # ---------- 4. x 轴排序 ----------
    if args.x_axis == "r0":
        key = fields_arrays["radius"][:, 0] if "radius" in fields_arrays else None
        order = np.argsort(-np.nan_to_num(key, nan=-np.inf))
        x_vals = np.nan_to_num(key[order], nan=0.0)
        x_label = r"Initial shard radius $r_0$ [m]"
    elif args.x_axis == "psi0":
        key = C_sel[:, 0] if args.color == "psin" else None
        if key is None or not np.isfinite(key).any():
            log("[x] WARNING: psi0 排序需要 --color psin 且 rectangle 数据可用，"
                "退化为 index 排序。")
            order = np.argsort(i_spi)
            x_vals = i_spi[order].astype(float)
            x_label = "Shard index i_spi"
        else:
            order = np.argsort(np.nan_to_num(key, nan=np.inf))
            x_vals = np.nan_to_num(key[order], nan=0.0)
            x_label = r"Initial $\Psi_N$"
    elif args.x_axis == "death":
        # 按"半径置零"（烧蚀完成）的先后排序：
        #   1) 时间窗内半径首次置零的碎片排在前（越早置零越靠前）；
        #   2) 窗口结束仍未置零的碎片排在其后，按窗口最后时刻的半径升序
        #      （越接近烧蚀完越靠前）；
        #   3) 键完全相同时按原始碎片编号 i_spi 保持稳定。
        # 注：排序在 --t_min_ms/--t_max_ms 裁剪后的时间窗内进行，窗口外才置零的
        # 碎片在本次绘图中按"仍存活"处理。
        rr = np.asarray(r_sel, dtype=float)
        fin = np.isfinite(rr)
        if not (fin & (rr > 0)).any():
            log("[x] WARNING: 时间窗内没有可用的半径数据，death 排序退化为 index 排序。")
            order = np.argsort(i_spi)
            x_vals = i_spi[order].astype(float)
            x_label = "Shard index i_spi"
        else:
            n_sh, nt = rr.shape
            idx_all = np.arange(nt)
            alive = fin & (rr > 0)
            # "出生" = 首次出现正半径的索引（从未存活记 nt，避免把注入前的 0 当成置零）
            birth = np.where(alive.any(axis=1), alive.argmax(axis=1), nt)
            # 出生之后首次出现 r<=0 的索引 = 该碎片置零（烧蚀完成）的时刻
            died_mat = fin & (rr <= 0) & (idx_all[None, :] > birth[:, None])
            died = died_mat.any(axis=1)
            death_idx = np.where(died, died_mat.argmax(axis=1), nt)
            # 末时刻半径（取窗口内最后一个有限值），仅作为仍存活碎片的二级排序键
            fin_rev = fin[:, ::-1]
            last_idx = nt - 1 - fin_rev.argmax(axis=1)
            final_r = rr[np.arange(n_sh), last_idx]
            final_r = np.where(np.isfinite(final_r) & ~died, final_r, 0.0)
            order = np.lexsort((i_spi, final_r, death_idx))
            x_vals = np.arange(1, n_sh + 1, dtype=float)
            x_label = "Shard order (1 = first fully ablated)"
            n_died = int(died.sum())
            if n_died:
                t_die = np.array([t_sel[int(k)] for k in death_idx[died]], dtype=float)
                log(f"[x] death 排序: {n_died}/{n_sh} 个碎片在时间窗内烧蚀完"
                    f"（置零于 {t_die.min():.4f} ~ {t_die.max():.4f} ms）；"
                    f"其余 {n_sh - n_died} 个按末时刻半径升序排在其后。")
            else:
                log(f"[x] death 排序: 时间窗内没有碎片烧蚀完，全部按末时刻半径升序排列。")
            if debug:
                for rank0, k in enumerate(order[:5], start=1):
                    when = f"{t_sel[int(death_idx[k])]:.4f} ms" if died[k] else "alive"
                    log(f"[x][debug] rank {rank0}: i_spi={i_spi[k]}, death={when}, "
                        f"r_end={final_r[k]:.4e} m")
    else:
        order = np.argsort(i_spi)
        x_vals = i_spi[order].astype(float)
        x_label = "Shard index i_spi"

    dV_sel = dV_sel[order, :]
    dr_sel = dr_sel[order, :]
    C_sel = C_sel[order, :]
    r_sel = r_sel[order, :]
    acc_sel = acc_sel[order, :]
    r0_vec = r0_vec[order]
    n_alive = n_alive[tsel]

    # ---------- 5. 抽稀 ----------
    if args.max_plotted_shards and 0 < args.max_plotted_shards < len(x_vals):
        xi = np.linspace(0, len(x_vals) - 1, args.max_plotted_shards).astype(int)
        xi = np.unique(xi)
        x_vals, dV_sel, dr_sel, C_sel, r_sel, acc_sel = (
            x_vals[xi], dV_sel[xi, :], dr_sel[xi, :], C_sel[xi, :], r_sel[xi, :], acc_sel[xi, :])
        r0_vec = r0_vec[xi]
        log(f"[plot] x 方向抽稀到 {len(x_vals)} 个碎片")

    # z 量：abs = |ΔV_R| = |V_R(t) - V_R(注入)|（默认，“速度衰减幅度”，朝上为正）
    if args.dvr_sign == "abs":
        dV_plot = np.abs(dV_sel)
        dvr_label, dvr_txt = r"$|\Delta V_R|$  [m/s]", "|ΔV_R|"
    elif args.dvr_sign == "minus":
        dV_plot, dvr_label, dvr_txt = -dV_sel, r"$-\Delta V_R$  [m/s]", "-ΔV_R"
    else:
        dV_plot, dvr_label, dvr_txt = dV_sel, r"$\Delta V_R$  [m/s]", "ΔV_R"
    log(f"[plot] Fig.1 z 轴 = {dvr_txt} [m/s]，"
        f"范围 [{np.nanmin(dV_plot):.3f}, {np.nanmax(dV_plot):.3f}]")

    # 加速度 z 量：a_R = d(ΔV_R)/dt（默认取幅度 |a_R|，与 |ΔV_R| 的约定一致）
    if args.acc_sign == "abs":
        acc_plot = np.abs(acc_sel)
        acc_label, acc_txt = r"$|a_R|$  [m/s$^2$]", "|a_R|"
    elif args.acc_sign == "minus":
        acc_plot, acc_label, acc_txt = -acc_sel, r"$-a_R$  [m/s$^2$]", "-a_R"
    else:
        acc_plot, acc_label, acc_txt = acc_sel, r"$a_R$  [m/s$^2$]", "a_R"
    if args.accel == "on":
        if np.isfinite(acc_plot).any():
            log(f"[plot] 加速度图 z 轴 = {acc_txt} [m/s²]，"
                f"范围 [{np.nanmin(acc_plot):.3e}, {np.nanmax(acc_plot):.3e}]")
        else:
            log("[plot] 加速度图无有效数据（时间轴退化或速度缺失），将跳过输出。")

    # x/y 抽稀步长只算一次，保证两张图的网格完全一致。
    # 注意：--time_stride 默认是 1，此时应当由 --max_t_points 自动决定步长，
    # 否则 600+ 个时间步会直接铺满整个曲面。
    nx, nt = len(x_vals), len(t_sel)
    x_stride_eff = max(1, int(np.ceil(nx / max(1, args.max_x_points))))
    if args.time_stride > 1:
        t_stride_eff = args.time_stride
    else:
        t_stride_eff = max(1, int(np.ceil(nt / max(1, args.max_t_points))))
    n_cells = int(np.ceil(nx / x_stride_eff)) * int(np.ceil(nt / t_stride_eff))
    if n_cells > args.max_cells:
        factor = np.sqrt(n_cells / float(args.max_cells)) * 1.05
        x_stride_eff = max(x_stride_eff, int(np.ceil(x_stride_eff * factor)))
        t_stride_eff = max(t_stride_eff, int(np.ceil(t_stride_eff * factor)))
        log(f"[plot] 网格 {n_cells} 超过 --max_cells={args.max_cells}，"
            f"自动加大抽稀步长")
    log(f"[plot] 抽稀步长: x_stride={x_stride_eff}, t_stride={t_stride_eff} "
        f"(原始 {nx} 碎片 x {nt} 时间步)")

    # 透明度：把"碎片还剩多少"编码成 alpha（r/r0），这样就不必再占一张 Δr 的图
    if args.alpha == "radius":
        alpha_full = _alpha_from_radius(r_sel, r0_vec[:, None], args.alpha_min, args.alpha_max)
        alpha_note = (f"Opacity = r/r0  (opaque = fresh shard r=r0, transparent = fully ablated; "
                      f"alpha in [{args.alpha_min:.2f}, {args.alpha_max:.2f}])")
        log(f"[plot] 透明度映射 = r/r0 (alpha {args.alpha_min:.2f}..{args.alpha_max:.2f}); "
            f"r0 范围 [{np.nanmin(r0_vec) * 1e3:.4f}, {np.nanmax(r0_vec) * 1e3:.4f}] mm")
    else:
        alpha_full = np.ones_like(r_sel)
        alpha_note = ""

    # 抽稀：x/y 用同一组索引，保证 z / 颜色 / 透明度三个数组严格对齐
    xi = _stride_indices(len(x_vals), x_stride_eff)
    ti = _stride_indices(len(t_sel), t_stride_eff)
    xs = x_vals[xi]
    ts = t_sel[ti]
    Z1 = dV_plot[np.ix_(xi, ti)]        # 速度衰减量
    Z2 = dr_sel[np.ix_(xi, ti)]         # 半径损失
    Z3 = acc_plot[np.ix_(xi, ti)]       # 火箭加速度
    C1 = C_sel[np.ix_(xi, ti)]          # 颜色场
    A1 = alpha_full[np.ix_(xi, ti)]     # 透明度
    log(f"[plot] 网格: {Z1.shape[0]} x {Z1.shape[1]} (碎片 x 时间)")

    # ---------- 6. 共享颜色标度（两张图/两个渲染共用同一 norm） ----------
    norm, scale_used = _make_norm(C1, args.color, args.vmin, args.vmax, args.color_log)
    cmap = plt.get_cmap("turbo")
    if color_label and norm is not None:
        log(f"[color] {color_label}，范围 [{norm.vmin:.4g}, {norm.vmax:.4g}]，标度 = {scale_used}")

    # ---------- 7. 出图 ----------
    prefix = os.path.join(args.out_dir, args.out_prefix)
    t_min_lbl = t_label
    want_html = args.html == "on" or (args.html == "auto" and _has_plotly())
    color_txt = color_label or "none"

    if args.png and args.layout in ("merged", "both"):
        plot_merged_fig(
            xs, ts, Z1, C1, A1, norm, cmap,
            title="Rocket effect: " + dvr_txt + " vs time, color = " + color_txt
                  + ", opacity = shard radius",
            z_label=dvr_label, x_label=x_label, t_label=t_min_lbl,
            color_label=color_label, alpha_note=alpha_note,
            out_png=f"{prefix}_merged.png", elev=args.elev, azim=args.azim,
            render=args.render,
        )
        if args.also2d:
            plot_merged_2d(xs, ts, Z1, C1, A1, norm, cmap,
                           title=dvr_label, z_label=dvr_label,
                           color_label=color_label, t_label=t_min_lbl,
                           out_png=f"{prefix}_merged_2d.png")

    if args.png and args.accel == "on" and np.isfinite(Z3).any():
        plot_merged_fig(
            xs, ts, Z3, C1, A1, norm, cmap,
            title="Rocket acceleration: " + acc_txt + " vs time, color = " + color_txt
                  + ", opacity = shard radius",
            z_label=acc_label, x_label=x_label, t_label=t_min_lbl,
            color_label=color_label, alpha_note=alpha_note,
            out_png=f"{prefix}_accel.png", elev=args.elev, azim=args.azim,
            render=args.render,
        )
        if args.also2d:
            plot_merged_2d(xs, ts, Z3, C1, A1, norm, cmap,
                           title=acc_label, z_label=acc_label,
                           color_label=color_label, t_label=t_min_lbl,
                           out_png=f"{prefix}_accel_2d.png")

    if args.png and args.layout in ("split", "both"):
        plot_surface_fig(
            xs, ts, Z1, C1, norm, cmap,
            title="Fig.1  Rocket effect: shard velocity change vs time "
                  f"(color = {color_txt})",
            z_label=dvr_label, x_label=x_label, t_label=t_min_lbl,
            color_label=color_label, out_png=f"{prefix}_fig1_dVR.png",
            elev=args.elev, azim=args.azim,
        )
        plot_surface_fig(
            xs, ts, Z2, C1, norm, cmap,
            title="Fig.2  Ablation: shard radius loss vs time (same color scale)",
            z_label=r"$\Delta r$  [mm]", x_label=x_label, t_label=t_min_lbl,
            color_label=color_label, out_png=f"{prefix}_fig2_dr.png",
            elev=args.elev, azim=args.azim,
        )
        if args.also2d:
            plot_2d_companion(xs, ts, Z1, C1, norm, cmap,
                              title=dvr_label, z_label=dvr_label,
                              color_label=color_label, t_label=t_min_lbl,
                              out_png=f"{prefix}_fig1_dVR_2d.png")
            plot_2d_companion(xs, ts, Z2, C1, norm, cmap,
                              title=r"$\Delta r$ [mm]", z_label=r"$\Delta r$ [mm]",
                              color_label=color_label, t_label=t_min_lbl,
                              out_png=f"{prefix}_fig2_dr_2d.png")

    if want_html:
        plot_html(xs, ts, Z1, C1, norm,
                  f"Merged  {dvr_txt} (color = {color_txt})", dvr_label, color_label,
                  t_min_lbl, f"{prefix}_merged.html")
        if args.accel == "on" and np.isfinite(Z3).any():
            plot_html(xs, ts, Z3, C1, norm,
                      f"Accel  {acc_txt} (color = {color_txt})", acc_label, color_label,
                      t_min_lbl, f"{prefix}_accel.html")
        if args.layout in ("split", "both"):
            plot_html(xs, ts, Z2, C1, norm,
                      "Fig.2  Δr (color = seen field)", "Δr [mm]", color_label,
                      t_min_lbl, f"{prefix}_fig2_dr.html")

    # ---------- 8. 文字摘要 ----------
    log("")
    log("=" * 74)
    log("摘要")
    log("=" * 74)
    log(f"  碎片数 (n_spi_tot)      : {len(i_spi)}")
    log(f"  时间步数                : {len(steps)}  ({steps[0]} .. {steps[-1]})")
    log(f"  参考步 (Δ 的零点)       : step {steps[0]}")
    log(f"  存活碎片数 (首/末)      : {int(n_alive[0])} / {int(n_alive[-1])}")
    log(f"  {dvr_txt} [m/s] 范围        : "
        f"[{np.nanmin(Z1) if np.isfinite(Z1).any() else float('nan'):.3f}, "
        f"{np.nanmax(Z1) if np.isfinite(Z1).any() else float('nan'):.3f}]")
    log(f"  Δr 范围                 : [{np.nanmin(Z2) if np.isfinite(Z2).any() else float('nan'):.4f}, "
        f"{np.nanmax(Z2) if np.isfinite(Z2).any() else float('nan'):.4f}] mm")
    if args.accel == "on":
        log(f"  {acc_txt} [m/s²] 范围       : "
            f"[{np.nanmin(Z3) if np.isfinite(Z3).any() else float('nan'):.3e}, "
            f"{np.nanmax(Z3) if np.isfinite(Z3).any() else float('nan'):.3e}]")
    if norm is not None:
        log(f"  颜色范围                : [{norm.vmin:.4g}, {norm.vmax:.4g}]  ({color_label})")
    log("=" * 74)
    return 0


def _has_plotly() -> bool:
    try:
        import plotly  # noqa: F401
        return True
    except Exception:  # noqa: BLE001
        return False


if __name__ == "__main__":
    sys.exit(main())
