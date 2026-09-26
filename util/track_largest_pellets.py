#!/usr/bin/env python3
import os
import re
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, LinearNDInterpolator

def load_time_mapping(live_data_file):
    """解析 macroscopic_vars.dat 以提取 step 到 time 的映射及单位换算常数"""
    time_map = {}
    time_scaling = None
    with open(live_data_file, 'r') as f:
        for line in f:
            if line.startswith('@sqrt_mu0_rho0:'):
                try:
                    time_scaling = float(line.strip().split()[1])
                except Exception:
                    pass
            elif line.startswith('@times:'):
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        step = int(parts[1])
                        time_val = float(parts[2])
                        time_map[step] = time_val
                    except ValueError:
                        continue
    if time_scaling is None:
        raise ValueError(
            f"FATAL: Could not find '@sqrt_mu0_rho0:' in {live_data_file}. "
            "This file must contain the time normalization factor. "
            "Ensure macroscopic_vars.dat is from the same simulation run."
        )
    return time_map, time_scaling

def load_abl_history(abl_file):
    """读取 abl_history.dat 提取第一次注入的消融率"""
    if not os.path.exists(abl_file):
        return None
    try:
        data = []
        with open(abl_file, 'r') as f:
            lines = f.readlines()
        
        start_idx = 1 if ('time' in lines[0] and 'abl_rate' in lines[0]) else 0
        
        for line in lines[start_idx:]:
            parts = line.split()
            if len(parts) >= 3:
                step = int(parts[0])
                time_val = float(parts[1])
                abl_rate_inj1 = float(parts[2])
                data.append((step, time_val, abl_rate_inj1))
        
        return pd.DataFrame(data, columns=['step', 'time', 'abl_rate_inj1'])
    except Exception as e:
        print(f"Failed to read {abl_file}: {e}")
        return None

def load_energies(live_data_file, time_scaling):
    """提取 macroscopic_vars.dat 中的 integrated_energies 数据
    返回 DataFrame: time_jorek, time_si, total_energy, magnetic, kin_par, kin_perp, thermal, e_thermal, i_thermal
    """
    if not os.path.exists(live_data_file):
        return None
    data = []
    with open(live_data_file, 'r') as f:
        for line in f:
            if line.startswith('@integrated_energies:') and not line.startswith('@integrated_energies_'):
                parts = line.strip().split()
                if len(parts) >= 9 and '%' not in parts[1]:
                    try:
                        t_jorek = float(parts[1])
                        total_energy = float(parts[2])
                        magnetic = float(parts[3])
                        kin_par = float(parts[4])
                        kin_perp = float(parts[5])
                        thermal = float(parts[6])
                        e_thermal = float(parts[7])
                        i_thermal = float(parts[8])
                        data.append((t_jorek, t_jorek * time_scaling, total_energy, magnetic, kin_par, kin_perp, thermal, e_thermal, i_thermal))
                    except ValueError:
                        pass
    if data:
        return pd.DataFrame(data, columns=['time_jorek', 'time_si', 'total_energy', 'magnetic', 'kin_par', 'kin_perp', 'thermal', 'e_thermal', 'i_thermal'])
    return None

def load_radiated_power(live_data_file, time_scaling):
    """提取 macroscopic_vars.dat 中的 radiated power 数据
    返回 DataFrame: time_si, radiated_power
    """
    if not os.path.exists(live_data_file):
        return None
    data = []
    with open(live_data_file, 'r') as f:
        for line in f:
            if line.startswith('@dissipative_terms:') and not line.startswith('@dissipative_terms_'):
                parts = line.strip().split()
                if len(parts) >= 7 and '%' not in parts[1]:
                    try:
                        t_jorek = float(parts[1])
                        # 0="@dissipative_terms:", 1="time", 2="Ohmic", 3="Frictional", 4="Perp_visc", 5="Par_visc", 6="Radiated"
                        rad_power = float(parts[6])
                        data.append((t_jorek * time_scaling, rad_power))
                    except ValueError:
                        pass
    if data:
        return pd.DataFrame(data, columns=['time_si', 'radiated_power'])
    return None

def get_target_q_psi_n(q_profile_file, target_q=-2.0):
    """(DEPRECATED) 从单个 q-profile 文件读取初始 q(Psi_N) 并插值固定 Psi_N。
    保留此函数仅用于向后兼容；新流程应使用 load_qprofile_evolution()。"""
    if not q_profile_file or not os.path.exists(q_profile_file):
        return None
    try:
        col_names = ["psi_n", "q"]
        df = pd.read_csv(q_profile_file, delim_whitespace=True, comment='#', names=col_names)
        
        df_sorted = df.drop_duplicates(subset=["q"]).sort_values(by="q")
        if target_q < df_sorted["q"].min() or target_q > df_sorted["q"].max():
            print(f"Target q={target_q} is out of bounds in {q_profile_file}.")
            return None
        
        f_psi_n = interp1d(df_sorted["q"], df_sorted["psi_n"])
        return float(f_psi_n(target_q))
    except Exception as e:
        print(f"Failed to read/interpolate {q_profile_file}: {e}")
        return None

def load_qprofile_evolution(q_profile_file, time_map, time_scaling, t0, target_q=-2.0):
    """解析合并的 q-profile 时间演化文件，返回每个时刻 q=target_q 所在的 Psi_N。
    
    输入文件格式（由 jorek2_postproc qprofile 命令在 for 循环中输出）：
        # time step #XXXXXX
        Psi_N   q
        ...     ...
        # time step #YYYYYY
        ...
    
    Parameters
    ----------
    q_profile_file : str
        合并后的 q-profile 文件路径
    time_map : dict
        step → jorek time 映射
    time_scaling : float
    t0 : float
        起始时间
    target_q : float
        目标 q 值
    
    Returns
    -------
    psi_n_at_q : np.ndarray (n_steps,) or None
        每个时间步的 Psi_N(q=target_q)
    time_si_at_q : np.ndarray (n_steps,) or None
    steps_at_q : np.ndarray (n_steps,) or None
    """
    if not q_profile_file or not os.path.exists(q_profile_file):
        return None, None, None
    
    blocks = []  # list of (step, time_si, psi_n_at_target)
    
    current_step = None
    psi_n_vals = []
    q_vals = []
    
    def process_block(step, pn, qv):
        if step is None or len(pn) < 2:
            return
        try:
            # 按 q 排序并去重
            order = np.argsort(qv)
            pn_sorted = np.array(pn)[order]
            q_sorted = np.array(qv)[order]
            unique_mask = np.concatenate([[True], np.diff(q_sorted) > 0])
            pn_unique = pn_sorted[unique_mask]
            q_unique = q_sorted[unique_mask]
            
            if len(q_unique) < 2:
                return
            if target_q < q_unique.min() or target_q > q_unique.max():
                return
            
            f_pn = interp1d(q_unique, pn_unique, kind='linear')
            psi_n_tgt = float(f_pn(target_q))
            real_time = time_map.get(step, step) * time_scaling
            if real_time >= t0:
                blocks.append((step, real_time, psi_n_tgt))
        except Exception:
            pass
    
    with open(q_profile_file, "r") as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if "time step" in line:
                    if current_step is not None:
                        process_block(current_step, psi_n_vals, q_vals)
                        psi_n_vals = []
                        q_vals = []
                    try:
                        step_str = line.split("#")[-1].strip()
                        current_step = int(step_str)
                    except ValueError:
                        pass
            else:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        psi_n_vals.append(float(parts[0]))
                        q_vals.append(float(parts[1]))
                    except ValueError:
                        pass
    
    if current_step is not None:
        process_block(current_step, psi_n_vals, q_vals)
    
    if not blocks:
        return None, None, None
    
    blocks.sort(key=lambda x: x[1])
    steps_arr = np.array([b[0] for b in blocks])
    times_arr = np.array([b[1] for b in blocks])
    psi_n_arr = np.array([b[2] for b in blocks])
    
    return psi_n_arr, times_arr, steps_arr

def build_rectangle_step_map(rectangle_dir):
    """扫描目录中由 jorek2_postproc rectangle 导出的 HDF5 文件，
    建立 step → file_path 映射。
    
    文件命名格式: exprs_Rmin*_Rmax*_Zmin*_Zmax*_phi*_sXXXXXX.h5
    """
    step_map = {}
    pattern = re.compile(r'_s(\d{6})\.h5$')
    for f in sorted(glob.glob(os.path.join(rectangle_dir, "exprs_Rmin*.h5"))):
        m = pattern.search(os.path.basename(f))
        if m:
            step = int(m.group(1))
            step_map[step] = f
    return step_map

def infer_bounds_from_fragment_csv(csv_dir, margin=0.1):
    """从 fragments_position_*.csv 自动推断 (RMIN, RMAX, ZMIN, ZMAX)。

    仅使用匹配文件名的数据，并对边界增加固定 margin。
    """
    all_files = glob.glob(os.path.join(csv_dir, "*.csv"))
    csv_files = sorted([f for f in all_files if re.match(r"fragments_position_\d+\.csv", os.path.basename(f))])
    if not csv_files:
        return None

    r_min = np.inf
    r_max = -np.inf
    z_min = np.inf
    z_max = -np.inf
    found = False

    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file, delim_whitespace=True, header=None, engine='python')
        except Exception:
            continue

        # 需要至少包含列: i_spi, R, Z
        if df.empty or len(df.columns) < 3:
            continue

        r_col = pd.to_numeric(df.iloc[:, 1], errors='coerce')
        z_col = pd.to_numeric(df.iloc[:, 2], errors='coerce')
        valid = r_col.notna() & z_col.notna()
        if not valid.any():
            continue

        r_vals = r_col[valid].values
        z_vals = z_col[valid].values
        r_min = min(r_min, float(np.min(r_vals)))
        r_max = max(r_max, float(np.max(r_vals)))
        z_min = min(z_min, float(np.min(z_vals)))
        z_max = max(z_max, float(np.max(z_vals)))
        found = True

    if not found:
        return None

    return (r_min - margin, r_max + margin, z_min - margin, z_max + margin)

def centers_to_edges(values):
    """将离散采样中心点转换为绘图边界。

    适用于非均匀网格：相邻中心点之间采用中点作为边界，
    首尾边界外推半个间隔。
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size == 0:
        return values
    if values.size == 1:
        delta = 0.5 if values[0] == 0 else abs(values[0]) * 0.5
        return np.array([values[0] - delta, values[0] + delta], dtype=float)

    edges = np.empty(values.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (values[:-1] + values[1:])
    edges[0] = values[0] - 0.5 * (values[1] - values[0])
    edges[-1] = values[-1] + 0.5 * (values[-1] - values[-2])
    return edges

# Cache for Psi_N interpolators to avoid reloading the same HDF5 repeatedly
_psi_interp_cache = {}  # step → LinearNDInterpolator
_MAX_CACHE_SIZE = 8

def get_psi_n_interpolator(step, step_map, bounds=None):
    """获取指定 step 的 Psi_N(R,Z) 插值器（带 LRU 缓存）。
    
    Parameters
    ----------
    step : int
    step_map : dict
        step → rectangle HDF5 file path
    bounds : tuple or None
        (RMIN, RMAX, ZMIN, ZMAX) 备用边界
    
    Returns
    -------
    LinearNDInterpolator or None
    """
    global _psi_interp_cache
    
    if step in _psi_interp_cache:
        return _psi_interp_cache[step]
    
    # 找最接近的步骤（前后各取一个候选，取最近的）
    if step not in step_map:
        available = sorted(step_map.keys())
        if not available:
            return None
        # 找最近的步骤
        idx = np.searchsorted(available, step)
        candidates = []
        if idx < len(available):
            candidates.append(available[idx])
        if idx > 0:
            candidates.append(available[idx - 1])
        best_step = min(candidates, key=lambda s: abs(s - step))
        # 只在差距合理时使用（差距 > 10 步则警告）
        if abs(best_step - step) > 10:
            print(f"Warning: step {step} has no rectangle, nearest is {best_step} (Δ={abs(best_step-step)})")
        step = best_step
    
    file_path = step_map[step]
    R_data, Z_data, Psi_N_data, _ = load_psi_n_h5(file_path, bounds=bounds)
    interpolator = LinearNDInterpolator(list(zip(R_data, Z_data)), Psi_N_data)
    
    # LRU eviction
    if len(_psi_interp_cache) >= _MAX_CACHE_SIZE:
        oldest = next(iter(_psi_interp_cache))
        del _psi_interp_cache[oldest]
    _psi_interp_cache[step] = interpolator
    
    return interpolator

def get_te_evolution_at_moving_psi(psi_n_lookup, te_profile_files, time_map, time_scaling, t0, limit_ev=10.0):
    """提取移动磁面（Psi_N 随时间变化）上的局部 Te 演变，并寻找局部 TQ 时间。
    
    Parameters
    ----------
    psi_n_lookup : dict or None
        step → Psi_N 映射（由 load_qprofile_evolution 构建）
    te_profile_files : list of str
        average 输出的 profile 文件路径列表
    time_map : dict
    time_scaling : float
    t0 : float
    limit_ev : float
        Te 低于此阈值 (eV) 视为局部 TQ
        
    Returns
    -------
    t_tq : float or None
    t_tq_step : int or None
    te_df : pd.DataFrame
    psi_n_track : list of (step, psi_n_used)
    """
    if psi_n_lookup is None or not te_profile_files:
        return None, None, pd.DataFrame(), []
    
    # 构建 step → Psi_N 快速查询
    if isinstance(psi_n_lookup, dict):
        psi_n_dict = psi_n_lookup
    else:
        # 假设是 (psi_n_arr, times_arr, steps_arr) 三元组
        psi_n_arr, times_arr, steps_arr = psi_n_lookup
        psi_n_dict = {int(s): float(p) for s, p in zip(steps_arr, psi_n_arr)}
    
    te_history = []
    psi_n_track = []
    t_tq = None
    t_tq_step = None
    
    for f in te_profile_files:
        current_step = None
        data_rows = []
        
        def process_block(step, rows):
            nonlocal t_tq, t_tq_step
            if step not in psi_n_dict:
                return
            psi_target = psi_n_dict[step]
            try:
                df = pd.DataFrame(rows, columns=['psi_n', 'Te'])
                df = df.dropna().sort_values(by="psi_n")
                if not df.empty and psi_target >= df["psi_n"].min() and psi_target <= df["psi_n"].max():
                    f_te = interp1d(df["psi_n"], df["Te"])
                    te_val = max(float(f_te(psi_target)), 1e-4)
                    real_time = time_map.get(step, step) * time_scaling
                    if real_time >= t0:
                        te_history.append((step, real_time, te_val))
                        psi_n_track.append((step, psi_target))
                        if te_val <= limit_ev and t_tq is None:
                            t_tq = float(real_time)
                            t_tq_step = step
            except Exception:
                pass

        with open(f, "r") as fp:
            for line in fp:
                line = line.strip()
                if not line: continue
                if line.startswith("#"):
                    if "time step" in line:
                        if current_step is not None and data_rows:
                            process_block(current_step, data_rows)
                            data_rows = []
                        try:
                            step_str = line.split("#")[-1].strip()
                            current_step = int(step_str)
                        except ValueError:
                            pass
                else:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            data_rows.append([float(parts[0]), float(parts[1])])
                        except ValueError:
                            pass
                        
        if current_step is not None and data_rows:
            process_block(current_step, data_rows)
            
    te_df = pd.DataFrame(te_history, columns=["step", "time_si", "Te"])
    if not te_df.empty:
        te_df = te_df.sort_values(by="time_si")
    return t_tq, t_tq_step, te_df, psi_n_track

def load_averaged_nimp(profile_files, time_map, time_scaling, t0):
    """
    解析 average 命令输出的径向 profile 文件，提取 nimp(Psi_N_current, t) 二维分布。
    
    输入文件包含 3 列: Psi_N_current, T_e, nimp
    直接使用当前时刻的 Psi_N，不做初始磁面映射。
    
    Parameters
    ----------
    profile_files : list of str
    time_map : dict
    time_scaling : float
    t0 : float
    
    Returns
    -------
    nimp_2d : np.ndarray (n_times, n_psi) or None
    psi_bins : np.ndarray (n_psi,) or None
    time_si : np.ndarray (n_times,) or None
    """
    if not profile_files:
        return None, None, None
    
    all_blocks = []  # list of (step, time_si, psi_arr, nimp_arr)
    
    for f in profile_files:
        current_step = None
        data_rows = []
        
        def commit_block(step, rows):
            if step is None or not rows:
                return
            try:
                cols = ['psi_n_cur', 'Te', 'nimp']
                df = pd.DataFrame(rows, columns=cols)
                df = df.dropna()
                if df.empty or len(df) < 2:
                    return
                if df['nimp'].max() == df['nimp'].min():
                    return
                real_time = time_map.get(step, step) * time_scaling
                if real_time >= t0:
                    all_blocks.append((step, real_time, df['psi_n_cur'].values, df['nimp'].values))
            except Exception:
                pass
        
        with open(f, "r") as fp:
            for line in fp:
                line = line.strip()
                if not line: continue
                if line.startswith("#"):
                    if "time step" in line:
                        if current_step is not None and data_rows:
                            commit_block(current_step, data_rows)
                            data_rows = []
                        try:
                            step_str = line.split("#")[-1].strip()
                            current_step = int(step_str)
                        except ValueError:
                            pass
                else:
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            data_rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
                        except ValueError:
                            pass
        if current_step is not None and data_rows:
            commit_block(current_step, data_rows)
    
    if not all_blocks:
        return None, None, None
    
    # 按时间排序
    all_blocks.sort(key=lambda x: x[1])
    
    # 统一 Psi_N 网格：强制升序 (0→1) 并插值到统一参考网格
    sorted_blocks = []
    for step, t_si, psi_arr, nimp_arr in all_blocks:
        order = np.argsort(psi_arr)
        sorted_blocks.append((step, t_si, psi_arr[order], nimp_arr[order]))
    all_blocks = sorted_blocks
    
    psi_all_min = min(b[2][0] for b in all_blocks)
    psi_all_max = max(b[2][-1] for b in all_blocks)
    n_psi_median = int(np.median([len(b[2]) for b in all_blocks]))
    psi_bins = np.linspace(max(psi_all_min, 0.0), min(psi_all_max, 1.0), n_psi_median)
    
    nimp_rows = []
    time_rows = []
    
    for step, t_si, psi_arr, nimp_arr in all_blocks:
        # 沿 Psi 轴线性外推缺失边界区域
        try:
            f_interp = interp1d(psi_arr, nimp_arr, kind='linear',
                               bounds_error=False, fill_value='extrapolate')
        except ValueError:
            # 旧版 scipy 不支持 'extrapolate'，回退到边界常值
            print('Warning: interp1d does not support "extrapolate", falling back to edge values for out-of-bounds.')
            f_interp = interp1d(psi_arr, nimp_arr, kind='linear',
                               bounds_error=False,
                               fill_value=(float(nimp_arr[0]), float(nimp_arr[-1])))
        nimp_interped = f_interp(psi_bins)
        # 裁剪线性外推可能产生的物理上不可能的负值
        nimp_interped = np.maximum(nimp_interped, 0.0)
        nimp_rows.append(nimp_interped)
        time_rows.append(t_si)
    
    nimp_2d = np.array(nimp_rows)
    time_si = np.array(time_rows)
    
    if nimp_2d.max() <= 0 or np.allclose(nimp_2d, nimp_2d[0], atol=1e-30):
        return None, None, None
    
    return nimp_2d, psi_bins, time_si

def load_psi_n_h5(filename, bounds=None):
    """
    专门读取 jorek2_postproc 导出的 rectangle (.h5) 2D网格，返回所有点集供插值。
    
    Returns
    -------
    R_pts, Z_pts, Psi_N_pts, Psi_pts : np.ndarray
        Psi_pts 在仅有 Psi_N 表达式时为 None，用于向后兼容。
    """
    if not filename.endswith(".h5"):
        raise ValueError(f"Expected a .h5 file containing the Psi_N grid, got: {filename}")
    
    import h5py
    with h5py.File(filename, 'r') as f:
        variables = [v.decode('utf-8') if isinstance(v, bytes) else str(v) for v in f['variables'][()]]
        print(f"HDF5 Variables detected: {variables}")
        v = f['values'][()]
        
        n_var = v.shape[0]
        psi_pts_extra = None  # Psi field if available
        
        # --- 提取 R, Z 坐标轴 ---
        if 'axis1' in f:
            ax1 = f['axis1'][()]
            ax2 = f['axis2'][()]
        elif 'grid1' in f:
            ax1 = f['grid1'][()]
            ax2 = f['grid2'][()]
        else:
            # 无坐标轴：用 bounds 或默认
            if bounds is not None and len(bounds) == 4:
                rmin, rmax, zmin, zmax = bounds
                ax1 = np.linspace(rmin, rmax, v.shape[1])
                ax2 = np.linspace(zmin, zmax, v.shape[2])
            else:
                ax1 = np.linspace(0, 1, v.shape[1])
                ax2 = np.linspace(0, 1, v.shape[2])
            # 无坐标轴时原样使用（无法判断R/Z方向）
            R_ax, Z_ax = ax1, ax2
            R_mesh, Z_mesh = np.meshgrid(R_ax, Z_ax, indexing='xy')
            R_pts = R_mesh.flatten()
            Z_pts = Z_mesh.flatten()
            psi_pts = v[0].flatten() if n_var >= 1 else np.zeros_like(R_pts)
            if n_var >= 2:
                psi_pts_extra = v[1].flatten()
            mask = (R_pts != 0) | (Z_pts != 0)
            return R_pts[mask], Z_pts[mask], psi_pts[mask], (psi_pts_extra[mask] if psi_pts_extra is not None else None)
        
        # 有坐标轴：判断哪个是R，哪个是Z
        if ax1.mean() > ax2.mean() and ax1.min() >= 0:
            R_ax, Z_ax = ax1, ax2
        else:
            R_ax, Z_ax = ax2, ax1
        
        # indexing='xy' so that R varies fastest when flattened
        R_mesh, Z_mesh = np.meshgrid(R_ax, Z_ax, indexing='xy')
        R_pts = R_mesh.flatten()
        Z_pts = Z_mesh.flatten()
        
        if n_var >= 1:
            psi_pts = v[0].flatten()
        else:
            psi_pts = np.zeros_like(R_pts)
        if n_var >= 2:
            psi_pts_extra = v[1].flatten()
        
        mask = (R_pts != 0) | (Z_pts != 0)
        return (R_pts[mask], Z_pts[mask], psi_pts[mask],
                psi_pts_extra[mask] if psi_pts_extra is not None else None)

def main():
    parser = argparse.ArgumentParser(description="Track the average penetration of the largest N pellet fragments.")
    parser.add_argument("--csv_dir", type=str, required=True, help="Directory containing fragments_position_*.csv from convert2diagno_spi.py")
    parser.add_argument("--rectangle_dir", type=str, required=True, help="Directory containing time-dependent Psi_N(R,Z) rectangle HDF5 files from jorek2_postproc")
    parser.add_argument("--live_data", type=str, default='macroscopic_vars.dat', help="Path to macroscopic_vars.dat to extract physical time")
    parser.add_argument("--abl_file", type=str, default="abl_history.dat", help="Path to abl_history.dat produced by jorek2_diagno_spi")
    parser.add_argument("--start_time", type=float, default=None, help="Optional start time (in SI seconds) to plot from")
    parser.add_argument("--max_time_ms", type=float, default=2.0, help="Maximum plotted time in ms from the start_time (default: 2.0 ms)")
    parser.add_argument("--q_profile", type=str, default=None, help="Path to merged q-profile evolution dat file (from jorek2_postproc qprofile in for loop)")
    parser.add_argument("--target_q", type=float, default=-2.0, help="Target q-value to track evolving rational surface")
    parser.add_argument("--out_csv", type=str, default="pellet_average_trajectory.csv", help="Output trajectory file")
    parser.add_argument("--out_plot", type=str, default="penetration_trajectory.png", help="Output visualization plot")
    parser.add_argument("--bounds", type=float, nargs=4, metavar=('RMIN', 'RMAX', 'ZMIN', 'ZMAX'), default=None, help="Bounding box of the rectangle grid if needed")
    parser.add_argument("--te_profile", type=str, default=None, help="Path to averaged Te profile .dat file (from jorek2_postproc average command)")
    parser.add_argument("--nimp_profile", type=str, default=None, help="Path to averaged nimp profile .dat file (defaults to --te_profile if not given)")
    parser.add_argument("--debug", type=int, default=0, help="Enable debug output when set to 1")
    
    args = parser.parse_args()
    debug_enabled = int(args.debug) == 1

    # 0. 载入时间映射表
    time_map = {}
    time_scaling = 1.0
    if not args.live_data or not os.path.exists(args.live_data):
        raise FileNotFoundError(
            f"FATAL: macroscopic_vars.dat ('{args.live_data}') not found. "
            "This file is required for physical time conversion. "
            "Ensure the JOREK simulation has produced it before running diagnostics."
        )
    time_map, time_scaling = load_time_mapping(args.live_data)
    print(f"Loaded {len(time_map)} time steps from {args.live_data} (Scaling factor: {time_scaling:.4e})")

    # 1. 构建时间依赖的 Psi_N(R,Z) 插值器映射
    step_map = build_rectangle_step_map(args.rectangle_dir)
    if not step_map:
        print(f"No rectangle HDF5 files (exprs_Rmin*.h5) found in {args.rectangle_dir}")
        return
    print(f"Found {len(step_map)} rectangle HDF5 files in {args.rectangle_dir} (steps {min(step_map.keys())}–{max(step_map.keys())})")

    # 若未显式提供 bounds，则自动从 fragment CSV 推断边界，
    # 以避免某些 HDF5 缺少坐标轴时回退到 [0,1] 导致插值失效。
    if args.bounds is None:
        auto_bounds = infer_bounds_from_fragment_csv(args.csv_dir, margin=0.1)
        if auto_bounds is not None:
            args.bounds = auto_bounds
            print("Auto-inferred bounds from fragment CSVs: "
                  f"R=[{args.bounds[0]:.4f}, {args.bounds[1]:.4f}], "
                  f"Z=[{args.bounds[2]:.4f}, {args.bounds[3]:.4f}]")
        else:
            print("Warning: failed to infer bounds from fragment CSVs; \
falling back to HDF5 axis (or [0,1] if axis missing).")
    
    # 预热第一个插值器以验证数据
    first_step = min(step_map.keys())
    test_interp = get_psi_n_interpolator(first_step, step_map, bounds=args.bounds)
    print(f"Psi_N interpolator for step {first_step} ready.")

    # 2. 遍历并严格过滤所需前缀的 CSV (或 dat 重命名) 文件
    all_files = glob.glob(os.path.join(args.csv_dir, "*.csv"))
    csv_files = sorted([f for f in all_files if re.match(r"fragments_position_\d+\.csv", os.path.basename(f))])
    
    if not csv_files:
        print(f"No fragments_position_XXXXXX.csv files found in {args.csv_dir}")
        return

    time_series = []
    # 按照 Fortran 中 write(20,'(i7,3f12.3,e14.6,f12.3)') 输出定义的各列数据标题
    col_names = ['i_spi', 'spi_R', 'spi_Z', 'spi_phi', 'spi_radius', 'spi_species']
    
    for csv_file in csv_files:
        # 尝试从文件名萃取 step
        step_str = os.path.basename(csv_file).replace("fragments_position_", "").replace(".csv", "")
        # print(f"Processing file: {csv_file} (step: {step_str})")
        try:
            step = int(step_str)
        except ValueError:
            step = -1

        # 读取由于直接从 .dat 重命名的空格分割文件，不存在 header
        try:
            # 兼容多空格对齐并强制赋标题
            df = pd.read_csv(csv_file, delim_whitespace=True, header=None, names=col_names, engine='python')
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            continue
            
        if df.empty or 'spi_radius' not in df.columns:
            print(f"Skipping file: {csv_file}, no valid data.")
            continue

        # 选取所有半径大于 0 的碎片
        df_active = df[df['spi_radius'] > 0]
        
        if df_active.empty:
            # 如果没有活着的碎块，设定半径为 0
            avg_R, avg_Z, avg_radius, total_mass = df['spi_R'].mean(), df['spi_Z'].mean(), 0.0, 0.0
            psi_n_val, psi_n_std = np.nan, 0.0
        else:
            # 以体积（半径的三次方）为权重计算质心位置
            volumes = df_active['spi_radius']**3
            total_mass = volumes.sum()
            avg_R = (df_active['spi_R'] * volumes).sum() / total_mass
            avg_Z = (df_active['spi_Z'] * volumes).sum() / total_mass
            # 使用总体积等效半径来指示消融进程（解决从小碎片先死导致的算术平均值假象）
            avg_radius = total_mass**(1/3)
        
            # 计算每个碎片的精确 Psi_N 并求质量加权标准差
            # 使用当前时间步的 Psi_N(R,Z) 映射（而非初始平衡）
            step_interp = get_psi_n_interpolator(step, step_map, bounds=args.bounds)
            fragment_R = df_active['spi_R'].values
            fragment_Z = df_active['spi_Z'].values
            fragment_psi = step_interp(fragment_R, fragment_Z) if step_interp is not None else np.full(len(fragment_R), np.nan)
            
            valid_mask = ~np.isnan(fragment_psi)
            if valid_mask.sum() == 0:
                psi_n_val = np.nan
                psi_n_std = 0.0
            else:
                frag_psi_valid = fragment_psi[valid_mask]
                frag_vol_valid = volumes.values[valid_mask]
                mass_valid = frag_vol_valid.sum()
                
                if mass_valid > 0:
                    psi_n_val = float(np.average(frag_psi_valid, weights=frag_vol_valid))
                    variance = np.average((frag_psi_valid - psi_n_val)**2, weights=frag_vol_valid)
                    psi_n_std = float(np.sqrt(variance))
                else:
                    psi_n_val = np.nan
                    psi_n_std = 0.0

        if time_map:
            jorek_time = time_map.get(step, None)
            real_time = jorek_time * time_scaling if jorek_time is not None else None
        else:
            real_time = step

        time_series.append({
            'step': step,
            'time_si': real_time,
            'avg_R': avg_R,
            'avg_Z': avg_Z,
            'avg_radius': avg_radius,
            'psi_n': psi_n_val,
            'psi_n_std': psi_n_std
        })

    if not time_series:
        print("No valid data processed.")
        return

    # 3. 输出追踪数据结果
    res_df = pd.DataFrame(time_series)
    
    # 剔除无法映射时间的行
    initial_len = len(res_df)
    res_df = res_df.dropna(subset=['time_si'])
    if len(res_df) < initial_len:
        print(f"Dropped {initial_len - len(res_df)} steps due to missing time_si mapping.")
        
    if res_df.empty:
        print("No valid data left after dropping missing time steps.")
        return

    res_df = res_df.sort_values('step')
    res_df.to_csv(args.out_csv, index=False)
    print(f"Averaged trajectory saved to {args.out_csv}")

    # 4. 可视化：时间 vs 距离(穿深) 以及 选项式的消融率
    abl_df = load_abl_history(args.abl_file)

    if abl_df is not None and not abl_df.empty:
        if 'time_si' not in abl_df.columns:
            abl_df['time_si'] = abl_df['time'] * time_scaling
            
        # 如果未指定起始时间，则根据消融率动态寻找起始点
        if args.start_time is None:
            active_idx = abl_df[abl_df['abl_rate_inj1'] > 0].index
            if not active_idx.empty:
                first_idx = active_idx.min()
                start_idx = max(0, first_idx - 1)  # 往前推一步保留刚等于0时的时刻
                args.start_time = float(abl_df.loc[start_idx, 'time_si'])
                print(f"Auto-detected ablation start time: {args.start_time:.6e} s")

    # 应用过滤起始时间 (剔除垃圾时间)
    if args.start_time is not None:
        res_df = res_df[res_df['time_si'] >= args.start_time]
        if abl_df is not None and not abl_df.empty:
            abl_df = abl_df[abl_df['time_si'] >= args.start_time]

    t0 = args.start_time if args.start_time is not None else res_df['time_si'].min()

    energy_df = load_energies(args.live_data, time_scaling)
    if energy_df is not None and not energy_df.empty:
        energy_df = energy_df.iloc[:-1]  # 取 0 到 -2 的数据区间（去掉最后一行）
        if args.start_time is not None:
            energy_df = energy_df[energy_df['time_si'] >= args.start_time]

    rad_power_df = load_radiated_power(args.live_data, time_scaling)
    if rad_power_df is not None and not rad_power_df.empty:
        if args.start_time is not None:
            rad_power_df = rad_power_df[rad_power_df['time_si'] >= args.start_time]

    has_abl = abl_df is not None and not abl_df.empty
    has_egy = energy_df is not None and not energy_df.empty
    has_rad = rad_power_df is not None and not rad_power_df.empty
    
    # 加载 q-profile 时间演化，获取每个时刻 q=target_q 所在的 Psi_N
    psi_n_at_q_arr, time_at_q_arr, steps_at_q_arr = load_qprofile_evolution(
        args.q_profile, time_map, time_scaling, t0, target_q=args.target_q)
    
    has_q_evo = psi_n_at_q_arr is not None and len(psi_n_at_q_arr) > 0
    if has_q_evo:
        print(f"Loaded q-profile evolution: {len(psi_n_at_q_arr)} time steps, "
              f"Psi_N(q={args.target_q}) from {psi_n_at_q_arr[0]:.4f} to {psi_n_at_q_arr[-1]:.4f}")
        # 构建 step → Psi_N(q=target_q) 查询字典
        psi_n_at_q_dict = {int(s): float(p) for s, p in zip(steps_at_q_arr, psi_n_at_q_arr)}
        # 初始时刻的 q=-2 Psi_N（用于 subplot 1 参考线）
        target_q_psi_initial = float(psi_n_at_q_arr[0])
    else:
        psi_n_at_q_dict = {}
        target_q_psi_initial = None
        # 降级: 尝试从 qprofile 文件获取固定 Psi_N
        target_q_psi_initial = get_target_q_psi_n(args.q_profile, args.target_q)
        if target_q_psi_initial is not None:
            print(f"Interpolated Psi_N for q={args.target_q} is {target_q_psi_initial:.6f} (static, no evolution data)")
    
    # 构建时间插值函数 Psi_N(q=target_q, t) 用于 fragment 穿深事件检测
    if has_q_evo:
        f_psin_at_q_vs_t = interp1d(time_at_q_arr, psi_n_at_q_arr, kind='linear',
                                     bounds_error=False, fill_value='extrapolate')
    else:
        f_psin_at_q_vs_t = None

    # 计算标记事件时间
    t_r0 = None 
    t_r0_step = None
    radius_drop_idx = res_df[res_df['avg_radius'] <= 1e-10].first_valid_index()
    if radius_drop_idx is not None:
        t_r0 = float(res_df.loc[radius_drop_idx, 'time_si'])
        t_r0_step = int(res_df.loc[radius_drop_idx, 'step'])

    t_10_percent_thermal = None
    t_10_percent_step = None
    if has_egy:
        initial_thermal = float(energy_df['thermal'].iloc[0])
        th_drop_idx = energy_df[energy_df['thermal'] <= 0.1 * initial_thermal].first_valid_index()
        if th_drop_idx is not None:
            t_10_percent_thermal = float(energy_df.loc[th_drop_idx, 'time_si'])
            _closest_idx = (res_df['time_si'] - t_10_percent_thermal).abs().idxmin()
            t_10_percent_step = int(res_df.loc[_closest_idx, 'step'])

    # 碎片穿深到达移动 q=-2 面的时刻
    t_target_q = None
    t_target_q_step = None
    if has_q_evo and f_psin_at_q_vs_t is not None:
        # 对每个时刻，比较碎片 Psi_N 与移动 q=-2 面的 Psi_N
        for _, row in res_df.iterrows():
            t_si = row['time_si']
            pn_frag = row['psi_n']
            if np.isnan(pn_frag):
                continue
            pn_q = float(f_psin_at_q_vs_t(t_si))
            if pn_frag <= pn_q:
                t_target_q = float(t_si)
                t_target_q_step = int(row['step'])
                break
    elif target_q_psi_initial is not None:
        target_q_idx = res_df[res_df['psi_n'] <= target_q_psi_initial].first_valid_index()
        if target_q_idx is not None:
            t_target_q = float(res_df.loc[target_q_idx, 'time_si'])
            t_target_q_step = int(res_df.loc[target_q_idx, 'step'])
            
    _t0_idx = (res_df['time_si'] - t0).abs().idxmin()
    t0_step = int(res_df.loc[_t0_idx, 'step'])
    
    # 查找 Te 文件路径
    if args.te_profile:
        te_files = [args.te_profile]
    else:
        _candidates = [
            "auto_diagno_spi/postproc/exprs_averaged_merged.dat",
            "postproc/exprs_averaged_merged.dat",
        ]
        te_files = []
        for _cand in _candidates:
            if os.path.exists(_cand):
                te_files = [_cand]
                break
        if not te_files:
            _combined = sorted(glob.glob("auto_diagno_spi/postproc/exprs_averaged_s*..*.dat"))
            if not _combined:
                _combined = sorted(glob.glob("postproc/exprs_averaged_s*..*.dat"))
            if _combined:
                te_files = [_combined[-1]]
            else:
                _all = sorted(glob.glob("auto_diagno_spi/postproc/exprs_averaged_s*.dat"))
                if not _all:
                    _all = sorted(glob.glob("postproc/exprs_averaged_s*.dat"))
                te_files = [f for f in _all if ".." not in f]
    
    # 查找局部温度 TQ —— 跟踪移动 q=-2 面上的 Te
    t_local_tq, t_local_tq_step, te_df, psi_n_track = get_te_evolution_at_moving_psi(
        psi_n_at_q_dict, te_files, time_map, time_scaling, t0, limit_ev=10.0)
    
    if psi_n_track:
        print(f"Te tracked at evolving q={args.target_q} surface: {len(te_df)} points")
    
    # 加载 nimp(Psi_N, t) 二维分布用于热力图（使用当前时刻 Psi_N，不做映射）
    if args.nimp_profile:
        print(f"Using specified nimp profile file: {args.nimp_profile}")
        nimp_files = [args.nimp_profile]
    else:
        nimp_files = te_files
    nimp_2d, psi_bins_nimp, time_si_nimp = load_averaged_nimp(
        nimp_files, time_map, time_scaling, t0)
    if nimp_2d is not None:
        nimp_raw_excess = nimp_2d
        nimp_raw_excess = np.maximum(nimp_raw_excess, 0.0)
        nimp_display = np.maximum(nimp_raw_excess, 1e-10) * 1e-20
        nimp_threshold = 0.05 * nimp_raw_excess.max() if nimp_raw_excess.max() > 0 else np.inf
        cloud_front_psi = np.array([
            psi_bins_nimp[np.where(row > nimp_threshold)[0][0]]
            if np.any(row > nimp_threshold) else np.nan
            for row in nimp_raw_excess
        ])
        cloud_peak_psi = np.array([
            psi_bins_nimp[np.argmax(row)] if row.max() > 0 else np.nan
            for row in nimp_raw_excess
        ])
        has_nimp = True
        print(f"Loaded nimp 2D distribution: {nimp_2d.shape[0]} time steps × {nimp_2d.shape[1]} Psi_N bins")
    else:
        has_nimp = False
        print("No valid nimp data found, falling back to ±1σ envelope.")

    # Determine output directory for text summaries (same as plot)
    out_dir = os.path.dirname(os.path.abspath(args.out_plot))

    # --- Build and write EVENT TIMINGS ---
    event_lines = []
    event_lines.append("="*50)
    event_lines.append("                EVENT TIMINGS")
    event_lines.append("="*50)
    def format_event(name, t_val, step_val):
        if t_val is not None:
            return f"- {name:<27}: t = {(t_val - t0)*1000:8.4f} ms (relative), Step ~ {step_val}"
        else:
            return f"- {name:<27}: Not Reached"
    def emit_event(name, t_val, step_val):
        line = format_event(name, t_val, step_val)
        print(line)
        event_lines.append(line)
    emit_event("Injection / Start", t0, t0_step)
    if has_q_evo or target_q_psi_initial is not None:
        emit_event(f"Reached evolving q={args.target_q}", t_target_q, t_target_q_step)
        emit_event(f"Local TQ (Te<10eV at q={args.target_q})", t_local_tq, t_local_tq_step)
    emit_event("Thermal Quench (~90% drop)", t_10_percent_thermal, t_10_percent_step)
    emit_event("Complete Ablation", t_r0, t_r0_step)
    event_lines.append("="*50)
    event_text = "\n".join(event_lines) + "\n"
    print()
    summary_path = os.path.join(out_dir, "diagnostic_summary.txt")
    with open(summary_path, "w") as f:
        f.write(event_text)
    print(f"[INFO]   Event timings saved to {summary_path}")

    if has_abl:
        abl_lines = []
        abl_lines.append("="*50)
        abl_lines.append("             INTEGRATED ABLATION")
        abl_lines.append("="*50)

        # 统一积分上限：优先使用局部 Te<10eV 的 TQ，其次回退到 90% 热量下降时刻
        integration_tq_time = t_local_tq if t_local_tq is not None else t_10_percent_thermal
        abl_lines.append(format_event("Ablation/Radiation TQ Upper Bound", integration_tq_time, t_local_tq_step if t_local_tq is not None else t_10_percent_step))
        
        def calculate_integrated_ablation(end_time):
            if end_time is None or end_time < t0:
                return None
            mask = (abl_df['time_si'] >= t0) & (abl_df['time_si'] <= end_time)
            df_sub = abl_df[mask]
            if len(df_sub) < 2:
                return 0.0
            return np.trapz(df_sub['abl_rate_inj1'], df_sub['time_si'])

        abl_to_tq = calculate_integrated_ablation(integration_tq_time)
        abl_to_r0 = calculate_integrated_ablation(t_r0)
        abl_total = calculate_integrated_ablation(abl_df['time_si'].max())
        
        if abl_to_tq is not None:
            abl_lines.append(f"- To Thermal Quench     : {abl_to_tq:.4e}")
        if abl_to_r0 is not None:
            abl_lines.append(f"- To Complete Ablation  : {abl_to_r0:.4e}")
        abl_lines.append(f"- Total Ablated (Max)   : {abl_total:.4e}")
        abl_lines.append("="*50)
        abl_text = "\n".join(abl_lines) + "\n"
        print(abl_text)
        with open(summary_path, "a") as f:
            f.write("\n" + abl_text)
        print(f"[INFO]   Integrated ablation appended to {summary_path}")

    if has_rad:
        rad_lines = []
        rad_lines.append("="*50)
        rad_lines.append("         INTEGRATED RADIATED ENERGY")
        rad_lines.append("="*50)
        
        def calculate_integrated_radiation(end_time):
            if end_time is None or end_time < t0:
                return None
            mask = (rad_power_df['time_si'] >= t0) & (rad_power_df['time_si'] <= end_time)
            df_sub = rad_power_df[mask]
            if len(df_sub) < 2:
                return 0.0
            return np.trapz(df_sub['radiated_power'], df_sub['time_si'])

        rad_to_tq = calculate_integrated_radiation(integration_tq_time)
        rad_to_r0 = calculate_integrated_radiation(t_r0)
        rad_total = calculate_integrated_radiation(rad_power_df['time_si'].max())
        
        if rad_to_tq is not None:
            rad_lines.append(f"- To Thermal Quench     : {rad_to_tq:.4e} J")
        if rad_to_r0 is not None:
            rad_lines.append(f"- To Complete Ablation  : {rad_to_r0:.4e} J")
        rad_lines.append(f"- Total Radiated (Max)  : {rad_total:.4e} J")
        rad_lines.append("="*50)
        rad_text = "\n".join(rad_lines) + "\n"
        print(rad_text)
        with open(summary_path, "a") as f:
            f.write("\n" + rad_text)
        print(f"[INFO]   Integrated radiated energy appended to {summary_path}")

    num_subplots = 2
    if has_abl:
        num_subplots += 1
    if has_egy:
        num_subplots += 1

    plt.rcParams.update({
        'font.size': 14,
        'axes.labelsize': 16,
        'axes.titlesize': 18,
        'legend.fontsize': 14,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'lines.linewidth': 3
    })

    fig, axes = plt.subplots(num_subplots, 1, figsize=(10, 4*num_subplots), sharex=True)
    if num_subplots == 1:
        axes = [axes]
        
    x_val = res_df['time_si'] - t0 # s
    x_label = "Simulation Time (ms)" if time_scaling != 1.0 else ("Simulation Time" if args.live_data else "Simulation Step")
    
    # 画事件竖线的辅助函数
    def draw_event_lines(ax):
        # if t_psi95 is not None and t_psi95 >= t0:
        #     tx = (t_psi95 - t0)*1e3
        #     ax.axvline(tx, color='black', linestyle='-.', alpha=0.8, linewidth=2.0)
        if t_target_q is not None and t_target_q >= t0:
            tx = (t_target_q - t0)*1e3
            ax.axvline(tx, color='#9467bd', linestyle='-.', alpha=0.8, linewidth=2.0)
        
        # 优先画出基于局部表面温度10eV标记的TQ线，如果没有数据退化回90%下降的
        t_tq_plot = t_local_tq if (t_local_tq is not None) else t_10_percent_thermal
        if t_tq_plot is not None and t_tq_plot >= t0:
            tx = (t_tq_plot - t0)*1e3
            ax.axvline(tx, color='darkorange', linestyle='-.', alpha=0.8, linewidth=2.0)
            
        if t_r0 is not None and t_r0 >= t0:
            tx = (t_r0 - t0)*1e3
            ax.axvline(tx, color='grey', linestyle='-.', alpha=0.8, linewidth=2.0)

    ax_ind = 0

    # 第一张子图: 穿深 + nimp 热力图
    ax_depth = axes[ax_ind]
    
    if has_nimp:
        # 使用线性色标（替换原先的对数色标）
        time_ms_nimp = (time_si_nimp - t0) * 1e3

        # 裁剪到绘图时间范围
        mask_t = (time_ms_nimp >= 0) & (time_ms_nimp <= args.max_time_ms)
        if mask_t.sum() >= 2:
            t_plot = time_ms_nimp[mask_t]
            nimp_plot = nimp_display[mask_t, :]
            time_edges = centers_to_edges(t_plot)
            psi_edges = centers_to_edges(psi_bins_nimp)

            dt = np.diff(t_plot)
            if debug_enabled and len(dt) > 0:
                dt_min = float(np.min(dt))
                dt_max = float(np.max(dt))
                if dt_min > 0 and dt_max / dt_min > 1.05:
                    print(f"nimp heatmap uses non-uniform time spacing: dt_min={dt_min:.6e} ms, dt_max={dt_max:.6e} ms")

            # 动态色标范围: vmin 取非零数据的 5th 百分位以避免被噪声拉低
            vmax = float(nimp_plot.max())
            if np.any(nimp_plot > 0):
                vmin = float(np.percentile(nimp_plot[nimp_plot > 0], 5))
            else:
                vmin = 0.0

            im = ax_depth.pcolormesh(
                time_edges,
                psi_edges,
                nimp_plot.T,
                cmap='inferno',
                vmin=vmin,
                vmax=vmax,
                shading='flat'
            )
            cbar = plt.colorbar(im, ax=ax_depth, label=r'$n_{imp}$ ($10^{20}\,m^{-3}$)', 
                               pad=0.015, fraction=0.046)
            cbar.ax.tick_params(labelsize=11)

            if debug_enabled:
                print('time_edges:' + str(time_edges[0]) + ', ' + str(time_edges[-1]) + ' , ' + str(psi_edges[0]) + ', ' + str(psi_edges[-1]))
                print(f"nimp_plot shape: {nimp_plot.shape}, t_plot length: {len(t_plot)}, psi_bins_nimp length: {len(psi_bins_nimp)}")
                print(f'nimp_plot.T min: {nimp_plot.T.min():.4e}, 5th percentile: {vmin:.4e}, max: {vmax:.4e}')
                print(f'nimp_plot.T raw data xrange and len: {t_plot[0]}-{t_plot[-1]}, {nimp_plot.T.shape[0]}')
                print(f'nimp_plot.T raw data yrange and len: {psi_bins_nimp[0]}-{psi_bins_nimp[-1]}, {nimp_plot.T.shape[1]}')
                print(f"cloud_front_ range: {np.nanmin(cloud_front_psi):.4f} - {np.nanmax(cloud_front_psi):.4f}")
                
        # 叠加云前沿和云峰线（虚线）
        cloud_front = cloud_front_psi[mask_t]
        cloud_peak  = cloud_peak_psi[mask_t]
        ax_depth.plot(t_plot, cloud_front, color='white', linestyle='--', linewidth=1.5, alpha=1, label='Cloud Front')
        ax_depth.plot(t_plot, cloud_peak, color='#ffdd44', linestyle='--', linewidth=1.5, alpha=1, label='Peak Density')
    
    # 叠加质量加权质心线（始终绘制，在热力图上清晰可见）
    ax_depth.plot(x_val*1e3, res_df['psi_n'], color='#1f77b4', linestyle='-', linewidth=2.8, label="Avg Psi_N")
    
    if not has_nimp and 'psi_n_std' in res_df.columns:
        # 无 nimp 数据时回退到 ±1σ 包络
        ax_depth.fill_between(
            x_val*1e3,
            res_df['psi_n'] - res_df['psi_n_std'],
            res_df['psi_n'] + res_df['psi_n_std'],
            color='#1f77b4', alpha=0.3, label="Cloud Thickness"
        )
    
    # 叠加移动 q=-2 面的 Psi_N 位置曲线
    if has_q_evo:
        t_q_ms = (time_at_q_arr - t0) * 1e3
        mask_q = (t_q_ms >= 0) & (t_q_ms <= args.max_time_ms)
        if mask_q.sum() >= 2:
            ax_depth.plot(t_q_ms[mask_q], psi_n_at_q_arr[mask_q],
                         color='#9467bd', linestyle='--', linewidth=2.5,
                         label=f"q={args.target_q} surface")
    elif target_q_psi_initial is not None:
        ax_depth.axhline(target_q_psi_initial, color='#9467bd', linestyle='--', linewidth=2.5,
                        label=f"q={args.target_q} (Psi_N={target_q_psi_initial:.3f})")
        
    ax_depth.set_ylabel("Psi_N")
    ax_depth.set_title("Pellet Penetration and Ablation Cloud Evolution")
    ax_depth.grid(True, alpha=0.3)
    ax_depth.legend(loc='upper right')
    ax_depth.set_xlim(left=0, right=args.max_time_ms)
    ax_depth.set_ylim(bottom=0, top=1.0)
    ax_ind += 1

    # 第二张子图: 移动 q=target_q 面的 Te 演化
    ax_te = axes[ax_ind]
    if not te_df.empty:
        te_x = (te_df['time_si'] - t0) * 1e3
        ax_te.plot(te_x, te_df['Te'], color='#ff7f0e', linestyle='-', linewidth=2.5,
                   label=f"Te at evolving q={args.target_q} surface")
        ax_te.axhline(10.0, color='darkorange', linestyle=':', linewidth=2, label="10 eV Limit")
    ax_te.set_ylabel("Te (eV)")
    ax_te.set_yscale('log')
    ax_te.grid(True, which="both", ls="--", alpha=0.6)
    ax_te.legend()
    ax_te.set_xlim(left=0, right=args.max_time_ms)
    ax_ind += 1

    # 第三张子图 (可选): 消融率
    if has_abl:
        ax_abl = axes[ax_ind]
        t_abl_ms = (abl_df['time_si'] - t0)*1e3
        
        # 定义需要打印到 Legend 的积分区间大小（优先填完全消融前的质量，其次TQ前的总量，最后全图）
        display_abl = abl_to_r0 if abl_to_r0 is not None else abl_to_tq
        if display_abl is None:
            display_abl = abl_total
        label_text = f"Ablation Rate (Inj 1)\nIntegrated: {display_abl:.2e}" if display_abl else "Ablation Rate (Inj 1)"
        
        ax_abl.plot(t_abl_ms, abl_df['abl_rate_inj1'], color='#d62728', linestyle='-', label=label_text)
        
        # 确定积分区间的上限时间（优先取 local TQ 即 Te<10eV 时刻）
        t_tq_for_fill = t_local_tq if t_local_tq is not None else t_10_percent_thermal
        t_end_val = t_r0 if t_r0 is not None else t_tq_for_fill
        if t_end_val is not None:
            mask = abl_df['time_si'] <= t_end_val
            t_fill_ms = t_abl_ms[mask]
            y_fill = abl_df.loc[mask, 'abl_rate_inj1']
        else:
            t_fill_ms = t_abl_ms
            y_fill = abl_df['abl_rate_inj1']

        # 用半透明红色填绘画积分范围
        ax_abl.fill_between(t_fill_ms, 0, y_fill, color='#d62728', alpha=0.2)
        
        ax_abl.set_ylabel("Ablation Rate ( /s)")
        ax_abl.grid(True)
        ax_abl.legend()
        ax_abl.set_xlim(left=0, right=args.max_time_ms)
        ax_abl.set_ylim(bottom=0)
        ax_ind += 1

    # 第四张子图 (可选): 热能演化
    if has_egy:
        ax_egy = axes[ax_ind]
        t_egy = (energy_df['time_si'] - t0)*1e3
        ax_egy.plot(t_egy, energy_df['thermal'], linestyle='-', linewidth=3, label="Thermal energy")
        ax_egy.plot(t_egy, energy_df['e_thermal'], linestyle='--', linewidth=3, label="Electron thermal energy")
        ax_egy.plot(t_egy, energy_df['i_thermal'], linestyle=':', linewidth=3, label="Ion thermal energy")
        ax_egy.set_ylabel("Energy (J)")
        ax_egy.grid(True)
        ax_egy.legend()
        ax_egy.set_xlim(left=0, right=args.max_time_ms)
        ax_egy.set_ylim(bottom=0)
        ax_ind += 1

    axes[-1].set_xlabel(x_label)
    plt.tight_layout()
    
    clean_out = args.out_plot.replace('.png', '_clean.png')
    if clean_out == args.out_plot:
        clean_out = args.out_plot + "_clean.png"
        
    plt.savefig(clean_out, dpi=300)
    print(f"Clean time-penetration diagram saved to {clean_out}")

    # --- 后期补充所有事件线和文字，用于第二张图 ---
    draw_event_lines(ax_depth)
    draw_event_lines(ax_te)
    if has_abl: draw_event_lines(ax_abl)
    if has_egy: draw_event_lines(ax_egy)
    
    if t_r0 is not None and t_r0 >= t0:
        tx = (t_r0 - t0)*1e3
        ax_depth.text(
            tx,
            0.98,
            f' Complete Ablation ({tx:.2f}ms) ',
            transform=ax_depth.get_xaxis_transform(),
            color='black',
            va='top',
            ha='right',
            rotation=90,
            fontweight='bold',
            fontsize=14,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, pad=2.0)
        )
        
    t_tq_plot = t_local_tq if (t_local_tq is not None) else t_10_percent_thermal
    if t_tq_plot is not None and t_tq_plot >= t0:
        tx = (t_tq_plot - t0)*1e3
        label_text = 'TQ ' if (t_local_tq is not None) else ' TQ (90%) '
        ax_te.text(tx+0.045, 0.5e2, f'{label_text}({tx:.2f}ms) ', color='darkorange', va='bottom', ha='right', rotation=90, fontweight='bold', fontsize=14)

    # if t_psi95 is not None and t_psi95 >= t0:
    #     tx = (t_psi95 - t0)*1e3
    #     ax_depth.text(tx, 0.8, f' q95 ({tx:.2f}ms) ', color='black', va='top', ha='right', rotation=90, fontweight='bold', fontsize=14)
    if t_target_q is not None and t_target_q >= t0:
        tx = (t_target_q - t0)*1e3
        ax_depth.text(tx, 0.6, f'q={args.target_q} ({tx:.2f}ms) ', color="#f8f0ff", va='top', ha='right', rotation=90, fontweight='bold', fontsize=14)

    annotated_out = args.out_plot.replace('.png', '_annotated.png')
    if annotated_out == args.out_plot:
        annotated_out = args.out_plot + "_annotated.png"
        
    plt.savefig(annotated_out, dpi=300)
    print(f"Annotated time-penetration diagram saved to {annotated_out}")

if __name__ == "__main__":
    main()
