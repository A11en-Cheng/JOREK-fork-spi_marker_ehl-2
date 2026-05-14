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
    time_scaling = 1.0
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
    if not q_profile_file or not os.path.exists(q_profile_file):
        return None
    try:
        col_names = ["psi_n", "q"]
        df = pd.read_csv(q_profile_file, delim_whitespace=True, comment='#', names=col_names)
        
        # 移除重复或者无用的数值
        df_sorted = df.drop_duplicates(subset=["q"]).sort_values(by="q")
        if target_q < df_sorted["q"].min() or target_q > df_sorted["q"].max():
            print(f"Target q={target_q} is out of bounds in {q_profile_file}.")
            return None
        
        f_psi_n = interp1d(df_sorted["q"], df_sorted["psi_n"])
        return float(f_psi_n(target_q))
    except Exception as e:
        print(f"Failed to read/interpolate {q_profile_file}: {e}")
        return None

def get_te_evolution_at_psi(psi_target, te_profile_files, time_map, time_scaling, t0, limit_ev=10.0):
    """提取目标面局部Te演变，并寻找局部热猝灭（Te < 10eV）时间
    
    Parameters
    ----------
    psi_target : float or None
        目标 Psi_N 值
    te_profile_files : list of str
        由外部传入的 Te profile 文件路径列表（已解析好的，不再内部搜索）
    time_map : dict
        step → jorek time 映射
    time_scaling : float
        时间缩放因子
    t0 : float
        起始时间
    limit_ev : float
        Te 低于此阈值 (eV) 视为局部 TQ
        
    Returns
    -------
    t_tq : float or None
    t_tq_step : int or None
    te_df : pd.DataFrame
    """
    if psi_target is None:
        return None, None, pd.DataFrame()
    if not te_profile_files:
        return None, None, pd.DataFrame()
    
    te_history = []
    t_tq = None
    t_tq_step = None
    
    for f in te_profile_files:
        current_step = None
        data_rows = []
        
        def process_block(step, rows):
            nonlocal t_tq, t_tq_step
            try:
                df = pd.DataFrame(rows, columns=['psi_n', 'Te'])
                df = df.dropna().sort_values(by="psi_n")
                if not df.empty and psi_target >= df["psi_n"].min() and psi_target <= df["psi_n"].max():
                    f_te = interp1d(df["psi_n"], df["Te"])
                    te_val = float(f_te(psi_target))
                    real_time = time_map.get(step, step) * time_scaling
                    if real_time >= t0:
                        te_history.append((step, real_time, te_val))
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
                            # Parse step like "# time step #002750"
                            step_str = line.split("#")[-1].strip()
                            current_step = int(step_str)
                        except ValueError:
                            pass
                else:
                    parts = line.split()
                    if len(parts) >= 2:
                        data_rows.append([float(parts[0]), float(parts[1])])
                        
        if current_step is not None and data_rows:
            process_block(current_step, data_rows)
            
    te_df = pd.DataFrame(te_history, columns=["step", "time_si", "Te"])
    if not te_df.empty:
        te_df = te_df.sort_values(by="time_si")
    return t_tq, t_tq_step, te_df

def load_averaged_nimp(profile_files, time_map, time_scaling, t0, psi_to_psin_init=None):
    """
    解析 average 命令输出的径向 profile 文件，提取 nimp(Psi_N_initial, t) 二维分布。
    
    期望输入文件包含 3 或 4 列:
      3 列: Psi_N_current, T_e, nimp
      4 列: Psi_N_current, T_e, nimp, Psi
    
    若提供 psi_to_psin_init 映射且文件含 4 列，则用 Psi 将坐标映射到初始 Psi_N。
    
    Parameters
    ----------
    profile_files : list of str
    time_map : dict
    time_scaling : float
    t0 : float
    psi_to_psin_init : callable or None
        Psi → Psi_N_initial 的插值函数
    
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
        n_cols = None  # 每个文件独立检测列数
        
        def commit_block(step, rows, ncols):
            if step is None or not rows:
                return
            try:
                if ncols == 4:
                    cols = ['psi_n_cur', 'Te', 'nimp', 'psi']
                else:
                    cols = ['psi_n_cur', 'Te', 'nimp']
                df = pd.DataFrame(rows, columns=cols)
                df = df.dropna()
                if df.empty or len(df) < 2:
                    return
                if df['nimp'].max() == df['nimp'].min():
                    return
                real_time = time_map.get(step, step) * time_scaling
                if real_time >= t0:
                    if ncols == 4 and psi_to_psin_init is not None:
                        # 用 Psi → Psi_N_initial 映射
                        psi_n_init = psi_to_psin_init(df['psi'].values)
                        # 过滤超出范围的插值结果
                        valid = ~np.isnan(psi_n_init)
                        if valid.sum() < 2:
                            return
                        all_blocks.append((step, real_time, psi_n_init[valid], df['nimp'].values[valid]))
                    else:
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
                            commit_block(current_step, data_rows, n_cols)
                            data_rows = []
                        try:
                            step_str = line.split("#")[-1].strip()
                            current_step = int(step_str)
                        except ValueError:
                            pass
                else:
                    parts = line.split()
                    if n_cols is None:
                        n_cols = len(parts)
                    if len(parts) >= 3:
                        try:
                            data_rows.append([float(p) for p in parts[:n_cols]])
                        except ValueError:
                            pass
        if current_step is not None and data_rows:
            commit_block(current_step, data_rows, n_cols or 3)
    
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
        f_interp = interp1d(psi_arr, nimp_arr, kind='linear', 
                           bounds_error=False, fill_value=0.0)
        nimp_rows.append(f_interp(psi_bins))
        time_rows.append(t_si)
    
    nimp_2d = np.array(nimp_rows)
    time_si = np.array(time_rows)
    
    nimp_2d = np.maximum(nimp_2d, 1e-10)
    
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
    parser.add_argument("--psi_file", type=str, required=True, help="Flux surface dat file exported from jorek2_postproc")
    parser.add_argument("--live_data", type=str, default='macroscopic_vars.dat', help="Path to macroscopic_vars.dat to extract physical time")
    parser.add_argument("--abl_file", type=str, default="abl_history.dat", help="Path to abl_history.dat produced by jorek2_diagno_spi")
    parser.add_argument("--start_time", type=float, default=None, help="Optional start time (in SI seconds) to plot from")
    parser.add_argument("--max_time_ms", type=float, default=3.0, help="Maximum plotted time in ms from the start_time (default: 3.0 ms)")
    parser.add_argument("--q_profile", type=str, default=None, help="Path to q_profile dat file")
    parser.add_argument("--target_q", type=float, default=-2.0, help="Target q-value to draw horizontal line in penetration plot")
    parser.add_argument("--out_csv", type=str, default="pellet_average_trajectory.csv", help="Output trajectory file")
    parser.add_argument("--out_plot", type=str, default="penetration_trajectory.png", help="Output visualization plot")
    parser.add_argument("--bounds", type=float, nargs=4, metavar=('RMIN', 'RMAX', 'ZMIN', 'ZMAX'), default=None, help="Bounding box of the rectangle grid if needed")
    parser.add_argument("--te_profile", type=str, default=None, help="Path to averaged Te profile .dat file (from jorek2_postproc average command)")
    parser.add_argument("--nimp_profile", type=str, default=None, help="Path to averaged nimp profile .dat file (defaults to --te_profile if not given)")
    
    args = parser.parse_args()

    # 0. 载入时间映射表
    time_map = {}
    time_scaling = 1.0
    if args.live_data and os.path.exists(args.live_data):
        time_map, time_scaling = load_time_mapping(args.live_data)
        print(f"Loaded {len(time_map)} time steps from {args.live_data} (Scaling factor: {time_scaling:.4e})")
    else:
        if args.live_data:
            print(f"Warning: {args.live_data} not found, will only use steps.")

    # 1. 载入磁场网格数据
    R_data, Z_data, Psi_N_data, Psi_data = load_psi_n_h5(args.psi_file, bounds=args.bounds)
    print(f"Loaded 2D grid data with {len(R_data)} valid points for Psi_N interpolation.")
    
    psi_interpolator = LinearNDInterpolator(list(zip(R_data, Z_data)), Psi_N_data)
    
    # 构建 Psi → Psi_N_initial 映射（用于将 average 输出的 Psi 翻回初始 Psi_N）
    has_psi_mapping = Psi_data is not None
    if has_psi_mapping:
        # 以 Psi 为自变量，Psi_N 为因变量，通过排序构建单调映射
        order = np.argsort(Psi_data)
        psi_sorted = Psi_data[order]
        psi_n_sorted = Psi_N_data[order]
        # 去重（Psi 值可能因网格分辨率有重复）
        unique_mask = np.concatenate([[True], np.diff(psi_sorted) > 0])
        psi_to_psin_init = interp1d(psi_sorted[unique_mask], psi_n_sorted[unique_mask],
                                    kind='linear', bounds_error=False, 
                                    fill_value=(psi_n_sorted[unique_mask][0], psi_n_sorted[unique_mask][-1]))
        psi_axis_init = float(psi_sorted[unique_mask][0])
        psi_bnd_init = float(psi_sorted[unique_mask][-1])
        print(f"Built Psi→Psi_N_initial mapping: Psi_axis={psi_axis_init:.4f}, Psi_bnd={psi_bnd_init:.4f}")
    else:
        psi_to_psin_init = None
        psi_axis_init = psi_bnd_init = None

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
            fragment_R = df_active['spi_R'].values
            fragment_Z = df_active['spi_Z'].values
            fragment_psi = psi_interpolator(fragment_R, fragment_Z)
            
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
        if args.start_time is not None:
            energy_df = energy_df[energy_df['time_si'] >= args.start_time]

    rad_power_df = load_radiated_power(args.live_data, time_scaling)
    if rad_power_df is not None and not rad_power_df.empty:
        if args.start_time is not None:
            rad_power_df = rad_power_df[rad_power_df['time_si'] >= args.start_time]

    has_abl = abl_df is not None and not abl_df.empty
    has_egy = energy_df is not None and not energy_df.empty
    has_rad = rad_power_df is not None and not rad_power_df.empty
    
    target_q_psi = get_target_q_psi_n(args.q_profile, args.target_q)
    if target_q_psi is not None:
        print(f"Interpolated Psi_N for q={args.target_q} is {target_q_psi:.6f}")

    # 计算标记事件时间
    t_r0 = None 
    t_r0_step = None
    # 半径小于等于 1e-10 就视作 0
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

    t_psi95 = None
    t_psi95_step = None
    psi95_idx = res_df[res_df['psi_n'] <= 0.95].first_valid_index()
    if psi95_idx is not None:
        t_psi95 = float(res_df.loc[psi95_idx, 'time_si'])
        t_psi95_step = int(res_df.loc[psi95_idx, 'step'])

    t_target_q = None
    t_target_q_step = None
    if target_q_psi is not None:
        target_q_idx = res_df[res_df['psi_n'] <= target_q_psi].first_valid_index()
        if target_q_idx is not None:
            t_target_q = float(res_df.loc[target_q_idx, 'time_si'])
            t_target_q_step = int(res_df.loc[target_q_idx, 'step'])
            
    _t0_idx = (res_df['time_si'] - t0).abs().idxmin()
    t0_step = int(res_df.loc[_t0_idx, 'step'])
    
    # 查找局部温度 TQ 和获取 Te 演化历史
    # 由外部传入 te_profile 文件路径（支持单独调用），若未传则自动搜索
    if args.te_profile:
        te_files = [args.te_profile]
    else:
        # 自动搜索：优先新目录结构，降级到旧路径
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
    t_local_tq, t_local_tq_step, te_df = get_te_evolution_at_psi(
        target_q_psi, te_files, time_map, time_scaling, t0, limit_ev=10.0)
    
    # 加载 nimp(Psi_N, t) 二维分布用于热力图
    if args.nimp_profile:
        nimp_files = [args.nimp_profile]
    else:
        nimp_files = te_files  # 默认复用 Te 文件（同一 average 输出含 nimp 列）
    nimp_2d, psi_bins_nimp, time_si_nimp = load_averaged_nimp(
        nimp_files, time_map, time_scaling, t0, psi_to_psin_init=psi_to_psin_init)
    if nimp_2d is not None:
        # 减去本底（首步的 nimp profile），只显示 SPI 贡献的杂质增量
        nimp_baseline = nimp_2d[0, :].copy()
        nimp_excess = nimp_2d - nimp_baseline[np.newaxis, :]
        nimp_excess = np.maximum(nimp_excess, 0.0)
        # 计算 Psi_N 坐标下的云扩散前沿（超过阈值的最小 Psi_N）和云峰
        nimp_threshold = 0.05 * nimp_excess.max() if nimp_excess.max() > 0 else 1.0
        cloud_front_psi = np.array([
            psi_bins_nimp[np.where(row > nimp_threshold)[0][0]] 
            if np.any(row > nimp_threshold) else np.nan
            for row in nimp_excess
        ])
        cloud_peak_psi = np.array([
            psi_bins_nimp[np.argmax(row)] if row.max() > 0 else np.nan
            for row in nimp_excess
        ])
        has_nimp = True
        psi_label = "Psi_N_initial" if has_psi_mapping else "Psi_N_current"
        print(f"Loaded nimp 2D distribution: {nimp_2d.shape[0]} time steps × {nimp_2d.shape[1]} {psi_label} bins")
    else:
        has_nimp = False
        print("No valid nimp data found, falling back to ±1σ envelope.")

    print("\n" + "="*50)
    print("                EVENT TIMINGS")
    print("="*50)
    def print_event(name, t_val, step_val):
        if t_val is not None:
            print(f"- {name:<27}: t = {(t_val - t0)*1000:8.4f} ms (relative), Step ~ {step_val}")
        else:
            print(f"- {name:<27}: Not Reached")
            
    print_event("Injection / Start", t0, t0_step)
    print_event("Reached q_95 (Psi_N=0.95)", t_psi95, t_psi95_step)
    if target_q_psi is not None:
        print_event(f"Reached Target q={args.target_q}", t_target_q, t_target_q_step)
        print_event(f"Local TQ (Te<10eV at q={args.target_q})", t_local_tq, t_local_tq_step)
    print_event("Thermal Quench (~90% drop)", t_10_percent_thermal, t_10_percent_step)
    print_event("Complete Ablation", t_r0, t_r0_step)
    print("="*50 + "\n")

    if has_abl:
        print("="*50)
        print("             INTEGRATED ABLATION")
        print("="*50)
        
        def calculate_integrated_ablation(end_time):
            if end_time is None or end_time < t0:
                return None
            mask = (abl_df['time_si'] >= t0) & (abl_df['time_si'] <= end_time)
            df_sub = abl_df[mask]
            if len(df_sub) < 2:
                return 0.0
            return np.trapz(df_sub['abl_rate_inj1'], df_sub['time_si'])

        abl_to_tq = calculate_integrated_ablation(t_10_percent_thermal)
        abl_to_r0 = calculate_integrated_ablation(t_r0)
        abl_total = calculate_integrated_ablation(abl_df['time_si'].max())
        
        if abl_to_tq is not None:
            print(f"- To Thermal Quench     : {abl_to_tq:.4e}")
        if abl_to_r0 is not None:
            print(f"- To Complete Ablation  : {abl_to_r0:.4e}")
        print(f"- Total Ablated (Max)   : {abl_total:.4e}")
        print("="*50 + "\n")

    if has_rad:
        print("="*50)
        print("         INTEGRATED RADIATED ENERGY")
        print("="*50)
        
        def calculate_integrated_radiation(end_time):
            if end_time is None or end_time < t0:
                return None
            mask = (rad_power_df['time_si'] >= t0) & (rad_power_df['time_si'] <= end_time)
            df_sub = rad_power_df[mask]
            if len(df_sub) < 2:
                return 0.0
            return np.trapz(df_sub['radiated_power'], df_sub['time_si'])

        rad_to_tq = calculate_integrated_radiation(t_10_percent_thermal)
        rad_to_r0 = calculate_integrated_radiation(t_r0)
        rad_total = calculate_integrated_radiation(rad_power_df['time_si'].max())
        
        if rad_to_tq is not None:
            print(f"- To Thermal Quench     : {rad_to_tq:.4e} J")
        if rad_to_r0 is not None:
            print(f"- To Complete Ablation  : {rad_to_r0:.4e} J")
        print(f"- Total Radiated (Max)  : {rad_total:.4e} J")
        print("="*50 + "\n")

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
        if t_psi95 is not None and t_psi95 >= t0:
            tx = (t_psi95 - t0)*1e3
            ax.axvline(tx, color='black', linestyle='-.', alpha=0.8, linewidth=2.0)
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
        from matplotlib.colors import LogNorm
        time_ms_nimp = (time_si_nimp - t0) * 1e3
        
        # 裁剪到绘图时间范围
        mask_t = (time_ms_nimp >= 0) & (time_ms_nimp <= args.max_time_ms)
        if mask_t.sum() >= 2:
            t_plot = time_ms_nimp[mask_t]
            nimp_plot = nimp_excess[mask_t, :]
            
            # 动态色标范围: vmin 取最大值的 1e-5 或数据中位数, 避免噪声
            vmax = nimp_plot.max()
            vmin = max(vmax * 1e-5, np.median(nimp_plot[nimp_plot > 0]) if np.any(nimp_plot > 0) else 1e15)
            
            im = ax_depth.imshow(
                nimp_plot.T,
                aspect='auto',
                origin='lower',
                extent=[t_plot[0], t_plot[-1], psi_bins_nimp[0], psi_bins_nimp[-1]],
                cmap='inferno',
                norm=LogNorm(vmin=vmin, vmax=vmax),
                interpolation='bilinear'
            )
            cbar = plt.colorbar(im, ax=ax_depth, label=r'n$_{imp}$ excess (m$^{-3}$)', 
                               pad=0.015, fraction=0.046)
            cbar.ax.tick_params(labelsize=11)
        
        # 叠加云前沿和云峰线（虚线）
        cloud_front = cloud_front_psi[mask_t]
        cloud_peak  = cloud_peak_psi[mask_t]
        ax_depth.plot(t_plot, cloud_front, color='white', linestyle='--', linewidth=1.5, alpha=1, label='Cloud Front (5% max)')
        ax_depth.plot(t_plot, cloud_peak, color='#ffdd44', linestyle='--', linewidth=1.5, alpha=1, label='Cloud Peak Density')
    
    # 叠加质量加权质心线（始终绘制，在热力图上清晰可见）
    ax_depth.plot(x_val*1e3, res_df['psi_n'], color='#1f77b4', linestyle='-', linewidth=2.8, label="Mass-Weighted Avg Psi_N")
    
    if not has_nimp and 'psi_n_std' in res_df.columns:
        # 无 nimp 数据时回退到 ±1σ 包络
        ax_depth.fill_between(
            x_val*1e3,
            res_df['psi_n'] - res_df['psi_n_std'],
            res_df['psi_n'] + res_df['psi_n_std'],
            color='#1f77b4', alpha=0.3, label="±1σ Cloud Thickness"
        )
    
    if target_q_psi is not None:
         ax_depth.axhline(target_q_psi, color='#9467bd', linestyle='--', linewidth=2.5, label=f"q={args.target_q} (Psi_N={target_q_psi:.3f})")
        
    ax_depth.set_ylabel("Psi_N (initial equilibrium)")
    ax_depth.set_title("Pellet Penetration and Ablation Cloud Evolution")
    ax_depth.grid(True, alpha=0.3)
    ax_depth.legend(loc='upper right')
    ax_depth.set_xlim(left=0, right=args.max_time_ms)
    ax_depth.set_ylim(bottom=0, top=1.0)
    ax_ind += 1

    # 第二张子图: q=target_q 面的 Te 演化
    ax_te = axes[ax_ind]
    if not te_df.empty:
        te_x = (te_df['time_si'] - t0) * 1e3
        ax_te.plot(te_x, te_df['Te'], color='#ff7f0e', linestyle='-', linewidth=2.5, label=f"Te at q={args.target_q}")
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
        
        # 确定积分区间的上限时间
        t_end_val = t_r0 if t_r0 is not None else (t_10_percent_thermal if t_10_percent_thermal is not None else None)
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
        y_loc = res_df['psi_n'].max() if not res_df.empty else 1.0
        ax_te.text(tx, 1e2, f' r->0 ({tx:.2f}ms) ', color='dimgrey', va='top', ha='right', rotation=90, fontweight='bold', fontsize=14)
        
    t_tq_plot = t_local_tq if (t_local_tq is not None) else t_10_percent_thermal
    if t_tq_plot is not None and t_tq_plot >= t0:
        tx = (t_tq_plot - t0)*1e3
        label_text = ' Local TQ ' if (t_local_tq is not None) else ' TQ (90%) '
        ax_te.text(tx, 0.5e6, f'{label_text}({tx:.2f}ms) ', color='darkorange', va='bottom', ha='right', rotation=90, fontweight='bold', fontsize=14)

    if t_psi95 is not None and t_psi95 >= t0:
        tx = (t_psi95 - t0)*1e3
        ax_depth.text(tx, 0.8, f' q95 ({tx:.2f}ms) ', color='black', va='top', ha='right', rotation=90, fontweight='bold', fontsize=14)
    if t_target_q is not None and t_target_q >= t0:
        tx = (t_target_q - t0)*1e3
        ax_depth.text(tx, 1.3, f' q={args.target_q} ({tx:.2f}ms) ', color='#9467bd', va='top', ha='right', rotation=90, fontweight='bold', fontsize=14)

    annotated_out = args.out_plot.replace('.png', '_annotated.png')
    if annotated_out == args.out_plot:
        annotated_out = args.out_plot + "_annotated.png"
        
    plt.savefig(annotated_out, dpi=300)
    print(f"Annotated time-penetration diagram saved to {annotated_out}")

if __name__ == "__main__":
    main()
