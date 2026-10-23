#!/bin/bash
# 
# Automated SPI Diagnostics Workflow script
# This script orchestrates convert2diagno_spi, jorek2_diagno_spi,
# jorek2_postproc, track_largest_pellets.py and auto_spi_shard_surfaces.py
# to fully automate SPI pellet injection tracking and diagnostics.
#
export OMP_NUM_THREADS=1
ulimit -s unlimited
set -e

# Option / Parameter parsing
#   --plot-only : skip Steps 1-3 (fragment extraction, ablation history and
#                 postproc diagnostics) and only produce plots from the
#                 existing outputs in auto_diagno_spi/
PLOT_ONLY=false
POSITIONAL_ARGS=()
for arg in "$@"; do
    case "$arg" in
        --plot-only)
            PLOT_ONLY=true
            ;;
        *)
            POSITIONAL_ARGS+=("$arg")
            ;;
    esac
done
set -- "${POSITIONAL_ARGS[@]}"

DIAGNO_BIN=${1:-./jorek2_diagno_spi}
POSTPROC_BIN=${2:-./jorek2_postproc}
MACRO_VARS=${3:-./macroscopic_vars.dat}
NAMELIST=${4:-jorek_namelist_ntor11}

# Number of parallel workers (adjust based on available cores)
N_PARALLEL=30

# postproc 在 Step 3 里使用的表达式列表（rectangle / average 共用）。
# 重要：这也决定了 exprs_Rmin*.h5 里 values 的**变量顺序** ——
# jorek2_postproc 把 h5 的变量名数组写成了 character(len=1)，变量名会退化成
# 首字母（P/T/n/n），所以 Step 5 依赖 --expr_order 按位置复原变量名。
RECT_EXPRS=${RECT_EXPRS:-"Psi_N T_e nimp ne"}

# All outputs go under this directory
OUTPUT_DIR="auto_diagno_spi"
POSTPROC_DIR="$OUTPUT_DIR/postproc"
CSV_FOLDER="$OUTPUT_DIR/csv_files"

mkdir -p "$POSTPROC_DIR" "$CSV_FOLDER"

# Use a unique temp directory per run to avoid conflicts with parallel diagnostics
UNIQUE_TEMP="/tmp/jorek_convert_$$_$(date +%s)"
mkdir -p "$UNIQUE_TEMP"
# Ensure cleanup of temp dir on script exit
cleanup_temp() { rm -rf "$UNIQUE_TEMP"; }
trap cleanup_temp EXIT

# Cleanup worker temp files
rm -f "$OUTPUT_DIR"/postproc_step3_worker*.in "$OUTPUT_DIR"/postproc_step3_worker*.log 2>/dev/null

echo "============================================================"
echo "       JOREK SPI Diagnostic Suite — Automated Workflow        "
echo "   (fragment centroid, ablation cloud, q-slice Te, energies)  "
echo "============================================================"
echo "Output directory: $OUTPUT_DIR"
if [ "$PLOT_ONLY" = true ]; then
    echo "Mode: plot-only (reusing existing diagnostics, Steps 1-3 skipped)"
fi

# Wrapper to run postproc (single-core, stdin-driven Fortran program)
run_postproc() {
    local in_file=$1
    local log_file=$2
    "$POSTPROC_BIN" < "$in_file" > "$log_file" 2>&1
}

# Shared output paths used by Step 4 (also referenced in plot-only mode)
QPROF_MERGED="$POSTPROC_DIR/qprofile_merged.dat"
MERGED_OUT="$POSTPROC_DIR/exprs_averaged_merged.dat"

if [ "$PLOT_ONLY" = false ]; then

# SPI start timestep: minimum step from part*.h5 names
START_STEP_RAW=$(ls part*.h5 2>/dev/null | sed -n 's/^.*part\([0-9]\+\)\.h5.*$/\1/p' | sort -n | head -1)
if [ -z "$START_STEP_RAW" ]; then
    echo "Error: No part*.h5 files found in current directory."
    exit 1
fi
# Remove leading zeros to treat as base-10 integer
START_STEP=$((10#$START_STEP_RAW))
echo "[INFO] Detected SPI start step: $START_STEP"

# ==============================================================================
# Step 1: Extract 3D fragment positions across all time steps (→ fragments_position_*.csv)
# ==============================================================================
echo "[STEP 1/5] Extracting fragment (R,Z,φ) trajectories from all part*.h5 …"
python3 util/convert2diagno_spi.py "$DIAGNO_BIN" "$NAMELIST" "${START_STEP}-99999" -t "$CSV_FOLDER" -j ${N_PARALLEL} --temp-dir "$UNIQUE_TEMP"

# ==============================================================================
# Step 2: Ablation-rate history from the final restart file
# ==============================================================================
echo "[STEP 2/5] Computing ablation-rate history (abl_history.dat) …"
LAST_STEP_RAW=$(ls part*.h5 2>/dev/null | sed -n 's/^.*part\([0-9]\+\)\.h5.*$/\1/p' | sort -n | tail -1)
LAST_STEP=$((10#$LAST_STEP_RAW))
echo "[INFO]   Last available step: $LAST_STEP"

PART_FILE=$(printf "part%05d.h5" $LAST_STEP)
JOREK_FILE=$(printf "jorek%06d.h5" $LAST_STEP)

if [ ! -f "$JOREK_FILE" ] || [ ! -f "$PART_FILE" ]; then
    echo "Error: Restart files for step $LAST_STEP ($JOREK_FILE, $PART_FILE) are missing."
    exit 1
fi

cp "$JOREK_FILE" jorek_restart.h5
cp "$PART_FILE" part_restart.h5

# Execute diagno_spi taking the namelist as input
"$DIAGNO_BIN" < "$NAMELIST" > "$OUTPUT_DIR/diagno_spi_step2.log" || true
if [ -f "abl_history.dat" ]; then
    cp abl_history.dat "$OUTPUT_DIR/abl_history.dat"
    echo "[INFO]   abl_history.dat → $OUTPUT_DIR/"
else
    echo "[WARN]   abl_history.dat was not produced; ablation integrals will be skipped."
fi

# Compute total steps (used by Step 3)
TOTAL_STEPS=$((LAST_STEP - START_STEP + 1))
STEPS_PER_WORKER=$(( (TOTAL_STEPS + N_PARALLEL - 1) / N_PARALLEL ))

# ==============================================================================
# Step 3: Unified parallel diagnostics — rectangle + qprofile + average
#         All three commands run in a single for-loop per worker,
#         sharing one restart load per time step (≈30–40% less wall time).
# ==============================================================================
echo "[STEP 3/5] Unified diagnostics (rectangle + qprofile + average) for all ${TOTAL_STEPS} steps …"
BOUNDS=$(python3 -c "
import glob, sys, pandas as pd, re, os
csvs = glob.glob('$CSV_FOLDER/*.csv')
csvs = [f for f in csvs if re.match('^fragments_position_[0-9]+\\\\.csv$', os.path.basename(f))]
if not csvs: sys.exit(1)
dfs = []
for f in csvs:
    try:
        df_tmp = pd.read_csv(f, delim_whitespace=True, header=None, engine='python')
        if not df_tmp.empty and len(df_tmp.columns) >= 3:
            dfs.append(df_tmp)
    except:
        pass
if not dfs: sys.exit(1)
df = pd.concat(dfs)
rmin, rmax = df.iloc[:, 1].min(), df.iloc[:, 1].max()
zmin, zmax = df.iloc[:, 2].min(), df.iloc[:, 2].max()
print(f'{rmin-0.1:.4f} {rmax+0.1:.4f} {zmin-0.1:.4f} {zmax+0.1:.4f}')
")

if [ $? -eq 0 ] && [ ! -z "$BOUNDS" ]; then
    read RMIN RMAX ZMIN ZMAX <<< "$BOUNDS"
    echo "[INFO]   Fragment (R,Z) envelope: R in [$RMIN, $RMAX]  Z in [$ZMIN, $ZMAX]"
else
    echo "[ERROR]  Failed to compute fragment (R,Z) boundaries from CSV files."
    exit 1
fi

# ---- Clean up all output files from previous runs ----
rm -f "$POSTPROC_DIR"/exprs_Rmin*.h5 2>/dev/null
rm -f "$POSTPROC_DIR"/qprofile_s*.dat 2>/dev/null
rm -f "$POSTPROC_DIR"/exprs_averaged_s*.dat 2>/dev/null
rm -f "$QPROF_MERGED" "$MERGED_OUT" 2>/dev/null
rm -f "$OUTPUT_DIR"/postproc_step3_worker*.in "$OUTPUT_DIR"/postproc_step3_worker*.log 2>/dev/null

echo "[INFO]   ${TOTAL_STEPS} steps ÷ ${N_PARALLEL} workers ≈ ${STEPS_PER_WORKER} steps/worker"

# ---- Launch unified parallel workers ----
unified_pids=()
for ((w=0; w<N_PARALLEL; w++)); do
    w_start=$((START_STEP + w * STEPS_PER_WORKER))
    w_end=$((w_start + STEPS_PER_WORKER - 1))
    if [ $w_end -gt $LAST_STEP ]; then
        w_end=$LAST_STEP
    fi
    if [ $w_start -gt $LAST_STEP ]; then
        break
    fi

    worker_in="$OUTPUT_DIR/postproc_step3_worker${w}.in"
    worker_log="$OUTPUT_DIR/postproc_step3_worker${w}.log"

    cat <<EOF > "$worker_in"
namelist $NAMELIST
set_postproc_dir $POSTPROC_DIR
si-units
expressions $RECT_EXPRS
mark_coords 1
set surfaces 120
for step $w_start to $w_end do
  rectangle $RMIN $RMAX 200 $ZMIN $ZMAX 200 0.
  qprofile
  average
done
EOF

    echo "[INFO]   Unified worker ${w}: steps ${w_start} → ${w_end}"
    "$POSTPROC_BIN" < "$worker_in" > "$worker_log" 2>&1 &
    unified_pids+=($!)
done

echo "[INFO]   Waiting for ${#unified_pids[@]} unified workers …"
wait "${unified_pids[@]}"
echo "[INFO]   All unified workers finished."

# ---- Verify rectangle HDF5 files ----
RECT_COUNT=$(ls "$POSTPROC_DIR"/exprs_Rmin*.h5 2>/dev/null | wc -l)
if [ "$RECT_COUNT" -gt 0 ]; then
    echo "[INFO]   ${RECT_COUNT} rectangle HDF5 files exported → $POSTPROC_DIR/"
else
    echo "[ERROR]  jorek2_postproc did not produce any exprs_Rmin*.h5 files in $POSTPROC_DIR/"
    exit 1
fi

# ---- Merge qprofile output files ----
first_file=true
for dat_file in "$POSTPROC_DIR"/qprofile_s*.dat; do
    [ -f "$dat_file" ] || continue
    if $first_file; then
        cp "$dat_file" "$QPROF_MERGED"
        first_file=false
    else
        awk 'BEGIN{found=0} /^# time step #/{found=1} found{print}' "$dat_file" >> "$QPROF_MERGED"
    fi
done

if [ -f "$QPROF_MERGED" ]; then
    cp "$QPROF_MERGED" "$OUTPUT_DIR/qprofile_merged.dat"
    echo "[INFO]   Merged q-profile → $OUTPUT_DIR/qprofile_merged.dat  ($(wc -l < "$QPROF_MERGED") lines)"
else
    echo "[WARN]   No qprofile output files produced; target-q tracking will be skipped."
fi

# ---- Merge average output files ----
first_file=true
for dat_file in "$POSTPROC_DIR"/exprs_averaged_s*.dat; do
    [ -f "$dat_file" ] || continue
    if $first_file; then
        cp "$dat_file" "$MERGED_OUT"
        first_file=false
    else
        awk 'BEGIN{found=0} /^# time step #/{found=1} found{print}' "$dat_file" >> "$MERGED_OUT"
    fi
done

if [ -f "$MERGED_OUT" ]; then
    echo "[INFO]   Merged profile → $MERGED_OUT  ($(wc -l < "$MERGED_OUT") lines)"
else
    echo "[WARN]   No average output files produced; heatmap and Te-evolution will be skipped."
fi

else   # PLOT_ONLY=true: reuse existing outputs, no Steps 1-3

if [ ! -d "$CSV_FOLDER" ] || [ ! -d "$POSTPROC_DIR" ]; then
    echo "Error: existing outputs ($CSV_FOLDER, $POSTPROC_DIR) not found."
    echo "       Run the full workflow once before using --plot-only."
    exit 1
fi
echo "[INFO] Reusing existing outputs from $CSV_FOLDER and $POSTPROC_DIR."

fi   # end of PLOT_ONLY check



# ==============================================================================
# Step 4: Python synthesis — centroid tracking, heatmap, event detection, PDF
# ==============================================================================
echo "[STEP 4/5] Running track_largest_pellets.py (centroid + cloud heatmap + events) …"
TE_PROFILE_ARG=""
if [ -f "$MERGED_OUT" ]; then
    TE_PROFILE_ARG="--te_profile $MERGED_OUT"
    echo "[INFO]   Feeding Te/nimp profile: $MERGED_OUT"
else
    echo "[INFO]   No merged profile available; Python will attempt auto-discovery."
fi
QPROF_ARG=""
if [ -f "$QPROF_MERGED" ]; then
    QPROF_ARG="--q_profile $QPROF_MERGED"
    echo "[INFO]   Feeding evolving q-profile: $QPROF_MERGED"
else
    echo "[INFO]   No merged q-profile available; target-q line will be omitted."
fi
BOUNDS_ARG=""
if [ -n "${RMIN:-}" ] && [ -n "${RMAX:-}" ] && [ -n "${ZMIN:-}" ] && [ -n "${ZMAX:-}" ]; then
    BOUNDS_ARG="--bounds $RMIN $RMAX $ZMIN $ZMAX"
elif [ "$PLOT_ONLY" = true ]; then
    echo "[INFO]   No stored bounds; Python will infer them from the fragment CSVs."
fi
python3 util/track_largest_pellets.py \
    --csv_dir "$CSV_FOLDER" \
    --rectangle_dir "$POSTPROC_DIR" \
    --abl_file "$OUTPUT_DIR/abl_history.dat" \
    --live_data "$MACRO_VARS" \
    --max_time_ms 0.5 \
    $BOUNDS_ARG \
    $QPROF_ARG \
    $TE_PROFILE_ARG \
    --out_plot "$OUTPUT_DIR/spi_diagnostic_overview.png" \
    --out_csv "$OUTPUT_DIR/pellet_trajectory_summary.csv"

# ==============================================================================
# Step 5: 3D shard surfaces（火箭效应）
#   一张合并图：x=碎片排序轴（默认按半径置零/烧蚀完成先后，见 SURF_XAXIS）,
#               y=自消融开始时间[ms], z=|ΔV_R|,
#               颜色 = 碎片见到的量, 透明度 = r/r0（碎片还剩多少）
#   默认还额外出"碎片加速度"合并图：z=|a_R|=|dV_R/dt|（见 SURF_ACCEL）
#   * 时间轴同样以 macroscopic_vars.dat 的 @sqrt_mu0_rho0 归一化，
#     并以 abl_history.dat 定位的"消融开始"为零点
#   * 环境变量：
#       SURF_COLOR=te|tn|psin|none   颜色场（默认 te）
#       SURF_LAYOUT=merged|split|both 合并图/分开两张/都要（默认 merged）
#       SURF_RENDER=surface|lines|both 曲面/逐碎片折线/叠加（默认 surface）
#       SURF_XAXIS=index|r0|psi0|death 碎片排序（默认 death=按半径置零/烧蚀完成的先后：
#                                     窗口内越早烧蚀完越靠前，窗口末仍存活的按
#                                     末时刻半径升序排在其后；index=恢复原始编号）
#       SURF_ACCEL=on|off            碎片加速度 a_R=dV_R/dt 合并图（默认 on）
#       SKIP_SHARD_SURFACES=1        跳过本步
# ==============================================================================
if [ "${SKIP_SHARD_SURFACES:-0}" != "1" ]; then
    SURF_COLOR=${SURF_COLOR:-te}
    SURF_LAYOUT=${SURF_LAYOUT:-merged}
    SURF_RENDER=${SURF_RENDER:-surface}
    SURF_XAXIS=${SURF_XAXIS:-death}
    SURF_ACCEL=${SURF_ACCEL:-on}
    echo "[STEP 5/5] Running auto_spi_shard_surfaces.py (color=$SURF_COLOR layout=$SURF_LAYOUT render=$SURF_RENDER x_axis=$SURF_XAXIS accel=$SURF_ACCEL) …"
    mkdir -p "$OUTPUT_DIR/surfaces"
    # 删除旧缓存，避免复用到上一次运行的数组（缓存仅供手动重画时通过 --cache_npz 复用）
    rm -f "$OUTPUT_DIR/surfaces/shard_surfaces_cache.npz"
    python3 util/auto_spi_shard_surfaces.py \
        --csv_dir "$CSV_FOLDER" \
        --rectangle_dir "$POSTPROC_DIR" \
        --abl_file "$OUTPUT_DIR/abl_history.dat" \
        --live_data "$MACRO_VARS" \
        --color "$SURF_COLOR" \
        --layout "$SURF_LAYOUT" \
        --render "$SURF_RENDER" \
        --x_axis "$SURF_XAXIS" \
        --accel "$SURF_ACCEL" \
        --expr_order "$RECT_EXPRS" \
        $BOUNDS_ARG \
        --out_dir "$OUTPUT_DIR/surfaces" \
        --out_prefix shard_surfaces \
        --cache_npz "$OUTPUT_DIR/surfaces/shard_surfaces_cache.npz" \
        || echo "[WARN]   auto_spi_shard_surfaces.py 失败或数据不足，已跳过 3D 曲面输出。"
else
    echo "[STEP 5/5] SKIP_SHARD_SURFACES=1, 跳过 3D 曲面输出。"
fi

echo "============================================================"
echo "[DONE]  SPI diagnostic suite finished successfully."
echo ""
echo "  Output directory        : $OUTPUT_DIR/"
echo "  Diagnostic plot         : $OUTPUT_DIR/spi_diagnostic_overview.png"
echo "  Trajectory CSV          : $OUTPUT_DIR/pellet_trajectory_summary.csv"
echo "  Diagnostic summary      : $OUTPUT_DIR/diagnostic_summary.txt"
echo "============================================================"
