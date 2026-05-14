#!/bin/bash
# 
# Automated SPI Diagnostics Workflow script
# This script orchestrates convert2diagno_spi, jorek2_diagno_spi, 
# jorek2_postproc, and track_largest_pellets.py to fully automate
# SPI pellet injection tracking and diagnostics.
#
export OMP_NUM_THREADS=1
ulimit -s unlimited
set -e

# Option / Parameter parsing could be added here, but for now we take defaults.
DIAGNO_BIN=${1:-./jorek2_diagno_spi}
POSTPROC_BIN=${2:-./jorek2_postproc}
MACRO_VARS=${3:-./macroscopic_vars.dat}
NAMELIST=${4:-input_EHL_t2}

# Number of parallel workers (adjust based on available cores)
N_PARALLEL=16

# All outputs go under this directory
OUTPUT_DIR="auto_diagno_spi"
POSTPROC_DIR="$OUTPUT_DIR/postproc"
CSV_FOLDER="$OUTPUT_DIR/csv_files"

mkdir -p "$POSTPROC_DIR" "$CSV_FOLDER"

# Cleanup worker temp files
rm -f "$OUTPUT_DIR"/postproc_step45_worker*.in "$OUTPUT_DIR"/postproc_step45_worker*.log 2>/dev/null

echo "================================================="
echo "      JOREK Automated SPI Diagnostics Workflow     "
echo "================================================="

# Wrapper to run postproc (single-core, stdin-driven Fortran program)
run_postproc() {
    local in_file=$1
    local log_file=$2
    "$POSTPROC_BIN" < "$in_file" > "$log_file" 2>&1
}

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
# Step 1: Create CSV folder mapping 3D positions over time
# ==============================================================================
echo "[STEP 1] Generating CSV files using convert2diagno_spi.py..."
python3 util/convert2diagno_spi.py "$DIAGNO_BIN" "$NAMELIST" "${START_STEP}-99999" -t "$CSV_FOLDER" -j ${N_PARALLEL}

# ==============================================================================
# Step 2: Extract Ablation history 
# ==============================================================================
echo "[STEP 2] Generating ablation history (abl_history.dat)..."
LAST_STEP_RAW=$(ls part*.h5 2>/dev/null | sed -n 's/^.*part\([0-9]\+\)\.h5.*$/\1/p' | sort -n | tail -1)
LAST_STEP=$((10#$LAST_STEP_RAW))
echo "[INFO] Last step found: $LAST_STEP"

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
    echo "[INFO] abl_history.dat copied to $OUTPUT_DIR/"
else
    echo "Warning: abl_history.dat was not produced!"
fi

# ==============================================================================
# Step 3: Produce average rectangle Psi_N mappings around fragments
# ==============================================================================
echo "[STEP 3] Generating R,Z to Psi_N mapping..."
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
# Add safety margin
print(f'{rmin-0.1:.4f} {rmax+0.1:.4f} {zmin-0.1:.4f} {zmax+0.1:.4f}')
")

if [ $? -eq 0 ] && [ ! -z "$BOUNDS" ]; then
    read RMIN RMAX ZMIN ZMAX <<< "$BOUNDS"
    echo "[INFO] Fragment trajectory bounding box: R[$RMIN, $RMAX], Z[$ZMIN, $ZMAX]"
    
    # 清理之前可能产生的同名残留文件以避免 HDF5 写入冲突
    rm -f "$POSTPROC_DIR"/exprs_*.h5 2>/dev/null
    
    cat <<EOF > "$OUTPUT_DIR/postproc_step3.in"
namelist $NAMELIST
set_postproc_dir $POSTPROC_DIR
expressions Psi_N
for step $START_STEP do
  rectangle $RMIN $RMAX 200 $ZMIN $ZMAX 200 0.
done
EOF
    run_postproc "$OUTPUT_DIR/postproc_step3.in" "$OUTPUT_DIR/postproc_step3.log"

    RECC_FILE=$(ls "$POSTPROC_DIR"/exprs_*.h5 2>/dev/null | tail -1)
    if [ -n "$RECC_FILE" ]; then
        cp "$RECC_FILE" "$OUTPUT_DIR/rectangle_psi.h5"
        echo "[INFO] Successfully copied $RECC_FILE to $OUTPUT_DIR/rectangle_psi.h5"
    else
        echo "Error: Failed to generate exprs_*.h5 file in $POSTPROC_DIR/"
        exit 1
    fi
else
    echo "Error: Failed to compute trajectory boundaries from CSVs."
    exit 1
fi

# ==============================================================================
# Step 4: Generate precise q-profile for TQ thresholds
# ==============================================================================
echo "[STEP 4] Generating q-profile..."
cat <<EOF > "$OUTPUT_DIR/postproc_step4.in"
namelist $NAMELIST
set_postproc_dir $POSTPROC_DIR
for step $START_STEP do
  qprofile
done
EOF
run_postproc "$OUTPUT_DIR/postproc_step4.in" "$OUTPUT_DIR/postproc_step4.log"

QPROF_FILE=$(ls "$POSTPROC_DIR"/qprofile_*.dat 2>/dev/null | tail -1)
if [ -n "$QPROF_FILE" ]; then
    cp "$QPROF_FILE" "$OUTPUT_DIR/qprofile.dat"
    echo "[INFO] Successfully copied $QPROF_FILE to $OUTPUT_DIR/qprofile.dat"
else
    echo "Warning: Failed to generate qprofile_*.dat file in $POSTPROC_DIR/"
fi

# ==============================================================================
# Step 4.5: Generate averaged Te profile for local TQ detection
# ==============================================================================
echo "[STEP 4.5] Generating averaged Te profiles (parallelized)..."



TOTAL_STEPS=$((LAST_STEP - START_STEP + 1))
STEPS_PER_WORKER=$(( (TOTAL_STEPS + N_PARALLEL - 1) / N_PARALLEL ))

echo "[INFO] Total steps: $TOTAL_STEPS, workers: $N_PARALLEL, steps/worker: ~$STEPS_PER_WORKER"

# Clean up previous temp files
rm -f "$POSTPROC_DIR"/exprs_averaged_*.dat 2>/dev/null

# Launch parallel postproc workers
worker_pids=()
for ((w=0; w<N_PARALLEL; w++)); do
    w_start=$((START_STEP + w * STEPS_PER_WORKER))
    w_end=$((w_start + STEPS_PER_WORKER - 1))
    if [ $w_end -gt $LAST_STEP ]; then
        w_end=$LAST_STEP
    fi
    if [ $w_start -gt $LAST_STEP ]; then
        break
    fi

    worker_in="$OUTPUT_DIR/postproc_step45_worker${w}.in"
    worker_log="$OUTPUT_DIR/postproc_step45_worker${w}.log"

    cat <<EOF > "$worker_in"
namelist $NAMELIST
set_postproc_dir $POSTPROC_DIR
si-units
expressions Psi_N T_e
mark_coords 1
set surfaces 120
for step $w_start to $w_end do
  average
done
EOF

    echo "[INFO] Launching worker $w: steps $w_start -> $w_end"
    "$POSTPROC_BIN" < "$worker_in" > "$worker_log" 2>&1 &
    worker_pids+=($!)
done

# Wait for all workers to finish
echo "[INFO] Waiting for ${#worker_pids[@]} postproc workers to complete..."
wait "${worker_pids[@]}"
echo "[INFO] All workers finished."

# Merge output files: collect all .dat files, keep step delimiters, deduplicate column headers
MERGED_OUT="$POSTPROC_DIR/exprs_averaged_merged.dat"
first_file=true
for dat_file in "$POSTPROC_DIR"/exprs_averaged_s*.dat; do
    [ -f "$dat_file" ] || continue
    if $first_file; then
        # Keep entire first file (header + all step blocks)
        cp "$dat_file" "$MERGED_OUT"
        first_file=false
    else
        # Skip column header block, keep everything from first "# time step" onward
        awk 'BEGIN{found=0} /^# time step #/{found=1} found{print}' "$dat_file" >> "$MERGED_OUT"
    fi
done

if [ -f "$MERGED_OUT" ]; then
    echo "[INFO] Merged averaged profiles into $MERGED_OUT ($(wc -l < "$MERGED_OUT") lines)"
else
    echo "Warning: No averaged profile output files found."
fi



# ==============================================================================
# Step 5: Final Python Diagnostic Plots
# ==============================================================================
echo "[STEP 5] Generating visualization and data synthesis..."
TE_PROFILE_ARG=""
if [ -f "$MERGED_OUT" ]; then
    TE_PROFILE_ARG="--te_profile $MERGED_OUT"
    echo "[INFO] Passing Te profile: $MERGED_OUT"
else
    echo "[INFO] No merged Te profile found, Python will auto-search."
fi
python3 util/track_largest_pellets.py \
    --csv_dir "$CSV_FOLDER" \
    --psi_file "$OUTPUT_DIR/rectangle_psi.h5" \
    --abl_file "$OUTPUT_DIR/abl_history.dat" \
    --q_profile "$OUTPUT_DIR/qprofile.dat" \
    --live_data "$MACRO_VARS" \
    --max_time_ms 2.0 \
    --bounds $RMIN $RMAX $ZMIN $ZMAX \
    $TE_PROFILE_ARG \
    --out_plot "$OUTPUT_DIR/spi_diagnostic_overview.png" \
    --out_csv "$OUTPUT_DIR/pellet_trajectory_summary.csv"

echo "================================================="
echo "[SUCCESS] SPI Diagnostics Complete."
echo "Please check $OUTPUT_DIR/spi_diagnostic_overview.png."
