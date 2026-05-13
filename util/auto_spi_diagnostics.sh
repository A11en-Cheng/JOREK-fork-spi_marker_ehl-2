#!/bin/bash
# 
# Automated SPI Diagnostics Workflow script
# This script orchestrates convert2diagno_spi, jorek2_diagno_spi, 
# jorek2_postproc, and track_largest_pellets.py to fully automate
# SPI pellet injection tracking and diagnostics.
#

set -e

# Option / Parameter parsing could be added here, but for now we take defaults.
DIAGNO_BIN=${1:-./jorek2_diagno_spi}
POSTPROC_BIN=${2:-./jorek2_postproc}
MACRO_VARS=${3:-./macroscopic_vars.dat}
NAMELIST=${4:-input_EHL_t2}
CSV_FOLDER="csv_file_automatic"

echo "================================================="
echo "      JOREK Automated SPI Diagnostics Workflow     "
echo "================================================="

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
mkdir -p "$CSV_FOLDER"
python3 util/convert2diagno_spi.py "$DIAGNO_BIN" "$NAMELIST" "${START_STEP}-99999" -t "$CSV_FOLDER" -j16

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
"$DIAGNO_BIN" < "$NAMELIST" > diagno_spi_step2.log || true
if [ ! -f "abl_history.dat" ]; then
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
    rm -f postproc/exprs_*.h5 2>/dev/null
    
    cat <<EOF | "$POSTPROC_BIN" > postproc_step3.log
namelist $NAMELIST
expressions Psi_N
for step $START_STEP do
  rectangle $RMIN $RMAX 200 $ZMIN $ZMAX 200 0.
done
EOF

    RECC_FILE=$(ls postproc/exprs_*.h5 2>/dev/null | tail -1)
    if [ -n "$RECC_FILE" ]; then
        cp "$RECC_FILE" rectangle_psi.h5
        echo "[INFO] Successfully copied $RECC_FILE to rectangle_psi.h5"
    else
        echo "Error: Failed to generate exprs_*.h5 file in postproc/"
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
cat <<EOF | "$POSTPROC_BIN" > postproc_step4.log
namelist $NAMELIST
for step $START_STEP do
  qprofile
done
EOF

QPROF_FILE=$(ls postproc/qprofile_*.dat 2>/dev/null | tail -1)
if [ -n "$QPROF_FILE" ]; then
    cp "$QPROF_FILE" qprofile.dat
    echo "[INFO] Successfully copied $QPROF_FILE to qprofile.dat"
else
    echo "Warning: Failed to generate qprofile_*.dat file in postproc/"
fi

# ==============================================================================
# Step 5: Final Python Diagnostic Plots
# ==============================================================================
echo "[STEP 5] Generating visualization and data synthesis..."
python3 util/track_largest_pellets.py \
    --csv_dir "$CSV_FOLDER" \
    --psi_file rectangle_psi.h5 \
    --top_n 100 \
    --abl_file abl_history.dat \
    --q_profile qprofile.dat \
    --live_data "$MACRO_VARS" \
    --max_time_ms 2.0 \
    --bounds $RMIN $RMAX $ZMIN $ZMAX \
    --out_plot spi_diagnostic_overview.png

echo "================================================="
echo "[SUCCESS] SPI Diagnostics Complete."
echo "Please check spi_diagnostic_overview.png."
