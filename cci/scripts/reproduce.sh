#!/usr/bin/env bash
# Recompute into a new result namespace. Does not republish the six frozen reports.
set -euo pipefail
run=${1:?Usage: bash scripts/reproduce.sh results/new_run_name}
if [[ -e "$run" ]]; then echo "Output already exists: $run" >&2; exit 1; fi
mkdir -p "$run"
export LD_LIBRARY_PATH=/opt/anaconda3/lib
export PYTHONPATH="$PWD/.deps"
export NUMBA_CACHE_DIR="$PWD/.cache/numba"
export MPLCONFIGDIR="$PWD/.cache/mpl"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python scripts/fetch_inputs.py > "$run/input_hashes.log"
Rscript scripts/export_cellchat_resources.R > "$run/resources.log" 2>&1
Rscript scripts/export_visium_annotation.R > "$run/annotation.log" 2>&1
python scripts/run_kang.py "$run/kang" > "$run/kang.log" 2>&1
python scripts/null_strata.py "$run/kang" > "$run/null_strata.log"
python scripts/run_ic_sanity.py "$run/kang" "$run/ic_sanity" > "$run/ic_sanity.log" 2>&1
python scripts/run_scaling.py "$run/scaling" > "$run/scaling.log" 2>&1
python scripts/run_scaling.py "$run/pooled_full" --full-only > "$run/pooled_full.log" 2>&1
python scripts/run_spatial.py "$run/spatial" > "$run/spatial.log" 2>&1
python scripts/finish_spatial.py "$run/spatial" > "$run/typed_kernel.log"
python scripts/prepare_baselines.py > "$run/prepare_baselines.log"
python scripts/prepare_cellchat_input.py > "$run/prepare_cellchat.log"
for method in CellPhoneDB FastCCC LIANA; do
 /usr/bin/time -v -o "$run/native_${method}.time.txt" python scripts/run_native.py "$method" "$run/native_${method}" > "$run/native_${method}.log" 2>&1
done
/usr/bin/time -v -o "$run/native_CellChat.time.txt" Rscript scripts/run_cellchat.R "$run/native_CellChat" > "$run/native_CellChat.log" 2>&1
python scripts/analyze_native.py "$run/native_analysis" "$run/native_CellPhoneDB" "$run/native_FastCCC" "$run/native_LIANA" "$run/native_CellChat" > "$run/native_analysis.log" 2>&1
python scripts/prepare_spatial_cellchat.py "$run/spatial/spots.tsv" > "$run/prepare_spatial_native.log"
/usr/bin/time -v -o "$run/native_SpatialCellChat.time.txt" Rscript scripts/run_cellchat_spatial.R "$run/native_SpatialCellChat" > "$run/native_SpatialCellChat.log" 2>&1
python scripts/index_data.py > "$run/index.log"
printf 'Completed new numerical outputs: %s\n' "$run"
