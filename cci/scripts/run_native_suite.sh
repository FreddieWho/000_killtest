#!/usr/bin/env bash
set -u
export NUMBA_CACHE_DIR="$PWD/.cache/numba"
export MPLCONFIGDIR="$PWD/.cache/mpl"
export LD_LIBRARY_PATH=/opt/anaconda3/lib
export PYTHONPATH="$PWD/.deps"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for method in CellPhoneDB FastCCC LIANA; do
  /usr/bin/time -v -o "results/native_${method}_01.time.txt" python scripts/run_native.py "$method" "results/native_${method}_01" > "results/native_${method}_01.log" 2>&1
  rc=$?
  printf '%s\t%s\n' "$method" "$rc" >> results/native_suite_exit.tsv
done
/usr/bin/time -v -o results/native_CellChat_01.time.txt Rscript scripts/run_cellchat.R results/native_CellChat_01 > results/native_CellChat_01.log 2>&1
rc=$?
printf 'CellChat\t%s\n' "$rc" >> results/native_suite_exit.tsv
