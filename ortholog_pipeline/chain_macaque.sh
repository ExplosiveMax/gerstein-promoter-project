#!/usr/bin/env bash
# Unattended chain: wait for the Task-1 offset grid to finish (so it does NOT
# share the GPU), then run the full macaque cross-species experiment.
#   1. TARGET=macaque N_PAIRS=1000 run_experiment.sh  (raw/repeat/repeat_gc x seeds)
#   2. classify_holdout.py on runs/macaque/repeat_gc for seeds 0 1 42
#   3. summarize_macaque.py -> summary table
# NOTE: run_experiment.sh cp's fastas to ./promoters.fasta / ./random_sequences.fa
# INSIDE ortholog_pipeline/ (its own scratch) — the repo TOP-LEVEL
# ~/promoter_project/promoters.fasta is never touched.
set -uo pipefail
cd "$(dirname "$0")"
# CRITICAL: run_experiment.sh calls bare `python`. Put the dnabert2 env FIRST on
# PATH so bare `python` == python 3.9 + transformers 4.36.2 (NOT base 3.12, which
# breaks DNABERT-2 with "BertConfig has no attribute pad_token_id").
export PATH="/opt/homebrew/Caskroom/miniforge/base/envs/dnabert2/bin:$PATH"
PY=/opt/homebrew/Caskroom/miniforge/base/envs/dnabert2/bin/python
DRIVER_LOG="$HOME/promoter_project/meeting_prep/task1_offset_run/driver_run.log"

echo "CHAIN_START $(date)"
# fail fast if the env is wrong, rather than burning the night on a broken run
ver=$(python -c "import transformers; print(transformers.__version__)" 2>/dev/null)
echo "python=$(which python)  transformers=$ver"
if [ "$ver" != "4.36.2" ]; then
    echo "FATAL: transformers is '$ver', expected 4.36.2 — aborting." ; exit 3
fi

# ---- gate: offset grid must be finished and no classify_ap.py (offset) running ----
# DONE marker prevents a false positive during the ~1min data-gen gap between
# offsets (when momentarily no classify_ap.py is running but the grid isn't done).
echo "Waiting for offset grid to finish (poll every 60s)..."
while true; do
    if ! pgrep -f classify_ap.py >/dev/null 2>&1; then
        if grep -q '^DONE' "$DRIVER_LOG" 2>/dev/null \
           || ! pgrep -f 'bash driver.sh' >/dev/null 2>&1; then
            break
        fi
    fi
    sleep 60
done
echo "OFFSET_GRID_CLEAR $(date) — GPU free, starting macaque"

# ---- 1. full decomposition (AUROC + AUPRC via classify_ap.py drop-in) ----
export CLASSIFY="$HOME/promoter_project/ortholog_pipeline/classify_ap.py"
TARGET=macaque N_PAIRS=1000 SEEDS="0 1 42" bash run_experiment.sh
rc_decomp=$?
echo "DECOMP_DONE rc=$rc_decomp $(date)"

# ---- 2. chromosome-holdout eval on repeat_gc, seeds 0 1 42 ----
for seed in 0 1 42; do
    echo "=== macaque holdout repeat_gc seed $seed  $(date +%H:%M:%S) ==="
    $PY classify_holdout.py --data-dir runs/macaque/repeat_gc --seed "$seed" \
        > "runs/macaque/repeat_gc/holdout_seed_${seed}.txt" 2>&1
    grep -h HOLDOUT_RESULT "runs/macaque/repeat_gc/holdout_seed_${seed}.txt" || echo "  (no HOLDOUT_RESULT — check log)"
done
echo "HOLDOUT_DONE $(date)"

# ---- 3. summary table ----
echo
$PY summarize_macaque.py
echo "CHAIN_ALL_DONE $(date)"
