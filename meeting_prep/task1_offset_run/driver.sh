#!/bin/bash
# Task 1: fine-scale chr1 offset scan.
# 500bp windows, stepped every 100bp across 0-2000bp from the TSS.
# Promoter-vs-GC-matched-random (chr1_offset.py: GC +/-0.10, repeats<0.25, N<0.05).
# Data generated ONCE per offset with a FIXED seed (same gene set at every
# offset); classifier trained with seeds 0 1 42 (MPS is not deterministic, so
# we report the spread). classify_ap.py reports AUROC + AUPRC.
set -u
cd "$(dirname "$0")"
PY=/opt/homebrew/Caskroom/miniforge/base/envs/dnabert2/bin/python
DATAGEN_SEED=42            # fixed -> identical genes + negatives across offsets

echo "START $(date)"
for offset in $(seq 0 100 2000); do
    echo "===== OFFSET $offset  $(date +%H:%M:%S) ====="
    # regenerate promoters.fasta / random_sequences.fa for THIS offset (local to run dir)
    $PY chr1_offset.py --offset "$offset" --seed "$DATAGEN_SEED" \
        > "gen_offset_${offset}.log" 2>&1
    for seed in 0 1 42; do
        $PY classify_ap.py --seed "$seed" \
            --promoters promoters.fasta --randoms random_sequences.fa \
            > "chr1_fine_offset_${offset}_seed_${seed}.txt" 2>&1
        line=$(grep FINAL_RESULT "chr1_fine_offset_${offset}_seed_${seed}.txt")
        echo "  offset=$offset seed=$seed -> ${line:-NO_RESULT}"
    done
done
echo "DONE $(date)"
