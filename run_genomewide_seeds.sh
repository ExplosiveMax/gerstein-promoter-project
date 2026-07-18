#!/bin/bash
# Run the remaining seeds for the genome-wide 3-way classifier sequentially.
# Uses the isolated venv (DNABERT-2 needs the older transformers pinned there).
# CPU-only, ~2h/seed, so this is an overnight job.
cd "$(dirname "$0")"
source venv_genomewide/bin/activate
for seed in 0 1 123 999; do
    echo "===== genome-wide seed $seed starting $(date) ====="
    PYTHONUNBUFFERED=1 python classify_3way_genomewide.py --seed "$seed" \
        > "threeway_genomewide_seed_${seed}.txt" 2>&1
    echo "===== genome-wide seed $seed done $(date) ====="
done
echo "All genome-wide seeds complete $(date)"
