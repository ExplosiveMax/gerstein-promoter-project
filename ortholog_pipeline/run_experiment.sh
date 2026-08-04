#!/usr/bin/env bash
set -euo pipefail

TARGET=${TARGET:-mouse}
N_PAIRS=${N_PAIRS:-2000}
SEEDS=${SEEDS:-"0 1 42"}
ORTHO="orthologs_${TARGET}.tsv"
ROOT="runs/${TARGET}"
CLASSIFY=${CLASSIFY:-"$HOME/promoter_project/classify.py"}

# --- env guard: this script calls bare `python`, so it MUST resolve to the
# dnabert2 env (python 3.9 + transformers 4.36.2). If you forgot to
# `conda activate dnabert2`, bare python falls through to base 3.12, whose newer
# transformers breaks DNABERT-2 with "BertConfig has no attribute pad_token_id"
# — and you only find out after the datasets are built. Fail fast instead.
ver=$(python -c "import transformers; print(transformers.__version__)" 2>/dev/null || true)
echo "env check: python=$(command -v python)  transformers=${ver:-<none>}"
if [ "$ver" != "4.36.2" ]; then
    echo "FATAL: transformers is '${ver:-<none>}', expected 4.36.2." >&2
    echo "       Run 'conda activate dnabert2' (or put its bin first on PATH) and retry." >&2
    exit 3
fi

[ -f "$ORTHO" ] || python get_orthologs.py --target "$TARGET" --out "$ORTHO"

python build_dataset.py --orthologs "$ORTHO" --target "$TARGET" \
    --out-dir "$ROOT/raw"       --n-pairs "$N_PAIRS" --no-repeat-filter --no-gc-match
python build_dataset.py --orthologs "$ORTHO" --target "$TARGET" \
    --out-dir "$ROOT/repeat"    --n-pairs "$N_PAIRS" --no-gc-match
python build_dataset.py --orthologs "$ORTHO" --target "$TARGET" \
    --out-dir "$ROOT/repeat_gc" --n-pairs "$N_PAIRS"

for level in raw repeat repeat_gc; do
  for seed in $SEEDS; do
    echo "=== $TARGET / $level / seed $seed ==="
    cp "$ROOT/$level/human_promoters.fasta"     promoters.fasta
    cp "$ROOT/$level/${TARGET}_promoters.fasta" random_sequences.fa
    python -u "$CLASSIFY" --seed "$seed" | tee "$ROOT/$level/seed_${seed}.txt"
  done
done

echo "--- AUC summary ---"
grep -H FINAL_RESULT "$ROOT"/*/seed_*.txt || true
