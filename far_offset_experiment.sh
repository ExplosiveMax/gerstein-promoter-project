#!/bin/bash
PYTHON=/opt/homebrew/Caskroom/miniforge/base/envs/dnabert2/bin/python

for offset in 0 1000 2000 4000 6000 8000 10000; do
    echo "===== Offset: $offset bp upstream ====="
    $PYTHON parse_gtf_far.py --offset $offset --window 500
    bedtools getfasta -fi chr22.fa -bed promoters.bed -fo promoters.fasta
    $PYTHON gcmatch_to_promoters.py
    bedtools getfasta -fi chr22.fa -bed random_regions.bed -fo random_sequences.fa
    for seed in 0 42 999; do
        $PYTHON classify.py --seed $seed | tee far_offset_${offset}_seed_${seed}.txt
    done
done
echo "All done!"