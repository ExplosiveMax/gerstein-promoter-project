#!/bin/bash
PYTHON=/opt/homebrew/Caskroom/miniforge/base/envs/dnabert2/bin/python

for offset in -1000 -500 0 500 1000; do
    echo "===== Offset: $offset ====="
    $PYTHON parse_gtf_offset.py --offset $offset
    bedtools getfasta -fi chr22.fa -bed promoters.bed -fo promoters.fasta

    # NEW: match randoms to THIS offset's promoter GC
    $PYTHON gcmatch_to_promoters.py

    bedtools getfasta -fi chr22.fa -bed random_regions.bed -fo random_sequences.fa

    for seed in 0 1 42 123 999; do
        $PYTHON classify.py --seed $seed | tee offset_${offset}_seed_${seed}.txt
    done
    echo "===== Done with offset $offset ====="
done
echo "All done!"