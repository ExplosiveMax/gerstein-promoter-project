#!/bin/bash

PYTHON=/opt/homebrew/Caskroom/miniforge/base/envs/dnabert2/bin/python

for window in 256 500 1000 2000; do
    echo "===== Window size: $window ====="
    
    # Generate promoter regions
    $PYTHON parse_gtf_window.py --window $window
    
    # Extract promoter sequences
    bedtools getfasta -fi chr22.fa -bed promoters.bed -fo promoters.fasta
    
    # Generate GC matched random regions
    $PYTHON gcmatch_window.py --window $window
    
    # Extract random sequences
    bedtools getfasta -fi chr22.fa -bed random_regions.bed -fo random_sequences.fa
    
    # Run 5 seeds
    for seed in 0 1 42 123 999; do
        $PYTHON classify.py --seed $seed | tee window_${window}_seed_${seed}.txt
    done
    
    echo "===== Done with window $window ====="
done

echo "All done!"