# Cross-species promoter classification (human vs mouse, extensible)

Question: can DNABERT-2 distinguish the promoter of a human gene from the
promoter of that gene's 1:1 mouse ortholog, from sequence alone — and does the
signal survive the controls, or is it just transposons / GC?

## What makes this defensible (the four controls)

1. **Matched gene set.** Every human promoter is paired with the promoter of its
   1:1 ortholog. Same functions on both sides → separability = species sequence
   divergence, not gene-composition differences.
2. **Repeat control (the cross-species trap).** Human and mouse carry *different*
   transposons (Alu vs B1/B2). Without masking, the model "detects species" by
   spotting lineage-specific repeats, not promoter biology. We use the
   soft-masked genome and drop pairs where either promoter exceeds a repeat
   fraction. **This is the single most important control here** — more so than in
   the within-human work.
3. **GC control.** Per-pair drop at |Δgc| > 0.10 (lab standard). We also print the
   raw pool GC gap so you can answer Joel's "is GC matching even necessary?"
4. **Pair-level / chromosome holdout.** Each ortholog pair → train/val/test by the
   **human** gene's chromosome (chr8+chr9 test, chr10 val, rest train). Both
   species' promoters of a pair share a split, so a gene's human promoter can't
   be in train while its mouse ortholog is in test (ortholog leakage).

## Downloads (mouse) — assembly MUST match

BioMart current orthologs are **GRCm39**, so use GENCODE mouse **GRCm39** (do NOT
grab an old GRCm38/mm10 FASTA — IDs would silently map to wrong coordinates).

```
# GENCODE mouse M39 (GRCm39) — genome + annotation (same assembly guaranteed)
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/GRCm39.primary_assembly.genome.fa.gz
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/gencode.vM39.basic.annotation.gtf.gz
# gunzip both, then: samtools faidx GRCm39.primary_assembly.genome.fa
```
Put the paths into `species_config.py` (mouse entry). Human already points at your
GRCh38 + GENCODE v47.

## Run order

```
conda activate dnabert2            # python 3.9, MPS, transformers pinned 4.36.2
pip install pyfaidx requests       # if not already present

python get_orthologs.py --target mouse --out orthologs_mouse.tsv
bash run_experiment.sh             # builds raw / repeat / repeat_gc, trains, evals
```

`run_experiment.sh` wraps the canonical DNABERT-2 finetune call; swap in your
`classify.py` invocation inside `train_and_eval()` if you prefer the FASTA path
(`human_promoters.fasta` + `mouse_promoters.fasta` are emitted too, with distinct
names to dodge the classify.py clobbering footgun).

## Reading the result (report all three)

| Level        | If AUROC is…                          | Interpretation |
|--------------|----------------------------------------|----------------|
| raw          | high                                   | expected — GC + transposons |
| repeat       | drops a lot                            | signal was mostly lineage-specific repeats |
| repeat_gc    | still clearly > 0.5                    | **real conserved-but-divergent promoter signal** |
| repeat_gc    | collapses to ~0.5                      | no promoter-specific cross-species signal after controls |

Either repeat_gc outcome is a clean, reportable result. Report AUROC **and**
AUPRC with across-seed spread (MPS variance is meaningful).

## Extending along the evolutionary timeline (Joel's follow-up)

Add a species to `TARGETS` in `species_config.py` (only `homolog_prefix` +
GRCm39-equivalent genome/gtf), then `TARGET=chimp bash run_experiment.sh`.
Expectation: closer species (chimp) → harder to separate → lower repeat_gc AUROC,
giving an AUROC-vs-divergence curve. For a shared gene set traceable across many
species, pull multi-species orthologs and intersect the human gene IDs across all
target TSVs before building.

## Files
- `species_config.py` — species registry (paths, BioMart prefixes)
- `get_orthologs.py` — high-confidence 1:1 orthologs from Ensembl BioMart
- `build_dataset.py` — extract, control, split, emit (fastas + CSVs + manifest)
- `eval_metrics.py` — AUROC + AUPRC + bootstrap CI
- `run_experiment.sh` — confound decomposition × seeds
- `manifest.csv` (per run) — full provenance; audit any pair
