# Overnight work — 2026-07-18

Everything below happened in the worktree
`.claude/worktrees/awesome-mendeleev-f66987`. **Nothing in `~/promoter_project`
was edited, run, or touched.** The reference inputs were symlinked in read-only.

---

## 1. Symlinks set up

Symlinked these read-only reference files from `~/promoter_project` into the
worktree (all are `.gitignore`d, so they don't get committed):

- `chr22.fa` (+ `chr22.fa.fai`)
- `gencode.v47.basic.annotation.gtf.gz`
- `gencode.v47.2wayconspseudos.gtf.gz`  ← the Yale-UCSC 2-way consensus set
- `uorfs_capped.fasta`

`chr22_gene_intervals.tsv` was **not** symlinked: it's already committed to git
in the worktree and is byte-identical to the copy in `~/promoter_project`, so I
left the tracked version in place.

I deliberately did **not** symlink `promoters.bed`, `random_regions.bed`,
`promoters.fasta`, `random_sequences.fa`, or any `.txt`/`.log` — those are your
running job's live outputs.

Note on your background job: by the time I looked, the fixed-offset sweep had
already **finished** (last write to `fixed_run.log` was 01:23, no
`caffeinate`/`classify` process was alive, and the log ends cleanly after
`offset_fixed_10000`). So my genome-wide run below was not actually competing
with it. Final numbers from your sweep, for the record:

| offset (bp upstream) | AUC (seed 42) |
|---|---|
| 0    | 0.772 |
| 2000 | 0.854 |
| 4000 | 0.817 |
| 6000 | 0.811 |
| 8000 | 0.775 |
| 10000| 0.734 |

---

## 2. Code review

### 2.1 Cross-script inconsistencies (the thing you specifically asked about)

**GC-matching tolerance is defined four different ways.** Same "GC-matched
random control" concept, four different windows:

| script | tolerance | notes |
|---|---|---|
| `gcmatch_window_v3.py` | **±0.10** per-sequence | this is the one your live pipeline uses |
| `gcmatch_to_promoters.py` | ±0.03 around promoter mean | `--tolerance` default |
| `gcmatch_3way.py` | ±0.03 around combined-class mean | hardcoded `target ± 0.03` |
| `gc_control_test.py` | absolute `0.40–0.60` band | not promoter-relative at all |
| `gcmatch_picker.py` / `gcmatch_window.py` | absolute `0.40–0.60` | older versions |

±0.10 (v3) is a **loose** match — a "GC-matched" negative can sit 10 percentage
points off the promoter it's matched to. Given the AUCs live around 0.77–0.85,
some of that signal is plausibly residual GC/compositional difference rather
than regulatory grammar. Worth deciding on one tolerance and using it everywhere;
I'd argue ±0.03–0.05. (For the genome-wide work I matched the value your live
pipeline actually uses, ±0.10, so the comparison is apples-to-apples — but see §3.)

**Repeat-mask threshold `0.25` is hardcoded in 8+ places** (`randompicker.py`,
`gcmatch_*.py`, `genic_windows.py`, `build_orf_decoys.py`,
`build_4way_classes.py`, `build_uorf_negatives.py` implicitly, …). It's at least
*consistent*, but it's a copy-pasted literal. If you ever revisit it you'll have
to change it in every file. Should be a shared constant / CLI arg.

**GC helper has two different definitions.** Most scripts use
`(G+C)/len(seq)` — which silently counts `N` in the denominator. Only
`gcmatch_window_v3.py` uses `(G+C)/(A+T+G+C)`, which correctly ignores `N`. On
windows with assembly gaps these disagree. Minor on chr22, but it's a real
inconsistency in what "GC content" means across the codebase.

### 2.2 Bugs / edge cases

- **`gcmatch_window_v3.py:71`** — the closing print reads
  `window if positives else args.window`, but `window` is the loop variable and
  only exists if `positives` is non-empty; if `promoters.bed` were empty this
  line would `NameError`. Cosmetic (the loop above would already have done
  nothing), but the guard is backwards.

- **Off-by-`max(0,...)` boundary handling is uneven.** `parse_gtf.py`,
  `parse_gtf_window.py`, `parse_gtf_offset.py` clamp `prom_start` with `max(0,…)`
  but never check `prom_end` against the chromosome length — a minus-strand gene
  near the end of the chromosome can produce a window that runs off the end.
  `parse_gtf_fixed.py` *does* check both ends (`win_end > CHR_LEN`) but hardcodes
  `CHR_LEN = 50818468` (chr22 only) — that constant is silently wrong for any
  other chromosome, which matters the moment you go genome-wide.

- **`overlaps_any_except` name-collision assumption.** `gene_overlap.py`'s
  filter skips windows overlapping a gene *other than the one you're centred on*,
  keyed on `gene_name`. Two distinct genes that share a `gene_name` (happens with
  some readthrough / PAR / duplicated annotations) would be treated as "self" and
  a real overlap would be let through. Rare, but it's a correctness gap in the
  isolation filter.

- **`build_uorf_negatives.py` vs `build_orf_decoys.py` model two different
  negatives** ("random sub-interval of non-uORF 5'UTR" vs "real decoy ORF with
  start+stop in non-uORF 5'UTR"). Both are wired into different 4-way runs. Not a
  bug, but the two negative definitions are easy to confuse — the 4-way results
  will mean different things depending on which negatives file was used. I
  verified the decoy coordinates are at least internally valid: all 711 entries
  in `uorf_decoys_chr22.bed` do start with a start-codon and end with a stop-codon
  in the annotated strand (spot-check passed, 711 OK / 0 bad).

- **`load_fasta` variants everywhere.** At least 4 near-identical copies (some in
  `classify*.py`, some open the file without a `with`, leaking the handle). Fine
  functionally; worth factoring into one shared module given how many classifiers
  now exist.

- **`classify.py` default is CPU-fragile on a shared machine.** All the
  `classify*.py` scripts hard-select `mps` if available. Two of them running at
  once share Apple-Silicon unified memory and will fight. For the genome-wide run
  I forced **CPU** specifically so I would not contend with anything on MPS.

### 2.3 Reproducibility nits

- `randompicker.py`/`shufflepicker.py` reseed **inside** the loop
  (`random.seed(i)` per sequence in `shufflepicker.py`) — intentional there (per-
  sequence shuffle) but easy to mistake for a global seed.
- The experiment shell scripts hardcode a python at
  `…/miniforge/base/envs/dnabert2/bin/python`. That env currently has
  `transformers 4.57.6` installed — which, see §4, **breaks DNABERT-2's
  trust_remote_code loader**. Your runs work because they were launched before the
  upgrade and the interpreter is still holding the old module in memory. A fresh
  launch of those scripts today would fail at model load. Flagging because it will
  bite the next time you start a sweep.

---

## 3. Task I picked: scaling the Yale-UCSC 2-way-consensus pseudogene
##    analysis genome-wide

### Why this one

You gave two candidate next-steps. The GO / GTEx multi-way classification needs
data files that **aren't downloaded** and, more importantly, needs a labelling
scheme you and Joel haven't pinned down yet (which GO slim? which GTEx tissues?
expression threshold?) — not a good fit for unattended overnight work where I
shouldn't be inventing the experimental design.

The genome-wide pseudogene scale-up, by contrast, was **fully specified by work
already in the repo**: the chr22 3-way pipeline
(`parse_promoters_3way.py` → `gcmatch_3way.py` → `classify_3way.py`) already
exists, and the genome-wide Yale-UCSC 2-way consensus GTF
(`gencode.v47.2wayconspseudos.gtf.gz`) was already sitting in your project
directory — you clearly downloaded it on 07-15 for exactly this. So this is
"finish the thing that's already started," which is what you asked for.

### What I built

New scripts (all committed):

- `gene_overlap_genomewide.py` — per-chromosome gene-interval loader with a
  pickle cache, so the 3.1 GB GTF is only parsed once across the sweep.
- `parse_promoters_genomewide.py` — extracts functional (protein-coding) and
  pseudogene (Yale-UCSC 2-way consensus) promoter windows across all standard
  chromosomes, **with the gene-overlap and boundary fixes from the code review
  folded in** (see below).
- `gcmatch_genomewide.py` — GC-matched, non-genic, repeat-filtered **and now
  N-gap-filtered** random controls, sampled genome-wide.
- `extract_fasta_genomewide.py` — pulls sequences via `pyfaidx` (indexed random
  access) instead of loading a 3.1 GB genome string into RAM.
- `classify_3way_genomewide.py` — the DNABERT-2 3-way classifier, CPU-pinned.
- `aggregate_genomewide.py` — mean±std across seeds.
- `run_genomewide_seeds.sh` — sequential multi-seed driver.

Fixes carried in from the review, so the genome-wide version doesn't repeat the
chr22 pipeline's shortcuts:

1. **Overlap filtering** — chr22 `parse_promoters_3way.py` / `gcmatch_3way.py`
   did *not* exclude promoter windows falling inside another gene body. The
   genome-wide version does (`overlaps_any_except`), across both the
   protein-coding and pseudogene interval sets.
2. **Assembly-gap filter** — the first negative-control pass pulled one window
   that was 36% `N` (a centromeric/telomeric gap). Added a `--max_n 0.05` reject.
   Genome-wide this matters; on chr22 it rarely triggered.
3. **Randomised selection** — genes/pseudogenes are shuffled before taking the
   first N, instead of taking whatever comes first in file (chromosome) order, so
   the sample isn't biased toward chr1.
4. **No hardcoded chromosome length** — boundary checks use the actual per-
   chromosome length from the FASTA index.

### Data built (committed as BED; FASTAs are gitignored)

- `functional_promoters_genomewide.bed` — 300 protein-coding promoters
- `pseudo_promoters_genomewide.bed` — 300 Yale-UCSC 2-way consensus pseudogene
  promoters
- `random_genomewide.bed` — 599 GC-matched non-genic controls (1 of 600 couldn't
  be placed within the attempt cap; acceptable)

Sanity: all windows 1000 bp, 0% N after the gap filter, mean GC functional 0.505
/ pseudo 0.409 / random 0.434.

### Result (seed 42, genome-wide)

```
Overall accuracy: 0.613
AUC functional vs rest: 0.814
AUC pseudogene vs rest: 0.722
AUC random  vs rest:    0.756
```

Compared to the chr22-only 3-way (5-seed mean): functional-vs-rest ≈ 0.87,
pseudogene-vs-rest ≈ 0.85, random-vs-rest ≈ 0.92.

**So the AUCs drop when we go genome-wide** — most visibly for pseudogene and
random. Two things to keep in mind before reading too much into that:

- **Different pseudogene population.** The chr22 3-way defined pseudogenes by
  `gene_type` string-matching against the *basic annotation* (~90 pseudogenes on
  chr22). The genome-wide run uses the **Yale-UCSC 2-way consensus** set, which is
  a different, more stringently-defined population. Lower pseudogene AUC may partly
  be "these are harder/cleaner pseudogenes," not just "more chromosomes."
- **Stricter negatives.** The genome-wide controls now exclude gene-body overlaps
  and assembly gaps that the chr22 pipeline didn't, which removes some of the easy
  compositional giveaways — expected to *lower* AUC and make the number more
  honest.

This is a **single seed** (each CPU run is ~2 h). I launched seeds 0/1/123/999 in
sequence (`run_genomewide_seeds.sh`, PID logged in `genomewide_seeds_driver.log`);
they should trickle in overnight. Run `python3 aggregate_genomewide.py` in the
morning to get the mean±std across whatever finished.

---

## 4. Environment issue you should know about (didn't touch it)

The shared miniforge base env has been upgraded to `transformers 5.12.1`, and the
`dnabert2` conda env to `4.57.6`. **Both are too new for DNABERT-2's
`trust_remote_code` model files** — a fresh model load now fails two different
ways (`BertConfig has no attribute pad_token_id`, then a `config_class`
registration mismatch). Your currently-running processes are fine because they
loaded the module before the upgrade, but the next fresh launch of any
`classify*.py` will break.

Rather than modify any shared environment (you asked me to keep things stable), I
created an **isolated venv** in the worktree, `venv_genomewide/`
(`--system-site-packages`, `transformers==4.36.2`), and added a small monkeypatch
in `classify_3way_genomewide.py` to relax the over-strict `register()` check. That
venv is gitignored. When you want to re-run the older chr22 classifiers, the
cleanest fix is to pin `transformers==4.36.2` (or similar) in the `dnabert2` env —
happy to do that with you, but I didn't want to change a shared env unprompted.

---

## 5. What I'd do next (for you + Joel)

1. Let the remaining genome-wide seeds finish, then `aggregate_genomewide.py` for
   error bars; add a bar plot mirroring `plot_3way.py`.
2. Decide on **one** GC tolerance and repeat threshold and thread them through as
   shared constants — the ±0.03 vs ±0.10 split makes cross-experiment comparison
   shaky.
3. If genome-wide pseudogene-vs-rest holds around ~0.72 across seeds, that's the
   interesting headline: DNABERT-2 separates functional from pseudogene promoters
   noticeably *less* well on the stringent 2-way-consensus set than the chr22
   number suggested.
