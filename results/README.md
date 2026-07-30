# results/

Terminal experiment outputs, moved out of the repo root (2026-07). These are
stdout captures, run logs, and plots — **not** inputs. Scripts still run from the
repo root and read/write their inputs and intermediates there; only these
after-the-fact artifacts live here.

| Folder | Contents |
|---|---|
| `breadth/` | Tissue-breadth classification results (`breadth_*`, `breadth_gcmatched_*`) |
| `offset/` | Offset / distance scans (`chr1_offset_*`, `offset_v*`, `offset_fixed_*`) |
| `onoff/` | Per-tissue on/off classification (`onoff_<Tissue>_*`) |
| `go_3way/` | GO 3-way classification (`go_3way_*`) |
| `go_pairwise/` | GO pairwise classification, incl. GC-matched (`pairwise_*`, `pairwise_gcm_*`) |
| `threeway/` | 3-way prom/genic/random (`threeway_*`, `threeway_genomewide_*`) |
| `fourway/` | 4-way classification (`fourway_*`, `fourway_trimmed_*`) |
| `uorf/` | uORF experiments (`uorf_*`, warmup / decoy variants) |
| `prom_genic/` | Promoter-vs-genic, GC controls (`prom_vs_genic_*`, `genicB*`, `gccontrol*`) |
| `logs/` | Run logs (`*.log`) |
| `plots/` | Figures (`*.png`) — plotting scripts now write here |

## What stayed at the repo root (and why)

- **Scripts** — `*.py`, `*.sh`.
- **Gene-list inputs** — `*_genes.txt`, `onoff_*_high/_low.txt`,
  `housekeeping_genes.txt`, `tissuespecific_genes.txt` (read by build/extract scripts).
- **Intermediates that scripts read back** — `*.bed`, `*.fasta`, `*.fa`,
  `ism_importance_*.npy` (read by `scan_motifs_hivlo.py`).
- **Pre-existing output dirs referenced by scripts** — `window_output/`,
  `offset_output/`, `jaspar_motifs/`, `onoff_liver_model_seed*/`.
- **Raw data** — `*.gz` / `.gct.gz` / `.gaf.gz`.
