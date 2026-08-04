"""Fetch high-confidence 1-to-1 orthologs between human and a target species.

Source: Ensembl BioMart (query on the human dataset, homolog attributes for the
target). We keep only orthology_type == 'ortholog_one2one' AND confidence == 1.
That confidence flag is Ensembl's high-confidence call (based on gene-order
conservation + tree consistency) and cuts a lot of dubious mappings.

Output: a 2-column TSV  <human_gene_id>\t<target_gene_id>  (unversioned Ensembl IDs).

Usage:
    python get_orthologs.py --target mouse --out orthologs_mouse.tsv
    # offline fallback if BioMart is flaky:
    python get_orthologs.py --target mouse --out orthologs_mouse.tsv \
        --from-file some_biomart_export.tsv
"""

import argparse
import sys
import time

from species_config import HUMAN, TARGETS

BIOMART_URL = "https://www.ensembl.org/biomart/martservice"

_XML = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="{dataset}" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="{prefix}_homolog_ensembl_gene"/>
    <Attribute name="{prefix}_homolog_orthology_type"/>
    <Attribute name="{prefix}_homolog_orthology_confidence"/>
  </Dataset>
</Query>"""


def _strip_version(gene_id: str) -> str:
    # ENSG00000123456.7 -> ENSG00000123456 ; also drop GENCODE _PAR_Y suffix
    return gene_id.split(".")[0].replace("_PAR_Y", "")


def fetch_biomart(dataset: str, prefix: str, retries: int = 4) -> str:
    import requests  # imported here so offline/--from-file path needs no network deps

    xml = _XML.format(dataset=dataset, prefix=prefix)
    last = None
    for attempt in range(retries):
        try:
            r = requests.post(BIOMART_URL, data={"query": xml}, timeout=180)
            r.raise_for_status()
            if r.text.startswith("Query ERROR"):
                raise RuntimeError(r.text[:400])
            return r.text
        except Exception as e:  # noqa: BLE001
            last = e
            wait = 5 * (attempt + 1)
            print(f"  BioMart attempt {attempt+1} failed ({e}); retry in {wait}s",
                  file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"BioMart failed after {retries} tries: {last}")


def parse_one2one(tsv_text: str) -> list[tuple[str, str]]:
    """Rows are: human_id, target_id, orthology_type, confidence."""
    pairs, seen_h, seen_t = [], set(), set()
    for line in tsv_text.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        h_id, t_id, otype, conf = f[0], f[1], f[2], f[3]
        if not h_id or not t_id:
            continue
        if otype != "ortholog_one2one":
            continue
        if conf.strip() not in {"1", "1.0"}:
            continue
        h, t = _strip_version(h_id), _strip_version(t_id)
        # enforce strict 1:1 even if BioMart mislabels: skip any id reused
        if h in seen_h or t in seen_t:
            continue
        seen_h.add(h)
        seen_t.add(t)
        pairs.append((h, t))
    return pairs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, choices=sorted(TARGETS))
    ap.add_argument("--out", required=True)
    ap.add_argument("--from-file", default=None,
                    help="Parse a local BioMart TSV export instead of querying online.")
    args = ap.parse_args()

    prefix = TARGETS[args.target].homolog_prefix
    if args.from_file:
        with open(args.from_file) as fh:
            text = fh.read()
    else:
        print(f"Querying BioMart: human x {args.target} ({prefix}) ...", file=sys.stderr)
        text = fetch_biomart(HUMAN.biomart_dataset, prefix)

    pairs = parse_one2one(text)
    with open(args.out, "w") as out:
        out.write("human_gene_id\ttarget_gene_id\n")
        for h, t in pairs:
            out.write(f"{h}\t{t}\n")
    print(f"Wrote {len(pairs):,} high-confidence 1:1 orthologs -> {args.out}")


if __name__ == "__main__":
    main()
