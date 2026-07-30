import torch
from transformers.models.auto.auto_factory import _BaseAutoModelClass

def _lenient_register(cls, config_class, model_class, exist_ok=False):
    cls._model_mapping.register(config_class, model_class, exist_ok=True)
_BaseAutoModelClass.register = classmethod(_lenient_register)

from transformers import AutoTokenizer, AutoModelForSequenceClassification
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--model_dir", default="onoff_liver_model_seed42")
parser.add_argument("--fasta", default="onoff_Liver_high.fasta")
parser.add_argument("--n_seqs", type=int, default=5)
parser.add_argument("--out_prefix", required=True)
args = parser.parse_args()

BASES = ["A", "C", "G", "T"]

def load_fasta(filepath, n):
    seqs, cur = [], ""
    for line in open(filepath):
        line = line.strip()
        if line.startswith(">"):
            if cur: seqs.append(cur.upper()); cur = ""
        else:
            cur += line
    if cur: seqs.append(cur.upper())
    return seqs[:n]

print("Loading model...")
tokenizer = AutoTokenizer.from_pretrained(args.model_dir, trust_remote_code=True)
model = AutoModelForSequenceClassification.from_pretrained(args.model_dir, trust_remote_code=True)
device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
model = model.to(device)
model.eval()

def score(seq):
    enc = tokenizer(seq, return_tensors="pt", truncation=True, max_length=512).to(device)
    with torch.no_grad():
        logits = model(**enc).logits
        prob = torch.softmax(logits, dim=1)[0, 1].item()  # P(high expression)
    return prob

print(f"Loading {args.n_seqs} sequences from {args.fasta}...")
seqs = load_fasta(args.fasta, args.n_seqs)

for si, seq in enumerate(seqs):
    print(f"\n=== Sequence {si} (len={len(seq)}) ===")
    baseline = score(seq)
    print(f"Baseline P(high): {baseline:.4f}")

    importance = np.zeros(len(seq))
    for pos in range(len(seq)):
        orig = seq[pos]
        deltas = []
        for b in BASES:
            if b == orig: continue
            mutated = seq[:pos] + b + seq[pos+1:]
            deltas.append(abs(score(mutated) - baseline))
        importance[pos] = np.mean(deltas)  # avg abs shift across 3 possible mutations

    np.save(f"ism_importance_{args.out_prefix}_seq{si}.npy", importance)
    top10 = np.argsort(importance)[::-1][:10]
    print("Top 10 most important positions (by avg |delta|):")
    for p in sorted(top10):
        print(f"  pos {p:>4}  base={seq[p]}  importance={importance[p]:.4f}")

print("\nSaved per-position importance arrays as ism_importance_seq*.npy")
