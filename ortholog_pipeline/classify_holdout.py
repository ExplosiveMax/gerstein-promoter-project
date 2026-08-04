import torch
from transformers.models.auto.auto_factory import _BaseAutoModelClass
def _lenient_register(cls, config_class, model_class, exist_ok=False):
    cls._model_mapping.register(config_class, model_class, exist_ok=True)
_BaseAutoModelClass.register = classmethod(_lenient_register)

from transformers import AutoTokenizer, AutoModelForSequenceClassification, AutoConfig
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import roc_auc_score, average_precision_score
import numpy as np, csv, argparse, os

ap = argparse.ArgumentParser()
ap.add_argument("--data-dir", required=True)   # e.g. runs/mouse/repeat_gc
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()
seed = args.seed
torch.manual_seed(seed); np.random.seed(seed)
if torch.backends.mps.is_available(): torch.mps.manual_seed(seed)

def load_csv(p):
    seqs, labs = [], []
    with open(p) as f:
        for row in csv.DictReader(f):
            seqs.append(row["sequence"].upper()); labs.append(int(row["label"]))
    return seqs, labs

class DS(Dataset):
    def __init__(self, enc, lab): self.enc, self.lab = enc, lab
    def __len__(self): return len(self.lab)
    def __getitem__(self, i):
        item = {k: v[i] for k, v in self.enc.items()}
        item["labels"] = torch.tensor(self.lab[i]); return item

tr_seqs, tr_lab = load_csv(os.path.join(args.data_dir, "train.csv"))
te_seqs, te_lab = load_csv(os.path.join(args.data_dir, "test.csv"))
print(f"[{args.data_dir}] seed={seed}  train={len(tr_seqs)} test={len(te_seqs)}")

tok = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
tr_enc = tok(tr_seqs, padding=True, truncation=True, max_length=512, return_tensors="pt")
te_enc = tok(te_seqs, padding=True, truncation=True, max_length=512, return_tensors="pt")
tr_loader = DataLoader(DS(tr_enc, tr_lab), batch_size=8, shuffle=True)
te_loader = DataLoader(DS(te_enc, te_lab), batch_size=8)

cfg = AutoConfig.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
cfg.use_flash_attn = False; cfg.num_labels = 2
model = AutoModelForSequenceClassification.from_pretrained(
    "zhihan1996/DNABERT-2-117M", config=cfg, trust_remote_code=True)
dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
model = model.to(dev)
opt = torch.optim.AdamW(model.parameters(), lr=2e-5)

for ep in range(3):
    model.train(); tot = 0
    for i, b in enumerate(tr_loader):
        b = {k: v.to(dev) for k, v in b.items()}
        out = model(**b); out.loss.backward(); opt.step(); opt.zero_grad()
        tot += out.loss.item()
        if (i+1) % 20 == 0: print(f"ep{ep+1} b{i+1}/{len(tr_loader)} loss {out.loss.item():.4f}")
    print(f"ep{ep+1} avg {tot/len(tr_loader):.4f}")

model.eval(); probs, labs = [], []
with torch.no_grad():
    for b in te_loader:
        lb = b.pop("labels"); b = {k: v.to(dev) for k, v in b.items()}
        p = torch.softmax(model(**b).logits, dim=1)[:, 1]
        probs.extend(p.cpu().numpy()); labs.extend(lb.numpy())
auc = roc_auc_score(labs, probs); ap_ = average_precision_score(labs, probs)
print(f"HOLDOUT_RESULT dir={args.data_dir} seed={seed} auroc={auc:.4f} auprc={ap_:.4f}")
