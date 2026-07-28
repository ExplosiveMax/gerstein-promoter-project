import torch
from transformers.models.auto.auto_factory import _BaseAutoModelClass

def _lenient_register(cls, config_class, model_class, exist_ok=False):
    cls._model_mapping.register(config_class, model_class, exist_ok=True)
_BaseAutoModelClass.register = classmethod(_lenient_register)

from transformers import AutoTokenizer, AutoModelForSequenceClassification, AutoConfig
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--cat_a", required=True)
parser.add_argument("--cat_b", required=True)
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()
seed = args.seed
torch.manual_seed(seed)
np.random.seed(seed)
if torch.backends.mps.is_available():
    torch.mps.manual_seed(seed)

class DNADataset(Dataset):
    def __init__(self, encodings, labels):
        self.encodings = encodings
        self.labels = labels
    def __len__(self):
        return len(self.labels)
    def __getitem__(self, idx):
        item = {k: v[idx] for k, v in self.encodings.items()}
        item['labels'] = torch.tensor(self.labels[idx])
        return item

def load_fasta(filepath):
    seqs, cur = [], ""
    for line in open(filepath):
        line = line.strip()
        if line.startswith(">"):
            if cur: seqs.append(cur.upper()); cur = ""
        else: cur += line
    if cur: seqs.append(cur.upper())
    return seqs

a = load_fasta(f"go_{args.cat_a}_pw.fasta")
b = load_fasta(f"go_{args.cat_b}_pw.fasta")
n = min(len(a), len(b))
a, b = a[:n], b[:n]
print(f"{args.cat_a} vs {args.cat_b}: {n} each")

sequences = a + b
labels = [1]*len(a) + [0]*len(b)

train_seqs, test_seqs, train_labels, test_labels = train_test_split(
    sequences, labels, test_size=0.2, random_state=seed, stratify=labels
)

tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
train_enc = tokenizer(train_seqs, padding=True, truncation=True, max_length=512, return_tensors="pt")
test_enc  = tokenizer(test_seqs, padding=True, truncation=True, max_length=512, return_tensors="pt")

train_loader = DataLoader(DNADataset(train_enc, train_labels), batch_size=8, shuffle=True)
test_loader  = DataLoader(DNADataset(test_enc, test_labels), batch_size=8)

config = AutoConfig.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
config.use_flash_attn = False
config.num_labels = 2
model = AutoModelForSequenceClassification.from_pretrained(
    "zhihan1996/DNABERT-2-117M", config=config, trust_remote_code=True)

device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
model = model.to(device)
optimizer = torch.optim.AdamW(model.parameters(), lr=2e-5)

for epoch in range(3):
    model.train()
    for batch in train_loader:
        batch = {k: v.to(device) for k, v in batch.items()}
        outputs = model(**batch)
        loss = outputs.loss
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()
    print(f"Epoch {epoch+1} done")

model.eval()
all_probs, all_labels = [], []
with torch.no_grad():
    for batch in test_loader:
        lb = batch.pop("labels")
        batch = {k: v.to(device) for k, v in batch.items()}
        out = model(**batch)
        probs = torch.softmax(out.logits, dim=1)[:, 1].cpu().numpy()
        all_probs.extend(probs)
        all_labels.extend(lb.numpy())

auc = roc_auc_score(all_labels, all_probs)
print(f"FINAL_PAIR {args.cat_a}_vs_{args.cat_b} seed={seed} auc={auc:.4f}")
