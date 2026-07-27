import torch
from transformers.models.auto.auto_factory import _BaseAutoModelClass

def _lenient_register(cls, config_class, model_class, exist_ok=False):
    cls._model_mapping.register(config_class, model_class, exist_ok=True)
_BaseAutoModelClass.register = classmethod(_lenient_register)

from transformers import AutoTokenizer, AutoModelForSequenceClassification, AutoConfig
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.model_selection import train_test_split
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--high", required=True)
parser.add_argument("--low", required=True)
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
        else:
            cur += line
    if cur: seqs.append(cur.upper())
    return seqs

print("Loading sequences...")
high = load_fasta(args.high)
low  = load_fasta(args.low)
print(f"High: {len(high)}, Low: {len(low)}")

sequences = high + low
labels = [1]*len(high) + [0]*len(low)

train_seqs, test_seqs, train_labels, test_labels = train_test_split(
    sequences, labels, test_size=0.2, random_state=seed, stratify=labels
)
print(f"Train: {len(train_seqs)}, Test: {len(test_seqs)}")

print("Loading tokenizer...")
tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
train_enc = tokenizer(train_seqs, padding=True, truncation=True, max_length=512, return_tensors="pt")
test_enc  = tokenizer(test_seqs, padding=True, truncation=True, max_length=512, return_tensors="pt")

train_loader = DataLoader(DNADataset(train_enc, train_labels), batch_size=8, shuffle=True)
test_loader  = DataLoader(DNADataset(test_enc, test_labels), batch_size=8)

print("Loading model...")
config = AutoConfig.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
config.use_flash_attn = False
config.num_labels = 2
model = AutoModelForSequenceClassification.from_pretrained(
    "zhihan1996/DNABERT-2-117M", config=config, trust_remote_code=True)

device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
print(f"Using device: {device}")
model = model.to(device)

optimizer = torch.optim.AdamW(model.parameters(), lr=2e-5)

print("\nFine-tuning...")
for epoch in range(3):
    model.train()
    total = 0
    for i, batch in enumerate(train_loader):
        batch = {k: v.to(device) for k, v in batch.items()}
        outputs = model(**batch)
        loss = outputs.loss
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()
        total += loss.item()
        if (i+1) % 10 == 0:
            print(f"Epoch {epoch+1}, batch {i+1}/{len(train_loader)}, loss: {loss.item():.4f}")
    print(f"Epoch {epoch+1} complete. Avg loss: {total/len(train_loader):.4f}")

print("\nEvaluating...")
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
print(f"\nAUC Score: {auc:.4f}")
print(f"FINAL_RESULT seed={seed} auc={auc:.4f}")
