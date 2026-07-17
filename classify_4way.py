import torch
from transformers import AutoTokenizer, AutoModelForSequenceClassification, AutoConfig
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import roc_auc_score, confusion_matrix, accuracy_score
from sklearn.model_selection import train_test_split
import numpy as np
import argparse

parser = argparse.ArgumentParser()
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
functional = load_fasta("uorfs_capped.fasta")
pseudo     = load_fasta("uorf_negatives_chr22.fasta")
random_seq = load_fasta("nonutr_4way.fasta")
fourth     = load_fasta("random_4way.fasta")
print(f"uORF: {len(functional)}, withinUTR: {len(pseudo)}, nonUTR: {len(random_seq)}, random: {len(fourth)}")

sequences = functional + pseudo + random_seq + fourth
labels = [0]*len(functional) + [1]*len(pseudo) + [2]*len(random_seq) + [3]*len(fourth)

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
config.num_labels = 4
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
all_logits, all_labels = [], []
with torch.no_grad():
    for batch in test_loader:
        lb = batch.pop("labels")
        batch = {k: v.to(device) for k, v in batch.items()}
        out = model(**batch)
        all_logits.append(out.logits.cpu().numpy())
        all_labels.extend(lb.numpy())

logits = np.concatenate(all_logits, axis=0)
probs = torch.softmax(torch.tensor(logits), dim=1).numpy()
preds = probs.argmax(axis=1)
all_labels = np.array(all_labels)

acc = accuracy_score(all_labels, preds)
# one-vs-rest AUC per class
auc_ovr = roc_auc_score(all_labels, probs, multi_class="ovr", average=None)

names = ["uORF", "withinUTR", "nonUTR", "random"]
print(f"\n=== RESULTS seed={seed} ===")
print(f"Overall accuracy: {acc:.4f}")
for i, n in enumerate(names):
    print(f"AUC {n} vs rest: {auc_ovr[i]:.4f}")
print("\nConfusion matrix (rows=true, cols=pred):")
print("           " + "  ".join(f"{n[:6]:>6}" for n in names))
cm = confusion_matrix(all_labels, preds)
for i, n in enumerate(names):
    print(f"{n[:10]:>10} " + "  ".join(f"{v:>6}" for v in cm[i]))
print(f"\nFINAL_RESULT seed={seed} acc={acc:.4f}")