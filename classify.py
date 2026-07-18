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
    sequences = []
    current_seq = ""
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    sequences.append(current_seq.upper())
                current_seq = ""
            else:
                current_seq += line
        if current_seq:
            sequences.append(current_seq.upper())
    return sequences

print("Loading sequences...")
promoters = load_fasta("promoters.fasta")
randoms = load_fasta("random_sequences.fa")
print(f"Random seed: {seed}")
print(f"Promoters: {len(promoters)}, Random: {len(randoms)}")

sequences = promoters + randoms
labels = [1] * len(promoters) + [0] * len(randoms)

# train test
train_seqs, test_seqs, train_labels, test_labels = train_test_split(
    sequences, labels, test_size=0.2, random_state=seed, stratify=labels
)
print(f"Train: {len(train_seqs)}, Test: {len(test_seqs)}")

print("Loading tokenizer...")
tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)

print("Tokenizing...")
train_encodings = tokenizer(train_seqs, padding=True, truncation=True, max_length=512, return_tensors="pt")
test_encodings = tokenizer(test_seqs, padding=True, truncation=True, max_length=512, return_tensors="pt")

train_dataset = DNADataset(train_encodings, train_labels)
test_dataset = DNADataset(test_encodings, test_labels)

train_loader = DataLoader(train_dataset, batch_size=8, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=8)

print("Loading model...")
config = AutoConfig.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
config.use_flash_attn = False
config.num_labels = 2
model = AutoModelForSequenceClassification.from_pretrained(
    "zhihan1996/DNABERT-2-117M", config=config, trust_remote_code=True
)

device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
print(f"Using device: {device}")
model = model.to(device)

optimizer = torch.optim.AdamW(model.parameters(), lr=2e-5)
num_epochs = 3

print("\nFine-tuning...")
for epoch in range(num_epochs):
    model.train()
    total_loss = 0
    for i, batch in enumerate(train_loader):
        batch = {k: v.to(device) for k, v in batch.items()}
        outputs = model(**batch)
        loss = outputs.loss
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()
        total_loss += loss.item()
        if (i + 1) % 10 == 0:
            print(f"Epoch {epoch+1}, batch {i+1}/{len(train_loader)}, loss: {loss.item():.4f}")
    print(f"Epoch {epoch+1} complete. Avg loss: {total_loss/len(train_loader):.4f}")

print("\nEvaluating...")
model.eval()
all_probs = []
all_labels = []

with torch.no_grad():
    for batch in test_loader:
        labels_batch = batch.pop('labels')
        batch = {k: v.to(device) for k, v in batch.items()}
        outputs = model(**batch)
        probs = torch.softmax(outputs.logits, dim=1)[:, 1]
        all_probs.extend(probs.cpu().numpy())
        all_labels.extend(labels_batch.numpy())

auc = roc_auc_score(all_labels, all_probs)
print(f"\nAUC Score: {auc:.4f}")
print(f"FINAL_RESULT seed={seed} auc={auc:.4f}")
print("0.5 = random chance, higher = model can distinguish promoters from random DNA")