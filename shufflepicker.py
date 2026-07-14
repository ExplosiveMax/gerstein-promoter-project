import random

def load_fasta(filepath):
    sequences = []
    headers = []
    current_seq = ""
    current_header = ""
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    sequences.append(current_seq.upper())
                    headers.append(current_header)
                current_seq = ""
                current_header = line
            else:
                current_seq += line
        if current_seq:
            sequences.append(current_seq.upper())
            headers.append(current_header)
    return headers, sequences

print("Loading promoter sequences...")
headers, sequences = load_fasta("promoters.fasta")
print(f"Loaded {len(sequences)} sequences")

print("Shuffling sequences...")
with open("random_sequences.fa", "w") as f:
    for i, (header, seq) in enumerate(zip(headers, sequences)):
        seq_list = list(seq)
        random.seed(i)
        random.shuffle(seq_list)
        shuffled = ''.join(seq_list)
        f.write(f">shuffled_{i}\n{shuffled}\n")

print("Done! Saved to random_sequences.fa")