import random
random.seed(42)

def gc_content(seq):
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq)

# Load chr22 into memory
print("Loading chr22...")
chr22_seq = ""
with open('chr22.fa') as f:
    for line in f:
        if not line.startswith('>'):
            chr22_seq += line.strip()

print(f"chr22 length: {len(chr22_seq)}")

region_size = 1000
num_regions = 500
min_gc = 0.40
max_gc = 0.60
max_repeat = 0.25

regions = []
attempts = 0

while len(regions) < num_regions:
    attempts += 1
    start = random.randint(0, len(chr22_seq) - region_size)
    end = start + region_size
    seq = chr22_seq[start:end]
    
    # Check repeat content
    lowercase = sum(1 for c in seq if c.islower())
    if lowercase / len(seq) > max_repeat:
        continue
    
    # Check GC content
    gc = gc_content(seq)
    if gc < min_gc or gc > max_gc:
        continue
    
    regions.append((start, end))

print(f"Found {len(regions)} GC-matched regions after {attempts} attempts")

with open('random_regions.bed', 'w') as f:
    for start, end in regions:
        f.write(f'chr22\t{start}\t{end}\n')

print("Saved to random_regions.bed")