import random
random.seed(42)

# Load chr22 sequence into memory so we can check for repeats
print("Loading chr22...")
chr22_seq = ""
with open('chr22.fa') as f:
    for line in f:
        if not line.startswith('>'):
            chr22_seq += line.strip()

print(f"chr22 length: {len(chr22_seq)}")

region_size = 1000
num_regions = 500
max_repeat_fraction = 0.25  # reject if >25% lowercase (repetitive)

regions = []
attempts = 0

while len(regions) < num_regions:
    attempts += 1
    start = random.randint(0, len(chr22_seq) - region_size)
    end = start + region_size
    seq = chr22_seq[start:end]
    lowercase = sum(1 for c in seq if c.islower())
    if lowercase / len(seq) <= max_repeat_fraction:
        regions.append((start, end))

print(f"Found {len(regions)} clean regions after {attempts} attempts")

with open('random_regions.bed', 'w') as f:
    for start, end in regions:
        f.write(f'chr22\t{start}\t{end}\n')

print("Saved to random_regions.bed")