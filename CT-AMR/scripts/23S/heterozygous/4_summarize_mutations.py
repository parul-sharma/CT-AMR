from Bio import AlignIO
from collections import defaultdict, Counter
import csv
import sys

# --- USER INPUT ---
alignment_file = "aligned_23S_within5pb.fasta"  # aligned protein FASTA file
reference_id = "ref_NC_000117.1_878039-880902_Chlamydia_trachomatis_D/UW-3/CX"           # replace with your reference sequence ID
# ------------------

# Load aligned sequences
alignment = AlignIO.read(alignment_file, "fasta")

# Get reference sequence
reference = next((rec for rec in alignment if rec.id == reference_id), None)
if reference is None:
    sys.exit(f"Reference ID '{reference_id}' not found in the alignment.")

ref_seq = str(reference.seq)
mutations_by_sample = defaultdict(list)
mutation_counter = Counter()

# Extract mutations with +1 position shift
for record in alignment:
    if record.id == reference_id:
        continue

    sample_id = record.id
    sample_seq = str(record.seq)

    for pos, (ref_aa, sample_aa) in enumerate(zip(ref_seq, sample_seq), start=1):
        real_pos = pos  # Shift position by +1

        if ref_aa != sample_aa:
            if sample_aa == "-":
                mutation = f"{ref_aa}{real_pos}del"
            elif ref_aa == "-":
                mutation = f"{real_pos}ins{sample_aa}"
            else:
                mutation = f"{ref_aa}{real_pos}{sample_aa}"

            mutations_by_sample[sample_id].append(mutation)
            mutation_counter[mutation] += 1

# Save mutations per sample to CSV
with open("mutations_per_sample.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Sample_ID", "Mutations"])
    for sample, muts in mutations_by_sample.items():
        writer.writerow([sample, "; ".join(muts) if muts else "No mutations"])

# Save mutation summary to CSV
with open("mutation_summary.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Mutation", "Count"])
    for mutation, count in mutation_counter.most_common():
        writer.writerow([mutation, count])

print("Mutation extraction complete!")
print("Results saved to 'mutations_per_sample.csv' and 'mutation_summary.csv'")
