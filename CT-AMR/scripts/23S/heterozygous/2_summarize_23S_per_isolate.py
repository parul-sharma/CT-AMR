from Bio import SeqIO
import pandas as pd

input_fasta = "filtered_23S_sequences.fasta"
output_table = "23S_copy_summary.tsv"

# Dictionary: isolate -> list of sequence lengths
isolate_seqs = {}

for record in SeqIO.parse(input_fasta, "fasta"):
    isolate = record.id.split("|")[0]
    seq_len = len(record.seq)
    isolate_seqs.setdefault(isolate, []).append(seq_len)

# Build summary table
rows = []
for isolate, lengths in isolate_seqs.items():
    row = {
        "isolate": isolate,
        "num_23S_copies": len(lengths)
    }
    # Add columns for each copy length
    for i, l in enumerate(lengths, start=1):
        row[f"copy{i}_length"] = l
    rows.append(row)

# Convert to DataFrame
df = pd.DataFrame(rows)

# Sort by isolate name
df = df.sort_values("isolate")

# Save
df.to_csv(output_table, sep="\t", index=False)

print(f"✅ Summary table written to {output_table}")
print(df.head())
