from Bio import SeqIO
import pandas as pd
import re

# --- INPUTS ---
fasta_file = "filtered_23S_sequences.fasta"       # your big fasta file
summary_table = "23S_copy_summary.tsv"       # from the previous step
output_table = "23S_copy_summary_filtered.tsv"
output_fasta = "23S_filtered_within5bp.fasta"

# --- PARAMETERS ---
min_len = 2910     # minimum acceptable full-length 23S
max_diff = 5       # maximum bp difference between full-length copies

# --- LOAD SUMMARY TABLE ---
df = pd.read_csv(summary_table, sep="\t")

# --- FUNCTION TO SELECT QUALIFYING ISOLATES ---
def select_within5bp(row):
    copy_cols = [c for c in row.index if c.startswith("copy")]
    lengths = [row[c] for c in copy_cols if not pd.isna(row[c])]
    if len(lengths) < 2:
        return False

    lengths = sorted(lengths, reverse=True)
    # Get indices of the top two
    top_two = lengths[:2]

    # keep if both are ≥ 2900 and within 5bp
    if top_two[0] >= min_len and top_two[1] >= min_len and abs(top_two[0] - top_two[1]) <= max_diff:
        return True
    return False

# --- FILTER DATAFRAME ---
filtered_df = df[df.apply(select_within5bp, axis=1)].copy()
print(f"✅ {len(filtered_df)} isolates retained (within 5bp full-length pairs).")

# --- SAVE FILTERED TABLE ---
filtered_df.to_csv(output_table, sep="\t", index=False)
print(f"💾 Filtered table saved to {output_table}")

# --- LOAD FASTA RECORDS ---
records = list(SeqIO.parse(fasta_file, "fasta"))

# --- BUILD DICTIONARY BY ISOLATE PREFIX ---
fasta_dict = {}
for record in records:
    isolate_id = record.id.split("|")[0]
    fasta_dict.setdefault(isolate_id, []).append(record)

# --- EXTRACT MATCHING SEQUENCES ---
selected_records = []
for isolate in filtered_df["isolate"]:
    if isolate not in fasta_dict:
        continue

    seqs = fasta_dict[isolate]
    # Get their lengths
    lens = [len(r.seq) for r in seqs]
    # Sort by length descending
    seqs_sorted = sorted(seqs, key=lambda r: len(r.seq), reverse=True)

    # Keep sequences within 5 bp of the longest one
    longest_len = len(seqs_sorted[0].seq)
    close_seqs = [r for r in seqs_sorted if abs(len(r.seq) - longest_len) <= max_diff and len(r.seq) >= min_len]

    selected_records.extend(close_seqs)

SeqIO.write(selected_records, output_fasta, "fasta")
print(f"💾 Extracted {len(selected_records)} sequences written to {output_fasta}")
