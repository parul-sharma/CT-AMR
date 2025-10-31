from Bio import SeqIO

# Input and output files
input_fasta = "all_23S_rRNA_sequences.fasta"
output_fasta = "filtered_23S_sequences.fasta"

# List of low-quality genomes to remove
remove_list = {
    "E112-66_wgs_batch3",
    "E112-67_wgs_batch3",
    "E112-13_wgs_batch3",
    "E112-10_wgs_batch3",
    "E112-75_wgs_batch3",
    "E112-57_wgs_batch3",
    "929C_wgs_batch1",
    "E112-72_wgs_batch3",
    "1270_wgs_batch2",
    "1077F1C_wgs_batch1",
    "SAMEA1398201",
    "Undetermined_wgs_batch3",
    "106F2C_wgs_batch1",
    "869F2V_wgs_batch1",
    "E112-22_wgs_batch3",
    "E112-28_wgs_batch3",
    "E112-40_wgs_batch3",
    "E112-48_wgs_batch3",
    "E112-43_wgs_batch3",
    "E112-64_wgs_batch3",
    "E112-44_wgs_batch3",
    "E112-73_wgs_batch3",
    "E112-61_wgs_batch3",
    "421_wgs_batch2",
    "869F2C_wgs_batch1",
    "E112-49_wgs_batch3",
    "E112-17_wgs_batch3",
    "E112-39_wgs_batch3",
    "E112-26_wgs_batch3",
    "E112-36_wgs_batch3",
    "E112-37_wgs_batch3",
    "SAMEA1398219",
    "SAMEA1527843",
    "SAMEA1527913",
    "SAMEA1398212",
    "SAMEA1559202",
    "SAMEA1398228"
}

kept_records = []
found_removed = set()

for record in SeqIO.parse(input_fasta, "fasta"):
    isolate = record.id.split("|")[0]
    if isolate in remove_list:
        found_removed.add(isolate)
        continue
    kept_records.append(record)

SeqIO.write(kept_records, output_fasta, "fasta")

missing = remove_list - found_removed

print(f"✅ Wrote {len(kept_records)} records to {output_fasta}")
print(f"❌ Removed {len(found_removed)} isolates present in FASTA.")
if missing:
    print(f"⚠️  {len(missing)} isolates from remove list not found:")
    for m in sorted(missing):
        print("   ", m)
