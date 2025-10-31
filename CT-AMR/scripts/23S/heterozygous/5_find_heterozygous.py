import pandas as pd

# --- INPUT ---
input_file = "mutations_per_sample.csv"
output_file = "23S_homo_hetero_summary.tsv"

# --- LOAD DATA ---
df = pd.read_csv(input_file)

# Extract isolate ID without gene/copy info
df['Isolate'] = df['Sample_ID'].apply(lambda x: x.split('|')[0])

# Group by isolate
summary_list = []

for isolate, group in df.groupby('Isolate'):
    # List of mutation sets per copy
    copy_mut_sets = [set(m.strip() for m in muts.split(';')) for muts in group['Mutations']]
    
    # Homozygous = mutations present in all copies
    homozygous = set.intersection(*copy_mut_sets)
    # Heterozygous = mutations present in some but not all copies
    heterozygous = set.union(*copy_mut_sets) - homozygous
    
    summary_list.append({
        'Isolate': isolate,
        'num_copies': len(copy_mut_sets),
        'homozygous_mutations': '; '.join(sorted(homozygous)) if homozygous else '',
        'heterozygous_mutations': '; '.join(sorted(heterozygous)) if heterozygous else ''
    })

summary_df = pd.DataFrame(summary_list)

# Save summary
summary_df.to_csv(output_file, sep='\t', index=False)
print(f"✅ Homozygous/heterozygous summary saved to {output_file}")
