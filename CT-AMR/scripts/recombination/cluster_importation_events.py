import pandas as pd

# Parameters
INPUT_FILE = "CF_Core_alignment.importation_status.txt"   # your file
OUTPUT_FILE = "core_alignment_recomb_hotspots.tsv"
MAX_DIST = 1000      # max cluster span
CUSHION = 100        # extra tolerance (bp)

# Load data
df = pd.read_csv(INPUT_FILE, sep="\t", header=0, names=["Node", "Beg", "End"])

# Sort by genomic position
df = df.sort_values("Beg").reset_index(drop=True)

clusters = []
current_start = df.loc[0, "Beg"]
current_end = df.loc[0, "End"]
current_nodes = [df.loc[0, "Node"]]

for i in range(1, len(df)):
    beg, end, node = df.loc[i, ["Beg", "End", "Node"]]
    
    # If adding this interval keeps cluster span <= (MAX_DIST + CUSHION)
    if (max(current_end, end) - min(current_start, beg)) <= (MAX_DIST + CUSHION):
        current_start = min(current_start, beg)
        current_end = max(current_end, end)
        current_nodes.append(node)
    else:
        # Save current cluster
        clusters.append([current_start, current_end, len(current_nodes), ",".join(set(current_nodes))])
        
        # Start new cluster
        current_start, current_end = beg, end
        current_nodes = [node]

# Save last cluster
clusters.append([current_start, current_end, len(current_nodes), ",".join(set(current_nodes))])

# Write output
out_df = pd.DataFrame(clusters, columns=["Cluster_Start", "Cluster_End", "Num_Events", "Nodes"])
out_df = out_df.sort_values("Num_Events", ascending=False)  # most frequent hotspots first
out_df.to_csv(OUTPUT_FILE, sep="\t", index=False)

print(f"Saved clustered recombination hotspots to {OUTPUT_FILE}")
