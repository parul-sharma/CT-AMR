import csv
from collections import Counter

input_file = "mutations_per_sample.csv"
output_file = "mutation_combination_summary.csv"

# Read mutations per sample
combination_counter = Counter()

with open(input_file, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        muts = row["Mutations"]
        if muts.strip().lower() in ("no mutations", "", "na", "none"):
            combo = ""  # Represent no mutation as empty string
        else:
            muts_list = [m.strip() for m in muts.split(";")]
            # Sort mutations alphabetically to normalize combinations
            muts_list.sort()
            combo = ", ".join(muts_list)
        combination_counter[combo] += 1

# Write summary of unique mutation combinations
with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Mutation_Combination", "Number_of_Samples"])
    for combo, count in combination_counter.most_common():
        display_combo = combo if combo else "(No mutations)"
        writer.writerow([display_combo, count])

print(f"Unique mutation combination summary saved to '{output_file}'.")
