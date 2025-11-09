import os
import sys
import glob
import pandas as pd
import time

files=sys.argv[1:-1]

print("Files found:", files)

all_counts = []

for file in files:
    start_time = time.time()
    df = pd.read_csv(file, sep="\t", comment="#")
    
    sample_name = os.path.basename(file).replace("_featurecounts.txt", "")
    df = df[["Geneid", df.columns[-1]]]
    df.rename(columns={df.columns[-1]: sample_name}, inplace=True)
    
    all_counts.append(df)
    
    elapsed = (time.time() - start_time) / 60  # minutes
    print(f"Completed {sample_name} | Rows: {df.shape[0]} | Time: {elapsed:.2f} min")

counts_matrix = all_counts[0]
for df in all_counts[1:]:
    counts_matrix = counts_matrix.merge(df, on="Geneid", how="outer")

output_file=sys.argv[-1]
counts_matrix.to_csv(output_file, index=False)

print("All files processed!")
print("Merged matrix shape:", counts_matrix.shape)