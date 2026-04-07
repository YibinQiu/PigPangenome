import sys
import pandas as pd

# Get input/output paths from command line
parfile = sys.argv[1]  # Input parquet file
outfile = sys.argv[2]  # Output compressed tsv file

# Read parquet and save as gzipped tsv
df = pd.read_parquet(parfile)
df.to_csv(outfile, compression="gzip", sep="\t", index=False)