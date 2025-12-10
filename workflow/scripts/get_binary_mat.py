import pandas as pd
import numpy as np
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("--counts", help="RSEM TPM filepath", required=True)
parser.add_argument("--bmatoutfile", help="output CSV path for binary matrix", required=True)

args = parser.parse_args()

counts = args.counts
bmatoutfile = args.bmatoutfile


########################################
# Load in data
########################################

print("Loading in data")

# load in expression matrix
mat = pd.read_table(counts)

########################################
# Make binary matrix
########################################

print("Making binary matrix")

# compute lower thirds
mat["lower_third"] = (
	mat[mat.columns[1:]]	# all but gene symbol column
	.apply(lambda row: np.quantile(row.values, 0.33), axis=1)
)

# binarize by lower third
sampleids = [c for c in mat.columns if c not in ["Hugo_Symbol", "lower_third"]]
mat[sampleids] = (mat[sampleids].ge(mat["lower_third"], axis=0)).astype(int)
mat = mat.drop(columns=["lower_third"])

# save binary matrix
print("Saving binary matrix")
mat.to_csv(bmatoutfile, index=False)

print("done")