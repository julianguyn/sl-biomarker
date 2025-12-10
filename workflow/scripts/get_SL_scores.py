import pandas as pd
import numpy as np
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("--bmat", help="output CSV path for binary matrix", required=True)
parser.add_argument("--SL_pairs", help="path to SL pairs (from get networks)", required=True)
parser.add_argument("--SLoutfile", help="Output CSV path for SL scores", required=True)

args = parser.parse_args()

bmat = args.bmat
SL_pairs = args.SL_pairs
SLoutfile = args.SLoutfile


########################################
# Load in data
########################################

print("Loading in data")

# load in known SL pairs
SL_pairs = pd.read_csv(SL_pairs)

# load in expression matrix
mat = pd.read_csv(bmat)

########################################
# Get SL Scores
########################################

print("Computing SL scores")

SL_scores = []
sampleids = [c for c in mat.columns if c not in ["Hugo_Symbol"]]
for geneA, geneB in zip(SL_pairs["geneA"], SL_pairs["geneB"]):
	#print(geneA, geneB)
	# iterate through each SL pair
	for sample in sampleids:
		#if geneA == 'KRAS' and geneB == 'LOC643802':
		#	print(sample)
		# iterate through each patient
		binA = mat.loc[mat["Hugo_Symbol"] == geneA, sample]
		binB = mat.loc[mat["Hugo_Symbol"] == geneB, sample]

		# annotate gene states
		if len(binA) == 0:
			gstateA = None
		else:
			gstateA = "knockout" if binA.iloc[0] == 0 else "normal"
		if len(binB) == 0:
			gstateB = None
		else:
			gstateB = "knockout" if binB.iloc[0] == 0 else "normal"

		# annotate cell fate
		if gstateA is None or gstateB is None:
			outcome = None
			group = None
		elif gstateA == "knockout" and gstateB == "knockout":
			outcome = "doubleKO"
			group = 0
		elif gstateA == "knockout" and gstateB == "normal":
			outcome = "singleKO"
			group = 1
		elif gstateA == "normal" and gstateB == "knockout":
			outcome = "singleKO"
			group = 2
		else:
			outcome = "noKO"
			group = 3

		SL_scores.append({
			"pair": geneA + "_" + geneB,
			"geneA": geneA,
			"geneB": geneB,
			"geneA_state": gstateA,
			"geneB_state": gstateB,
			"group": group,
			"outcome": outcome,
			"sample": sample,
		})

result = pd.DataFrame(SL_scores)

# save SL scores
result.to_csv(SLoutfile, index=False)
print('done')
