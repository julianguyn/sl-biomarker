import pandas as pd
import numpy as np

def score_SL(filename):
	"""
	Compute SL scores from expression matrix

	:param filename: MOHCCN cohort name
	:export: binary matrix of gene KO to data/procdata/binarymat
	:export: SL scores to data/procdata/SLscores
	"""
	print("Starting %s" % filename)
	# load in known SL pairs
	SL_pairs = pd.read_csv("data/procdata/SL_Pairs/human_venn_all_regions.csv")

	# load in expression matrix
	indir = "data/rawdata/" + filename + "_RSEM_TPM_reheadered.tsv"
	mat = pd.read_table(indir)

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
	outdir = "data/procdata/binarymat/" + filename + ".csv"
	mat.to_csv(outdir, index=False)

	# compute SL scores
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
				"geneA": geneA,
				"geneB": geneB,
				"geneA_state": gstateA,
				"geneB_state": gstateB,
				"group": group,
				"outcome": outcome,
				"sample": sample,
			})
	
	result = pd.DataFrame(SL_scores)
	
	# save binary matrix
	outdir = "data/procdata/SLscores/" + filename + ".csv"
	result.to_csv(outdir, index=False)
	print('done')
	

score_SL("INSPIRE")
