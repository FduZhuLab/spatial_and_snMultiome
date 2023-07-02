# packages
import os
import gc
import dill
import pickle
import argparse

import warnings

import numpy as np
import pandas as pd

## post-anlys
from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
from scenicplus.eregulon_enrichment import score_eRegulons, binarize_AUC
from scenicplus.cistromes import generate_pseudobulks, TF_cistrome_correlation
from scenicplus.RSS import regulon_specificity_scores
from scenicplus.dimensionality_reduction import run_eRegulons_umap, harmony

os.chdir(os.getenv("HOME"))

# arguments
## io
parser = argparse.ArgumentParser(description="SCENICplus Post-Anlys")
parser.add_argument(
	'--input_dir', '-i', type=str, required=True
)
parser.add_argument(
	'--output_dir', '-o', type=str, default=None
)
parser.add_argument(
	'--raw_output', '-raw', type=str, default=None
)

## filter and post-Anlys
parser.add_argument(
	'--groups', '-g', type=str, required=True, nargs='+',
)
parser.add_argument(
	'--suffix', '-suf', type=str, default="_filtered"
)
parser.add_argument(
	'--auc_key', '-k', type=str, default=None
)
parser.add_argument(
	'--gene_auc_thres', '-gthres', type=float, default=0.05
)
parser.add_argument(
	'--region_auc_thres', '-rthres', type=float, default=0.05
)
parser.add_argument(
	'--gene_auc_norm', '-gnorm', default=False, action='store_true'
)
parser.add_argument(
	'--region_auc_norm', '-rnorm', default=False, action='store_true'
)

## anlys optional
parser.add_argument(
	'--batch', '-batch', type=str, default=None
)
parser.add_argument(
	'--bin', '-bin', default=True, action='store_true'
)
parser.add_argument(
	'--no_bin', '-no_bin', dest="bin", action='store_false'
)
parser.add_argument(
	'--rss', '-rss', default=True, action='store_true'
)
parser.add_argument(
	'--no_rss', '-no_rss', dest="rss", action='store_false'
)
parser.add_argument(
	'--umap', '-umap', default=True, action='store_true'
)
parser.add_argument(
	'--no_umap', '-no_umap', dest="umap", action='store_false'
)
parser.add_argument(
	'--rumap', '-rumap', default=True, action='store_true'
)
parser.add_argument(
	'--no_rumap', '-no_rumap', dest="rumap", action='store_false'
)

## save optional
parser.add_argument(
	'--save_raw', '-save_raw', default=True, action='store_true'
)
parser.add_argument(
	'--no_raw', '-no_raw', dest="save_raw", action='store_false'
)
parser.add_argument(
	'--save_filter', '-save_filter', default=True, action='store_true'
)
parser.add_argument(
	'--no_filter', '-no_filter', dest="save_filter", action='store_false'
)

## others
parser.add_argument(
	'--ncpu', '-n', type=int, required=True
)

## parse
args = parser.parse_args()

input_dir = args.input_dir
output_dir = args.output_dir

if output_dir is None:
	output_dir = os.path.join(input_dir, "PostAnlys")

os.makedirs(output_dir, exist_ok=True)
print(f"The output dir is {output_dir}")	

auc_key = args.auc_key
suffix = args.suffix
if auc_key is None:
	auc_key = f"eGRN_AUC{suffix}"

print(f"The filtered AUC key is {auc_key}")

# Load
print("Loading raw SCENIC+ object...")

with open(os.path.join(input_dir, "scplus_obj.pkl"), "rb") as f:
    scplus_obj = dill.load(f)
	
## test groups
groups = args.groups
if type(groups) == list:
	print(f"The groups are: {groups}")
	foo = scplus_obj.metadata_cell.loc[:, groups]
else:
	raise(ValueError("The groups argument must be a list"))

## test batch
batch = args.batch
if batch is not None:
	print(f"The batch is: {batch}")
	foo = scplus_obj.metadata_cell.loc[:, batch]
	if batch in groups:
		raise(ValueError("The batch cannot be in groups"))

# Save raw
if args.save_raw:
	if args.raw_output is None:
		raw_output = os.path.join(output_dir, "raw")

	print(f"The raw output dir is {raw_output}")
	print("Saving raw SCENIC+ eRegulon...")
	os.makedirs(raw_output, exist_ok=True)

	with open(os.path.join(raw_output, "scplus_uns.dill"), "wb") as f:
		dill.dump(scplus_obj.uns, f)

	scplus_obj.uns['eRegulon_metadata'].reset_index(drop=True).to_feather(
		os.path.join(raw_output, "eRegulon_metadata.feather")
	)

	with open(os.path.join(raw_output, "eRegulon_signatures.pkl"), "wb") as f:
		pickle.dump(scplus_obj.uns["eRegulon_signatures"], f)

	## Save raw post-anlys data
	print("Saving raw SCENIC+ post-anlys data...")

	AUC = scplus_obj.uns['eRegulon_AUC']
	for k in AUC:
		pd.DataFrame(AUC[k]).reset_index().to_feather(
			os.path.join(raw_output, f"eRegulon_AUC_{k}.feather")
		)

	AUC_T = scplus_obj.uns['eRegulon_AUC_thresholds']
	for k in AUC_T:
		X = pd.DataFrame(AUC_T[k])
		X.columns = ["thres"]
		X.to_csv(
			os.path.join(raw_output, f"eRegulon_AUC_thresholds_{k}.csv")
		)

	TFCorr = scplus_obj.uns['TF_cistrome_correlation']
	for k in TFCorr:
		TFCorr[k].to_csv(
			os.path.join(raw_output, f"TFCorr_{k}.csv"), index=False
		)

	RSS = scplus_obj.uns["RSS"]
	for k in RSS:
		RSS[k].to_csv(
			os.path.join(raw_output, f"eRegulon_RSS_{k}.csv")
		)


# Filter
print("Filterimg SCENIC+ eRegulon...")

apply_std_filtering_to_eRegulons(scplus_obj)

## save
if args.save_filter:
	print("Saving filtered SCENIC+ eRegulon...")

	scplus_obj.uns['eRegulon_metadata_filtered'].reset_index(drop=True).to_feather(
		os.path.join(output_dir, "eRegulon_metadata_filtered.feather")
	)

	with open(os.path.join(output_dir, "eRegulon_signatures_flitered.pkl"), "wb") as f:
		pickle.dump(scplus_obj.uns["eRegulon_signatures_filtered"], f)

	# Jaccard
	exec(open("project/multiomics/SCENIC/script/PostAnlys/my_scplus_jaccard.py").read())

	for k in scplus_obj.uns['eRegulon_signatures_filtered']:
		print(f"Jaccard for {k}...")

		jaccard = my_scplus_jaccard(
			scplus_obj, gene_or_region_based = k,
			signature_key = 'eRegulon_signatures_filtered'
		)

		jaccard.to_csv(
			os.path.join(output_dir, f"eGRN_jaccard_{k}.csv")
		)
	
## AUC enrich
with open(os.path.join(input_dir, "region_ranking.pkl"), "rb") as f:
    region_ranking = dill.load(f)

with open(os.path.join(input_dir, "gene_ranking.pkl"), "rb") as f:
    gene_ranking = dill.load(f)

params = {
	"gene": {
		"rank": gene_ranking, "thres": args.gene_auc_thres, "norm": args.gene_auc_norm
	},
    "region": {
		"rank": region_ranking, "thres": args.region_auc_thres, "norm": args.region_auc_norm
	}
}

print("AUC enrichment...")

for t in params:
    score_eRegulons(
        scplus_obj,
        ranking = params[t]["rank"],
        eRegulon_signatures_key = 'eRegulon_signatures_filtered',
        key_added = auc_key,
        enrichment_type= t,
        auc_threshold = params[t]["thres"],
        normalize = params[t]["norm"],
        n_cpu = args.ncpu
    )

del gene_ranking, region_ranking, params
gc.collect()

## save AUC
AUC = scplus_obj.uns[auc_key]
for k in AUC:
    pd.DataFrame(AUC[k]).reset_index().to_feather(
        os.path.join(output_dir, f"{auc_key}_{k}.feather")
    )

# Harmony
if batch is not None:
	print("Harmony correcting...")
	auc_key_c = f"{auc_key}_correct"

	harmony(
		scplus_obj, variable = batch,
		auc_key = auc_key,
		out_key = auc_key_c,
		signature_keys = list(scplus_obj.uns[auc_key].keys())
	)
	
	## save correct AUC
	AUC_c = scplus_obj.uns[auc_key_c]
	for k in AUC_c:
		pd.DataFrame(AUC_c[k]).reset_index().to_feather(
			os.path.join(output_dir, f"{auc_key_c}_{k}.feather")
		)
else:
	args.rumap = False

# Binarise
if args.bin:
	print("Binarising...")
	auc_t_key = f"{auc_key}_thresholds"

	binarize_AUC(
		scplus_obj,
		auc_key = auc_key,
		out_key = auc_t_key,
		signature_keys = list(scplus_obj.uns[auc_key].keys()),
		n_cpu = args.ncpu
	)

	AUC_T = scplus_obj.uns[auc_t_key]
	for k in AUC_T:
		X = pd.DataFrame(AUC_T[k])
		X.columns = ["thres"]
		X.to_csv(
			os.path.join(output_dir, f"{auc_t_key}_{k}.csv")
		)

# For each group
corr_dir = f"{output_dir}/TFCorr"
rss_dir = f"{output_dir}/RSS"
umap_dir = f"{output_dir}/UMAP"

os.makedirs(corr_dir, exist_ok=True)
os.makedirs(rss_dir, exist_ok=True)
os.makedirs(umap_dir, exist_ok=True)

for g in groups:
	print(f"Post-Anlys for {g}")

	for k in AUC:
		print(f"TF-Cis Correlation for {k}...")
		# TF-Cis Correlation
		k2 = f"{g}{suffix}_{k}"
		generate_pseudobulks(
			scplus_obj, variable = g, auc_key = auc_key, signature_key = k
		)

		with warnings.catch_warnings():
			warnings.simplefilter(action='ignore', category=FutureWarning)
			TF_cistrome_correlation(
				scplus_obj, variable = g, auc_key = auc_key, signature_key = k,
				use_pseudobulk = True, out_key = k2
			)

		## save
		scplus_obj.uns['TF_cistrome_correlation'][k2].to_csv(
			f"{corr_dir}/TFCorr_{k2}.csv", index=False
		)

		# RSS
		if args.rss:
			print(f"RSS for {k}...")
			regulon_specificity_scores(
				scplus_obj, variable = g, auc_key = auc_key, signature_keys = [k],
				out_key_suffix = f"{suffix}_{k}"
			)

			## save
			scplus_obj.uns["RSS"][k2].to_csv(f"{rss_dir}/eGRN_RSS_{k2}.csv")

	# RSS for all
	if args.rss:
		print("RSS for all...")
		regulon_specificity_scores(
			scplus_obj, variable = g, auc_key = auc_key,
			signature_keys = list(AUC.keys()), out_key_suffix = f"{suffix}_merge"
		)

		## save
		x = f"{g}{suffix}_merge"
		scplus_obj.uns["RSS"][x].to_csv(f"{rss_dir}/eGRN_RSS_{x}.csv")


# RunUMAP
if args.rumap:
	auc_key = auc_key_c
	AUC = AUC_c

if args.umap:
	for k in AUC:
		print(f"RunUMAP for {k}...")
		k2 = f'eRegulons{suffix}_{k}_UMAP'
		run_eRegulons_umap(
			scplus_obj = scplus_obj, auc_key = auc_key,
			reduction_name = k2, signature_keys = [k]
		)
		scplus_obj.dr_cell[k2].to_csv(f"{umap_dir}/{k2}.csv")

	print("RunUMAP for all...")
	k2 = f'eRegulons{suffix}_merge_UMAP'
	run_eRegulons_umap(
		scplus_obj = scplus_obj, auc_key = auc_key,
		reduction_name = k2,
		signature_keys = list(AUC.keys())
	)
	scplus_obj.dr_cell[k2].to_csv(f"{umap_dir}/{k2}.csv")

# Save
print("Saving SCENIC+ eRegulon and post-anlys data")

with open(os.path.join(output_dir, "scplus_uns.dill"), "wb") as f:
    dill.dump(scplus_obj.uns, f)

