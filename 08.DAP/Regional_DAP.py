import os
import gc
import argparse
from warnings import warn

import numpy as np
import pandas as pd
import scanpy as sc

os.chdir(os.getenv("HOME"))

# arguments
parser = argparse.ArgumentParser(description="Regional DAP")
parser.add_argument(
	'--adata', '-adata', type=str, required=True
)
parser.add_argument(
	'--celltype_label', '-cl', type=str, required=True
)
parser.add_argument(
	'--subset_label', '-subset_label', type=str, default=None
)
parser.add_argument(
	'--subset_by', '-subset_by', type=str, default=None, nargs='+'
)
parser.add_argument(
	'--region_label', '-rl', type=str, required=True
)
parser.add_argument(
	'--selected_region', '-r', type=str, required=True
)
parser.add_argument(
	'--out_dir', '-out', type=str, required=True
)
## optional
parser.add_argument(
	'--method', '-m', type=str, default='wilcoxon'
)
parser.add_argument(
	'--pval_cutoff', '-pv', type=float, default=0.05
)
parser.add_argument(
	'--pct_cutoff', '-pct', type=float, default=0.1
)
parser.add_argument(
	'--log2fc_cutoff', '-log2fc', type=float, default=0.
)
## norm
parser.add_argument(
	'--norm', '-norm', default=False, action="store_true"
)

args = parser.parse_args()

# parse
out_dir = args.out_dir
os.makedirs(out_dir, exist_ok=True)

subset_label = args.subset_label
subset_list = args.subset_by

clabel = args.celltype_label

# load
adata = sc.read_h5ad(args.adata)

# norm
if args.norm:
	print("Running normalize and log1p(base=2) for adata")
	sc.pp.normalize_total(adata, target_sum=10**4)
	sc.pp.log1p(adata, base=2)

if subset_label is not None:
	adata = adata[adata.obs[subset_label].isin(subset_list), ].copy()
	gc.collect()
	print(f"Subsetted adata for {subset_list} by {subset_label}")

# run
print("This is the adata:")
print(adata)

celltypes = adata.obs[clabel].unique()

# Run
all_results = []
exec(open("project/multiomics/SCENIC/script/my_runDAP.py").read())

for cl in celltypes:
	gc.collect()
	print(f"\nRunning for {cl}:")
	sub_adata = adata[adata.obs[clabel] == cl, ]
	
	try:
		codap_results = my_runDAP(
			sub_adata,
			group_by=args.region_label,
			selected_group=args.selected_region,
			method=args.method,
			pval_cutoff=args.pval_cutoff,
			pct_cutoff=args.pct_cutoff,
			log2fc_cutoff=args.log2fc_cutoff,
			save_dir=os.path.join(out_dir, cl)
		)

		if codap_results is not None:
			for k in codap_results:
				df = codap_results[k]
				if df is None:
					pass
				elif df.shape[0] > 0:
					df["sign"] = k
					df["cluster"] = cl
					all_results.append(df)

					res_dir = os.path.join(out_dir, cl, k)
					os.makedirs(res_dir, exist_ok=True)
					df.to_csv(
						os.path.join(res_dir, "cons.markers_df.csv"), index=False
					)
					print(f"There are {df.shape[0]} DAPs for {k}")
				else:
					warn(f"There is no DAP for {k}")


	except Exception as e:
		warn(f"Error for {cl}: {e}")
		continue


# merge
if len(all_results) > 0:
	mdf = pd.concat(all_results)

	# mdf.to_csv(
	# 	os.path.join(out_dir, "agg.cons.mdf.csv.gz"), index=False
	# )
	mdf.reset_index(drop=True).to_feather(
		os.path.join(out_dir, "agg.cons.mdf.feather")
	)
else:
	warn(f"There is no DAP")

