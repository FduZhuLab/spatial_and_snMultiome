# packages
import os
import gc
import pickle
import argparse

import anndata
from scipy.sparse import csr_matrix
import scanpy as sc

from pycisTopic.diff_features import impute_accessibility

os.chdir(os.getenv("HOME"))

# arguments
parser = argparse.ArgumentParser(description="Imputed ACC")
parser.add_argument(
	'--home_dir', '-i', type=str, required=True
)
parser.add_argument(
	'--save_path', '-o', type=str, default=None
)
parser.add_argument(
	'--norm', '-norm', default=False, action='store_true'
)
parser.add_argument(
	'--sparse', '-sparse', default=False, action='store_true'
)
parser.add_argument(
	'--scale_factor', '-scale', type=int, default=10**4
)

args = parser.parse_args()

# Load
home_dir = args.home_dir
save_path = args.save_path
if save_path is None:
	save_path = os.path.join(home_dir, "imputed_ACC.h5ad")

os.makedirs(os.path.dirname(save_path), exist_ok=True)

with open(os.path.join(home_dir, "topic_obj.pkl"), "rb") as f:
    topic_obj = pickle.load(f)

# Imputed
imputed_acc_obj = impute_accessibility(
	topic_obj, selected_cells=None, selected_regions=None
)

# create
X = imputed_acc_obj.mtx.T
if args.sparse:
	X = csr_matrix(X)

adata = anndata.AnnData(
    X=X,
	obs=topic_obj.cell_data.loc[imputed_acc_obj.cell_names, :],
	var=topic_obj.region_data.loc[imputed_acc_obj.feature_names, :],
	dtype="int32"
)

adata.obs["cisTopic_nr_frag"] = adata.obs["cisTopic_nr_frag"].astype('int')
adata.obs["cisTopic_log_nr_frag"] = adata.obs["cisTopic_nr_frag"].astype('float')
adata.obs["cisTopic_nr_acc"] = adata.obs["cisTopic_nr_acc"].astype('int')
adata.obs["cisTopic_log_nr_acc"] = adata.obs["cisTopic_log_nr_acc"].astype('float')

del topic_obj
foo = gc.collect()

if args.norm:
	adata.raw = adata
	sc.pp.normalize_total(adata, target_sum=args.scale_factor)
	sc.pp.log1p(adata, base=2)


# save
adata.write_h5ad(save_path)
