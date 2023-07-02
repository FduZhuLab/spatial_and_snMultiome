# Load packages
import os
import gc
import argparse

import numpy as np
import pandas as pd

import scanpy as sc

from cell2location.models import RegressionModel
from cell2location.models import Cell2location

os.chdir(os.getenv("HOME"))

# Argements
parser = argparse.ArgumentParser()
parser.add_argument('--donor', '-d', type=str, required=True)
parser.add_argument('--region', '-r', type=str, required=True)
parser.add_argument('--spt_input', '-spt', type=str, required=True)
parser.add_argument('--out_dir', '-out', type=str, required=True)
parser.add_argument('--group_by', '-grp', type=str, required=True)

parser.add_argument('--Ncell', '-N', type=int, required=True)
parser.add_argument('--alpha', '-a', type=float, default=20)

parser.add_argument('--ref_epochs', '-ep1', type=int, default=1000)
parser.add_argument('--map_epochs', '-ep2', type=int, default=15000)

args = parser.parse_args()

donor = args.donor
region = args.region
group_by = args.group_by

# Load Data
spt_ada = sc.read(args.spt_input)
ind = ~(
	spt_ada.obs["is_outlier"] | (spt_ada.obs["Spatial_Layer"] == "_Outlier")
)
spt_ada = spt_ada[ind, :].copy()
spt_ada.X = spt_ada.layers["counts"].astype("int32")

sc_ada = sc.read("project/multiomics/GeneAnlys/Merge/merge_rna.h5ad")
sc_ada.obs[group_by] = sc_ada.obs[group_by].astype("str")
ind = (
	(sc_ada.obs["Region"] == region) &
	(sc_ada.obs["Sample"] == donor)
)
sc_ada = sc_ada[ind, :].copy()
sc_ada.X = sc_ada.layers["counts"].astype("int32")

gc.collect()

print('\nData was loaded')

out_dir = args.out_dir
os.makedirs(out_dir, exist_ok=True)

# Prepare data

## avail clusters
grp_num = sc_ada.obs[group_by].value_counts()
groups = grp_num.index[grp_num >= 3].sort_values()

all_cluster_metric = pd.read_csv(
	"project/multiomics/SplitObj/Seurat/eval_cluster/all_metric_df.csv"
)
ind = (
	(all_cluster_metric["region"] == region) &
	(all_cluster_metric["donor"] == donor) &
	(all_cluster_metric["is_validate_f05"] == 1)
)
valid_grp = all_cluster_metric.loc[ind, "group"].unique()

groups = groups[groups.isin(valid_grp)]

print(f'The trained {group_by}: {np.asarray(groups)}')

## trained genes
deg_dir = "project/multiomics/GeneAnlys/Merge/celltype_marker/cluster_in_batch/top200"
mdf = pd.read_csv(f"{deg_dir}/{donor}_{region}_cluster_label_top200_mdf.csv")
# mdf = mdf[mdf['group'].isin(groups)]

train_genes = mdf['feature']
train_genes = train_genes[train_genes.isin(spt_ada.var_names)].unique()

print(f'The length of trained genes is: {len(train_genes)}')

## prepare data
sc_ada = sc_ada[
	sc_ada.obs[group_by].isin(groups), train_genes
].copy()
sc_ada.obs[group_by] = sc_ada.obs[group_by].astype("str")

spt_ada = spt_ada[:, train_genes].copy()

gc.collect()

print(f"\nThis is scRNA data: {sc_ada}")
print(f"\nThis is spatial data: {spt_ada}")

# Reference signatures
RegressionModel.setup_anndata(adata=sc_ada, labels_key=group_by)
ref_model = RegressionModel(sc_ada)
ref_model.train(max_epochs=args.ref_epochs, use_gpu=True)

sc_ada = ref_model.export_posterior(
	sc_ada, sample_kwargs={'use_gpu': True}
)

## save
sc_ada.write(f"{out_dir}/sc.h5ad")

inf_aver = sc_ada.varm['means_per_cluster_mu_fg'].copy()
inf_aver.columns = inf_aver.columns.str.replace("means_per_cluster_mu_fg_", "")

inf_aver.to_csv(f"{out_dir}/ref_profile.csv")

ref_model.save(f"{out_dir}/ref_signatures", overwrite=True)

print("\nSaved reference model")

# Mapping to Spatial
print(f'\nThe N_cells_per_location is {args.Ncell}')
print(f'\nThe detection_alpha is {args.alpha}')
gc.collect()

Cell2location.setup_anndata(adata=spt_ada)
map_model = Cell2location(
	spt_ada, cell_state_df=inf_aver,
	N_cells_per_location=args.Ncell,
	detection_alpha=args.alpha
)
map_model.train(max_epochs=args.map_epochs, batch_size=None, use_gpu=True)

spt_ada = map_model.export_posterior(
	spt_ada, sample_kwargs={'use_gpu': True, 'batch_size': map_model.adata.n_obs}
)

## save
spt_ada.write(f"{out_dir}/spt.h5ad")

map_model.save(f"{out_dir}/spt_map", overwrite=True)

print("\nSaved mapping model")

# Save abundance
map_annot = spt_ada.obsm['q05_cell_abundance_w_sf'].copy()
map_annot.columns = map_annot.columns.str.replace("q05cell_abundance_w_sf_", "")

map_annot = spt_ada.obs[["Spatial_Layer"]].join(map_annot)
map_annot.index.name = "spot"

map_annot.to_csv(f"{out_dir}/map_annot.csv")

print("\nSaved mapping annotation abundance")

