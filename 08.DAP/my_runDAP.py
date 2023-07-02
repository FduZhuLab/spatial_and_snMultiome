import os
import pickle
from warnings import warn
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from typing import List, Optional


def my_runDAP(
	adata: anndata.AnnData,
	group_by: str,
	selected_group: str,
	min_cells: Optional[int] = 3,
	method: Optional[str] = "wilcoxon",
	pval_cutoff: Optional[float] = 0.05,
	pct_cutoff: Optional[float] = 0.1,
	log2fc_cutoff: Optional[float] = 0.,
	save_dir: Optional[str] = None
):
	grp_num = adata.obs.groupby(group_by).size()
	avail_grp = grp_num.index[grp_num >= min_cells]
	
	if len(avail_grp) < 2:
		warn(f"Few available groups with min cells {min_cells}")
		return None

	if selected_group not in avail_grp:
		warn(f"{selected_group} is not in the available groups")
		return None

	vs_groups = avail_grp[avail_grp != selected_group]
	
	dap_results = {"pos": {}, "neg": {}}

	for vsg in vs_groups:
		k = f"{group_by}_{selected_group}"
		sc.tl.rank_genes_groups(
			adata, groupby=group_by,
			groups=[selected_group], reference=vsg,
			method=method, pts=True, key_added=k
		)
		dap_df = sc.get.rank_genes_groups_df(adata, group=selected_group, key=k)

		dap_df = dap_df.rename(columns={"names": "gene", "pct_nz_group": "pct.1"})
		dap_df["pct.2"] = dap_df["gene"].map(adata.uns[k]['pts'][vsg].to_dict())
		pct_max = dap_df[["pct.1", "pct.2"]].max(axis=1)
		dap_df["pct.qdiff"] = (dap_df["pct.1"] - dap_df["pct.2"]).abs() / pct_max
		dap_df.set_index("gene", inplace=True)

		dap_pos = dap_df[dap_df["logfoldchanges"] > 0]
		dap_neg = dap_df[dap_df["logfoldchanges"] < 0]
		if dap_pos.shape[0] > 0:
			dap_results["pos"][vsg] = dap_pos
		if dap_neg.shape[0] > 0:
			dap_results["neg"][vsg] = dap_neg


	# save
	if save_dir is not None:
		os.makedirs(save_dir, exist_ok=True)
		with open(f"{save_dir}/dap_results.pkl", "wb") as f:
			pickle.dump(dap_results, f)


	# intersection
	codap_results = dict()
	for s in ["pos", "neg"]:
		dap_results_sign = dap_results[s]
		dap_list = [dap_results_sign[k] for k in dap_results_sign]
		if len(dap_list) == 0:
			continue

		co_dap = pd.Index(set.intersection(*[set(df.index) for df in dap_list]))
		if len(co_dap) == 0:
			continue

		for i, df in enumerate(dap_list):
			df = df.loc[co_dap, :]
			for k in ["scores", "logfoldchanges"]:
				df[k] = df[k].abs()

			dap_list[i] = df

		codap_df = {}
		for k in dap_list[0].columns:
			if k != "pct.1":
				X = np.array([df[k].to_numpy() for df in dap_list])
				codap_df[f"max_{k}"] = X.max(axis=0)
				codap_df[f"min_{k}"] = X.min(axis=0)

		codap_df = pd.DataFrame(codap_df)
		codap_df["pct.1"] = dap_list[0]["pct.1"].to_numpy()
		codap_df.index = co_dap
		codap_df.sort_values(by="min_scores", ascending=False, inplace=True)

		codap_df.index.name = "gene"
		codap_df.reset_index(inplace=True)
		codap_df["region"] = selected_group

		# filter
		ind = (
			(codap_df["min_logfoldchanges"] > log2fc_cutoff) &
			(codap_df["max_pvals_adj"] < pval_cutoff) &
			((codap_df["min_pct.2"] > pct_cutoff) | (codap_df["pct.1"] > pct_cutoff))
		)
		codap_df = codap_df[ind]

		if codap_df.shape[0] > 0:
			codap_results[s] = codap_df
		else:
			warn(f"No conserved DAPs for {s}")


	if len(codap_results) > 0:
		return codap_results
	else:
		warn("No conserved DAPs")
		return None
