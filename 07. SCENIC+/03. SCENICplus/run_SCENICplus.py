import os
import pickle
import argparse

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import scanpy as sc
import dill
from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.wrappers.run_scenicplus import run_scenicplus


# arguments
parser = argparse.ArgumentParser(description="Run SCENICplus.")
parser.add_argument(
	'--gex_adata', '-gex', type=str, required=True
)
parser.add_argument(
	'--cis_topic', '-acc', type=str, required=True
)
parser.add_argument(
	'--menr', '-menr', type=str, required=True
)
parser.add_argument(
	'--object', '-obj', type=str, required=False, default=None
)
parser.add_argument(
	'--groups', '-g', type=str, required=True, nargs='+'
)
parser.add_argument(
	'--tf_file', '-tf', type=str, required=True
)
parser.add_argument(
	'--outdir', '-out', type=str, required=True
)
parser.add_argument(
	'--ncpu', '-n', type=int, required=True
)
parser.add_argument(
	'--tooldir', '-tooldir', type=str, required=False, default="/home/whe/Programs/miniconda3/envs/scenicplus/bin"
)
parser.add_argument(
	'--tmpdir', '-tmp', type=str, required=False, default=None
)
parser.add_argument(
	'--species', '-species', type=str, required=False, default="hsapiens"
)
parser.add_argument(
	'--assembly', '-assembly', type=str, required=False, default="hg38"
)
parser.add_argument(
	'--biomart', '-biomart', type=str, required=False, default="http://www.ensembl.org"
)

args = parser.parse_args()


# load object
if args.object is not None:
	with open(args.object, "rb") as f:
		scplus_obj = dill.load(f)

else:
	gex_adata = sc.read(args.gex_adata)

	with open(args.cis_topic, "rb") as f:
		topic_obj = pickle.load(f)
	
	n = topic_obj.selected_model.cell_topic.shape[0]
	print(f"The topic object has {n} topics")

	gex_adata = gex_adata[topic_obj.cell_names, ].copy()

	with open(args.menr, "rb") as f:
		menr = dill.load(f)

	scplus_obj = create_SCENICPLUS_object(
		GEX_anndata = gex_adata,
		cisTopic_obj = topic_obj,
		menr = menr
	)


# test TF
tf_list = pd.read_csv(args.tf_file, header=None).iloc[:, 0]
if not tf_list.isin(scplus_obj.metadata_genes.index).all():
	raise(ValueError("Some TFs are not in the gene list")) 


# test input variables
groups = args.groups
if type(groups) == list:
	print("The variables are:", groups)
	foo = scplus_obj.metadata_cell.loc[:, groups]
else:
	raise(ValueError("The groups argument must be a list"))

# test reduction
for k in scplus_obj.dr_cell:
	reduct = scplus_obj.dr_cell[k]
	if reduct.shape[1] > 2:
		print("Modify reduction:", k)
		reduct = reduct.iloc[:, 0:2]
		scplus_obj.dr_cell[k] = reduct

# to dense array
if type(scplus_obj.X_EXP) == csr_matrix:
	print("Making sparse array to dense array:")
	scplus_obj.X_EXP = scplus_obj.X_EXP.toarray()

# print
print("This is the object:")
print(scplus_obj)

# dir
outdir = args.outdir
os.makedirs(outdir, exist_ok=True)

tmpdir = args.tmpdir
if tmpdir is not None:
	os.makedirs(tmpdir, exist_ok=True)


print("The out directory is:")
print(outdir)
print("The ray temp directory is:")
print(tmpdir)

# Run
print("\nStart running:")

try:
	run_scenicplus(
		scplus_obj = scplus_obj,
		variable = groups,
		tf_file = args.tf_file,
		save_path = outdir,
		upstream = [1000, 150000],
		downstream = [1000, 150000],
		calculate_TF_eGRN_correlation = True,
		calculate_DEGs_DARs = False,
		export_to_loom_file = True,
		export_to_UCSC_file = False,
		species = args.species,
		assembly = args.assembly,
		biomart_host = args.biomart,
		path_bedToBigBed = args.tooldir,
		n_cpu = args.ncpu,
		_temp_dir = tmpdir
	)
except Exception as e:
	t = type(scplus_obj.X_EXP)
	if (t == np.ndarray) | (t == np.matrix):
		scplus_obj.X_EXP = csr_matrix(scplus_obj.X_EXP)

	#in case of failure, still save the object
	with open(os.path.join(outdir, 'scplus_obj.pkl'), 'wb') as f:
		dill.dump(scplus_obj, f, protocol=-1)

	raise(e)
