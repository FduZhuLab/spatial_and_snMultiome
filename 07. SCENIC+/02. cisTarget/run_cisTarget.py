import os
import pickle
import argparse

from scenicplus.wrappers.run_pycistarget import run_pycistarget


# arguments
parser = argparse.ArgumentParser(description="Run cisTarget enrichment.")
parser.add_argument(
	'--region_set_path', '-sets', type=str, required=True
)
parser.add_argument(
	'--rankings_db', '-r', type=str, required=True
)
parser.add_argument(
	'--scores_db', '-s', type=str, required=True
)
parser.add_argument(
	'--motif_annot', '-m', type=str, required=True
)
parser.add_argument(
	'--verion', '-v', type=str, required=True, help="motif version, like 'v10nr_clust'"
)
parser.add_argument(
	'--outdir', '-out', type=str, required=True
)
parser.add_argument(
	'--ncpu', '-n', type=int, required=True
)
parser.add_argument(
	'--tmpdir', '-tmp', type=str, required=False, default=None
)
parser.add_argument(
	'--species', '-species', type=str, required=False, default="homo_sapiens"
)
parser.add_argument(
	'--biomart', '-biomart', type=str, required=False, default="http://www.ensembl.org"
)
parser.add_argument(
	'--no_promoter', '-nopro', default=False, action='store_true',
)
parser.add_argument(
	'--with_pro', '-with_pro', dest='no_promoter', action='store_false',
)

args = parser.parse_args()

# Load args
with open(args.region_set_path, "rb") as f:
	region_sets = pickle.load(f)

outdir = args.outdir
os.makedirs(outdir, exist_ok=True)

tmpdir = args.tmpdir
if tmpdir is not None:
	os.makedirs(tmpdir, exist_ok=True)


# Run
run_pycistarget(
    region_sets = region_sets,
	ctx_db_path = args.rankings_db,
    dem_db_path = args.scores_db,
    path_to_motif_annotations = args.motif_annot,
	annotation_version = args.verion,
    save_path = outdir,
	_temp_dir = tmpdir,
	save_partial = True,
    species = args.species,
    run_without_promoters = args.no_promoter,
	biomart_host=args.biomart,
    n_cpu = args.ncpu
)

