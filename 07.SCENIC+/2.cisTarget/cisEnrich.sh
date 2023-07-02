#!/bin/bash
#PBS -q fat
#PBS -l walltime=72:00:00 -l nodes=1:ppn=16
#PBS -N cisTarget
#PBS -j oe
#PBS -o /home/whe/qsub_opt/${PBS_JOBID}.${PBS_JOBNAME}.log

conda_env="scenicplus"
_CONDA_ROOT="${HOME}/Programs/miniconda3"

source ${_CONDA_ROOT}/bin/activate ${conda_env}

export OMP_NUM_THREADS=${PBS_NP}
export MKL_NUM_THREADS=${PBS_NP}
export OPENBLAS_NUM_THREADS=${PBS_NP}
export VECLIB_MAXIMUM_THREADS=${PBS_NP}
export NUMEXPR_NUM_THREADS=${PBS_NP}

echo -e "the script is running\n"

# database
db_base="${HOME}/project/SCENIC/data/human"
rank_db="${db_base}/hg38/new/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
score_db="${db_base}/hg38/new/hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
motif_db="${db_base}/Motif/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
version="v10nr_clust"
biomart="http://sep2019.archive.ensembl.org"

# read mapping file
map_file="${HOME}/script/multiomics/map_group.csv"
group=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})

echo -e "Running for ${group}\n"

# Run
TMPDIR="/var/tmp"
tmp_dir=$(mktemp -dp $TMPDIR)
base_dir="${HOME}/project/multiomics/SCENIC/${group}"

script="${HOME}/project/multiomics/SCENIC/script/cisTarget/run_cisTarget.py"
python ${script} \
	-sets "${base_dir}/Evaluation/region_sets.pkl" \
	-r $rank_db -s $score_db -m $motif_db -v $version -biomart $biomart \
	-out "${base_dir}/cisTarget" \
	-tmp $tmp_dir \
	-n $PBS_NP


echo -e "\nthe script is over"
