#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=32
#PBS -N scenicplus
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
tf_file="${HOME}/project/multiomics/SCENIC/Data/avail_TF.txt"
biomart="http://sep2019.archive.ensembl.org"

# read mapping file
map_file="${HOME}/script/multiomics/map_group.csv"
celltype=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})

echo -e "Running for ${celltype}\n"

# Run
TMPDIR="/var/tmp"
tmp_dir=$(mktemp -dp $TMPDIR)

base_dir="${HOME}/project/multiomics/SCENIC/${celltype}"

script="${HOME}/project/multiomics/SCENIC/script/SCENICplus/run_SCENICplus.py"
python ${script} \
	-gex "${base_dir}/../full_gex_raw.h5ad" \
	-acc "${base_dir}/topic_obj.pkl" \
	-menr "${base_dir}/cisTarget/menr.pkl" \
	-g "cluster_label" \
	-tf $tf_file \
	-out "${base_dir}/SCENICplus" \
	-tmp $tmp_dir \
	-n $PBS_NP \
	-biomart $biomart


echo -e "\nthe script is over"
