#!/bin/bash
#PBS -q fat
#PBS -l walltime=72:00:00 -l nodes=1:ppn=16 -l mem=100G
#PBS -N RegDAP
#PBS -j oe
#PBS -o /home/whe/qsub_opt/${PBS_JOBID}.${PBS_JOBNAME}.log

# activate env
conda_env="py38_bio"
_CONDA_ROOT="${HOME}/Programs/miniconda3"

source ${_CONDA_ROOT}/bin/activate ${conda_env}

export OMP_NUM_THREADS=${PBS_NP}
export MKL_NUM_THREADS=${PBS_NP}
export OPENBLAS_NUM_THREADS=${PBS_NP}
export VECLIB_MAXIMUM_THREADS=${PBS_NP}
export NUMEXPR_NUM_THREADS=${PBS_NP}

# read mapping file
map_file="${HOME}/script/multiomics/map_group+region.csv"
group=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})
region=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $2}' ${map_file})

echo -e "Running for ${group}, ${region}\n"

# run
base_dir="${HOME}/project/multiomics/SCENIC"

script="${base_dir}/script/Regional_DAP.py"
python ${script} \
	-adata "${base_dir}/${group}/imputed_ACC.h5ad" \
	-out "${base_dir}/Data/Regional_DA/${group}/${region}" \
	-cl "subclass_label" \
	-rl "Region" -r $region \
	-norm


echo -e "\nthe script is over"
