#!/bin/bash
#PBS -q fat
#PBS -l walltime=72:00:00 -l nodes=1:ppn=30 -l mem=1000G
#PBS -N postanlys
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

# read mapping file
map_file="${HOME}/script/multiomics/map_group.csv"
group=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})

echo -e "Running for ${group}\n"

# Run
base_dir="${HOME}/project/multiomics/SCENIC"
home_dir="${base_dir}/${group}/SCENICplus"

script="${base_dir}/script/PostAnlys/post_SCENICplus.py"
python ${script} \
	-i $home_dir \
	-o "${home_dir}/PostAnlys" \
	-suf "" \
	-g "cluster_label" \
	-batch "Sample" \
	-n ${PBS_NP}


echo -e "\nthe script is over"
