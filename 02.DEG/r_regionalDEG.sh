#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=16
#PBS -N RegionalDEG
#PBS -j oe
#PBS -o /home/whe/qsub_opt/${PBS_JOBID}.${PBS_JOBNAME}.log

# activate env
conda_env="r4_bio"
_CONDA_ROOT="${HOME}/Programs/miniconda3"

source ${_CONDA_ROOT}/bin/activate ${conda_env}

export OMP_NUM_THREADS=${PBS_NP}
export MKL_NUM_THREADS=${PBS_NP}
export OPENBLAS_NUM_THREADS=${PBS_NP}

# read mapping file
map_file="${HOME}/script/multiomics/map_region.csv"
region=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})
echo -e "Running for ${region}\n"

# run
echo -e "the script is running\n"

base_dir="${HOME}/project/multiomics/GeneAnlys"

script="${HOME}/project/multiomics/GeneAnlys/script/Regional_DEG.R"
Rscript ${script} \
  -r $region \
  -i "${base_dir}/Merge/merge_rna_intg.qs" \
  -o "${base_dir}/Regional_DEG/${region}" \
  -grp "subclass_label" -rgrp "Region"


echo -e "\nthe script is over"
