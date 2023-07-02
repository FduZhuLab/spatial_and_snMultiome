#!/bin/bash
#PBS -q batch
#PBS -l walltime=24:00:00 -l nodes=1:ppn=32 -l mem=200G
#PBS -N preprocess_each
#PBS -j oe
#PBS -o /home/whe/qsub_opt/${PBS_JOBID}.${PBS_JOBNAME}.log

script="/home/whe/script/multiomics/withPars/preprocess.R"
conda_env="r4_bio"
_CONDA_ROOT="/home/whe/Programs/miniconda3"

source ${_CONDA_ROOT}/bin/activate ${conda_env}

export OMP_NUM_THREADS=${PBS_NP}
export OPENBLAS_NUM_THREADS=${PBS_NP}

echo -e 'the script is running\n'

Rscript ${script} --iter --sample="No_53,No_61" --region="M1C,V1C,ACC,dlPFC"

echo -e '\nthe script is over'
