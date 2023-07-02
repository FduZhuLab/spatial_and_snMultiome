#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=1
#PBS -N align_tnf5
#PBS -j oe
#PBS -o /home/whe/qsub_opt/${PBS_JOBID}.${PBS_JOBNAME}.log

# activate env
conda_env="py38_bio"
_CONDA_ROOT="${HOME}/Programs/miniconda3"

source ${_CONDA_ROOT}/bin/activate ${conda_env}

export OMP_NUM_THREADS=${PBS_NP}
export MKL_NUM_THREADS=${PBS_NP}
export OPENBLAS_NUM_THREADS=${PBS_NP}

# read mapping file
map_file="${HOME}/project/multiomics/CallPeak/data/map_replicate.csv"
region=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})
rep=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $3}' ${map_file})
echo -e "Running for ${region}, ${rep}\n"

# Run
echo -e "the script is running\n"

python "${HOME}/project/multiomics/CallPeak/script/align_tnf5.py" \
  -in "${HOME}/project/multiomics/CallPeak/bedpe/split_byannot/${region}" \
  -out "project/multiomics/CallPeak/tnf5_bed/rep/${region}" \
  -rep $rep

echo -e "\nthe script is over"
