#!/bin/bash
#SBATCH --partition=organ
#SBATCH --nodelist=gpu1
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --signal=USR2
#SBATCH --job-name=Seg_Calc
#SBATCH --output=/home/%u/slurm_log/%j.%x.out

map_file="${HOME}/script/map_batch.csv"

donor=$(awk -F',' -v x=${SLURM_ARRAY_TASK_ID} 'NR==x {print $1}' ${map_file})
region=$(awk -F',' -v x=${SLURM_ARRAY_TASK_ID} 'NR==x {print $2}' ${map_file})

conda_env="squid"
_CONDA_ROOT="${HOME}/programs/miniconda3"

source ${_CONDA_ROOT}/bin/activate ${conda_env}

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export VECLIB_MAXIMUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo -e "Running for ${donor}, ${region}\n"

img_base="${HOME}/project/Spatial/Squid/Segement/output"

python "${HOME}/project/Spatial/Squid/Segement/pipline_segment.py" \
  -d $donor -r $region -spt "${HOME}/project/Spatial/Object/Scanpy" \
  -img "${HOME}/project/Spatial/Images" \
  -o $img_base

echo -e "the Segment is over\n"

python "${HOME}/project/Spatial/Squid/Segement/pipline_calc_seg.py" \
  -d $donor -r $region -spt "${HOME}/project/Spatial/Object/Scanpy" \
  -img $img_base

echo -e "\nthe script is over"
