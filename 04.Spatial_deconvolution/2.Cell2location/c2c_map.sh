#!/bin/bash
#PBS -q gpu
#PBS -l walltime=200:00:00 -l nodes=gpu03:ppn=8 -l mem=50G
#PBS -N C2C
#PBS -j oe
#PBS -o /home/whe/qsub_opt/${PBS_JOBID}.${PBS_JOBNAME}.log

# activate env
conda_env="cell2location"
_CONDA_ROOT="${HOME}/Programs/miniconda3"

source ${_CONDA_ROOT}/bin/activate ${conda_env}

export OMP_NUM_THREADS=${PBS_NP}
export MKL_NUM_THREADS=${PBS_NP}
export OPENBLAS_NUM_THREADS=${PBS_NP}
export VECLIB_MAXIMUM_THREADS=${PBS_NP}
export NUMEXPR_NUM_THREADS=${PBS_NP}

# read map file
map_file="${HOME}/project/Spatial/Cell2loc/Nucleus_num.csv"

donor=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})
region=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $2}' ${map_file})
NN=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $5}' ${map_file})

echo -e "Running for ${donor}, ${region}\n"

# run
out_base="${HOME}/project/Spatial/Cell2loc"
spt_base="${HOME}/project/Spatial/Object/Scanpy"

python -u "${out_base}/script/cell2loc_pipeline.py" \
  -d $donor -r $region \
  -spt "${spt_base}/${donor}/${region}/spt_adata.h5ad" \
  -out "${out_base}/output-3/${donor}/${region}" \
  -N $NN \
  -grp "cluster_label"


echo -e "\nthe script is over"
