#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=16
#PBS -N cci
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

label="cluster_label"
input="${HOME}/project/multiomics/SplitObj/Seurat/perRegion/${region}"
output="${HOME}/project/multiomics/CellChat/${label}/${region}"

customDB="${HOME}/project/multiomics/CellChat/Data/liana_curated_ccdb.rds"

script="${HOME}/project/multiomics/CellChat/script/run_cellchat.R"
Rscript ${script} \
	-i "${input}/rna.intg.qs" \
	-o ${output} \
	-db $customDB \
	-grp $label


echo -e "\nthe script is over"
