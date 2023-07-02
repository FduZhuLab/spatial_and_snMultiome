#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=12
#PBS -N Intg-ssv4
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
map_file="${HOME}/script/multiomics/map_region+ref.csv"
region=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})
ref_region=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $2}' ${map_file})

echo -e "Running for ${region}, ${ref_region}\n"

# run
qry_base="${HOME}/project/multiomics/SplitObj/Seurat/perRegion"
ref_file="${HOME}/project/reference/Processed/Allen_SS_human_Multiple/SplitRegion/${ref_region}/cortex.obj.qs"
out_base="${HOME}/project/multiomics/Integration/Regional-Allen_SSv4/${region}"

script="${HOME}/R_func/seurat/Integrate/my.Integrate.R"
Rscript ${script} \
	-obj_list "${qry_base}/${region}/rna.intg.qs" $ref_file \
	-out $out_base \
	-anchor_type "hvg" \
	-ref_id 2 -use_ref


echo -e "\nthe script is over"
