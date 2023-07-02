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
celltype="GABAergic"

echo -e "Running for ${celltype}\n"

# run
qry_file="${HOME}/project/multiomics/GeneAnlys/${celltype}/merge_rna_intg.qs"
ref_file="/home/whe/project/reference/Processed/Allen_SS_human_Multiple/subclass/${celltype}/human_obj.qs"
out_base="${HOME}/project/multiomics/Integration/Allen_SSv4/${celltype}"

script="${HOME}/R_func/seurat/Integrate/my.Integrate.R"
Rscript ${script} \
	-obj_list $qry_file $ref_file \
	-out $out_base \
	-anchor_type "hvg" \
	-ref_id 2


echo -e "\nthe script is over"
