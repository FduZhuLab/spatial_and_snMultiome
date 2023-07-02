#!/bin/bash
#PBS -q fat
#PBS -l walltime=72:00:00 -l nodes=1:ppn=16
#PBS -N Intg-spe-EXC
#PBS -j oe
#PBS -o /home/whe/qsub_opt/${PBS_JOBID}.${PBS_JOBNAME}.log

# activate env
conda_env="r4_bio"
_CONDA_ROOT="${HOME}/Programs/miniconda3"

source ${_CONDA_ROOT}/bin/activate ${conda_env}

export OMP_NUM_THREADS=${PBS_NP}
export MKL_NUM_THREADS=${PBS_NP}
export OPENBLAS_NUM_THREADS=${PBS_NP}

celltype="Glutamatergic"
echo -e "Running for ${celltype}\n"

# run
hs_base="${HOME}/project/multiomics/GeneAnlys"
hs_file="${hs_base}/${celltype}/merge_rna_intg.qs"

mm_base="${HOME}/project/reference/Processed/Allen_Mouse_CTX_HPF/downspl"
mm_file="${mm_base}/spl.mouse_${celltype}.qs"

mq_base="${HOME}/project/reference/Processed/Macaque_cortex"
mq_file="${mq_base}/${celltype}/spl.mq.obj.qs"

out_base="${HOME}/project/multiomics/Integration/Full_Intg_mm_mq"
out_dir="${out_base}/${celltype}"

script="${HOME}/R_func/seurat/Integrate/my.Integrate.R"
Rscript ${script} \
	-obj_list $hs_file $mq_file $mm_file \
	-out $out_dir \
	-anchor_type "hvg" \
	-ref_id 3 \
	-do_sct \
	-num_hvg 5000 \
	-new_select \
	-tf_reduc 'cca'


echo -e "\nthe script is over"
