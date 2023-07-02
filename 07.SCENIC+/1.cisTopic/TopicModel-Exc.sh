#!/bin/bash
#PBS -q fat
#PBS -l walltime=72:00:00 -l nodes=1:ppn=30
#PBS -N TopicModel
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

echo -e "the script is running\n"

base_dir="${HOME}/project/multiomics/newSCENIC/Glutamatergic"
n_topics=$(echo {75..100..5} | sed 's/ /,/g')

script="${HOME}/project/multiomics/SCENIC/script/cisTopic/runMultiModels_lda_mallet.py"
python ${script} \
  -i "${base_dir}/topic_obj.pkl" \
  -o "${base_dir}/models.pkl" \
  -nt ${n_topics} \
  -c ${PBS_NP} \
  -it 300 \
  -s 42 \
  -sp "${base_dir}/Model/intermediate" \
  -td "${base_dir}/Model/tmp"


echo -e "\nthe script is over"
