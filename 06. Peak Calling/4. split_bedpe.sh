#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=2
#PBS -N split_bed
#PBS -j oe
#PBS -o /home/whe/qsub_opt/${PBS_JOBID}.${PBS_JOBNAME}.log

# activate env
conda_env="py38_bio"
_CONDA_ROOT="${HOME}/Programs/miniconda3"

source ${_CONDA_ROOT}/bin/activate ${conda_env}

export OMP_NUM_THREADS=${PBS_NP}
export MKL_NUM_THREADS=${PBS_NP}
export OPENBLAS_NUM_THREADS=${PBS_NP}

# read variables for PBS jobs
var_list="${HOME}/project/multiomics/CallPeak/data/batch_list.csv"

donor=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${var_list})
region=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $2}' ${var_list})
batch=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $3}' ${var_list})

echo -e "Running for ${donor}, ${region}, ${batch}\n"

# Run
echo -e "the script is running\n"

input_file="${HOME}/project/multiomics/CallPeak/bedpe/${donor}_${region}_bedpe.gz"
annot_file="${HOME}/project/multiomics/CallPeak/cell_annot/split/${donor}_${region}_cell_annot.csv"
out_dir="${HOME}/project/multiomics/CallPeak/bedpe/split_byannot/${region}"

python "${HOME}/project/multiomics/CallPeak/script/split_bedpe.py" \
  -d $donor -b $batch \
  -in $input_file -ano $annot_file -out $out_dir \
  -celltype `cat ${HOME}/project/multiomics/CallPeak/data/celltype.txt`

echo -e "\nthe script is over"
