#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=1
#PBS -N bam2bedpe
#PBS -j oe
#PBS -o /home/whe/qsub_opt/${PBS_JOBID}.${PBS_JOBNAME}.log

# activate env
conda_env="bulk-seq"
_CONDA_ROOT="${HOME}/Programs/miniconda3"

source ${_CONDA_ROOT}/bin/activate ${conda_env}

export OMP_NUM_THREADS=${PBS_NP}
export MKL_NUM_THREADS=${PBS_NP}
export OPENBLAS_NUM_THREADS=${PBS_NP}

# read mapping file
map_file="${HOME}/script/multiomics/map_batch.csv"

donor=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})
region=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $2}' ${map_file})

echo -e "Running for ${donor}, ${region}\n"

# Run
echo -e "the script is running\n"

out_dir="${HOME}/project/multiomics/CallPeak/filter_bam/${donor}/${region}"
bam_file="${out_dir}/filtered_dedup_nsrt.bam"
bedpe_file="${HOME}/project/multiomics/CallPeak/bedpe/${donor}_${region}_bedpe.gz"

mkdir -p $out_dir && cd $out_dir

bedtools bamtobed -bedpe -mate1 -i ${bam_file} | gzip > ${bedpe_file}

echo -e "\nthe script is over"
