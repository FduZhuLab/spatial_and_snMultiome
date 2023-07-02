#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=8
#PBS -N filter_bam
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
bam_in="${out_dir}/possorted.CB.bam"
bam_out="${out_dir}/filtered_dedup_nsrt.bam"

cd $out_dir

samtools view -@ ${PBS_NP} -bu -f 2 -F 0x70c -q 30 ${bam_in} | samtools collate - -@ ${PBS_NP} -u | samtools fixmate -@ ${PBS_NP} -rmu - - | samtools view -@ ${PBS_NP} -bu -f 2 -F 1804 - | samtools sort - -@ ${PBS_NP} -u | samtools markdup -ru -@ ${PBS_NP} - - | samtools sort - -@ ${PBS_NP} -n -o ${bam_out}

rm -f $bam_in # to save disk space

echo -e "\nthe script is over"
