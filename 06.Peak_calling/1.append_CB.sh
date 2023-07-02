#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=8
#PBS -N AppendCB
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

bam_file="${HOME}/project/multiomics/cellranger_out/${donor}/${region}/atac_possorted_bam.bam"
out_dir="${HOME}/project/multiomics/CallPeak/filter_bam/${donor}/${region}"

mkdir -p $out_dir && cd $out_dir

### extract the header file
samtools view ${bam_file} -H > possorted.header.sam

### create a bam file with the barcode embedded into the read name
cat <( cat possorted.header.sam ) \
<( samtools view -@ ${PBS_NP} -u ${bam_file} | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) | samtools view -@ ${PBS_NP} -bS - -o ${out_dir}/possorted.CB.bam

echo -e "\nthe script is over"
