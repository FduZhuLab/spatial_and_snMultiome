#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=4
#PBS -N RenameFrag
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
map_file="${HOME}/project/multiomics/CallPeak/data/map_batch.csv"

donor=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})
region=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $2}' ${map_file})
batch=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $3}' ${map_file})

echo -e "Running for ${donor}, ${region}, ${batch}\n"

# Run
echo -e "the script is running\n"
input="${HOME}/project/multiomics/cellranger_out/${donor}/${region}/atac_fragments.tsv.gz"
out_dir="${HOME}/project/multiomics/RenameFrag/${donor}/${region}"
mkdir -p $out_dir && cd $out_dir

zcat $input \
| awk -v bb="${batch}_" 'BEGIN{FS=OFS="\t"} {if ($1~/^#/) {print $0} else {$4=bb$4; print $0}}' \
| bgzip -@ $PBS_NP -c > atac_fragments.tsv.gz && tabix -p bed atac_fragments.tsv.gz

echo -e "\nthe script is over"
