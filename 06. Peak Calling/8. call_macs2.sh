#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=2
#PBS -N macs2
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
map_file="${HOME}/project/multiomics/CallPeak/data/map_replicate.csv"
region=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})
subdir=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $2}' ${map_file})
rep=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $3}' ${map_file})
echo -e "Running for ${region}, ${rep}\n"

# Run
input_dir="${HOME}/project/multiomics/CallPeak/tnf5_bed/rep/${region}"
out_dir="${HOME}/project/multiomics/CallPeak/peaks/split/${region}/${subdir}"
mkdir -p $out_dir

## for each celltype
annot="${HOME}/project/multiomics/CallPeak/data/celltype.txt"
for clu in `cat $annot`; do

tnf5_bed="${input_dir}/${clu}_${rep}.tnf5.bed"
out_name="${clu}_${rep}"

if [ -f $tnf5_bed ]; then
    printf "Running for ${clu}, output name contains ${out_name}..."
    macs2 callpeak -t $tnf5_bed -f BED -n $out_name -g hs -q 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits --outdir $out_dir --seed 42
    echo "done"
else
    echo "For ${clu}, the file is missing."
fi

done

echo -e "\nthe script is over"
