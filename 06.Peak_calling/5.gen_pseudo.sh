#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=1
#PBS -N pseudo
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
map_file="${HOME}/project/multiomics/CallPeak/data/map_region.csv"
region=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})
echo -e "Running for ${region}\n"

input_dir="${HOME}/project/multiomics/CallPeak/bedpe/split_byannot/${region}"
out_dir=$input_dir
mkdir -p $out_dir && cd $out_dir

# define random seed function
my_get_seed()
{
  seed="$1";
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null;
}

# for each celltype
annot="${HOME}/project/multiomics/CallPeak/data/celltype.txt"
for clu in `cat $annot`; do

## shuffle and split equally
rep1="${input_dir}/${clu}_No_53.bedpe"
rep2="${input_dir}/${clu}_No_61.bedpe"

if [ -f $rep1 ] && [ -f $rep2 ]; then
    printf "Running for ${clu}..."
    nlines=$(cat $rep1 $rep2 | wc -l)
    nlines=$(( (nlines + 1) / 2 ))
    cat $rep1 $rep2 | shuf --random-source=<(my_get_seed 42) | split -d -l $nlines - "${clu}_pseudo"
    echo "done"
else
    echo "For ${clu}, at least one file is missing."
fi

done

# rename files with .bedpe suffix
find . -type f -name "*_pseudo*" -exec mv {} {}.bedpe \;

echo -e "\nthe script is over"
