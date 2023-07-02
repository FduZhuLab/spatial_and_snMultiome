#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=1
#PBS -N pooled
#PBS -j oe
#PBS -o /home/whe/qsub_opt/${PBS_JOBID}.${PBS_JOBNAME}.log

# read mapping file
map_file="${HOME}/project/multiomics/CallPeak/data/map_region.csv"
region=$(awk -F',' -v x=${PBS_ARRAYID} 'NR==x {print $1}' ${map_file})
echo -e "Running for ${region}\n"

in_dir="${HOME}/project/multiomics/CallPeak/tnf5_bed/rep/${region}"
out_dir="${HOME}/project/multiomics/CallPeak/tnf5_bed/Pool/${region}"
mkdir -p $out_dir

## for each celltype
annot="${HOME}/project/multiomics/CallPeak/data/celltype.txt"
for clu in `cat $annot`; do

tnf5_no53="${in_dir}/${clu}_No_53.tnf5.bed"
tnf5_no61="${in_dir}/${clu}_No_61.tnf5.bed"
pool_tnf5="${out_dir}/${clu}.tnf5.bed"

if [ -f $tnf5_no53 ] && [ -f $tnf5_no61 ]; then
    printf "Running for ${clu}..."
    cat $tnf5_no53 $tnf5_no61 > $pool_tnf5
    echo "done"
else
	echo "For ${clu}, at least one file is missing."
fi

done

# over
echo -e "\nthe script is over"
