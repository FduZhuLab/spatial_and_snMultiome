#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=2
#PBS -N overlap
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

# you'd better remove the list file if it is exist, but do not do it here because this script would be run for each region.
summit_list="${HOME}/project/multiomics/CallPeak/data/naiveSummit.tsv"
peak_list="${HOME}/project/multiomics/CallPeak/data/naivePeak.tsv"

echo -e "Running for ${region}\n"

# for each celltype
annot="${HOME}/project/multiomics/CallPeak/data/celltype.txt"
for clu in `cat $annot`; do

peak_dir="${HOME}/project/multiomics/CallPeak/peaks/split/${region}"
pooled_peak="${peak_dir}/Pool/${clu}_peaks.narrowPeak"
pooled_summit="${peak_dir}/Pool/${clu}_summits.bed"
rep1_peak="${peak_dir}/replicate/${clu}_No_53_peaks.narrowPeak"
rep2_peak="${peak_dir}/replicate/${clu}_No_61_peaks.narrowPeak"
pseudo1_peak="${peak_dir}/pseudorep/${clu}_pseudo00_peaks.narrowPeak"
pseudo2_peak="${peak_dir}/pseudorep/${clu}_pseudo01_peaks.narrowPeak"

# overlapped peak
overlap_dir="${HOME}/project/multiomics/CallPeak/peaks/overlap/${region}"
mkdir -p $overlap_dir
PooledInRep="${overlap_dir}/${clu}_pooled+rep.narrowPeak"
PooledInPseudo="${overlap_dir}/${clu}_pooled+pseudo.narrowPeak"

# parse peak summit
parse_dir="${HOME}/project/multiomics/CallPeak/peaks/parse/${region}"
mkdir -p $parse_dir
naive_peak="${parse_dir}/${clu}_narrowPeak"
naive_summit="${parse_dir}/${clu}_summit.bed"

if [ -f $pooled_peak ] && [ -f $pooled_summit ] && [ -f $rep1_peak ] && [ -f $rep2_peak ] && [ -f $pseudo1_peak ] && [ -f $pseudo2_peak ]
then
printf "Running for ${clu}..."

## overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5

## overlap with replicate
intersectBed -nonamecheck -wo -a ${pooled_peak} -b ${rep1_peak} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' \
| cut -f 1-10 | sort -k1,1 -k2,2n \
| intersectBed -nonamecheck -wo -a stdin -b ${rep2_peak} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' \
| cut -f 1-10 | sort -k1,1 -k2,2n > ${PooledInRep}


## overlap with pseudorep
intersectBed -nonamecheck -wo -a ${pooled_peak} -b ${pseudo1_peak} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' \
| cut -f 1-10 | sort -k1,1 -k2,2n \
| intersectBed -nonamecheck -wo -a stdin -b ${pseudo2_peak} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' \
| cut -f 1-10 | sort -k1,1 -k2,2n > ${PooledInPseudo}


## merge overlapped replicate and overlapped pseudorep
cat ${PooledInRep} ${PooledInPseudo} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){if ($5>1000) $5=1000; print $0}' \
| sort -k1,1 -k2,2n | uniq > ${naive_peak}


## get summit 
join -1 1 -2 4 <(cat ${naive_peak} | cut -f 4 | sort) <(sort -k4,4 ${pooled_summit}) -t$'\t' \
| awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$1,$5}' \
| sort -k1,1 -k2,2n | uniq > ${naive_summit}


## write file path
if [ -s $naive_peak ]; then
echo -e "${region}_${clu}\t${naive_peak}" >> $peak_list
echo -e "${region}_${clu}\t${naive_summit}" >> $summit_list
echo "done."
else
echo "no overlapped peaks, skip."
fi

else
echo "Skip for ${clu}. At least one required file is missing."
fi

done

echo -e "\nthe script is over"
