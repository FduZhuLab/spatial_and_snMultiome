#!/bin/bash
#PBS -q batch
#PBS -l walltime=72:00:00 -l nodes=1:ppn=8
#PBS -N MergePeak
#PBS -j oe
#PBS -o /home/whe/qsub_opt/${PBS_JOBID}.${PBS_JOBNAME}.log

# activate env
conda_env="r4_bio"
_CONDA_ROOT="${HOME}/Programs/miniconda3"

source ${_CONDA_ROOT}/bin/activate ${conda_env}

export OMP_NUM_THREADS=${PBS_NP}
export MKL_NUM_THREADS=${PBS_NP}
export OPENBLAS_NUM_THREADS=${PBS_NP}

# Run
blacklist="${HOME}/project/multiomics/CallPeak/data/ENCFF356LFX_GRCh38_unified_blacklist.bed.gz"
chromSize="${HOME}/project/multiomics/CallPeak/data/chrom_size.txt"

input_list="${HOME}/project/multiomics/CallPeak/data/naiveSummit.tsv"
out_dir="${HOME}/project/multiomics/CallPeak/peaks/merge"
mkdir -p $out_dir

echo -e "the script is running\n"

Rscript "${HOME}/project/multiomics/CallPeak/script/merge_peak.R" \
  -i $input_list \
  -g hg38 \
  --blacklist $blacklist \
  --chromSize $chromSize \
  -o $out_dir \
  -t ${PBS_NP}


## filter by "Score Per Million" (SPM) >= 5 
sed '1d' "${out_dir}/filteredNfixed.union.peakSet" \
| awk 'BEGIN{FS=OFS="\t"}($11>=5){print $1,$2,$3,$7,$6}' \
| sort -k1,1 -k2,2n | uniq > "${out_dir}/filteredNfixed.union.bed"

awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"peak"NR,$5}' "${out_dir}/filteredNfixed.union.bed" > "${out_dir}/final.union.bed"

echo -e "\nthe script is over"
