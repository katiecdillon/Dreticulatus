#!/bin/bash
#SBATCH --job-name=align
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=900gb
#SBATCH --export=NONE
#SBATCH --time=5-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=kcd88651@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS
#SBATCH --array=1-8

# ---------------------------------------------------------------------------- #
# CREATE CONDA ENVIRONMENTS
# ---------------------------------------------------------------------------- #
source ~/.bashrc

conda create -n align_env -c bioconda minimap2=2.29 samtools=1.21 -y

# ---------------------------------------------------------------------------- #
# FILE PATHS
# ---------------------------------------------------------------------------- #
config=/scratch/kcd88651/ticks/D_reticulatus/STEP03_config_align.txt

dir=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $2}' "$config")
assembly=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $3}' "$config")
raw=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $4}' "$config")
header=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $5}' "$config")
ref=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $6}' "$config")

# ---------------------------------------------------------------------------- #
# ALIGN, CONVERT TO BAM, SORT BAM, AND INDEX BAM
# ---------------------------------------------------------------------------- #
conda activate align_env

cd "$dir"

minimap2 -t 32 -ax map-ont "$assembly" "$raw" | \
    samtools view -@ 32 -b | \
    samtools sort -@ 32 -o "${header}.${ref}.sorted.bam"

samtools index \
    -@ 32 \
    "${header}.${ref}.sorted.bam"

conda deactivate