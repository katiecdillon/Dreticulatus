#!/bin/bash
#SBATCH --job-name=dorado
#SBATCH --partition=gpu_p
#SBATCH --gres=gpu:A100:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=900gb
#SBATCH --export=NONE
#SBATCH --time=10-
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=kcd88651@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-4


# ---------------------------------------------------------------------------- #
# DIRECTORY PATHS
# ---------------------------------------------------------------------------- #
wd='/scratch/kcd88651/ticks/D_reticulatus'

# ---------------------------------------------------------------------------- #
# FILE PATHS
# ---------------------------------------------------------------------------- #
config=$wd/config_dorado.txt

dir=$(awk -v Array_ID=$SLURM_ARRAY_TASK_ID '$1==Array_ID {print $2}' $config)
fastqgz=$(awk -v Array_ID=$SLURM_ARRAY_TASK_ID '$1==Array_ID {print $3}' $config)
fastq=$(awk -v Array_ID=$SLURM_ARRAY_TASK_ID '$1==Array_ID {print $4}' $config)
correct=$(awk -v Array_ID=$SLURM_ARRAY_TASK_ID '$1==Array_ID {print $5}' $config)

# ---------------------------------------------------------------------------- #
# DORADO
# ---------------------------------------------------------------------------- #
ml dorado/2.0.0-foss-2024a-CUDA-12.6.0

# Make sure output directory exists
mkdir -p "$dir"

# unzip the compressed FASTQ files
cd "$dir"
gunzip -c "$fastqgz" > "$fastq"

# Run 'dorado correct' for Elska, Louise, Penny, and DretUK
dorado correct --threads 64 "$fastq" > "$correct"
