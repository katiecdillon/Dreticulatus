#!/bin/bash
#SBATCH --job-name=SNP
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=900gb
#SBATCH --export=NONE
#SBATCH --time=10-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=kcd88651@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS
#SBATCH --array=1-11


# ---------------------------------------------------------------------------- #
# CREATE CONDA ENVIRONMENTS
# ---------------------------------------------------------------------------- #
source ~/.bashrc

conda create -n snp_env -c bioconda bcftools=1.23.1 samtools=1.21 -y

# ---------------------------------------------------------------------------- #
# FILE PATHS
# ---------------------------------------------------------------------------- #
DretUK_asmbly='scratch/kcd88651/ticks/D_reticulatus/DretUK/DretUK_assembly.fna'
wd='/scratch/kcd88651/ticks/D_reticulatus/pixy/DretUKref'

config=${wd}/STEP04_config_SNP_DretUKref.txt
numb=$SLURM_ARRAY_TASK_ID
chr=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $2}' "$config")

# ---------------------------------------------------------------------------- #
# VARIANT CALLING - DRETUK REFERENCE
# ---------------------------------------------------------------------------- #
conda activate snp_env
cd "$wd"

# make chromosome directories
for i in {1..11}
do
    mkdir -p chr${i}
done

# make a list of bam files
if [[ -f DretUKref.bamlist.txt ]]; then rm DretUKref.bamlist.txt; fi ## avoid adding to list twice
for i in Elska Louise Penny
do
    echo "$i.DretUK.sorted.bam" >> DretUKref.bamlist.txt
done
bamlist="DretUKref.bamlist.txt"

# generate SNPs and indels
cd "$wd/chr${numb}" || exit 1
bcftools mpileup --threads 32 -f "$DretUK_a" -b $bamlist -r "$chr" \
    | bcftools call --threads 32 -m -O z -a GQ -o chr${numb}.DretUKref.vcf.gz

echo "Done: chr${numb}.DretUKref.vcf.gz"

conda deactivate
