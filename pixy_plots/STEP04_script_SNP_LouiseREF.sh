#!/bin/bash
#SBATCH --job-name=SNP
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200gb
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
Louise_asmbly='scratch/kcd88651/ticks/D_reticulatus/Louise/Louise_assembly.fna'
wd='/scratch/kcd88651/ticks/D_reticulatus/pixy/LouiseREF'
LvD_bam='/scratch/kcd88651/ticks/D_reticulatus/pixy/DretUK_ref/LouiseFNA.DretUKFNA.sorted.bam'

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
if [[ -f LouiseREF.bamlist.txt ]]; then rm LouiseREF.bamlist.txt; fi ## avoid adding to list twice
for i in Elska Louise Penny DretUK
do
    echo "$i.Louise.sorted.bam" >> LouiseREF.bamlist.txt
done
bamlist="LouiseREF.bamlist.txt"

# regions file contains the Louise assembly contigs that aligned...
# to the 11 chromosome-scale scaffolds of the DretUK assembly
regions="Louise.DretUK_regions_chr${numb}.txt"
samtools view -q 20 -F 0x900 "$LvD_bam" "$chr" | cut -f1 | sort -u > "$regions"

n_contigs=$(wc -l < "$regions")
echo "chr $numb ($chr): $n_contigs Louise contigs"
if [[ "$n_contigs" -eq 0 ]]; then
    echo "WARNING: no Louise contigs mapped to $chr; skipping VCF for chr $numb." >&2
    exit 0
fi

# generate SNPs and indels
bcftools mpileup --threads 32 -f "$Louise_asmbly" -b "$bamlist" -R "$regions" \
    | bcftools call --threads 32 -m -O z -a GQ -o "chr${numb}.Louise.vcf.gz"

echo "Done: chr${numb}.Louise.vcf.gz"

conda deactivate
