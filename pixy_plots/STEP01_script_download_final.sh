#!/bin/bash
#SBATCH --job-name=SRA_access
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500gb
#SBATCH --export=NONE
#SBATCH --time=1-
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=kcd88651@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ---------------------------------------------------------------------------- #
# CREATE CONDA ENVIRONMENTS
# ---------------------------------------------------------------------------- #
source ~/.bashrc
conda create -n ncbi_tools_env -c bioconda entrez-direct=25.3 sra-tools=3.4.1 -y
conda create -n nanoplot_env -c bioconda nanoplot=1.46.2 -y

# ---------------------------------------------------------------------------- #
# DIRECTORY PATHS
# ---------------------------------------------------------------------------- #
Elska_dir='/scratch/kcd88651/ticks/D_reticulatus/Elska/raw_reads'
Louise_dir='/scratch/kcd88651/ticks/D_reticulatus/Louise/raw_reads'
Penny_dir='/scratch/kcd88651/ticks/D_reticulatus/Penny/raw_reads'
DretUK_dir='/scratch/kcd88651/ticks/D_reticulatus/DretUK/raw_reads'

THREADS=16

# ---------------------------------------------------------------------------- #
# DOWNLOAD, CONVERT, CONCATENATE, AND COMPRESS FOR ELSKA, LOUISE, AND PENNY
# ---------------------------------------------------------------------------- #
conda activate ncbi_tools_env

for i in Elska Louise Penny
do
    dir_var="${i}_dir"
    sample_dir="${!dir_var}"

    cd "$sample_dir"

    # Get accessions
    esearch -db sra -query "PRJNA1130541[BioProject]" | efetch -format runinfo > "${i}_runinfo.csv"
    cut -d',' -f1 "${i}_runinfo.csv" | tail -n +2 > "${i}_sra_accessions.txt"
    echo "Accessions for ${i}:"
    cat "${i}_sra_accessions.txt"

    # Prefetch
    prefetch --option-file "${i}_sra_accessions.txt" --output-directory "$sample_dir"

    # fasterq-dump
    echo "Converting SRA to FASTQ for ${i}..."
    while read acc; do
        fasterq-dump "$sample_dir/$acc/$acc.sra" \
            --outdir "$sample_dir" \
            --threads $THREADS
    done < "${i}_sra_accessions.txt"

    # Concatenate
    echo "Concatenating FASTQ files for ${i}..."
    cat "$sample_dir"/*.fastq > "$sample_dir/${i}_combined.fastq"

    # Compress
    echo "Compressing FASTQ files for ${i}..."
    gzip "$sample_dir/${i}_combined.fastq"

    # Cleanup
    echo "Cleaning up intermediate files for ${i}..."
    find "$sample_dir" -name "*.fastq" -delete

    echo "Done with ${i}! Final files:"
    ls -lh "$sample_dir/${i}_combined.fastq.gz"
done

# ---------------------------------------------------------------------------- #
# DOWNLOAD, CONVERT, CONCATENATE, AND COMPRESS FOR DRETUK
# ---------------------------------------------------------------------------- #
cd "$DretUK_dir"

# Get accessions
esearch -db sra -query "PRJNA1418530[BioProject]" | efetch -format runinfo > DretUK_runinfo.csv
cut -d',' -f1 DretUK_runinfo.csv | tail -n +2 > DretUK_sra_accessions.txt
echo "Accessions for DretUK:"
cat DretUK_sra_accessions.txt

# Prefetch
prefetch --option-file DretUK_sra_accessions.txt --output-directory "$DretUK_dir"

# fasterq-dump
echo "Converting SRA to FASTQ for DretUK..."
while read acc; do
    fasterq-dump "$DretUK_dir/$acc/$acc.sra" \
        --outdir "$DretUK_dir" \
        --threads $THREADS
done < DretUK_sra_accessions.txt

# Concatenate
echo "Concatenating FASTQ files for DretUK..."
cat "$DretUK_dir"/*.fastq > "$DretUK_dir/DretUK_combined.fastq"

# Compress
echo "Compressing FASTQ files for DretUK..."
gzip "$DretUK_dir/DretUK_combined.fastq"

# Cleanup
echo "Cleaning up intermediate files for DretUK..."
find "$DretUK_dir" -name "*.fastq" -delete

echo "Done with DretUK! Final files:"
ls -lh "$DretUK_dir/DretUK_combined.fastq.gz"

conda deactivate

# ---------------------------------------------------------------------------- #
# NANOPLOT QC
# ---------------------------------------------------------------------------- #
conda activate nanoplot_env

for i in Elska Louise Penny DretUK
do
    dir_var="${i}_dir"
    sample_dir="${!dir_var}"

    cd "$sample_dir"

    NanoPlot --fastq "${i}_combined.fastq.gz" --info_in_report --N50 --outdir "nanoplot_raw_${i}"
done

conda deactivate
