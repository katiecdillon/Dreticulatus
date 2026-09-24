#!/bin/bash
#SBATCH --job-name=gffread
#SBATCH --output=gffread_%j.out
#SBATCH --error=gffread_%j.err
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cag74075@uga.edu

module load gffread/0.12.7-GCCcore-12.3.0

# Change these paths to your file paths~

DEVON_GFF="/lustre2/scratch/cag74075/Tick_Compare_2026/Devon_Guided.gff3"

DEVON_FASTA="/lustre2/scratch/cag74075/Tick_Compare_2026/Dretic_Racon_Polish_Guided_Clean.renamed.fasta"


LOUISE_GFF="/lustre2/scratch/cag74075/Tick_Compare_2026/Louise.augustus.hints.gff3"

LOUISE_FASTA="/lustre2/scratch/cag74075/Tick_Compare_2026/Louise.fna"


# Generate your Output folder and port your outputs to that Output folder~

OUT="/lustre2/scratch/cag74075/Tick_Compare_2026/Proteomes"

mkdir -p "$OUT"



echo "Extracting Devon proteins..."

gffread \
"$DEVON_GFF" \
-g "$DEVON_FASTA" \
-y "$OUT/Devon_Dretic_Racon.faa" \
-x "$OUT/Devon_Dretic_Racon.cds.fna"



echo "Extracting Louise proteins..."

gffread \
"$LOUISE_GFF" \
-g "$LOUISE_FASTA" \
-y "$OUT/Louise_Dretic_Elska.faa" \
-x "$OUT/Louise_Dretic_Elska.cds.fna"



echo "Protein counts"

echo "Devon:"
grep -c "^>" "$OUT/Devon_Dretic_Racon.faa"


echo "Louise:"
grep -c "^>" "$OUT/Louise_Dretic_Elska.faa"



echo "Finished"