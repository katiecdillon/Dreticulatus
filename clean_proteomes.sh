#!/bin/bash
#SBATCH --job-name=clean_proteomes
#SBATCH --output=clean_proteomes_%j.out
#SBATCH --error=clean_proteomes_%j.err
#SBATCH --ntasks=1
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

cd /lustre2/scratch/cag74075/Tick_Compare_2026

rm -rf OrthoFinder_input
mkdir OrthoFinder_input

for file in Proteomes/*.faa
do
    name=$(basename "$file")

    awk '
    /^>/ {
        sub(/\r$/, "")
        print
        next
    }
    {
        sub(/\r$/, "")
        gsub(/\*/, "")
        gsub(/\./, "X")
        gsub(/[[:space:]]/, "")
        if (length($0) > 0) print
    }
    ' "$file" > "OrthoFinder_input/$name"
done

echo "Protein counts:"
grep -c "^>" OrthoFinder_input/*.faa

echo "Checking for invalid characters:"
grep -v "^>" OrthoFinder_input/*.faa |
grep -n '[^ABCDEFGHIKLMNPQRSTVWXYZOUJ]' |
head || echo "No invalid amino-acid characters found."

echo "Finished."