#!/bin/bash
#SBATCH --job-name=orthofinder_3way
#SBATCH --output=orthofinder_3way_%j.out
#SBATCH --error=orthofinder_3way_%j.err
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cag74075@uga.edu


module load OrthoFinder/3.1.0-foss-2023a

cd /lustre2/scratch/cag74075/Tick_Compare_2026

orthofinder \
-f OrthoFinder_input \
-t 16 \
-a 16
-o Tick_Compare_2026/ Orthofinder_Threeway