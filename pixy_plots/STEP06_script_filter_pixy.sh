#!/bin/bash
#SBATCH --job-name=pixyFilter
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200gb
#SBATCH --export=NONE
#SBATCH --time=3-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=kcd88651@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS


source ~/.bashrc
conda activate pixy_env
ml BCFtools/1.23.1-GCC-13.3.0

Louise_wd='/scratch/kcd88651/ticks/Dermacentor_Reticulatus/pixy/Louise_merge'
DretUK_wd='/scratch/kcd88651/ticks/Dermacentor_Reticulatus/pixy/DretUK_merge'

for k in Louise DretUK
do
    if [[ "$k" == "Louise" ]]; then
    wd="$Louise_wd"
    pop="$Louise_wd/pixy_input/Louise_pop.txt"
    else
    wd="$DretUK_wd"
    pop="$DretUK_wd/pixy_input/DretUK_pop.txt"
    fi
for i in {1..11}
do
    cd "$wd/chr$i" || exit 1

    bcftools view \
      -e 'TYPE="indel" || F_MISSING > 0 || INFO/DP < 20 || INFO/DP > 300' \
      "chr$i.${k}merge.vcf.gz" \
      -Oz \
      -o "chr${i}.${k}merge_filtered.vcf.gz"

    bcftools index --csi -f "chr${i}.${k}merge_filtered.vcf.gz"

    cd "$wd/pixy_input/" || exit 1
    pixy \
    --stats pi watterson_theta tajima_d dxy fst \
    --vcf "$wd/chr$i/chr${i}.${k}merge_filtered.vcf.gz" \
    --populations "$pop" \
    --window_size 50000 \
    --n_cores 32 \
    --fst_type hudson \
    --fst_components \
    --tajima_components \
    --output_prefix "Chr${i}${k}" \
    --output_folder "$wd/pixy_input"
done
done

conda deactivate