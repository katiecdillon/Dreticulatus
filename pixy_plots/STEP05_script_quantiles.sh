#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200gb
#SBATCH --export=NONE
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=kcd88651@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

ml BCFtools/1.23.1-GCC-13.3.0 R/4.5.1-gfbf-2025a

Louise_wd='/scratch/kcd88651/ticks/D_reticulatus/pixy/LouiseREF'
DretUK_wd='/scratch/kcd88651/ticks/D_reticulatus/pixy/DretUKref'

for k in Louise DretUK
do
    # choose working directory
    if [[ "$k" == "Louise" ]]; then
        wd="$Louise_wd"
        sample="Louise"
    else
        wd="$DretUK_wd"
        sample="DretUK"
    fi

    cd "$wd"

for i in {1..11}
do
    cd chr$i

    bcftools query -f '%INFO/DP\n' "${wd}/chr${i}/chr${i}.${sample}.vcf.gz" | \
    Rscript -e "
    depths <- scan(file='stdin')

    sink('${wd}/chr${i}_depth_quantiles.txt')

    cat('Upper quantiles\n')
    print(quantile(depths, probs=c(0.95,0.99,0.995), na.rm=TRUE))

    cat('\nLower quantiles\n')
    print(quantile(depths, probs=c(0.05,0.01,0.005), na.rm=TRUE))

    sink()
    "

    cp chr${i}_depth_quantiles.txt $wd
    cd ..
done

(
echo -e "chromosome\tquantile\tvalue"

for f in chr*_depth_quantiles.txt
do
    awk -v file="$f" '
    BEGIN {
        sub(/_depth_quantiles\.txt$/, "", file)
    }

    /Upper quantiles|Lower quantiles/ {next}

    /%/ {
        p1=$1
        p2=$2
        p3=$3

        getline

        print file "\t" p1 "\t" $1
        print file "\t" p2 "\t" $2
        print file "\t" p3 "\t" $3
    }
    ' "$f"

done
done
) > "${k}_depth_quantiles.txt"
