#!/bin/bash
#SBATCH --job-name=longest_isoform
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=02:00:00
#SBATCH --output=longest_isoform_%j.out
#SBATCH --error=longest_isoform_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=skc49482@uga.edu

start_time=`date +%s`

export TMPDIR=/scratch/skc49482/Rtmp
mkdir -p $TMPDIR

module load AGAT/1.4.2

srcdir="/scratch/skc49482/Tick"
cd /scratch/skc49482/Katie_genespace/Tick_isoform

for species in D_albipictus D_andersoni D_silvarum D_variabilis; do
  echo "Processing ${species}..."

  # Step 1: Extract protein-coding gene IDs from GTF
  awk -F'\t' '$3=="gene" && /gene_biotype "protein_coding"/' ${srcdir}/${species}.gtf | \
    grep -oP 'gene_id "[^"]+"' | sort -u > ${species}_pc_ids.txt

  # Step 2: Pull all lines for those genes
  grep -Ff ${species}_pc_ids.txt ${srcdir}/${species}.gtf | grep -v '^#' > ${species}_pc.gtf

  # Step 3: AGAT longest isoform
  agat_sp_keep_longest_isoform.pl -gff ${species}_pc.gtf -o ${species}_longest.gtf

  # Step 4: Remove non-mRNA transcripts and their children
  awk -F'\t' '$3=="transcript" && !/transcript_biotype=mRNA/' ${species}_longest.gtf | \
    grep -oP 'ID=[^;]+' > ${species}_remove.txt
  grep -vFf ${species}_remove.txt ${species}_longest.gtf > ${species}_final.gtf

  # Verify
  echo "Stats for ${species}_final.gtf"
  cut -f3 ${species}_final.gtf | sort | uniq -c
done

## D_reticulatus and D_retic_UK (renamed for GENESPACE compatibility)
D_reticulatus="/scratch/skc49482/Tick/ElskaRT.augustus.hints.gtf"
D_retic_UK="/scratch/skc49482/Tick/D_retic_UK.gtf"

agat_sp_keep_longest_isoform.pl -gff "$D_reticulatus" -o Dret_EU_isoform.gtf
agat_sp_keep_longest_isoform.pl -gff "$D_retic_UK" -o Dret_UK_isoform.gtf

echo "Stats for Dret_EU"
cut -f3 Dret_EU_isoform.gtf | sort | uniq -c
echo "Stats for Dret_UK"
cut -f3 Dret_UK_isoform.gtf | sort | uniq -c

echo "Done"

end_time=`date +%s`
runtime=$((end_time - start_time))
runtimeH=$((runtime / 3600))
runtimeM=$(((runtime % 3600) / 60))
runtimeS=$((runtime % 60))

echo "Duration: $runtime seconds"
echo "Duration: $runtimeH hours, $runtimeM minutes, $runtimeS seconds"
