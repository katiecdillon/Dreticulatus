#!/bin/bash

gtfdir="/scratch/skc49482/Katie_genespace/Tick_isoform/"
outdir="/scratch/skc49482/Katie_genespace/Tick_isoform/genespace/bed"
mkdir -p $outdir

# NCBI species
for species in D_albipictus D_andersoni D_silvarum D_variabilis; do
  echo "Making BED for ${species}..."
  awk -F'\t' '$3=="transcript" {
    match($0, /ID=[^;]+/); tid=substr($0, RSTART+3, RLENGTH-3)
    print $1"\t"$4"\t"$5"\t"tid
  }' ${gtfdir}/${species}_final.gtf > ${outdir}/${species}.bed
done

# Dret_EU
echo "Making BED for Dret_EU..."
awk -F'\t' '$3=="transcript" {
  match($0, /ID=[^;]+/); tid=substr($0, RSTART+3, RLENGTH-3)
  print $1"\t"$4"\t"$5"\t"tid
}' ${gtfdir}/Dret_EU_isoform.gtf > ${outdir}/Dret_EU.bed

# Dret_UK
echo "Making BED for Dret_UK..."
awk -F'\t' '$3=="transcript" {
  match($0, /ID=[^;]+/); tid=substr($0, RSTART+3, RLENGTH-3)
  print $1"\t"$4"\t"$5"\t"tid
}' ${gtfdir}/Dret_UK_isoform.gtf > ${outdir}/Dret_UK.bed

# Verify
echo "BED file stats:"
for f in ${outdir}/*.bed; do
  echo "$(basename $f): $(wc -l < $f) entries"
  head -2 $f
  echo ""
done

# Check BED-peptide ID concordance
pepdir="/scratch/skc49482/Katie_genespace/Tick_isoform/genespace/peptide"
echo "Checking BED-peptide ID matches:"
for f in ${outdir}/*.bed; do
  sp=$(basename $f .bed)
  mismatches=$(diff <(cut -f4 $f | sort) <(grep "^>" ${pepdir}/${sp}.fa | sed 's/>//' | sort) | wc -l)
  if [ "$mismatches" -eq 0 ]; then
    echo "  ${sp}: PASS"
  else
    echo "  ${sp}: MISMATCH ($mismatches differences)"
  fi
done

echo "Done"
