#!/bin/bash

module load gffread/0.12.7-GCCcore-12.3.0
module load SeqKit/2.8.2

srcdir="/scratch/skc49482/Tick"
gtfdir="/scratch/skc49482/Katie_genespace/Tick_isoform"
outdir="/scratch/skc49482/Katie_genespace/Tick_isoform/genespace/peptide"
mkdir -p $outdir

# NCBI species
for species in D_albipictus D_andersoni D_silvarum D_variabilis; do
  echo "Extracting peptides for ${species}..."
  gffread ${gtfdir}/${species}_final.gtf \
    -g ${srcdir}/${species}.fasta \
    -S -y ${outdir}/${species}.fa
done

# Dret_EU
echo "Extracting peptides for Dret_EU..."
gffread ${gtfdir}/Dret_EU_isoform.gtf \
  -g ${srcdir}/Elska.ragtag.scaffold.fasta \
  -S -y ${outdir}/Dret_EU.fa

# Dret_UK
echo "Extracting peptides for Dret_UK..."
gffread ${gtfdir}/Dret_UK_isoform.gtf \
  -g ${srcdir}/D_retic_UK.fasta \
  -S -y ${outdir}/Dret_UK.fa

# Stats
echo "Peptide sequence stats:"
seqkit stats ${outdir}/*.fa

echo "Done"
