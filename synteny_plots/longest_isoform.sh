#!/bin/bash
# prepare 6 genomes for genespace run (LONGEST-ISOFORM)

input="/scratch/skc49482/Katie_genespace/genespace_aug/prep_input"
mkdir -p genespace_isoform/bed genespace_isoform/peptide

module load AGAT/1.4.2
module load gffread/0.12.7-GCCcore-12.3.0
module load SeqKit/2.8.2

## D_albipictus 
sp=D_albipictus
awk -F'\t' '$3=="gene" && /gene_biotype "protein_coding"/' $input/$sp.gtf | grep -oP 'gene_id "[^"]+"' | sort -u > ${sp}_pc.txt
grep -Ff ${sp}_pc.txt $input/$sp.gtf | grep -v '^#' > ${sp}_pc.gtf
agat_sp_keep_longest_isoform.pl -gff ${sp}_pc.gtf -o ${sp}_longest.gtf
awk -F'\t' '$3=="transcript" && !/transcript_biotype=mRNA/ {match($0,/ID=[^;]+/);print substr($0,RSTART+3,RLENGTH-3)}' ${sp}_longest.gtf > ${sp}_rm.txt
if [ -s ${sp}_rm.txt ]; then grep -vFf ${sp}_rm.txt ${sp}_longest.gtf > ${sp}_final.gtf; else cp ${sp}_longest.gtf ${sp}_final.gtf; fi
awk -F'\t' '$3=="transcript"{match($0,/ID=[^;]+/);print $1"\t"($4-1)"\t"$5"\t"substr($0,RSTART+3,RLENGTH-3)}' ${sp}_final.gtf > genespace_isoform/bed/$sp.bed
sed -i 's/^NC_091821.1\t/1\t/;s/^NC_091822.1\t/2\t/;s/^NC_091823.1\t/3\t/;s/^NC_091824.1\t/4\t/;s/^NC_091825.1\t/5\t/;s/^NC_091826.1\t/6\t/;s/^NC_091827.1\t/7\t/;s/^NC_091828.1\t/8\t/;s/^NC_091829.1\t/9\t/;s/^NC_091830.1\t/10\t/' genespace_isoform/bed/$sp.bed
awk -F'\t' '$1 ~ /^[0-9]+$/' genespace_isoform/bed/$sp.bed > t && mv t genespace_isoform/bed/$sp.bed
awk -F'\t' '$3=="transcript"||$3=="exon"||$3=="CDS"' ${sp}_final.gtf > ${sp}_gff.gtf
gffread ${sp}_gff.gtf -g $input/$sp.fasta -S -y genespace_isoform/peptide/$sp.fa
cut -f4 genespace_isoform/bed/$sp.bed | sort > ${sp}.ids; seqkit grep -f ${sp}.ids genespace_isoform/peptide/$sp.fa > t && mv t genespace_isoform/peptide/$sp.fa
echo "$sp: bed=$(wc -l <genespace_isoform/bed/$sp.bed) pep=$(grep -c '>' genespace_isoform/peptide/$sp.fa)"

## D_andersoni 
sp=D_andersoni
awk -F'\t' '$3=="gene" && /gene_biotype "protein_coding"/' $input/$sp.gtf | grep -oP 'gene_id "[^"]+"' | sort -u > ${sp}_pc.txt
grep -Ff ${sp}_pc.txt $input/$sp.gtf | grep -v '^#' > ${sp}_pc.gtf
agat_sp_keep_longest_isoform.pl -gff ${sp}_pc.gtf -o ${sp}_longest.gtf
awk -F'\t' '$3=="transcript" && !/transcript_biotype=mRNA/ {match($0,/ID=[^;]+/);print substr($0,RSTART+3,RLENGTH-3)}' ${sp}_longest.gtf > ${sp}_rm.txt
if [ -s ${sp}_rm.txt ]; then grep -vFf ${sp}_rm.txt ${sp}_longest.gtf > ${sp}_final.gtf; else cp ${sp}_longest.gtf ${sp}_final.gtf; fi
awk -F'\t' '$3=="transcript"{match($0,/ID=[^;]+/);print $1"\t"($4-1)"\t"$5"\t"substr($0,RSTART+3,RLENGTH-3)}' ${sp}_final.gtf > genespace_isoform/bed/$sp.bed
sed -i 's/^NC_092814.1\t/1\t/;s/^NC_092815.1\t/2\t/;s/^NC_092816.1\t/3\t/;s/^NC_092817.1\t/4\t/;s/^NC_092818.1\t/5\t/;s/^NC_092819.1\t/6\t/;s/^NC_092820.1\t/7\t/;s/^NC_092821.1\t/8\t/;s/^NC_092822.1\t/9\t/;s/^NC_092823.1\t/10\t/;s/^NC_092824.1\t/11\t/' genespace_isoform/bed/$sp.bed
awk -F'\t' '$1 ~ /^[0-9]+$/' genespace_isoform/bed/$sp.bed > t && mv t genespace_isoform/bed/$sp.bed
awk -F'\t' '$3=="transcript"||$3=="exon"||$3=="CDS"' ${sp}_final.gtf > ${sp}_gff.gtf
gffread ${sp}_gff.gtf -g $input/$sp.fasta -S -y genespace_isoform/peptide/$sp.fa
cut -f4 genespace_isoform/bed/$sp.bed | sort > ${sp}.ids; seqkit grep -f ${sp}.ids genespace_isoform/peptide/$sp.fa > t && mv t genespace_isoform/peptide/$sp.fa
echo "$sp: bed=$(wc -l <genespace_isoform/bed/$sp.bed) pep=$(grep -c '>' genespace_isoform/peptide/$sp.fa)"

## D_silvarum
sp=D_silvarum
awk -F'\t' '$3=="gene" && /gene_biotype "protein_coding"/' $input/$sp.gtf | grep -oP 'gene_id "[^"]+"' | sort -u > ${sp}_pc.txt
grep -Ff ${sp}_pc.txt $input/$sp.gtf | grep -v '^#' > ${sp}_pc.gtf
agat_sp_keep_longest_isoform.pl -gff ${sp}_pc.gtf -o ${sp}_longest.gtf
awk -F'\t' '$3=="transcript" && !/transcript_biotype=mRNA/ {match($0,/ID=[^;]+/);print substr($0,RSTART+3,RLENGTH-3)}' ${sp}_longest.gtf > ${sp}_rm.txt
if [ -s ${sp}_rm.txt ]; then grep -vFf ${sp}_rm.txt ${sp}_longest.gtf > ${sp}_final.gtf; else cp ${sp}_longest.gtf ${sp}_final.gtf; fi
awk -F'\t' '$3=="transcript"{match($0,/ID=[^;]+/);print $1"\t"($4-1)"\t"$5"\t"substr($0,RSTART+3,RLENGTH-3)}' ${sp}_final.gtf > genespace_isoform/bed/$sp.bed
sed -i 's/^NC_051154.1\t/1\t/;s/^NC_051155.1\t/2\t/;s/^NC_051156.1\t/3\t/;s/^NC_051157.2\t/4\t/;s/^NC_051158.1\t/5\t/;s/^NC_051159.1\t/6\t/;s/^NC_051160.1\t/7\t/;s/^NC_051161.1\t/8\t/;s/^NC_051162.1\t/9\t/;s/^NC_051163.1\t/10\t/;s/^NC_051164.1\t/11\t/' genespace_isoform/bed/$sp.bed
awk -F'\t' '$1 ~ /^[0-9]+$/' genespace_isoform/bed/$sp.bed > t && mv t genespace_isoform/bed/$sp.bed
awk -F'\t' '$3=="transcript"||$3=="exon"||$3=="CDS"' ${sp}_final.gtf > ${sp}_gff.gtf
gffread ${sp}_gff.gtf -g $input/$sp.fasta -S -y genespace_isoform/peptide/$sp.fa
cut -f4 genespace_isoform/bed/$sp.bed | sort > ${sp}.ids; seqkit grep -f ${sp}.ids genespace_isoform/peptide/$sp.fa > t && mv t genespace_isoform/peptide/$sp.fa
echo "$sp: bed=$(wc -l <genespace_isoform/bed/$sp.bed) pep=$(grep -c '>' genespace_isoform/peptide/$sp.fa)"

## D_variabilis 
sp=D_variabilis
awk -F'\t' '$3=="gene" && /gene_biotype "protein_coding"/' $input/$sp.gtf | grep -oP 'gene_id "[^"]+"' | sort -u > ${sp}_pc.txt
grep -Ff ${sp}_pc.txt $input/$sp.gtf | grep -v '^#' > ${sp}_pc.gtf
agat_sp_keep_longest_isoform.pl -gff ${sp}_pc.gtf -o ${sp}_longest.gtf
awk -F'\t' '$3=="transcript" && !/transcript_biotype=mRNA/ {match($0,/ID=[^;]+/);print substr($0,RSTART+3,RLENGTH-3)}' ${sp}_longest.gtf > ${sp}_rm.txt
if [ -s ${sp}_rm.txt ]; then grep -vFf ${sp}_rm.txt ${sp}_longest.gtf > ${sp}_final.gtf; else cp ${sp}_longest.gtf ${sp}_final.gtf; fi
awk -F'\t' '$3=="transcript"{match($0,/ID=[^;]+/);print $1"\t"($4-1)"\t"$5"\t"substr($0,RSTART+3,RLENGTH-3)}' ${sp}_final.gtf > genespace_isoform/bed/$sp.bed
sed -i 's/^NC_134568.1\t/1\t/;s/^NC_134569.1\t/2\t/;s/^NC_134570.1\t/3\t/;s/^NC_134571.1\t/4\t/;s/^NC_134572.1\t/5\t/;s/^NC_134573.1\t/6\t/;s/^NC_134574.1\t/7\t/;s/^NC_134575.1\t/8\t/;s/^NC_134576.1\t/9\t/;s/^NC_134577.1\t/10\t/;s/^NC_134578.1\t/11\t/' genespace_isoform/bed/$sp.bed
awk -F'\t' '$1 ~ /^[0-9]+$/' genespace_isoform/bed/$sp.bed > t && mv t genespace_isoform/bed/$sp.bed
awk -F'\t' '$3=="transcript"||$3=="exon"||$3=="CDS"' ${sp}_final.gtf > ${sp}_gff.gtf
gffread ${sp}_gff.gtf -g $input/$sp.fasta -S -y genespace_isoform/peptide/$sp.fa
cut -f4 genespace_isoform/bed/$sp.bed | sort > ${sp}.ids; seqkit grep -f ${sp}.ids genespace_isoform/peptide/$sp.fa > t && mv t genespace_isoform/peptide/$sp.fa
echo "$sp: bed=$(wc -l <genespace_isoform/bed/$sp.bed) pep=$(grep -c '>' genespace_isoform/peptide/$sp.fa)"

## DretUK 
sp=DretUK
agat_sp_keep_longest_isoform.pl -gff $input/D_retic_UK.gtf -o ${sp}_longest.gtf
awk -F'\t' '$3=="transcript"{match($0,/ID=[^;]+/);print $1"\t"($4-1)"\t"$5"\t"substr($0,RSTART+3,RLENGTH-3)}' ${sp}_longest.gtf > genespace_isoform/bed/$sp.bed
sed -i 's/_RagTag//; s/^Dretic_Chr_//' genespace_isoform/bed/$sp.bed
awk -F'\t' '$1 ~ /^[0-9]+$/' genespace_isoform/bed/$sp.bed > t && mv t genespace_isoform/bed/$sp.bed
awk -F'\t' '$3=="transcript"||$3=="exon"||$3=="CDS"' ${sp}_longest.gtf > ${sp}_gff.gtf
gffread ${sp}_gff.gtf -g $input/D_retic_UK.fasta -S -y genespace_isoform/peptide/$sp.fa
cut -f4 genespace_isoform/bed/$sp.bed | sort > ${sp}.ids; seqkit grep -f ${sp}.ids genespace_isoform/peptide/$sp.fa > t && mv t genespace_isoform/peptide/$sp.fa
echo "$sp: bed=$(wc -l <genespace_isoform/bed/$sp.bed) pep=$(grep -c '>' genespace_isoform/peptide/$sp.fa)"

## DretLouise (AUGUSTUS, AGAT longest, use .aa file, CM accessions -> 1-11)
sp=DretLouise
agat_sp_keep_longest_isoform.pl -gff $input/Louise_ragtag_augustus.hints.gtf -o ${sp}_longest.gtf
awk -F'\t' '$3=="transcript"{match($0,/ID=[^;]+/);print $1"\t"($4-1)"\t"$5"\t"substr($0,RSTART+3,RLENGTH-3)}' ${sp}_longest.gtf > genespace_isoform/bed/$sp.bed
sed -i 's/_RagTag//' genespace_isoform/bed/$sp.bed
sed -i 's/^CM170536.1\t/1\t/;s/^CM170537.1\t/2\t/;s/^CM170538.1\t/3\t/;s/^CM170539.1\t/4\t/;s/^CM170540.1\t/5\t/;s/^CM170541.1\t/6\t/;s/^CM170542.1\t/7\t/;s/^CM170543.1\t/8\t/;s/^CM170544.1\t/9\t/;s/^CM170545.1\t/10\t/;s/^CM170546.1\t/11\t/' genespace_isoform/bed/$sp.bed
awk -F'\t' '$1 ~ /^[0-9]+$/' genespace_isoform/bed/$sp.bed > t && mv t genespace_isoform/bed/$sp.bed
# peptide: filter the AUGUSTUS .aa to the longest-isoform IDs
cp $input/Louise_ragtag_augustus.hints.aa genespace_isoform/peptide/$sp.fa
cut -f4 genespace_isoform/bed/$sp.bed | sort > ${sp}.ids; seqkit grep -f ${sp}.ids genespace_isoform/peptide/$sp.fa > t && mv t genespace_isoform/peptide/$sp.fa
echo "$sp: bed=$(wc -l <genespace_isoform/bed/$sp.bed) pep=$(grep -c '>' genespace_isoform/peptide/$sp.fa)"

echo "=== DONE ==="
