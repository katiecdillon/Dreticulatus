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

PROT_DIR="/lustre2/scratch/cag74075/Derm_compare/Proteomes"

#Protein Fastas
P1="$PROT_DIR/Dretic_Racon_Polish_Guided.faa"
P2="$PROT_DIR/Dermacentor_silvarum.faa"
P3="$PROT_DIR/Dermacentor_Reticulatus.faa"

# Id-tags
S1="DRETIC"
S2="DSILV"
S3="DRET"

# Output base attachment
OUTBASE="$SLURM_SUBMIT_DIR/OF_3way_${SLURM_JOB_ID}"
IN="$OUTBASE/input"
mkdir -p "$IN"

echo "=== Inputs ==="
ls -lh "$P1" "$P2" "$P3"

# check for file inputs
for f in "$P1" "$P2" "$P3"; do
  [[ -s "$f" ]] || { echo "[ERROR] Missing/empty: $f" >&2; exit 1; }
done

command -v orthofinder >/dev/null 2>&1 || { echo "[ERROR] orthofinder not in PATH"; exit 1; }

# -------------------------
# Clean FASTA:
# - rewrite headers to TAG|0000001 
# - write header_map.tsv for ease of use
# - remove '*' and CRLF and any unkowns.
# -------------------------
clean_fa () {
  local infa="$1"
  local outfa="$2"
  local map="$3"
  local prefix="$4"

  rm -f "$map"

  awk -v pref="$prefix" -v MAP="$map" '
    BEGIN{n=0}
    /^>/{
      n++
      old=substr($0,2)
      gsub(/\r/,"",old)
      new=sprintf("%s|%07d", pref, n)
      print new "\t" old >> MAP
      print ">" new
      next
    }
    {
      gsub(/\r/,"")
      gsub(/\*/,"")
      if(length($0)>0) print
    }
  ' "$infa" > "$outfa"
}

echo "[1/3] Cleaning proteomes into OrthoFinder input folder..."
clean_fa "$P1" "$IN/${S1}.faa" "$OUTBASE/${S1}.header_map.tsv" "$S1"
clean_fa "$P2" "$IN/${S2}.faa" "$OUTBASE/${S2}.header_map.tsv" "$S2"
clean_fa "$P3" "$IN/${S3}.faa" "$OUTBASE/${S3}.header_map.tsv" "$S3"

echo "Sequence counts after cleaning:"
for f in "$IN/"*.faa; do
  echo -n "$(basename "$f") : "
  grep -c '^>' "$f"
done

# -------------------------
# Run OrthoFinder
# -------------------------
echo "[2/3] Running OrthoFinder..."
THREADS="${SLURM_CPUS_PER_TASK:-16}"
orthofinder -f "$IN" -t "$THREADS" -a "$THREADS"

echo "[3/3] Done."
echo "OUTBASE:      $OUTBASE"
echo "INPUT:        $IN"
echo "HEADER MAPS:  $OUTBASE/*header_map.tsv"
echo "RESULTS:      $IN/OrthoFinder/Results_*"
echo "Key file for Colab (shared/unique orthogroups):"
echo "  $IN/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv"