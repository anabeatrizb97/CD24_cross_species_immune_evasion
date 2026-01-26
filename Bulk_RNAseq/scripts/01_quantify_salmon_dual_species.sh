#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
#  Dual-species transcript quantification with Salmon
# - Human tumor (GRCh38) and zebrafish host/TME (GRCz11) quantified separately.
#
# Inputs:
# - Raw FASTQ files (lane-split, single-end, 76bp) 
#    - All files in one folder, e.g.: B_SW480_L001.fastq.gz, B_SW480_L002.fastq.gz, ...
# - Reference transcriptomes FASTA files (Ensemble FTP; release 115)
#    - Human "Homo_sapiens.GRCh38.cdna.all.fa.gz"
#      (from https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/)
#    - Zebrafish "Danio_rerio.GRCz11.cdna.all.fa.gz"
#      (from https://ftp.ensembl.org/pub/release-115/fasta/danio_rerio/cdna/)
#
# Output:
# - Salmon quant directories per sample per species
#
# Path structure:
# <PROJECT_ROOT>/
#     data/fastq/
#     refs/
#       Homo_sapiens.GRCh38.cdna.all.fa.gz
#       Danio_rerio.GRCz11.cdna.all.fa.gz
#       salmon_index/
#         human_k31/
#         zebrafish_k31/
#     results/salmon/
#       human/<SAMPLE>/
#       zebrafish/<SAMPLE>/
#
# Software:
# - salmon v1.10.3
# ==============================================================================

# ------------------------------
# User-configurable parameters
# ------------------------------

# Project root (default - directory containing this script)
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input FASTQ folder (raw FASTQ files, lane-split single-end .fastq.gz)
FASTQ_DIR="${PROJECT_ROOT}/data/fastq"

# Reference folder (FASTA files for each species reference transcriptome)
REF_DIR="${PROJECT_ROOT}/refs"

# Output folder
OUT_DIR="${PROJECT_ROOT}/results/salmon"

# Reference FASTA files
HUMAN_FA="${REF_DIR}/Homo_sapiens.GRCh38.cdna.all.fa.gz"
ZEBRA_FA="${REF_DIR}/Danio_rerio.GRCz11.cdna.all.fa.gz"

# Salmon index folders (created if missing)
INDEX_DIR="${REF_DIR}/salmon_index"
KMER=31
HUMAN_INDEX="${INDEX_DIR}/human_k${KMER}"
ZEBRA_INDEX="${INDEX_DIR}/zebrafish_k${KMER}"

THREADS=8

# Samples (must match FASTQ name prefix)
SAMPLES=("B_SW480" "C_SW620" "I_SW620" "L_SW620" "N_SW480" "O_SW620" "P_SW480")

# ------------------------------
# Sanity checks
# ------------------------------
command -v salmon >/dev/null 2>&1 || { echo "ERROR: salmon not found in PATH"; exit 1; }

[[ -d "${FASTQ_DIR}" ]] || { echo "ERROR: FASTQ_DIR not found: ${FASTQ_DIR}"; exit 1; }
[[ -f "${HUMAN_FA}" ]]  || { echo "ERROR: Missing human FASTA: ${HUMAN_FA}"; exit 1; }
[[ -f "${ZEBRA_FA}" ]]  || { echo "ERROR: Missing zebrafish FASTA: ${ZEBRA_FA}"; exit 1; }

mkdir -p "${OUT_DIR}/human" "${OUT_DIR}/zebrafish" "${INDEX_DIR}"

# Summary log
echo "Salmon version:"
salmon --version || true
echo "K-mer: ${KMER} | Threads: ${THREADS}"
echo "Human FASTA: ${HUMAN_FA}"
echo "Zebrafish FASTA: ${ZEBRA_FA}"
echo "FASTQ_DIR: ${FASTQ_DIR}"
echo "OUT_DIR: ${OUT_DIR}"
echo "INDEX_DIR: ${INDEX_DIR}"

# ------------------------------
# Build indices (if missing)
# ------------------------------
if [[ ! -d "${HUMAN_INDEX}" ]]; then
  echo "Building human Salmon index..."
  salmon index -t "${HUMAN_FA}" -i "${HUMAN_INDEX}" -k "${KMER}"
else
  echo "Human index exists: ${HUMAN_INDEX} (skipping)"
fi

if [[ ! -d "${ZEBRA_INDEX}" ]]; then
  echo "Building zebrafish Salmon index..."
  salmon index -t "${ZEBRA_FA}" -i "${ZEBRA_INDEX}" -k "${KMER}"
else
  echo "Zebrafish index exists: ${ZEBRA_INDEX} (skipping)"
fi

# ------------------------------
# Quantification
# ------------------------------
for SAMPLE in "${SAMPLES[@]}"; do
  echo "=============================================================================="
  echo "Sample: ${SAMPLE}"

  # Find all lane FASTQs for this sample (common Illumina naming)
  mapfile -t LANES < <(
    find "${FASTQ_DIR}" -maxdepth 1 -type f \
      \( -name "${SAMPLE}*_L00*.fastq.gz" -o -name "${SAMPLE}*_L00*.fq.gz" \) \
      | sort
  )

  if [[ "${#LANES[@]}" -eq 0 ]]; then
    echo "ERROR: No lane FASTQs found for ${SAMPLE} in ${FASTQ_DIR}"
    echo "Tip: check naming with: ls \"${FASTQ_DIR}/${SAMPLE}*\""
    exit 1
  fi

  echo "  Found lanes:"
  printf "    %s\n" "${LANES[@]}"

  OUT_H="${OUT_DIR}/human/${SAMPLE}"
  OUT_Z="${OUT_DIR}/zebrafish/${SAMPLE}"
  mkdir -p "${OUT_H}" "${OUT_Z}"

  echo "  Quantifying against human..."
  salmon quant \
    -i "${HUMAN_INDEX}" \
    -l A \
    -r "${LANES[@]}" \
    -p "${THREADS}" \
    --validateMappings \
    -o "${OUT_H}"

  echo "  Quantifying against zebrafish..."
  salmon quant \
    -i "${ZEBRA_INDEX}" \
    -l A \
    -r "${LANES[@]}" \
    -p "${THREADS}" \
    --validateMappings \
    -o "${OUT_Z}"

  echo "  Done: ${SAMPLE}"
done

echo "=============================================================================="
echo "All samples processed."
echo "Outputs:"
echo "  Human:     ${OUT_DIR}/human/<SAMPLE>/"
echo "  Zebrafish: ${OUT_DIR}/zebrafish/<SAMPLE>/"