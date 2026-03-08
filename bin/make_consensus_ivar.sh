#!/usr/bin/env bash
set -euo pipefail

BAM="$1"
REF="$2"
BED="$3"
OUT_PREFIX="$4"      # 예: barcode/GENE/barcode_GENE_ivar
T="${5:-0.6}"        # allele freq threshold for ambiguity
M="${6:-10}"         # min depth for calling
OUTFA="${7:-${OUT_PREFIX}.consensus.fasta}"

# ivar consensus 역시 -p prefix로 파일 생성(prefix.fa)
samtools mpileup -A -d 0 -Q 0 \
  -f "${REF}" \
  -l "${BED}" \
  "${BAM}" \
| ivar consensus \
    -p "${OUT_PREFIX}" \
    -t "${T}" \
    -m "${M}" \
> /dev/null

# 결과 파일명은 보통 ${OUT_PREFIX}.fa
if [[ -f "${OUT_PREFIX}.fa" ]]; then
  mv -f "${OUT_PREFIX}.fa" "${OUTFA}"
else
  echo "ERROR: iVar consensus output FASTA not found: ${OUT_PREFIX}.fa" >&2
  exit 1
fi

