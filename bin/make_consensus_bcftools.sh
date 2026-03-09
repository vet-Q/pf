#!/usr/bin/env bash
set -euo pipefail

VCF="$1"      # vcf.gz
REF="$2"      # genome fasta
BED="$3"      # region.bed (chr start end)
OUTFA="$4"    # consensus fasta
PREFIX="${5:-tmp}"

region=$(awk 'NR==1{print $1":"($2+1)"-"$3; exit}' "${BED}")

REFSLICE="${PREFIX}.ref.fasta"
samtools faidx "${REF}" "${region}" > "${REFSLICE}"

TRIMVCF="${PREFIX}.trim.vcf.gz"
bcftools view -R "${BED}" -Oz -o "${TRIMVCF}" "${VCF}"

tabix -f -p vcf "${TRIMVCF}"

bcftools consensus -f "${REFSLICE}" "${TRIMVCF}" > "${OUTFA}"

