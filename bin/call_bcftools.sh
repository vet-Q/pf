#!/usr/bin/env bash
set -euo pipefail

BAM="$1"
REF="$2"
BED="$3"
OUT="$4"            # *.vcf.gz
CALLMODE="${5:-mv}" # mv|cv
SNP_ONLY="${6:-true}"
MAXD="${7:-10000}"

if [[ "${CALLMODE}" == "cv" ]]; then
  CALL_ARGS="-cv"
else
  CALL_ARGS="-mv"
fi

# NOTE:
# - region 제한은 mpileup의 -R 로 이미 적용됨
# - 파이프 중간에 view -R 를 또 넣으면 stdin 처리 꼬이는 경우가 있어 제거
# - 최종 압축은 bcftools + -Oz 로 바로

if [[ "${SNP_ONLY}" == "true" ]]; then
  bcftools mpileup -Ou -f "${REF}" -R "${BED}" \
    --annotate FORMAT/DP,FORMAT/AD --max-depth "${MAXD}" "${BAM}" \
  | bcftools call ${CALL_ARGS} -Ou \
  | bcftools view -v snps -Oz -o "${OUT}"
else
  bcftools mpileup -Ou -f "${REF}" -R "${BED}" \
    --annotate FORMAT/DP,FORMAT/AD --max-depth "${MAXD}" "${BAM}" \
  | bcftools call ${CALL_ARGS} -Oz -o "${OUT}"
fi

tabix -f -p vcf "${OUT}"
