#!/usr/bin/env bash
set -euo pipefail

BAM="$1"
REF="$2"
BED="$3"
OUT_PREFIX="$4"          # 예: .../barcode..._GENE_ivar  (여기서 .all 붙일지 말지 우리가 결정)
MIN_ALT_FRAC="${5:-0}"   # 0
MIN_DP="${6:-1}"         # 1
SNP_ONLY="${7:-true}"    # true|false

RAW_TSV="${OUT_PREFIX}.all.tsv"
SNP_TSV="${OUT_PREFIX}.snp.tsv"

# iVar variants: -p prefix 로 prefix.tsv 생성
# 여기서는 "원래 iVar가 만드는 파일이 RAW_TSV와 동일"하도록 prefix를 설계하는 게 제일 깔끔함.
# 즉, -p에 "${OUT_PREFIX}.all" 을 주면 결과가 "${OUT_PREFIX}.all.tsv"가 된다.
samtools mpileup -A -d 0 -Q 0 \
  -f "${REF}" \
  -l "${BED}" \
  "${BAM}" \
| ivar variants \
    -p "${OUT_PREFIX}.all" \
    -q 0 \
    -t "${MIN_ALT_FRAC}" \
    -m "${MIN_DP}" \
> /dev/null

# 여기서 iVar가 만든 파일은 "${OUT_PREFIX}.all.tsv" = RAW_TSV
# 따라서 mv 필요 없음. (있어도 같은 파일이라 에러남)

if [[ ! -s "${RAW_TSV}" ]]; then
  echo "ERROR: iVar variants output TSV not found or empty: ${RAW_TSV}" >&2
  exit 1
fi

if [[ "${SNP_ONLY}" == "true" ]]; then
  # 헤더 유지 + SNP만 (REF/ALT 길이 1)
  awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1 {print; next}
    {
      ref=$3; alt=$4
      if(length(ref)==1 && length(alt)==1 && ref!="N" && alt!="N") print
    }' "${RAW_TSV}" > "${SNP_TSV}"
else
  cp -f "${RAW_TSV}" "${SNP_TSV}"
fi

