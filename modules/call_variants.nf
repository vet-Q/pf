process CALL_VARIANTS {

  tag { "${barcode}" }

  input:
    tuple val(barcode), path(bam)
    path bed
    path ref
    val genes_str
    val outroot

    // from main.nf
    val method        // "bcftools" | "ivar" | "both" (지금은 bcftools만 써도 됨)
    val bcft_mode     // "mv" | "cv"
    val snp_only      // true | false
    val threads       // int
    val max_depth     // int

  output:
    path("${barcode}")

  // ✅ 결과를 outdir/results 아래로 통일
  publishDir { "${outroot}/calling/${method}" }, mode: 'copy', overwrite: true


  script:
  """
  set -euo pipefail

  mkdir -p ${barcode}
  GENES="${genes_str}"

  for GENE in \$GENES; do
    mkdir -p ${barcode}/\$GENE

    awk -v g="\$GENE" '\$4==g{print \$1"\\t"\$2"\\t"\$3}' ${bed} > ${barcode}/\$GENE/region.bed

    # ✅ gene별 BAM (필요하면 유지)
    samtools view -b -L ${barcode}/\$GENE/region.bed ${bam} > ${barcode}/\$GENE/${barcode}.\$GENE.bam
    samtools index ${barcode}/\$GENE/${barcode}.\$GENE.bam

    if [[ "${method}" == "bcftools" || "${method}" == "both" ]]; then
      # 네가 만들어둔 스크립트 사용 (이스케이프 최소)
      ${projectDir}/bin/call_bcftools.sh \
        ${barcode}/\$GENE/${barcode}.\$GENE.bam \
        ${ref} \
        ${barcode}/\$GENE/region.bed \
        ${barcode}/\$GENE/${barcode}.\$GENE.bcftools.vcf.gz \
        ${bcft_mode} \
        ${snp_only} \
        ${max_depth}

      ${projectDir}/bin/make_consensus_bcftools.sh \
        ${barcode}/\$GENE/${barcode}.\$GENE.bcftools.vcf.gz \
        ${ref} \
        ${barcode}/\$GENE/region.bed \
        ${barcode}/\$GENE/${barcode}.\$GENE.bcftools.consensus.fasta \
        ${barcode}/\$GENE/${barcode}.\$GENE
    fi

    # (optional) ivar도 같이 돌릴 거면 여기 블록 추가 가능
  done
  """
}
