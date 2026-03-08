process IVAR_CONSENSUS {

  tag { barcode }

  publishDir "${outdir}/calling/ivar", mode: 'copy', overwrite: true

  input:
    tuple val(barcode), path(bam)
    path bed
    path ref
    val genes_str
    val outdir
    val ivar_t
    val ivar_m_cons
    val ivar_m_var

  output:
    path("${barcode}")

  script:
  """
  set -euo pipefail

  mkdir -p ${barcode}
  GENES="${genes_str}"

  for GENE in \$GENES; do
    echo "▶ ${barcode} - \$GENE"
    mkdir -p ${barcode}/\$GENE

    # gene bed (chr start end) 3컬럼
    awk -v g="\$GENE" '\$4==g{print \$1"\\t"\$2"\\t"\$3}' ${bed} > ${barcode}/\$GENE/region.bed

    # per-position depth
    samtools depth -aa -b ${barcode}/\$GENE/region.bed ${bam} > ${barcode}/\$GENE/depth.tsv

    # (1) iVar consensus
    ${projectDir}/bin/make_consensus_ivar.sh \\
      ${bam} \\
      ${ref} \\
      ${barcode}/\$GENE/region.bed \\
      ${barcode}/\$GENE/${barcode}_\${GENE}_ivar \\
      ${ivar_t} \\
      ${ivar_m_cons} \\
      ${barcode}/\$GENE/${barcode}.\$GENE.ivar.consensus.fasta

    # (2) iVar variants (SNP only)
    ${projectDir}/bin/call_ivar.sh \\
      ${bam} \\
      ${ref} \\
      ${barcode}/\$GENE/region.bed \\
      ${barcode}/\$GENE/${barcode}_\${GENE}_ivar \\
      0 \\
      ${ivar_m_var} \\
      true

    # 산출물 정리(이름 통일)
    # RAW: ${barcode}_\${GENE}_ivar.all.tsv
    # SNP: ${barcode}_\${GENE}_ivar.snp.tsv
  done
  """
}

