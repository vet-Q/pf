/**
 * SPLICE_TRANSLATE
 *
 * Per-sample process: reads per-gene consensus FASTAs produced by CALL_VARIANTS,
 * splices out intronic regions using CDS coordinates from cds_coords.tsv,
 * and translates to protein.
 *
 * Input:
 *   barcode      : sample name (may include suffixes like ".Pf3D7" from BAM filename)
 *   variant_dir  : per-barcode directory from CALL_VARIANTS (contains gene subdirs)
 *   cds_tsv      : resources/cds_coords.tsv
 *   genes_str    : space-separated gene list (e.g. "CRT DHFR DHPS K13 MDR1 MSP2 PMI PMIII")
 *   outroot      : base output directory
 *   method       : variant caller used ("bcftools" | "ivar")
 *   af_threshold : AF threshold for IUPAC resolution (e.g. 0.6)
 *
 * Output per sample (clean barcode label, e.g. barcode01):
 *   {clean_barcode}.{GENE}.nt.fasta   spliced nucleotide CDS
 *   {clean_barcode}.{GENE}.aa.fasta   translated protein
 *   {clean_barcode}.{GENE}.af.tsv     IUPAC resolution table (AF per ambiguous CDS site)
 *
 * Note: barcode is cleaned to the first dot-separated field so that
 *   "barcode09.Pf3D7" → "barcode09"
 *   This ensures gene_msa_variants.R can parse barcode names from filenames.
 */
process SPLICE_TRANSLATE {

  tag { "${barcode}" }

  input:
    tuple val(barcode), path(variant_dir)
    path cds_tsv
    val  genes_str
    val  outroot
    val  method
    val  af_threshold

  output:
    path("*.nt.fasta"), emit: nt_fastas
    path("*.aa.fasta"), emit: aa_fastas
    path("*.af.tsv"),   emit: af_tsvs

  publishDir "${outroot}/translation/per_sample", mode: 'copy', overwrite: true

  script:
  """
  set -euo pipefail

  # barcode 이름에서 첫 번째 점 이전 부분만 사용
  # (예: barcode09.Pf3D7 → barcode09)
  CLEAN_BARCODE="\$(echo '${barcode}' | cut -d. -f1)"

  GENES="${genes_str}"

  for GENE in \$GENES; do
    # CALL_VARIANTS 가 만든 consensus FASTA 경로
    FA="${variant_dir}/\${GENE}/${barcode}.\${GENE}.${method}.consensus.fasta"
    # 같은 레이어의 bcftools VCF (IUPAC 해소용 AF 정보)
    VCF="${variant_dir}/\${GENE}/${barcode}.\${GENE}.bcftools.vcf.gz"
    if [ ! -f "\$FA" ]; then
      echo "  [SPLICE_TRANSLATE] WARNING: \$FA not found, skipping \${GENE}" >&2
      continue
    fi

    python3 ${projectDir}/bin/splice_and_translate.py \\
      --fasta        "\$FA"                                          \\
      --vcf          "\$VCF"                                         \\
      --af-threshold ${af_threshold}                                  \\
      --cds          ${cds_tsv}                                      \\
      --gene         "\$GENE"                                        \\
      --sample       "\${CLEAN_BARCODE}"                             \\
      --out-nt       "\${CLEAN_BARCODE}.\${GENE}.nt.fasta"          \\
      --out-aa       "\${CLEAN_BARCODE}.\${GENE}.aa.fasta"          \\
      --out-af       "\${CLEAN_BARCODE}.\${GENE}.af.tsv"
  done

  # 출력 파일이 하나도 없을 경우 placeholder 생성 (Nextflow output 검증 통과용)
  ls ./*.nt.fasta 2>/dev/null || touch "EMPTY.nt.fasta"
  ls ./*.aa.fasta 2>/dev/null || touch "EMPTY.aa.fasta"
  ls ./*.af.tsv   2>/dev/null || touch "EMPTY.af.tsv"
  """
}
