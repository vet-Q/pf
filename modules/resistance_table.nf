process RESISTANCE_TABLE {

  publishDir "${outdir}/resistance", mode: 'copy', overwrite: true

  input:
    path variants_tsv   // variant_AF_all.tsv (all samples merged)
    path cds_tsv        // resources/cds_coords.tsv
    path r_script       // bin/make_resistance_table.R
    val  af_threshold   // e.g. 0.5 — 이 이상이면 변이 콜로 표시
    val  dp_min         // e.g. 20 — 이 미만이면 lowcov(DP)로 표시
    val  outdir

  output:
    path "resistance_markers.tsv"

  script:
  """
  Rscript ${r_script} \
    ${variants_tsv} \
    ${cds_tsv} \
    resistance_markers.tsv \
    ${af_threshold} \
    ${dp_min}
  """
}
