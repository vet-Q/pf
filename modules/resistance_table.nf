process RESISTANCE_TABLE {

  publishDir "${outdir}/resistance", mode: 'copy', overwrite: true

  input:
    path variants_tsv   // variant_AF_all.tsv (all samples merged)
    path cds_tsv        // resources/cds_coords.tsv
    path r_script       // bin/make_resistance_table.R
    val  af_threshold   // e.g. 0.5
    val  outdir

  output:
    path "resistance_markers.tsv"

  script:
  """
  Rscript ${r_script} \
    ${variants_tsv} \
    ${cds_tsv} \
    resistance_markers.tsv \
    ${af_threshold}
  """
}
