process PLOT_GENE_1x8 {

  tag { barcode }
  publishDir { "${outdir}/gene_bins" }, mode: 'copy', overwrite: true

  input:
    tuple val(barcode), path(barcode_dir), val(outdir)
    val genes_csv
    path plot_script

  output:
    path("${barcode}_genes_1x8.pdf")

  script:
  """
  set -euo pipefail
  Rscript ${plot_script} ${barcode_dir} "${genes_csv}" ${barcode}_genes_1x8.pdf
  """
}
