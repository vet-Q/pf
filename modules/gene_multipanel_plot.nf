process PLOT_GENE_MULTIPANEL {

  tag { barcode }

  publishDir { "${outdir}/gene_bins" }, mode: 'copy', overwrite: true

  input:
    tuple val(barcode), path(barcode_dir), val(outdir)
    path plot_script

  output:
    path("${barcode}_gene_multipanel.pdf")

  script:
  """
  set -euo pipefail
  Rscript ${plot_script} ${barcode_dir} ${barcode}_gene_multipanel.pdf
  """
}

