process PF_DEPTH {

  tag { barcode }

  input:
    tuple val(barcode), path(bam)
    path bed
    path plot_script
    val outroot

  output:
    path("${barcode}")

  publishDir { "${outroot}/pf_depth" }, mode: 'copy', overwrite: true

  script:
  """
  set -euo pipefail

  mkdir -p ${barcode}/depth_by_gene
  mkdir -p ${barcode}/plots_by_gene

  samtools depth -b ${bed} ${bam} > ${barcode}/all_amplicons.depth.txt

  Rscript ${plot_script} \
    ${barcode}/all_amplicons.depth.txt \
    ${bed} \
    ${barcode}
  """
}

