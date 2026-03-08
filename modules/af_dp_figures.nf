process AF_DP_FIGURES {

  tag "AFDP"

  publishDir "${params.outdir}/figures/afdp", mode: 'copy'

  input:
    path vcf_files
    path afdp_r

  output:
    path "Figure1_by_gene.pdf"
    path "Figure2_overall.pdf"

  script:
  """
  # vcf_files는 Nextflow가 stage-in한 파일들이고,
  # R은 그 목록만 있으면 됨.
  ls -1 *.vcf.gz > vcf_list.txt

  Rscript ${afdp_r} \
    --vcf_list vcf_list.txt \
    --out_fig1 Figure1_by_gene.pdf \
    --out_fig2 Figure2_overall.pdf \
    --genes "${params.genes.tokenize(' ').join(',')}" \
    --min_dp ${params.afdp_min_dp} \
    --af_cut ${params.afdp_af_cut} \
    --shade_xmin ${params.afdp_shade_xmin} \
    --shade_xmax ${params.afdp_shade_xmax}
  """
}
