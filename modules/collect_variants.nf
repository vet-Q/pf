process COLLECT_VARIANTS {

  tag "COLLECT_VARIANTS"

  publishDir "${outdir}/calling", mode: 'copy', overwrite: true

  input:
    path vcf_files   // collect()로 모은 모든 .vcf.gz 파일들
    path r_script    // bin/collect_variants.R
    val  outdir

  output:
    path "variant_AF_all.tsv", emit: variants_tsv

  script:
  """
  ls -1 *.vcf.gz > vcf_list.txt

  Rscript ${r_script} \
    --vcf_list vcf_list.txt \
    --out_tsv  variant_AF_all.tsv \
    --genes    "${params.genes.tokenize(' ').join(',')}"
  """
}
