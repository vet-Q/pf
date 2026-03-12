nextflow.enable.dsl = 2

/**
 * Plasmodium falciparum (Pf) Analysis Pipeline
 *
 * Input : FASTQ (Nanopore) or pre-aligned BAM files
 * Output: VCF, consensus FASTA, depth plots, AF/DP figures,
 *         gene-level MSA + amino-acid variant tables
 *
 * Usage:
 *   nextflow run main.nf -profile standard
 *   nextflow run main.nf -profile standard --do_alignment true --fastq_root ./data
 *   All parameters are defined in nextflow.config (params block).
 */

// ==================== MODULE INCLUDES ====================
include { FASTQ_QC }             from "./modules/fastq_qc.nf"
include { MINIMAP2_ALIGN }       from "./modules/minimap2_align.nf"
include { PF_DEPTH }             from "./modules/pf_depth.nf"
include { IVAR_CONSENSUS }       from "./modules/ivar_consensus.nf"
include { GENE_BIN_DEPTH }       from "./modules/gene_bins_nf.nf"
include { CALL_VARIANTS }        from "./modules/call_variants.nf"
include { PLOT_GENE_MULTIPANEL } from "./modules/gene_bins_plot_nf.nf"
include { PLOT_GENE_1x8 }        from "./modules/gene_1to8_plot.nf"
include { AF_DP_FIGURES }        from "./modules/af_dp_figures.nf"
include { SPLICE_TRANSLATE }     from "./modules/splice_translate.nf"
include { ALIGN_VARIANTS }       from "./modules/align_variants.nf"

// ==================== MAIN WORKFLOW ====================
workflow {

  // ── 경로 확정 (null 이면 resources/ 기본값 사용) ──────────────────────────
  def bed_abs = params.bed ?: "${projectDir}/resources/amplicons.bed"
  def ref_abs = params.ref ?: "${projectDir}/resources/PlasmoDB-67_Pfalciparum3D7_Genome.fasta"

  log.info """
  ========================================
    pfNextflow Pipeline
  ========================================
    FASTQ root  : ${params.fastq_root}
    BAM root    : ${params.bam_root}
    Output dir  : ${params.outdir}
    Reference   : ${ref_abs}
    BED file    : ${bed_abs}
    Alignment   : ${params.do_alignment}
    Preset      : ${params.minimap2_preset}
  ========================================
  """

  // ── 입력 채널: FASTQ 정렬 또는 기존 BAM 로드 ────────────────────────────
  if (params.do_alignment || params.do_fastqc) {
    channel
      .fromPath("${params.fastq_root}/${params.fastq_glob}", checkIfExists: true)
      .ifEmpty { error "❌ No FASTQ files at: ${params.fastq_root}/${params.fastq_glob}" }
      .map     { fq -> tuple(fq.parent.name, fq) }
      .groupTuple()
      .set { fastq_ch }
  } else if (params.do_depth || params.do_ivar || params.do_gene_bins || params.do_variants) {
    // BAM 이 실제로 필요한 경우에만 로드
    channel
      .fromPath("${params.bam_root}/${params.bam_glob}", checkIfExists: true)
      .ifEmpty { error "❌ No BAM files at: ${params.bam_root}/${params.bam_glob}" }
      .map     { bam -> tuple(bam.getBaseName().replaceAll(/\.sorted$/, ""), bam) }
      .set { bam_ch }
  }

  // ── 공통 리소스 채널 ─────────────────────────────────────────────────────
  bed_ch = channel.value(file(bed_abs))
  ref_ch = channel.value(file(ref_abs))

  // ── R 시각화 스크립트 채널 ───────────────────────────────────────────────
  r_afdp      = channel.value(file("${projectDir}/bin/make_af_dp_figures.R"))
  r_depth     = channel.value(file("${projectDir}/bin/plot_depth_per_gene.R"))
  r_multipanel = channel.value(file("${projectDir}/bin/plot_gene_multipanel.R"))
  r_gene1x8   = channel.value(file("${projectDir}/bin/plot_gene_1x8.R"))

  // ==================== FASTQ QC ====================
  if (params.do_fastqc) {
    log.info "Running: FastQ Quality Control"
    FASTQ_QC(fastq_ch)
    fastq_qc_pass = FASTQ_QC.out.fastq_pass
  } else if (params.do_alignment) {
    fastq_qc_pass = fastq_ch
  }

  // ==================== ALIGNMENT (MINIMAP2) ====================
  if (params.do_alignment) {
    log.info "Running: Minimap2 Alignment (preset: ${params.minimap2_preset})"
    MINIMAP2_ALIGN(fastq_qc_pass, ref_ch, params.minimap2_preset)
    bam_ch = MINIMAP2_ALIGN.out.bam
  }

  // ==================== DEPTH ANALYSIS ====================
  if (params.do_depth) {
    log.info "Running: Depth Analysis"
    PF_DEPTH(bam_ch, bed_ch, r_depth, params.outdir)
  }

  // ==================== IVAR CONSENSUS ====================
  if (params.do_ivar) {
    log.info "Running: iVar Consensus Generation"
    IVAR_CONSENSUS(
      bam_ch, bed_ch, ref_ch, params.genes,
      params.outdir,
      params.ivar_threshold,
      params.ivar_min_cons_cov,
      params.ivar_min_var_cov
    )
  }

  // ==================== GENE BIN DEPTH ====================
  def gene_bins_ch = null
  if (params.do_gene_bins || params.do_multipanel || params.do_gene_1x8) {
    log.info "Running: Gene Bin Depth Analysis"
    gene_bins_ch = GENE_BIN_DEPTH(bam_ch, bed_ch, params.outdir)
  }

  // ==================== VARIANT CALLING ====================
  def vcf_ch = null
  if (params.do_variants) {
    log.info "Running: Variant Calling (${params.variant_method})"
    vcf_ch = CALL_VARIANTS(
      bam_ch, bed_ch, ref_ch, params.genes,
      params.outdir,
      params.variant_method,
      params.bcftools_mode,
      params.snp_only,
      params.bcftools_threads,
      params.max_depth
    )
  }

  // ==================== GENE MULTIPANEL PLOT ====================
  if (params.do_multipanel) {
    log.info "Running: Gene Multipanel Plot"
    if (gene_bins_ch == null) error "❌ do_multipanel requires do_gene_bins=true"
    PLOT_GENE_MULTIPANEL(gene_bins_ch, r_multipanel)
  }

  // ==================== GENE 1x8 PLOT ====================
  if (params.do_gene_1x8) {
    log.info "Running: Gene 1x8 Plot"
    if (gene_bins_ch == null) error "❌ do_gene_1x8 requires do_gene_bins=true"
    PLOT_GENE_1x8(gene_bins_ch, params.genes.tokenize(' ').join(','), r_gene1x8)
  }

  // ==================== AF/DP FIGURES ====================
  if (params.do_afdp) {
    log.info "Running: AF/DP Figures"
    if (vcf_ch == null) error "❌ do_afdp requires do_variants=true"

    vcf_ch
      .flatMap { dir -> file("${dir}/**/*.vcf.gz") }
      .collect()
      .set { all_vcfs_ch }

    AF_DP_FIGURES(all_vcfs_ch, r_afdp)
  }

  // ==================== CDS SPLICE + TRANSLATE + VARIANT TABLE ====================
  if (params.do_translate) {
    log.info "Running: CDS Splice, Translate & Variant Table"

    def cds_tsv_ch = channel.value(file("${projectDir}/resources/cds_coords.tsv"))

    // vcf_ch 가 없으면 (단독 do_translate 실행) 디스크의 기존 calling 결과에서 로드
    if (vcf_ch == null) {
      log.info "  기존 calling 결과를 디스크에서 로드: ${params.outdir}/calling/${params.variant_method}/"
      vcf_ch = channel
        .fromPath("${params.outdir}/calling/${params.variant_method}/*", type: 'dir', checkIfExists: true)
        .ifEmpty { error "❌ do_translate 에는 do_variants=true 또는 ${params.outdir}/calling/${params.variant_method}/ 하위 디렉토리가 필요합니다." }
    }

    // 샘플별 변이 디렉토리에 barcode 이름 부여
    vcf_ch
      .map { dir -> tuple(dir.name, dir) }
      .set { variants_with_barcode_ch }

    // 샘플별: consensus FASTA → exon 접합 → 번역
    SPLICE_TRANSLATE(
      variants_with_barcode_ch,
      cds_tsv_ch,
      params.genes,
      params.outdir,
      params.variant_method,
      params.ivar_threshold   // IUPAC 해소 AF threshold (bcftools consensus 결과용)
    )

    // NT: 유전자별로 전체 샘플 모음 (barcode01.CRT.nt.fasta → gene key = CRT)
    SPLICE_TRANSLATE.out.nt_fastas
      .flatten()
      .filter { f -> !f.name.startsWith('EMPTY') }
      .map    { f -> tuple(f.name.tokenize('.')[1], f) }
      .groupTuple()
      .set { nt_by_gene_ch }

    // AA: 동일 구조
    SPLICE_TRANSLATE.out.aa_fastas
      .flatten()
      .filter { f -> !f.name.startsWith('EMPTY') }
      .map    { f -> tuple(f.name.tokenize('.')[1], f) }
      .groupTuple()
      .set { aa_by_gene_ch }

    // 유전자별 NT+AA 합쳐서 MAFFT 정렬 + 변이표 생성
    nt_by_gene_ch
      .join(aa_by_gene_ch)
      .set { by_gene_ch }

    ALIGN_VARIANTS(by_gene_ch, ref_ch, cds_tsv_ch, params.outdir)
  }

  log.info "Pipeline completed successfully!"
}

