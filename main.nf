nextflow.enable.dsl = 2

/**
 * Plasmodium falciparum (Pf) Analysis Pipeline
 * 
 * Purpose: Process aligned/sorted BAM files to generate variant calls, 
 *          consensus sequences, depth analysis, and gene-level visualizations
 * 
 * Input: Pre-aligned and sorted BAM files (BAM/BAI format) or FASTQ files
 * Output: VCF files, consensus sequences, depth plots, and AF/DP figures
 * 
 * Usage:
 *   nextflow run main.nf -profile standard
 *   nextflow run main.nf -profile standard --do_alignment true --fastq_root ./data
 *   All parameters are configured in nextflow.config (params block)
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

// ==================== MAIN WORKFLOW ====================
workflow {

  // ==================== RESOLVE PATHS ====================
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

  // Load BAM files OR align from FASTQ
  if (params.do_alignment || params.do_fastqc) {
    // Load FASTQ files if alignment is requested
    channel
      .fromPath("${params.fastq_root}/${params.fastq_glob}", checkIfExists: true)
      .ifEmpty { error "❌ ERROR: No FASTQ files found at: ${params.fastq_root}/${params.fastq_glob}\n" +
                       "   Please check params.fastq_root and params.fastq_glob settings." }
      .map { fastq ->
        def sample = fastq.parent.name
        tuple(sample, fastq)
      }
      .groupTuple()
      .set { fastq_ch }
  } else {
    // Load pre-aligned BAM files
    channel
      .fromPath("${params.bam_root}/${params.bam_glob}", checkIfExists: true)
      .ifEmpty { error "❌ ERROR: No BAM files found at: ${params.bam_root}/${params.bam_glob}\n" +
                       "   Please check params.bam_root and params.bam_glob settings." }
      .map { bam ->
        def barcode = bam.getBaseName().replaceAll(/\.sorted$/, "")
        tuple(barcode, bam)
      }
      .set { bam_ch }
  }

  // Load reference resources
  bed_ch = channel.value(file(bed_abs))
  ref_ch = channel.value(file(ref_abs))
  
  // Load R scripts for visualization
  r_afdp = channel.value(file("${projectDir}/bin/make_af_dp_figures.R"))
  r_depth = channel.value(file("${projectDir}/bin/plot_depth_per_gene.R"))
  r_multipanel = channel.value(file("${projectDir}/bin/plot_gene_multipanel.R"))
  r_gene1x8 = channel.value(file("${projectDir}/bin/plot_gene_1x8.R"))

  // ==================== FASTQ QUALITY CONTROL ====================
  if (params.do_fastqc) {
    log.info "Running: FastQ Quality Control"
    FASTQ_QC(fastq_ch)
    fastq_qc_pass = FASTQ_QC.out.fastq_pass
  } else if (params.do_alignment) {
    fastq_qc_pass = fastq_ch
  }

  // ==================== READ ALIGNMENT (MINIMAP2) ====================
  if (params.do_alignment) {
    log.info "Running: Minimap2 Alignment (preset: ${params.minimap2_preset})"
    
    MINIMAP2_ALIGN(
      fastq_qc_pass, 
      ref_ch,
      params.minimap2_preset
    )
    
    // Output from minimap2 becomes our BAM channel
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
  // Note: gene_bins output is required by multipanel and gene_1x8 plots
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
    if (gene_bins_ch == null) {
      error "❌ ERROR: do_multipanel=true requires gene_bins output.\n" +
            "   Please enable: do_gene_bins=true OR do_multipanel=true OR do_gene_1x8=true"
    }
    PLOT_GENE_MULTIPANEL(gene_bins_ch, r_multipanel)
  }

  // ==================== GENE 1x8 PLOT ====================
  if (params.do_gene_1x8) {
    log.info "Running: Gene 1x8 Plot"
    if (gene_bins_ch == null) {
      error "❌ ERROR: do_gene_1x8=true requires gene_bins output.\n" +
            "   Please enable: do_gene_bins=true OR do_multipanel=true OR do_gene_1x8=true"
    }
    def genes_csv = params.genes.tokenize(' ').join(',')
    PLOT_GENE_1x8(gene_bins_ch, genes_csv, r_gene1x8)
  }

  // ==================== AF/DP FIGURES ====================
  if (params.do_afdp) {
    log.info "Running: AF/DP Figures"
    
    if (vcf_ch == null) {
      error "❌ ERROR: do_afdp=true requires VCF files.\n" +
            "   Please enable: do_variants=true"
    }

    vcf_ch
      .flatMap { dir -> file("${dir}/**/*.vcf.gz") }
      .collect()
      .set { all_vcfs_ch }

    AF_DP_FIGURES(
      all_vcfs_ch, 
      r_afdp
    )
  }

  log.info "Pipeline completed successfully!"
}
