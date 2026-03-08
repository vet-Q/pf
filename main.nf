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
 *   nextflow run main.nf -params-file resources/params.yaml
 *   OR modify resources/params.yaml and run: nextflow run main.nf
 */

// ==================== LOAD PARAMETERS FROM YAML ====================
def loadParamsFromYaml(yamlFile) {
    def yamlPath = file(yamlFile)
    if (!yamlPath.exists()) {
        error "❌ ERROR: YAML 파일을 찾을 수 없습니다: ${yamlFile}\n" +
              "   Please check: ${projectDir}/configs/params.yaml"
    }
    
    // YAML 파싱 간단한 버전 (대안: snakeyaml 라이브러리 사용 가능)
    def lines = yamlPath.text.split('\n')
    def config = [:]
    def section = null
    
    lines.each { line ->
        line = line.trim()
        if (line.isEmpty() || line.startsWith('#')) return
        
        if (line.endsWith(':') && !line.startsWith(' ')) {
            section = line.replaceAll(':$', '')
        } else if (line.contains(':')) {
            def parts = line.split(':', 2)
            def key = parts[0].trim()
            def value = parts[1].trim()
            if (section) key = section + '_' + key
            
            // 값 파싱 (문자열, boolean, 숫자)
            if (value.equalsIgnoreCase('true')) config[key] = true
            else if (value.equalsIgnoreCase('false')) config[key] = false
            else if (value.isEmpty()) config[key] = null
            else if (value.matches(/-?\\d+(\\.\\d+)?/)) config[key] = value.toDouble()
            else config[key] = value.replaceAll('^["\']|["\']$', '')
        }
    }
    
    return config
}

// Load YAML configuration
def yamlConfig = loadParamsFromYaml("${projectDir}/configs/params.yaml")
log.info "✓ Loaded parameters from: configs/params.yaml"

// ==================== SET PARAMETERS FROM YAML ====================
// Paths
params.bam_root    = params.bam_root    ?: yamlConfig.get('paths_bam_root', ".")
params.bam_glob    = params.bam_glob    ?: yamlConfig.get('paths_bam_glob', "*.bam")
params.fastq_root  = params.fastq_root  ?: yamlConfig.get('paths_fastq_root', ".")
params.fastq_glob  = params.fastq_glob  ?: yamlConfig.get('paths_fastq_glob', "**/*.{fastq,fq,fastq.gz,fq.gz}")
params.outdir      = params.outdir      ?: yamlConfig.get('paths_outdir', "./results")

// References
params.bed         = params.bed         ?: "${projectDir}/resources/" + yamlConfig.get('references_bed', "amplicons.bed")
params.ref         = params.ref         ?: "${projectDir}/resources/" + yamlConfig.get('references_ref', "PlasmoDB-67_Pfalciparum3D7_Genome.fasta")

// iVar parameters
params.ivar_threshold    = params.ivar_threshold    ?: yamlConfig.get('ivar_threshold', 0.6)
params.ivar_min_cons_cov = params.ivar_min_cons_cov ?: yamlConfig.get('ivar_min_cons_cov', 10)
params.ivar_min_var_cov  = params.ivar_min_var_cov  ?: yamlConfig.get('ivar_min_var_cov', 1)
params.genes             = params.genes             ?: yamlConfig.get('ivar_genes', "CRT DHFR DHPS K13 MDR1 MSP2 PMI PMIII")

// Minimap2 parameters
params.minimap2_preset = params.minimap2_preset ?: yamlConfig.get('minimap2_preset', "sr")
params.minimap2_threads = params.minimap2_threads ?: yamlConfig.get('minimap2_threads', 4)

// Bcftools parameters
params.variant_method   = params.variant_method   ?: "bcftools"
params.bcftools_mode    = params.bcftools_mode    ?: yamlConfig.get('bcftools_mode', "mv")
params.snp_only         = params.snp_only != null ? params.snp_only : yamlConfig.get('bcftools_snp_only', true)
params.bcftools_threads = params.bcftools_threads ?: yamlConfig.get('bcftools_threads', 4)
params.max_depth        = params.max_depth        ?: yamlConfig.get('bcftools_max_depth', 10000)

// AF/DP Figure parameters
params.afdp_min_depth    = params.afdp_min_depth    ?: yamlConfig.get('afdp_figures_min_depth', 15)
params.afdp_af_threshold = params.afdp_af_threshold ?: yamlConfig.get('afdp_figures_af_threshold', 0.50)
params.afdp_shade_xmin   = params.afdp_shade_xmin   ?: yamlConfig.get('afdp_figures_shade_xmin', 0.48)
params.afdp_shade_xmax   = params.afdp_shade_xmax   ?: yamlConfig.get('afdp_figures_shade_xmax', 0.52)

// VCF pattern
params.vcf_glob = params.vcf_glob ?: "${params.outdir}/calling/bcftools/**/*.vcf.gz"

// Analysis toggles (pipeline steps)
params.do_fastqc     = params.do_fastqc     != null ? params.do_fastqc     : yamlConfig.get('pipeline_do_fastqc', false)
params.do_alignment  = params.do_alignment  != null ? params.do_alignment  : yamlConfig.get('pipeline_do_alignment', false)
params.do_depth      = params.do_depth      != null ? params.do_depth      : yamlConfig.get('pipeline_do_depth', true)
params.do_ivar       = params.do_ivar       != null ? params.do_ivar       : yamlConfig.get('pipeline_do_ivar', true)
params.do_gene_bins  = params.do_gene_bins  != null ? params.do_gene_bins  : yamlConfig.get('pipeline_do_gene_bins', true)
params.do_variants   = params.do_variants   != null ? params.do_variants   : yamlConfig.get('pipeline_do_variants', true)
params.do_multipanel = params.do_multipanel != null ? params.do_multipanel : yamlConfig.get('pipeline_do_multipanel', false)
params.do_gene_1x8   = params.do_gene_1x8   != null ? params.do_gene_1x8   : yamlConfig.get('pipeline_do_gene_1x8', false)
params.do_afdp       = params.do_afdp       != null ? params.do_afdp       : yamlConfig.get('pipeline_do_afdp', true)

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

  // Load BAM files OR align from FASTQ
  if (params.do_alignment || params.do_fastqc) {
    // Load FASTQ files if alignment is requested
    Channel
      .fromPath("${params.fastq_root}/${params.fastq_glob}", checkIfExists: true)
      .ifEmpty { error "❌ ERROR: No FASTQ files found at: ${params.fastq_root}/${params.fastq_glob}\n" +
                       "   Please check params.fastq_root and params.fastq_glob settings." }
      .map { fastq ->
        def sample = fastq.getBaseName().replaceAll(/\.(fastq|fq|fastq\.gz|fq\.gz).*/, "")
        tuple(sample, fastq)
      }
      .groupTuple()
      .set { fastq_ch }
  } else {
    // Load pre-aligned BAM files
    Channel
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
  bed_ch = Channel.value(file(params.bed))
  ref_ch = Channel.value(file(params.ref))
  
  // Load R scripts for visualization
  r_afdp = Channel.value(file("${projectDir}/bin/make_af_dp_figures.R"))
  r_depth = Channel.value(file("${projectDir}/bin/plot_depth_per_gene.R"))
  r_multipanel = Channel.value(file("${projectDir}/bin/plot_gene_multipanel.R"))
  r_gene1x8 = Channel.value(file("${projectDir}/bin/plot_gene_1x8.R"))

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
      ref_ch
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
      r_afdp,
      params.afdp_min_depth,
      params.afdp_af_threshold,
      params.afdp_shade_xmin,
      params.afdp_shade_xmax
    )
  }

  log.info "Pipeline completed successfully!"
}
