/**
 * ALIGN_VARIANTS
 *
 * Per-gene process (runs ONCE per gene, collecting ALL samples).
 *
 * Steps:
 *   1. Extract the 3D7 reference CDS sequence from the reference genome
 *      (splice_and_translate.py --ref-fasta, Python)
 *   2. Call gene_msa_variants.R which:
 *        a. Reads all per-sample nt/aa FASTAs for this gene
 *        b. Prepends the 3D7 reference
 *        c. Runs MAFFT to produce a multi-sample MSA per gene
 *        d. Produces a variant TSV (3D7 as reference row)
 *
 * Requires: mafft, python3, Rscript (all in PATH)
 *
 * Input:
 *   gene       : gene name (e.g. "CRT")
 *   nt_fastas  : ALL per-sample {barcode}.{GENE}.nt.fasta files (collected)
 *   aa_fastas  : ALL per-sample {barcode}.{GENE}.aa.fasta files (collected)
 *   ref        : reference genome FASTA
 *   cds_tsv    : resources/cds_coords.tsv
 *   outroot    : base output directory
 *
 * Output: {gene}/
 *   {GENE}.nt.combined.fasta    3D7 + all samples (raw)
 *   {GENE}.nt.msa.fasta         MAFFT-aligned nucleotide MSA
 *   {GENE}.nt.variants.tsv      nucleotide variant table
 *   {GENE}.aa.combined.fasta    3D7 + all samples (raw)
 *   {GENE}.aa.msa.fasta         MAFFT-aligned amino-acid MSA
 *   {GENE}.aa.variants.tsv      amino-acid variant table  ← main output
 */
process ALIGN_VARIANTS {

  tag { "${gene}" }

  input:
    tuple val(gene), path(nt_fastas), path(aa_fastas)
    path ref
    path cds_tsv
    val  outroot

  output:
    path("${gene}"), emit: gene_dir

  publishDir "${outroot}/translation", mode: 'copy', overwrite: true

  script:
  """
  set -euo pipefail

  # ── 1. 3D7 reference CDS 서열 추출 ──────────────────────────────────────────
  # splice_and_translate.py 를 --ref-fasta 모드로 실행하여
  # 레퍼런스 게놈에서 exon을 접합·번역한 3D7 기준 서열을 만든다.
  mkdir -p ref_seqs

  python3 ${projectDir}/bin/splice_and_translate.py \\
    --ref-fasta ${ref}                  \\
    --cds       ${cds_tsv}             \\
    --gene      "${gene}"              \\
    --sample    "3D7"                  \\
    --out-nt    "ref_seqs/3D7.${gene}.nt.fasta" \\
    --out-aa    "ref_seqs/3D7.${gene}.aa.fasta"

  # ── 2. EMPTY placeholder 파일을 작업 디렉토리에서 제거 ──────────────────────
  # SPLICE_TRANSLATE 가 모든 유전자가 실패했을 때 생성하는 EMPTY.nt.fasta 제거
  mkdir -p nt_dir aa_dir
  for f in *.nt.fasta; do
    [[ "\$f" == EMPTY* ]] && continue
    cp "\$f" nt_dir/
  done
  for f in *.aa.fasta; do
    [[ "\$f" == EMPTY* ]] && continue
    cp "\$f" aa_dir/
  done

  # ── 3. R 스크립트로 MSA + 변이표 생성 ──────────────────────────────────────
  # gene_msa_variants.R 이 다음을 수행한다:
  #   - nt_dir / aa_dir 에서 이 유전자의 파일만 수집
  #   - 3D7 reference 를 첫 번째 서열로 붙여 combined FASTA 작성
  #   - MAFFT 정렬 수행
  #   - 3D7 대비 variant 위치를 TSV 로 저장
  Rscript ${projectDir}/bin/gene_msa_variants.R \\
    --gene   "${gene}"                           \\
    --nt_dir "nt_dir"                            \\
    --aa_dir "aa_dir"                            \\
    --ref_nt "ref_seqs/3D7.${gene}.nt.fasta"    \\
    --ref_aa "ref_seqs/3D7.${gene}.aa.fasta"    \\
    --outdir "."
  """
}
