#!/usr/bin/env Rscript
# =============================================================================
# gene_msa_variants.R
#
# 목적:
#   단일 유전자(GENE)에 대해 여러 샘플의 spliced CDS FASTA 파일을 모아서
#   ① 3D7 reference 서열을 맨 앞에 붙인 합본 FASTA를 만들고
#   ② MAFFT로 다중서열정렬(MSA)을 수행하고
#   ③ 3D7 대비 변이 위치만 추출한 TSV 변이표를 생성한다.
#   NT(핵산)과 AA(아미노산) 두 가지 모두 동일 흐름으로 처리한다.
#
# 입력:
#   --gene      유전자 이름 (예: CRT)
#   --nt_dir    샘플별 *.nt.fasta 파일들이 들어있는 디렉토리
#   --aa_dir    샘플별 *.aa.fasta 파일들이 들어있는 디렉토리
#   --ref_nt    3D7 reference 핵산 FASTA (splice_and_translate.py가 생성)
#   --ref_aa    3D7 reference 아미노산 FASTA
#   --outdir    결과를 저장할 디렉토리
#
# 출력: outdir/{GENE}/
#   {GENE}.nt.combined.fasta   3D7 + 전체 샘플 NT 합본
#   {GENE}.nt.msa.fasta        MAFFT NT 정렬 결과
#   {GENE}.nt.variants.tsv     NT 변이표
#   {GENE}.aa.combined.fasta   3D7 + 전체 샘플 AA 합본
#   {GENE}.aa.msa.fasta        MAFFT AA 정렬 결과
#   {GENE}.aa.variants.tsv     AA 변이표  ← 핵심 결과
#
# 의존성:
#   - mafft (PATH에 있어야 함)
#   - R 패키지: optparse (CRAN)
#     * Biostrings 없이 base R + 직접 파싱으로 구현
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(optparse)
})


# ── 0. 커맨드라인 인수 정의 ─────────────────────────────────────────────────
option_list <- list(
  make_option("--gene",   type = "character", help = "유전자 이름 (예: CRT)"),
  make_option("--nt_dir", type = "character", help = "샘플 NT FASTA 디렉토리"),
  make_option("--aa_dir", type = "character", help = "샘플 AA FASTA 디렉토리"),
  make_option("--ref_nt", type = "character", help = "3D7 reference NT FASTA"),
  make_option("--ref_aa", type = "character", help = "3D7 reference AA FASTA"),
  make_option("--outdir", type = "character", default = ".", help = "결과 디렉토리")
)

opt <- parse_args(OptionParser(option_list = option_list))

# 필수 인수 확인
for (arg in c("gene", "nt_dir", "aa_dir", "ref_nt", "ref_aa")) {
  if (is.null(opt[[arg]])) stop(paste("필수 인수 누락:", arg))
}

GENE   <- opt$gene
OUTDIR <- file.path(opt$outdir, GENE)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

message(sprintf("[%s] 유전자별 MSA + 변이표 생성 시작", GENE))


# ── 1. FASTA 읽기 유틸 함수 ─────────────────────────────────────────────────
# FASTA 파일을 읽어서 named character vector (이름=header, 값=서열)로 반환
# 여러 레코드가 있을 수 있으므로 전부 읽음
read_fasta <- function(path) {
  if (!file.exists(path)) stop(paste("파일 없음:", path))
  lines <- readLines(path, warn = FALSE)

  headers  <- c()
  seqs_raw <- list()
  current  <- NULL

  for (line in lines) {
    line <- trimws(line)
    if (nchar(line) == 0) next   # 빈 줄 무시

    if (startsWith(line, ">")) {
      # '>' 이후 부분이 헤더 (공백 이후 설명 포함)
      current <- substr(line, 2, nchar(line))
      headers  <- c(headers, current)
      seqs_raw[[current]] <- character(0)
    } else {
      # 서열 라인: 대문자로 통일하고 기존 서열에 이어 붙임
      seqs_raw[[current]] <- c(seqs_raw[[current]], toupper(line))
    }
  }

  # 각 레코드 서열을 하나의 문자열로 합침
  seqs <- sapply(headers, function(h) paste(seqs_raw[[h]], collapse = ""))
  return(seqs)
}


# ── 2. FASTA 쓰기 유틸 함수 ─────────────────────────────────────────────────
# named character vector → FASTA 파일
# line_width: 한 줄에 쓸 문자 수 (0 이면 줄바꿈 없이 한 줄)
write_fasta <- function(seqs, path, line_width = 80) {
  con <- file(path, "w")
  for (nm in names(seqs)) {
    writeLines(paste0(">", nm), con)
    seq_str <- seqs[[nm]]
    if (line_width > 0 && nchar(seq_str) > line_width) {
      # 지정된 너비로 래핑
      chunks <- sapply(
        seq(1, nchar(seq_str), by = line_width),
        function(i) substr(seq_str, i, min(i + line_width - 1, nchar(seq_str)))
      )
      writeLines(chunks, con)
    } else {
      writeLines(seq_str, con)
    }
  }
  close(con)
}


# ── 3. MAFFT 실행 유틸 함수 ─────────────────────────────────────────────────
# input_fasta → output_fasta
# 서열 수가 많을 때는 --auto 가 적합한 알고리즘을 자동 선택
run_mafft <- function(input_fasta, output_fasta) {
  cmd <- sprintf("mafft --auto --quiet '%s' > '%s'", input_fasta, output_fasta)
  message(sprintf("  MAFFT 실행: %s", basename(input_fasta)))
  ret <- system(cmd)
  if (ret != 0) stop(sprintf("MAFFT 실패 (exit code %d): %s", ret, cmd))
  if (!file.exists(output_fasta) || file.size(output_fasta) == 0) {
    stop("MAFFT 출력 파일이 비어있거나 생성되지 않았음: ", output_fasta)
  }
}


# ── 4. 변이표 생성 유틸 함수 ─────────────────────────────────────────────────
# MSA FASTA (named vector) + reference 헤더 식별자 → 변이 TSV data.frame
#
# 동작 방식:
#   - reference 서열에서 gap(-)인 열(column)은 insertion → 건너뜀
#   - reference 기준 1-based 위치를 pos 로 기록
#   - 그 열에서 하나 이상의 샘플이 reference와 다른 문자를 가지면 variant row
#   - gap(-), N, X, ? 는 "불확정"으로 취급해 variant 판정에서 제외
#     (단, 변이표에는 실제 문자 그대로 기록)
make_variant_table <- function(msa_seqs, ref_tag, gene) {
  # ref_tag 를 포함하는 첫 번째 헤더를 reference 로 사용
  ref_idx <- which(grepl(ref_tag, names(msa_seqs), fixed = TRUE))
  if (length(ref_idx) == 0) {
    warning(sprintf("  [%s] reference '%s' 를 MSA에서 찾지 못함", gene, ref_tag))
    ref_idx <- 1
  }
  ref_key <- names(msa_seqs)[ref_idx[1]]
  ref_seq <- msa_seqs[[ref_key]]

  # 샘플 목록 (reference 제외)
  sample_keys <- names(msa_seqs)[names(msa_seqs) != ref_key]
  if (length(sample_keys) == 0) {
    warning(sprintf("  [%s] 샘플 서열이 없음", gene))
    return(data.frame(gene = character(), pos = integer(), ref = character()))
  }

  # 정렬 길이 (모든 서열이 같아야 함 - MAFFT 출력 보장)
  aln_len <- nchar(ref_seq)

  # 서열을 문자 벡터로 분해 (strsplit 후 첫 번째 원소)
  ref_chars  <- strsplit(ref_seq, "")[[1]]
  samp_chars <- lapply(msa_seqs[sample_keys], function(s) {
    chars <- strsplit(s, "")[[1]]
    # 길이 불일치 보정 (이상 케이스 대비)
    if (length(chars) < aln_len) chars <- c(chars, rep("-", aln_len - length(chars)))
    chars[seq_len(aln_len)]
  })

  # 컬럼별 순회하며 variant 위치 수집
  variant_rows <- list()
  ref_pos      <- 0L   # reference 기준 1-based 위치 (gap 열은 카운트 제외)

  for (col in seq_len(aln_len)) {
    r <- ref_chars[col]

    if (r == "-") next   # reference gap = insertion → 건너뜀
    ref_pos <- ref_pos + 1L

    # 각 샘플의 이 위치 문자
    calls <- sapply(samp_chars, function(ch) ch[col])

    # 불확정 문자 제외 후, reference와 다른 것이 하나라도 있으면 variant
    is_variable <- any(
      calls != r & !calls %in% c("-", "N", "X", "?", "n", "x")
    )

    if (is_variable) {
      row <- c(list(gene = gene, pos = ref_pos, ref = r), as.list(calls))
      variant_rows[[length(variant_rows) + 1]] <- row
    }
  }

  if (length(variant_rows) == 0) {
    message(sprintf("  [%s] variant 없음 (모든 샘플이 reference와 동일)", gene))
    empty <- as.data.frame(
      matrix(ncol = 3 + length(sample_keys), nrow = 0)
    )
    colnames(empty) <- c("gene", "pos", "ref", sample_keys)
    return(empty)
  }

  # list-of-lists → data.frame
  df <- do.call(rbind, lapply(variant_rows, as.data.frame, stringsAsFactors = FALSE))
  df$pos <- as.integer(df$pos)

  message(sprintf("  [%s] variant %d 개 위치 발견", gene, nrow(df)))
  return(df)
}


# ── 5. 샘플 FASTA 파일 수집 유틸 ────────────────────────────────────────────
# 지정된 디렉토리에서 "{barcode}.*.{GENE}.{type}.fasta" 패턴 파일 목록 반환
# 파일명 패턴: barcode01.Pf3D7.CRT.nt.fasta  (GENE 은 3번째 점 구분 필드)
collect_sample_fastas <- function(dir_path, gene, type) {
  # type = "nt" or "aa"
  if (!dir.exists(dir_path)) stop(paste("디렉토리 없음:", dir_path))

  pattern <- sprintf("\\.%s\\.%s\\.fasta$", gene, type)
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)

  # placeholder(EMPTY) 파일 제거
  files <- files[!grepl("EMPTY", basename(files))]

  if (length(files) == 0) {
    warning(sprintf("  [%s] %s 샘플 FASTA 파일 없음 (패턴: %s)", gene, type, pattern))
  } else {
    message(sprintf("  [%s] %s 파일 %d 개 발견", gene, type, length(files)))
  }
  return(files)
}


# ── 6. 합본 FASTA 생성 및 헤더 정규화 ──────────────────────────────────────
# reference FASTA 한 개 + 샘플 FASTA 여러 개를 합쳐서 combined FASTA 저장
# 헤더는 "barcode01|GENE" 형태로 통일 (공백·특수문자 방지)
build_combined_fasta <- function(ref_path, sample_files, gene, out_path) {
  combined <- c()

  # 3D7 reference 먼저 추가
  ref_seqs <- read_fasta(ref_path)
  # 헤더를 "3D7|GENE" 으로 통일
  names(ref_seqs) <- paste0("3D7|", gene)
  combined <- c(combined, ref_seqs)

  # 샘플별 서열 추가
  for (f in sample_files) {
    seqs <- read_fasta(f)
    # 헤더에서 barcode 이름 추출: "barcode01.Pf3D7|CRT" → "barcode01"
    # 파일명에서 barcode 추출이 더 안정적
    barcode <- sub("^(barcode[0-9]+)\\..*$", "\\1", basename(f))

    # 동일한 헤더가 있으면 덮어쓰기 방지 (중복 샘플 파일 경우 경고)
    new_header <- paste0(barcode, "|", gene)
    if (new_header %in% names(combined)) {
      warning(sprintf("  [%s] 중복 샘플 헤더 무시: %s", gene, new_header))
      next
    }
    names(seqs) <- new_header
    combined <- c(combined, seqs[1])  # 첫 번째 레코드만 사용
  }

  write_fasta(combined, out_path)
  message(sprintf("  [%s] combined FASTA 저장: %s (%d 개 서열)", gene, basename(out_path), length(combined)))
  return(combined)
}


# ── 7. 실제 처리 함수 (NT 또는 AA) ──────────────────────────────────────────
# type = "nt" | "aa"
process_type <- function(type) {
  message(sprintf("\n[%s] === %s 처리 시작 ===", GENE, toupper(type)))

  # 7-1. 파일 경로 설정
  dir_path    <- if (type == "nt") opt$nt_dir else opt$aa_dir
  ref_path    <- if (type == "nt") opt$ref_nt  else opt$ref_aa
  combined_fa <- file.path(OUTDIR, sprintf("%s.%s.combined.fasta", GENE, type))
  msa_fa      <- file.path(OUTDIR, sprintf("%s.%s.msa.fasta",      GENE, type))
  variant_tsv <- file.path(OUTDIR, sprintf("%s.%s.variants.tsv",   GENE, type))

  # 7-2. 샘플 파일 수집
  sample_files <- collect_sample_fastas(dir_path, GENE, type)
  if (length(sample_files) == 0) {
    message(sprintf("  [%s] 샘플 없음 → 건너뜀", GENE))
    return(invisible(NULL))
  }

  # 7-3. reference 확인
  if (!file.exists(ref_path)) {
    stop(sprintf("3D7 reference 파일 없음: %s", ref_path))
  }

  # 7-4. 합본 FASTA 생성 (3D7 + 모든 샘플)
  combined_seqs <- build_combined_fasta(ref_path, sample_files, GENE, combined_fa)

  # 7-5. MAFFT 다중서열정렬
  run_mafft(combined_fa, msa_fa)

  # 7-6. 정렬 결과 읽기
  msa_seqs <- read_fasta(msa_fa)
  if (length(msa_seqs) == 0) stop(sprintf("MSA 결과가 비어있음: %s", msa_fa))

  # 7-7. 변이표 생성 (3D7 기준)
  vt <- make_variant_table(msa_seqs, ref_tag = "3D7", gene = GENE)

  # 7-8. 변이표 저장
  write.table(vt, variant_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("  [%s] %s 변이표 저장: %s", GENE, toupper(type), basename(variant_tsv)))
}


# ── 8. NT / AA 순서로 처리 실행 ─────────────────────────────────────────────
process_type("nt")
process_type("aa")

message(sprintf("\n[%s] 완료 → %s/", GENE, OUTDIR))
