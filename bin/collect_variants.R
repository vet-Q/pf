#!/usr/bin/env Rscript
# collect_variants.R
#
# 파이프라인이 생성한 per-sample, per-gene VCF 파일들을 하나의
# variant_AF_all.tsv 로 병합합니다.
#
# Usage:
#   Rscript collect_variants.R \
#     --vcf_list vcf_list.txt \
#     --out_tsv  variant_AF_all.tsv \
#     --genes    "CRT,DHFR,DHPS,K13,MDR1,MSP2,PMI,PMIII"

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(optparse)
})

# ── CLI ─────────────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("--vcf_list"), type = "character",
              help = "Text file with one .vcf.gz path per line"),
  make_option(c("--out_tsv"),  type = "character",
              default = "variant_AF_all.tsv",
              help = "Output path for merged TSV [default: %default]"),
  make_option(c("--genes"),    type = "character",
              default = "CRT,DHFR,DHPS,K13,MDR1,MSP2,PMI,PMIII",
              help = "Comma-separated gene names (used for gene extraction from filename)")
)

opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$vcf_list))

GENES <- str_split(opt$genes, ",", simplify = TRUE)[1, ]
GENES <- GENES[GENES != ""]

# ── 헬퍼 ─────────────────────────────────────────────────────────────────────
extract_barcode <- function(path) {
  bc <- str_match(basename(path), "(barcode\\d+)")[, 2]
  if (is.na(bc)) bc <- "UNKNOWN"
  bc
}

extract_gene_from_path <- function(path, genes = GENES) {
  # 파일명에서 .GENE. 패턴 우선 (스테이징 후 barcode01.CRT.bcftools.vcf.gz 형태)
  m <- vapply(genes, function(g) {
    grepl(paste0("\\.", g, "\\."), basename(path))
  }, logical(1))
  if (any(m)) return(genes[which(m)[1]])

  # 폴더 구조에서 추출 (원본 경로 사용 시 폴백)
  m2 <- vapply(genes, function(g) {
    grepl(paste0("(^|/)", g, "(/|$)"), path)
  }, logical(1))
  if (any(m2)) return(genes[which(m2)[1]])

  "GENE"
}

has_format_field <- function(vcf, field) {
  gt_hdr <- colnames(vcf@gt)
  if (!("FORMAT" %in% gt_hdr)) return(FALSE)
  fmt <- vcf@gt[, "FORMAT"]
  any(grepl(paste0("(^|:)", field, "(:|$)"), fmt))
}

read_one_vcf <- function(vcf_path) {
  vcf <- tryCatch(
    read.vcfR(vcf_path, verbose = FALSE),
    error = function(e) {
      message("[ERR] cannot read: ", vcf_path, "\n", e$message)
      return(NULL)
    }
  )
  if (is.null(vcf) || nrow(vcf@fix) == 0) return(NULL)

  if (!has_format_field(vcf, "AD") || !has_format_field(vcf, "DP")) {
    message("[SKIP] no AD/DP in FORMAT: ", vcf_path)
    return(NULL)
  }

  fix <- as.data.frame(vcf@fix) %>%
    transmute(
      CHROM = CHROM,
      POS   = as.integer(POS),
      REF   = REF,
      ALT   = ALT
    )

  ad <- extract.gt(vcf, element = "AD", as.numeric = FALSE)
  dp <- extract.gt(vcf, element = "DP", as.numeric = FALSE)
  gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

  sample_name <- colnames(ad)[1]

  out <- fix %>%
    mutate(
      barcode = extract_barcode(vcf_path),
      gene    = extract_gene_from_path(vcf_path),
      sample  = sample_name,
      GT      = as.character(gt[, 1]),
      AD_raw  = as.character(ad[, 1]),
      DP      = suppressWarnings(as.integer(dp[, 1]))
    ) %>%
    mutate(
      AD_raw = ifelse(is.na(AD_raw), ".", AD_raw),

      AD_REF = ifelse(
        AD_raw == "." | AD_raw == "",
        NA_integer_,
        suppressWarnings(as.integer(str_split_fixed(AD_raw, ",", 2)[, 1]))
      ),

      AD_ALT = ifelse(
        AD_raw == "." | AD_raw == "",
        NA_integer_,
        {
          parts <- str_split(AD_raw, ",", simplify = FALSE)
          vapply(parts, function(p) {
            nums <- suppressWarnings(as.integer(p))
            if (length(nums) < 2) return(NA_integer_)
            sum(nums[-1], na.rm = TRUE)   # multi-ALT는 합산
          }, integer(1))
        }
      )
    ) %>%
    mutate(
      AD_DP = AD_REF + AD_ALT,
      AF    = ifelse(!is.na(AD_DP) & AD_DP > 0, AD_ALT / AD_DP, NA_real_)
    ) %>%
    filter(!is.na(AF)) %>%
    select(barcode, gene, sample, CHROM, POS, REF, ALT, GT, DP, AD_REF, AD_ALT, AD_DP, AF)

  if (nrow(out) == 0) return(NULL)
  out
}

# ── 메인 ─────────────────────────────────────────────────────────────────────
vcf_files <- readLines(opt$vcf_list, warn = FALSE)
vcf_files <- vcf_files[nzchar(vcf_files)]

if (length(vcf_files) == 0) stop("vcf_list is empty.")

message("Processing ", length(vcf_files), " VCF files...")

df <- bind_rows(map(vcf_files, read_one_vcf))

if (is.null(df) || nrow(df) == 0) {
  message("[WARN] No variants extracted — writing empty TSV.")
  df <- data.frame(
    barcode = character(), gene = character(), sample = character(),
    CHROM = character(), POS = integer(), REF = character(), ALT = character(),
    GT = character(), DP = integer(), AD_REF = integer(), AD_ALT = integer(),
    AD_DP = integer(), AF = numeric()
  )
}

write.table(df, opt$out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
message("Written: ", opt$out_tsv, " (", nrow(df), " rows)")
