#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(optparse)
})

# -------------------------
# CLI
# -------------------------
option_list <- list(
  make_option(c("--vcf_list"), type="character", help="Text file with one .vcf.gz path per line"),
  make_option(c("--out_fig1"), type="character", default="Figure1_by_gene.pdf"),
  make_option(c("--out_fig2"), type="character", default="Figure2_overall.pdf"),
  make_option(c("--min_dp"),  type="integer",   default=15),
  make_option(c("--af_cut"),  type="double",    default=0.50),
  make_option(c("--shade_xmin"), type="double", default=0.48),
  make_option(c("--shade_xmax"), type="double", default=0.52),
  make_option(c("--genes"), type="character", default="CRT,DHFR,DHPS,K13,MDR1,MSP2,PMI,PMIII",
              help="Comma-separated gene folder names (order controls facet order)")
)

opt <- parse_args(OptionParser(option_list=option_list))
stopifnot(!is.null(opt$vcf_list))

GENES <- str_split(opt$genes, ",", simplify = TRUE)[1,]
GENES <- GENES[GENES != ""]

# -------------------------
# Helpers (same logic, but more Nextflow-friendly)
# -------------------------
extract_barcode <- function(path){
  bc <- str_match(path, "(barcode\\d+)\\.Pf3D7\\.final\\.sorted")[,2]
  if (is.na(bc)) bc <- str_match(path, "(barcode\\d+)")[,2]
  if (is.na(bc)) bc <- "UNKNOWN"
  bc
}

extract_gene_from_path <- function(path, genes = GENES){
  m <- vapply(genes, function(g){
    grepl(paste0("(^|/)", g, "(/|$)"), path)
  }, logical(1))
  if (any(m)) return(genes[which(m)[1]])
  "GENE"
}

has_format_field <- function(vcf, field){
  gt_hdr <- colnames(vcf@gt)
  if (!("FORMAT" %in% gt_hdr)) return(FALSE)
  fmt <- vcf@gt[,"FORMAT"]
  any(grepl(paste0("(^|:)", field, "(:|$)"), fmt))
}

read_one_vcf <- function(vcf_path, MIN_DP){
  vcf <- tryCatch(read.vcfR(vcf_path, verbose = FALSE),
                  error = function(e) return(NULL))
  if (is.null(vcf)) return(NULL)

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

  # single sample 가정
  sample <- colnames(ad)[1]

  out <- fix %>%
    mutate(
      file    = basename(vcf_path),
      barcode = extract_barcode(vcf_path),
      gene    = extract_gene_from_path(vcf_path),
      sample  = sample,
      AD_raw  = as.character(ad[,1]),
      DP      = suppressWarnings(as.integer(dp[,1]))
    ) %>%
    mutate(
      AD_raw = ifelse(is.na(AD_raw), ".", AD_raw),

      refAD = ifelse(
        AD_raw == "." | AD_raw == "",
        NA_integer_,
        suppressWarnings(as.integer(str_split_fixed(AD_raw, ",", 2)[,1]))
      ),

      altAD = ifelse(
        AD_raw == "." | AD_raw == "",
        NA_integer_,
        {
          parts <- str_split(AD_raw, ",", simplify = FALSE)
          vapply(parts, function(p){
            nums <- suppressWarnings(as.integer(p))
            if (length(nums) < 2) return(NA_integer_)
            sum(nums[-1], na.rm = TRUE)  # multi-ALT면 ALT들 합산
          }, integer(1))
        }
      )
    ) %>%
    mutate(
      totalAD = refAD + altAD,
      AF = ifelse(!is.na(totalAD) & totalAD > 0, altAD / totalAD, NA_real_)
    ) %>%
    filter(!is.na(AF), !is.na(DP), DP >= MIN_DP)

  if (nrow(out) == 0) return(NULL)
  out
}

theme_fig <- function(base_size = 12){
  theme_bw(base_size = base_size) +
    theme(
      strip.text = element_text(color = "#7A1F2B", face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.grid.minor = element_blank()
    )
}

# -------------------------
# Main
# -------------------------
vcf_files <- readLines(opt$vcf_list, warn = FALSE)
vcf_files <- vcf_files[nzchar(vcf_files)]
if (length(vcf_files) == 0) stop("vcf_list is empty.")

df <- bind_rows(map(vcf_files, ~read_one_vcf(.x, MIN_DP = opt$min_dp)))
if (nrow(df) == 0) stop("No valid variants after filtering (AD/DP missing or MIN_DP too strict).")

df$gene <- factor(df$gene, levels = GENES)

# Figure 1 (facet)
p1_scatter <- ggplot(df, aes(x = AF, y = log10(DP))) +
  annotate("rect", xmin = opt$shade_xmin, xmax = opt$shade_xmax, ymin = -Inf, ymax = Inf, alpha = 0.15) +
  geom_vline(xintercept = opt$af_cut, linetype = "dashed") +
  geom_point(alpha = 0.20, size = 1) +
  facet_wrap(~gene, ncol = 3) +
  labs(x = "Allele fraction (ALT / (REF+ALT))",
       y = expression(log[10]*"(Depth)")) +
  theme_fig()

p1_hist <- ggplot(df, aes(x = AF)) +
  annotate("rect", xmin = opt$shade_xmin, xmax = opt$shade_xmax, ymin = -Inf, ymax = Inf, alpha = 0.15) +
  geom_vline(xintercept = opt$af_cut, linetype = "dashed") +
  geom_histogram(bins = 40) +
  labs(x = NULL, y = "# of variants") +
  theme_fig() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

FIG1 <- p1_hist / p1_scatter + plot_layout(heights = c(1, 4))
ggsave(opt$out_fig1, FIG1, width = 12, height = 8, dpi = 300, limitsize = FALSE)

# Figure 2 (overall)
p2_scatter <- ggplot(df, aes(x = AF, y = log10(DP))) +
  annotate("rect", xmin = opt$shade_xmin, xmax = opt$shade_xmax, ymin = -Inf, ymax = Inf, alpha = 0.15) +
  geom_vline(xintercept = opt$af_cut, linetype = "dashed") +
  geom_point(alpha = 0.20, size = 1) +
  labs(x = "Allele fraction (ALT / (REF+ALT))",
       y = expression(log[10]*"(Depth)")) +
  theme_fig()

p2_hist <- ggplot(df, aes(x = AF)) +
  annotate("rect", xmin = opt$shade_xmin, xmax = opt$shade_xmax, ymin = -Inf, ymax = Inf, alpha = 0.15) +
  geom_vline(xintercept = opt$af_cut, linetype = "dashed") +
  geom_histogram(bins = 40) +
  labs(x = NULL, y = "# of variants") +
  theme_fig() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

FIG2 <- p2_hist / p2_scatter + plot_layout(heights = c(1, 4))
ggsave(opt$out_fig2, FIG2, width = 10, height = 7, dpi = 300, limitsize = FALSE)

message("[OK] wrote: ", opt$out_fig1)
message("[OK] wrote: ", opt$out_fig2)
