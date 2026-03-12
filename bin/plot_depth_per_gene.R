#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
depth_file <- args[1]   # ex) barcode01/all_amplicons.depth.txt
bed_file   <- args[2]   # ex) amplicons.bed
out_dir    <- args[3]   # ex) barcode01

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

dir.create(file.path(out_dir, "depth_by_gene"), showWarnings = FALSE)
dir.create(file.path(out_dir, "plots_by_gene"), showWarnings = FALSE)

## ========= 1) gene별 depth 텍스트 + 기본 플롯 =========

# 1) depth 파일 읽기 (samtools depth 결과)
# chr pos depth
depth <- read.table(depth_file,
                    header = FALSE,
                    col.names = c("chr", "pos", "depth"))

# 2) BED 읽기 (amplicon 정보)
# chr start end gene
bed <- read.table(bed_file,
                  header = FALSE,
                  col.names = c("chr", "start", "end", "gene"))

genes <- unique(bed$gene)

for (g in genes) {

  region <- bed %>% filter(gene == g)

  if (nrow(region) != 1) {
    message("[WARN] ", g, " 에 대해 BED가 1행이 아님 (무시)")
    next
  }

  sub <- depth %>%
    filter(chr == region$chr[1],
           pos >= region$start[1],
           pos <= region$end[1])

  if (nrow(sub) == 0) {
    message("[WARN] ", g, " 에서 depth 행이 없음")
    next
  }

  ## --- 텍스트 파일 저장 (예: MSP2_depth.txt)
  out_txt <- file.path(out_dir, "depth_by_gene",
                       paste0(g, "_depth.txt"))
  write.table(sub,
              file      = out_txt,
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE,
              sep       = "\t")

  ## --- 스무스 depth 플롯 (loess)
  p <- ggplot(sub, aes(x = pos, y = depth)) +
    geom_point(alpha = 0.3, size = 0.4) +
    geom_smooth(method = "loess", se = FALSE, span = 0.1) +
    theme_bw(base_size = 11) +
    labs(title = paste0("Depth across ", g),
         x = paste0(region$chr[1], ":", region$start[1], "-", region$end[1]),
         y = "Depth") +
    theme(plot.title = element_text(hjust = 0.5))

  out_pdf <- file.path(out_dir, "plots_by_gene",
                       paste0(g, "_depth.pdf"))
  ggsave(out_pdf, plot = p, width = 7, height = 4)
}

## ========= 2) MSP2 멀티패널 (옵션) =========
## out_dir 안에
##   - msp2_0-1k_depth.txt ... msp2_4kplus_depth.txt
##   - MSP2.fa
## 가 있을 때만 수행

msp2_depth_files <- list.files(
  path       = out_dir,
  pattern    = "^msp2_.*_depth.txt$",
  full.names = TRUE
)

msp2_fa_path <- file.path(out_dir, "MSP2.fa")

if (length(msp2_depth_files) > 0 && file.exists(msp2_fa_path)) {

  message("[INFO] MSP2 멀티패널 플롯을 생성합니다...")

  ## --- 2-1) read-length bin stacked area ---

  depth_bins_list <- lapply(msp2_depth_files, function(f) {

    # 파일 사이즈 0이면 스킵
    if (file.info(f)$size == 0) {
      message("Skipping empty file: ", basename(f))
      return(NULL)
    }

    df <- read.table(f, col.names = c("chr", "pos", "depth"))

    # bin 이름: msp2_0-1k_depth.txt -> 0-1k
    bin <- gsub("msp2_|_depth.txt", "", basename(f))
    df$bin <- bin
    df
  })

  # NULL 제거
  depth_bins_list <- depth_bins_list[!sapply(depth_bins_list, is.null)]

  if (length(depth_bins_list) > 0) {

    depth_bins <- dplyr::bind_rows(depth_bins_list)

    bin_levels_all <- c("0-1k", "1-2k", "2-3k", "3-4k", "4kplus")
    bin_levels_use <- bin_levels_all[bin_levels_all %in% unique(depth_bins$bin)]
    depth_bins$bin <- factor(depth_bins$bin, levels = bin_levels_use)

    p_cov <- ggplot(depth_bins, aes(x = pos, y = depth, fill = bin)) +
      geom_area(position = "stack") +
      theme_bw(base_size = 11) +
      labs(title = "MSP2 coverage by read length",
           x = NULL,
           y = "Coverage") +
      theme(
        legend.position = "top",
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
      )

    ## --- 2-2) A+T% (슬라이딩 윈도우) ---

    fa_lines <- readLines(msp2_fa_path)
    seq_str  <- toupper(paste(fa_lines[-1], collapse = ""))
    bases    <- strsplit(seq_str, "")[[1]]
    n        <- length(bases)

    window <- 10  # 10bp 창

    at_vals <- sapply(1:(n - window + 1), function(i) {
      w <- bases[i:(i + window - 1)]
      at <- sum(w %in% c("A", "T"))
      at / window * 100
    })

    msp2_start <- 271576  # MSP2 시작 좌표(Pf3D7_02_v3:271576-274917)
    pos_at <- msp2_start + (window / 2) + seq_along(at_vals) - 1

    df_at <- data.frame(pos = pos_at, at = at_vals)

    p_at <- ggplot(df_at, aes(x = pos, y = at)) +
      geom_line(color = "steelblue", linewidth = 0.5) +
      theme_bw(base_size = 11) +
      labs(title = "A+T% across MSP2",
           x = NULL,
           y = "A+T (%)") +
      theme(
        plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
      )

    ## --- 2-3) Homopolymer length per position ---

    hp_len <- integer(n)
    i <- 1
    while (i <= n) {
      base <- bases[i]
      j <- i
      while (j <= n && bases[j] == base) {
        j <- j + 1
      }
      run_len <- j - i
      hp_len[i:(j - 1)] <- run_len
      i <- j
    }

    pos_hp <- msp2_start + 0:(n - 1)
    df_hp <- data.frame(pos = pos_hp, hp = hp_len)

    p_hp <- ggplot(df_hp, aes(x = pos, y = hp)) +
      geom_line(color = "firebrick", linewidth = 0.5) +
      theme_bw(base_size = 11) +
      labs(title = "Homopolymer length across MSP2",
           x = "Genomic position (bp)",
           y = "Homopolymer length (bp)") +
      theme(
        plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
      )

    ## --- 2-4) 세 패널 결합 & 저장 ---

    p_multi <- p_cov / p_at / p_hp + plot_layout(heights = c(3, 2, 2))

    out_multi <- file.path(out_dir,
                           "MSP2_multipanel_coverage_AT_homopolymer.pdf")

    ggsave(out_multi, p_multi, width = 8, height = 7)
    message("[INFO] MSP2 멀티패널 저장: ", out_multi)
  } else {
    message("[WARN] MSP2 depth bin 파일이 모두 비어있음. 멀티패널 생략.")
  }
} else {
  message("[INFO] MSP2 멀티패널용 파일이 없어 생략합니다.")
}

