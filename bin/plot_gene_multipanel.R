#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
barcode_dir <- args[1]   # e.g. barcode01
out_pdf     <- args[2]   # e.g. barcode01_gene_multipanel.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

genes <- list.dirs(barcode_dir, recursive = FALSE, full.names = FALSE)

bin_levels <- c("0-1k","1-2k","2-3k","3-4k","4kplus")

read_depth_bins <- function(gene_path){
  files <- list.files(gene_path, pattern="^bin_.*_depth.txt$", full.names=TRUE)
  if(length(files) == 0) return(NULL)

  dfs <- lapply(files, function(f){
    if(file.info(f)$size == 0) return(NULL)
    df <- read.table(f, col.names=c("chr","pos","depth"))
    bin <- basename(f) |>
      str_replace("^bin_", "") |>
      str_replace("_depth.txt$", "")
    df$bin <- bin
    df
  })
  dfs <- dfs[!sapply(dfs, is.null)]
  if(length(dfs) == 0) return(NULL)
  dplyr::bind_rows(dfs)
}

read_readlen <- function(gene_path){
  f <- file.path(gene_path, "readlen.txt")
  if(!file.exists(f) || file.info(f)$size == 0) return(NULL)
  x <- scan(f, quiet=TRUE)
  data.frame(read_len = x)
}

plots <- list()

for(g in genes){
  gp <- file.path(barcode_dir, g)

  depth_bins <- read_depth_bins(gp)
  if(is.null(depth_bins)) next

  depth_bins$bin <- factor(depth_bins$bin, levels=bin_levels)

  p_cov <- ggplot(depth_bins, aes(x=pos, y=depth, fill=bin)) +
    geom_area(position="stack") +
    theme_bw(base_size = 10) +
    labs(title = paste0(g, " coverage by read length bin"),
         x=NULL, y="Coverage") +
    theme(legend.position="top")

  rl <- read_readlen(gp)
  if(!is.null(rl)){
    rl$gene <- g
    p_len <- ggplot(rl, aes(x=gene, y=read_len/1000)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.15, alpha=0.25, size=0.6) +
      theme_bw(base_size = 10) +
      labs(title = paste0(g, " read length distribution"),
           x=NULL, y="Read length (kbp)")

    p <- p_cov / p_len + plot_layout(heights=c(3,2))
  } else {
    p <- p_cov
  }

  plots[[g]] <- p
}

if(length(plots) == 0){
  message("[INFO] No plots generated (no depth bin files).")
  quit(save="no", status=0)
}

pdf(out_pdf, width=8, height=7)
for(nm in names(plots)){
  print(plots[[nm]])
}
dev.off()

message("[INFO] saved: ", out_pdf)

