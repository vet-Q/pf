#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

gene_bins_root <- args[1]  # e.g. /.../results/gene_bins
genes_arg      <- args[2]  # e.g. "CRT,DHFR,DHPS,K13,MDR1,MSP2"
out_pdf        <- args[3]

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(readr)
})

genes <- str_split(genes_arg, ",", simplify = TRUE) |> as.character()
genes <- genes[genes != ""]

# gene_bins_root 아래 barcode 디렉토리들
barcodes <- list.dirs(gene_bins_root, recursive = FALSE, full.names = TRUE)

read_depth_one <- function(depth_file, barcode, gene){
  # 파일이 비면 skip
  if (!file.exists(depth_file) || file.info(depth_file)$size == 0) return(NULL)
  df <- read.table(depth_file, col.names = c("chr","pos","depth"))
  df$barcode <- barcode
  df$gene <- gene
  df
}

all <- list()

for (bc_dir in barcodes) {
  barcode <- basename(bc_dir)
  
  for (g in genes) {
    gp <- file.path(bc_dir, g)
    if (!dir.exists(gp)) next
    
    depth_files <- list.files(gp, pattern="^bin_.*_depth\\.txt$", full.names=TRUE)
    if (length(depth_files) == 0) next
    
    # bin들 전부 합쳐 “포지션 depth 분포”로 사용(원하면 bin별 facet도 가능)
    dfs <- lapply(depth_files, read_depth_one, barcode=barcode, gene=g)
    dfs <- dfs[!sapply(dfs, is.null)]
    if (length(dfs) == 0) next
    
    all[[length(all)+1]] <- bind_rows(dfs)
  }
}

if (length(all) == 0) {
  message("[INFO] No depth files found.")
  quit(save="no", status=0)
}

dat <- bind_rows(all) %>%
  mutate(
    barcode = factor(barcode, levels = sort(unique(barcode))),
    gene = factor(gene, levels = genes)
  )

# barcode×gene 별 평균 계산(점으로 표시)
means <- dat %>%
  group_by(gene, barcode) %>%
  summarise(mean_depth = mean(depth, na.rm=TRUE), .groups="drop")

p <- ggplot(dat, aes(x=barcode, y=depth)) +
  geom_boxplot(outlier.size=0.2) +
  geom_point(data=means, aes(y=mean_depth), size=1.2) +
  facet_wrap(~gene, nrow=1, scales="free_y") +
  theme_bw(base_size=9) +
  labs(x=NULL, y="Depth (per-position distribution)", title="Depth distribution by barcode (mean shown as dot)") +
  theme(
    axis.text.x = element_text(angle=60, hjust=1, vjust=1),
    panel.grid.minor = element_blank()
  )

ggsave(out_pdf, p, width = 2.2*length(genes), height=4.0, limitsize=FALSE)
message("[INFO] saved: ", out_pdf)
