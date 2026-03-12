#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
barcode_dir <- args[1]   # e.g. .../gene_bins/barcode01
genes_arg   <- args[2]   # "CRT,DHFR,DHPS,K13,MDR1,MSP2"
out_pdf     <- args[3]

# 패키지를 시작할때 뜨는 내용을 안뜨도록 조용히 부르는 메서드
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

# GENE을 벡터로 만들고, 딕셔너리로 만드는거?
genes <- str_split(genes_arg, ",", simplify = TRUE) |> as.character()
genes <- genes[genes != ""]

#바코드 이름 추출
barcode <- basename(barcode_dir)

# depth bin의 순서를 정하기 위한 기준벡터. 나중에 factor로 쓰기 위함
bin_levels <- c("0-1k","1-2k","2-3k","3-4k","4kplus")

# “gene 하나의 디렉토리”를 입력으로 받아서 bin별 depth 데이터를 하나의 data.frame으로 합쳐주는 함수
read_depth_bins <- function(gene_path){
  files <- list.files(gene_path, pattern="^bin_.*_depth\\.txt$", full.names=TRUE)
  if(length(files) == 0) return(NULL)
  
  dfs <- lapply(files, function(f){
    if(file.info(f)$size == 0) return(NULL)
    df <- read.table(f, col.names=c("chr","pos","depth"))
    bin <- basename(f) |>
      str_replace("^bin_", "") |>
      str_replace("_depth\\.txt$", "")
    df$bin <- bin
    df
  })
  dfs <- dfs[!sapply(dfs, is.null)]
  if(length(dfs) == 0) return(NULL)
  bind_rows(dfs)
}

# ------------------------------------------------------------------------------
# read_readlen(): gene 폴더(gene_path) 안의 readlen.txt를 읽어서
# read length 값들을 data.frame 형태로 반환하는 함수
# - readlen.txt가 없거나(파일 없음) 비어있으면(NULL) 반환
# - scan()은 텍스트 파일에서 숫자 벡터를 빠르게 읽는 함수(한 줄/여러 줄 상관없이 숫자만 쭉)
# ------------------------------------------------------------------------------
read_readlen <- function(gene_path){
  f <- file.path(gene_path, "readlen.txt")
  if(!file.exists(f) || file.info(f)$size == 0) return(NULL)
  x <- scan(f, quiet=TRUE)
  data.frame(read_len = x)
}


# ------------------------------------------------------------------------------
# make_gene_plot(g): gene 이름(g)을 받아서
# 1) depth bin별 coverage plot(p_cov)을 만들고
# 2) readlen.txt가 있으면 read length boxplot(p_len)도 만들어
# 3) patchwork로 위/아래로 합쳐서 반환
# - 중간에 필요한 파일/폴더가 없으면(NULL) 반환해서 자동 스킵
# ------------------------------------------------------------------------------
make_gene_plot <- function(g){
  gp <- file.path(barcode_dir, g) # barcode_dir/gene 경로
  if(!dir.exists(gp)) return(NULL)  # 해당 gene 디렉토리가 없으면 스킵(NULL)
  
  depth_bins <- read_depth_bins(gp)   # bin_XXX_depth.txt 들 읽어서 한 df로 합치기(앞에서 정의한 함수)
  if(is.null(depth_bins)) return(NULL) # depth 파일이 하나도 없거나 다 비었으면 스킵(NULL)
  
  # bin을 factor로 바꿔서 "0-1k -> 1-2k -> ..." 순서가 고정되게 함(그래프/범례 순서 안정화)
  depth_bins$bin <- factor(depth_bins$bin, levels=bin_levels)
  
  
  # --------------------------------------------------------------------------
  # p_cov: coverage(깊이) 영역 그래프
  # - x축: position(pos)
  # - y축: depth
  # - fill: bin (bin별로 색이 달라지고, position="stack"으로 누적 영역이 됨)
  # --------------------------------------------------------------------------
  p_cov <- ggplot(depth_bins, aes(x=pos, y=depth, fill=bin)) +
    geom_area(position="stack") +
    theme_bw(base_size = 7) +  # 흰 배경 테마(글자 기본 7)
    labs(title=g, x=NULL, y="Cov") +  # 제목=gene, x축 라벨 제거, y축=Cov
    theme(
      legend.position="none",  # gene별 플롯에서는 범례 숨김(많아지면 지저분해짐)
      plot.title = element_text(size=8, face="bold", hjust=0.5),  # 제목 폰트/가운데 정렬
      axis.text.x = element_blank(), # x축 글자 숨김(위쪽 플롯은 보통 공유하므로)
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # --------------------------------------------------------------------------
  # read length 데이터가 있으면 boxplot 추가
  # readlen.txt는 보통 read 길이(bp)들이라서 /1000 해서 kb 단위로 표시
  # --------------------------------------------------------------------------
  
  rl <- read_readlen(gp)  # readlen 읽기(없으면 NULL)
  if(!is.null(rl)){
    p_len <- ggplot(rl, aes(x=1, y=read_len/1000)) + 
      geom_jitter(
        width = 0.15,              # x축으로 살짝 흩뿌리기
        size = 0.8,                # 점 크기
        alpha = 0.4,               # 투명도 (연하게)
        color = "grey70"           # 연한 회색
      ) +
      geom_boxplot(outlier.shape=NA) +   # 극단치 점(outlier) 숨김(스케일 망가지는 것 방지)
      theme_bw(base_size = 7) +
      labs(x=NULL, y="Len(kb)") +  # y축 라벨: Len(kb)
      theme(
        axis.text.x = element_blank(),  # x축 글자 숨김(의미 없는 더미축이므로)
        axis.ticks.x = element_blank(),    # x축 눈금 숨김
        panel.grid.minor = element_blank()
      )
    # patchwork 문법:
    # p_cov / p_len  => 위(p_cov), 아래(p_len)로 배치
    # heights=c(4,2) => 위:아래 높이 비율 4:2
    p_cov / p_len + plot_layout(heights=c(4,2))
  } else {
    # readlen.txt가 없으면 coverage plot만 반환
    p_cov
  }
}

# ------------------------------------------------------------------------------
# genes 벡터의 각 gene에 대해 make_gene_plot을 적용
# - 결과는 list로 모이고, NULL(스킵된 gene)은 제거
# ------------------------------------------------------------------------------
plots <- lapply(genes, make_gene_plot)
plots <- plots[!sapply(plots, is.null)]


# ------------------------------------------------------------------------------
# 플롯이 하나도 없으면 메시지 출력 후 정상 종료(status=0)
# (실행 환경이 Rscript라면 quit()로 종료하는 패턴)
# ------------------------------------------------------------------------------
if(length(plots) == 0){
  message("[INFO] No plots generated.")
  quit(save="no", status=0)
}

# ------------------------------------------------------------------------------
# wrap_plots: 플롯 리스트를 한 화면에 배치
# nrow=1 => 한 줄에 가로로 쭉 나열
# ------------------------------------------------------------------------------
p_all <- wrap_plots(plots, nrow=1)

# 가로로 길어지니 width를 충분히
pdf(out_pdf, width=2.2*length(plots), height=4.0)
print(p_all)  # patchwork 결과 출력(장치에 그리기)
dev.off()  # PDF 장치 닫기(파일 저장 완료)

message("[INFO] saved: ", out_pdf)
