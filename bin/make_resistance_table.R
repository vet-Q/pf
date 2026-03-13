#!/usr/bin/env Rscript
# =============================================================================
# make_resistance_table.R
#
# 목적:
#   variant_AF_all.tsv (bcftools 변이표) + cds_coords.tsv를 이용해
#   P. falciparum 주요 약제내성 마커 위치만 추출하고
#   샘플 × 마커 형태의 요약 표를 TSV로 저장한다.
#
# 입력:
#   args[1] : variant_AF_all.tsv (파이프라인 결과)
#   args[2] : cds_coords.tsv (resources/)
#   args[3] : 출력 TSV 경로
#   args[4] : AF 임계값 (기본 0.5 — 이 이상이면 변이로 표시)
#
# 출력:
#   샘플 × 마커 wide-format TSV
#   예) barcode | CRT_K76T | DHFR_S108N | K13_C580Y | ...
#       barcode01 | T(1.00) | N(0.99) | WT | ...
#
# 마커 명명 규칙:
#   {GENE}_{ref_aa}{aa_pos}{alt_aa}  ← 문헌 표준 표기 (예: CRT_K76T)
#
# 의존성:
#   R 패키지: dplyr, tidyr (CRAN)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript make_resistance_table.R <variants.tsv> <cds_coords.tsv> <output.tsv> [af_threshold]")
}

variants_path <- args[1]
cds_path      <- args[2]
out_path      <- args[3]
af_threshold  <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.5)

message("[INFO] AF threshold: ", af_threshold)

# =============================================================================
# 1. 알려진 약제내성 마커 정의
#    출처: WHO 2023, MalariaGEN, 문헌 기반 주요 validated/candidate 마커
#    aa_pos: 1-based 아미노산 번호 (PlasmoDB-67 기준)
#    ref_aa: 3D7 참조 아미노산 / alt_aa: 내성 아미노산
# =============================================================================
MARKERS <- tribble(
  ~gene,   ~aa_pos, ~ref_aa, ~alt_aa, ~drug,
  # ── CRT (PF3D7_0709000) — chloroquine ───────────────────────────────────
  # 참조 하플로타입 (3D7): C72-V73-M74-N75-K76 (CVMNK)
  # 내성 하플로타입: CVIET — 72번(C)과 73번(V)은 동일하나 하플로타입 맥락으로 포함
  # C72S: 동남아 일부 분리주에서 보고됨
  # V73I: 드물게 보고됨 (보통 V로 보존)
  "CRT",    72, "C", "S", "CQ",
  "CRT",    73, "V", "I", "CQ",
  "CRT",    74, "M", "I", "CQ",
  "CRT",    75, "N", "E", "CQ",
  "CRT",    76, "K", "T", "CQ",      # 핵심 CQ 내성 마커
  "CRT",   220, "A", "S", "CQ",
  "CRT",   271, "Q", "E", "CQ",
  "CRT",   326, "N", "S", "CQ",
  "CRT",   356, "I", "T", "CQ",
  "CRT",   371, "R", "I", "CQ",

  # ── DHFR (PF3D7_0417200) — pyrimethamine / trimethoprim ─────────────────
  "DHFR",   51, "N", "I", "PYR",
  "DHFR",   59, "C", "R", "PYR",
  "DHFR",  108, "S", "N", "PYR",    # 핵심 PYR 내성 마커
  "DHFR",  164, "I", "L", "PYR",

  # ── DHPS (PF3D7_0810800) — sulfadoxine ──────────────────────────────────
  "DHPS",  436, "S", "A", "SDX",
  "DHPS",  437, "A", "G", "SDX",    # 핵심 SDX 내성 마커
  "DHPS",  540, "K", "E", "SDX",
  "DHPS",  581, "A", "G", "SDX",
  "DHPS",  613, "A", "S", "SDX",

  # ── K13 (PF3D7_1343700) — artemisinin — WHO validated & candidate ────────
  "K13",   446, "F", "I", "ART",
  "K13",   476, "M", "I", "ART",
  "K13",   493, "Y", "H", "ART",
  "K13",   533, "G", "V", "ART",
  "K13",   539, "R", "T", "ART",
  "K13",   543, "I", "T", "ART",
  "K13",   553, "P", "L", "ART",
  "K13",   561, "R", "H", "ART",    # 동남아 주요 마커
  "K13",   574, "P", "L", "ART",
  "K13",   580, "C", "Y", "ART",    # 동아프리카 주요 마커
  "K13",   675, "A", "V", "ART",

  # ── MDR1 (PF3D7_0523000) — lumefantrine / amodiaquine ───────────────────
  "MDR1",   86, "N", "Y", "LUM/AQ",
  "MDR1",  184, "Y", "F", "LUM/AQ",
  "MDR1", 1034, "S", "C", "LUM/AQ",
  "MDR1", 1042, "N", "D", "LUM/AQ",
  "MDR1", 1246, "D", "Y", "LUM/AQ",
)

# =============================================================================
# 2. CDS 좌표 로드
# =============================================================================
cds <- read.table(
  cds_path, header = TRUE, sep = "\t",
  comment.char = "#", stringsAsFactors = FALSE
)
# 컬럼: gene, chrom, cds_start(0-based), cds_end(0-based excl), strand, exon_rank

# =============================================================================
# 3. 게놈 위치(1-based) → AA 위치(1-based) 매핑 테이블 생성
#    + strand : 엑손 오름차순, 엑손 내부도 오름차순
#    - strand : 엑손 내림차순, 엑손 내부도 내림차순 (역방향이 5'→3')
# =============================================================================
build_pos_table <- function(gene_name, cds_df) {
  exons  <- cds_df[cds_df$gene == gene_name, ]
  if (nrow(exons) == 0) return(NULL)
  exons  <- exons[order(exons$cds_start), ]
  strand <- exons$strand[1]

  if (strand == "+") {
    rows <- lapply(seq_len(nrow(exons)), function(i) {
      ex <- exons[i, ]
      # cds_start는 0-based → +1 해서 1-based genomic pos
      data.frame(gpos = seq(ex$cds_start + 1L, ex$cds_end),
                 stringsAsFactors = FALSE)
    })
  } else {
    # - strand: 역순 엑손, 역순 위치
    rows <- lapply(nrow(exons):1, function(i) {
      ex <- exons[i, ]
      data.frame(gpos = seq(ex$cds_end, ex$cds_start + 1L),
                 stringsAsFactors = FALSE)
    })
  }

  positions <- do.call(rbind, rows)
  positions$cds_0based <- seq(0L, nrow(positions) - 1L)
  positions$aa_pos     <- positions$cds_0based %/% 3L + 1L  # 1-based
  positions$gene       <- gene_name
  positions
}

genes_needed <- unique(MARKERS$gene)
pos_tables   <- lapply(genes_needed, build_pos_table, cds_df = cds)
names(pos_tables) <- genes_needed
pos_df <- do.call(rbind, Filter(Negate(is.null), pos_tables))
rownames(pos_df) <- NULL

message("[INFO] Position map built: ", nrow(pos_df), " CDS positions across ",
        length(genes_needed), " genes")

# =============================================================================
# 4. 변이표 로드
# =============================================================================
vars <- read.table(
  variants_path, header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, quote = ""
)
message("[INFO] Variants loaded: ", nrow(vars), " rows, ",
        length(unique(vars$barcode)), " samples")

# =============================================================================
# 5. 게놈 위치 → AA 위치 합류, 마커 위치 변이만 추출
# =============================================================================
vars_at_markers <- vars %>%
  left_join(
    pos_df %>% select(gene, gpos, aa_pos),
    by = c("gene" = "gene", "POS" = "gpos")
  ) %>%
  filter(!is.na(aa_pos)) %>%
  inner_join(
    MARKERS %>% select(gene, aa_pos, ref_aa, alt_aa, drug),
    by = c("gene", "aa_pos")
  ) %>%
  mutate(marker_id = paste0(gene, "_", ref_aa, aa_pos, alt_aa))

# =============================================================================
# 6. 전체 샘플 × 전체 마커 그리드 생성 후 변이 정보 합류
#    변이가 없는 위치는 ref_aa(레퍼런스 아미노산)로 채운다 — "WT" 문자열 사용 안 함
# =============================================================================
all_barcodes <- unique(vars$barcode)
all_markers  <- MARKERS %>%
  mutate(marker_id = paste0(gene, "_", ref_aa, aa_pos, alt_aa))

# 모든 (barcode, marker) 조합의 기준 그리드
grid <- expand.grid(
  barcode   = all_barcodes,
  marker_id = all_markers$marker_id,
  stringsAsFactors = FALSE
) %>%
  left_join(all_markers %>% select(marker_id, ref_aa, alt_aa),
            by = "marker_id")

# 마커별 코돈 내 최대 AF 변이 하나만 (같은 코돈에 여러 SNP이 있을 경우)
best_variant <- vars_at_markers %>%
  group_by(barcode, marker_id) %>%
  slice_max(order_by = AF, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(barcode, marker_id, AF)

# 그리드에 합류
result <- grid %>%
  left_join(best_variant, by = c("barcode", "marker_id")) %>%
  mutate(
    # 변이 없음(NA) → ref_aa 그대로 표시
    # 변이 있고 AF >= threshold → alt_aa(AF)
    # 변이 있고 0.2 <= AF < threshold → ref_aa/alt_aa(AF) (혼합감염)
    # 변이 있고 AF < 0.2 → ref_aa (노이즈 수준)
    call = case_when(
      is.na(AF)              ~ ref_aa,
      AF >= af_threshold     ~ paste0(alt_aa, "(", sprintf("%.2f", AF), ")"),
      AF >= 0.2              ~ paste0(ref_aa, "/", alt_aa, "(", sprintf("%.2f", AF), ")"),
      TRUE                   ~ ref_aa
    )
  )

# =============================================================================
# 7. Wide-format 표 생성 (마커 컬럼 순서는 MARKERS 정의 순서 고정)
# =============================================================================
marker_order <- all_markers$marker_id

wide <- result %>%
  select(barcode, marker_id, call) %>%
  pivot_wider(
    names_from  = marker_id,
    values_from = call
  ) %>%
  select(barcode, all_of(marker_order)) %>%
  arrange(barcode)

# =============================================================================
# 8. 저장
# =============================================================================
write.table(wide, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
message("[INFO] Resistance table saved: ", out_path)
message("[INFO] Samples: ", nrow(wide), "  Markers: ", ncol(wide) - 1)
