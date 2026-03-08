library(ape)
library(ggtree)
library(ggrepel)
library(ggplot2)

TREE_FILE <- "./consensus_Pf3D7_02_v3/aligned_MSP2.nwk"
N_PERM <- 999
K_CLADE <- 6
set.seed(123)

tree <- read.tree(TREE_FILE)
labels <- tree$tip.label

# ---- robust parser ----
extract_meta <- function(label){
  # 숨은 문자 제거 (윈도우 CR 등)
  label <- gsub("\r", "", label)
  
  # 마지막 4자리 연도 추출 (끝에 고정하지 않음)
  m <- regmatches(label, regexpr("[0-9]{4}(?!.*[0-9]{4})", label, perl=TRUE))
  if (length(m)==0 || nchar(m)==0) stop("Cannot find 4-digit year in label: ", label)
  year <- as.numeric(m)
  
  # year 및 앞의 공백/탭 제거해서 core 만들기
  core <- sub("[ \t]*[0-9]{4}(?!.*[0-9]{4}).*$", "", label, perl=TRUE)
  
  parts <- strsplit(core, "\\|")[[1]]
  sample <- parts[1]
  country <- paste(parts[-1], collapse="|")
  
  # country에 붙은 탭/공백 정리
  country <- gsub("[ \t]+$", "", country)
  
  c(sample=sample, country=country, year=year)
}

meta_mat <- t(sapply(labels, extract_meta))
meta <- data.frame(
  Sample  = meta_mat[, "sample"],
  Country = meta_mat[, "country"],
  Year    = as.numeric(meta_mat[, "year"]),
  stringsAsFactors = FALSE
)
rownames(meta) <- labels

# ---- distance matrices ----
genetic_dist <- cophenetic.phylo(tree)
meta <- meta[rownames(genetic_dist), , drop=FALSE]

country_dist <- as.dist(1 - outer(meta$Country, meta$Country, FUN="=="))
year_dist    <- as.dist(abs(outer(meta$Year, meta$Year, FUN="-")))

# ---- Mantel tests ----
mantel_country <- mantel(genetic_dist, country_dist, permutations=N_PERM)
mantel_year    <- mantel(genetic_dist, year_dist, permutations=N_PERM)
mantel_year_partial <- mantel.partial(genetic_dist, year_dist, country_dist, permutations=N_PERM)

mantel_country
mantel_year
mantel_year_partial

# ---- clade permutation test ----
hc <- hclust(as.dist(cophenetic.phylo(tree)), method = "average")
clade <- cutree(hc, k = K_CLADE)

clade_country_stat <- function(clade, country){
  vals <- c()
  for(k in unique(clade)){
    idx <- which(clade==k)
    if(length(idx) < 3) next
    same <- outer(country[idx], country[idx], FUN="==")
    vals <- c(vals, mean(same[upper.tri(same)]))
  }
  mean(vals)
}

obs <- clade_country_stat(clade, meta$Country)
perm <- replicate(2000, clade_country_stat(clade, sample(meta$Country)))
p_perm <- mean(perm >= obs)

p_perm

# ---- save outputs ----
sink("mantel_results.txt")
cat("Mantel: genetic vs country\n"); print(mantel_country); cat("\n")
cat("Mantel: genetic vs year\n"); print(mantel_year); cat("\n")
cat("Partial Mantel: year | country\n"); print(mantel_year_partial); cat("\n")
cat("Clade-country permutation p =", p_perm, "\n")
sink()

pdf("MSP2_tree_country.pdf", width=10, height=10)
country_factor <- factor(meta$Country)
country_factor <- factor(meta$Country)

plot(
  tree,
  type = "phylogram",
  show.tip.label = FALSE,   # 🔥 라벨 제거
  tip.color = as.integer(country_factor),
  cex = 0.5
)

legend(
  "topleft",
  legend = levels(country_factor),
  col = seq_along(levels(country_factor)),
  pch = 19,
  cex = 0.6,
  bty = "n"
)

p2 <- ggtree(tree, layout = "rectangular") +
  geom_tiplab_repel(size = 1.6, max.overlaps = Inf)

ggsave("MSP2_tree_A4_rectangular.pdf", p2, width = 11.69, height = 8.27, units = "in")



dev.off()
