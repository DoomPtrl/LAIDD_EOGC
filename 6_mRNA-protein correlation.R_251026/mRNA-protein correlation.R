# =========================
# 1) Libraries
# =========================
library(readr)
library(dplyr)
library(biomaRt)
library(preprocessCore)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(SummarizedExperiment)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 3) preprocessCore 설치 (업데이트/질문 끔)
BiocManager::install("preprocessCore", update = FALSE, ask = FALSE)

# 4) 로드
library(preprocessCore)

# 5) 동작 확인
exists("normalize.quantiles")


options(repos = c(CRAN = "https://cran.rstudio.com"))
dir.create("C:/Rlibs", showWarnings = FALSE)
.libPaths(c("C:/Rlibs", .libPaths()))

# 2) 의존성까지 설치
install.packages("survminer", dependencies = TRUE)

# 3) 로드
library(survminer)

# 4) 확인
packageVersion("survminer")
# =========================
# 2) mRNA data: load → log2(TPM+1) → quantile norm → T/N FC
# =========================
rna_dir      <- "D:/LAIDD/RSEM_Gene/RSEM_Gene"  # folder with *results files
protein_file <- "D:/LAIDD/protein_expression_change.txt"

# list and read
rna_files <- list.files(rna_dir, pattern = "results", full.names = TRUE)
stopifnot(length(rna_files) > 1)

# read into named list
rna_list <- setNames(
  lapply(rna_files, \(f){
    df <- readr::read_delim(f, delim = "\t", show_col_types = FALSE)
    # keep minimal cols
    stopifnot(all(c("gene_id","TPM") %in% names(df)))
    df <- df[, c("gene_id","TPM")]
    df$gene_id <- as.character(df$gene_id)
    df
  }),
  nm = gsub("_rsem_genes_original_results.*$|\\.tsv$|\\.txt$", "", basename(rna_files))
)

# 2-1 Filter: keep genes with TPM >= 1 in >= 30% of samples
long <- dplyr::bind_rows(
  lapply(names(rna_list), \(nm){
    df <- rna_list[[nm]]
    df$sample <- nm
    df
  })
)

keep_ids <- long %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarize(frac = mean(TPM >= 1, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(frac >= 0.30) %>%
  dplyr::pull(gene_id)

rna_list <- lapply(rna_list, \(df) df[df$gene_id %in% keep_ids, c("gene_id","TPM")])

# 2-2 log2(TPM+1)
rna_list <- lapply(
  rna_list,
  \(df){
    df$log2TPM1 <- log2(df$TPM + 1)
    df[, c("gene_id","log2TPM1")]
  }
)

# 2-3 Map Ensembl -> Symbol for the union once
all_gene_ids <- unique(unlist(lapply(rna_list, \(df) df$gene_id)))
ens <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map <- biomaRt::getBM(
  attributes = c("ensembl_gene_id","external_gene_name"),
  filters    = "ensembl_gene_id",
  values     = all_gene_ids,
  mart       = ens
)

# 2-4 Attach gene_name and align rows by common gene set
common_ids <- Reduce(intersect, lapply(rna_list, \(df) df$gene_id))

rna_list <- lapply(
  rna_list,
  \(df){
    df <- df[df$gene_id %in% common_ids, ]
    df <- df[match(common_ids, df$gene_id), ]
    df$gene_name <- map$external_gene_name[match(df$gene_id, map$ensembl_gene_id)]
    df$gene_name[is.na(df$gene_name)] <- NA_character_
    df[, c("gene_id","gene_name","log2TPM1")]
  }
)

# 2-5 Quantile normalization across samples (columns) - 필요에 따라 생략 가능
mat    <- do.call(cbind, lapply(rna_list, \(df) df$log2TPM1))
mat_qn <- preprocessCore::normalize.quantiles(as.matrix(mat))
for (i in seq_along(rna_list)) rna_list[[i]]$log2TPM1 <- mat_qn[, i]

# 2-6 Build paired N/T and compute FC = Tumor - Normal
samps <- names(rna_list)
isN   <- grepl("N$", samps)
isT   <- grepl("T$", samps)
nums  <- suppressWarnings(as.integer(gsub("[^0-9]", "", samps)))

pairs <- data.frame(N = samps[isN], Nnum = nums[isN])
pairs$T <- paste0(pairs$Nnum + 1, "T")
pairs   <- pairs[pairs$T %in% samps, ]
stopifnot(nrow(pairs) > 0)

# function for pair label
pair_label <- \(nnum) paste0("N", nnum, "T", nnum + 1)

# compute FC for each pair
fc_list <- lapply(seq_len(nrow(pairs)), function(k){
  ndf <- rna_list[[pairs$N[k]]]
  tdf <- rna_list[[pairs$T[k]]]
  stopifnot(identical(ndf$gene_id, tdf$gene_id))
  
  lbl      <- pair_label(pairs$Nnum[k])  # ex) "N15T16"
  diff_vec <- tdf$log2TPM1 - ndf$log2TPM1
  
  out <- data.frame(
    gene_id     = ndf$gene_id,
    gene_name   = ndf$gene_name,
    check.names = FALSE
  )
  out[[lbl]] <- diff_vec  # 동적 컬럼 추가
  out
})

# merge into one wide table
merged_df <- Reduce(function(x, y) dplyr::left_join(x, y, by = c("gene_id","gene_name")), fc_list)

# ensure no duplicated rows
merged_df <- merged_df[!duplicated(merged_df$gene_id), ]

# Inspect
dim(merged_df)
head(merged_df[, 1:6])

# =========================
# 3) Proteomics data
# =========================
prot_raw <- read.csv(protein_file, sep = "\t", header = TRUE, check.names = FALSE)
str(prot_raw)

prot_clean <- prot_raw %>%
  dplyr::mutate(
    Symbol = trimws(Symbol),                # Symbol 앞뒤 공백 제거
    Symbol = dplyr::na_if(Symbol, "#N/A")   # "#N/A" 문자열을 NA로 변환
  ) %>%
  dplyr::filter(!is.na(Symbol), Symbol != "")  # 유효한 Symbol만 남김

# log2FC 컬럼 자동 탐지: "N숫자T숫자" 패턴(예: N15T16, N23T24)
fc_cols <- grep("^N\\d+T\\d+$", names(prot_clean), value = TRUE)
stopifnot(length(fc_cols) > 0)  # log2FC 열이 하나도 없으면 에러

# 숫자형으로 강제 변환(문자 섞임 방지)
prot_clean[fc_cols] <- lapply(prot_clean[fc_cols], function(x) suppressWarnings(as.numeric(x)))

# 중복 Symbol 통합(각 페어별 중앙값)
prot_agg <- prot_clean %>%
  dplyr::group_by(Symbol) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(fc_cols),
      ~ if (all(is.na(.x))) NA_real_ else stats::median(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# 커버리지 필터: 환자 30% 이상에서 값 존재하는 단백질만 유지
n_pairs       <- length(fc_cols)
target_frac   <- 0.30
min_present_n <- ceiling(target_frac * n_pairs)

prot_keep <- prot_agg %>%
  dplyr::mutate(
    present_cnt = rowSums(!is.na(dplyr::across(dplyr::all_of(fc_cols))))
  ) %>%
  dplyr::filter(present_cnt >= min_present_n) %>%
  dplyr::select(-present_cnt)

# 최종 결과(예시 행 수 주석): 중복 통합(중앙값) + 커버리지 조건 만족 테이블  # ~7329
nrow(prot_keep)


library(dplyr)

## 3.5) Proteomics wide (Symbol + log2FC 페어 컬럼들만)
prot_wide <- prot_keep %>%
  dplyr::arrange(Symbol) %>%
  dplyr::select(Symbol, dplyr::all_of(fc_cols))

## 4) RNA: gene_name 기준 중복을 중앙값으로 집계
rna_pairs_all <- setdiff(names(merged_df), c("gene_id","gene_name"))

# 숫자 강제(문자형 섞였을 가능성 방지)
tmp_rna <- merged_df
tmp_rna[rna_pairs_all] <- lapply(tmp_rna[rna_pairs_all], function(x) suppressWarnings(as.numeric(x)))

rna_fc_by_gene <- tmp_rna %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(rna_pairs_all),
      ~ if (all(is.na(.))) NA_real_ else stats::median(., na.rm = TRUE)
    ),
    .groups = "drop"
  )

## 5) 유전자 및 페어 교집합 구하기
genes_intersect <- intersect(rna_fc_by_gene$gene_name, prot_wide$Symbol)
if (length(genes_intersect) == 0) stop("No overlapping genes between RNA and Proteomics after filtering.")

prot_pairs <- setdiff(names(prot_wide), "Symbol")
rna_pairs  <- setdiff(names(rna_fc_by_gene), "gene_name")
common_pairs <- intersect(rna_pairs, prot_pairs)
if (length(common_pairs) == 0) {
  stop("No common sample-pair columns between RNA and Proteomics. Check pair naming (e.g., N15T16).")
}

## 6) 교집합으로 정렬·정합 (행: 유전자, 열: 공통 페어)
prot_wide2 <- prot_wide %>%
  dplyr::filter(Symbol %in% genes_intersect) %>%
  dplyr::arrange(Symbol) %>%
  dplyr::select(Symbol, dplyr::all_of(common_pairs))

rna_fc2 <- rna_fc_by_gene %>%
  dplyr::filter(gene_name %in% genes_intersect) %>%
  dplyr::arrange(gene_name) %>%
  dplyr::select(gene_name, dplyr::all_of(common_pairs))

stopifnot(identical(rna_fc2$gene_name, prot_wide2$Symbol))

## 7) 행렬 변환 (유전자 × 페어), 숫자화 보장
rna_mat  <- as.matrix(rna_fc2[, common_pairs, drop = FALSE])
prot_mat <- as.matrix(prot_wide2[, common_pairs, drop = FALSE])
mode(rna_mat) <- "numeric"
mode(prot_mat) <- "numeric"

## 8) 유전자 벡터
gene_symbols <- rna_fc2$gene_name
stopifnot(length(gene_symbols) == nrow(rna_mat))

## 9) Spearman 상관 (유전자별, 공통 페어 축)
ng <- nrow(rna_mat)
rho  <- numeric(ng); rho[]  <- NA_real_
pval <- numeric(ng); pval[] <- NA_real_

for (i in seq_len(ng)) {
  xi <- as.numeric(rna_mat[i, ])
  yi <- as.numeric(prot_mat[i, ])
  ok <- is.finite(xi) & is.finite(yi)  # 완전쌍만 사용
  
  # 유효쌍이 최소 3개, 양쪽 모두 변동이 있어야 테스트 가능
  if (sum(ok) >= 3 && length(unique(xi[ok])) > 1 && length(unique(yi[ok])) > 1) {
    ct <- suppressWarnings(stats::cor.test(xi[ok], yi[ok], method = "spearman"))
    rho[i]  <- unname(ct$estimate)
    pval[i] <- ct$p.value
  }
}

## 10) 요약 통계 + FDR
overall_mean <- mean(rho, na.rm = TRUE)
pct_pos      <- mean(rho > 0, na.rm = TRUE) * 100
fdr          <- p.adjust(pval, method = "BH")
pct_sig_pos  <- mean((rho > 0) & (fdr < 0.01), na.rm = TRUE) * 100

summary_list <- list(
  n_genes                     = ng,
  n_pairs                     = ncol(rna_mat),
  overall_mean_rho            = overall_mean,
  percent_positive            = pct_pos,
  percent_positive_FDR_lt_0_01 = pct_sig_pos,
  n_sig_pos_FDR_lt_0_01       = sum((rho > 0) & (fdr < 0.01), na.rm = TRUE)
)
print(summary_list)

## 11) 결과 테이블(유전자별 상관)
cor_df <- data.frame(
  gene = gene_symbols,
  rho  = rho,
  p    = pval,
  fdr  = fdr,
  stringsAsFactors = FALSE
)

# 참고: 상위/하위 몇 개 확인
head(cor_df[order(cor_df$rho, decreasing = TRUE), ], 10)
tail(cor_df[order(cor_df$rho, decreasing = TRUE), ], 10)

# --- Final gene symbols vector (rna_fc2 기준) ---
gene_symbols <- rna_fc2$gene_name
stopifnot(length(gene_symbols) == nrow(rna_mat),
          identical(gene_symbols, prot_wide2$Symbol))

# --- 4) Spearman correlation (per gene across common pairs) ---
ng <- nrow(rna_mat)
rho  <- rep(NA_real_, ng)
pval <- rep(NA_real_, ng)

for (i in seq_len(ng)) {
  xi <- as.numeric(rna_mat[i, ])
  yi <- as.numeric(prot_mat[i, ])
  ok <- is.finite(xi) & is.finite(yi)  # 완전쌍만 사용
  
  # 유효쌍 ≥ 3, 양쪽 모두 변동(상수 아님)일 때만 테스트
  if (sum(ok) >= 3 && length(unique(xi[ok])) > 1 && length(unique(yi[ok])) > 1) {
    ct <- suppressWarnings(stats::cor.test(xi[ok], yi[ok], method = "spearman"))
    rho[i]  <- unname(ct$estimate)
    pval[i] <- ct$p.value
  }
}

overall_mean <- mean(rho, na.rm = TRUE)
pct_pos      <- mean(rho > 0, na.rm = TRUE) * 100
fdr          <- p.adjust(pval, method = "BH")
pct_sig_pos  <- mean((rho > 0) & (fdr < 0.01), na.rm = TRUE) * 100

list(
  n_genes                      = ng,
  n_pairs                      = ncol(rna_mat),
  overall_mean_rho             = overall_mean,
  percent_positive             = pct_pos,
  percent_positive_FDR_lt_0_01 = pct_sig_pos
)

# --- Final gene symbols vector (rna_fc2 기준) ---
gene_symbols <- rna_fc2$gene_name
stopifnot(length(gene_symbols) == nrow(rna_mat),
          identical(gene_symbols, prot_wide2$Symbol))

# --- 4) Spearman correlation (per gene across common pairs) ---
ng <- nrow(rna_mat)
rho  <- rep(NA_real_, ng)
pval <- rep(NA_real_, ng)

for (i in seq_len(ng)) {
  xi <- as.numeric(rna_mat[i, ])
  yi <- as.numeric(prot_mat[i, ])
  ok <- is.finite(xi) & is.finite(yi)  # 완전쌍만 사용
  
  # 유효쌍 ≥ 3, 양쪽 모두 변동(상수 아님)일 때만 테스트
  if (sum(ok) >= 3 && length(unique(xi[ok])) > 1 && length(unique(yi[ok])) > 1) {
    ct <- suppressWarnings(stats::cor.test(xi[ok], yi[ok], method = "spearman"))
    rho[i]  <- unname(ct$estimate)
    pval[i] <- ct$p.value
  }
}

overall_mean <- mean(rho, na.rm = TRUE)
pct_pos      <- mean(rho > 0, na.rm = TRUE) * 100
fdr          <- p.adjust(pval, method = "BH")
pct_sig_pos  <- mean((rho > 0) & (fdr < 0.01), na.rm = TRUE) * 100

## 10) Summary stats
summary_list <- list(
  n_genes                      = ng,
  n_pairs                      = ncol(rna_mat),
  overall_mean_rho             = overall_mean,
  percent_positive             = pct_pos,
  percent_positive_FDR_lt_0_01 = pct_sig_pos,
  n_sig_pos_FDR_lt_0_01        = sum((rho > 0) & (fdr < 0.01), na.rm = TRUE)
)
print(summary_list)

## 11) Gene-level correlation table
cor_df <- data.frame(
  gene = gene_symbols,
  rho  = rho,
  p    = pval,
  fdr  = fdr,
  stringsAsFactors = FALSE
)

# quick peek: top/bottom 10 by rho
head(cor_df[order(cor_df$rho, decreasing = TRUE), ], 10)
tail(cor_df[order(cor_df$rho, decreasing = TRUE), ], 10)

## Figure 2A — Distribution of Spearman's ρ (robust)

stopifnot(exists("rho"))  # rho가 있어야 함

# 1) 유한값만 사용
rho_fin <- rho[is.finite(rho)]
if (length(rho_fin) < 2) {
  stop(sprintf("Not enough finite rho values to plot (n=%d).", length(rho_fin)))
}

# 2) bin 설정 (엣지케이스 방어)
bw <- 0.035
rng <- range(rho_fin)
start <- floor(rng[1] / bw) * bw
end   <- ceiling(rng[2] / bw) * bw
if (start == end) {  # 모든 값이 동일할 때 폭 확보
  start <- start - bw
  end   <- end + bw
}
brks <- seq(start, end, by = bw)
if (length(brks) < 2L) {
  brks <- c(start, start + bw)  # 최소 2개 보장
}

# 3) 히스토그램(계산만)
h <- hist(rho_fin, breaks = brks, plot = FALSE)

# 4) 히스토그램 데이터 프레임
dfh <- data.frame(
  mid     = h$mids,
  density = h$density
)
dfh$fill <- ifelse(dfh$mid < 0, "neg", "pos")

# 5) 레이블 텍스트 (앞서 계산한 overall_mean, pct_pos, pct_sig_pos 활용)
ymax <- max(dfh$density, na.rm = TRUE)
mean_txt <- sprintf("Mean = %.3f", overall_mean)
sub_txt  <- sprintf("Positive = %.2f%%   |   Positive & FDR<0.01 = %.2f%%",
                    pct_pos, pct_sig_pos)

# 6) ggplot (버전 호환을 위해 ..density.. 사용)
library(ggplot2)

p <- ggplot(dfh, aes(x = mid, y = density)) +
  # 막대 (음수=주황, 양수=파랑)
  geom_col(aes(fill = fill), width = bw, color = "black", alpha = 0.85) +
  scale_fill_manual(values = c(neg = "orange3", pos = "darkblue"), guide = "none") +
  # 밀도 곡선 (linewidth / after_stat 사용)
  geom_density(
    data = data.frame(rho = rho_fin),
    aes(x = rho, y = after_stat(density)),
    color = "red3", linewidth = 0.9
  ) +
  # 평균선 (linewidth 사용)
  geom_vline(xintercept = overall_mean, linetype = "dashed",
             color = "green3", linewidth = 1) +
  # 평균 텍스트
  annotate("text",
           x = overall_mean, y = ymax * 0.95,
           label = mean_txt, color = "green3", fontface = "bold", vjust = 1) +
  labs(x = "Spearman's correlation coefficient",
       y = "Probability density",
       subtitle = sub_txt) +
  coord_cartesian(xlim = range(dfh$mid), ylim = c(0, ymax * 1.05)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

print(p)

