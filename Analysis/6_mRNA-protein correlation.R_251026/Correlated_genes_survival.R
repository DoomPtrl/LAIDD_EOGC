# =========================================================
# TCGA-STAD survival analysis (Q1 vs Q4) + exports + KM plots (top 5)
# =========================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(biomaRt)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(cowplot)   # 중앙 정렬용
})

## ---------------- I/O ----------------
out_dir    <- "D:/LAIDD"                     # ← 필요시 변경
tables_dir <- file.path(out_dir, "tables")
plots_dir  <- file.path(out_dir, "plots")
dir.create(out_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")

## ---------------- 0) 입력 존재 확인 ----------------
if (!exists("cor_df") || !all(c("gene","rho","fdr") %in% names(cor_df))) {
  stop("cor_df(gene, rho, fdr)가 메모리에 필요합니다. (mRNA–protein 상관 분석 표)")
}

## ---------------- 1) TCGA-STAD: expression (HTSeq - FPKM-UQ) ----------------
# NOTE: 'STAR - Counts'가 아니라 FPKM-UQ를 다운받아 log2(FPKM-UQ+1) 분석합니다.
query <- GDCquery(
  project       = "TCGA-STAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "HTSeq - FPKM-UQ",
  sample.type   = "Primary Tumor"
)
GDCdownload(query)
se <- GDCprepare(query)  # SummarizedExperiment

# rows: Ensembl IDs (with version); cols: barcodes
expr <- assay(se) %>% as.data.frame()
expr$ensembl <- rownames(expr)
expr$ensembl_base <- sub("\\..*$", "", expr$ensembl)

# Ensembl → SYMBOL
ens <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapE <- biomaRt::getBM(
  attributes = c("ensembl_gene_id","external_gene_name"),
  filters    = "ensembl_gene_id",
  values     = unique(expr$ensembl_base),
  mart       = ens
)
expr$SYMBOL <- mapE$external_gene_name[match(expr$ensembl_base, mapE$ensembl_gene_id)]
expr <- expr[!is.na(expr$SYMBOL) & nzchar(expr$SYMBOL), , drop = FALSE]

val_cols <- setdiff(colnames(expr), c("ensembl","ensembl_base","SYMBOL"))
mat <- as.matrix(expr[, val_cols, drop = FALSE])
rownames(mat) <- expr$SYMBOL

# 동일 SYMBOL 다수 → 중앙값 집계 후 log2(FPKM-UQ+1)
mat <- apply(mat, 2, function(x) tapply(x, INDEX = rownames(mat), FUN = median, na.rm = TRUE))
mat <- t(mat); mat <- t(mat)   # 행=SYMBOL, 열=샘플 보장
mode(mat) <- "numeric"
mat <- log2(mat + 1)

# 환자 ID(바코드 앞 12자)로 평균 → patient-level matrix
col_patient <- substr(colnames(mat), 1, 12)
mat_df <- as.data.frame(
  t(apply(mat, 1, function(v) tapply(as.numeric(v), col_patient, mean, na.rm = TRUE)))
)

## ---------------- 2) Clinical (OS) from colData(se) ----------------
cd <- as.data.frame(SummarizedExperiment::colData(se))
to_char <- function(x) as.character(x)
to_num  <- function(x) suppressWarnings(as.numeric(to_char(x)))
pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) df[[hit[1]]] else rep(NA, nrow(df))
}

patient <- if ("submitter_id" %in% names(cd)) {
  to_char(cd$submitter_id)
} else if ("case_submitter_id" %in% names(cd)) {
  to_char(cd$case_submitter_id)
} else if ("barcode" %in% names(cd)) {
  substr(to_char(cd$barcode), 1, 12)
} else {
  substr(colnames(se), 1, 12)
}

d_death <- to_num(pick_col(cd, c("days_to_death","days_to_death.x","days_to_death.y")))
d_fu    <- to_num(pick_col(cd, c("days_to_last_follow_up","days_to_lastfollowup",
                                 "days_to_last_follow_up.x","days_to_last_follow_up.y")))
vital   <- to_char(pick_col(cd, c("vital_status","vital_status.x","vital_status.y")))

OS_time  <- ifelse(!is.na(d_death), d_death, d_fu)
OS_event <- ifelse(tolower(vital) == "dead", 1, 0)

surv_df <- data.frame(
  patient  = patient,
  OS_time  = OS_time,
  OS_event = OS_event,
  stringsAsFactors = FALSE
) %>% filter(is.finite(OS_time), is.finite(OS_event))

# 교집합 환자만 정렬
common_pat <- intersect(colnames(mat_df), surv_df$patient)
mat_df <- mat_df[, common_pat, drop = FALSE]
surv_df <- surv_df[match(common_pat, surv_df$patient), , drop = FALSE]
stopifnot(identical(colnames(mat_df), surv_df$patient))

## ---------------- 3) 유전자별 생존분석 (Top25% vs Bottom25%) ----------------
genes_test <- intersect(unique(cor_df$gene), rownames(mat_df))
cat("Genes to test:", length(genes_test), "\n")

surv_by_gene <- function(exp_vec, surv_df) {
  q25 <- stats::quantile(exp_vec, 0.25, na.rm = TRUE)
  q75 <- stats::quantile(exp_vec, 0.75, na.rm = TRUE)
  grp <- rep(NA_integer_, length(exp_vec))
  grp[exp_vec <= q25] <- 0  # Low
  grp[exp_vec >= q75] <- 1  # High
  ok <- which(!is.na(grp) & is.finite(surv_df$OS_time) & is.finite(surv_df$OS_event))
  if (length(ok) < 10 || length(unique(grp[ok])) < 2) return(c(NA, NA, NA, NA))
  s <- survival::survdiff(Surv(OS_time, OS_event) ~ grp, data = data.frame(surv_df, grp = grp)[ok, ])
  chisq <- unname(s$chisq)
  p_lr  <- 1 - pchisq(chisq, df = 1)
  fit <- try(survival::coxph(Surv(OS_time, OS_event) ~ grp, data = data.frame(surv_df, grp = grp)[ok, ]), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(p_lr, chisq, NA, NA))
  hr <- suppressWarnings(exp(coef(fit)))
  p_hr <- suppressWarnings(summary(fit)$coefficients[,"Pr(>|z|)"])
  c(p_lr, chisq, hr, p_hr)
}

res_mat <- matrix(NA_real_, nrow = length(genes_test), ncol = 4,
                  dimnames = list(genes_test, c("p_logrank","chisq","HR_high_vs_low","p_cox")))
for (i in seq_along(genes_test)) {
  g <- genes_test[i]
  res_mat[i, ] <- surv_by_gene(as.numeric(mat_df[g, ]), surv_df)
}
surv_res <- as.data.frame(res_mat)
surv_res$gene <- rownames(surv_res)

# 라벨링
surv_res <- surv_res %>%
  mutate(
    sig_surv  = ifelse(p_logrank < 0.05, TRUE, FALSE),
    direction = dplyr::case_when(
      sig_surv & HR_high_vs_low < 1 ~ "positive",  # High expr → 낮은 hazard
      sig_surv & HR_high_vs_low > 1 ~ "negative",  # High expr → 높은 hazard
      TRUE ~ "ns"
    )
  )

## ---------------- 4) 저장 (테이블 & 그림) ----------------
# 전체 테이블
fn_all <- file.path(tables_dir, paste0("TCGA_STAD_survival_by_gene_", timestamp, ".csv"))
write.csv(surv_res, fn_all, row.names = FALSE, fileEncoding = "UTF-8-BOM")

# Fig 2C: rho by survival significance
df2c <- cor_df %>%
  dplyr::inner_join(surv_res[, c("gene","sig_surv")], by = "gene") %>%
  dplyr::mutate(group = ifelse(sig_surv, "Survival sig (p<0.05)", "Not sig"))
mean_rho_sig <- mean(df2c$rho[df2c$sig_surv], na.rm = TRUE)
mean_rho_nos <- mean(df2c$rho[!df2c$sig_surv], na.rm = TRUE)
wil_p <- wilcox.test(rho ~ sig_surv, data = df2c)$p.value
p_2c <- ggplot(df2c, aes(x = group, y = rho, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  labs(x = NULL, y = "mRNA–protein correlation (rho)",
       subtitle = sprintf("Wilcoxon p = %.2e", wil_p)) +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out_dir, paste0("Fig2C_rho_by_survivalSig_", timestamp, ".png")),
       p_2c, width = 5.5, height = 4.5, dpi = 300)

# Fig S2A: log-rank chi-square by MP significance
df_s2a <- cor_df %>%
  dplyr::mutate(mp_grp = dplyr::case_when(
    fdr < 0.01 ~ "FDR<0.01",
    fdr > 0.10 ~ "FDR>0.10",
    TRUE ~ NA_character_
  )) %>%
  dplyr::inner_join(surv_res[, c("gene","chisq")], by = "gene") %>%
  dplyr::filter(!is.na(mp_grp))
wil_p2 <- wilcox.test(chisq ~ mp_grp, data = df_s2a)$p.value
p_s2a <- ggplot(df_s2a, aes(x = mp_grp, y = chisq, fill = mp_grp)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  labs(x = NULL, y = "Log-rank chi-square",
       subtitle = sprintf("Wilcoxon p = %.2e", wil_p2)) +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out_dir, paste0("FigS2A_chisq_by_MPsignif_", timestamp, ".png")),
       p_s2a, width = 5.5, height = 4.5, dpi = 300)

# 유의 유전자 리스트(논문식): MP FDR<0.01 ∩ survival p<0.05 (HR 방향에 따라 분리)
sig_mp <- subset(cor_df, fdr < 0.01, select = c(gene))
sel <- surv_res %>% dplyr::inner_join(sig_mp, by = "gene")
pos_genes <- sel$gene[sel$direction == "positive" & sel$sig_surv]
neg_genes <- sel$gene[sel$direction == "negative" & sel$sig_surv]
write.table(pos_genes,
            file = file.path(tables_dir, paste0("genes_survival_positive_", timestamp, ".txt")),
            row.names = FALSE, col.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
write.table(neg_genes,
            file = file.path(tables_dir, paste0("genes_survival_negative_", timestamp, ".txt")),
            row.names = FALSE, col.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

## ---------------- 5) KM helper + 상위 5개 자동 저장 ----------------
# 중앙 정렬/오른쪽 여백 제거 버전
km_plot_gene_pretty2 <- function(gene_symbol, mat_df, surv_df,
                                 out_dir = getwd(), timestamp = format(Sys.time(), "%Y%m%d_%H%M"),
                                 q = 0.25, width = 7, height = 6, dpi = 320,
                                 table_height = 0.22, title_y = "") {
  stopifnot(gene_symbol %in% rownames(mat_df))
  ex  <- as.numeric(mat_df[gene_symbol, ])
  q25 <- stats::quantile(ex, q,     na.rm = TRUE)
  q75 <- stats::quantile(ex, 1 - q, na.rm = TRUE)
  grp <- rep(NA_integer_, length(ex)); grp[ex <= q25] <- 0; grp[ex >= q75] <- 1
  
  dat <- data.frame(OS_time = surv_df$OS_time, OS_event = surv_df$OS_event,
                    grp = factor(grp, levels = c(0,1), labels = c("Low (≤25%)","High (≥75%)")))
  dat <- dat[!is.na(dat$grp) & is.finite(dat$OS_time) & is.finite(dat$OS_event), , drop = FALSE]
  if (nrow(dat) < 10 || length(unique(dat$grp)) < 2) return(invisible(NULL))
  
  fit   <- survival::survfit(Surv(OS_time, OS_event) ~ grp, data = dat)
  sdiff <- survival::survdiff(Surv(OS_time, OS_event) ~ grp, data = dat)
  p_lr  <- 1 - pchisq(unname(sdiff$chisq), df = 1)
  cfit  <- survival::coxph(Surv(OS_time, OS_event) ~ grp, data = dat)
  hr    <- unname(exp(coef(cfit))); hr_ci <- unname(exp(confint(cfit)))
  p_cox <- unname(summary(cfit)$coefficients[, "Pr(>|z|)"])
  sub_txt <- sprintf("HR (High vs Low) = %.2f [%.2f, %.2f]; Cox p = %.2e",
                     hr, hr_ci[1], hr_ci[2], p_cox)
  
  tmax <- max(dat$OS_time, na.rm = TRUE)
  brk  <- if (tmax > 3000) 365 else 180
  
  gp <- survminer::ggsurvplot(
    fit, data = dat,
    risk.table = TRUE, risk.table.height = table_height,
    risk.table.y.text = TRUE, risk.table.y.text.col = TRUE,
    tables.theme = survminer::theme_cleantable() %+replace% ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 9, margin = ggplot2::margin(r = 4))
    ),
    conf.int = FALSE,
    legend = "top", legend.title = NULL,
    legend.labs = levels(dat$grp),
    break.time.by = brk,                         # xlim은 아래 scale_x에서 통일
    pval = sprintf("Log-rank p = %.2e", p_lr),
    pval.coord = c(tmax * 0.02, 0.10),
    palette = c("#E69F00", "#0072B2"),
    ggtheme = ggplot2::theme_classic(base_size = 11)
  )
  
  # 오른쪽 공백 제거 + 좌우 대칭 마진
  gp$plot <- gp$plot +
    ggplot2::labs(title = gene_symbol, subtitle = sub_txt, x = "Days", y = title_y) +
    ggplot2::scale_x_continuous(limits = c(0, tmax), expand = c(0, 0)) +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 10, 2, 10), legend.position = "top")
  gp$table <- gp$table +
    ggplot2::scale_x_continuous(limits = c(0, tmax), expand = c(0, 0)) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 10, 6, 10))
  
  # 두 패널 폭을 맞춰 중앙 정렬
  combo <- cowplot::plot_grid(
    gp$plot, gp$table,
    ncol = 1,
    rel_heights = c(1 - table_height, table_height),
    align = "v", axis = "lr"
  )
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fn <- file.path(out_dir, sprintf("KM_%s_%s.png", make.names(gene_symbol), timestamp))
  ggplot2::ggsave(fn, combo, width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", fn)
  invisible(fn)
}

# 상위 5개 (p_logrank 오름차순) 저장
top_km <- surv_res %>%
  dplyr::filter(is.finite(p_logrank)) %>%
  dplyr::arrange(p_logrank) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::pull(gene)

invisible(lapply(
  top_km,
  function(g) try(km_plot_gene_pretty2(g, mat_df, surv_df,
                                       out_dir = plots_dir,
                                       timestamp = timestamp),
                  silent = TRUE)
))

message("\nSaved tables in: ", normalizePath(tables_dir, winslash = "/"))
message("Saved figures in: ", normalizePath(out_dir, winslash = "/"))
message("Saved KM plots (top 5 by p_logrank) in: ", normalizePath(plots_dir, winslash = "/"))
