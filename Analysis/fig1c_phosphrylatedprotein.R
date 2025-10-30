---
title: "fig1c" - nonsynonymous mutated gene and phsophorylated protein relation
author: "Yujung Ahn"
date: "2025-09-09"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(data.table)
# ========== 준비 ========== #
# 1. 데이터 불러오기
pprotein <- fread("pprotein.txt", na.strings = c("NaN","NA"))

# 2. NaN → 0
pprotein[is.na(pprotein)] <- 0

# 3. mut/wt 그룹 설정 - ARDI1A gene
mut_samples <- c("N63T64", "N111T112", "N117T118")

all_samples <- setdiff(colnames(pprotein), c("Peptide","Symbol"))
wt_samples  <- setdiff(all_samples, mut_samples)

# 4. T-test .. data.table이라서 .. 붙여야 함
t_results <- apply(pprotein[, ..all_samples], 1, function(x){
  mut <- as.numeric(x[mut_samples])
  wt  <- as.numeric(x[wt_samples])
  if (length(unique(mut)) > 1 | length(unique(wt)) > 1) {
    t_res <- t.test(mut, wt)
    return(c(mean_mut = mean(mut), mean_wt = mean(wt),
             pval = t_res$p.value))
  } else {
    return(c(mean_mut = mean(mut), mean_wt = mean(wt), pval = NA))
  }
})

# 5. t_results를 data.frame으로 변환
t_results <- as.data.frame(t(t_results))
t_results$Symbol <- pprotein$Symbol  # 유전자 이름 붙이기

# 6. p < 0.05인 심볼만 추출
sig_symbols <- t_results$Symbol[!is.na(t_results$pval) & t_results$pval < 0.05]

# 7. 해당 심볼 행만 뽑아서 매트릭스 만들기
heatmap_mat <- pprotein[pprotein$Symbol %in% sig_symbols, ..all_samples]
heatmap_mat <- as.matrix(heatmap_mat)
rownames(heatmap_mat) <- pprotein$Symbol[pprotein$Symbol %in% sig_symbols]

# 8. Mut / WT annotation
ann_col <- data.frame(
  Group = ifelse(colnames(heatmap_mat) %in% mut_samples, "Mut", "WT")
)
rownames(ann_col) <- colnames(heatmap_mat)



# Mut / WT samples 명시적 정의
wt_samples <- setdiff(colnames(heatmap_mat), mut_samples)

# Mut → WT 순서로 column 정렬
col_order <- c(mut_samples, wt_samples)
heatmap_mat <- heatmap_mat[, col_order]   # BUG?

# annotation_col을 col_order 기준으로 새로 생성
ann_col <- data.frame(
  Group = ifelse(col_order %in% mut_samples, "Mut", "WT")
)
rownames(ann_col) <- col_order

# 색상 팔레트
my_colors <- colorRampPalette(c("green", "black", "red"))(100)

# Mut/WT annotation 색상 지정
ann_colors <- list(
  Group = c("Mut" = "black", "WT" = "white")
)

# Mut, WT 평균 구하기
mut_mean <- rowMeans(heatmap_mat[, mut_samples, drop = FALSE])
wt_mean  <- rowMeans(heatmap_mat[, wt_samples, drop = FALSE])

# 조건 만족하는 row만 추출
keep_rows <- (mut_mean > 0) & (wt_mean < 0)
heatmap_mat <- heatmap_mat[keep_rows, ]

# breakpoints 정의 (-1.5 ~ 1.5)
breaks_list <- seq(-1.5, 1.5, length.out = 101)


# Heatmap
library(pheatmap)
         
pheatmap(heatmap_mat,
         scale = "row",
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         use_raster = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         legend = FALSE,
         main = "RHOA",
         color = my_colors,
         breaks = breaks_list)


## 드디어 마침내 해내다... 농담곰의 축복덕분이야...
```

