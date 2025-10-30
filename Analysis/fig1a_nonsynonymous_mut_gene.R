---
title: "mut in wes - fig1a. nonsynonymous mutated gene"
author: "Yujung Ahn"
date: "2025-09-05"
output: html_document
---

```{r, echo=FALSE}
library(data.table)
library(dplyr)
library(tidyr)
library(dndscv)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

```{r, echo=FALSE}

wes <- fread("/home/yujungahn/laidd/wes.txt", stringsAsFactors = FALSE)

```


```{r}
wes <- fread("/home/yujungahn/laidd/wes.txt")

# End_Position 칼럼 생성 (Locus 복사)
wes$End_Position <- wes$Locus

# 컬럼 이름 MutSigCV 맞춤
colnames(wes)[colnames(wes) == "Official Symbol"] <- "Hugo_Symbol"
colnames(wes)[colnames(wes) == "Patient"] <- "Tumor_Sample_Barcode"
colnames(wes)[colnames(wes) == "Chr"] <- "Chromosome"
colnames(wes)[colnames(wes) == "Locus"] <- "Start_Position"
colnames(wes)[colnames(wes) == "From"] <- "Reference_Allele"
colnames(wes)[colnames(wes) == "To"] <- "Tumor_Seq_Allele2"

# Variant_Classification 생성
wes$Variant_Classification <- ifelse(grepl("nonsynonymous", wes$Type, ignore.case = TRUE), "Missense_Mutation",
                             ifelse(grepl("stopgain", wes$Type, ignore.case = TRUE), "Nonsense_Mutation",
                             ifelse(grepl("stoploss", wes$Type, ignore.case = TRUE), "Nonstop_Mutation",
                             ifelse(grepl("frameshift", wes$Type, ignore.case = TRUE) & grepl("ins", wes$Type, ignore.case = TRUE), "Frame_Shift_Ins",
                             ifelse(grepl("frameshift", wes$Type, ignore.case = TRUE) & grepl("del", wes$Type, ignore.case = TRUE), "Frame_Shift_Del",
                             "In_Frame_Ins")))))

# MutSigCV 필수 칼럼만 선택 & 순서 맞춤
wes_maf <- wes %>%
  select(Hugo_Symbol,
         Tumor_Sample_Barcode,
         Chromosome,
         Start_Position,
         End_Position,
         Reference_Allele,
         Tumor_Seq_Allele2,
         Variant_Classification)

# 저장
fwrite(wes_maf, "/home/yujungahn/laidd/wes_for_mutsig_maf.maf",
       sep = "\t", quote = FALSE, na = "")


```



0. mut count 가 높은 gene만 추리기
```{r}
library(data.table)
library(dplyr)

# 파일 불러오기
wes <- fread("/home/yujungahn/laidd/wes.txt")

# nonsynonymous / stopgain / frameshift만 선택
wes_sig <- wes %>% 
  filter(grepl("nonsynonymous|stopgain|stoploss|frameshift", Type, ignore.case = TRUE))

# gene별 mutation 수 계산
gene_count <- wes_sig %>%
  group_by(`Official Symbol`) %>%
  summarise(mut_count = n()) %>%
  arrange(desc(mut_count))

# 예: mutation 상위 20개
top_genes <- head(gene_count, 20)
top_genes

```

wes_sig <- wes %>%
  filter(grepl("nonsynonymous|stopgain|stoploss|frameshift", Type, ignore.case = TRUE))


```{r}
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
# =====================
# 1. Mutations per megabase per patient (mutation load)
# =====================
patient_load <- wes_sig %>%
  group_by(Patient) %>%
  summarise(mutations = n())

ggplot(patient_load, aes(x = Patient, y = mutations)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Mutations per patient", y = "Mutation count", x = "Patient") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

```

```{r}
# =====================
# 2. Mutation types of significantly mutated genes per patient
# =====================

library(data.table)

# data.table로 변환
wes_dt <- as.data.table(wes_sig)

# mutation_flag 0/1 matrix 생성
mutation_matrix <- dcast(
  wes_dt[, .(mutation_flag = 1), by = .(`Official Symbol`, Patient)],
  `Official Symbol` ~ Patient,
  value.var = "mutation_flag",
  fill = 0
)

# rownames 설정 (gene 이름)
mat <- as.matrix(mutation_matrix[, -1, with = FALSE])
rownames(mat) <- mutation_matrix$`Official Symbol`

# heatmap 그리기
library(pheatmap)
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = c("white", "red"),
         main = "Mutation type per patient")

# top 20 gene heatmap
top_genes <- rowSums(mat) %>% sort(decreasing = TRUE) %>% head(30) %>% names()
mat_top <- mat[top_genes, ]
pheatmap(mat_top, cluster_rows=TRUE, cluster_cols=TRUE,
         color=c("white", "red"), main="Top 30 mutated genes")

#mute 빈도 heatmap
freq_genes <- rowSums(mat > 0)/ncol(mat)
mat_filtered <- mat[freq_genes > 0.05, ]  # 5% 이상 mutation gene만
pheatmap(mat_filtered, cluster_rows=TRUE, cluster_cols=TRUE,
         color=c("white", "red"), main=">0.05freq genes")

```

```{r}
# =====================
# 3. Mutation frequency per gene across patients
# =====================
gene_freq_top50 <- gene_freq %>%
  slice_max(freq, n = 50)  # 상위 50개

ggplot(gene_freq_top50, aes(x = reorder(`Official Symbol`, -freq), y = freq)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  theme_minimal() +
  labs(title = "Top 50 mutated genes", y = "Number of patients", x = "Gene") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

```



```{r}
library(dplyr)

# 관심 유전자 벡터
genes_of_interest <- c("CDH1", "ARID1A", "RHOA")

# 각 유전자별 mutation 환자 번호
patients_per_gene <- wes_sig %>%
  filter(`Official Symbol` %in% genes_of_interest) %>%
  group_by(`Official Symbol`) %>%
  summarise(Patients = list(unique(Patient)))

# 결과 확인
patients_per_gene

```


# defining significantly muatated gene 7 "sig_"
```{r}
library(dplyr)
library(dndscv)

# 1. 데이터 불러오기
wes <- read.table("wes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2. dndscv용 포맷 만들기
# 필요한 컬럼: sampleID, chr, pos, ref, mut, gene
sig_muts <- wes %>%
  transmute(
    sampleID = Patient,
    chr = gsub("chr", "", Chr),   # chr6 → 6
    pos = Locus,                  # 단일염기 위치
    ref = From,
    mut = To,
    gene = `Official.Symbol`
  )

# 3. dNdScv 실행
sig_res <- dndscv(sig_muts)

# 4. 결과에서 significant gene 추출 (q < 0.1)
siggenes <- sig_res$sel_cv %>%
  filter(qallsubs_cv < 0.1)
print(siggenes)

# 5. 결과 txt 저장
write.table(siggenes,
            file = "fig1a_significant_genes_q0.1.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

```


# significantly mutated genes list to plot - bar plot
```{r}
library(ggplot2)
library(dplyr)

# mutation count 계산 (시각화용)
siggenes <- siggenes %>%
  mutate(mutation_count = n_mis + n_non + n_spl)

# 상위 20개만 선택 (dNdScv 결과 순서 유지)
top_genes <- siggenes %>% slice_head(n = 20)

# factor 수준을 dNdScv 순서 그대로 설정
top_genes$gene_name <- factor(top_genes$gene_name, levels = top_genes$gene_name)

# bar plot
ggplot(top_genes, aes(x = gene_name, y = mutation_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal(base_size = 10) +
  labs(x = "Gene", y = "Mutation count",
       title = "Top 20 Significantly Mutated Genes (q < 0.1, dNdScv order)") +
  coord_flip()

```


### 아래는 실제로 result에 포함되지는 않았으나 데이터 cut off 범위 설정 용 top 40 도 확인하였음.
# translation effect - synonymous / non-synonymous
## siggenes top (20)
## bar size proportional to mutation frequency of each gene across patients 
```{r}
library(dplyr)
library(ggplot2)

# 상위 40개 유전자 목록
top_genes_list <- top_genes$gene_name

# 데이터 필터링: 상위 유전자만
wes_top <- wes %>%
  filter(`Official.Symbol` %in% top_genes_list)

# Type 간 분류
wes_top <- wes_top %>%
  mutate(effect = case_when(
    grepl("synonymous", Type, ignore.case = TRUE) ~ "synonymous",
    grepl("nonsynonymous|stopgain|frameshift|nonsense", Type, ignore.case = TRUE) ~ "non-synonymous",
    TRUE ~ "other"
  ))

# 환자 수 확인
n_patients <- length(unique(wes$Patient))

# 유전자별 mutation frequency 계산
effect_summary <- wes_top %>%
  group_by(`Official.Symbol`, effect) %>%
  summarise(muts = n(), .groups = "drop") %>%
  mutate(mutation_frequency = muts / n_patients)

# stacked bar plot
ggplot(effect_summary, aes(x = `Official.Symbol`, y = mutation_frequency, fill = effect)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal(base_size = 10) +
  labs(x = "Gene", y = "Mutation frequency across patients",
       fill = "Translation Effect",
       title = "Synonymous & Non-synonymous Mutation Frequency (Top Genes)")

```


# mutation frquency across patients - about mut type!
```{r}
library(dplyr)
library(ggplot2)

# 상위 40개 유전자 목록
top_genes_list <- top_genes$gene_name

# 상위 유전자 데이터만 필터링
wes_top <- wes %>%
  filter(`Official.Symbol` %in% top_genes_list)

# mutation type 분류
wes_top <- wes_top %>%
  mutate(mutation_type = case_when(
    grepl("stopgain|nonsense", Type, ignore.case = TRUE) ~ "nonsense",
    grepl("frameshift", Type, ignore.case = TRUE) ~ "frameshift",
    grepl("inframe", Type, ignore.case = TRUE) ~ "inframe indel",
    grepl("splice", Type, ignore.case = TRUE) ~ "splice site",
    grepl("nonsynonymous|missense", Type, ignore.case = TRUE) ~ "missense",
    grepl("synonymous", Type, ignore.case = TRUE) ~ "synonymous",
    TRUE ~ "other"
  ))

# 환자 수
n_patients <- length(unique(wes$Patient))

# 유전자별 mutation type frequency 계산
mut_summary <- wes_top %>%
  group_by(`Official.Symbol`, mutation_type) %>%
  summarise(muts = n(), .groups = "drop") %>%
  mutate(frequency = muts / n_patients)

# stacked bar plot
ggplot(mut_summary, aes(x = `Official.Symbol`, y = frequency, fill = mutation_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal(base_size = 10) +
  labs(
    x = "Gene",
    y = "Mutation frequency across patients",
    fill = "Mutation Type",
    title = "Mutation Types in Top Genes"
  ) +
  scale_fill_brewer(palette = "Set2")  # 색깔 조합 깔끔하게

```







