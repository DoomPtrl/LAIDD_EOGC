.libPaths()
# R, 터미널에서 R 들어가기
setwd("/home/Data_Drive_8TB/tyws0213/0507/") 
getwd()
source('./renv/activate.R')
# Test 폴더 내에 renv 파일이 생성되었는지 확인
.libPaths()
library(dplyr)
library(tidyr)
library(preprocessCore)
library(parallel)
library(pheatmap)

# ---------------------------
# 1. 데이터 불러오기
# ---------------------------
phospho <- read.delim(
  "~/LAIDD/phosphopeptide_change.txt",
  header=TRUE, sep="\t", check.names=FALSE
)

mut <- read.delim(
  "~/LAIDD/SM_nonsyn_80_3level.txt",
  header=TRUE, sep="\t", stringsAsFactors=FALSE
)

# ---------------------------
# 2. phospho 전처리
# ---------------------------
gene_symbols <- phospho$Symbol
peptide_ids  <- make.unique(phospho[,1])   # peptide ID
phospho_mat <- phospho[, -(1:2)]
rownames(phospho_mat) <- peptide_ids
phospho_mat <- apply(phospho_mat, 2, as.numeric)
rownames(phospho_mat) <- peptide_ids
phospho_mat[is.nan(phospho_mat)] <- NA

# ---------------------------
# 3. 단일 환자 케이스만 선택 (T = N+1)
# ---------------------------
valid_cols <- sapply(colnames(phospho_mat), function(nm) {
  N_id <- as.numeric(sub("N","", stringr::str_extract(nm, "N[0-9]+")))
  T_id <- as.numeric(sub("T","", stringr::str_extract(nm, "T[0-9]+")))
  !is.na(N_id) && !is.na(T_id) && (T_id == N_id + 1)
})
phospho_single <- phospho_mat[, valid_cols, drop=FALSE]

# ---------------------------
# 4. 열 확장 (rep1/rep2 Normal, rep1/rep2 Tumor)
# ---------------------------
expand_single <- function(mat) {
  new_list <- list()
  for(nm in colnames(mat)) {
    values <- mat[, nm]
    N_id <- as.numeric(sub("N","", stringr::str_extract(nm, "N[0-9]+")))
    T_id <- as.numeric(sub("T","", stringr::str_extract(nm, "T[0-9]+")))
    nlab <- paste0(N_id,"N")
    tlab <- paste0(T_id,"T")
    new_list[[paste0(nlab,"_rep1")]] <- values
    new_list[[paste0(nlab,"_rep2")]] <- values
    new_list[[paste0(tlab,"_rep1")]] <- values
    new_list[[paste0(tlab,"_rep2")]] <- values
  }
  expanded <- as.data.frame(new_list)
  rownames(expanded) <- rownames(mat)
  return(expanded)
}
phospho_expanded <- expand_single(phospho_single)
colnames(phospho_expanded) <- sub("_rep[12]","", colnames(phospho_expanded))
colnames(phospho_expanded) <- sub("^X", "", colnames(phospho_expanded))

# ---------------------------
# 5. mutation 데이터 wide-format 변환
# ---------------------------
mut_status <- mut %>%
  dplyr::select(Patient, Official.Symbol) %>%
  dplyr::distinct() %>%
  dplyr::mutate(MutStatus="Mut") %>%
  tidyr::pivot_wider(
    id_cols=Patient,
    names_from=Official.Symbol,
    values_from=MutStatus,
    values_fill="WT"
  ) %>%
  as.data.frame()
rownames(mut_status) <- mut_status$Patient
mut_status <- mut_status[, -1]

# ---------------------------
# 6. 샘플 매칭
# ---------------------------
common_samples <- intersect(colnames(phospho_expanded), rownames(mut_status))
phospho_expanded <- phospho_expanded[, common_samples, drop=FALSE]
mut_status <- mut_status[common_samples, , drop=FALSE]

# ---------------------------
# 7. filtering + normalization
# ---------------------------
min_detected <- ceiling(ncol(phospho_expanded) * 0.5)
valid_rows <- rowSums(!is.na(phospho_expanded)) >= min_detected
phospho_filt <- phospho_expanded[valid_rows, ]

phospho_norm <- preprocessCore::normalize.quantiles(as.matrix(phospho_filt))
rownames(phospho_norm) <- rownames(phospho_filt)
colnames(phospho_norm) <- colnames(phospho_filt)

########################## 유의미 gene 뽑기 ############################

test_mut_wt_lm <- function(gene, expr, mut_expr) {
  group <- mut_expr[, gene]
  
  # Mut, WT 샘플 수 확인
  n_mut <- sum(group == "Mut", na.rm = TRUE)
  n_wt  <- sum(group == "WT", na.rm = TRUE)
  
  # Mut이 최소 3개 이상, WT도 최소 3개 이상이어야
  if (n_mut < 3 | n_wt < 3) return(NULL)
  
  design <- model.matrix(~0 + factor(group))
  colnames(design) <- c("Mut","WT")
  
  fit <- lmFit(expr, design)
  cont <- makeContrasts(Mut - WT, levels=design)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)
  
  res <- topTable(fit2, number=Inf, adjust.method="BH")
  res$Gene <- gene
  res$PhosphoID <- rownames(res)
  res
}

# 필터링된 유전자만 병렬 실행
library(limma)
library(BiocParallel)

# 1. valid_genes 재정의
valid_genes <- colnames(mut_status)[
  colSums(mut_status == "Mut", na.rm = TRUE) >= 3 &
    colSums(mut_status == "WT", na.rm = TRUE) >= 3
]

# 2. 샘플 순서 확인
all(colnames(phospho_norm) == rownames(mut_status))  # TRUE여야 함

# 3. 병렬 실행 (50코어)
param <- MulticoreParam(workers = 50)
results_list <- bplapply(
  valid_genes,
  FUN = function(g) test_mut_wt_lm(g, phospho_norm, mut_status),
  BPPARAM = param
)

all_results <- bind_rows(results_list)

gene_summary <- all_results %>%
  group_by(Gene) %>%
  dplyr::summarise(
    min_p = min(P.Value, na.rm=TRUE),
    min_FDR = min(adj.P.Val, na.rm=TRUE),
    n_sig = sum(adj.P.Val < 0.05, na.rm=TRUE)   # 유의한 phosphopeptide 개수
  ) %>%
  arrange(min_p)

head(gene_summary, 20)

# 패키지 설치
install.packages("writexl")

# 패키지 로드
library(writexl)

# 엑셀로 저장
write_xlsx(gene_summary, path = "~/LAIDD/gene_summary.xlsx")


library(pheatmap)
library(RColorBrewer)

plot_mut_cistrans_heatmap <- function(
    gene,
    phospho_norm,
    mut_status,
    p_cutoff = 0.05,
    diff_cutoff = 0.5
) {
  if (!(gene %in% colnames(mut_status))) {
    message("Gene not in mutation data: ", gene)
    return(NULL)
  }
  
  mut_samples <- rownames(mut_status[mut_status[[gene]] == "Mut", , drop=FALSE])
  wt_samples  <- rownames(mut_status[mut_status[[gene]] == "WT", , drop=FALSE])
  
  common_samples <- intersect(colnames(phospho_norm), c(mut_samples, wt_samples))
  if (length(common_samples) == 0) {
    message("No overlapping samples for ", gene)
    return(NULL)
  }
  
  sub_mat <- phospho_norm[, common_samples, drop=FALSE]
  
  res <- apply(sub_mat, 1, function(x) {
    x_mut <- x[colnames(sub_mat) %in% mut_samples]
    x_wt  <- x[colnames(sub_mat) %in% wt_samples]
    
    if (sum(!is.na(x_mut)) < 2 || sum(!is.na(x_wt)) < 2) {
      return(c(pval=NA, diff=NA, mean_mut=NA, mean_wt=NA))
    }
    
    suppressWarnings({
      p <- wilcox.test(x_mut, x_wt)$p.value
      diff <- median(x_mut, na.rm=TRUE) - median(x_wt, na.rm=TRUE)
      mean_mut <- mean(x_mut, na.rm=TRUE)
      mean_wt  <- mean(x_wt, na.rm=TRUE)
    })
    return(c(pval = p, diff = diff, mean_mut = mean_mut, mean_wt = mean_wt))
  })
  
  res <- t(res)
  res <- as.data.frame(res)
  res <- res[complete.cases(res), ]
  
  sig_res <- res[
    res$pval < p_cutoff &
      res$diff > diff_cutoff &
      res$mean_mut > 0 & res$mean_wt < 0, 
  ]
  
  sig_peps <- rownames(sig_res)
  heat_mat <- sub_mat[sig_peps, , drop=FALSE]
  
  ordered_samples <- c(intersect(mut_samples, colnames(heat_mat)),
                       intersect(wt_samples, colnames(heat_mat)))
  heat_mat <- heat_mat[, ordered_samples, drop=FALSE]
  
  ann_col <- data.frame(
    Mutation = ifelse(colnames(heat_mat) %in% mut_samples, "Mut", "WT")
  )
  rownames(ann_col) <- colnames(heat_mat)
  
  library(pheatmap)
  
  # heatmap 그릴 때만 색상 조정
  my_palette <- colorRampPalette(c("limegreen", "black", "red3"))(50)
  
  # 🔹 강제적으로 색상 대비를 주기 위해 breaks 설정
  # (z-score -2 ~ +2를 중점적으로 색상에 배분)
  my_breaks <- seq(-2, 2, length.out = 51)
  
  pheatmap::pheatmap(
    heat_mat,
    scale="row",
    color=my_palette,
    breaks=my_breaks,
    annotation_col = ann_col,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = FALSE,
    na_col = "black",
    border_color=NA,
    main = paste0(gene, " mutation: cis+trans phosphopeptides (p<", p_cutoff, ")")
  )
}

options(bitmapType = 'cairo')
plot_mut_cistrans_heatmap("ARID1A", phospho_norm, mut_status)
peplot_mut_cistrans_heatmap("MUC5B", phospho_norm, mut_status)
plot_mut_cistrans_heatmap("TP53", phospho_norm, mut_status)
plot_mut_cistrans_heatmap("TSHZ1", phospho_norm, mut_status)

get_sig_peptides <- function(
    gene,
    phospho_norm,
    mut_status,
    p_cutoff = 0.05,
    diff_cutoff = 0.5
) {
  if (!(gene %in% colnames(mut_status))) {
    message("Gene not in mutation data: ", gene)
    return(NULL)
  }
  
  # Mut / WT 샘플 나누기
  mut_samples <- rownames(mut_status[mut_status[[gene]] == "Mut", , drop=FALSE])
  wt_samples  <- rownames(mut_status[mut_status[[gene]] == "WT", , drop=FALSE])
  
  common_samples <- intersect(colnames(phospho_norm), c(mut_samples, wt_samples))
  if (length(common_samples) == 0) {
    message("No overlapping samples for ", gene)
    return(NULL)
  }
  
  sub_mat <- phospho_norm[, common_samples, drop=FALSE]
  
  # 각 peptide 별 Wilcoxon test
  res <- apply(sub_mat, 1, function(x) {
    x_mut <- x[colnames(sub_mat) %in% mut_samples]
    x_wt  <- x[colnames(sub_mat) %in% wt_samples]
    
    if (sum(!is.na(x_mut)) < 2 || sum(!is.na(x_wt)) < 2) {
      return(c(pval=NA, diff=NA, mean_mut=NA, mean_wt=NA))
    }
    
    suppressWarnings({
      p <- wilcox.test(x_mut, x_wt)$p.value
      diff <- median(x_mut, na.rm=TRUE) - median(x_wt, na.rm=TRUE)
      mean_mut <- mean(x_mut, na.rm=TRUE)
      mean_wt  <- mean(x_wt, na.rm=TRUE)
    })
    return(c(pval = p, diff = diff, mean_mut = mean_mut, mean_wt = mean_wt))
  })
  
  res <- t(res)
  res <- as.data.frame(res)
  res <- res[complete.cases(res), ]
  
  if (nrow(res) == 0) {
    message("No valid peptides for ", gene)
    return(NULL)
  }
  
  # 🔹 조건: p < cutoff, diff > cutoff, mean_mut > 0, mean_wt < 0
  sig_res <- res[
    res$pval < p_cutoff &
      res$diff > diff_cutoff &
      res$mean_mut > 0 &
      res$mean_wt < 0, 
  ]
  
  if (nrow(sig_res) == 0) {
    message("No significant peptides for ", gene, " under thresholds.")
    return(NULL)
  }
  
  sig_peps <- rownames(sig_res)
  return(sig_peps)
}

# 예시 실행
sig_peptides_ARID1A <- get_sig_peptides(
  gene = "ARID1A",
  phospho_norm = phospho_norm,
  mut_status = mut_status,
  p_cutoff = 0.05,
  diff_cutoff = 0.5
)
sig_peptides_CDH1 <- get_sig_peptides(
  gene = "CDH1",
  phospho_norm = phospho_norm,
  mut_status = mut_status,
  p_cutoff = 0.05,
  diff_cutoff = 0.5
)
sig_peptides_RHOA <- get_sig_peptides(
  gene = "RHOA",
  phospho_norm = phospho_norm,
  mut_status = mut_status,
  p_cutoff = 0.05,
  diff_cutoff = 0.5
)
sig_peptides_TP53

######################################## Gene Set 분석 #########################

# peptide → gene 매핑 테이블
peptide_to_gene <- data.frame(
  peptide_id = make.unique(phospho[,1]),  # 첫 번째 열이 peptide ID
  gene_symbol = phospho$Symbol
)

library(dplyr)

# phospho_norm → data.frame 변환
phospho_df <- as.data.frame(phospho_norm)
phospho_df$peptide_id <- rownames(phospho_norm)

# gene symbol 붙이기
phospho_df <- left_join(phospho_df, peptide_to_gene, by="peptide_id")

# gene 단위 평균/합계로 집계 (여기서는 평균 예시)
phospho_gene <- phospho_df %>%
  group_by(gene_symbol) %>%
  dplyr::summarise(across(where(is.numeric), mean, na.rm=TRUE)) %>%
  as.data.frame()

rownames(phospho_gene) <- phospho_gene$gene_symbol
phospho_gene <- phospho_gene[,-1]


phospho_gene

library(dplyr)

# 1. peptide ID 컬럼 만들기
phospho <- phospho %>%
  mutate(PeptideID = make.unique(phospho[,1]))

# 2. PeptideID와 Symbol만 선택
peptide_to_gene <- phospho %>%
  dplyr::select(PeptideID, Symbol)

# 3. sig_peptides와 mapping
sig_peptides_gene_ARID1A <- peptide_to_gene %>%
  dplyr::filter(PeptideID %in% sig_peptides_ARID1A)
sig_peptides_gene_CDH1 <- peptide_to_gene %>%
  dplyr::filter(PeptideID %in% sig_peptides_CDH1)
sig_peptides_gene_RHOA <- peptide_to_gene %>%
  dplyr::filter(PeptideID %in% sig_peptides_RHOA)

sig_peptides_gene_TP53

# ---------------------------
# 0. 필요한 패키지 로드
# ---------------------------
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)   # 인간 유전자 annotation
library(dplyr)

# ---------------------------
# 1. sig_peptides에서 gene 목록 추출
# ---------------------------
sig_genes_ARID1A <- unique(sig_peptides_gene_ARID1A$Symbol)
sig_genes_CDH1 <- unique(sig_peptides_gene_CDH1$Symbol)
sig_genes_RHOA <- unique(sig_peptides_gene_RHOA$Symbol)

sig_genes <- sig_genes[!is.na(sig_genes)]  # NA 제거

# ---------------------------
# 2. gene symbol → Entrez ID 변환
# ---------------------------
gene_entrez_ARID1A <- bitr(sig_genes_ARID1A, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
gene_entrez_CDH1 <- bitr(sig_genes_CDH1, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
gene_entrez_RHOA <- bitr(sig_genes_RHOA, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
# ---------------------------
# 3. ORA 수행 (KEGG 예시)
# ---------------------------
ora_kegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",       # 인간
  pvalueCutoff = 0.05
)

# ---------------------------
# 4. ORA 수행 (GO 예시)
# ---------------------------
ora_go_ARID1A <- enrichGO(
  gene = gene_entrez_ARID1A$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",             # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE         # EntrezID → gene symbol
)
ora_go_CDH1 <- enrichGO(
  gene = gene_entrez_CDH1$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",             # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE         # EntrezID → gene symbol
)
ora_go_RHOA <- enrichGO(
  gene = gene_entrez_RHOA$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",             # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE         # EntrezID → gene symbol
)

# ---------------------------
# 5. 결과 확인
# ---------------------------
head(ora_kegg)
head(ora_go)

# ---------------------------
# 6. 시각화 (Dotplot 예시)
# ---------------------------
library(enrichplot)

dotplot(ora_kegg, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")
dotplot(ora_go_CDH1, showCategory = 20) + ggtitle("KRAS GO Biological Process Enrichment")

install.packages("ggplot2")
library(ggplot2)
install.packages("RCPA")
library(RCPA)

install.packages("CePa")
devtools::install_github("jokergoo/igraph")

# 패키지 로드
library(CePa)
library(igraph)

# ---------------------------
# 0. 필요한 패키지 로드
# ---------------------------
library(CePa)
library(ggplot2)

# ---------------------------
# 1. CPDB 경로 데이터 불러오기 (미리 다운로드한 .gmt 파일 필요)
# ---------------------------
# CPDB 웹사이트에서 "Pathways → Export" → GMT 포맷 파일 다운로드

# ---------------------------
# 1. 파일 불러오기
# ---------------------------
cpdb <- read.delim(
  "~/LAIDD/CPDB_pathways_genes.tab",
  header = FALSE,
  stringsAsFactors = FALSE
)

# cpdb[,1] = Pathway 이름
# cpdb[,2] = Gene 목록 (쉼표로 구분되어 있다고 가정)

# ---------------------------
# 2. ORA용 list 변환
# ---------------------------
cpdb <- cpdb[-1, ]

cpdb_list <- setNames(
  lapply(cpdb[,4], function(x) strsplit(x, ",")[[1]]),
  cpdb[,1]
)

# cpdb_list : list 형태 (이름 = 경로, 값 = gene 벡터)

cpdb_df <- do.call(rbind, lapply(names(cpdb_list), function(pw) {
  data.frame(term = pw, gene = cpdb_list[[pw]], stringsAsFactors = FALSE)
}))

# 확인
head(cpdb_df)


# ---------------------------
# 3. 확인
# ---------------------------
str(cpdb_list[[1]])  # 첫 번째 경로 확인
head(names(cpdb_list))  # 경로 이름 확인


# ---------------------------
# 0. 필요한 패키지 로드
# ---------------------------
library(clusterProfiler)
library(dplyr)

# ---------------------------
# 1. sig_peptides_gene에서 유전자 목록 추출
# ---------------------------
sig_genes <- unique(sig_peptides_gene$Symbol)
sig_genes <- sig_genes[!is.na(sig_genes)]  # NA 제거

# ---------------------------
# 2. ORA 수행 (ConsensusPathDB 경로)
# ---------------------------
ora_cpdb_ARID1A <- enricher(
  gene         = sig_genes_ARID1A,
  TERM2GENE    = cpdb_df,   # cpdb_list는 경로별 유전자 리스트
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
ora_cpdb_CDH1 <- enricher(
  gene         = sig_genes_CDH1,
  TERM2GENE    = cpdb_df,   # cpdb_list는 경로별 유전자 리스트
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
ora_cpdb_RHOA <- enricher(
  gene         = sig_genes_RHOA,
  TERM2GENE    = cpdb_df,   # cpdb_list는 경로별 유전자 리스트
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# ---------------------------
# 3. 결과 확인
# ---------------------------
head(ora_cpdb)

# ---------------------------
# 4. 시각화 (Dotplot)
# ---------------------------
library(enrichplot)
dotplot(ora_cpdb_CDH1, showCategory = 20) + ggtitle("TTN ORA: CPDB Pathways")
dotplot(ora_cpdb_RHOA, showCategory = 20) + ggtitle("MUC5B ORA: CPDB Pathways")
dotplot(ora_cpdb_ARID1A, showCategory = 20) + ggtitle("TP53 ORA: CPDB Pathways")



library(clusterProfiler)
library(dplyr)
library(reshape2)
library(pheatmap)

# 1. 각 ORA 결과 데이터프레임화
df_TSHZ1 <- as.data.frame(ora_cpdb_TSHZ1) %>% mutate(gene_set="TSHZ1")
df_TLR4  <- as.data.frame(ora_cpdb_TLR4)  %>% mutate(gene_set="TLR4")
df_TTN  <- as.data.frame(ora_cpdb_TTN)  %>% mutate(gene_set="TTN")

# 2. 필요한 컬럼만 선택
df_all <- bind_rows(df_TSHZ1, df_TLR4, df_NRG1) %>%
  dplyr::select(gene_set, ID, Description, p.adjust)

# 3. heatmap용 데이터 준비 (-log10 p.adjust, NA는 0으로 대체 가능)
heat_df <- df_all %>%
  mutate(logp = -log10(p.adjust)) %>%
  dcast(ID + Description ~ gene_set, value.var="logp")

rownames(heat_df) <- heat_df$Description   # pathway 이름을 rowname으로
heat_mat_overlap <- as.matrix(heat_df[ , -(1:2) ]) # ID, Description 제외하고 매트릭스화

library(pheatmap)

# 1. NA는 0처럼 표시되도록 NA를 0으로 치환 (선택)
heat_mat_overlap[is.na(heat_mat_overlap)] <- 0

# 2. 사용자 정의 색상 팔레트
my_colors <- colorRampPalette(c("white", "yellow", "red"))(100)
rownames(heat_mat_overlap)[3] <- "Diseases of signal transduction by growth factor receptors\nand second messengers"
rownames(heat_mat_overlap)[24] <- "Thyroid hormones production and their peripheral downstream\nsignaling effects"
# 3. heatmap 그리기
pheatmap(
  heat_mat_overlap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = my_colors,
  main = "Shared Pathways",
  na_col = "white",
  angle_col = "45",
  fontsize_row = 6,   # 행 이름 글자 크기 조정
  legend = 
) 


BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

# 색상 함수 정의
col_fun <- colorRamp2(
  breaks = c(0, max(heat_mat_overlap, na.rm=TRUE)/2, max(heat_mat_overlap, na.rm=TRUE)),
  colors = c("white", "yellow", "red")
)

ht <- Heatmap(
  heat_mat_overlap,
  column_title = "Shared Pathways",
  column_title_gp = gpar(fontsize = 12),
  name = "-log10(p.adjust)",   # 범례 제목
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 9),
  
  row_names_gp = gpar(fontsize = 6),   # 원하는 크기로 줄이기
  # 🔹 각 셀마다 네모 테두리 추가
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = "gray45", fill = NA, lwd = 0.5))
  },
  # 🔹 Heatmap 테두리 넣기
  heatmap_legend_param = list(
    direction = "horizontal",  # 범례 수평
    title_position = "topcenter",
    legend_width = unit(6, "cm")
    
  )
)

draw(ht, heatmap_legend_side = "bottom")  # ✅ 범례를 위로 이동
