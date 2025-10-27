df_RNA<- read.csv("~/LAIDD/RNAseq_data.csv", row.names = 1)
df_TCGA <- read.csv("~/LAIDD/TCGA_STAD_RNA_count_matrix.csv", row.names = 1)

# 0. 라이브러리 불러오기
BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)
library(DESeq2)
library(ggplot2)

############################################################################
# ✨ STEP 0: 데이터 준비 (이 부분만 실제 데이터 파일로 교체하세요)
############################################################################

df_RNA<- read.csv("~/LAIDD/RNAseq_data.csv", row.names = 1)
df_TCGA <- read.csv("~/LAIDD/TCGA_STAD_RNA_count_matrix.csv", row.names = 1)

############################################################################
# ✨ STEP 1: 유전자 ID 통일 (ENSG -> Gene Symbol)
############################################################################

cat("\n--------- 1. 유전자 ID 변환 시작 ---------\n")

# 1.1. biomaRt로 ENSG ID와 Gene Symbol 매핑 정보 가져오기
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = rownames(df_RNA),
  mart = ensembl
)

# 1.2. 데이터셋1에 매핑 정보 결합 및 정리
dataset1_mapped <- df_RNA %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  inner_join(gene_map, by = "ensembl_gene_id") %>%
  dplyr::select(-ensembl_gene_id) %>% # 기존 ENSG ID 열 제거
  dplyr::filter(hgnc_symbol != "" & !is.na(hgnc_symbol)) # Gene Symbol이 없는 경우 제외

dataset1_mapped <- df_RNA %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_map, by = "ensembl_gene_id") %>%   # inner_join → left_join으로 변경
  dplyr::filter(hgnc_symbol != "" & !is.na(hgnc_symbol)) %>%
  dplyr::group_by(hgnc_symbol) %>%                 # symbol별로 묶기
  dplyr::summarise(across(where(is.numeric), mean)) %>%  # 중복된 유전자 평균 처리
  ungroup() %>%                    
  tibble::column_to_rownames("hgnc_symbol")

cat("✅ 매핑 완료! 남은 유전자 수:", nrow(dataset1_mapped), "\n")



############################################################################
# ✨ STEP 2: 데이터셋 통합
############################################################################

cat("\n--------- 2. 데이터셋 통합 시작 ---------\n")

# 2.1. 두 데이터셋에 공통으로 존재하는 유전자 찾기
common_genes <- intersect(rownames(dataset1_mapped), rownames(df_TCGA))
cat("공통 유전자 개수:", length(common_genes), "\n")

# 2.2. 공통 유전자를 기준으로 데이터 필터링 및 결합
combined_counts <- cbind(
  dataset1_mapped[common_genes, ],
  df_TCGA[common_genes, ]
)

# NA 값이 없는지 확인 (만약 있다면 0으로 대체)
combined_counts[is.na(combined_counts)] <- 0

cat("--------- 통합된 Count Matrix 확인 ---------\n")
print(head(combined_counts))

# 예시 메타데이터 생성
meta <- data.frame(
  sample = colnames(combined_counts),
  batch = ifelse(grepl("^SRR", colnames(combined_counts)), "SRR", "TCGA"),
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$sample

# install.packages("sva")  # 필요시
BiocManager::install("sva")
library(sva)

# meta$batch이 있어야 함
batch <- meta[colnames(combined_counts), "batch"]

# ComBat-Seq 실행 (정수 matrix 사용)
combat_seq_counts <- ComBat_seq(counts = as.matrix(combined_counts),
                                batch = batch,
                                group = NULL)   # group(조건)을 넣어 생물학적 차이 보존 가능

# 결과는 정수 행렬
dim(combat_seq_counts)

# 라이브러리
library(edgeR)
library(ggplot2)
library(limma)

# ComBat-Seq 결과가 "combat_seq_counts"
# meta에 batch, condition 정보가 있다고 가정
# meta의 rownames == colnames(combat_seq_counts) 이어야 함

library(edgeR)

# edgeR DGEList 생성
dge_raw <- DGEList(counts = as.matrix(combined_counts))  # matrix 형태
dge_raw <- calcNormFactors(dge_raw)

# logCPM 변환
logCPM_before <- cpm(dge_raw, log = TRUE, prior.count = 1)

dge_combat <- DGEList(counts = as.matrix(combat_seq_counts))
dge_combat <- calcNormFactors(dge_combat)
logCPM_combat <- cpm(dge_combat, log = TRUE, prior.count = 1)

library(ggplot2)

pca_plot <- function(mat, meta, title, color_by="batch", shape_by=NULL){
  pca <- prcomp(t(mat), scale. = TRUE)
  df <- data.frame(PC1 = pca$x[,1],
                   PC2 = pca$x[,2],
                   batch = meta[colnames(mat), color_by])
  
  if(!is.null(shape_by)){
    df$shape <- meta[colnames(mat), shape_by]
    ggplot(df, aes(x = PC1, y = PC2, color = batch, shape = shape)) +
      geom_point(size = 2.5, alpha = 0.8) +
      theme_bw() + ggtitle(title)
  } else {
    ggplot(df, aes(x = PC1, y = PC2, color = batch)) +
      geom_point(size = 2.5, alpha = 0.8) +
      theme_bw() + ggtitle(title)
  }
}

p1 <- pca_plot(logCPM_before, meta, "Before ComBat-Seq")
p2 <- pca_plot(logCPM_combat, meta, "After ComBat-Seq")

p1
p2

write.table(combat_seq_counts, "~/LAIDD/bulk_counts_CIBERSORTx.txt",
            sep = "\t", quote = FALSE, col.names = NA)

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
singlecell.data <- Read10X(data.dir = "~/LAIDD/10X_Counts/")
# Initialize the Seurat object with the raw (non-normalized data).
singlecell <- CreateSeuratObject(counts = singlecell.data, project = "singlecell", min.cells = 3, min.features = 200)
singlecell

singlecell <- SCTransform(singlecell, verbose = FALSE)

celltype_anno <- read.table("~/LAIDD/major_celltypes.txt", header = TRUE, row.names = 1)
# rownames = cell barcodes

# Seurat object에 metadata로 추가
singlecell$celltype <- celltype_anno[colnames(singlecell), "major_celltype"]

library(dplyr)

# counts → log-normalized
sc_counts <- as.matrix(singlecell@assays$SCT@counts)  # SCTransform 기준
sc_counts <- log1p(sc_counts)  # log1p transform

# cell type별 평균 발현
signature <- sc_counts %>%
  t() %>%  # cells in rows
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  left_join(celltype_anno %>% tibble::rownames_to_column("cell"), by="cell") %>%
  dplyr::group_by(major_celltype) %>%
  dplyr::summarise(across(where(is.numeric), mean)) %>%
  tibble::column_to_rownames("major_celltype") %>%
  t()  # genes in rows, cell types in columns

# 파일로 저장
write.table(signature, "~/LAIDD/sc_signature_CIBERSORTx.txt", sep="\t", quote = FALSE, col.names = NA)

