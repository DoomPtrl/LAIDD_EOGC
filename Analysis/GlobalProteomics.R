library(ConsensusClusterPlus)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(CancerSubtypes)
library(cluster)
library(NMF)
library(doParallel)
library(foreach)
set.seed(17)

GlobalP <- read_tsv('./txt/protein_expression_change.txt')

# Norm
GlobalP <- as.data.frame(GlobalP)
GlobalP <- GlobalP[rowSums(is.na(GlobalP)) < 40, ]
rownames(GlobalP) <- GlobalP[, 2]
GlobalP[, 7:86] <- apply(GlobalP[, 7:86], 2, as.numeric)
GlobalP[, 7:86] <- normalizeBetweenArrays(GlobalP[, 7:86], method="quantile")

# MAD
Mads <- apply(GlobalP[, 7:86], 1, mad, na.rm = T)
Th <- quantile(Mads, 0.9, na.rm = T)
Top <- which(Mads >= Th)
GlobalP_Mads <- GlobalP[, 7:86][Top, ]
GlobalP_Mads[is.na(GlobalP_Mads)] <- 0

# NMF
GlobalP_NMF <- as.matrix(GlobalP_Mads)
GlobalP_NMF_pos <- GlobalP_NMF
GlobalP_NMF_pos[GlobalP_NMF_pos < 0] <- 0
GlobalP_NMF_neg <- - GlobalP_NMF
GlobalP_NMF_neg[GlobalP_NMF_neg < 0] <- 0
GlobalP_NMF <- rbind(GlobalP_NMF_pos, GlobalP_NMF_neg)


# cNMF cluster
k_values <- 2:6
results_list <- list()  # 전체 결과 저장용

num_cores <- 6
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results_list <- foreach(k = k_values, .packages = "CancerSubtypes") %dopar% {
  ExecuteCNMF(GlobalP_NMF, clusterNum = k, nrun = 200)
}

# names 할당
names(results_list) <- paste0("k_", k_values)

stopCluster(cl)


for (name in names(results_list)) {
  dist_mat <- results_list[[name]]$distanceMatrix
  # pheatmap은 기본적으로 바로 출력하지만, knit에서는 print()를 함께 써야 확실
  print(
    pheatmap(
      dist_mat,
      main = paste0("Distance Matrix Heatmap (", name, ")"),
      border_color = NA,
      labels_row = rep('', nrow(dist_mat)),
      labels_col = rep('', ncol(dist_mat))
    )
  )
}


k_values <- as.numeric(gsub("k_", "", names(results_list)))

for (k in k_values) {
  result <- results_list[[paste0("k_", k)]]
  
  # PCA 수행 (scale. = TRUE 권장)
  pca_res <- prcomp(t(GlobalP_NMF), scale. = TRUE)
  
  # PC1, PC2 좌표와 클러스터 할당 결합
  pca_df <- data.frame(PC1 = pca_res$x[,1],
                       PC2 = pca_res$x[,2],
                       cluster = factor(result$group))
  
  # 시각화
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(size = 4, alpha = 0.8) +
    theme_minimal() +
    labs(title = paste0("PCA Plot (k = ", k, ")"),
         x = "Principal Component 1",
         y = "Principal Component 2")
  
  print(p)
}



saveRDS(results_list, './rds/GlobalP.rds')
results_list <- readRDS('./rds/GlobalP.rds')



group <- factor(results_list[['k_3']]$group)
group <- paste0("Group", group)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(factor(group))

fit <- lmFit(GlobalP[, 7:86], design)
fit2 <- eBayes(fit)

top_table_all <- topTable(fit2, adjust = "BH", number = Inf)
head(top_table_all)

top_table_all2 <- top_table_all %>%
  filter(adj.P.Val < 0.00001)
sig_P <- rownames(top_table_all2)[!(top_table_all2$Group1 > 0 & top_table_all2$Group2 > 0 & top_table_all2$Group3 > 0 | top_table_all2$Group1 < 0 & top_table_all2$Group2 < 0 & top_table_all2$Group3 < 0)]



# signature 시각화
GlobalP_Sel <- GlobalP %>%
  filter(.$Prot_ID %in% sig_P)

Grp <- factor(results_list[['k_3']]$group, levels = c(1,2,3,4,5), labels = c("Cl1", "Cl2", "Cl3", "Cl4", "Cl5"))
group_colors <- c("Cl1" = "#E41A1C", "Cl2" = "#377EB8", "Cl3" = "#4DAF4A", "Cl4" = "#984EA3", "Cl5" = "#f4d505")
GlobalP_Sel[is.na(GlobalP_Sel)] <- 0

Heatmap(as.matrix(GlobalP_Sel[, 7:86]),
        name = "Global protein",
        col = colorRamp2(c(-0.5, 0, 0.5), c("blue2", "gray10", "red2")),
        show_row_names = FALSE,
        show_column_names = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        show_column_dend = F,
        cluster_columns = F,
        cluster_rows = T,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete",
        show_row_dend = F,
        column_split  = Grp,                 # 그룹별 행 분할
        row_title = NULL,               # 기본 그룹명 대신 없앰
        row_title_side = "left",
        row_gap = unit(5, "mm"),
        top_annotation = HeatmapAnnotation(Group = Grp, col = list(Group = group_colors)),
        column_title = "Hierarchical Clustering Heatmap of RNA"
)
