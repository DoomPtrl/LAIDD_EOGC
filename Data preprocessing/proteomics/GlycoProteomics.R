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

GlycoP <- read_tsv('./txt/glycopeptide_change.txt')

# Norm
GlycoP <- as.data.frame(GlycoP)
GlycoP <- GlycoP[rowSums(is.na(GlycoP)) < 40, ]
rownames(GlycoP) <- GlycoP[, 1]
GlycoP[, 3:82] <- apply(GlycoP[, 3:82], 2, as.numeric)
GlycoP[, 3:82] <- normalizeBetweenArrays(GlycoP[, 3:82], method="quantile")

# MAD
Mads <- apply(GlycoP[, 3:82], 1, mad, na.rm = T)
Th <- quantile(Mads, 0.8, na.rm = T)
Top <- which(Mads >= Th)
GlycoP_Mads <- GlycoP[, 3:82][Top, ]
GlycoP_Mads[GlycoP_Mads == 'NaN'] <- 0

# NMF
GlycoP_NMF <- as.matrix(GlycoP_Mads)
GlycoP_NMF_pos <- GlycoP_NMF
GlycoP_NMF_pos[GlycoP_NMF_pos < 0] <- 0
GlycoP_NMF_neg <- - GlycoP_NMF
GlycoP_NMF_neg[GlycoP_NMF_neg < 0] <- 0
GlycoP_NMF <- rbind(GlycoP_NMF_pos, GlycoP_NMF_neg)


# cNMF cluster
k_values <- 2:6
results_list <- list()  # 전체 결과 저장용

num_cores <- 10
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results_list <- foreach(k = k_values, .packages = "CancerSubtypes") %dopar% {
  ExecuteCNMF(GlycoP_NMF, clusterNum = k, nrun = 200)
}

# names 할당
names(results_list) <- paste0("k_", k_values)

stopCluster(cl)


for (name in names(results_list)) {
  dist_mat <- results_list[[name]]$distanceMatrix
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
  pca_res <- prcomp(t(GlycoP_NMF), scale. = TRUE)
  
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


k_values <- as.numeric(gsub("k_", "", names(results_list)))  # k값 추출
cophenetic_values <- numeric(length(k_values))

for (i in seq_along(k_values)) {
  result <- results_list[[paste0("k_", k_values[i])]]
  dist_mat <- as.dist(result$distanceMatrix)
  hc <- hclust(dist_mat)
  cophenetic_mat <- cophenetic(hc)
  cophenetic_values[i] <- cor(dist_mat, cophenetic_mat)
}

df <- data.frame(k = k_values, CopheneticCorrelation = cophenetic_values)

ggplot(df, aes(x = k, y = CopheneticCorrelation)) +
  geom_line() +
  geom_point(color = "darkred", size = 3) +
  labs(title = "Cophenetic Correlation by Number of Clusters",
       x = "Number of Clusters (k)", y = "Cophenetic Correlation") +
  theme_minimal()




saveRDS(results_list, './rds/GlycoP.rds')
results_list <- readRDS('./rds/GlycoP.rds')








group <- factor(results_list[['k_2']]$group)
group <- paste0("Group", group)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(factor(group))

fit <- lmFit(GlycoP[, 3:82], design)
contrast <- makeContrasts(Group1vsGroup2 = Group1 - Group2, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "BH", number = Inf)
head(results)

sig_GlycoP <- rownames(results)[results$adj.P.Val < 0.00001 & abs(results$logFC) > 0.5]



# signature 시각화
GlycoP_Sel <- GlycoP %>%
  filter(.$Peptide %in% sig_GlycoP)

Grp <- factor(results_list[['k_2']]$group, levels = c(1,2,3,4,5), labels = c("Cl1", "Cl2", "Cl3", "Cl4", "Cl5"))
group_colors <- c("Cl1" = "#E41A1C", "Cl2" = "#377EB8", "Cl3" = "#4DAF4A", "Cl4" = "#984EA3", "Cl5" = "#f4d505")
GlycoP_Sel[is.na(GlycoP_Sel)] <- 0

Heatmap(as.matrix(GlycoP_Sel[, 3:82]),
        name = "Glyco",
        col = colorRamp2(c(-1, 0, 1), c("blue2", "gray10", "red2")),
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
        column_title = "Hierarchical Clustering Heatmap of Glyco"
)
