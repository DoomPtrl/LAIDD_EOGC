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
library(reshape2)
set.seed(17)

PhosphoP <- read_tsv('./txt/Phosphopeptide_change.txt')

# Norm
PhosphoP <- as.data.frame(PhosphoP)
PhosphoP <- PhosphoP[rowSums(is.na(PhosphoP)) < 40, ]
rownames(PhosphoP) <- PhosphoP[, 1]
PhosphoP[, 3:82] <- apply(PhosphoP[, 3:82], 2, as.numeric)
PhosphoP[, 3:82] <- normalizeBetweenArrays(PhosphoP[, 3:82], method="quantile")

# MAD
Mads <- apply(PhosphoP[, 3:82], 1, mad, na.rm = T)
Th <- quantile(Mads, 0.9, na.rm = T)
Top <- which(Mads >= Th)
PhosphoP_Mads <- PhosphoP[, 3:82][Top, ]
PhosphoP_Mads[PhosphoP_Mads == 'NaN'] <- 0

# NMF
PhosphoP_NMF <- as.matrix(PhosphoP_Mads)
PhosphoP_NMF_pos <- PhosphoP_NMF
PhosphoP_NMF_pos[PhosphoP_NMF_pos < 0] <- 0
PhosphoP_NMF_neg <- - PhosphoP_NMF
PhosphoP_NMF_neg[PhosphoP_NMF_neg < 0] <- 0
PhosphoP_NMF <- rbind(PhosphoP_NMF_pos, PhosphoP_NMF_neg)


# cNMF cluster
k_values <- 2:6
results_list <- list()  # 전체 결과 저장용

num_cores <- 10
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results_list <- foreach(k = k_values, .packages = "CancerSubtypes") %dopar% {
  ExecuteCNMF(PhosphoP_NMF, clusterNum = k, nrun = 200)
}

# names 할당
names(results_list) <- paste0("k_", k_values)

stopCluster(cl)











for (name in names(results_list)) {
  dist_mat <- results_list[[name]]$originalResult@consensus
  print(
    pheatmap(
      dist_mat,
      main = paste0("Consensus value (", name, ")"),
      border_color = NA,
      labels_row = rep('', nrow(dist_mat)),
      labels_col = rep('', ncol(dist_mat))
    )
  )
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




saveRDS(results_list, './rds/PhosphoP.rds')
results_list <- readRDS('./rds/PhosphoP.rds')



group <- factor(results_list[['k_4']]$group)
group <- paste0("Group", group)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(factor(group))

fit <- lmFit(PhosphoP[, 3:82], design)
fit2 <- eBayes(fit)

top_table_all <- topTable(fit2, adjust = "BH", number = Inf)
head(top_table_all)

top_table_all2 <- top_table_all %>%
  filter(adj.P.Val < 0.000001)
sig_P <- rownames(top_table_all2)[!(top_table_all2$Group1 > 0 & top_table_all2$Group2 > 0 & top_table_all2$Group3 > 0 | top_table_all2$Group1 < 0 & top_table_all2$Group2 < 0 & top_table_all2$Group3 < 0)]



# signature 시각화
PhosphoP_Sel <- PhosphoP %>%
  filter(.$Peptide %in% sig_P)

Grp <- factor(results_list[['k_4']]$group, levels = c(1,2,3,4,5), labels = c("Cl1", "Cl2", "Cl3", "Cl4", "Cl5"))
group_colors <- c("Cl1" = "#E41A1C", "Cl2" = "#377EB8", "Cl3" = "#4DAF4A", "Cl4" = "#984EA3", "Cl5" = "#f4d505")
PhosphoP_Sel[is.na(PhosphoP_Sel)] <- 0

Heatmap(as.matrix(PhosphoP_Sel[, 3:82]),
        name = "Global protein",
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
        column_title = "Hierarchical Clustering Heatmap of RNA"
)
