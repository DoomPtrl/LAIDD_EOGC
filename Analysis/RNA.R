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


Files <- list.files('./txt/RNA/')

RNA_list <- list()

setwd('./txt/RNA/')
RNA_list <- lapply(Files, read_table)
setwd('X:/home/LAIDD_DATA/subtype_all/')

RNA_df <- RNA_list[[1]][, 1:4]
TPMs <- lapply(RNA_list, function(x) x[[6]])
RNA_df <- bind_cols(RNA_df, as.data.frame(TPMs))
colnames(RNA_df) <- c('geneId', 'transcriptId', 'length', 'eff_length', sub('_.*', '', Files))

RNA_df <- RNA_df[rowSums(RNA_df == 0) < 80, ]
RNA_df[, 5:164] <- log2(RNA_df[, 5:164] + 1)
RNA_df <- as.data.frame(RNA_df)

RNA_T <- RNA_df %>% dplyr::select(-matches('N$'))
RNA_N <- RNA_df %>% dplyr::select(-matches('T$'))

RNA_TN <- RNA_df[, 1:4]
RNA_TN[, 5:84] <- - RNA_N[, 5:84] + RNA_T[, 5:84]
RNA_TN[, 5:84] <- normalizeBetweenArrays(RNA_TN[, 5:84], method="quantile")

# MAD
Mads <- apply(RNA_TN[, 5:84], 1, mad, na.rm = T)
Th <- quantile(Mads, 0.9, na.rm = T)
Top <- which(Mads >= Th)
RNA_TN_Mads <- RNA_TN[, 5:84][Top, ]
RNA_TN_Mads <- RNA_TN_Mads[, order(colnames(RNA_TN_Mads))]

# NMF
RNA_TN_NMF <- as.matrix(RNA_TN_Mads)
RNA_TN_NMF_pos <- RNA_TN_NMF
RNA_TN_NMF_pos[RNA_TN_NMF_pos < 0] <- 0
RNA_TN_NMF_neg <- - RNA_TN_NMF
RNA_TN_NMF_neg[RNA_TN_NMF_neg < 0] <- 0
RNA_TN_NMF <- rbind(RNA_TN_NMF_pos, RNA_TN_NMF_neg)


# cNMF cluster
k_values <- 2:6
results_list <- list()  # 전체 결과 저장용

num_cores <- 10
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results_list <- foreach(k = k_values, .packages = "CancerSubtypes") %dopar% {
  ExecuteCNMF(RNA_TN_NMF, clusterNum = k, nrun = 200)
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




saveRDS(results_list, './rds/RNA.rds')
results_list <- readRDS('./rds/RNA.rds')
rownames(RNA_TN) <- RNA_TN$geneId

group <- factor(results_list[['k_2']]$group)
group <- paste0("Group", group)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(factor(group))

fit <- lmFit(RNA_TN[, 5:84], design)
contrast <- makeContrasts(Group1vsGroup2 = Group1 - Group2, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "BH", number = Inf)
head(results)

sig_gene <- rownames(results)[results$adj.P.Val < 0.00001]



# signature 시각화
RNA_TN_Sel <- RNA_TN %>%
  filter(.$geneId %in% sig_gene)

Grp <- factor(results_list[['k_2']]$group, levels = c(1,2,3,4,5), labels = c("Cl1", "Cl2", "Cl3", "Cl4", "Cl5"))
group_colors <- c("Cl1" = "#E41A1C", "Cl2" = "#377EB8", "Cl3" = "#4DAF4A", "Cl4" = "#984EA3", "Cl5" = "#f4d505")

Heatmap(as.matrix(RNA_TN_Sel[, 5:84]),
        name = "RNA",
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
