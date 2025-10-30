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
set.seed(17)

GlobalP <- read_tsv('./txt/protein_expression_change.txt')

# Norm
GlobalP <- as.data.frame(GlobalP)
GlobalP <- GlobalP[rowSums(is.na(GlobalP)) < 40, ]
rownames(GlobalP) <- GlobalP[, 2]
GlobalP[, 7:86] <- apply(GlobalP[, 7:86], 2, as.numeric)
GlobalP[, 7:86] <- normalizeBetweenArrays(GlobalP[, 7:86], method="quantile")

# NMF
GlobalP_NMF <- as.matrix(GlobalP[, 7:86])
GlobalP_NMF_pos <- GlobalP_NMF
GlobalP_NMF_pos[GlobalP_NMF_pos < 0] <- 0
GlobalP_NMF_neg <- - GlobalP_NMF
GlobalP_NMF_neg[GlobalP_NMF_neg < 0] <- 0
GlobalP_NMF <- rbind(GlobalP_NMF_pos, GlobalP_NMF_neg)

# MAD
Mads <- apply(GlobalP_NMF, 1, mad, na.rm = T)
Th <- quantile(Mads, 0.9, na.rm = T)
Top <- which(Mads >= Th)
GlobalP_Mads <- GlobalP_NMF[Top, ]
GlobalP_Mads[is.na(GlobalP_Mads)] <- 0


# Consensus cluster
GlobalPRes <- ConsensusClusterPlus(
  GlobalP_Mads,
  maxK = 6, reps = 100,
  pItem = 0.8, pFeature = 1,
  clusterAlg = 'km', distance = 'euclidean',
  title = 'GlobalProteomics/10p', plot = 'png',
  seed = 17
)

# Consensus plot
GlobalPRes$CopheneticValues <- sapply(2:6, function(k) {
  Consensus_mat <- GlobalPRes[[k]]$consensusMatrix
  Dist_mat <- as.dist(1 - Consensus_mat)
  hc <- hclust(Dist_mat, method = "complete")
  CopheneticDist <- cophenetic(hc)
  cor(Dist_mat, CopheneticDist)
})

GlobalPRes$CV_df <- data.frame(k = 2:6, cophenetic_correlation = GlobalPRes$CopheneticValues)

GlobalPRes$CVPlot <- ggplot(GlobalPRes$CV_df, aes(x = k, y = cophenetic_correlation)) +
  geom_line() +
  geom_point(size = 3) +
  labs(title = "Global Proteomics",
       x = "k (Number of Clusters)",
       y = "Cophenetic Correlation") +
  theme_minimal() ; GlobalPRes$CVPlot

# Save
saveRDS(GlobalPRes, './GlobalProteomics/GlobalP_10p.rds')



# visualization (Heatmap)
Mads <- apply(GlobalP[, 7:86], 1, mad, na.rm = T)
Th <- quantile(Mads, 0.9, na.rm = T)
Top <- which(Mads >= Th)
GlobalP_Vis <- GlobalP[Top, ]
GlobalP_Vis[is.na(GlobalP_Vis)] <- 0

k = 3
Grp <- factor(GlobalPRes[[k]]$consensusClass, levels = c(1,2,3,4,5,6), labels = c("Cl1", "Cl2", "Cl3", "Cl4", "Cl5", "Cl6"))
GroupColors <- c("Cl1" = "#E41A1C", "Cl2" = "#377EB8", "Cl3" = "#4DAF4A", "Cl4" = "#984EA3", "Cl5" = "#f4d505", "Cl6" = "#24f505")

Heatmap(as.matrix(GlobalP_Vis[, 7:86]),
        name = "Global Proteomics",
        col = colorRamp2(c(-1.5, 0, 1.5), c("blue2", "gray10", "red2")),
        show_row_names = FALSE,
        show_column_names = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        show_column_dend = F, 
        cluster_columns = TRUE,
        cluster_rows = T,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete",
        show_row_dend = F,
        column_split  = Grp,
        row_title = NULL,
        row_title_side = "left",
        row_gap = unit(5, "mm"),
        top_annotation = HeatmapAnnotation(Group = Grp, col = list(Group = GroupColors)),
        column_title = "Hierarchical Clustering Heatmap of Global Proteins"
)