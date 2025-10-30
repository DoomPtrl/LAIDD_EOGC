library(dplyr)
library(readr)
library(circlize)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
set.seed(17)

GlobalP <- read_tsv('./txt/protein_expression_change.txt')
SampleName <- colnames(GlobalP[, 7:86])

RNA <- readRDS('./rds/RNA.rds')
Pro <- readRDS('./rds/GlobalP.rds')
Gly <- readRDS('./rds/GlycoP.rds')
Pho <- readRDS('./rds/PhosphoP.rds')


RNACl <- RNA[['k_2']]$group
ProCl <- Pro[['k_3']]$group
GlyCl <- Gly[['k_2']]$group
PhoCl <- Pho[['k_4']]$group

names(RNACl) <- SampleName
names(ProCl) <- SampleName
names(GlyCl) <- SampleName
names(PhoCl) <- SampleName

RNACl <- as.factor(RNACl)
ProCl <- as.factor(ProCl)
GlyCl <- as.factor(GlyCl)
PhoCl <- as.factor(PhoCl)

RNACl.Mat <- model.matrix(~ RNACl - 1)
ProCl.Mat <- model.matrix(~ ProCl - 1)
GlyCl.Mat <- model.matrix(~ GlyCl - 1)
PhoCl.Mat <- model.matrix(~ PhoCl - 1)

AllCl <- data.frame()
AllCl <- cbind(RNACl.Mat, ProCl.Mat, GlyCl.Mat, PhoCl.Mat)
rownames(AllCl) <- SampleName



results <- ConsensusClusterPlus(t(AllCl),
                                maxK = 6,
                                reps = 1000,
                                pItem = 0.8,
                                pFeature = 1,
                                title = './All',
                                clusterAlg = "pam",
                                distance = "binary",
                                plot = 'png',
                                seed = 17)

col_cluster_labels <- factor(results[[3]]$consensusClass)
cluster_levels <- levels(col_cluster_labels)
cluster_colors <- rainbow(length(cluster_levels))
names(cluster_colors) <- cluster_levels
col_anno <- HeatmapAnnotation(Cluster = col_cluster_labels,
                              col = list(Cluster = cluster_colors),
                              show_annotation_name = TRUE)


results$heatmap <- Heatmap(t(AllCl),
        name = "Binary Matrix",
        col = c("gray5", "red2"),
        show_row_names = T,
        show_column_names = F,
        top_annotation = col_anno,
        cluster_rows = F,
        show_column_dend = F,
        row_order = c('ProCl1', 'PhoCl4', 'PhoCl1', 'GlyCl2', 'RNACl1', 'PhoCl2', 'ProCl2', 'RNACl2', 'GlyCl1', 'ProCl3', 'PhoCl3'),
        column_split = results[[3]]$consensusClass) ; results$heatmap

saveRDS(results, 'clustering2.rds')
