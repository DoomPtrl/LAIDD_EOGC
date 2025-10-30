## =========================================================
suppressPackageStartupMessages({
  library(dplyr); library(tibble)
  library(clusterProfiler); library(org.Hs.eg.db)
  library(msigdbr); library(enrichplot); library(ggplot2); library(patchwork)
})

# -------- 기본 입력 --------
stopifnot(exists("cor_df"), all(c("gene","rho") %in% names(cor_df)))
if (!exists("gene_symbols")) gene_symbols <- unique(cor_df$gene)

# 0) 랭크 정렬 & 상/하위 1000 (양/음 모두 포함)
df_all <- cor_df %>% dplyr::filter(is.finite(rho)) %>% dplyr::arrange(dplyr::desc(rho))
N_use  <- min(1000, nrow(df_all) %/% 2)
top_genes <- head(df_all$gene, N_use)
bot_genes <- tail(df_all$gene, N_use)
cat(sprintf("Top=%d, Bottom=%d (of %d ranked genes)\n", length(top_genes), length(bot_genes), nrow(df_all)))

# 1) SYMBOL -> ENTREZ + Universe (실제 테스트 ∩ KEGG)
mp <- bitr(unique(gene_symbols), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

## --- msigdbr: KEGG LEGACY 우선 선택 (버전 호환) ---
coll <- msigdbr::msigdbr_collections(); cols <- names(coll)
ccol <- if ("collection"       %in% cols) "collection"       else if ("gs_collection" %in% cols) "gs_collection" else if ("gs_cat" %in% cols) "gs_cat" else stop("no collection col")
scol <- if ("subcollection"    %in% cols) "subcollection"    else if ("gs_subcollection" %in% cols) "gs_subcollection" else if ("gs_subcat" %in% cols) "gs_subcat" else stop("no subcollection col")

subvals   <- unique(coll[[scol]][coll[[ccol]] == "C2"])
kegg_order <- c("CP:KEGG_LEGACY","CP:KEGG_MEDICUS","CP:KEGG")
kegg_sub  <- intersect(kegg_order, subvals); stopifnot(length(kegg_sub) > 0)

m_kegg <- msigdbr::msigdbr(species = "human", collection = "C2", subcollection = kegg_sub[1])
gcol <- intersect(c("entrez_gene","ncbi_gene","entrez_gene_id"), names(m_kegg))[1]
stopifnot(length(gcol) == 1)

term2g <- m_kegg[, c("gs_name", gcol)]
names(term2g) <- c("gs_name","entrez_gene")

# --- 질병 관련 KEGG 제외 (요청 유지) ---
disease_regex <- "(?i)(DISEASE|PATHWAYS?_IN_CANCER|CANCER|CARCINOMA|SARCOMA|NEOPLAS|TUMOR|INFECTION|INFECTIOUS|VIR(AL|US)|BACTER|FUNGAL|PARASIT|PROTOZOA|HIV|AIDS|INFLUENZA|MEASLES|MUMPS|RUBELLA|EBOLA|DENGUE|ZIKA|COVID|CORONAVIRUS|HEPATITIS|HPV|H\\.?_?PYLORI|HELICOBACTER|SALMONELLA|SHIGELLA|STAPHYLOCOCCUS|STREPTOCOCCUS|CHOLERA|TUBERCULOSIS|TOXOPLASMOSIS|ALZHEIMER|PARKINSON|HUNTINGTON|EPILEPSY|SCHIZOPHRENIA|CARDIOMYOPATHY|ARRHYTHMIA|ATHEROSCLEROSIS|THROMBOSIS|HYPERTENS|OBESITY|DIABETES|LUPUS|ARTHRITIS|PSORIASIS|NEPHROPATHY|RETINOPATHY|ENCEPHAL)"
term2g <- term2g %>% dplyr::filter(!grepl(disease_regex, gs_name, perl = TRUE))

# Universe = 테스트 전체 ∩ (질병 제외된) KEGG
universe <- intersect(unique(mp$ENTREZID), unique(term2g$entrez_gene))

to_entrez <- function(syms){
  out <- unique(mp$ENTREZID[match(syms, mp$SYMBOL)]); out[!is.na(out)]
}
g_top <- to_entrez(top_genes)
g_bot <- to_entrez(bot_genes)

# 2) ORA (네 방식 그대로)
run_enrich <- function(gset){
  if (length(gset) < 10) return(NULL)
  ek <- enricher(gene = gset, TERM2GENE = term2g, universe = universe,
                 pAdjustMethod = "BH", pvalueCutoff = 0.05)
  if (is.null(ek)) return(NULL)
  setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}
ek_top <- run_enrich(g_top)
ek_bot <- run_enrich(g_bot)

# 3) 이름 예쁘게(접두어 제거, Title Case) — 네 함수 유지
clean_kegg <- function(x){
  x <- sub("^KEGG(\\s+LEGACY|\\s+MEDICUS)?(\\s+REFERENCE)?\\s*", "", x)
  x <- gsub("_", " ", x); x <- trimws(x)
  tools::toTitleCase(tolower(x))
}
fix_desc <- function(ek){
  if (is.null(ek)) return(NULL)
  ek@result$Description <- clean_kegg(ek@result$Description)
  ek
}
ek_top <- fix_desc(ek_top); ek_bot <- fix_desc(ek_bot)

# 4) Strong/Weak 각각 유의성 높은 순서로 4개만 선택 (결과 테이블 순서 그대로)
pick_topk_by_result_order <- function(ek, k=4){
  if (is.null(ek)) return(data.frame())
  df <- as.data.frame(ek@result)
  if (!nrow(df)) return(df[0,])
  df[seq_len(min(k, nrow(df))), , drop = FALSE]
}
top4_df <- pick_topk_by_result_order(ek_top, 4)
bot4_df <- pick_topk_by_result_order(ek_bot, 4)

ids_top <- as.character(top4_df$ID)
ids_bot <- as.character(bot4_df$ID)
labs_top <- setNames(as.character(top4_df$Description), ids_top)
labs_bot <- setNames(as.character(bot4_df$Description), ids_bot)

# 5) 랭크/색 데이터 (A모드용)
rank_df <- df_all %>% dplyr::mutate(rank = dplyr::row_number(),
                                    pos  = pmax(rho, 0),
                                    neg  = pmin(rho, 0))
n_tot <- nrow(rank_df)
rho_by_gene <- setNames(rank_df$rho, rank_df$gene)

# ENTREZ -> SYMBOL, set 멤버 전체
ent2sym  <- mp %>% dplyr::distinct(ENTREZID, SYMBOL)
term2sym <- term2g %>% dplyr::left_join(ent2sym, by = c("entrez_gene" = "ENTREZID")) %>%
  dplyr::filter(!is.na(SYMBOL)) %>% dplyr::distinct(gs_name, SYMBOL)

# 6) A모드 바코드 데이터(멤버 전부, 색은 각 유전자 rho 부호)
build_ticks_all <- function(ids, labels){
  if (!length(ids)) return(tibble())
  dplyr::bind_rows(Map(function(id_one, label_one){
    syms <- term2sym$SYMBOL[term2sym$gs_name == id_one]
    pos  <- match(syms, rank_df$gene); pos <- pos[!is.na(pos)]
    if (!length(pos)) return(tibble())
    side <- ifelse(rho_by_gene[rank_df$gene[pos]] >= 0, "Top", "Bottom")
    tibble(x = pos, term = label_one, side = side)
  }, ids, labels))
}
ticks_top <- build_ticks_all(ids_top, labs_top)
ticks_bot <- build_ticks_all(ids_bot, labs_bot)
ticks_all <- dplyr::bind_rows(ticks_top, ticks_bot)

# 7) y축 순서(위: strong 4, 아래: weak 4)
order_terms <- c(unname(labs_top), unname(labs_bot))
term_index  <- setNames(seq_along(order_terms), order_terms)
if (nrow(ticks_all)) ticks_all$y <- unname(term_index[ticks_all$term])

# 8) 플롯 (곡선 + 바코드)
p_curve <- ggplot(rank_df, aes(x = rank)) +
  geom_area(aes(y = pos), fill = "#2C7FB8") +
  geom_area(aes(y = neg), fill = "#C07B33") +
  labs(y = "Spearman's correlation coefficient", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank())

p_barcode <- if (nrow(ticks_all)) {
  ggplot(ticks_all, aes(x = x, y = y, color = side)) +
    geom_linerange(aes(ymin = y - 0.45, ymax = y + 0.45),
                   size = 0.35, alpha = 0.95) +
    scale_x_continuous(limits = c(1, n_tot), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(order_terms),
                       labels = order_terms,
                       trans = "reverse",
                       position = "right") +
    scale_color_manual(values = c(Top = "#2C7FB8", Bottom = "#C07B33")) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(t = 0, r = 5, b = 0, l = 0))
} else {
  ggplot() + annotate("text", 0.5, 0.5,
                      label = "No KEGG pathways (disease-excluded) passed the cutoff.",
                      size = 5) + theme_void()
}

p_full <- p_curve / p_barcode + plot_layout(heights = c(3, 1.6))
print(p_full)

# 9) 저장(그림 + 표)
ggsave("kegg_barcode_top4_bot4.png", p_full, width = 11, height = 6.2, dpi = 300)
ggsave("kegg_barcode_top4_bot4.pdf",  p_full, width = 11, height = 6.2)

write.csv(top4_df, "pathways_selected_strong_top4.csv", row.names = FALSE)
write.csv(bot4_df, "pathways_selected_weak_top4.csv",   row.names = FALSE)

if (!is.null(ek_top)) write.csv(as.data.frame(ek_top@result), "enrichment_strong_all.csv", row.names = FALSE)
if (!is.null(ek_bot)) write.csv(as.data.frame(ek_bot@result), "enrichment_weak_all.csv",   row.names = FALSE)

# 10) 유전자 ↔ 경로 매핑 저장 (선택 8개 / 전체 유의)
rank_vec <- setNames(rank_df$rho, rank_df$gene)

build_gene_map <- function(ids, labels, group_name){
  if (!length(ids)) return(tibble())
  dplyr::bind_rows(Map(function(id_one, label_one){
    entrez <- term2g$entrez_gene[term2g$gs_name == id_one]
    syms_all <- unique(mp$SYMBOL[match(entrez, mp$ENTREZID)])
    syms_all <- syms_all[!is.na(syms_all)]
    syms_in_rank <- syms_all[syms_all %in% names(rank_vec)]
    if (!length(syms_in_rank)) return(tibble())
    rho <- as.numeric(rank_vec[syms_in_rank])
    side <- ifelse(rho >= 0, "Top", "Bottom")
    in_input <- ifelse(syms_in_rank %in% top_genes, "TopSet",
                       ifelse(syms_in_rank %in% bot_genes, "BottomSet", "None"))
    tibble(pathway_id = id_one,
           pathway    = label_one,
           group      = group_name,      # Strong/Weak 선택군
           gene       = syms_in_rank,
           rho        = rho,
           side       = side,
           in_input_set = in_input)
  }, ids, labels))
}
gene_map_selected <- dplyr::bind_rows(
  build_gene_map(ids_top, labs_top, "Strong"),
  build_gene_map(ids_bot, labs_bot, "Weak")
)
write.csv(gene_map_selected, "gene_to_pathways_selected8.csv", row.names = FALSE)

build_gene_map_from_result <- function(ek, group_name){
  if (is.null(ek)) return(tibble())
  df <- as.data.frame(ek@result); if (!nrow(df)) return(tibble())
  ids  <- as.character(df$ID)
  labs <- setNames(as.character(df$Description), ids)
  build_gene_map(ids, labs, group_name)
}
gene_map_all <- dplyr::bind_rows(
  build_gene_map_from_result(ek_top, "Strong"),
  build_gene_map_from_result(ek_bot, "Weak")
)
write.csv(gene_map_all, "gene_to_pathways_all_enriched.csv", row.names = FALSE)

message("Done. Files saved:
- kegg_barcode_top4_bot4.png / .pdf
- pathways_selected_strong_top4.csv, pathways_selected_weak_top4.csv
- enrichment_strong_all.csv, enrichment_weak_all.csv
- gene_to_pathways_selected8.csv
- gene_to_pathways_all_enriched.csv")
