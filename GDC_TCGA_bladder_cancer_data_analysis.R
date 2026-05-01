#================================ssGSEA==========================
# Load libraries
library(GSVA)
library(GSEABase)
library(Biobase)
library(pheatmap)
library(ggplot2)
library(dplyr)

#====================================== Load and preprocess data ===========

#All data used here are downloaded from https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Bladder%20Cancer%20(BLCA)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

# Step 1: Load expression data 
count <- read.delim("TCGA-BLCA.star_counts.tsv.gz")


ann <- read.delim("gencode.v36.annotation.gtf.gene.probemap")

# Strip version numbers from Ensembl_IDs in both count and ann
count$Ensembl_ID <- sub("\\..*", "", count$Ensembl_ID)

ann$id_clean <- sub("\\..*", "", ann$id)

# Merge annotation with counts by Ensembl ID
count <- merge(ann, count, by.x = "id_clean", by.y = "Ensembl_ID")

rownames(count) <- make.names(count$gene, unique = TRUE)

write.csv(count, "exp_rownames_counts.csv")

exp_data <- read.csv("exp_rownames_counts.csv", row.names = 1)
exp_data <- exp_data[!grepl("^ENSG", rownames(exp_data)), ]
colnames(exp_data) <- gsub("\\.", "-", colnames(exp_data))

# Step 2: Load metadata
meta <- read.delim(gzfile("TCGA-BLCA.clinical.tsv.gz"), sep = "\t", stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
primary_meta <- meta[meta$tissue_type.samples == "Tumor", ]
colnames(primary_meta)
table(primary_meta$ajcc_pathologic_t.diagnoses)

# Step 3: Subset expression and metadata to matched primary tumor samples
common_samples <- intersect(colnames(exp_data), rownames(primary_meta))
exp_data <- exp_data[, common_samples]
primary_meta <- primary_meta[common_samples, ]

all(colnames(exp_data) == rownames(primary_meta))

# Step 4: Prepare expression matrix
expr_mat <- as.matrix(exp_data)
stopifnot(mode(expr_mat) == "numeric", is.matrix(expr_mat), !any(is.na(expr_mat)))

# Step 5: Define hypoxia hallmark gene set
# Source: MSigDB Hallmark gene set: HALLMARK_HYPOXIA

hypoxia_genes <- read.csv("hypoxia_geneset.csv", header = TRUE)
gene_list <- as.character(hypoxia_genes[[1]])
hypoxia_list <- list(Hypoxia = gene_list)

# Step 6: Run ssGSEA

param <- ssgseaParam(expr_mat, hypoxia_list)
ssgsea_scores <- gsva(param)
ssgsea_t <- t(ssgsea_scores)

# Step 7: Merge hypoxia score with metadata
primary_meta$HypoxiaScore <- ssgsea_t[, "Hypoxia"]

summary(primary_meta$HypoxiaScore)

# Optional: Define groups by score (median split)
primary_meta$HypoxiaGroup <- ifelse(
  primary_meta$HypoxiaScore > median(primary_meta$HypoxiaScore),
  "High",
  "Low"
)
table(primary_meta$HypoxiaGroup)


write.csv(primary_meta, "Hypoxia_score_final_metadata.csv")


# Order samples by hypoxia score
primary_meta_ordered <- primary_meta[order(primary_meta$HypoxiaScore), ]

# Make sample a factor so ggplot respects the order
primary_meta_ordered$sample <- factor(
  primary_meta_ordered$sample,
  levels = primary_meta_ordered$sample
)

primary_meta_ordered$HypoxiaScore_z <-
  scale(primary_meta_ordered$HypoxiaScore)

ggplot(primary_meta_ordered,
       aes(x = sample, y = HypoxiaScore_z, fill = HypoxiaGroup)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Low" = "#1F77B4", "High" = "#D62728")) +
  theme_classic() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = "Sample",
    y = "Hypoxia score (z-score)",
    fill = "Group",
    title = "Hypoxia score across samples"
  )



#================================== Box plot ===================================


hypoxia_scores <- data.frame(
  SampleID = colnames(ssgsea_scores),
  HypoxiaScore = as.numeric(ssgsea_scores[1, ])
)

# Add group information from metadata
hypoxia_scores$RiskGroup <- annotation_col$RiskGroup[
  match(hypoxia_scores$SampleID, rownames(annotation_col))
]

# Check the structure
str(hypoxia_scores)
head(hypoxia_scores)

# Boxplot comparing Hypoxia Score by Risk Group
library(ggplot2)

ggplot(primary_meta, aes(x = HypoxiaGroup, y = HypoxiaScore, fill = HypoxiaGroup)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = c("High" = "#D62728", "Low" = "#1F77B4")) +
  labs(
    title = "Hypoxia ssGSEA Score",
    x = "Risk Group",
    y = "Hypoxia ssGSEA Score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# Perform statistical test (Wilcoxon rank-sum test)
test_result <- wilcox.test(HypoxiaScore ~ HypoxiaGroup, data = primary_meta)
test_result$p.value  # show p-value

# Optional: Add p-value annotation to the plot
library(ggpubr)
ggplot(primary_meta, aes(x = HypoxiaGroup, y = HypoxiaScore, fill = HypoxiaGroup)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = c("High" = "#D62728", "Low" = "#1F77B4")) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(
    title = "Hypoxia ssGSEA Score by Hypoxia Score Group",
    x = "Hypoxia Group",
    y = "Hypoxia ssGSEA Score"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.4)
  )

#================================survival=======================================

surv_data <- read.delim("TCGA-BLCA.survival.tsv.gz")
surv_data$sample <- gsub("\\.", "-", surv_data$sample)


# Add sample ID to metadata

primary_meta$sample <- rownames(primary_meta)

# Merge survival data with hypoxia group
merged_surv <- merge(surv_data, primary_meta[, c("sample", "HypoxiaScore" ,"HypoxiaGroup")], by = "sample")

# Check structure
head(merged_surv)

write.csv(merged_surv, "merged_surv.csv")

library(survival)
library(survminer)


cutp_hypoxia <- surv_cutpoint(
  merged_surv,
  time = "OS.time",
  event = "OS",
  variables = "HypoxiaScore",
  minprop = 0.1  # Ensures groups have at least 10% of samples (adjustable)
)

# Extract the optimal threshold
opt_thr_hypoxia <- cutp_hypoxia$cutpoint["HypoxiaScore", "cutpoint"]
print(opt_thr_hypoxia)


# Assign group based on optimal cutoff
merged_surv$hypoxia_group_opt <- ifelse(
  merged_surv$HypoxiaScore > median(merged_surv$HypoxiaScore), "High", "Low"
)



# Convert to factor for KM plot
merged_surv$hypoxia_group_opt <- factor(merged_surv$hypoxia_group_opt, levels = c("Low", "High"))


# Create survival object
surv_obj_hypoxia <- Surv(time = merged_surv$OS.time, event = merged_surv$OS)

# Fit KM model
fit_hypoxia_opt <- survfit(surv_obj_hypoxia ~ hypoxia_group_opt, data = merged_surv)

# Plot KM curve

ggsurvplot(
  fit_hypoxia_opt,
  data = merged_surv,
  pval = TRUE,
  risk.table = FALSE,
  palette = c("#D91656", "#3ABEF9"),
  title = "Survival Based on Optimized Hypoxia Score ",
  legend.title = "Hypoxia Group",
  legend.labs = c("High", "Low")
)


#=================================DEG===========================================

library(limma)
library(DESeq2)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)

# Ensure expression matrix matches metadata
primary_meta <- primary_meta[colnames(exp_data), ]
stopifnot(all(colnames(exp_data) == rownames(primary_meta)))

# Log2 transform expression
expr_log2 <- voom(exp_data)


# Create design matrix for ssGSEA Group (e.g., High vs Low Hypoxia)
primary_meta$HypoxiaGroup <- factor(primary_meta$HypoxiaGroup, levels = c("Low", "High"))
design <- model.matrix(~ 0 + HypoxiaGroup, data = primary_meta)
colnames(design) <- gsub("HypoxiaGroup", "", colnames(design))  # Columns: "Low", "High"

# Fit model
fit <- lmFit(expr_log2, design)

# Contrast: High vs Low
contrast.matrix <- makeContrasts(HighvsLow = High - Low, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract DEGs
deg_results <- topTable(fit2, adjust.method = "BH", number = Inf, sort.by = "P")
deg_results <- deg_results %>% rownames_to_column(var = "SYMBOL") %>% na.omit()

# Add expression status
deg_results <- deg_results %>%
  mutate(
    nlog10 = -log10(P.Value),
    expression = case_when(
      P.Value < 0.05 & logFC >= 0.5  ~ "Up-regulated",
      P.Value < 0.05 & logFC <= -0.5 ~ "Down-regulated",
      TRUE                           ~ "Stable"
    )
  )

# Check DEG summary
table(deg_results$expression)

# Get top genes for labeling
top_n <- 20
top_deg_genes <- bind_rows(
  deg_results %>%
    filter(expression == "Up-regulated") %>%
    arrange(P.Value, desc(abs(logFC))) %>%
    head(top_n),
  deg_results %>%
    filter(expression == "Down-regulated") %>%
    arrange(P.Value, desc(abs(logFC))) %>%
    head(top_n)
)

# Volcano Plot
volcano_plot <- ggplot(deg_results, aes(x = logFC, y = nlog10, fill = expression)) +
  geom_point(alpha = 0.8, size = 1.8, shape = 21, stroke = 0.05) +
  scale_fill_manual(values = c("Up-regulated" = "#ff6347", "Stable" = "#c4c4c4", "Down-regulated" = "#4682b4")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_deg_genes, aes(label = SYMBOL), size = 2.5, max.overlaps = 100) +
  theme_bw() +
  xlab("logFC (High vs Low ssGSEA group)") +
  ylab("-log10(P-value)") +
  ggtitle("DEGs between High vs Low Hypoxia ssGSEA Group")

print(volcano_plot)

# Save results if needed
write.csv(deg_results, "DEGs_High_vs_Low_ssGSEA.csv", row.names = FALSE)

#================================Enrichment=====================================

# Load libraries
# Ensure we call clusterProfiler::GSEA, not fgsea::gsea
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(forcats)
library(stringr)


# Read DE results and map SYMBOL → ENTREZ
topTable_results <- read.csv("DEGs_High_vs_Low_ssGSEA.csv", row.names = 1)
topTable_results$SYMBOL <- rownames(topTable_results)

# Check existence of SYMBOL and logFC
if(!all(c("SYMBOL","logFC") %in% colnames(topTable_results))) stop("CSV must contain SYMBOL and logFC columns")

# Map SYMBOL to ENTREZID and build ranked list
symbol2entrez <- bitr(
  topTable_results$SYMBOL,
  fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db
)
gene_rank_entrez <- setNames(
  topTable_results$logFC[match(symbol2entrez$SYMBOL, topTable_results$SYMBOL)],
  symbol2entrez$ENTREZID
)


#rank genes based on fold change
gene_rank_entrez <- sort(gene_rank_entrez, decreasing = TRUE)
# Optional: write .rnk
write.table(gene_rank_entrez, file = "gene_rank_all_entrez.rnk", sep = "\t", quote = FALSE, col.names = FALSE)




# 1) Retrieve MSigDB gene sets
msig_all <- msigdbr(species = "Homo sapiens")

# Helper to build TERM2GENE with Entrez IDs from gene_symbol
build_term2gene <- function(df) {
  # df must have columns: gs_name, gene_symbol
  term2symbol <- df %>% select(term = gs_name, symbol = gene_symbol) %>% distinct()
  mapping <- bitr(
    term2symbol$symbol,
    fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db
  )
  term2entrez <- term2symbol %>%
    inner_join(mapping, by = c("symbol" = "SYMBOL")) %>%
    select(term, gene = ENTREZID) %>%
    distinct() %>%
    as.data.frame()
  return(term2entrez)
}

# 2) Build collections
# Hallmark
homoH <- msig_all %>% filter(gs_collection == "H")
homoH_term2gene <- build_term2gene(homoH)
# GO Biological Process
homoGOBP <- msig_all %>% filter(gs_collection == "C5", gs_subcollection == "GO:BP")
homoGOBP_term2gene <- build_term2gene(homoGOBP)
# C2 collections
homoC2 <- msig_all %>% filter(gs_collection == "C2")
homoKEGG <- homoC2 %>% filter(gs_subcollection %in% c("CP:KEGG_LEGACY","CP:KEGG_MEDICUS"))
homoKEGG_term2gene <- build_term2gene(homoKEGG)
homoWiki <- homoC2 %>% filter(gs_subcollection == "CP:WIKIPATHWAYS")
homoWiki_term2gene <- build_term2gene(homoWiki)
homoReactome <- homoC2 %>% filter(gs_subcollection == "CP:REACTOME")
homoReactome_term2gene <- build_term2gene(homoReactome)


# 4) GSEA runner using clusterProfiler::GSEA
run_gsea <- function(geneList, term2gene, prefix, db) {
  res <- clusterProfiler::GSEA(
    geneList     = geneList,
    TERM2GENE    = term2gene,
    pvalueCutoff = 1,
    verbose      = FALSE
  )
  df <- as_tibble(res@result) %>% mutate(database = db)
  write.csv(df, file = paste0(prefix, ".csv"), row.names = FALSE)
  return(df)
}

# Run analyses
gsea_H        <- run_gsea(gene_rank_entrez, homoH_term2gene,        "GSEA_Hallmark_entrez", "Hallmark")
gsea_GOBP     <- run_gsea(gene_rank_entrez, homoGOBP_term2gene,     "GSEA_GOBP_entrez",     "GO:BP")
gsea_KEGG     <- run_gsea(gene_rank_entrez, homoKEGG_term2gene,     "GSEA_KEGG_entrez",     "KEGG")
gsea_Wiki     <- run_gsea(gene_rank_entrez, homoWiki_term2gene,     "GSEA_Wiki_entrez",     "WikiPathways")
gsea_Reactome <- run_gsea(gene_rank_entrez, homoReactome_term2gene, "GSEA_Reactome_entrez", "Reactome")


# 5) Plot top pathways example for each result
theme_set(DOSE::theme_dose(12))
plot_top <- function(df, title, n = 15) {
  top_df <- df %>% top_n(n, NES)
  ggplot(top_df, aes(x = enrichmentScore, y = fct_reorder(Description, enrichmentScore), fill = NES, size = -log10(p.adjust))) +
    geom_point(pch = 21) +
    scale_fill_gradient() +
    scale_size_continuous(range = c(2, 10)) +
    labs(title = title, x = "Enrichment Score", fill = "NES", size = "-log10(adj.p)") +
    scale_y_discrete(labels = function(x) str_wrap(x, 40)) +
    theme_minimal()
}

# Correct variable names and calls:
plot_top(gsea_H,        "Top Hallmark Pathways")
#ggsave("Top_Hallmark_entrez_dotplot.pdf", width = 10, height = 8)
plot_top(gsea_GOBP,     "Top GO Biological Processes")
#ggsave("Top_GOBP_entrez_dotplot.pdf", width = 10, height = 8)
plot_top(gsea_KEGG,     "Top KEGG Pathways")
#ggsave("Top_KEGG_entrez_dotplot.pdf", width = 10, height = 8)
plot_top(gsea_Wiki,     "Top WikiPathways")
#ggsave("Top_Wiki_entrez_dotplot.pdf", width = 10, height = 8)
plot_top(gsea_Reactome, "Top Reactome Pathways")
#ggsave("Top_Reactome_entrez_dotplot.pdf", width = 10, height = 8)
