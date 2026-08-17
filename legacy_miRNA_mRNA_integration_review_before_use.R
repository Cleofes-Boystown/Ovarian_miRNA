# --- BEGIN SCRIPT: OVARIAN-ONLY miRNA–mRNA INTEGRATION ---

# 1) Load libraries
library(readr);   library(dplyr)
library(tidyr);   library(purrr)
library(tibble);  library(multiMiR)
library(igraph);  library(ggraph)
library(clusterProfiler); library(org.Hs.eg.db)
library(ggplot2)

# 2) Define paths
proj_dir   <- "/media/cleofessarmiento/DATA1/Ov_integration"
miRNA_csv  <- file.path(proj_dir, "log2_normalised_miRNA_only.csv")
mrna_csv   <- file.path(proj_dir, "../Normalized_Log2_Count_Matrix.csv")
sample_csv <- file.path(proj_dir, "samplesheet.csv")
de_miRNA   <- file.path(proj_dir, "DE_Ov_vs_Ctrl_miRNA_Combined.csv")
de_mRNA    <- file.path(proj_dir, "DESeq2_Ovarian_vs_Control_Annotated.csv")

out_edges  <- file.path(proj_dir, "edges_ovarian_only.csv")
out_full   <- file.path(proj_dir, "network_ovarian_only.png")
out_core   <- file.path(proj_dir, "core_network_ovarian_only.png")
out_go     <- file.path(proj_dir, "GO_BP_ovarian_only.png")

# 3) Load & filter sample sheet to ovarian tumors only
samplesheet <- read_csv(sample_csv, show_col_types=FALSE) %>%
  filter(phenotype == "ovarian_cancer") %>%
  mutate(clean_prefix = toupper(gsub("\\d+$","", gsub("_","", sampleID))))

ov_samples  <- samplesheet$sampleID

# 4) Build miRNA matrix only for ovarian samples
mirna_mat <- read_csv(miRNA_csv, show_col_types=FALSE) %>%
  rename(miRNA = Name) %>%
  pivot_longer(-miRNA, names_to="RCC_Filename", values_to="expr") %>%
  left_join(select(samplesheet, RCC_Filename=IDFILE, SampleID=sampleID),
            by="RCC_Filename") %>%
  filter(SampleID %in% ov_samples) %>%
  select(miRNA, SampleID, expr) %>%
  pivot_wider(names_from=SampleID, values_from=expr, values_fn=mean) %>%
  filter(!is.na(miRNA)) %>%                                       # drop NA rows
  column_to_rownames("miRNA")

# 5) Load DE lists (computed vs. SRA controls)
sig_mirnas <- read_csv(de_miRNA, show_col_types=FALSE) %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  pull(miRNA)

de_mrna    <- read_csv(de_mRNA, show_col_types=FALSE)
sig_mrnas  <- de_mrna %>%
  filter(abs(log2FoldChange) > 1, padj < 0.05) %>%
  pull(gene_name)

# 6) Build mRNA matrix only for ovarian samples
mrna_mat <- read_csv(mrna_csv, show_col_types=FALSE) %>%
  rename(gene_id = `...1`) %>%
  left_join(select(de_mrna, gene_id, gene_name), by="gene_id") %>%
  filter(gene_name %in% sig_mrnas) %>%
  pivot_longer(-c(gene_id, gene_name), names_to="BAM_Filename", values_to="expr") %>%
  mutate(prefix = sub("_Aligned.*$", "", BAM_Filename)) %>%
  left_join(select(samplesheet, clean_prefix, SampleID=sampleID),
            by=c("prefix"="clean_prefix")) %>%
  filter(SampleID %in% ov_samples) %>%
  select(gene_name, SampleID, expr) %>%
  pivot_wider(names_from=SampleID, values_from=expr, values_fn=mean) %>%
  filter(!is.na(gene_name)) %>%                                   # drop NA rows
  column_to_rownames("gene_name")

# 7) Fetch validated targets between your DE miRNAs & DE mRNAs
validated_targets <- get_multimir(
  mirna      = sig_mirnas,
  table      = "validated",
  legacy.out = FALSE
)
targets_df <- as_tibble(validated_targets@data) %>%
  filter(mature_mirna_id %in% rownames(mirna_mat),
         target_symbol   %in% rownames(mrna_mat)) %>%
  select(mature_mirna_id, target_symbol) %>%
  distinct()

# 8) Compute Spearman correlations across ovarian samples
cor_results <- purrr::pmap_dfr(
  targets_df,
  function(mature_mirna_id, target_symbol) {
    x <- as.numeric(mirna_mat[mature_mirna_id, ])
    y <- as.numeric(mrna_mat[target_symbol, ])
    idx <- which(!is.na(x) & !is.na(y))
    if (length(idx) < 3) {
      return(tibble(miRNA=mature_mirna_id, mRNA=target_symbol,
                    rho=NA_real_, pval=NA_real_))
    }
    ct <- cor.test(x[idx], y[idx], method="spearman")
    tibble(miRNA = mature_mirna_id,
           mRNA   = target_symbol,
           rho    = ct$estimate,
           pval   = ct$p.value)
  }
) %>%
  filter(!is.na(rho), rho < -0.3, pval < 0.05)

# 9) Save edges table
write_csv(cor_results, out_edges)
message("Edges written to: ", out_edges)

# 10) Full network plot
png(out_full, width=2400, height=1600, res=300)
ggraph(
  graph_from_data_frame(select(cor_results, from=miRNA, to=mRNA),
                        directed=TRUE),
  layout="fr"
) +
  geom_edge_link(alpha=0.3) +
  geom_node_point(aes(color=ifelse(name %in% sig_mirnas,"miRNA","mRNA")),
                  size=2) +
  theme_void() +
  ggtitle("Ovarian-only Network (ρ < –0.3)")
dev.off()
message("Full network plot: ", out_full)

# 11) Core network (ρ < -0.7, top-10 miRNA hubs)
core_edges <- cor_results %>% filter(rho < -0.7)
g_core     <- graph_from_data_frame(select(core_edges, from=miRNA, to=mRNA),
                                    directed=FALSE)
vertex_attr(g_core,"type")   <- ifelse(vertex_attr(g_core,"name") %in% sig_mirnas,
                                       "miRNA","mRNA")
vertex_attr(g_core,"degree") <- degree(g_core, mode="all")

miinds <- which(vertex_attr(g_core,"type")=="miRNA")
namesv <- vertex_attr(g_core,"name")
degv   <- vertex_attr(g_core,"degree")
top_mi <- namesv[miinds][order(degv[miinds],decreasing=TRUE)[1:10]]

nodes  <- base::unique(c(top_mi,
              unlist(lapply(top_mi, function(x) namesv[neighbors(g_core,x)]))))
subg   <- induced_subgraph(g_core, vids=nodes)
vertex_attr(subg,"label")  <- ifelse(vertex_attr(subg,"name") %in% top_mi,
                                     vertex_attr(subg,"name"), NA)
vertex_attr(subg,"degree") <- degree(subg, mode="all")

png(out_core, width=2400, height=1600, res=300)
ggraph(subg, layout="fr") +
  geom_edge_link(alpha=0.4) +
  geom_node_point(aes(size=degree, color=type)) +
  geom_node_text(aes(label=label), repel=TRUE) +
  scale_color_manual(values=c(miRNA="tomato", mRNA="steelblue")) +
  scale_size(range=c(4,10)) +
  theme_void() +
  labs(title="Core Network (ρ < -0.7, Top-10 miRNA Hubs)",
       color="", size="Degree")
dev.off()
message("Core network plot: ", out_core)

# 12) GO-BP enrichment & dotplot
entrez_ids <- bitr(unique(cor_results$mRNA), fromType="SYMBOL", toType="ENTREZID",
                   OrgDb=org.Hs.eg.db)
go_enrich <- enrichGO(gene=entrez_ids$ENTREZID, OrgDb=org.Hs.eg.db,
                      keyType="ENTREZID", ont="BP", pAdjustMethod="BH",
                      qvalueCutoff=0.05, readable=TRUE)
png(out_go, width=2000, height=1200, res=300)
dotplot(go_enrich, showCategory=15) +
  ggtitle("GO Biological Processes (Ovarian Only)")
dev.off()
message("GO–BP plot: ", out_go)

# --- END SCRIPT ---
