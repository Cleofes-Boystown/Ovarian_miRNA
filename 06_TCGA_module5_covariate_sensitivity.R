suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

base <- "/media/cleofessarmiento/DATA1/Ovarian_cancer_paper all files/ova_paper/supplementary/tcga_validation_gdc/module5_tcga_validation_results"
rdata <- file.path(base, "module5_tcga_validation_results.RData")
expr_fp <- "/media/cleofessarmiento/DATA1/Ovarian_cancer_paper all files/ova_paper/supplementary/tcga_validation/TcgaTargetGtex_RSEM_gene_tpm.gz"
out <- file.path(base, "module5_tcga_nonoverlap_composition_adjusted_results.csv")

e <- new.env()
load(rdata, envir=e)

dt <- as.data.table(e$matched_expr_dt)
dt[, log2_mir15a := log2(mir15a_5p_rpm + 1)]

marker_sets <- list(
  endothelial = c("PECAM1","VWF","CDH5","KDR","FLT1","ESAM","CLDN5","ENG"),
  stromal = c("COL1A1","COL1A2","DCN","LUM","FAP","ACTA2","PDGFRB"),
  immune = c("PTPRC","CD3D","CD3E","CD8A","CD68","MS4A1","NKG7")
)

all_markers <- unique(unlist(marker_sets))

map <- as.data.table(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_markers,
  keytype = "SYMBOL",
  columns = c("SYMBOL","ENSEMBL")
))
map <- unique(map[!is.na(ENSEMBL)])

hdr <- names(fread(expr_fp, nrows=0))
patients <- dt$patient

matched_cols <- sapply(patients, function(p) {
  hits <- grep(paste0("^", p, "-"), hdr, value=TRUE)
  hits <- hits[grepl("-01($|-)", hits)]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
})

cat("Patient column matching:\n")
cat("Patients:", length(patients), "\n")
cat("Matched expression columns:", sum(!is.na(matched_cols)), "\n")
cat("Unmatched:", sum(is.na(matched_cols)), "\n")

keep <- !is.na(matched_cols)
dt2 <- copy(dt[keep])
matched_cols <- matched_cols[keep]

x <- fread(expr_fp, select=c("sample", matched_cols), nThread=2)

gene_col <- names(x)[1]
x[, ensembl_base := sub("\\..*$", "", x[[gene_col]])]

needed_ens <- unique(map$ENSEMBL)
subx <- x[ensembl_base %in% needed_ens]

cat("\nExtracted marker rows:\n")
print(dim(subx))
print(subx[, .(gene_id=get(gene_col), ensembl_base)])

if (nrow(subx) == 0) {
  stop("No marker rows extracted. Ensembl ID matching failed.")
}

mat <- as.matrix(subx[, ..matched_cols])
storage.mode(mat) <- "numeric"

ens_to_sym <- unique(map[, .(ENSEMBL, SYMBOL)])
sym_vec <- ens_to_sym$SYMBOL[match(subx$ensembl_base, ens_to_sym$ENSEMBL)]
rownames(mat) <- sym_vec

symbols <- unique(rownames(mat))
collapsed <- sapply(symbols, function(s) {
  colMeans(mat[rownames(mat) == s, , drop=FALSE], na.rm=TRUE)
})
mat2 <- t(collapsed)
rownames(mat2) <- symbols

z <- t(scale(t(mat2)))

for (nm in names(marker_sets)) {
  present <- intersect(marker_sets[[nm]], rownames(z))
  cat("\n", nm, " markers present: ", paste(present, collapse=", "), "\n", sep="")
  dt2[[paste0(nm, "_score")]] <- if (length(present) > 0) {
    colMeans(z[present, , drop=FALSE], na.rm=TRUE)
  } else {
    NA_real_
  }
}

run_model <- function(formula, data) {
  fit <- lm(formula, data=data)
  co <- summary(fit)$coefficients
  data.frame(
    model = deparse(formula),
    term = rownames(co),
    estimate = co[, "Estimate"],
    std_error = co[, "Std. Error"],
    t_value = co[, "t value"],
    p_value = co[, "Pr(>|t|)"],
    r_squared = summary(fit)$r.squared,
    adj_r_squared = summary(fit)$adj.r.squared,
    stringsAsFactors=FALSE
  )
}

res <- rbindlist(list(
  unadjusted = run_model(module5_score ~ log2_mir15a, dt2),
  endothelial_adjusted = run_model(module5_score ~ log2_mir15a + endothelial_score, dt2),
  stromal_adjusted = run_model(module5_score ~ log2_mir15a + stromal_score, dt2),
  immune_adjusted = run_model(module5_score ~ log2_mir15a + immune_score, dt2),
  full_proxy_adjusted = run_model(module5_score ~ log2_mir15a + endothelial_score + stromal_score + immune_score, dt2)
), idcol="analysis")

write.csv(res, out, row.names=FALSE)

cat("\nWrote:\n", out, "\n\n")

cat("log2_mir15a terms:\n")
print(res[term == "log2_mir15a"])

cat("\nComposition score correlations with module5_score:\n")
print(cor(
  dt2[, .(module5_score, endothelial_score, stromal_score, immune_score)],
  use="pairwise.complete.obs",
  method="spearman"
))
