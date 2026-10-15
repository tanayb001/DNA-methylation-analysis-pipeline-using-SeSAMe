# STEP 8 - Pathway analysis: over-representation (GO BP + MSigDB) and GSEA
source("scripts/utils.R")
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(msigdbr); library(dplyr); library(ggplot2)
})

cfg <- load_config(); d <- out_dirs(cfg); pw <- cfg$pathway
pdir <- ensure_dir(d$fig, "pathways")
ann <- read.csv(file.path(d$tab, "dmp_annotated.csv.gz"), stringsAsFactors = FALSE)
ann <- ann[!is.na(ann$gene) & ann$region %in% pw$use_regions, ]
universe <- unique(ann$gene)
message(length(universe), " genes in the tested universe (regions: ", paste(pw$use_regions, collapse = ", "), ")")

gene_sets <- lapply(list(all   = ann$gene[ann$Significant],
                         hyper = ann$gene[ann$Significant & ann$Direction == "Hyper"],
                         hypo  = ann$gene[ann$Significant & ann$Direction == "Hypo"]), unique)

load_msig <- function(s) {
  m <- tryCatch(msigdbr(species = "Homo sapiens", collection = s$collection, subcollection = s$subcollection),
                error = function(e) msigdbr(species = "Homo sapiens", category = s$collection, subcategory = s$subcollection))
  unique(as.data.frame(m)[, c("gs_name", "gene_symbol")])
}
msig <- setNames(lapply(pw$msigdb, load_msig), vapply(pw$msigdb, function(s) s$name, ""))

save_enrich <- function(obj, name) {
  df <- if (is.null(obj)) data.frame() else as.data.frame(obj)
  if (!nrow(df)) { message("  [", name, "] no terms tested"); return(invisible()) }
  write.csv(df, file.path(d$tab, paste0("pathway_", name, ".csv")), row.names = FALSE)
  n_sig <- sum(df$p.adjust < pw$fdr_cutoff, na.rm = TRUE)
  message("  [", name, "] ", n_sig, " terms with FDR < ", pw$fdr_cutoff)
  if (n_sig == 0) return(invisible())
  p <- tryCatch({
    if (inherits(obj, "gseaResult"))
      dotplot(obj, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)
    else dotplot(obj, showCategory = pw$top_n_plot)
  }, error = function(e) NULL)
  if (!is.null(p)) save_plot(p + ggtitle(name), file.path(pdir, paste0("dotplot_", name)), 9, 7)
  if (inherits(obj, "gseaResult")) {
    ids <- head(df$ID[order(df$pvalue)], 3)
    tryCatch({
      pdf(file.path(pdir, paste0("gseaplot_", name, ".pdf")), 8, 6)
      print(gseaplot2(obj, geneSetID = ids, pvalue_table = TRUE)); invisible(dev.off())
    }, error = function(e) invisible(dev.off()))
  }
}

# ---- Over-representation -------------------------------------------------------
for (nm in names(gene_sets)) {
  g <- gene_sets[[nm]]
  message("ORA set '", nm, "': ", length(g), " genes")
  if (length(g) < pw$min_genes) { message("  skipped (< ", pw$min_genes, " genes)"); next }
  save_enrich(tryCatch(enrichGO(g, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                universe = universe, pAdjustMethod = "BH", pvalueCutoff = 1,
                                qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500),
                       error = function(e) NULL), paste0("ORA_GO_BP_", nm))
  for (s in names(msig))
    save_enrich(tryCatch(enricher(g, universe = universe, TERM2GENE = msig[[s]], pAdjustMethod = "BH",
                                  pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500),
                         error = function(e) NULL), paste0("ORA_", s, "_", nm))
}

# ---- GSEA (all tested genes ranked by signed -log10 P, averaged per gene) -------
rk <- ann %>% group_by(gene) %>%
  summarise(score = mean(sign(delta_beta) * -log10(pmax(P, 1e-300))), .groups = "drop")
ranks <- sort(setNames(rk$score, rk$gene), decreasing = TRUE)
write.csv(data.frame(gene = names(ranks), rank_score = ranks), file.path(d$tab, "gsea_ranked_genes.csv"), row.names = FALSE)
message("GSEA on ", length(ranks), " ranked genes")
set.seed(cfg$seed)
save_enrich(tryCatch(gseGO(ranks, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", minGSSize = 15,
                           maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0,
                           seed = TRUE, verbose = FALSE), error = function(e) { message(e$message); NULL }),
            "GSEA_GO_BP")
for (s in names(msig))
  save_enrich(tryCatch(GSEA(ranks, TERM2GENE = msig[[s]], minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1,
                            pAdjustMethod = "BH", eps = 0, seed = TRUE, verbose = FALSE),
                       error = function(e) { message(e$message); NULL }), paste0("GSEA_", s))
write_session_info(cfg, "08_pathways")
