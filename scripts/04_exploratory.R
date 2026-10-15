# STEP 4 - PCA + hierarchical-clustering heatmap of the most variable probes
source("scripts/utils.R")
suppressPackageStartupMessages(library(pheatmap))

cfg <- load_config(); d <- out_dirs(cfg)
obj <- readRDS(file.path(d$obj, "betas_analysis.rds"))
betas <- obj$betas; pheno <- obj$pheno
cols <- group_colors(pheno$group)

# impute remaining NAs with the probe mean (exploratory plots only; NOT used for testing)
imp <- betas
na <- which(is.na(imp), arr.ind = TRUE)
if (nrow(na)) imp[na] <- rowMeans(betas, na.rm = TRUE)[na[, 1]]
sds <- sqrt(rowSums((imp - rowMeans(imp))^2) / (ncol(imp) - 1))
n_top <- min(cfg$exploratory$n_top_variable, nrow(imp))
top <- names(sort(sds, decreasing = TRUE))[seq_len(n_top)]
write.csv(data.frame(Probe_ID = top, sd_beta = sds[top]), file.path(d$tab, "top_variable_probes.csv"), row.names = FALSE)
write.csv(imp[top, ], file.path(d$tab, "top_variable_betas.csv"))   # input for UMAP (python)

# ---- PCA ----------------------------------------------------------------------
pca <- prcomp(t(imp[top, ]), center = TRUE, scale. = FALSE)
ve  <- 100 * pca$sdev^2 / sum(pca$sdev^2)
sc  <- data.frame(pca$x[, 1:min(5, ncol(pca$x)), drop = FALSE], group = pheno$group, label = pheno$label)
write.csv(data.frame(PC = paste0("PC", seq_along(ve)), variance_pct = ve),
          file.path(d$tab, "pca_variance_explained.csv"), row.names = FALSE)
write.csv(sc, file.path(d$tab, "pca_scores.csv"), row.names = FALSE)
p <- ggplot(sc, aes(PC1, PC2, colour = group)) +
  geom_point(size = 3.5) +
  ggrepel::geom_text_repel(aes(label = label), size = 3.2, show.legend = FALSE) +
  scale_colour_manual(values = cols) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), colour = "Group",
       title = paste0("PCA - top ", n_top, " most variable probes")) +
  theme_bw(base_size = 12)
save_plot(p, file.path(d$fig, "04_PCA"), 6.5, 5)

# ---- Heatmap ------------------------------------------------------------------
ann <- data.frame(Group = pheno$group, row.names = rownames(pheno))
ph <- pheatmap(imp[top, ], color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
               clustering_method = "ward.D2", clustering_distance_cols = "correlation",
               clustering_distance_rows = "euclidean", show_rownames = FALSE,
               annotation_col = ann, annotation_colors = list(Group = cols),
               main = paste0("Top ", n_top, " most variable probes (beta)"), silent = TRUE)
pdf(file.path(d$fig, "05_heatmap_top_variable_probes.pdf"), 7, 8); grid::grid.newpage(); grid::grid.draw(ph$gtable); invisible(dev.off())
png(file.path(d$fig, "05_heatmap_top_variable_probes.png"), 7, 8, units = "in", res = 300); grid::grid.newpage(); grid::grid.draw(ph$gtable); invisible(dev.off())
write_session_info(cfg, "04_exploratory")
