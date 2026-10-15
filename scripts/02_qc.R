# STEP 2 - Initial QC: detection p-value per array
#  * primary metric  : negative-control based detection p-value (SeSAMe detectionPnegEcdf)
#                      = same concept as the minfi / RnBeads detection p-value used in the paper
#  * also reported   : pOOBAH (SeSAMe's own, stricter for FFPE) mean p and fraction detected
source("scripts/utils.R")
suppressPackageStartupMessages({ library(sesame); library(BiocParallel) })

cfg   <- load_config(); d <- out_dirs(cfg); q <- cfg$qc
sheet <- read_sample_sheet(cfg)
prefixes <- resolve_prefixes(sheet, cfg$paths$work_idat_dir)
bp <- make_bp(cfg)

message("Reading ", length(prefixes), " IDAT pairs ...")
sdfs <- bplapply(unname(prefixes), readIDATpair, BPPARAM = bp)
names(sdfs) <- sheet$sample_id

platform <- if (identical(cfg$platform, "auto")) {
  tryCatch(sesame::sdfPlatform(sdfs[[1]]),
           error = function(e) sesame::inferPlatformFromProbeIDs(sdfs[[1]]$Probe_ID))
} else cfg$platform
writeLines(platform, file.path(d$tab, "platform.txt"))
message("Platform: ", platform)

one_sample <- function(id) {
  sdf   <- sdfs[[id]]
  p_oob <- sesame::pOOBAH(sdf, return.pval = TRUE)
  p_neg <- tryCatch(sesame::detectionPnegEcdf(sdf, return.pval = TRUE), error = function(e) NULL)
  data.frame(sample_id = id, n_probes = length(p_oob),
             mean_detP_negctl     = if (is.null(p_neg)) NA_real_ else mean(p_neg, na.rm = TRUE),
             frac_detected_negctl = if (is.null(p_neg)) NA_real_ else mean(p_neg < 0.05, na.rm = TRUE),
             mean_detP_pOOBAH     = mean(p_oob, na.rm = TRUE),
             frac_detected_pOOBAH = mean(p_oob < 0.05, na.rm = TRUE))
}
qc <- do.call(rbind, bplapply(names(sdfs), one_sample, BPPARAM = bp))
qc <- merge(sheet[, c("sample_id", "label", "patient", "group")], qc, by = "sample_id")

metric <- q$metric %||% "negctl"
if (metric == "negctl" && all(is.na(qc$mean_detP_negctl))) {
  message("WARNING: negative-control detection p-values unavailable; using pOOBAH fraction detected.")
  metric <- "poobah_frac"
}
if (metric == "negctl") {
  qc$mean_detP <- qc$mean_detP_negctl
  qc$pass <- qc$mean_detP <= q$max_mean_detp
  rule <- paste0("mean negative-control detection p <= ", q$max_mean_detp)
} else {
  qc$mean_detP <- qc$mean_detP_pOOBAH
  qc$pass <- qc$frac_detected_pOOBAH >= q$min_frac_detected_poobah
  rule <- paste0("pOOBAH fraction detected >= ", q$min_frac_detected_poobah)
}
write.csv(qc, file.path(d$tab, "qc_detection_pvalues.csv"), row.names = FALSE)
print(qc[, c("sample_id", "label", "mean_detP_negctl", "mean_detP_pOOBAH", "frac_detected_pOOBAH", "pass")])
message(sum(qc$pass), "/", nrow(qc), " samples pass QC (", rule, ")")

# Safety net: never continue silently when most arrays 'fail'
if (isTRUE(q$stop_if_majority_fail %||% TRUE) && sum(qc$pass) < max(4, ceiling(nrow(qc) / 2)))
  stop("Most arrays fail the QC rule (", rule, "). Inspect results/tables/qc_detection_pvalues.csv ",
       "and adjust qc$metric / thresholds in config.yaml deliberately, or set qc$stop_if_majority_fail: false.")

qc$label <- factor(qc$label, levels = qc$label[order(qc$mean_detP)])
p <- ggplot(qc, aes(x = mean_detP, y = label, fill = group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = group_colors(qc$group)) +
  labs(x = if (metric == "negctl") "Mean detection p-value (negative controls)" else "Mean detection p-value (pOOBAH)",
       y = "Sample ID", fill = "Group", title = "Mean detection p-value per array",
       subtitle = paste("QC rule:", rule)) +
  theme_bw(base_size = 12)
if (metric == "negctl")
  p <- p + geom_vline(xintercept = q$max_mean_detp, linetype = "dashed", colour = "red")
save_plot(p, file.path(d$fig, "01_QC_mean_detection_pvalue_barplot"), 7, max(4, 0.32 * nrow(qc) + 1.5))

p2 <- ggplot(qc, aes(x = frac_detected_pOOBAH, y = label, fill = group)) +
  geom_col(width = 0.7) + xlim(0, 1) +
  scale_fill_manual(values = group_colors(qc$group)) +
  labs(x = "Fraction of probes detected (pOOBAH p < 0.05)", y = "Sample ID", fill = "Group",
       title = "Probe detection rate per array") + theme_bw(base_size = 12)
save_plot(p2, file.path(d$fig, "01b_QC_fraction_detected_pOOBAH"), 7, max(4, 0.32 * nrow(qc) + 1.5))
write_session_info(cfg, "02_qc")
