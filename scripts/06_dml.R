# STEP 6 - Differential methylation with SeSAMe DML (linear model per probe)
source("scripts/utils.R")
suppressPackageStartupMessages({ library(sesame); library(SummarizedExperiment); library(BiocParallel) })

cfg <- load_config(); d <- out_dirs(cfg)
obj <- readRDS(file.path(d$obj, "betas_analysis.rds"))
case <- cfg$contrast$case; ref <- cfg$contrast$reference
ph <- obj$pheno
keep <- ph$group %in% c(case, ref)
ph <- ph[keep, , drop = FALSE]; betas <- obj$betas[, rownames(ph), drop = FALSE]
tb <- table(ph$group); print(tb)
if (any(tb[c(case, ref)] < 2, na.rm = TRUE) || length(tb) < 2)
  stop("Need >= 2 independent units in each group (case='", case, "', reference='", ref, "')")
ph$group <- factor(ph$group, levels = c(ref, case))

covs <- unlist(cfg$covariates)
for (cv in covs) { stopifnot(cv %in% names(ph)); if (is.character(ph[[cv]])) ph[[cv]] <- factor(ph[[cv]]) }
fm <- as.formula(paste("~", paste(c("group", covs), collapse = " + ")))
message("Model: ", deparse(fm), "   (", case, " vs ", ref, ")")

se <- SummarizedExperiment(assays = list(betas = betas),
                           colData = S4Vectors::DataFrame(ph[, c("group", covs), drop = FALSE]))
res  <- DML(se, fm, BPPARAM = make_bp(cfg))
smry <- summaryExtractTest(res)
if (is.null(smry$Probe_ID)) smry$Probe_ID <- rownames(smry)
write_gz(smry, file.path(d$tab, "dml_sesame_raw.csv.gz"))

est_col <- grep("^Est_group", names(smry), value = TRUE)
p_col   <- grep("^Pval_group", names(smry), value = TRUE)
stopifnot(length(est_col) == 1, length(p_col) == 1)

dmp <- data.frame(Probe_ID = smry$Probe_ID, delta_beta = smry[[est_col]], P = smry[[p_col]],
                  stringsAsFactors = FALSE)
dmp$FDR <- p.adjust(dmp$P, "BH")
in_case <- ph$group == case
dmp$mean_beta_case <- rowMeans(betas[dmp$Probe_ID, in_case, drop = FALSE], na.rm = TRUE)
dmp$mean_beta_ref  <- rowMeans(betas[dmp$Probe_ID, !in_case, drop = FALSE], na.rm = TRUE)
chk <- cor(dmp$delta_beta, dmp$mean_beta_case - dmp$mean_beta_ref, use = "complete.obs")
message(sprintf("Sanity check: cor(model estimate, raw mean difference) = %.3f (should be ~1; sign = %s minus %s)",
                chk, case, ref))
dmp <- dmp[!is.na(dmp$P), ]
dmp <- flag_significant(dmp, cfg)
dmp <- dmp[order(dmp$P), ]
dmp$contrast <- paste(case, "vs", ref)

write_gz(dmp, file.path(d$tab, "dmp_all.csv.gz"))
summ <- as.data.frame(table(dmp$Direction)); names(summ) <- c("Direction", "n_probes")
summ$criteria <- dmp$sig_basis[1]
write.csv(summ, file.path(d$tab, "dmp_summary.csv"), row.names = FALSE)
print(summ)
write_session_info(cfg, "06_dml")
