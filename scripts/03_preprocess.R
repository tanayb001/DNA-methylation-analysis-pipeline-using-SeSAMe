# STEP 3 - SeSAMe preprocessing, probe filtering, beta-value distribution
source("scripts/utils.R")
suppressPackageStartupMessages({ library(sesame); library(BiocParallel) })

cfg <- load_config(); d <- out_dirs(cfg); pp <- cfg$preprocess
sheet <- read_sample_sheet(cfg)
qc <- read.csv(file.path(d$tab, "qc_detection_pvalues.csv"), stringsAsFactors = FALSE)
failed <- qc$sample_id[!qc$pass]
if (length(failed)) message("Excluding ", length(failed), " sample(s) failing QC: ", paste(failed, collapse = ", "))
sheet <- sheet[sheet$sample_id %in% qc$sample_id[qc$pass], , drop = FALSE]
if (nrow(sheet) < 4) stop("Only ", nrow(sheet), " sample(s) left after QC - check results/tables/qc_detection_pvalues.csv")
prefixes <- resolve_prefixes(sheet, cfg$paths$work_idat_dir)
platform <- read_platform(cfg)
bp <- make_bp(cfg)

steps <- list()
log_step <- function(label, n) steps[[length(steps) + 1]] <<- data.frame(step = label, n_probes = n)

# ---- SeSAMe: mask -> infer channel -> dye-bias -> pOOBAH -> noob -> betas ----
betas <- openSesame(unname(prefixes), prep = pp$sesame_prep, BPPARAM = bp)
if (is.list(betas)) betas <- do.call(cbind, betas)
if (is.null(dim(betas))) betas <- matrix(betas, ncol = 1, dimnames = list(names(betas), NULL))
stopifnot(ncol(betas) == nrow(sheet))
colnames(betas) <- sheet$sample_id          # openSesame keeps input order
log_step("SeSAMe output", nrow(betas))

if (platform == "EPICv2") {                  # EPICv2 has replicate probes (cgXXXX_TC21 ...)
  betas <- sesame::betasCollapseToPfx(betas)
  log_step("EPICv2 replicate probes collapsed to prefix", nrow(betas))
}
if (isTRUE(pp$keep_cg_only)) {
  betas <- betas[grepl("^cg", rownames(betas)), , drop = FALSE]
  log_step("CpG (cg) probes only", nrow(betas))
}
na_frac <- rowMeans(is.na(betas))
betas <- betas[na_frac <= pp$max_na_fraction, , drop = FALSE]
log_step(paste0("NA fraction <= ", pp$max_na_fraction, " (masked/undetected)"), nrow(betas))

gr <- get_manifest(platform, cfg$genome)
betas <- betas[rownames(betas) %in% names(gr), , drop = FALSE]
log_step("Present in manifest", nrow(betas))
if (isTRUE(pp$remove_sex_chr)) {
  chr <- as.character(GenomicRanges::seqnames(gr[rownames(betas)]))
  betas <- betas[!chr %in% c("chrX", "chrY", "X", "Y"), , drop = FALSE]
  log_step("Sex-chromosome probes removed", nrow(betas))
}
write.csv(do.call(rbind, steps), file.path(d$tab, "probe_filtering_summary.csv"), row.names = FALSE)

# ---- Beta distribution (all high-quality arrays) ------------------------------
set.seed(cfg$seed)
idx <- sample(nrow(betas), min(200000, nrow(betas)))
m   <- match(colnames(betas), sheet$sample_id)
long <- data.frame(sample = rep(sheet$label[m], each = length(idx)),
                   group  = rep(sheet$group[m], each = length(idx)),
                   beta   = as.vector(betas[idx, ]))
long <- long[!is.na(long$beta), ]
p <- ggplot(long, aes(beta, colour = group, group = sample)) +
  geom_density(linewidth = 0.6) +
  scale_colour_manual(values = group_colors(sheet$group)) +
  labs(x = "Beta value", y = "Density", colour = "Group",
       title = paste0("Global methylation (beta) distribution - ", ncol(betas), " high-quality arrays"),
       subtitle = paste0(format(nrow(betas), big.mark = ","), " probes after filtering (200k random probes plotted)")) +
  theme_bw(base_size = 12)
save_plot(p, file.path(d$fig, "02_beta_value_distribution"), 7.5, 5)

# ---- Technical replicate concordance (Spearman) --------------------------------
reps <- names(which(table(sheet$patient) > 1))
if (length(reps)) {
  rc <- list(); sc <- list()
  for (pt in reps) {
    ids <- sheet$sample_id[sheet$patient == pt][1:2]
    ok  <- complete.cases(betas[, ids[1]], betas[, ids[2]])
    rc[[pt]] <- data.frame(patient = pt, rep_a = ids[1], rep_b = ids[2], n_probes = sum(ok),
                           spearman = cor(betas[ok, ids[1]], betas[ok, ids[2]], method = "spearman"))
    i <- sample(which(ok), min(30000, sum(ok)))
    sc[[pt]] <- data.frame(patient = pt, x = betas[i, ids[1]], y = betas[i, ids[2]])
  }
  rc <- do.call(rbind, rc)
  write.csv(rc, file.path(d$tab, "replicate_spearman.csv"), row.names = FALSE)
  sc <- do.call(rbind, sc)
  p <- ggplot(sc, aes(x, y)) + geom_point(size = 0.2, alpha = 0.15) +
    facet_wrap(~patient) + coord_equal() +
    labs(x = "Replicate 1 (beta)", y = "Replicate 2 (beta)", title = "Technical replicate concordance") +
    theme_bw()
  save_plot(p, file.path(d$fig, "03_replicate_concordance"), 6, 4.5)
}

# ---- Analysis unit = patient (replicates averaged) or sample -------------------
if (isTRUE(pp$collapse_replicates)) {
  pat  <- sheet$patient[match(colnames(betas), sheet$sample_id)]
  upat <- unique(pat)
  grp_per_pat <- tapply(sheet$group, sheet$patient, function(x) length(unique(x)))
  if (any(grp_per_pat > 1)) stop("A patient has replicates assigned to different groups")
  beta_a <- sapply(upat, function(p) {
    r <- rowMeans(betas[, pat == p, drop = FALSE], na.rm = TRUE); r[is.nan(r)] <- NA; r })
  pheno_a <- sheet[!duplicated(sheet$patient), , drop = FALSE]
  pheno_a <- pheno_a[match(upat, pheno_a$patient), , drop = FALSE]
  rownames(pheno_a) <- pheno_a$patient
  colnames(beta_a) <- pheno_a$patient
  pheno_a$label <- pheno_a$patient
} else {
  beta_a <- betas; pheno_a <- sheet[colnames(betas), , drop = FALSE]
}
message("Analysis units: ", ncol(beta_a), "  (", paste(names(table(pheno_a$group)), table(pheno_a$group), collapse = "; "), ")")

saveRDS(list(betas = betas, pheno = sheet[colnames(betas), ], platform = platform),
        file.path(d$obj, "betas_samples.rds"))
saveRDS(list(betas = beta_a, pheno = pheno_a, platform = platform),
        file.path(d$obj, "betas_analysis.rds"))
write.csv(pheno_a, file.path(d$tab, "pheno_analysis.csv"))
write_session_info(cfg, "03_preprocess")
