# Shared helpers, sourced by every R step -------------------------------------
suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

args_config <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) >= 1) a[1] else "config/config.yaml"
}

load_config <- function(path = args_config()) {
  stopifnot(file.exists(path))
  cfg <- yaml::read_yaml(path)
  set.seed(cfg$seed %||% 42)
  cfg
}

# Parallel backend: serial by default (robust on shared servers with CPU-time limits)
make_bp <- function(cfg) {
  n <- as.integer(cfg$n_cores %||% 1)
  if (n <= 1) BiocParallel::SerialParam(progressbar = FALSE)
  else BiocParallel::MulticoreParam(workers = n, timeout = 86400, stop.on.error = TRUE)
}

ensure_dir <- function(...) {
  d <- file.path(...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

out_dirs <- function(cfg) {
  r <- cfg$paths$results_dir
  list(fig = ensure_dir(r, "figures"), tab = ensure_dir(r, "tables"),
       obj = ensure_dir(r, "objects"), log = ensure_dir(r, "logs"))
}

read_sample_sheet <- function(cfg) {
  sh <- read.csv(cfg$paths$sample_sheet, stringsAsFactors = FALSE, check.names = FALSE)
  sc <- cfg$sample_cols
  miss <- setdiff(c(sc$sample_id, sc$patient, sc$group), names(sh))
  if (length(miss)) stop("Sample sheet is missing column(s): ", paste(miss, collapse = ", "))
  sh$sample_id <- as.character(sh[[sc$sample_id]])
  sh$patient   <- as.character(sh[[sc$patient]])
  sh$group     <- as.character(sh[[sc$group]])
  if (!"label" %in% names(sh)) sh$label <- sh$sample_id
  sh <- sh[!is.na(sh$group) & nzchar(sh$group), , drop = FALSE]
  if (anyDuplicated(sh$sample_id)) stop("Duplicated sample_id in sample sheet")
  rownames(sh) <- sh$sample_id
  sh
}

# Map each sample_id to exactly one IDAT prefix (path without _Grn.idat / _Red.idat)
resolve_prefixes <- function(sheet, idat_dir) {
  files <- list.files(idat_dir, pattern = "_Grn\\.idat(\\.gz)?$")
  pref  <- unique(sub("_Grn\\.idat(\\.gz)?$", "", files))
  out <- vapply(sheet$sample_id, function(id) {
    hit <- pref[pref == id]
    if (!length(hit)) hit <- pref[startsWith(pref, id)]
    if (length(hit) != 1)
      stop("sample_id '", id, "' matches ", length(hit), " IDAT prefixes (need exactly 1)")
    if (!file.exists(file.path(idat_dir, paste0(hit, "_Red.idat"))) &&
        !file.exists(file.path(idat_dir, paste0(hit, "_Red.idat.gz"))))
      stop("Red channel file missing for ", hit)
    file.path(idat_dir, hit)
  }, character(1))
  names(out) <- sheet$sample_id
  out
}

read_platform <- function(cfg) {
  if (!identical(cfg$platform, "auto")) return(cfg$platform)
  readLines(file.path(cfg$paths$results_dir, "tables", "platform.txt"), n = 1)
}

save_plot <- function(p, path_no_ext, width = 7, height = 5) {
  ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height)
  ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 300)
}

group_colors <- function(groups) {
  g <- unique(as.character(groups))
  pal <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
  setNames(rep(pal, length.out = length(g)), g)
}

write_gz <- function(df, path) write.csv(df, gzfile(path), row.names = FALSE)

# SeSAMe probe manifest as GRanges named by the 10-char CpG prefix (cg########)
get_manifest <- function(platform, genome = "hg38") {
  gr <- sesameData::sesameData_getManifestGRanges(platform, genome = genome)
  nm <- names(gr)
  if (is.null(nm)) nm <- as.character(S4Vectors::mcols(gr)$Probe_ID)
  names(gr) <- sub("_.*$", "", nm)
  gr[!duplicated(names(gr))]
}

# Call significant probes according to config (with transparent fallback)
flag_significant <- function(dmp, cfg) {
  dm   <- cfg$dm
  pcol <- if (toupper(dm$p_column) == "FDR") "FDR" else "P"
  sig  <- !is.na(dmp[[pcol]]) & dmp[[pcol]] < dm$p_cutoff &
          abs(dmp$delta_beta) >= dm$delta_beta_cutoff
  basis <- paste0(pcol, " < ", dm$p_cutoff, " & |delta beta| >= ", dm$delta_beta_cutoff)
  if (sum(sig) == 0 && isTRUE(dm$fallback_to_nominal) && pcol == "FDR") {
    sig <- !is.na(dmp$P) & dmp$P < dm$p_cutoff & abs(dmp$delta_beta) >= dm$delta_beta_cutoff
    basis <- paste0("EXPLORATORY: nominal P < ", dm$p_cutoff, " & |delta beta| >= ",
                    dm$delta_beta_cutoff, " (no probe passed FDR)")
    message("WARNING: no probe passed FDR; using NOMINAL P (exploratory only).")
  }
  dmp$Significant <- sig
  dmp$sig_basis   <- basis
  dmp$Direction   <- ifelse(sig, ifelse(dmp$delta_beta > 0, "Hyper", "Hypo"), "NS")
  dmp
}

# Package versions only (no paths / host names) -> safe to publish
write_session_info <- function(cfg, step) {
  d  <- out_dirs(cfg)$log
  si <- sessionInfo()
  pk <- c(si$otherPkgs, si$loadedOnly)
  df <- data.frame(package = names(pk),
                   version = vapply(pk, function(x) as.character(x$Version), ""),
                   row.names = NULL)
  df <- df[order(df$package), ]
  f  <- file.path(d, paste0("R_packages_", step, ".txt"))
  writeLines(c(R.version.string, paste("step:", step), "", "package\tversion",
               paste(df$package, df$version, sep = "\t")), f)
}
