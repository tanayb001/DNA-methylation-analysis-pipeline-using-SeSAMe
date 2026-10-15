# STEP 7 - Probe annotation and volcano plots
# Uses existing dmp_all.csv.gz; does not rerun preprocessing or DML.

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(GenomeInfoDb)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

cfg <- load_config()
d <- out_dirs(cfg)
an <- cfg$annotation

if (!identical(cfg$genome, "hg38")) {
  stop("This annotation script is configured for hg38.")
}

dmp <- read.csv(
  file.path(d$tab, "dmp_all.csv.gz"),
  stringsAsFactors = FALSE
)

required <- c("Probe_ID", "delta_beta", "P", "FDR", "Direction", "Significant")
missing <- setdiff(required, names(dmp))
if (length(missing)) {
  stop("Missing required columns in dmp_all.csv.gz: ",
       paste(missing, collapse = ", "))
}

# ---- Probe coordinates and hg38 gene models ----
platform <- read_platform(cfg)
manifest <- get_manifest(platform, cfg$genome)

if (is.null(names(manifest)) || anyDuplicated(names(manifest))) {
  stop("The probe manifest must have unique probe IDs as names.")
}

idx <- match(dmp$Probe_ID, names(manifest))
has_coord <- !is.na(idx)
gr <- manifest[idx[has_coord]]

GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")

# Re-match after chromosome filtering so excluded probes remain in the output,
# with missing coordinate/annotation fields.
coord_ids <- names(gr)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_gr <- suppressMessages(genes(txdb))
gene_gr <- GenomeInfoDb::keepStandardChromosomes(
  gene_gr, pruning.mode = "coarse"
)

gene_gr$symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = names(gene_gr),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

valid_gene <- !is.na(gene_gr$symbol) & nzchar(gene_gr$symbol)
gene_gr <- gene_gr[valid_gene]

prom_gr <- promoters(
  gene_gr,
  upstream = an$promoter_upstream,
  downstream = an$promoter_downstream
)
prom_gr$symbol <- gene_gr$symbol

# Set up a one-row-per-probe annotation table. If a probe has multiple
# qualifying hits, the selection order is: promoter, gene body, nearest gene.
anno <- data.frame(
  Probe_ID = dmp$Probe_ID,
  chr = NA_character_,
  pos = NA_integer_,
  gene = NA_character_,
  region = "Intergenic",
  dist_to_gene_bp = NA_integer_,
  stringsAsFactors = FALSE
)

if (length(gr)) {
  gr_ids <- names(gr)
  m <- match(gr_ids, anno$Probe_ID)
  anno$chr[m] <- as.character(seqnames(gr))
  anno$pos[m] <- start(gr)

  # Nearest-gene assignment is limited to the configured distance.
  nearest <- distanceToNearest(gr, gene_gr, ignore.strand = TRUE)
  if (length(nearest)) {
    q <- queryHits(nearest)
    s <- subjectHits(nearest)
    dd <- mcols(nearest)$distance
    within_limit <- dd <= an$max_nearest_distance
    q <- q[within_limit]
    s <- s[within_limit]
    dd <- dd[within_limit]

    if (length(q)) {
      nearest_rows <- data.frame(
        probe = gr_ids[q],
        gene = gene_gr$symbol[s],
        distance = as.integer(dd),
        stringsAsFactors = FALSE
      )
      nearest_rows <- nearest_rows[
        order(nearest_rows$probe, nearest_rows$distance,
              nearest_rows$gene),
      ]
      nearest_rows <- nearest_rows[!duplicated(nearest_rows$probe), ]
      k <- match(nearest_rows$probe, anno$Probe_ID)
      anno$gene[k] <- nearest_rows$gene
      anno$dist_to_gene_bp[k] <- nearest_rows$distance
    }
  }

  # Assign gene-body hits, deterministically choosing the alphabetically
  # first symbol if a probe overlaps multiple gene bodies.
  body_hits <- findOverlaps(gr, gene_gr, ignore.strand = TRUE)
  if (length(body_hits)) {
    body <- data.frame(
      probe = gr_ids[queryHits(body_hits)],
      gene = gene_gr$symbol[subjectHits(body_hits)],
      stringsAsFactors = FALSE
    )
    body <- body[order(body$probe, body$gene), ]
    body <- body[!duplicated(body$probe), ]
    k <- match(body$probe, anno$Probe_ID)
    anno$gene[k] <- body$gene
    anno$region[k] <- "Gene body"
    anno$dist_to_gene_bp[k] <- 0L
  }

  # Promoter assignment has highest priority, including when promoters
  # overlap one another.
  promoter_hits <- findOverlaps(gr, prom_gr, ignore.strand = TRUE)
  if (length(promoter_hits)) {
    promoter <- data.frame(
      probe = gr_ids[queryHits(promoter_hits)],
      gene = prom_gr$symbol[subjectHits(promoter_hits)],
      stringsAsFactors = FALSE
    )
    promoter <- promoter[order(promoter$probe, promoter$gene), ]
    promoter <- promoter[!duplicated(promoter$probe), ]
    k <- match(promoter$probe, anno$Probe_ID)
    anno$gene[k] <- promoter$gene
    anno$region[k] <- "Promoter"
    anno$dist_to_gene_bp[k] <- 0L
  }
}

ann <- left_join(dmp, anno, by = "Probe_ID")
write_gz(ann, file.path(d$tab, "dmp_annotated.csv.gz"))
write.csv(
  ann[ann$Significant %in% TRUE, ],
  file.path(d$tab, "dmp_significant_annotated.csv"),
  row.names = FALSE
)

message(
  "Annotated ", sum(!is.na(ann$gene)), " / ", nrow(ann),
  " probes with a gene; ",
  sum(ann$Significant %in% TRUE), " significant"
)

# ---- Gene-level descriptive summary ----
# Each assigned probe contributes once to its selected gene. The minimum-P
# adjustment within each gene accounts for the number of assigned probes.
gl <- ann %>%
  filter(!is.na(gene), is.finite(P), !is.na(delta_beta)) %>%
  group_by(gene) %>%
  summarise(
    n_probes = n(),
    n_sig_probes = sum(Significant %in% TRUE),
    mean_delta_beta = mean(delta_beta),
    min_P = min(P),
    .groups = "drop"
  ) %>%
  mutate(
    P_gene = pmin(1, min_P * n_probes),
    FDR_gene = p.adjust(P_gene, method = "BH"),
    Significant = FDR_gene < 0.05,
    Direction = case_when(
      mean_delta_beta > 0 ~ "Hyper",
      mean_delta_beta < 0 ~ "Hypo",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(P_gene)

write.csv(
  gl,
  file.path(d$tab, "gene_level_summary.csv"),
  row.names = FALSE
)

# ---- Volcano plots ----
dbc <- cfg$dm$delta_beta_cutoff
n_labels <- cfg$dm$n_labels

make_volcano <- function(df, x_col, p_col, label_col, x_label,
                         title, subtitle, caption = NULL) {
  df <- df[is.finite(df[[x_col]]) & is.finite(df[[p_col]]), , drop = FALSE]
  if (!nrow(df)) stop("No finite data available for plot: ", title)

  df$plot_y <- -log10(pmax(df[[p_col]], .Machine$double.xmin))
  df$plot_direction <- factor(
    df$Direction,
    levels = c("Hyper", "Hypo", "NS")
  )

  # Thin only the non-significant points for display; all data remain in files.
  ns <- which(df$plot_direction == "NS")
  if (length(ns) > 150000) {
    set.seed(1)
    keep_ns <- sample(ns, 150000)
    df <- df[sort(c(which(df$plot_direction != "NS"), keep_ns)), ,
             drop = FALSE]
  }

  labeled <- df[df$plot_direction != "NS", , drop = FALSE]
  labeled <- labeled[order(labeled[[p_col]]), , drop = FALSE]
  labeled <- head(labeled, n_labels)

  ggplot(df, aes(x = .data[[x_col]], y = plot_y,
                 colour = plot_direction)) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_colour_manual(
      values = c(Hyper = "#D55E00", Hypo = "#0072B2", NS = "grey75"),
      drop = FALSE
    ) +
    geom_vline(
      xintercept = c(-dbc, dbc),
      linetype = "dashed",
      colour = "grey40"
    ) +
    geom_text_repel(
      data = labeled,
      aes(label = .data[[label_col]]),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    labs(
      x = x_label,
      y = expression(-log[10](P)),
      title = title,
      subtitle = subtitle,
      caption = caption,
      colour = "Direction"
    ) +
    theme_bw(base_size = 12)
}

basis <- if ("sig_basis" %in% names(ann)) ann$sig_basis[1] else ""
probe_counts <- table(factor(ann$Direction, levels = c("Hyper", "Hypo")))
ann$lab <- ifelse(
  is.na(ann$gene),
  ann$Probe_ID,
  paste0(ann$gene, " (", ann$Probe_ID, ")")
)

p_probe <- make_volcano(
  ann,
  x_col = "delta_beta",
  p_col = "P",
  label_col = "lab",
  x_label = paste0(
    "SeSAMe model estimate: delta beta (",
    cfg$contrast$case, " - ", cfg$contrast$reference, ")"
  ),
  title = "Differentially methylated probes",
  subtitle = paste0(
    "Hyper: ", probe_counts["Hyper"],
    " | Hypo: ", probe_counts["Hypo"],
    if (nzchar(basis)) paste0("\n", basis) else ""
  ),
  caption = "Non-significant probes are randomly down-sampled for display only."
)

save_plot(p_probe, file.path(d$fig, "07_volcano_probes"), 7.5, 6)

gene_counts <- table(factor(gl$Direction, levels = c("Hyper", "Hypo")))
p_gene <- make_volcano(
  gl,
  x_col = "mean_delta_beta",
  p_col = "P_gene",
  label_col = "gene",
  x_label = "Mean assigned-probe delta beta per gene",
  title = "Gene-level summary volcano plot",
  subtitle = paste0(
    "Hyper: ", gene_counts["Hyper"],
    " | Hypo: ", gene_counts["Hypo"],
    "\nP_gene = min probe P multiplied by assigned probe count; BH-adjusted across genes"
  )
)

save_plot(p_gene, file.path(d$fig, "08_volcano_genes"), 7.5, 6)

write_session_info(cfg, "07_annotate_volcano")
