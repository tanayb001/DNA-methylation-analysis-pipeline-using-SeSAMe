# Install all R dependencies (run once):  Rscript scripts/00_install_dependencies.R
# Needs R >= 4.4 and Bioconductor >= 3.19 (first release with EPICv2 support in SeSAMe).
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran <- c("yaml", "ggplot2", "ggrepel", "pheatmap", "dplyr", "msigdbr")
bioc <- c("sesame", "sesameData", "SummarizedExperiment", "S4Vectors", "BiocParallel",
          "GenomicRanges", "GenomicFeatures", "GenomeInfoDb", "AnnotationDbi",
          "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db",
          "clusterProfiler", "enrichplot")

need <- c(cran, bioc)[!vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) BiocManager::install(need, ask = FALSE, update = FALSE)

# SeSAMe downloads its probe annotations once into the ExperimentHub cache
sesameData::sesameDataCache()
message("All dependencies installed.")
