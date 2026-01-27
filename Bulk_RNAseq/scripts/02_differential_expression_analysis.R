#!/usr/bin/env Rscript

# ==============================================================================
# Bulk RNA-seq differential expression analysis (SW620 vs SW480) from Salmon quantifications
# - Dual-species quantification performed separately against:
#     Human tumor transcriptome (GRCh38)
#     Zebrafish host/TME transcriptome (GRCz11)
#
# Inputs:
# - Salmon quantification files
#     results/salmon/human/<SAMPLE>/quant.sf
#     results/salmon/zebrafish/<SAMPLE>/quant.sf
# - Ensembl GTF annotation files (for transcript-to-gene mapping and gene_name labels):
#     refs/Homo_sapiens.GRCh38.115.gtf(.gz)
#     (from https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/)
#     refs/"Danio_rerio.GRCz11.115.gtf(.gz)
#     (from https://ftp.ensembl.org/pub/release-115/gtf/danio_rerio/)
#
# Outputs (written to processed/):
#   - Human_allgenes_stats.csv
#   - Human_FULL_SW620_vs_SW480.csv
#   - Human_DEG_FDR0.05_LFC1_SW620_vs_SW480.csv
#   - Human_normalized_expression.csv
#   (same for Zebrafish)
#   - processed_data.RData
#   - sessionInfo.txt
#
# Contrast: SW620 - SW480 (positive = higher in SW620)
# ==============================================================================

suppressPackageStartupMessages({
  library(tximport)
  library(GenomicFeatures)
  library(edgeR)
  library(limma)
  library(here)
  library(rtracklayer)
})

cat("=== Session info (key packages) ===\n")
cat("R:", R.version.string, "\n")
pkgs <- c("tximport","GenomicFeatures","AnnotationDbi","edgeR","limma","here","rtracklayer")
for (p in pkgs) cat(sprintf("%s: %s\n", p, as.character(packageVersion(p))))
cat("===================================\n\n")

# -----------------------------
# User-configurable parameters
# -----------------------------
FDR_CUTOFF <- 0.05
LFC_TREAT  <- 1    # treat threshold in log2 units

PROJECT_ROOT <- here::here()

SALMON_HUMAN_DIR <- file.path(PROJECT_ROOT, "results", "salmon", "human")
SALMON_ZFISH_DIR <- file.path(PROJECT_ROOT, "results", "salmon", "zebrafish")

REF_DIR       <- file.path(PROJECT_ROOT, "refs")
PROCESSED_DIR <- file.path(PROJECT_ROOT, "processed")
CACHE_DIR     <- file.path(PROCESSED_DIR, "cache")

dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR,     showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Reference annotation (GTF)
# -----------------------------
pick_gtf <- function(path_no_gz) {
  if (file.exists(path_no_gz)) return(path_no_gz)
  if (file.exists(paste0(path_no_gz, ".gz"))) return(paste0(path_no_gz, ".gz"))
  stop("Missing GTF: ", path_no_gz, " (.gz also not found)")
}

human_gtf <- pick_gtf(file.path(REF_DIR, "Homo_sapiens.GRCh38.115.gtf"))
zfish_gtf <- pick_gtf(file.path(REF_DIR, "Danio_rerio.GRCz11.115.gtf"))

# -----------------------------
# Caching helpers
# -----------------------------
cache_read <- function(path) if (file.exists(path)) readRDS(path) else NULL
cache_write <- function(obj, path) { saveRDS(obj, path); obj }

tx2gene_cache_h <- file.path(CACHE_DIR, "tx2gene_human.rds")
tx2gene_cache_z <- file.path(CACHE_DIR, "tx2gene_zfish.rds")
txdb_cache_h    <- file.path(CACHE_DIR, "txdb_human.sqlite")
txdb_cache_z    <- file.path(CACHE_DIR, "txdb_zfish.sqlite")
genemap_cache_h <- file.path(CACHE_DIR, "gene_map_human.rds")
genemap_cache_z <- file.path(CACHE_DIR, "gene_map_zfish.rds")

# -----------------------------
# Build/load tx2gene mapping from GTF (via TxDb)
# -----------------------------
cat("Creating/loading transcript-to-gene mappings...\n")

make_or_load_tx2gene <- function(gtf, txdb_sqlite_path, tx2gene_rds_path) {
  tx2gene <- cache_read(tx2gene_rds_path)
  if (!is.null(tx2gene)) return(tx2gene)

  txdb <- NULL
  if (file.exists(txdb_sqlite_path)) {
    txdb <- loadDb(txdb_sqlite_path)
  } else {
    txdb <- txdbmaker::makeTxDbFromGFF(gtf)
    saveDb(txdb, txdb_sqlite_path)
  }

  df <- select(txdb,
               keys = keys(txdb, "GENEID"),
               columns = "TXNAME",
               keytype = "GENEID")
  tx2gene <- df[, c("TXNAME", "GENEID")]
  cache_write(tx2gene, tx2gene_rds_path)
}

tx2gene_human <- make_or_load_tx2gene(human_gtf, txdb_cache_h, tx2gene_cache_h)
tx2gene_zfish <- make_or_load_tx2gene(zfish_gtf, txdb_cache_z, tx2gene_cache_z)

# -----------------------------
# Find Salmon quant.sf files
# -----------------------------
cat("Importing Salmon quantifications...\n")

get_quant_files <- function(parent_dir) {
  if (!dir.exists(parent_dir)) stop("Missing Salmon directory: ", parent_dir)

  sample_dirs <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
  files <- file.path(sample_dirs, "quant.sf")
  ok <- file.exists(files)
  files <- files[ok]
  names(files) <- basename(dirname(files))
  files
}

files_human <- get_quant_files(SALMON_HUMAN_DIR)
files_zfish <- get_quant_files(SALMON_ZFISH_DIR)

if (length(files_human) == 0) stop("No human quant.sf found in: ", SALMON_HUMAN_DIR)
if (length(files_zfish) == 0) stop("No zebrafish quant.sf found in: ", SALMON_ZFISH_DIR)

samples <- intersect(names(files_human), names(files_zfish))
if (length(samples) == 0) {
  stop("No matching sample folders between human and zebrafish Salmon outputs.\n",
       "Human samples: ", paste(names(files_human), collapse = ", "), "\n",
       "Zfish samples: ", paste(names(files_zfish), collapse = ", "))
}

samples <- sort(samples)
files_human <- files_human[samples]
files_zfish <- files_zfish[samples]

cat("Samples detected:\n")
print(samples)

# -----------------------------
# tximport (gene-level counts from lengthScaledTPM)
# -----------------------------
txi_human <- tximport(files_human, type = "salmon", tx2gene = tx2gene_human,
                      countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)

txi_zfish <- tximport(files_zfish, type = "salmon", tx2gene = tx2gene_zfish,
                      countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)

counts_human <- txi_human$counts
counts_zfish <- txi_zfish$counts

# -----------------------------
# Gene name mapping from GTF (gene_id -> gene_name)
# -----------------------------
gtf_gene_map <- function(gtf_file) {
  gtf <- rtracklayer::import(gtf_file)
  gtf <- gtf[gtf$type == "gene"]
  data.frame(
    ensembl_gene_id = as.character(mcols(gtf)$gene_id),
    external_gene_name = as.character(mcols(gtf)$gene_name),
    stringsAsFactors = FALSE
  )
}

cat("Creating/loading gene name maps...\n")

gene_info_human <- cache_read(genemap_cache_h)
if (is.null(gene_info_human)) gene_info_human <- cache_write(gtf_gene_map(human_gtf), genemap_cache_h)

gene_info_zfish <- cache_read(genemap_cache_z)
if (is.null(gene_info_zfish)) gene_info_zfish <- cache_write(gtf_gene_map(zfish_gtf), genemap_cache_z)

make_gene_names <- function(count_matrix, gene_info_df) {
  m <- match(rownames(count_matrix), gene_info_df$ensembl_gene_id)
  gene_names <- gene_info_df$external_gene_name[m]
  gene_names[is.na(gene_names) | gene_names == ""] <- rownames(count_matrix)[is.na(gene_names) | gene_names == ""]
  make.unique(gene_names, sep = "_")
}

gene_names_human <- make_gene_names(counts_human, gene_info_human)
gene_names_zfish <- make_gene_names(counts_zfish, gene_info_zfish)

# -----------------------------
# Groups + design matrix
# -----------------------------
group <- ifelse(grepl("SW480", samples), "SW480",
                ifelse(grepl("SW620", samples), "SW620", NA))

if (any(is.na(group))) {
  stop("Could not infer group for: ", paste(samples[is.na(group)], collapse = ", "),
       "\nRename samples or set group manually.")
}

group <- factor(group, levels = c("SW480", "SW620"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
contrast <- makeContrasts(SW620_vs_SW480 = SW620 - SW480, levels = design)

# -----------------------------
# Differential expression helper
# -----------------------------
run_de <- function(counts_mat, gene_names, label) {
  cat("\n-----------------------------------\n")
  cat("Processing", label, "differential expression...\n")
  cat("-----------------------------------\n")

  y <- DGEList(counts = counts_mat, group = group)
  rownames(y$counts) <- gene_names

  keep <- filterByExpr(y, group = group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, contrast)

  # 1) "All-genes stats" for GSEA ranking (limma eBayes t-statistics)
  fit_eb <- eBayes(fit)
  all_stats <- topTable(fit_eb, coef = 1, n = Inf, sort.by = "none")
  all_stats$gene <- rownames(all_stats)

  write.csv(all_stats,
            file.path(PROCESSED_DIR, paste0(label, "_allgenes_stats.csv")),
            row.names = FALSE)

  # 2) DE testing with treat threshold
  tfit <- treat(fit, lfc = LFC_TREAT)
  full <- topTreat(tfit, coef = 1, n = Inf)  # all tested genes, sorted

  write.csv(full,
            file.path(PROCESSED_DIR, paste0(label, "_FULL_SW620_vs_SW480.csv")),
            row.names = TRUE)

  deg <- full[full$adj.P.Val < FDR_CUTOFF & abs(full$logFC) > LFC_TREAT, , drop = FALSE]
  write.csv(deg,
            file.path(PROCESSED_DIR, paste0(label, "_DEG_FDR", FDR_CUTOFF, "_LFC", LFC_TREAT, "_SW620_vs_SW480.csv")),
            row.names = TRUE)

  cat(label, ": tested genes after filtering =", nrow(full), "\n")
  cat(label, ": DEGs (FDR <", FDR_CUTOFF, "and |logFC| >", LFC_TREAT, ") =", nrow(deg), "\n")

  list(y = y, v = v, fit = fit, fit_eb = fit_eb, tfit = tfit,
       full = full, deg = deg, all_stats = all_stats)
}

# -----------------------------
# Run DE for Human + Zebrafish
# -----------------------------
res_human <- run_de(counts_human, gene_names_human, "Human")
res_zfish <- run_de(counts_zfish, gene_names_zfish, "Zebrafish")

# -----------------------------
# Save processed data
# -----------------------------
write.csv(res_human$v$E, file.path(PROCESSED_DIR, "Human_normalized_expression.csv"))
write.csv(res_zfish$v$E, file.path(PROCESSED_DIR, "Zebrafish_normalized_expression.csv"))

gene_names_human_filtered <- rownames(res_human$v$E)
gene_names_zfish_filtered <- rownames(res_zfish$v$E)

save(counts_human, counts_zfish,
     gene_names_human, gene_names_zfish,
     gene_names_human_filtered, gene_names_zfish_filtered,
     group, samples,
     txi_human, txi_zfish,
     res_human, res_zfish,
     file = file.path(PROCESSED_DIR, "processed_data.RData"))

cat("\nAnalysis complete. Results saved in:", PROCESSED_DIR, "\n")

# -----------------------------
# Save R session information
# -----------------------------
sink(file.path(PROCESSED_DIR, "sessionInfo.txt"))
sessionInfo()
sink()