# Analyse différentielle (DESeq2) + MA-plots.
# Usage : Rscript deseq2_ma.R <counts.txt> <translation_genes.txt> <SRR1:cond,SRR2:cond,...> <dossier_sortie>
#
# - counts.txt : sortie de featureCounts (--extraAttributes Name)
# - contraste  : persister vs control (log2FC > 0 = plus exprimé chez les persisters)
# - seuil      : padj < 0.05 (Benjamini-Hochberg, comme dans l'article)
suppressMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage : Rscript deseq2_ma.R counts translation samples outdir")
counts_file <- args[1]; translation_file <- args[2]
samples_arg <- args[3]; outdir <- args[4]

dir.create(file.path(outdir, "deseq2"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "figures"), showWarnings = FALSE, recursive = TRUE)

# --- 1. Table des échantillons (SRR -> condition)
pairs <- strsplit(strsplit(samples_arg, ",")[[1]], ":")
coldata <- data.frame(sample = sapply(pairs, `[`, 1),
                      condition = sapply(pairs, `[`, 2),
                      stringsAsFactors = FALSE)
rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition, levels = c("control", "persister"))
if (any(is.na(coldata$condition)))
  stop("Les conditions doivent être 'control' ou 'persister' (voir config.yaml)")

# --- 2. Lecture des comptages featureCounts
raw <- read.table(counts_file, header = TRUE, sep = "\t", comment.char = "#",
                  check.names = FALSE, stringsAsFactors = FALSE)
# Les colonnes de comptage sont nommées d'après les chemins des BAM : results/mapping/SRRxxx.bam
colnames(raw) <- sub("\\.bam$", "", sub(".*/", "", colnames(raw)))
missing <- setdiff(coldata$sample, colnames(raw))
if (length(missing) > 0) stop(paste("Échantillons absents de counts.txt :", paste(missing, collapse = ", ")))

counts <- as.matrix(raw[, coldata$sample, drop = FALSE])
# Identifiant du gène : on enlève un éventuel préfixe "gene-" (GFF NCBI) pour obtenir le locus tag
locus <- sub("^gene-", "", raw$Geneid)
rownames(counts) <- locus
gene_name <- if ("Name" %in% colnames(raw)) raw$Name else locus
names(gene_name) <- locus

# --- 3. DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "persister", "control"), alpha = 0.05)
res_df <- as.data.frame(res)
res_df$locus <- rownames(res_df)
res_df$name  <- gene_name[res_df$locus]
write.csv(res_df, file.path(outdir, "deseq2", "deseq2_results.csv"), row.names = FALSE)

n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
cat(sprintf("Gènes testés : %d ; différentiellement exprimés (padj < 0.05) : %d (up : %d, down : %d)\n",
            sum(!is.na(res_df$padj)), n_sig,
            sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE),
            sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)))

# --- 4. MA-plots
plot_df <- res_df[!is.na(res_df$log2FoldChange) & res_df$baseMean > 0, ]

ma_plot <- function(df, title, label = FALSE) {
  sig <- !is.na(df$padj) & df$padj < 0.05
  plot(df$baseMean, df$log2FoldChange, log = "x", pch = 20,
       col = ifelse(sig, "red", "grey60"),
       xlab = "Mean of normalized counts", ylab = "log2 fold change (persister / control)",
       main = title)
  abline(h = 0, lty = 2)
  if (label && nrow(df) > 0)
    text(df$baseMean, df$log2FoldChange, labels = df$name, pos = 3, cex = 0.5)
}

# Supp. Fig. 3 : jeu de données complet
png(file.path(outdir, "figures", "MA_plot_all.png"), width = 1600, height = 1200, res = 200)
ma_plot(plot_df, "MA-plot - all genes")
dev.off()

# Fig. 3c : gènes liés à la traduction (liste KEGG)
translation <- readLines(translation_file)
tr_df <- plot_df[plot_df$locus %in% translation, ]
cat(sprintf("Gènes 'traduction' retrouvés dans les résultats : %d / %d\n",
            nrow(tr_df), length(translation)))
png(file.path(outdir, "figures", "MA_plot_translation.png"), width = 1600, height = 1200, res = 200)
ma_plot(tr_df, "MA-plot - translation-related genes", label = TRUE)
dev.off()
