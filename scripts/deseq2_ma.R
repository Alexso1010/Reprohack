# Analyse différentielle (DESeq2) + MA-plots.
## Usage : Rscript deseq2_ma.R <counts.txt> <translation_genes.txt> <annotation.gff> <SRR1:cond,SRR2:cond,...> <dossier_sortie>
suppressMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage : Rscript deseq2_ma.R counts translation gff samples outdir")
counts_file <- args[1]; translation_file <- args[2]; gff_file <- args[3]
samples_arg <- args[4]; outdir <- args[5]

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

# Noms de gènes lus dans le GFF (colonne 9 : "ID=...;Name=...;...")
gff <- read.table(gff_file, sep = "\t", quote = "", comment.char = "#",
                  stringsAsFactors = FALSE, fill = TRUE)
gff_all <- gff 
gff <- gff[gff[, 3] == "gene", ]
get_attr <- function(attrs, key) {
  ifelse(grepl(paste0("(^|;)", key, "="), attrs),
         sub(paste0(".*(^|;)", key, "=([^;]*).*"), "\\2", attrs), NA)
}
gff_id   <- get_attr(gff[, 9], "ID")
gff_name <- get_attr(gff[, 9], "Name")
names_map <- setNames(ifelse(is.na(gff_name), gff_id, gff_name), sub("^gene-", "", gff_id))
gene_name <- unname(names_map[locus])
gene_name[is.na(gene_name)] <- locus[is.na(gene_name)]   # à défaut, on garde le locus tag
names(gene_name) <- locus

# Aminoacyl-tRNA synthétases : gènes dont le produit contient "tRNA ligase" ou "tRNA synthetase"
# (hypothèse à noter dans le rapport : l'article ne dit pas comment ces gènes ont été choisis)
cds <- gff_all[gff_all[, 3] == "CDS", ]
cds_parent  <- sub("^gene-", "", get_attr(cds[, 9], "Parent"))
cds_product <- get_attr(cds[, 9], "product")
aa_trna_ids <- unique(cds_parent[grepl("tRNA (ligase|synthetase)", cds_product, ignore.case = TRUE)])
cat(sprintf("Aminoacyl-tRNA synthétases repérées dans le GFF : %d\n", length(aa_trna_ids)))

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
# Gènes nommés sur la figure de l'article, avec le décalage (x, y) de leur étiquette,
# en unités du graphique. Ces décalages sont à ajuster à l'œil après un premier rendu.
label_offsets <- list(frr = c(-1.5, 1.8), infA = c(0.5, 2.6), infB = c(3.0, -0.8),
                      infC = c(3.0, 0.5),  tsf  = c(1.5, 1.0), pth  = c(-2.5, -0.8))

pdf(file.path(outdir, "figures", "MA_plot_translation.pdf"), width = 5.5, height = 5.5)

x <- log2(tr_df$baseMean)                       # (1) abscisse en log2 base mean
y <- tr_df$log2FoldChange
sig <- !is.na(tr_df$padj) & tr_df$padj < 0.05
fill <- ifelse(sig, "red", "grey50")

plot(x, y, pch = 20, col = fill, xlim = c(0, 20), ylim = c(-6, 5), las = 1,
     xlab = expression(Log[2]~base~Mean), ylab = expression(Log[2]~Fold~Change))
abline(h = 0, lty = 2, lwd = 2)

# (3) contour noir autour des aminoacyl-tRNA synthétases (on redessine ces points, remplis)
aa <- tr_df$locus %in% aa_trna_ids
points(x[aa], y[aa], pch = 21, bg = fill[aa], col = "black", cex = 1.2, lwd = 1.5)

# (2) étiquettes des gènes nommés, reliées à leur point par un trait
for (g in names(label_offsets)) {
  i <- which(tr_df$name == g)
  if (length(i) == 0) { message("Gène à étiqueter absent des résultats : ", g); next }
  i <- i[1]
  tx <- x[i] + label_offsets[[g]][1]
  ty <- y[i] + label_offsets[[g]][2]
  segments(x[i], y[i], tx, ty, lwd = 1.5)
  text(tx, ty, labels = g, font = 3, cex = 1.1)
}

legend("bottomleft", bty = "n", legend = c("Significant", "Non-significant", "AA-tRNA synthetases"),
       pch = c(20, 20, 1), col = c("red", "grey50", "black"), cex = 0.8)
dev.off()

# règle de sélection des synthétases (« produit contenant tRNA ligase ou tRNA synthetase dans le GFF NCBI »)
