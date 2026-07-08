library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(GenomicRanges)
library(SummarizedExperiment)
library(rtracklayer)
library(BiocParallel)
library(DiffBind)
library(DESeq2)
library(ggrepel)
library(data.table)
library(AnnotationDbi)
library(AnnotationHub)
library(dplyr)
library(ggplot2)
library(circlize)


### ---- Differential Expression Analysis ---- ###
# Replace this path with the directory containing the final sorted BAM files produced by preprocessing.sh. 
# This script assumes preprocessing.sh has completed successfully and generated a ${PREFIX}.sorted.bam file for each sample.
rnaseq_dir <- "/path/to/output/final/"

rnaseq_list <- list.files(
  rnaseq_dir,
  pattern = "*.sorted.bam$",
  full.names = TRUE
)
bams <- BamFileList(rnaseq_list)

ah <- AnnotationHub()
edb <- ah[[names(query(ah, "EnsDb.Mmusculus.v102"))]]

txdb <- makeTxDbFromGFF("Mus_musculus.GRCm38.102.gtf.gz", format = "gtf")
autosomal_chr <- as.character(1:19)
sex_chr <- c("X", "Y")
mitochondrial <- "MT"
desired_chr <- c(autosomal_chr, sex_chr, mitochondrial)
txdb_filtered <- keepSeqlevels(txdb, desired_chr, pruning.mode = "coarse")
ebg <- exonsBy(txdb_filtered, by = "gene")
seqlevelsStyle(ebg) <- "UCSC"

rnaseq <- summarizeOverlaps(
    features = ebg, reads = bams,
    mode = "Union",
    singleEnd = FALSE,
    ignore.strand = FALSE,
    fragments = TRUE,
    BPPARAM = MulticoreParam(workers = 16)
)

coldata <- data.frame(
  File = c("S1_ctrl.bam", "S2_ctrl.bam", "S3_ctrl.bam",
           "S4_NR.bam", "S5_NR.bam"),
  Name = c("S1_ctrl", "S2_ctrl", "S3_ctrl",
           "S4_NR", "S5_NR"),
  Treatment = c("Control", "Control", "Control",
                "NR", "NR"),
  stringsAsFactors = FALSE
)
colData(rnaseq) <- DataFrame(coldata)
rnaseq$Treatment <- relevel(factor(rnaseq$Treatment), "Control")
colnames(rnaseq) <- coldata$Name

# Function to process DESeq2 results and return the ordered results
process_deseq_results <- function(dds, contrast, edb, output_file = NULL) {
  res <- results(dds, contrast = contrast)
  res_ordered <- res[order(res$pvalue), ]
  
  res_df <- as.data.frame(res_ordered)
  res_df$ensembl <- rownames(res_df)
  
  annot <- AnnotationDbi::select(edb, 
                                 keys = res_df$ensembl, 
                                 keytype = "GENEID", 
                                 columns = c("GENENAME", "ENTREZID"))
  res_df <- merge(res_df, annot, by.x = "ensembl", by.y = "GENEID", all.x = TRUE)
  
  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts$ensembl <- rownames(norm_counts)
  res_df <- merge(res_df, norm_counts, by = "ensembl", all.x = TRUE)
  
  res_df$negLogPValue <- -log10(res_df$pvalue)
  res_df$negLogPAdjValue <- -log10(res_df$padj)
  
  res_df <- res_df %>%
    relocate(ensembl, GENENAME, ENTREZID, .before = baseMean) %>%
    rename(symbol = GENENAME, entrez = ENTREZID)
  
  if (!is.null(output_file)) {
    fwrite(res_df, output_file, sep = "\t")
  }
  
  return(res_df)
}

rnaseq <- rnaseq[rowSums(assay(rnaseq)) >= 5, ]
dss <- DESeqDataSet(rnaseq, design = ~Treatment)
dss <- DESeq(dss)
NR_ctrl_res_ordered <- process_deseq_results(dss, c("Treatment", "NR", "Control"), edb, "RNAseq.tsv")
vsd <- vst(dss)

### ---- Figure 2f ---- ###

fig2f_df <- NR_ctrl_res_ordered %>%
  mutate(
    log2fc = log2FoldChange,
    neglog10p = -log10(pvalue),
    regulation = case_when(
      pvalue < 0.05 & log2fc >=  0.5 ~ "Up in NR (S4/S5)",
      pvalue < 0.05 & log2fc <= -0.5 ~ "Down in NR (S4/S5)",
      TRUE ~ "Not significant"
    )
  )

fig2f <- ggplot(fig2f_df, aes(x = log2fc, y = neglog10p, color = regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c(
      "Up in NR (S4/S5)"   = "red",
      "Down in NR (S4/S5)" = "blue",
      "Not significant"    = "grey70"
    )
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    x = "log2 fold change (NR vs Ctrl)",
    y = "-log10(p-value)",
    color = "Group",
    title = "Volcano plot: NR vs Ctrl"
  ) +
  theme_bw()

### ---- Figure 2g ---- ###

genes <- c(
  "Neo1","Abca4","Rorb","Vwf","Thbd","Rps6ka3",
  "Plxnb2","Klhl4","Itgb3","Matn4","Muc1","Il15",
  "Socs3","Lsp1","Ndufaf6","Ngp","Ndufaf4"
)

mat_df <- NR_ctrl_res_ordered %>%
  filter(symbol %in% genes) %>%
  select(symbol, S1_ctrl, S2_ctrl, S3_ctrl, S4_NR, S5_NR) %>%
  column_to_rownames("symbol")

mat <- as.matrix(mat_df)
mat <- t(apply(mat, 1, scale))

col_anno <- HeatmapAnnotation(
  Treatment = c("Old_Ctrl","Old_Ctrl","Old_Ctrl","Old_NR","Old_NR"),
  col = list(Treatment = c(
    "Old_Ctrl" = pal_startrek("uniform")(2)[[2]],
    "Old_NR"   = pal_startrek("uniform")(1)
  ))
)

ht <- Heatmap(
  mat,
  name = "Z-score",
  top_annotation = col_anno,
  cluster_columns = FALSE,
  show_row_names = TRUE
)
