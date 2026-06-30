library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
library(DiffBind)
library(DESeq2)
library(ggrepel)
library(data.table)
library(ggplot2)
library(ggalt)
library(dplyr)
library(tibble)
library(trekcolors)
library(stringr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(GenomicRanges)
library(matrixStats)
library(SummarizedExperiment)
library(rtracklayer)

### ---- Differential Expression Analysis ---- ###
rnaseq_dir <- "/path/to/output/final/${PREFIX}.sorted.bam"
rnaseq_cd_list <- system(paste0("ls ", rnaseq_dir, "*.sorted.bam"), intern = TRUE)
bams_cd <- BamFileList(rnaseq_cd_list)

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


rnaseq_cd <- summarizeOverlaps(
    features = ebg, reads = bams_cd,
    mode = "Union",
    singleEnd = FALSE,
    ignore.strand = FALSE,
    fragments = TRUE,
    BPPARAM = MulticoreParam(workers = 16)
)

coldata_cd <- data.frame(
  File = c(
    "L143639_Track-187082Aligned.sortedByCoord.out.bam",
    "L143640_Track-187083Aligned.sortedByCoord.out.bam",
    "L143641_Track-187084Aligned.sortedByCoord.out.bam",
    "L143642_Track-187085Aligned.sortedByCoord.out.bam",
    "L143643_Track-187086Aligned.sortedByCoord.out.bam",
    "L143644_Track-187087Aligned.sortedByCoord.out.bam",
    "L143645_Track-187094Aligned.sortedByCoord.out.bam",
    "L143646_Track-187095Aligned.sortedByCoord.out.bam",
    "L143647_Track-187096Aligned.sortedByCoord.out.bam",
    "L143648_Track-187097Aligned.sortedByCoord.out.bam",
    "L143649_Track-187098Aligned.sortedByCoord.out.bam",
    "L143650_Track-187099Aligned.sortedByCoord.out.bam"
  ),
  Name = c(
    "S1_CD4_ctrl",
    "S2_CD4_ctrl",
    "S3_CD4_ctrl",
    "S4_CD4_NR",
    "S5_CD4_NR",
    "S6_CD4_NR",
    "S7_CD8_ctrl",
    "S8_CD8_ctrl",
    "S9_CD8_ctrl",
    "S10_CD8_NR",
    "S11_CD8_NR",
    "S12_CD8_NR"
  ),
  Treatment = c(
    "Control","Control","Control",
    "NR","NR","NR",
    "Control","Control","Control",
    "NR","NR","NR"
  ),
  stringsAsFactors = FALSE
)
colData(rnaseq_cd) <- DataFrame(coldata_cd)
rnaseq_cd$Treatment <- relevel(factor(rnaseq_cd$Treatment, ordered = FALSE), "Control")
colnames(rnaseq_cd) <- colData(rnaseq_cd)$Name

rnaseq_cd4 <- rnaseq_cd[, 1:6]
rnaseq_cd8 <- rnaseq_cd[, 7:12]


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

# Process CD4 data
rnaseq_cd4 <- rnaseq_cd4[rowSums(assay(rnaseq_cd4)) >= 5, ]
dss_cd4 <- DESeqDataSet(rnaseq_cd4, design = ~Treatment)
dss_cd4 <- DESeq(dss_cd4)
NR_ctrl_res_ordered_cd4 <- process_deseq_results(dss_cd4, c("Treatment", "NR", "Control"), edb, "RNAseq_CD4.tsv")
vsd_cd4 <- vst(dss_cd4)

# Process CD8 data
rnaseq_cd8 <- rnaseq_cd8[rowSums(assay(rnaseq_cd8)) >= 5, ]
dss_cd8 <- DESeqDataSet(rnaseq_cd8, design = ~Treatment)
dss_cd8 <- DESeq(dss_cd8)
NR_ctrl_res_ordered_cd8 <- process_deseq_results(dss_cd8, c("Treatment", "NR", "Control"), edb, "RNAseq_CD8.tsv")
vsd_cd8 <- vst(dss_cd8)



### ---- Figure 6b ---- ###
getPCAs <- function(object, intgroup = "condition", ntop = 500, returnData = FALSE) {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2 / sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    } else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], 
                    group = group, intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:4]
        return(d)
    }
}

# Reusable plotting function
plot_pca_pairs <- function(vsd_obj, prefix, pc_percentages) {
    pca_data <- getPCAs(vsd_obj, "Treatment", returnData = TRUE) %>%
        mutate(name = sub("_.*", "", name)) %>%
        mutate(group = as.character(group), Treatment = as.character(Treatment)) %>%
        mutate(group = ifelse(group == "Control", "Old_Ctrl", group),
               group = ifelse(group == "NR", "Old_NR", group),
               Treatment = ifelse(Treatment == "Control", "Old_Ctrl", Treatment),
               Treatment = ifelse(Treatment == "NR", "Old_NR", Treatment))
    
    pcs <- combn(c("PC1", "PC2", "PC3", "PC4"), 2, simplify = FALSE)
    plots <- lapply(pcs, function(pc) {
        ggplot(pca_data, aes(x = !!sym(pc[1]), y = !!sym(pc[2]), color = Treatment)) +
            geom_encircle(aes(fill = Treatment), alpha = 0.3) +
            geom_point() +
            geom_text_repel(aes(label = paste0(name)), show.legend = FALSE, size = 5) +
            theme_classic() +
            scale_color_manual(values = c("Old_Ctrl" = "red", "Old_NR" = "blue")) +
            scale_fill_manual(values = c("Old_Ctrl" = "red", "Old_NR" = "blue"))+
            guides(color = guide_legend(override.aes = list(linetype = 0))) +
            xlab(paste0(pc[1], " (", signif(pc_percentages[as.numeric(str_extract(pc[1], "\\d+"))], 3) * 100, "%)")) +
            ylab(paste0(pc[2], " (", signif(pc_percentages[as.numeric(str_extract(pc[2], "\\d+"))], 3) * 100, "%)")) +
            theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16),
                  legend.text = element_text(size = 16), legend.title = element_text(size = 16),
                  legend.key.size = unit(1.5, "lines"))
    })
    
    # Save all 6 pairs as PDF (multi-page)
    pdf_file <- paste0(prefix, "_PCA_pairs.pdf")
    pdf(pdf_file, height = 4, width = 5.5)
    lapply(plots, print)
    dev.off()
    
    # Individual PNGs too (optional)
    for (i in seq_along(pcs)) {
        ggsave(paste0(prefix, "_PCA_", pcs[[i]][1], "_", pcs[[i]][2], ".png"), 
               plots[[i]], height = 4, width = 5.5)
    }
    message("Saved: ", pdf_file)
}


pc_percentages_cd4 <- attr(getPCAs(vsd_cd4, "Treatment", returnData = TRUE), "percentVar")
plot_pca_pairs(vsd_cd4, "CD4", pc_percentages_cd4)

pc_percentages_cd8 <- attr(getPCAs(vsd_cd8, "Treatment", returnData = TRUE), "percentVar")
plot_pca_pairs(vsd_cd8, "CD8", pc_percentages_cd8)


### ---- Figure 6c ---- ###
mat <- counts(dss_cd4, normalized = TRUE)

mat <- as.data.frame(mat) %>%
  rownames_to_column("ensembl_raw")

# Keep only significant DEGs (padj filter)
deg_keep <- NR_ctrl_res_ordered_cd4 %>%
  filter(padj <= 0.05) %>%
  pull(ensembl)

mat <- mat %>%
  filter(ensembl_raw %in% deg_keep)

mat <- mat %>%
  mutate(ensembl = sub("\\..*$", "", ensembl_raw))
colnames(mat) <- c(
  "ensembl_raw",
  sub("_.*$", "", colnames(mat)[2:(ncol(mat)-1)]),
  "ensembl"
)

mat <- mat %>%
  left_join(
    AnnotationDbi::select(
      edb,
      keys = mat$ensembl,
      keytype = "GENEID",
      columns = "GENENAME"
    ) %>%
      rename(
        ensembl = GENEID,
        external_gene_name = GENENAME
      ),
    by = "ensembl"
  ) %>%
  mutate(
    external_gene_name = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl,
      external_gene_name
    )
  ) %>%
  relocate(external_gene_name, .before = ensembl_raw) %>%
  relocate(ensembl, .after = ensembl_raw)

expr_cols <- setdiff(
  colnames(mat),
  c("external_gene_name", "ensembl_raw", "ensembl")
)

zmat <- mat %>%
  mutate(across(all_of(expr_cols), as.numeric)) %>%
  rowwise() %>%
  mutate(across(
    all_of(expr_cols),
    ~ (. - mean(c_across(all_of(expr_cols)))) /
      sd(c_across(all_of(expr_cols)))
  )) %>%
  ungroup()

zmat <- zmat %>%
  select(-c(ensembl, ensembl_raw)) %>%
  group_by(external_gene_name) %>%
  mutate(
    external_gene_name = ifelse(
      row_number() > 1,
      paste0(external_gene_name, "_", row_number() - 1),
      external_gene_name
    )
  ) %>%
  ungroup() %>%
  column_to_rownames("external_gene_name")


genes_to_label <- c("Stat1", "Irgm1", "Ly6a", "Slc15a2")

annotation_col <- data.frame(
  Condition = factor(rep(c("Ctrl", "NR"), each = 3))
)

rownames(annotation_col) <- colnames(zmat)

labeled_data <- zmat
rownames(labeled_data) <- ifelse(
  rownames(zmat) %in% genes_to_label,
  rownames(zmat),
  ""
)

png(
  filename = "RNAseq_CD4_Sample_Heatmap.png",
  width = 2000,
  height = 3000,
  res = 200
)

pheatmap(
  labeled_data,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_cols = FALSE,
  show_rownames = TRUE,
  fontsize_row = 8,
  annotation_col = annotation_col,
  cellheight = 0.6,
  gaps_col = 3,
  main = "CD4+ T Cell Expression\nHighlighted Genes : Stat1, Irgm1, Ly6a, Slc15a2"
)

dev.off()


### ---- Figure 6d ---- ###

mat <- counts(dss_cd8, normalized = TRUE)

mat <- as.data.frame(mat) %>%
  rownames_to_column("ensembl_raw")

# Keep only significant DEGs (padj filter)
deg_keep <- NR_ctrl_res_ordered_cd8 %>%
  filter(padj <= 0.05) %>%
  pull(ensembl)

mat <- mat %>%
  filter(ensembl_raw %in% deg_keep)

mat <- mat %>%
  mutate(ensembl = sub("\\..*$", "", ensembl_raw))

colnames(mat) <- c(
  "ensembl_raw",
  sub("_.*$", "", colnames(mat)[2:(ncol(mat)-1)]),
  "ensembl"
)

mat <- mat %>%
  left_join(
    AnnotationDbi::select(
      edb,
      keys = mat$ensembl,
      keytype = "GENEID",
      columns = "GENENAME"
    ) %>%
      rename(
        ensembl = GENEID,
        external_gene_name = GENENAME
      ),
    by = "ensembl"
  ) %>%
  mutate(
    external_gene_name = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl,
      external_gene_name
    )
  ) %>%
  relocate(external_gene_name, .before = ensembl_raw) %>%
  relocate(ensembl, .after = ensembl_raw)

expr_cols <- setdiff(
  colnames(mat),
  c("external_gene_name", "ensembl_raw", "ensembl")
)

zmat <- mat %>%
  mutate(across(all_of(expr_cols), as.numeric)) %>%
  rowwise() %>%
  mutate(across(
    all_of(expr_cols),
    ~ (. - mean(c_across(all_of(expr_cols)))) /
      sd(c_across(all_of(expr_cols)))
  )) %>%
  ungroup()

zmat <- zmat %>%
  select(-c(ensembl, ensembl_raw)) %>%
  group_by(external_gene_name) %>%
  mutate(
    external_gene_name = ifelse(
      row_number() > 1,
      paste0(external_gene_name, "_", row_number() - 1),
      external_gene_name
    )
  ) %>%
  ungroup() %>%
  column_to_rownames("external_gene_name")

annotation_col <- data.frame(
  Condition = factor(rep(c("Ctrl", "NR"), each = 3))
)

rownames(annotation_col) <- colnames(zmat)

png(
  filename = "RNAseq_CD8_Sample_Heatmap.png",
  width = 2000,
  height = 3000,
  res = 200
)

pheatmap(
  zmat,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_cols = FALSE,
  show_rownames = TRUE,
  fontsize_row = 8,
  annotation_col = annotation_col,
  cellheight = 0.6,
  gaps_col = 3,
  main = "CD8+ T Cell Expression"
)

dev.off()
