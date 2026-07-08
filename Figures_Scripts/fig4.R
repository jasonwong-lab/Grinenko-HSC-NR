library(DiffBind)
library(GenomicFeatures)
library(ChIPseeker)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(forcats)
library(readxl)
library(GenomicRanges)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(stringr)
library(grid)
library(readr)
library(tidyr)
library(stringr)

### ---- Calling DARs ---- ###

## DiffBind processing ##
# Replace the example file names below with the corresponding outputs from preprocessing_commands.sh.
samples <- data.frame(
  SampleID = c("S1", "S2", "S3", "S4", "S5"),
  Condition = c("CTRL", "CTRL", "CTRL", "NR", "NR"),
  Replicate = c(1, 2, 3, 1, 2),
  bamReads = c(
    "S1_cont.markdup.bam",
    "S1_ctrl_2.markdup.bam",
    "s2_ctrl_3.markdup.bam",
    "S2_NR_tr.markdup.bam",
    "S3_NR_tr.markdup.bam"
  ),
  Peaks = c(
    "macs2_outdir/S1_cont_peaks.xls",
    "macs2_outdir/S1_ctrl_2_peaks.xls",
    "macs2_outdir/s2_ctrl_3_peaks.xls",
    "macs2_outdir/S2_NR_tr_peaks.xls",
    "macs2_outdir/S3_NR_tr_peaks.xls"
  ),
  PeakCaller = "macs"
)

atac <- dba(sampleSheet = samples)
atac_count <- dba.count(atac)

atac_clean <- dba.blacklist(
  atac_count,
  blacklist = DBA_BLACKLIST_MM10,
  greylist = FALSE
)

atac_norm <- dba.normalize(
  atac_clean,
  background = TRUE,
  normalize = DBA_NORM_NATIVE,
  method = DBA_DESEQ2,
  library = DBA_LIBSIZE_BACKGROUND,
  bParallel = FALSE
)

## Differential analysis ##
atac_contrast <- dba.contrast(atac_norm, minMembers = 2)
atac_diff <- dba.analyze(atac_contrast)

atac_stats <- dba.report(atac_diff, bUsePval = TRUE, th = 1)
seqlevelsStyle(atac_stats) <- "UCSC"

dars_df <- as.data.frame(atac_stats)
sig_dars <- dars_df %>%
  filter(`p.value` <= 0.01)

open_dars <- dars_df %>%
  filter(`p.value` <= 0.01, Fold > 0)

closed_dars <- dars_df %>%
  filter(`p.value` <= 0.01, Fold < 0)


### ---- Annotating DARs ---- ###

## Load gene model ##
txdb <- makeTxDbFromGFF(
    "ncbiRefSeqSelect.noBin.gtf",
    format = "gtf"
)

## Annotate open and closed DARs ##
peakAnno_open <- annotatePeak(
  GRanges(open_dars),
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  verbose = FALSE
)

peakAnno_closed <- annotatePeak(
  GRanges(closed_dars),
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  verbose = FALSE
)

open_df <- as.data.frame(peakAnno_open)
closed_df <- as.data.frame(peakAnno_closed)

## Load enhancer set ##
standard_chr <- paste0("chr", c(1:19, "X", "Y"))
enh_raw <- read_xls("Enhancers.xls", sheet = 1)

enh_filt <- enh_raw %>%
  filter(`LT_HSC Tag Count in 2000 bp (5185076.0 Total, normalization factor = 1.93, effective total = 10000000)` != 0) %>%
  dplyr::select(Chr, Start, End, Strand)

gr_mm9 <- makeGRangesFromDataFrame(
  enh_filt,
  seqnames.field = "Chr",
  start.field = "Start",
  end.field = "End",
  strand.field = "Strand"
)
chain <- import.chain("mm9ToMm10.over.chain")
gr_mm10 <- liftOver(gr_mm9, chain) %>% unlist()

seqlevelsStyle(gr_mm10) <- "UCSC"
enhancers <- gr_mm10[seqnames(gr_mm10) %in% standard_chr]

## Re-label intergenic peaks that overlap enhancers ##
mark_enhancers <- function(df, peaks, enhancers) {
  ov <- findOverlaps(peaks, enhancers, ignore.strand = TRUE)
  idx <- intersect(queryHits(ov), grep("intergenic", df$annotation, ignore.case = TRUE))
  df$annotation[idx] <- "Enhancer"
  df
}
open_df <- mark_enhancers(open_df, as.GRanges(peakAnno_open), enhancers)
closed_df <- mark_enhancers(closed_df, as.GRanges(peakAnno_closed), enhancers)

### ---- Figure 4a ---- ###
simplify_annotation <- function(x) {
  x <- as.character(x)

  case_when(
    grepl("^Promoter", x)                            ~ "Promoter-TSS",
    x == "Promoter-TSS"                              ~ "Promoter-TSS",
    x == "TTS"                                       ~ "TTS",
    x == "Enhancer"                                  ~ "Enhancer",
    grepl("intergenic", x, ignore.case = TRUE)       ~ "Intergenic",
    grepl("^Intron", x)                              ~ "Intron",
    grepl("^Exon", x)                                ~ "Exon",
    grepl("5' UTR", x)                               ~ "5' UTR",
    grepl("3' UTR", x)                               ~ "3' UTR",
    TRUE                                             ~ "Other"
  )
}

open_df$annotation_simple   <- simplify_annotation(open_df$annotation)
closed_df$annotation_simple <- simplify_annotation(closed_df$annotation)

dar_df <- bind_rows(
  open_df %>% mutate(DAR_status = "Open"),
  closed_df %>% mutate(DAR_status = "Closed")
) %>%
  mutate(
    DAR_status = factor(DAR_status, levels = c("Open", "Closed")),
    annotation_simple = if_else(
      is.na(annotation_simple) | annotation_simple == "",
      "Other",
      annotation_simple
    ),
    annotation_simple = factor(
      annotation_simple,
      levels = c(
        "Promoter-TSS", "Intergenic", "Intron", "Exon",
        "Enhancer", "TTS", "5' UTR", "3' UTR", "Other"
      )
    )
  )

dar_count <- dar_df %>%
  count(DAR_status, annotation_simple, name = "Count")

col_anno <- c(
  "Promoter-TSS" = "red",
  "Intergenic" = "blue",
  "Intron" = "purple",
  "Exon" = "#F781BF",
  "Enhancer" = "orange",
  "TTS" = "grey50",
  "5' UTR" = "darkgreen",
  "3' UTR" = "darkolivegreen3",
  "Other" = "grey80"
)

p <- ggplot(dar_count, aes(x = DAR_status, y = Count, fill = annotation_simple)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = col_anno, drop = FALSE, name = "Annotation") +
  labs(x = NULL, y = "Number of DARs Old_NR vs Old_Ctrl") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 16, colour = "black"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.7, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

print(p)

### ---- Figure 4b ---- ###

norm_gr <- dba.peakset(atac_norm, bRetrieve = TRUE)

open_gr <- makeGRangesFromDataFrame(
  open_dars,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand",
  keep.extra.columns = TRUE
)

closed_gr <- makeGRangesFromDataFrame(
  closed_dars,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand",
  keep.extra.columns = TRUE
)

my_dars <- c(open_gr, closed_gr)
my_annotations <- c(open_df$annotation_simple, closed_df$annotation_simple)

hits <- findOverlaps(my_dars, norm_gr)
keep <- !duplicated(subjectHits(hits))
qidx <- queryHits(hits)[keep]
sidx <- subjectHits(hits)[keep]

dar_gr <- norm_gr[sidx]

count_df <- as.data.frame(mcols(dar_gr))
count_df <- count_df[, vapply(count_df, is.numeric, logical(1)), drop = FALSE]
count_mat <- as.matrix(count_df)

z_mat <- t(scale(t(count_mat)))
z_mat[is.na(z_mat)] <- 0

dar_status <- factor(
  ifelse(my_dars$Fold[qidx] > 0, "Open", "Closed"),
  levels = c("Open", "Closed")
)

anno <- my_annotations[qidx]
anno <- as.character(anno)
anno[is.na(anno)] <- "Other"
anno <- factor(
  anno,
  levels = c("Promoter-TSS", "Intergenic", "Intron", "Exon", "Enhancer", "TTS", "5' UTR", "3' UTR", "Other")
)

col_anno <- c(
  "Promoter-TSS" = "red",
  "Intergenic" = "lightblue",
  "Intron" = "purple",
  "Exon" = "black",
  "Enhancer" = "orange",
  "TTS" = "grey50",
  "5' UTR" = "darkgreen",
  "3' UTR" = "darkolivegreen3",
  "Other" = "grey80"
)

treat <- factor(c("Ctrl", "Ctrl", "Ctrl", "NR", "NR"), levels = c("Ctrl", "NR"))

ha_top <- HeatmapAnnotation(
  Treatment = treat,
  col = list(Treatment = c("Ctrl" = "blue", "NR" = "red"))
)

ha_row <- rowAnnotation(
  Annotation = anno,
  col = list(Annotation = col_anno),
  width = unit(6, "mm")
)

ht <- Heatmap(
  z_mat,
  name = "Z score",
  top_annotation = ha_top,
  left_annotation = ha_row,
  row_split = dar_status,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 18),
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

print(ht)

#### ---- Figure 4c ---- ###

## Read HOMER known motif results
read_homer_known <- function(path) {
  raw <- read_tsv(
    file.path(path, "knownResults.txt"),
    col_names = TRUE,
    show_col_types = FALSE
  )
  
  raw %>%
    transmute(
      motif_name = str_replace(`Motif Name`, "/.*$", ""),
      p_value = as.numeric(`P-value`),
      tgt_pct = as.numeric(str_remove(`% of Target Sequences with Motif`, "%")) / 100,
      bgd_pct = as.numeric(str_remove(`% of Background Sequences with Motif`, "%")) / 100,
      fold_enrichment = tgt_pct / bgd_pct,
      log_p_value = -log10(p_value)
    ) %>%
    filter(
      !is.na(motif_name),
      !is.na(p_value),
      is.finite(fold_enrichment),
      is.finite(log_p_value)
    )
}

# Replace the example paths below with the corresponding outputs from preprocessing_commands.sh.
known_ctrl <- read_homer_known("path/to/homer_output/ctrl")
known_nr <- read_homer_known("path/to/homer_output/nr")

motif_combined <- full_join(
  known_ctrl %>%
    select(motif_name, fold_enrichment, log_p_value, tgt_pct) %>%
    rename(
      fold_enrich_ctrl = fold_enrichment,
      pval_ctrl = log_p_value,
      tgt_pct_ctrl = tgt_pct
    ),
  known_nr %>%
    select(motif_name, fold_enrichment, log_p_value, tgt_pct) %>%
    rename(
      fold_enrich_nr = fold_enrichment,
      pval_nr = log_p_value,
      tgt_pct_nr = tgt_pct
    ),
  by = "motif_name"
) %>%
  mutate(
    log2_fc = log2(fold_enrich_nr / fold_enrich_ctrl),
    pval_max = pmax(pval_ctrl, pval_nr, na.rm = TRUE),
    direction = if_else(log2_fc > 0, "Up in NR", "Down in NR")
  ) %>%
  filter(is.finite(log2_fc), pval_max >= -log10(0.05))

motif_combined <- motif_combined %>% 
  group_by(motif_name) %>% 
  slice_head(n = 1) %>% 
  ungroup()

## Select top motifs in each direction ## 
top_motifs <- bind_rows(
  motif_combined %>%
    filter(direction == "Up in NR") %>%
    arrange(desc(abs(log2_fc)), desc(pval_max)) %>%
    slice_head(n = 30),
  motif_combined %>%
    filter(direction == "Down in NR") %>%
    arrange(desc(abs(log2_fc)), desc(pval_max)) %>%
    slice_head(n = 30)
) %>%
  group_by(motif_name, direction) %>%
  summarise(
    log2_fc = mean(log2_fc, na.rm = TRUE),
    pval_max = mean(pval_max, na.rm = TRUE),
    .groups = "drop"
  )


wrap_at_parentheses <- function(x) {
  x <- gsub("\\(", "\n\\(", x)
  x <- gsub("\\)", "\\)\n", x)
  str_wrap(x, width = 30)
}

fig4c <- top_motifs %>%
  arrange(desc(log2_fc)) %>% 
  mutate(
    Motif = factor(motif_name, levels = motif_name),
    direction = factor(direction, levels = c("Up in NR", "Down in NR"))
  )

wrapped_labels <- wrap_at_parentheses(levels(fig4c$Motif))

fig4c_plot <- ggplot(fig4c, aes(x = log2_fc, y = Motif, size = pval_max, colour = direction)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(1.5, 5)) +
  scale_colour_manual(values = c("Up in NR" = "red", "Down in NR" = "blue")) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  scale_y_discrete(labels = wrapped_labels) +
  labs(
    x = "log2(Old NR/Old)",
    y = "Motif",
    size = "-log10(P)",
    colour = "Direction"
  ) +
  theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    axis.text.y = element_text(lineheight = 0.9)
  )

print(fig4c_plot)
