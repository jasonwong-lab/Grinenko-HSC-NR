library(dplyr)
library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(plyranges)
library(DSS)
library(purrr)
library(rtracklayer)
library(readxl)
library(AnnotationDbi)
library(ggplot2)
library(ggsignif)
library(ggsci)
library(extrafont)
library(stringr)
library(tidyr)
library(tibble)

theme_set(theme_classic(base_family='Arial'))

### ---- DSS and Annotation ---- ###
## Input Reference files ##

standard_chr <- paste0("chr", c(1:19, "X", "Y", "M"))

# CpG islands
cpgis <- fread("cpgIslandExt.txt.gz") %>%
  select(V2,V3,V4) %>%
  setNames(c("seqnames","start","end")) %>%
  makeGRangesFromDataFrame() %>%
  {.[seqnames(.) %in% standard_chr]}

# Enhancers
enhancers <- read_xls("Enhancers.xls", sheet = 1) %>%
  filter(`LT_HSC Tag Count in 2000 bp (5185076.0 Total, normalization factor = 1.93, effective total = 10000000)` != 0) %>%
  select(Chr, Start, End, Strand) %>%
  makeGRangesFromDataFrame() %>%
  liftOver(import.chain("mm9ToMm10.over.chain")) %>% unlist() %>%
  {.[seqnames(.) %in% standard_chr]}

# Promoters
txdb <- makeTxDbFromGFF(
    "ncbiRefSeqSelect.noBin.gtf",
    format = "gtf"
)

promoters <- transcripts(txdb) %>% promoters(3000,3000) %>%
  as.data.frame() %>% select(seqnames,start,end,strand) %>%
  makeGRangesFromDataFrame() %>%
  {.[seqnames(.) %in% standard_chr]}

# Gene Information
tx <- transcripts(txdb, columns = c("tx_name"))
tss_gr <- promoters(tx, upstream = 0, downstream = 1)
ref <- read.table(
  "ncbiRefSeqSelect.noBin.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(ref) <- c(
  "tx_name", "chrom", "strand", "txStart", "txEnd",
  "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
  "score", "SYMBOL", "cdsStartStat", "cdsEndStat", "exonFrames"
)
gene_annot <- unique(ref[, c("tx_name", "SYMBOL")])
mcols(tss_gr) <- as.data.frame(mcols(tss_gr)) %>%
  left_join(gene_annot, by = "tx_name")
tss_gr$gene_id <- tss_gr$tx_name
tss_gr <- tss_gr[seqnames(tss_gr) %in% standard_chr]

# TFs
tfs <- rtracklayer::import("reMap2022.bb")
tfs <- tfs %>% plyranges::filter(seqnames %in% standard_chr)

## Overlaps and coverage filtering ##
methyl_files <- list.files("path/to/bismarck/", pattern="\\.CpG_report\\.txt\\.gz$", full.names=TRUE)

methyl <- lapply(methyl_files, function(f) {
  fread(f) %>%
    mutate(end=V2) %>%
    rename(chr=V1,start=V2,strand=V3,methylated=V4,unmethylated=V5,cpg=V6,motif=V7) %>%
    filter(chr %in% standard_chr) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
})

features <- list(cpgis, promoters, enhancers)

methyl_regions <- lapply(methyl, function(m) {
  lapply(features, function(f) join_overlap_inner(f,m))
})

methyl_regions_granges <- lapply(methyl_regions, function(meth) {
  lapply(meth, function(region) {
    region %>%
      group_by(seqnames,start,end,strand) %>%
      summarize(methylated=sum(methylated,na.rm=TRUE), unmethylated=sum(unmethylated,na.rm=TRUE)) %>%
      makeGRangesFromDataFrame(keep.extra.columns=TRUE)
  })
})

methyl_regions <- lapply(methyl_regions_granges, function(meth) {
  lapply(meth, function(region) {
    region %>%
      plyranges::mutate(N = methylated + unmethylated) %>%
      plyranges::select(-unmethylated) %>%
      as.data.frame() %>%
      select(-c(strand,width)) %>%
      relocate(end,.after=N) %>%
      relocate(N,.after=start) %>%
      setNames(c("chr","pos","N","X","end"))
  })
})

methyl_regions <- list_transpose(methyl_regions)
methyl_regions <- lapply(methyl_regions, function(meth) lapply(meth, function(sample) filter(sample,N>=5)))

### ---- DSS using helper functions ---- ###

as_plain_df <- function(x) {
  df <- as.data.frame(x)
  class(df) <- "data.frame"
  df
}

format_dml_results <- function(dml_results, comparison_regions, filter_type = c("pval", "fdr")) {
  filter_type <- match.arg(filter_type)

  lapply(seq_along(dml_results), function(index) {
    df <- as_plain_df(dml_results[[index]])

    if (filter_type == "pval") {
      df <- df %>% filter(pval <= 0.05)
    }

    if (filter_type == "fdr") {
      df <- df %>%
        mutate(fdr = p.adjust(pval, method = "BH")) %>%
        filter(fdr <= 0.05)
    }

    annot <- do.call("rbind", comparison_regions[[index]])

    df %>%
      left_join(annot, by = c("chr", "pos")) %>%
      select(-c(N, X)) %>%
      relocate(end, .after = pos) %>%
      rename(start = pos) %>%
      distinct(.keep_all = TRUE)
  })
}

sample_names <- gsub("_S1.*", "", basename(methyl_files))

# Old Ctrl vs Old NR
old_ctrl_old_nr <- lapply(methyl_regions, function(regions) {
  list(regions[[1]], regions[[2]], regions[[3]], regions[[4]])
})
old_ctrl_old_nr_names <- sample_names[1:4]
diff_methyl_regions_old_ctrl_old_nr <- lapply(old_ctrl_old_nr, function(region) {
  makeBSseqData(region, old_ctrl_old_nr_names) %>%
    DMLtest(
      group1 = old_ctrl_old_nr_names[grep("NR", old_ctrl_old_nr_names)],
      group2 = old_ctrl_old_nr_names[grep("PBS", old_ctrl_old_nr_names)]
    ) %>%
    callDML(delta = 0, p.threshold = 1)
})
diff_methyl_regions_old_ctrl_old_nr_p <- format_dml_results(
  diff_methyl_regions_old_ctrl_old_nr,
  old_ctrl_old_nr,
  filter_type = "pval"
)

# Young Ctrl vs Old Ctrl
young_ctrl_old_ctrl <- lapply(methyl_regions, function(regions) {
  list(regions[[1]], regions[[2]], regions[[5]], regions[[6]])
})
young_ctrl_old_ctrl_names <- c(sample_names[1:2], sample_names[5:6])
diff_methyl_regions_young_ctrl_old_ctrl <- lapply(young_ctrl_old_ctrl, function(region) {
  makeBSseqData(region, young_ctrl_old_ctrl_names) %>%
    DMLtest(
      group1 = young_ctrl_old_ctrl_names[grep("old", young_ctrl_old_ctrl_names)],
      group2 = young_ctrl_old_ctrl_names[grep("young", young_ctrl_old_ctrl_names)]
    ) %>%
    callDML(delta = 0, p.threshold = 1)
})
diff_methyl_regions_young_ctrl_old_ctrl_p <- format_dml_results(
  diff_methyl_regions_young_ctrl_old_ctrl,
  young_ctrl_old_ctrl,
  filter_type = "pval"
)

# Young Ctrl vs Old NR
young_ctrl_old_nr <- lapply(methyl_regions, function(regions) {
  list(regions[[3]], regions[[4]], regions[[5]], regions[[6]])
})
young_ctrl_old_nr_names <- sample_names[3:6]
diff_methyl_regions_young_ctrl_old_nr <- lapply(young_ctrl_old_nr, function(region) {
  makeBSseqData(region, young_ctrl_old_nr_names) %>%
    DMLtest(
      group1 = young_ctrl_old_nr_names[grep("old", young_ctrl_old_nr_names)],
      group2 = young_ctrl_old_nr_names[grep("young", young_ctrl_old_nr_names)]
    ) %>%
    callDML(delta = 0, p.threshold = 1)
})
diff_methyl_regions_young_ctrl_old_nr_p <- format_dml_results(
  diff_methyl_regions_young_ctrl_old_nr,
  young_ctrl_old_nr,
  filter_type = "pval"
)


### ---- Figure 5b-d ---- ###

desired_order <- c("Old_Ctrl_NR", "Young_Old_Ctrl", "Young_Old_NR")
plot_fill <- pal_startrek("uniform")(2)[2]

make_violin_plot <- function(feature_index, feature_label, output_file) {
  plot_df <- Map(function(meth, name) {
    meth %>%
      mutate(percentage = diff / mu2) %>%
      as.data.frame() %>%
      mutate(Category = name) %>%
      select(Category, percentage)
  },
  list(
    diff_methyl_regions_old_ctrl_old_nr_p[[feature_index]],
    diff_methyl_regions_young_ctrl_old_ctrl_p[[feature_index]],
    diff_methyl_regions_young_ctrl_old_nr_p[[feature_index]]
  ),
  desired_order
  ) %>%
    bind_rows()

  plot_df$Category <- factor(plot_df$Category, levels = desired_order)

  p <- ggplot(plot_df, aes(x = Category, y = percentage)) +
    geom_violin(fill = plot_fill, alpha = 0.3) +
    geom_boxplot(
      width = 0.1,
      fill = plot_fill,
      outliers = FALSE,
      show.legend = FALSE
    ) +
    stat_summary(
      fun = mean,
      colour = "black",
      geom = "point",
      size = 3,
      alpha = 0.5
    ) +
    stat_summary(
      fun = mean,
      geom = "text",
      aes(label = round(..y.., 2)),
      vjust = -1.5,
      hjust = -0.5,
      size = 5
    ) +
    geom_signif(
      comparisons = list(
        c("Young_Old_Ctrl", "Young_Old_NR"),
        c("Young_Old_Ctrl", "Old_Ctrl_NR"),
        c("Young_Old_NR", "Old_Ctrl_NR")
      ),
      map_signif_level = TRUE,
      step_increase = 0.1,
      test = "t.test",
      textsize = 5
    ) +
    labs(
      x = NULL,
      y = "% methylation difference",
      title = feature_label
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 0.5, vjust = 0.5)
    )

  ggsave(output_file, p, device = cairo_pdf, width = 8, height = 8, units = "in")
  return(p)
}

cpgis_comp <- make_violin_plot(
  feature_index = 1,
  feature_label = "CpG islands",
  output_file = "Methylation_Difference_CpGIs_Violin.pdf"
)

promoters_comp <- make_violin_plot(
  feature_index = 2,
  feature_label = "Promoters",
  output_file = "Methylation_Difference_Promoters_Violin.pdf"
)

enhancers_comp <- make_violin_plot(
  feature_index = 3,
  feature_label = "Enhancers",
  output_file = "Methylation_Difference_Enhancers_Violin.pdf"
)


### ---- Figure 5e ---- ###
feature_totals <- c(
  CpG_islands = length(methyl_regions_granges[[1]][[1]]),
  Promoters   = length(methyl_regions_granges[[2]][[1]]),
  Enhancers   = length(methyl_regions_granges[[3]][[1]])
)

count_hyper_hypo <- function(feature_df) {
  c(
    Hyper = feature_df %>% filter(diff / mu2 >= 0.05) %>% nrow(),
    Hypo  = feature_df %>% filter(diff / mu2 <= -0.05) %>% nrow()
  )
}

make_bar_plot <- function(feature_index, feature_label, total_features, output_file) {

  counts <- list(
    Young_Ctrl_Old_Ctrl = count_hyper_hypo(diff_methyl_regions_young_ctrl_old_ctrl_p[[feature_index]]),
    Young_Ctrl_Old_NR   = count_hyper_hypo(diff_methyl_regions_young_ctrl_old_nr_p[[feature_index]]),
    Old_Ctrl_Old_NR     = count_hyper_hypo(diff_methyl_regions_old_ctrl_old_nr_p[[feature_index]])
  )

  count_df <- as.data.frame(counts) %>%
    rownames_to_column("Direction") %>%
    pivot_longer(
      cols = c(Young_Ctrl_Old_Ctrl, Young_Ctrl_Old_NR, Old_Ctrl_Old_NR),
      names_to = "Comparison",
      values_to = "Counts"
    )

  percent_df <- count_df %>%
    mutate(
      Percent = (Counts / total_features) * 100,
      Comparison = recode(
        Comparison,
        "Young_Ctrl_Old_Ctrl" = "Old/Young",
        "Young_Ctrl_Old_NR"   = "Old NR/Young",
        "Old_Ctrl_Old_NR"     = "Old NR/Old"
      ),
      Comparison = factor(
        Comparison,
        levels = c("Old/Young", "Old NR/Young", "Old NR/Old")
      ),
      Direction = factor(Direction, levels = c("Hyper", "Hypo"))
    )

  p <- ggplot(percent_df, aes(x = Comparison, y = Percent, fill = Direction)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8) +
    geom_text(
      aes(label = Counts),
      position = position_dodge(width = 0.9),
      vjust = -0.4,
      size = 5
    ) +
    scale_fill_manual(
      values = c(
        "Hyper" = "#A62A2A",
        "Hypo"  = "#556B2F"
      )
    ) +
    labs(
      x = NULL,
      y = "% of features",
      title = feature_label
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 0.5, vjust = 0.5),
      legend.title = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

  ggsave(output_file, p, device = cairo_pdf, width = 6, height = 5, units = "in")
  return(p)
}

cpgis_bar <- make_bar_plot(
  feature_index = 1,
  feature_label = "CpGIs",
  total_features = feature_totals["CpG_islands"],
  output_file = "Methylation_Diff_CpGIs_Bar_Plot.pdf"
)

promoters_bar <- make_bar_plot(
  feature_index = 2,
  feature_label = "Promoters",
  total_features = feature_totals["Promoters"],
  output_file = "Methylation_Diff_Promoters_Bar_Plot.pdf"
)

enhancers_bar <- make_bar_plot(
  feature_index = 3,
  feature_label = "Enhancers",
  total_features = feature_totals["Enhancers"],
  output_file = "Methylation_Diff_Enhancers_Bar_Plot.pdf"
)

### ---- Figure 5f-h ---- ###
# Annotating DML regions with gene information
add_gene_and_tss <- function(dml_df, tss_gr) {
  if (nrow(dml_df) == 0) return(dml_df)
  
  dml_gr <- makeGRangesFromDataFrame(
    dml_df,
    seqnames.field = "chr",
    start.field    = "start",
    end.field      = "end",
    keep.extra.columns = TRUE
  )
  
  hits <- distanceToNearest(dml_gr, tss_gr)
  
  tss_hits <- tss_gr[subjectHits(hits)]
  dml_hits <- dml_gr[queryHits(hits)]
  
  # centre of DML region
  dml_center <- round((start(dml_hits) + end(dml_hits)) / 2)
  
  # signed distance based off of TSS strand
  signed_dist <- ifelse(
    strand(tss_hits) == "+",
    dml_center - start(tss_hits),
    start(tss_hits) - dml_center
  )
  
  dml_gr$nearest_gene_id     <- tss_hits$gene_id
  dml_gr$nearest_gene_symbol <- tss_hits$SYMBOL
  dml_gr$dist_to_TSS_bp      <- signed_dist
  
  as.data.frame(dml_gr)
}

annotated_results <- list(
  old_ctrl_old_nr = lapply(
    diff_methyl_regions_old_ctrl_old_nr_p,
    add_gene_and_tss,
    tss_gr = tss_gr
  ),

  young_ctrl_old_ctrl = lapply(
    diff_methyl_regions_young_ctrl_old_ctrl_p,
    add_gene_and_tss,
    tss_gr = tss_gr
  ),

  young_ctrl_old_nr = lapply(
    diff_methyl_regions_young_ctrl_old_nr_p,
    add_gene_and_tss,
    tss_gr = tss_gr
  )
)

# Exporting gene lists to be used in Metascape (online)
feature_names <- c("CpGI", "Promoter", "Enhancer")

for (comp in names(annotated_results)) {
  names(annotated_results[[comp]]) <- feature_names
}

extract_promoter_genes <- function(annotated_results, comp, direction) {
  
  promoters_df <- annotated_results[[comp]][["Promoter"]] %>%
    mutate(methylation_difference = diff / mu2)
  
  if (direction == "hyper") {
    promoters_df <- promoters_df %>%
      filter(methylation_difference >= 0.05)
  }
  
  if (direction == "hypo") {
    promoters_df <- promoters_df %>%
      filter(methylation_difference <= -0.05)
  }
  
  promoters_df %>%
    pull(nearest_gene_symbol) %>%
    na.omit() %>%
    unique() %>%
    sort()
}

write_metascape_list <- function(genes, output_file) {
  write.table(
    genes,
    file = output_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

hyper_promoters_old_young_genes <- extract_promoter_genes(
  annotated_results,
  comp = "young_ctrl_old_ctrl",
  direction = "hyper"
)

hypo_promoters_old_young_genes <- extract_promoter_genes(
  annotated_results,
  comp = "young_ctrl_old_ctrl",
  direction = "hypo"
)

hypo_promoters_old_nr_old_genes <- extract_promoter_genes(
  annotated_results,
  comp = "old_ctrl_old_nr",
  direction = "hypo"
)

write_metascape_list(
  hyper_promoters_old_young_genes,
  "hyper_promoters_old_vs_young.txt"
)

write_metascape_list(
  hypo_promoters_old_young_genes,
  "hypo_promoters_old_vs_young.txt"
)

write_metascape_list(
  hypo_promoters_old_nr_old_genes,
  "hypo_promoters_old_nr_vs_old.txt"
)
# Use Metascape online at https://metascape.org with these gene lists for ORA results. 

### ---- Figure 5i-p ---- ###

tf_list <- c("BACH2", "SPIB", "CTCF", "NFAT", "SPI1", "NFE2", "ATF3", "TCF3")

nr_color   <- "#E64B35"
ctrl_color <- "#4DBBD5"

for (tf_in in tf_list) {
  if (tf_in == "NFAT") {
    tfs_filtered <- tfs %>% filter(str_detect(TF, "^NFAT"))
  } else {
    tfs_filtered <- tfs %>% filter(TF == tf_in)
  }
  # overlapping TF with methylation
  methyl_tfs <- lapply(seq_along(methyl), function(i) {

    meth <- methyl[[i]]
    fname <- methyl_files[i]

    tf_overlap <- join_overlap_inner(tfs_filtered, meth)

    tf_overlap %>%
      as_tibble() %>%
      group_by(seqnames, start, end, strand, TF) %>%
      summarise(
        methylated   = sum(methylated, na.rm = TRUE),
        unmethylated = sum(unmethylated, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        beta = methylated / (methylated + unmethylated),
        treatment = ifelse(str_detect(fname, "PBS"), "Ctrl", "NR"),
        age       = ifelse(str_detect(fname, "young"), "Young", "Old"),
        comparison = paste(age, treatment)
      )
  })

  methyl_tfs_plot <- bind_rows(methyl_tfs)


  methyl_tfs_plot <- methyl_tfs_plot %>%
    filter(comparison %in% c("Young Ctrl", "Old Ctrl", "Old NR"))

  methyl_tfs_plot$comparison <- factor(
    methyl_tfs_plot$comparison,
    levels = c("Young Ctrl", "Old Ctrl", "Old NR")
  )
  comparisons_all <- list(
    c("Young Ctrl", "Old Ctrl"),
    c("Young Ctrl", "Old NR"),
    c("Old Ctrl", "Old NR")
  )

  p <- ggplot(methyl_tfs_plot,
              aes(x = comparison, y = beta, fill = treatment)) +

    geom_violin(alpha = 0.3, trim = TRUE) +

    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      colour = "black",
      fill = NA,
      show.legend = FALSE
    ) +

    stat_summary(
      fun = mean,
      geom = "text",
      aes(label = round(..y.., 2)),
      vjust = 0.5,
      hjust = -0.3,
      size = 5
    ) +

    geom_signif(
      comparisons = comparisons_all,
      test = "t.test",
      map_signif_level = TRUE,
      step_increase = 0.1,
      textsize = 4
    ) +

    scale_y_continuous(
    limits = c(0, 1.4),
    expand = expansion(mult = c(0, 0.1))
    ) +

    scale_fill_manual(
      values = c("NR" = nr_color, "Ctrl" = ctrl_color),
      name = "Treatment"
    ) +

    theme_classic() +
    labs(
      x = NULL,
      y = "Methylation Beta"
    ) +
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.text.x = element_text(angle = 25, hjust = 0.5),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16)
    )

  ggsave(
    filename = paste0(tf_in, "_methylation_violin.pdf"),
    plot = p,
    width = 6,
    height = 6,
    units = "in"
  )
}
