library(readxl)
library(dplyr)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stringr)
library(patchwork)

### ---- Fig 1 a-b, Supplementary Fig 1 a-b ---- ###

xlsx_path <- "mmc7.xlsx"
outdir <- "GSEA_plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## GSEA Functions ##
get_gene_set <- function(pathway_name) {
  
  msig <- bind_rows(
    
    # KEGG + WikiPathways
    msigdbr(species = "Homo sapiens", collection = "C2") %>%
      select(gs_name, gene_symbol),
    
    # GO BP
    msigdbr(species = "Homo sapiens", collection = "C5") %>%
      select(gs_name, gene_symbol)
    
  ) %>%
    rename(term = gs_name, gene = gene_symbol) %>%
    distinct()
  
  msig %>%
    filter(term == pathway_name)
}


run_gsea <- function(sheet_name, gene_sets) {
  
  df <- read_excel(xlsx_path, sheet = sheet_name) %>%
    filter(!is.na(Gene), !is.na(stat)) %>%
    mutate(
      Gene = as.character(Gene),
      stat = as.numeric(stat)
    ) %>%
    filter(!is.na(stat)) %>%
    group_by(Gene) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ranks <- df$stat
  names(ranks) <- df$Gene
  ranks <- sort(ranks, decreasing = TRUE)
  
  GSEA(
    geneList = ranks,
    TERM2GENE = gene_sets,
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    seed = TRUE,
    by = "fgsea"
  )
}

plot_one_pathway <- function(gsea_obj, pathway_id, title_label) {
  
  res <- as.data.frame(gsea_obj)
  
  if (!(pathway_id %in% res$ID)) {
    stop(paste("Pathway not found:", pathway_id))
  }
  
  nes <- res$NES[res$ID == pathway_id]
  fdr <- res$p.adjust[res$ID == pathway_id]
  
  label_text <- paste0(
    "NES = ", signif(nes, 3),
    "\nFDR = ", signif(fdr, 3)
  )
  
  p <- enrichplot::gseaplot2(
    gsea_obj,
    geneSetID = pathway_id,
    title = title_label,
    base_size = 12,
    color = "darkgreen",
    rel_heights = c(1.5, 0.4, 0.8),
    pvalue_table = FALSE
  )
  
  p[[1]] <- p[[1]] +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = label_text,
      hjust = 1.1,
      vjust = 1.2,
      size = 5
    )
  
  return(p)
}

## USER INPUT NEEDED: change this to the desired pathway and title ##

pathway <- "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY"
title <- "Kegg T cell receptor signalling"

# To recreate figures, use the following pathways:
# "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY"
# "WP_TCELL_RECEPTOR_SIGNALING"  
# "GOBP_T_CELL_ACTIVATION" for supplementary figures 1 a-b

gene_set <- get_gene_set(pathway)

# Run GSEA
gsea_cd4 <- run_gsea("CD4", gene_set)
gsea_cd8 <- run_gsea("CD8", gene_set)

# Plot
p_cd4 <- plot_one_pathway(
  gsea_cd4,
  pathway,
  paste("CD4:", title)
)

p_cd8 <- plot_one_pathway(
  gsea_cd8,
  pathway,
  paste("CD8:", title)
)

### ---- Fig 1 c-d ---- ###
## Make Metascape Input ##
make_metascape_input <- function(sheet_name,
                                  padj_cutoff = 0.05) {

  df <- read_excel(xlsx_path, sheet = sheet_name)

  df <- df %>%
    filter(!is.na(Gene),
           !is.na(log2FoldChange),
           !is.na(padj))

  up <- df %>%
    filter(log2FoldChange > 0,
           padj < padj_cutoff) %>%
    pull(Gene) %>%
    unique()

  down <- df %>%
    filter(log2FoldChange < 0,
           padj < padj_cutoff) %>%
    pull(Gene) %>%
    unique()

  writeLines(up, file.path(outdir, paste0(sheet_name, "_UP_genes.txt")))
  writeLines(down, file.path(outdir, paste0(sheet_name, "_DOWN_genes.txt")))

  invisible(list(up = up, down = down))
}

make_metascape_input("CD4")
make_metascape_input("CD8")

# Use Metascape online at https://metascape.org with the gene lists for ORA results. 

### ---- Fig 1 e-f, Supplementary Fig 1 c-d ---- ###
## This section generates pathway-specific BED files to be used in fig1.sh ##

file <- "21598290cd181474-sup-213895_2_supp_5514125_pr2m22.xlsx"

write_pathway_beds <- function(
  sheet,
  collection,
  gene_set,
  subcollection = NULL
) {
  
  df <- read_excel(file, sheet = sheet)
  
  genes <- msigdbr(
    species = "Homo sapiens",
    collection = collection,
    subcollection = subcollection
  ) %>%
    filter(gs_name == gene_set) %>%
    pull(gene_symbol) %>%
    unique()
  
  peaks <- df %>%
    filter(!is.na(gene_symbol)) %>%
    filter(gene_symbol %in% genes) %>%
    filter(!is.na(dist_to_tss)) %>%
    filter(Direction == "Decreased")
  
  enhancers <- peaks %>%
    filter(abs(dist_to_tss) >= 3000) %>%
    select(chr, start, end)
  
  tss <- peaks %>%
    filter(abs(dist_to_tss) < 3000) %>%
    select(chr, start, end)
  
  clean_gene_set <- gsub("[^A-Za-z0-9_]+", "_", gene_set)
  enhancer_file <- paste0(sheet, "_", clean_gene_set, "_enhancers.bed")
  tss_file <- paste0(sheet, "_", clean_gene_set, "_TSS.bed")
  
  write.table(
    enhancers,
    file = enhancer_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  write.table(
    tss,
    file = tss_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  invisible(list(
    genes = genes,
    peaks = peaks,
    enhancers = enhancers,
    tss = tss
  ))
}


# Example usage: 

write_pathway_beds(
  sheet = "H3K27ac",
  collection = "C5",
  gene_set = "GOBP_T_CELL_ACTIVATION"
)

# To reproduce the manuscript figures, use the following combinations:
# Histone marks:
#   - H3K4me3
#   - H3K4me1
#   - H3K27ac
#   - H3K27me1
#
# Gene sets:
#   - KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY
#   - WP_TCELL_RECEPTOR_SIGNALING
#   - GOBP_T_CELL_ACTIVATION (Supplementary Figure 1c–d.)
