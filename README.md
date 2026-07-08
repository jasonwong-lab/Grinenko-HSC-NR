# Grinenko-HSC-NR
Code repository for Liu et al. Epigenetic memory in aged hematopoietic stem cells governs anti-tumor T cell immunity.

### preprocessing_commands.sh
Important commands use to preprocess the following: 
1. RNA-seq data for HSC/CD4/CD8 
2. ATACseq data for HSCs
3. EMseq data for HSCs

All ATAC-seq, EM-seq, and RNA-seq raw data are deposited in Gene Expression Omnibus under the accession number GSE335321. 

### Figures_Scripts
Scripts used to generate main figures.

### Data sources for fig2.R — fig6.R
**All input files** required by the figure-generation R scripts are produced by preprocessing_commands.sh from the raw data in GSE335321. Run preprocessing_commands.sh first, then execute the R scripts using the generated outputs.
All reference files needed to run these scripts are highlighted in the section below. 

### Data Sources for Figure 1
This script uses publicly available datasets from the following studies:

#### Dataset 1 – Human Naïve CD4 and CD8 cell 

- **Publication:** *Multidimensional profiling of human T cells reveals high CD38 expression, marking recent thymic emigrants and age-related naive T cell remodeling* (Bohacova et al., 2024). 
- **Required supplementary file:** `NIHMS2160335-supplement-mmc7.xlsx`.
- Bohacova P, Terekhova M, Tsurinov P, Mullins R, Husarcikova K, Shchukina I, Antonova AU, Echalar B, Kossl J, Saidu A, Francis T, Mannie C, Arthur L, Harridge SDR, Kreisel D, Mudd PA, Taylor AM, McNamara CA, Cella M, Puram SV, van den Broek T, van Wijk F, Eghtesady P, Artyomov MN. (2024). Multidimensional profiling of human T cells reveals high CD38 expression, marking recent thymic emigrants and age-related naive T cell remodeling. Immunity, 57(10), 2362–2379.e10. https://doi.org/10.1016/j.immuni.2024.08.019

#### Dataset 2 – Chromatin Profiles of Purified Human HSCs

- **Publication:** *Aging Human Hematopoietic Stem Cells Manifest Profound Epigenetic Reprogramming of Enhancers That May Predispose to Leukemia* (Adelman et al., 2019).
- **Required supplementary file:** Supplementary Table S2 (differential peak calling results) `21598290cd181474-sup-213895_2_supp_5514125_pr2m22.xlsx`
- **ChIP-seq data:** GEO Series **GSE104408**.
- Adelman ER, Huang HT, Roisman A, Olsson A, Colaprico A, Qin T, Lindsley RC, Bejar R, Salomonis N, Grimes HL, Figueroa ME. (2019). Aging Human Hematopoietic Stem Cells Manifest Profound Epigenetic Reprogramming of Enhancers That May Predispose to Leukemia. Cancer Discovery, 9(8), 1080–1101. https://doi.org/10.1158/2159-8290.CD-18-1474

### Reference files
The following reference files are required by one or more scripts in this repository.

| File | Description | Source |
|------|-------------|--------|
| `Mus_musculus.GRCm38.102.gtf` | Ensembl gene annotation (release 102) | Ensembl |
| `mm10.fa` | Mouse reference genome sequence (mm10/GRCm38) | UCSC Genome Browser |
| `mm9ToMm10.over.chain` | UCSC liftOver chain file for converting mm9 coordinates to mm10 | UCSC Genome Browser |
| `cpgIslandExt.txt.gz` | UCSC CpG island annotation for mm10 | UCSC Table Browser | 
| `hg38.chrom.sizes` | UCSC utilities | UCSC Genome Browser |
| `reMap2022.bb` | reMap 2022 transcriptional regulator binding site annotations | reMap 2022 database |
| `ncbiRefSeqSelect.noBin.gtf` | UCSC RefSeq Select table for genomic feature annotation | UCSC Table Browser |
| `Enhancers.xls` | Enhancer annotation derived from published HSC chromatin state maps | Lara-Astiaso D, Weiner A, Lorenzo-Vivas E, *et al.* (2014). **Chromatin state dynamics during blood formation.** *Science* 345(6199):943–949. DOI: 10.1126/science.1256271. Enhancer annotations were obtained from Supplementary Table `1256271tables1.xls` |

### Notes on Usage for fig1.R and fig1.sh

`fig1.R` is divided into three independent sections.

#### Part 1 – Gene Set Enrichment Analysis (GSEA)

The script contains functions for running GSEA. To reproduce the figures, specify the desired pathway:

```r
gene_set <- get_gene_set(pathway)
```

Available pathways:

| Pathway                                  | Figure                    |
| ---------------------------------------- | ------------------------- |
| `KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY` | Figure 1a–b               |
| `WP_TCELL_RECEPTOR_SIGNALING`            | Figure 1a–b               |
| `GOBP_T_CELL_ACTIVATION`                 | Supplementary Figure 1a–b |

---

#### Part 2 – Over-representation analysis (ORA)

The script generates gene lists for over-representation analysis. These lists should be uploaded to **Metascape** (https://metascape.org) to perform ORA and generate the corresponding enrichment results.

---

#### Part 3 – Generation of pathway BED files

The script generates BED files for downstream ChIP-seq signal profiling.

Example:

```r
write_pathway_beds(
  sheet = "H3K27ac",
  collection = "C5",
  gene_set = "GOBP_T_CELL_ACTIVATION"
)
```

Supported histone marks:

* `H3K4me3`
* `H3K4me1`
* `H3K27ac`
* `H3K27me1`

Supported pathways:

* `GOBP_T_CELL_ACTIVATION`
* `KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY`
* `WP_TCELL_RECEPTOR_SIGNALING`

The generated BED files are used as input for `fig1.sh` to generate the profile plots (Supplementary Figure 1c–d).
