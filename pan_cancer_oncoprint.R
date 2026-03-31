# =============================================================================
# Pan-Cancer UBL Oncoprint: E1 Genes
# =============================================================================
# Description:
#   Processes TCGA mutation data across all cancer types and generates
#   oncoplots for E1 ubiquitin-activating enzyme genes using non-overlapping
#   gene lists. Outputs include TMB summaries, gene-level mutation summaries,
#   and publication-quality oncoprint figures (PDF + TIFF).
#
# Author:    Syeda Mahnoor Ahmed
#            Research Assistant, Epigenetics and Genome Integrity Lab (EaGIL)
#            Lahore University of Management Sciences (LUMS)
# Project:   Pan-Cancer UBL Pathway Mutational Landscape
# Date:      2026
# Contact:   ahmedmahnoor818@gmail.com
#
# Input:
#   - TCGA MAF files (Excel format) organised by cancer type folder
#   - Gene list Excel file with sheet "COMBINED GENES", column "E1"
#
# Output:
#   - combined_pan_cancer_E1.maf
#   - pan_cancer_tmb_per_sample_E1.tsv
#   - Pan_Cancer_Patient_Mutation_Summary_E1.xlsx
#   - UBL_E1_Genes_Gene_Summary_PanCanceR_E1.xlsx
#   - PanCancer_Oncoprint_UBL_E1_Genes_E1_Genes.pdf
#   - PanCancer_Oncoprint_UBL_E1_Genes_E1_Genes.tiff
#
# Dependencies:
#   CRAN:   dplyr, data.table, purrr, tidyr, openxlsx, readr
#   Bioc:   maftools, ComplexHeatmap
# =============================================================================


# ---- 0) Dependencies --------------------------------------------------------

quiet_install <- function(pkgs) {
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
}
quiet_install(c("dplyr", "data.table", "purrr", "tidyr", "openxlsx", "readr"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("maftools",      quietly = TRUE)) BiocManager::install("maftools",      ask = FALSE)
if (!requireNamespace("ComplexHeatmap",quietly = TRUE)) BiocManager::install("ComplexHeatmap", ask = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(purrr)
  library(tidyr)
  library(maftools)
  library(ComplexHeatmap)
  library(grid)
  library(openxlsx)
  library(readr)
})


# ---- 1) User-defined Settings -----------------------------------------------
# UPDATE THESE to match your local directory structure and target gene category.

base_dir      <- "path/to/TCGA_clean_data"       # folder containing per-cancer-type subdirs
gene_dir      <- "path/to/Gene_List"             # folder containing the gene list Excel file
out_dir       <- "path/to/Oncoprint_Outputs"     # folder where all outputs will be saved
GENE_CATEGORY <- "E1"                            # column to pull from gene list: "E1", "E2", "E3", "DUB", etc.
capture_mb    <- 38                              # WES capture size in megabases for TMB calculation

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- 2) Load Non-Overlapping Gene Lists -------------------------------------

cat("\n=== Loading Non-Overlapping Gene Lists ===\n")

gene_file  <- file.path(gene_dir, "Ub Gene Verified List (w_ AI) (1).xlsx")
gene_sheet <- read.xlsx(gene_file, sheet = "COMBINED GENES")

gene_list <- unique(na.omit(trimws(gene_sheet[[GENE_CATEGORY]])))
gene_list <- gene_list[gene_list != ""]
cat("Total unique", GENE_CATEGORY, "genes loaded:", length(gene_list), "\n")


# ---- 3) Scan Cancer Type Folders --------------------------------------------

cat("\n=== Scanning Cancer Type Folders ===\n")

all_folders    <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
cancer_folders <- all_folders[dir.exists(file.path(all_folders, "Extracted_XLSX_Files"))]

cat("Found", length(cancer_folders), "cancer type folders:\n")
for (f in cancer_folders) cat("  -", basename(f), "\n")


# ---- 4) Read All Excel MAF Files Across Cancer Types ------------------------

cat("\n=== Reading MAF Files from All Cancer Types ===\n")

min_cols <- c(
  "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build",
  "Chromosome", "Start_Position", "End_Position", "Strand",
  "Variant_Classification", "Variant_Type",
  "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
  "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"
)

read_one_excel <- function(fp) {
  cat("  Reading:", basename(fp), "\n")
  df <- tryCatch(
    read.xlsx(fp, sheet = 1),
    error = function(e) {
      cat("    Could not read", basename(fp), "\n")
      return(NULL)
    }
  )
  if (is.null(df)) return(NULL)
  if (!all(min_cols %in% names(df))) {
    cat("    Skipping - missing required columns\n")
    return(NULL)
  }
  df[, min_cols]
}

all_maf_data <- list()

for (cancer_folder in cancer_folders) {
  cancer_type <- basename(cancer_folder)
  cat("\nProcessing:", cancer_type, "\n")

  xlsx_dir    <- file.path(cancer_folder, "Extracted_XLSX_Files")
  excel_files <- list.files(xlsx_dir, pattern = "\\.xlsx$|\\.xls$",
                            full.names = TRUE, recursive = FALSE)

  if (length(excel_files) == 0) {
    cat("  No Excel files found in", xlsx_dir, "\n")
    next
  }

  cat("  Found", length(excel_files), "Excel files\n")
  raw_df_list_clean <- list()

  for (fp in excel_files) {
    df <- read_one_excel(fp)
    if (!is.null(df)) {
      df[]           <- lapply(df, as.character)
      df$Cancer_Type <- cancer_type
      raw_df_list_clean[[length(raw_df_list_clean) + 1]] <- df
    }
    rm(df); gc()
  }

  if (length(raw_df_list_clean) > 0) {
    all_maf_data[[cancer_type]] <- bind_rows(raw_df_list_clean)
    cat("  Loaded", nrow(all_maf_data[[cancer_type]]), "mutations from", cancer_type, "\n")
  }
}

raw_df <- bind_rows(all_maf_data)
rm(all_maf_data); gc()

cat("\n=== Combined Data Summary ===\n")
cat("Total mutations loaded:", nrow(raw_df), "\n")
cat("Unique samples:",         length(unique(raw_df$Tumor_Sample_Barcode)), "\n")
cat("Cancer types:",           length(unique(raw_df$Cancer_Type)), "\n")


# ---- 5) Sample-to-Cancer Mapping --------------------------------------------

cat("\n=== Creating Sample-Cancer Type Mapping ===\n")

sample_cancer_map <- raw_df %>%
  select(Tumor_Sample_Barcode, Cancer_Type) %>%
  distinct()

cat("Sample-Cancer mapping created:", nrow(sample_cancer_map), "samples\n")


# ---- 6) Build MAF Object ----------------------------------------------------

cat("\n=== Creating MAF Object ===\n")

maf_path <- file.path(out_dir, paste0("combined_pan_cancer_", GENE_CATEGORY, ".maf"))
write.table(raw_df, maf_path, sep = "\t", quote = FALSE, row.names = FALSE)

m <- read.maf(maf = maf_path)
print(m)


# ---- 7) Compute Tumour Mutational Burden (TMB) ------------------------------

cat("\n=== Computing TMB ===\n")

tmb_df     <- tmb(maf = m, captureSize = capture_mb, logScale = FALSE)

if ("total_perMB" %in% names(tmb_df)) {
  tmb_vec <- setNames(tmb_df$total_perMB, tmb_df$Tumor_Sample_Barcode)
} else if ("total" %in% names(tmb_df)) {
  tmb_vec <- setNames(tmb_df$total, tmb_df$Tumor_Sample_Barcode)
} else {
  stop("tmb() output missing both 'total_perMB' and 'total' columns.")
}

write.table(tmb_df, file.path(out_dir, paste0("pan_cancer_tmb_per_sample_", GENE_CATEGORY, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("TMB file saved\n")

patient_summary <- m@data %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    Total_Mutations = n(),
    Unique_Genes    = n_distinct(Hugo_Symbol),
    UBL_Mutations   = sum(Hugo_Symbol %in% gene_list),
    .groups = "drop"
  ) %>%
  left_join(sample_cancer_map, by = "Tumor_Sample_Barcode") %>%
  arrange(desc(Total_Mutations))

write.xlsx(patient_summary, file.path(out_dir, paste0("Pan_Cancer_Patient_Mutation_Summary_", GENE_CATEGORY, ".xlsx")))
cat("Patient mutation summary saved\n")


# ---- 8) Gene-Level Mutation Summary -----------------------------------------

cat("\n=== Creating Gene Mutation Summary ===\n")

create_gene_summary <- function(gene_list, category_name) {
  maf_df <- as_tibble(m@data)

  # Drop Cancer_Type if already present to avoid .x/.y suffix collision
  if ("Cancer_Type" %in% colnames(maf_df)) maf_df <- maf_df %>% select(-Cancer_Type)

  cancer_map  <- as_tibble(sample_cancer_map)
  joined_data <- maf_df %>%
    filter(Hugo_Symbol %in% gene_list) %>%
    left_join(cancer_map, by = "Tumor_Sample_Barcode")

  gene_summary <- joined_data %>%
    group_by(Hugo_Symbol) %>%
    summarise(
      Total_Samples_Mutated = n_distinct(Tumor_Sample_Barcode),
      Total_Mutations       = n(),
      Cancer_Types_Affected = n_distinct(Cancer_Type, na.rm = TRUE),
      Cancer_Types          = paste(sort(unique(Cancer_Type[!is.na(Cancer_Type)])), collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(desc(Total_Samples_Mutated))

  write.xlsx(gene_summary,
             file.path(out_dir, paste0(category_name, "_Gene_Summary_PanCanceR_E1.xlsx")))
  cat(category_name, "gene summary saved\n")
  return(gene_summary)
}


# ---- 9) Oncoprint Generation ------------------------------------------------

create_oncoplot <- function(gene_list, category_name, topN = 50) {

  cat("\n=== Creating", category_name, "Oncoplot ===\n")

  # Mutation classification mapping
  mut_map <- c(
    Missense_Mutation                          = "missense",
    Nonsense_Mutation                          = "nonsense",
    Frame_Shift_Del                            = "frameshift",
    Frame_Shift_Ins                            = "frameshift",
    Splice_Site                                = "splicing",
    Nonstop_Mutation                           = "nonstop",
    In_Frame_Del                               = "inframe",
    In_Frame_Ins                               = "inframe"
  )

  d_filtered <- as_tibble(m@data) %>%
    filter(Hugo_Symbol %in% gene_list) %>%
    mutate(alt = recode(Variant_Classification, !!!mut_map, .default = NA_character_)) %>%
    filter(!is.na(alt))

  # Rank genes by number of mutated samples
  gene_rank <- d_filtered %>%
    group_by(Hugo_Symbol) %>%
    summarise(MutatedSamples = n_distinct(Tumor_Sample_Barcode), .groups = "drop") %>%
    arrange(desc(MutatedSamples))

  top_genes <- head(gene_rank$Hugo_Symbol, min(topN, nrow(gene_rank)))
  cat("Top genes selected:", length(top_genes), "\n")
  cat("Most mutated gene:", top_genes[1], "—", gene_rank$MutatedSamples[1], "samples\n")

  # Build oncoprint matrix
  all_samples <- sort(unique(m@data$Tumor_Sample_Barcode))
  n_total     <- length(all_samples)

  oncomat <- matrix("", nrow = length(top_genes), ncol = n_total,
                    dimnames = list(top_genes, all_samples))

  d2   <- d_filtered %>% filter(Hugo_Symbol %in% top_genes)
  cell <- d2 %>%
    group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    summarise(state = paste(sort(unique(alt)), collapse = ";"), .groups = "drop")

  oncomat[cbind(
    match(cell$Hugo_Symbol,        rownames(oncomat)),
    match(cell$Tumor_Sample_Barcode, colnames(oncomat))
  )] <- cell$state

  idx_multi          <- which(oncomat != "" & grepl(";", oncomat, fixed = TRUE))
  oncomat[idx_multi] <- paste0(oncomat[idx_multi], ";multiple")

  # Keep only samples with at least one mutation
  mutated_samples <- colnames(oncomat)[colSums(oncomat != "") > 0]
  oncomat         <- oncomat[, mutated_samples, drop = FALSE]
  n_mutated       <- length(mutated_samples)
  cat("Samples with mutations:", n_mutated, "out of", n_total, "\n")

  tmb_vec_filtered                         <- tmb_vec[mutated_samples]
  tmb_vec_filtered[is.na(tmb_vec_filtered)] <- 0

  # Waterfall sort by total mutations per sample
  waterfall_order <- names(sort(colSums(oncomat != ""), decreasing = TRUE))
  oncomat         <- oncomat[, waterfall_order, drop = FALSE]
  tmb_for_plot    <- tmb_vec_filtered[waterfall_order]

  # Sort genes by mutated sample count
  gene_order <- names(sort(rowSums(oncomat != ""), decreasing = TRUE))
  oncomat    <- oncomat[gene_order, , drop = FALSE]

  # Colour scheme
  stacked_bar_colors <- c(
    missense   = "#97BF65",
    nonsense   = "#E41A1C",
    frameshift = "#377EB8",
    splicing   = "#FFB74D",
    nonstop    = "#FFF176",
    inframe    = "#A65628",
    multiple   = "#999999",
    WT         = "#FFFFFF"
  )
  col_alter <- stacked_bar_colors

  # Top annotation: TMB barplot
  top_anno <- HeatmapAnnotation(
    TMB = anno_barplot(
      tmb_for_plot,
      axis = FALSE, border = FALSE,
      height = unit(2, "cm"),
      gp = gpar(fill = "grey70", col = NA)
    ),
    show_annotation_name = FALSE
  )

  # Right annotation: stacked mutation-type barplot + mutation %
  gene_type_counts <- d2 %>%
    distinct(Hugo_Symbol, Tumor_Sample_Barcode, alt) %>%
    filter(Tumor_Sample_Barcode %in% colnames(oncomat)) %>%
    count(Hugo_Symbol, alt, name = "n") %>%
    pivot_wider(names_from = alt, values_from = n, values_fill = 0) %>%
    as.data.frame()

  rownames(gene_type_counts) <- gene_type_counts$Hugo_Symbol
  gene_type_counts <- gene_type_counts[
    match(rownames(oncomat), rownames(gene_type_counts)),
    setdiff(colnames(gene_type_counts), "Hugo_Symbol"),
    drop = FALSE
  ]
  gene_type_counts[is.na(gene_type_counts)] <- 0

  mut_pct    <- round(100 * rowSums(oncomat != "") / n_total)
  xlim_ticks <- pretty(c(0, ceiling(max(rowSums(gene_type_counts)) * 1.05)), n = 5)
  xlim_max   <- max(xlim_ticks)

  right_anno <- rowAnnotation(
    `# samples` = anno_barplot(
      as.matrix(gene_type_counts),
      gp = gpar(fill = stacked_bar_colors[colnames(gene_type_counts)], col = NA),
      width  = unit(50, "mm"), border = FALSE, axis = TRUE,
      axis_param = list(side = "top", at = xlim_ticks,
                        labels = as.character(xlim_ticks)),
      xlim = c(0, xlim_max)
    ),
    `%` = anno_text(paste0(mut_pct, "%"),
                    gp = gpar(fontsize = 8), just = "left", location = 0.05),
    annotation_name_side = "top",
    annotation_name_rot  = 0,
    annotation_name_gp   = gpar(fontsize = 10)
  )

  # Alteration drawing functions
  alter_fun <- list(
    background = function(x, y, w, h) grid.rect(x, y, w, h,
                                                  gp = gpar(fill = "#FFFFFF", col = NA)),
    missense   = function(x, y, w, h) grid.rect(x, y, w*0.95, h*0.95,
                                                  gp = gpar(fill = col_alter["missense"],   col = NA)),
    nonsense   = function(x, y, w, h) grid.rect(x, y, w*0.95, h*0.95,
                                                  gp = gpar(fill = col_alter["nonsense"],   col = NA)),
    frameshift = function(x, y, w, h) grid.rect(x, y, w*0.95, h*0.95,
                                                  gp = gpar(fill = col_alter["frameshift"], col = NA)),
    splicing   = function(x, y, w, h) grid.rect(x, y, w*0.95, h*0.95,
                                                  gp = gpar(fill = col_alter["splicing"],   col = NA)),
    nonstop    = function(x, y, w, h) grid.rect(x, y, w*0.95, h*0.95,
                                                  gp = gpar(fill = col_alter["nonstop"],    col = NA)),
    inframe    = function(x, y, w, h) grid.rect(x, y, w*0.95, h*0.95,
                                                  gp = gpar(fill = col_alter["inframe"],    col = NA)),
    multiple   = function(x, y, w, h) grid.rect(x + w*0.32, y, w*0.18, h*0.18,
                                                  gp = gpar(fill = col_alter["multiple"],   col = NA))
  )

  alter_levels <- c("missense","nonsense","frameshift","splicing","nonstop","inframe","multiple")
  lgd_alter    <- Legend(title = "Alterations", at = alter_levels,
                         legend_gp = gpar(fill = col_alter[alter_levels]))
  lgd_wt       <- Legend(title = NULL, at = "WT",
                         legend_gp = gpar(fill = "#FFFFFF"), labels = "WT")

  plot_title <- paste0(
    "UBL Pan-Cancer Oncoprint: ", GENE_CATEGORY, " Genes (n = ", n_total, ", ", n_mutated, " mutated)"
  )

  draw_tmb_axis <- function() {
    tmb_max <- max(tmb_for_plot, na.rm = TRUE)
    if (is.infinite(tmb_max) || is.na(tmb_max)) tmb_max <- 0
    ticks  <- pretty(c(0, tmb_max))
    ticks  <- ticks[ticks >= 0 & ticks <= tmb_max]
    if (length(ticks) == 0) ticks <- c(0, tmb_max)
    tick_y <- if (tmb_max > 0) ticks / tmb_max else seq(0, 1, length.out = length(ticks))

    decorate_annotation("TMB", {
      grid::grid.lines(x = unit(rep(-1.5, 2), "mm"), y = unit(c(0, 1), "npc"),
                       gp = gpar(lwd = 0.5))
      grid::grid.text(label = ticks, x = unit(-6, "mm"), y = unit(tick_y, "npc"),
                      gp = gpar(fontsize = 8))
      grid::grid.text("TMB", x = unit(-10.5, "mm"), y = unit(0.5, "npc"),
                      rot = 90, gp = gpar(fontsize = 10, fontface = "bold"))
    })
  }

  ht <- oncoPrint(
    oncomat,
    alter_fun         = alter_fun,
    col               = col_alter,
    column_order      = colnames(oncomat),
    row_order         = rownames(oncomat),
    top_annotation    = top_anno,
    right_annotation  = right_anno,
    show_pct          = FALSE,
    remove_empty_columns     = FALSE,
    show_column_names = FALSE,
    show_row_names    = TRUE,
    row_names_side    = "left",
    row_names_gp      = gpar(fontsize = 9),
    column_title      = plot_title,
    alter_fun_is_vectorized = FALSE,
    show_heatmap_legend     = FALSE
  )

  # Save PDF
  pdf_file <- file.path(out_dir, paste0("PanCancer_Oncoprint_", category_name, "_E1_Genes.pdf"))
  pdf(pdf_file, width = 16, height = 10)
  draw(ht, merge_legends = TRUE, heatmap_legend_list = list(lgd_alter, lgd_wt))
  draw_tmb_axis()
  dev.off()
  cat("PDF saved:", pdf_file, "\n")

  # Save TIFF
  tiff_file <- file.path(out_dir, paste0("PanCancer_Oncoprint_", category_name, "_E1_Genes.tiff"))
  tiff(tiff_file, width = 16, height = 10, units = "in", res = 600, compression = "lzw")
  draw(ht, merge_legends = TRUE, heatmap_legend_list = list(lgd_alter, lgd_wt))
  draw_tmb_axis()
  dev.off()
  cat("TIFF saved:", tiff_file, "\n")
}


# ---- 10) Run Pipeline -------------------------------------------------------

ubl_summary <- create_gene_summary(gene_list, paste0("UBL_", GENE_CATEGORY, "_Genes"))
create_oncoplot(gene_list, paste0("UBL_", GENE_CATEGORY, "_Genes"))
