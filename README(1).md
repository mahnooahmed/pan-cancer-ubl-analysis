# Pan-Cancer UBL Pathway Mutational Landscape

**Author:** Syeda Mahnoor Ahmed  
**Affiliation:** Epigenetics and Genome Integrity Lab (EaGIL), LUMS  
**Contact:** ahmedmahnoor818@gmail.com

---

## Overview

This repository contains R scripts developed during my Research Assistantship at EaGIL, LUMS, as part of a pan-cancer project mapping the mutational landscape of the **ubiquitin and ubiquitin-like (UBL) pathway**, specifically the information flow from **E1 → E2 → E3 → Ubiquitin → UBLs** across 28 TCGA cancer types.

This script (`pan_cancer_oncoprint.R`) processes TCGA mutation data (MAF format) for any UBL gene category, **E1, E2, E3, DUBs**, or others by setting a single variable at the top of the script. It generates:
- Pan-cancer oncoplots with TMB annotation
- Gene-level mutation summaries across cancer types
- Patient-level mutation burden summaries

---

## Repository Structure

```
├── pan_cancer_oncoprint_E1.R     # Main oncoprint pipeline (this script)
├── README.md
└── .gitignore
```

---

## Requirements

### R version
R >= 4.1.0

### CRAN packages
```r
install.packages(c("dplyr", "data.table", "purrr", "tidyr", "openxlsx", "readr"))
```

### Bioconductor packages
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("maftools", "ComplexHeatmap"))
```

---

## Data Requirements

This script requires **TCGA MAF files** which are not included in this repository due to size and data use restrictions. Data can be accessed via:
- [TCGA GDC Data Portal](https://portal.gdc.cancer.gov/)
- [cBioPortal](https://www.cbioportal.org/)

### Expected folder structure
```
TCGA_clean_data/
├── BRCA/
│   └── Extracted_XLSX_Files/
│       ├── sample1.xlsx
│       └── sample2.xlsx
├── LUAD/
│   └── Extracted_XLSX_Files/
│       └── ...
└── ...
```

Each Excel file must contain at minimum these MAF columns:
`Hugo_Symbol`, `Entrez_Gene_Id`, `Center`, `NCBI_Build`, `Chromosome`,
`Start_Position`, `End_Position`, `Strand`, `Variant_Classification`,
`Variant_Type`, `Reference_Allele`, `Tumor_Seq_Allele1`, `Tumor_Seq_Allele2`,
`Tumor_Sample_Barcode`, `Matched_Norm_Sample_Barcode`

---

## Usage

1. Clone this repository
2. Install dependencies (see above)
3. Update the settings at the top of the script:

```r
base_dir      <- "path/to/TCGA_clean_data"
gene_dir      <- "path/to/Gene_List"
out_dir       <- "path/to/Oncoprint_Outputs"
GENE_CATEGORY <- "E1"   # change to "E2", "E3", "DUB", etc.
capture_mb    <- 38
```

4. Run the script:
```r
source("pan_cancer_oncoprint_E1.R")
```

---

## Output Files

| File | Description |
|------|-------------|
| `combined_pan_cancer_E1.maf` | Combined MAF file across all cancer types |
| `pan_cancer_tmb_per_sample_E1.tsv` | Per-sample tumour mutational burden |
| `Pan_Cancer_Patient_Mutation_Summary_E1.xlsx` | Patient-level mutation counts |
| `UBL_E1_Genes_Gene_Summary_PanCanceR_E1.xlsx` | Gene-level summary across cancers |
| `PanCancer_Oncoprint_UBL_E1_Genes_E1_Genes.pdf` | Publication-quality oncoprint (PDF) |
| `PanCancer_Oncoprint_UBL_E1_Genes_E1_Genes.tiff` | High-resolution oncoprint (600 DPI TIFF) |

---

## Notes

- TMB is calculated assuming a WES capture size of **38 Mb** (adjustable via `capture_mb`)
- Oncoprint displays top 50 most-mutated genes by default (adjustable via `topN` parameter)
- Samples with no mutations in the gene list are excluded from the oncoprint
- Multi-hit samples are flagged as `multiple`

---

## Citation / Acknowledgements

If you use or adapt this code, please acknowledge:

> Syeda Mahnoor Ahmed, EaGIL, LUMS (2026). Pan-Cancer UBL Pathway Mutational Landscape Pipeline. GitHub: https://github.com/mahnooahmed

Data sourced from [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/tcga).

---

## License

This code is shared for academic and research purposes. Please contact the author before use in publications or derivative works.
