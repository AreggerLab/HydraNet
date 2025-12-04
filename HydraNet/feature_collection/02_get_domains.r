#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: 02_get_domains.R <input_file> <output_file>")

INPUT_CSV  <- args[1]
OUTPUT_CSV <- args[2]

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(biomaRt)
  library(tidyr)
})

df <- read_csv(INPUT_CSV, show_col_types = FALSE, guess_max = 15000)

# Track original row index
df <- df %>% mutate(original_row_idx = dplyr::row_number())

# Exonic guides = those with a valid protein_id and aa_index > 0
exonic_guides <- df %>%
  dplyr::filter(!is.na(protein_id),
                !is.na(aa_index),
                aa_index > 0)

if (nrow(exonic_guides) > 0) {

  # Initialize as 0 for exonic guides (NA elsewhere)
  df$Binary_Domain <- NA_integer_
  df$Binary_Domain[exonic_guides$original_row_idx] <- 0L

  # Connect to Ensembl mouse, mm10 / GRCm38 archive
  mart <- useMart(
    "ENSEMBL_MART_ENSEMBL",
    dataset = "mmusculus_gene_ensembl",
    host    = "https://may2021.archive.ensembl.org"
  )
  #========For Human========#
  #mart <- useMart(
  #"ENSEMBL_MART_ENSEMBL",
  #dataset = "hsapiens_gene_ensembl",
  #host    = "https://may2021.archive.ensembl.org"
  #)

  # Pull PFAM + InterPro domain ranges for relevant proteins
  domains <- getBM(
    attributes = c(
      "ensembl_peptide_id",
      "pfam", "pfam_start", "pfam_end",
      "interpro", "interpro_start", "interpro_end"
    ),
    filters  = "ensembl_peptide_id",
    values   = unique(exonic_guides$protein_id),
    mart     = mart
  )

  # Build PFAM table
  pfam_long <- domains %>%
    dplyr::filter(!is.na(pfam),
                  !is.na(pfam_start),
                  !is.na(pfam_end)) %>%
    dplyr::transmute(
      ensembl_peptide_id,
      domain_source  = "pfam",
      domain_id      = pfam,
      prot_dom_start = as.numeric(pfam_start),
      prot_dom_end   = as.numeric(pfam_end)
    )

  # Build InterPro table
  interpro_long <- domains %>%
    dplyr::filter(!is.na(interpro),
                  !is.na(interpro_start),
                  !is.na(interpro_end)) %>%
    dplyr::transmute(
      ensembl_peptide_id,
      domain_source  = "interpro",
      domain_id      = interpro,
      prot_dom_start = as.numeric(interpro_start),
      prot_dom_end   = as.numeric(interpro_end)
    )

  # Combined domain table
  dom_long <- dplyr::bind_rows(pfam_long, interpro_long)

  if (nrow(dom_long) > 0) {
    # Merge with exonic guides and check AA overlap (±2 AA around aa_index)
    dom_df <- merge(
      dplyr::select(exonic_guides, original_row_idx, protein_id, aa_index),
      dom_long,
      by.x = "protein_id",
      by.y = "ensembl_peptide_id"
    ) %>%
      dplyr::filter(
        !is.na(prot_dom_start),
        !is.na(prot_dom_end),
        aa_index >= (prot_dom_start - 2),
        aa_index <= (prot_dom_end + 2)
      )

    if (nrow(dom_df) > 0) {
      df$Binary_Domain[unique(dom_df$original_row_idx)] <- 1L
    }
  }
}

message(paste0("🐉🐉🐉 DONE: Writing ", OUTPUT_CSV))
write_csv(df, OUTPUT_CSV)

