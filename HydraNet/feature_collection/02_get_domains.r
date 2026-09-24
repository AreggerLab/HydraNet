#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: 02_get_domains.R <input_file> <output_file>")

INPUT_CSV  <- args[1]
OUTPUT_CSV <- args[2]

suppressPackageStartupMessages({
  library(AnnotationHub)
  library(ensembldb)
  library(dplyr)
  library(readr)
})

message("🔹 Loading EnsDb v104...")
# --- TRAVER ANNOTATIONHUB CACHE FIX ---
AH_CACHE <- Sys.getenv(
    "ANNOTATION_HUB_CACHE",
    unset = "path/AnnotationHub_cache"
)

dir.create(
    AH_CACHE,
    recursive = TRUE,
    showWarnings = FALSE
)

AnnotationHub::setAnnotationHubOption(
    "CACHE",
    AH_CACHE
)

message(paste0("AnnotationHub cache: ", AH_CACHE))
# --- END TRAVER ANNOTATIONHUB CACHE FIX ---

ah <- AnnotationHub()
edb <- query(ah, c("EnsDb", "Homo sapiens", "GRCh38", "104"))[[1]]
message("  - ✅ EnsDb loaded.")

df <- read_csv(INPUT_CSV, show_col_types = FALSE, guess_max = 15000)

# Identify exonic guides
df <- df %>% mutate(original_row_idx = row_number())
exonic_guides <- df %>% filter(!is.na(protein_id))

if (nrow(exonic_guides) > 0) {
  message(paste0("🔹 Found ", nrow(exonic_guides), " exonic guides."))

  df$Binary_Domain[exonic_guides$original_row_idx] <- 0L

  prot_ids <- unique(exonic_guides$protein_id)
  dom_all <- proteins(edb, filter = ProteinIdFilter(prot_ids),
                      columns = c("protein_id", "protein_domain_source", "prot_dom_start", "prot_dom_end"))

  if (length(dom_all) > 0) {
    dom_all <- dom_all[dom_all$protein_domain_source %in%
                         c("smart", "pfam", "scanprosite", "gene3d", "mobidb"), ]
    dom_df <- merge(select(exonic_guides, original_row_idx, protein_id, aa_index),
                    as.data.frame(dom_all), by = "protein_id") %>%
              filter(!is.na(aa_index) &
                     prot_dom_start <= aa_index + 2 &
                     prot_dom_end >= aa_index - 2)

    if (nrow(dom_df) > 0) {
      df$Binary_Domain[unique(dom_df$original_row_idx)] <- 1L
      message(paste0("  - ✅ Domains assigned to ", length(unique(dom_df$original_row_idx)), " guides."))
    } else {
      message("  - ⚠️ No domain hits found.")
    }
  } else {
    message("  - ⚠️ No protein domain records found.")
  }
} else {
  message("  - ⚠️ No exonic guides found.")
}

write_csv(df, OUTPUT_CSV)
message(paste0("🎉 DONE: Wrote output to ", OUTPUT_CSV))


