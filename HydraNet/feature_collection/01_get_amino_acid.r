#!/usr/bin/env Rscript

# --- ARGS ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 01_get_amino_acid_duncan.r <input_guide_file> <output_csv>", call. = FALSE)
}
INPUT_FILE <- args[1]
OUTPUT_CSV <- args[2]

# --- LIBS ---
suppressPackageStartupMessages({
  library(AnnotationHub)
  library(ensembldb)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(dplyr)
  library(readr)
  library(purrr)
  library(memoise)
  library(tidyr)
  library(stringr)
})

message("🐉🐉🐉 Loading EnsDb for Mus musculus (GRCm38/mm10) from AnnotationHub...")#====GRCh38 for human=====#
ah <- AnnotationHub()

qry <- query(ah, c("EnsDb", "Mus musculus", "GRCm38"))
#qry <- query(ah, c("EnsDb", "Homo sapiens", "GRCh38"))#=====For Human====#
if (length(qry) == 0) {
  stop("No Mus musculus (GRCm38) EnsDb found in AnnotationHub. Try updating Bioconductor/AnnotationHub.", call. = FALSE)
}

if (!is.null(mcols(qry)$rdatadateadded)) {
  qry <- qry[order(mcols(qry)$rdatadateadded, decreasing = TRUE)]
}

edb <- qry[[1]]
message(sprintf("  -  EnsDb loaded: %s", mcols(qry)$title[1]))



# --- READ ---
message(paste("Reading:", basename(INPUT_FILE)))

df <- read_csv(INPUT_FILE, show_col_types = FALSE)


# --- Prepare data for annotation ---
message("🐉🐉🐉 Preparing data for annotation...")
df <- df %>%
  # Ensure the chromosome column is named correctly for GRanges
  rename(chromosome = seqnames) %>%
  # Ensure essential columns have the right type
  mutate(
    chromosome = as.character(chromosome),
    cut_site = as.integer(cut_site)
  )

# --- ANNOTATE BY EXON OVERLAP USING cut_site ---
message("🐉🐉🐉 Annotating guides with gene/exon information...")
# Filter out rows where cut_site is NA, as they cannot be used for annotation
annot_df <- df %>%
  filter(!is.na(cut_site) & !is.na(chromosome)) %>%
  mutate(original_row_idx = row_number()) # Use a temporary index for joining

if (nrow(annot_df) > 0) {
    guides_gr <- GRanges(
      seqnames = annot_df$chromosome,
      ranges = IRanges(start = annot_df$cut_site, end = annot_df$cut_site),
      original_row_idx = annot_df$original_row_idx
    )
    seqlevelsStyle(guides_gr) <- "Ensembl"
    all_exons_gr <- exons(
      edb,
      columns = c("gene_id", "exon_idx", "gene_name", "exon_id")
    )

    overlaps <- findOverlaps(guides_gr, all_exons_gr)

    if (length(overlaps) > 0) {
      raw_annotations <- data.frame(
        original_row_idx = guides_gr$original_row_idx[queryHits(overlaps)],
        gene_id = all_exons_gr$gene_id[subjectHits(overlaps)],
        gene_name = all_exons_gr$gene_name[subjectHits(overlaps)],
        exon_idx = all_exons_gr$exon_idx[subjectHits(overlaps)],
        exon_id = all_exons_gr$exon_id[subjectHits(overlaps)]
      )

      ambiguous_guides <- raw_annotations %>%
        group_by(original_row_idx) %>%
        summarise(n_genes = n_distinct(gene_id), .groups = 'drop') %>%
        filter(n_genes > 1) %>%
        pull(original_row_idx)

      if (length(ambiguous_guides) > 0) {
        message(paste("  - Excluding", length(ambiguous_guides), "guides that map to multiple genes."))
      }

      clean_annotations <- raw_annotations %>%
        filter(!original_row_idx %in% ambiguous_guides) %>%
        group_by(original_row_idx) %>%
        summarise(
          gene_id = first(gene_id), gene_name = first(gene_name),
          exon_idx = first(exon_idx), exon_id = first(exon_id),
          .groups = 'drop'
        )

      # Join annotations back to the temporary data frame
      annot_df <- left_join(annot_df, clean_annotations, by = "original_row_idx") %>%
        mutate(gene_type = if_else(!is.na(gene_id), "protein_coding", NA_character_))

      # Now, join these new columns back to the original full dataframe `df`
      df <- df %>%
        mutate(original_row_idx_main = row_number()) %>%
        left_join(
            select(annot_df, original_row_idx, gene_id, gene_name, exon_idx, exon_id, gene_type),
            by = c("original_row_idx_main" = "original_row_idx")
        ) %>%
        select(-original_row_idx_main) # Clean up index

    } else {
      message("⚠️ No guides overlap protein-coding exons.")
      df$gene_id <- NA_character_; df$gene_name <- NA_character_
      df$exon_idx <- NA_integer_;  df$exon_id <- NA_character_
      df$gene_type <- NA_character_
    }
} else {
    message("⚠️ No valid guides with coordinates to annotate.")
    df$gene_id <- NA_character_; df$gene_name <- NA_character_
    df$exon_idx <- NA_integer_;  df$exon_id <- NA_character_
    df$gene_type <- NA_character_
}

message("✅ Annotation complete.")


# --- INIT COLUMNS & FILTER EXONIC GUIDES ---
df <- df %>%
  mutate(
    gene_id_clean = sub("\\..*", "", gene_id),
    original_row_idx = row_number(),
    Binary_Domain = NA_integer_,
    PhyloP_AA_Up1 = NA_real_, PhyloP_AA_Cut = NA_real_,
    PhyloP_AA_Down1 = NA_real_, PhyloP_GuidePlusFlank = NA_real_
  )

exonic_df <- df %>% filter(!is.na(gene_id))

# --- AA CONTEXT (canonical transcript) ---
message("🐉🐉🐉 Projecting cut sites to protein space...")

get_gene_info   <- memoise(function(gene_id) genes(edb, filter = GeneIdFilter(gene_id)))
get_protein_seq <- memoise(function(tx_id)  proteins(edb, filter = TxIdFilter(tx_id))$protein_sequence[[1]])
get_protein_id  <- memoise(function(tx_id)  proteins(edb, filter = TxIdFilter(tx_id))$protein_id[1])

process_guide <- function(guide) {
  tryCatch({
    gene_info <- get_gene_info(guide$gene_id_clean)
    if (length(gene_info) == 0) return(NULL)

    # 1) Try canonical transcript
    tx_id <- gene_info$canonical_transcript[1]

    # 2) If canonical is missing, fall back to first transcript with protein
    if (is.na(tx_id)) {
      prot_info <- proteins(edb, filter = GeneIdFilter(guide$gene_id_clean))
      if (nrow(prot_info) == 0) return(NULL)
      tx_id <- prot_info$tx_id[1]
    }

    # Use cut_site from the input data
    cut_gr <- GRanges(
      seqnames = guide$chromosome,
      ranges   = IRanges(guide$cut_site, guide$cut_site)
    )
    seqlevelsStyle(cut_gr) <- "Ensembl"

    tx_coord_all <- genomeToTranscript(cut_gr, edb)
    tx_ranges <- tx_coord_all[[1]][mcols(tx_coord_all[[1]])$tx_id == tx_id]
    if (length(tx_ranges) == 0) return(NULL)

    prot_coord <- transcriptToProtein(tx_ranges, edb)
    if (length(prot_coord) == 0) return(NULL)

    aa_index <- start(prot_coord)
    prot_seq <- get_protein_seq(tx_id)
    if (is.na(prot_seq) || aa_index > nchar(prot_seq)) return(NULL)

    gene_strand <- as.character(strand(gene_info))[1]

    list(
      original_row_idx = guide$original_row_idx,
      protein_id = get_protein_id(tx_id),
      aa_index = aa_index,
      gene_strand = gene_strand,
      aa_cut  = substr(prot_seq, aa_index, aa_index),
      aa_up   = if (aa_index > 1) substr(prot_seq, max(1, aa_index - 16), aa_index - 1) else "",
      aa_down = if (aa_index < nchar(prot_seq)) substr(prot_seq, aa_index + 1, min(nchar(prot_seq), aa_index + 16)) else ""
    )
  }, error = function(e) NULL)
}


amino_acid_results <- if (nrow(exonic_df) > 0) {
  message(paste("🔹 Processing", nrow(exonic_df), "exonic guides..."))
  res <- purrr::map_dfr(split(exonic_df, seq_len(nrow(exonic_df))), process_guide)
  message(paste("✅ Amino-acid context found for", nrow(res), "guides."))
  res
} else {
  message("🐉🐉🐉 No exonic guides to process for AA features.")
  tibble(
    original_row_idx = integer(), protein_id = character(),
    aa_index = integer(), gene_strand = character(),
    aa_cut = character(), aa_up = character(), aa_down = character()
  )
}

# --- JOIN RESULTS ---
df_base <- df %>% select(-any_of(c("protein_id", "aa_index", "gene_strand", "aa_cut", "aa_up", "aa_down")))
df_final <- df_base %>% left_join(amino_acid_results, by = "original_row_idx")

# --- WRITE ---
message(paste0("🐉🐉🐉 DONE: Writing ", OUTPUT_CSV))
write_csv(df_final, OUTPUT_CSV, na = "NA")
