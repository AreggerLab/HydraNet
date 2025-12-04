#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: 03_get_phylop_and_impute.R <input_csv> <phyloP_bigwig> <output_csv>")

INPUT_CSV      <- args[1]
PHYLOP_BW_FILE <- args[2]
OUTPUT_CSV     <- args[3]
GUIDE_LENGTH   <- 23

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(readr)
  library(purrr)
})

#SEQINFO_HG38 <- readRDS("/home/nahar2/Cas12_Human_GuideOptimization/hg38_seqinfo.rds")

message(paste0("🐉🐉🐉 Reading: ", basename(INPUT_CSV)))
df <- read_csv(INPUT_CSV, show_col_types = FALSE, guess_max = 50000)

# =================================================================================
# FIX: The input file ALREADY has a 'chromosome' column from Script 1.
# We REMOVE the redundant rename(chromosome = seqnames) call.
# We KEEP the rename(cut_pos = cut_site) call, which is still needed.
# =================================================================================
message("🐉🐉🐉 Standardizing column names...")
df <- df %>%
  # The 'chromosome' column already exists. No action needed for it.
  # Rename 'cut_site' to 'cut_pos' for codon-level calculations.
  mutate(cut_pos = cut_site) %>%
  # Ensure chromosome names have 'chr' prefix for UCSC compatibility (bigWig).
  mutate(
    chromosome = if_else(!is.na(chromosome) & !startsWith(chromosome, "chr"), 
                         paste0("chr", chromosome), 
                         chromosome)
  )

df$guide_idx <- seq_len(nrow(df))


exonic_guides <- df %>% filter(!is.na(gene_strand))
if (nrow(exonic_guides) > 0) {
  message("🐉🐉🐉 Computing PhyloP Guide+Flank (±4 bp around protospacer) ...")

  flank_start <- pmax(exonic_guides$protospacer_start - 4L, 1L)
  flank_end   <- exonic_guides$protospacer_end + 4L

  gr_flank <- GRanges(
    seqnames = exonic_guides$chromosome,
    ranges   = IRanges(start = flank_start, end = flank_end),
    guide_idx = exonic_guides$guide_idx
  )
  seqlevelsStyle(gr_flank) <- "UCSC"
  #seqinfo(gr_flank) <- SEQINFO_HG38[seqlevels(gr_flank)]

  bw_flank  <- suppressWarnings(import(PHYLOP_BW_FILE, which = gr_flank))
  overlaps  <- findOverlaps(gr_flank, bw_flank)

  if (length(overlaps) > 0) {
    agg <- aggregate(score(bw_flank[subjectHits(overlaps)]),
                     by = list(guide_idx = gr_flank$guide_idx[queryHits(overlaps)]),
                     FUN = mean, na.rm = TRUE)
    df$PhyloP_GuidePlusFlank[agg$guide_idx] <- agg$x
  }
}

# CODON-LEVEL (EXONIC)
get_bulk_phyloP_scores <- function(guides_df, pos_vector, offset_vector = 0) {
  start_pos <- pos_vector + offset_vector
  valid <- !is.na(start_pos) & start_pos > 0 & !is.na(guides_df$chromosome)
  if (sum(valid) == 0) return(data.frame(guide_idx = integer(), score = numeric()))
  
  gr <- GRanges(
    seqnames  = guides_df$chromosome[valid],
    ranges    = IRanges(start = start_pos[valid], width = 3),
    guide_idx = guides_df$guide_idx[valid]
  )
  seqlevelsStyle(gr) <- "UCSC"
  #seqinfo(gr) <- SEQINFO_HG38[seqlevels(gr)]
  scores_list <- suppressWarnings(import(PHYLOP_BW_FILE, which = gr, as = "NumericList"))
  mean_scores <- vapply(scores_list, mean, numeric(1), na.rm = TRUE)
  data.frame(guide_idx = gr$guide_idx, score = mean_scores)
}

if (nrow(exonic_guides) > 0) {
  message("🐉🐉🐉 Computing codon-context PhyloP (Cut/Up/Down) ...")
  get_codon_start <- function(cut) cut - ((cut - 1) %% 3)
  
  exonic_guides_with_pos <- exonic_guides %>% filter(!is.na(cut_pos))
  
  if(nrow(exonic_guides_with_pos) > 0) {
    exonic_guides_with_pos$codon_start <- get_codon_start(exonic_guides_with_pos$cut_pos)
    offsets <- ifelse(exonic_guides_with_pos$gene_strand == "+", 3, -3)

    scores_cut  <- get_bulk_phyloP_scores(exonic_guides_with_pos, exonic_guides_with_pos$codon_start)
    scores_up   <- get_bulk_phyloP_scores(exonic_guides_with_pos, exonic_guides_with_pos$codon_start, -offsets)
    scores_down <- get_bulk_phyloP_scores(exonic_guides_with_pos, exonic_guides_with_pos$codon_start,  offsets)

    if (nrow(scores_cut)  > 0) df$PhyloP_AA_Cut[ scores_cut$guide_idx]  <- scores_cut$score
    if (nrow(scores_up)   > 0) df$PhyloP_AA_Up1[ scores_up$guide_idx]   <- scores_up$score
    if (nrow(scores_down) > 0) df$PhyloP_AA_Down1[scores_down$guide_idx] <- scores_down$score
  } else {
    message("  - ⚠️ No exonic guides with valid cut positions to process.")
  }
}

# NUCLEOTIDE-LEVEL (EXONIC)
message("🐉🐉🐉 Computing nucleotide-level PhyloP (guide orientation) ...")
nuc_col_names <- paste0("PhyloP_Nuc_", 1:GUIDE_LENGTH)
df[, nuc_col_names] <- NA_real_

if (nrow(exonic_guides) > 0) {
  gr_nuc <- GRanges(
    seqnames  = exonic_guides$chromosome,
    ranges    = IRanges(start = pmax(exonic_guides$protospacer_start, 1L), width = GUIDE_LENGTH),
    guide_idx = exonic_guides$guide_idx,
    strand    = exonic_guides$gene_strand
  )
  seqlevelsStyle(gr_nuc) <- "UCSC"
  #seqinfo(gr_nuc) <- SEQINFO_HG38[seqlevels(gr_nuc)]
  scores_list <- suppressWarnings(import(PHYLOP_BW_FILE, which = gr_nuc, as = "NumericList"))
  
  fix_len <- function(vec, len) {
    v <- as.numeric(vec); n <- length(v)
    if (n == len) return(v)
    if (n >  len) return(v[1:len])
    c(v, rep(NA_real_, len - n))
  }

  if (length(scores_list) > 0) {
    mat <- map2(scores_list, as.character(strand(gr_nuc)), function(v, s) {
      vv <- fix_len(v, GUIDE_LENGTH)
      if (identical(s, "-")) rev(vv) else vv
    }) %>% do.call(rbind, .)
    df[gr_nuc$guide_idx, nuc_col_names] <- mat
  }
}

# FINAL
df <- df %>%
  select(-any_of(c("guide_idx")))


write_csv(df, OUTPUT_CSV, na = "NA")
message(paste0("🐉🐉🐉 DONE: Wrote output to ", OUTPUT_CSV))
