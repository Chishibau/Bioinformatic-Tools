#!/usr/bin/env Rscript

# Get the directory where this script is located
script_path <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", script_path[grep("--file=", script_path)])
script_dir <- dirname(script_path)

# Source the helper file from the same directory
helper_file <- file.path(script_dir, "seq_stats_helpers.R")

if (!file.exists(helper_file)) {
  stop("Error: seq_stats_helpers.R not found in script directory: ", script_dir)
}

source(helper_file)
cat("Loaded helpers from:", helper_file, "\n\n")

## load in variables
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript parse_seq_stats.R <working_dir> [target_bool]\n")
  cat("\nExpected directory structure:\n")
  cat("  working_dir/\n")
  cat("    mosdepth/           (optional: *.mosdepth.global.dist.txt, *.mosdepth.summary.txt, *.regions.bed.gz)\n")
  cat("    samblaster/         (optional: *.err files with duplication rates)\n")
  cat("    alignment_metrics/  (optional: *.alignment_metrics.txt)\n")
  cat("    gc_bias/            (optional: *.gc_summary_metrics.txt)\n")
  cat("    wgs_metrics/        (optional: *.wgs_metrics.txt)\n")
  cat("\nAlternatively, all files can be in the working_dir root.\n")
  quit(status = 1)
}

working_dir <- args[1]
target_bool <- if (length(args) >= 2) args[2] else "false"

if (!dir.exists(working_dir)) {
  stop("Error: Working directory does not exist: ", working_dir)
}

cat("Working directory:", working_dir, "\n")
cat("Target mode:", target_bool, "\n")

## Define possible subdirectories
subdirs <- c(
  "mosdepth",
  "samblaster", 
  "alignment_metrics",
  "gc_bias",
  "wgs_metrics"
)

## Recursively find all files in working_dir and subdirectories
all_files <- list.files(
  working_dir,
  recursive = TRUE,
  full.names = TRUE
)

cat("Found", length(all_files), "files total\n")

## Alternatively, if files are in root or specific subdirs
files <- c()
for (subdir in subdirs) {
  subdir_path <- file.path(working_dir, subdir)
  if (dir.exists(subdir_path)) {
    subdir_files <- list.files(subdir_path, full.names = TRUE, recursive = TRUE)
    files <- c(files, subdir_files)
    cat("Found", length(subdir_files), "files in", subdir, "\n")
  }
}

# Also check root directory
root_files <- list.files(working_dir, full.names = TRUE, recursive = FALSE)
root_files <- root_files[!dir.exists(root_files)]  # exclude directories
files <- unique(c(files, root_files))

# If no subdirectories found, use all files from recursive search
if (length(files) == 0) {
  files <- all_files
  cat("No subdirectories found, using all files from recursive search\n")
}

cat("Total files to process:", length(files), "\n")

## depth information from mosdepth
global_files <- files[grepl("\\.mosdepth\\.global\\.dist\\.txt$", basename(files))]
cat("Mosdepth global files:", length(global_files), "\n")

if (length(global_files) > 0) {
  medians <- map_dfr(global_files, function(f) {
    dat <- read_tsv(
      f,
      col_names = c("region","depth","frac"),
      col_types = "cid",
      show_col_types = FALSE
    )
    
    total <- dat %>% filter(region == "total")
    
    m <- total %>%
      filter(frac >= 0.5) %>%
      summarise(MEDIAN_COVERAGE = max(depth, na.rm = TRUE)) %>%
      pull(MEDIAN_COVERAGE)
    
    thresholds <- c(40L, 90L, 100L, 120L)
    thr_tbl <- tibble(depth = thresholds) %>%
      left_join(total %>% select(depth, frac), by = "depth") %>%
      mutate(frac = coalesce(frac, 0))
    
    frac_cols <- setNames(as.list(thr_tbl$frac), paste0("FRAC_GE_", thr_tbl$depth))
    
    sum_file <- sub("\\.mosdepth\\.global\\.dist\\.txt$", ".mosdepth.summary.txt", f)

    tot_cov_100m <- NA_real_
    mean_cov_total <- NA_real_

    if (file.exists(sum_file)) {
    sdat <- read_tsv(sum_file, show_col_types = FALSE, col_types = cols())

    chrom_col <- intersect(c("chrom", "region", "chromosome"), names(sdat))
    chrom_col <- if (length(chrom_col)) chrom_col[1] else names(sdat)[1]

    bases_col <- if ("bases" %in% names(sdat)) "bases" else names(sdat)[3]
    mean_col  <- if ("mean"  %in% names(sdat)) "mean"  else names(sdat)[4]

    total_row <- sdat %>% filter(.data[[chrom_col]] == "total")

    ## total coverage in 1e8 bases (your existing metric)
    if (nrow(total_row) > 0) {
        tot_cov_100m <- as.numeric(total_row[[bases_col]][1]) / 1e8
        mean_cov_total <- as.numeric(total_row[[mean_col]][1])
    } else {
        tot_cov_100m <- sum(sdat[[bases_col]], na.rm = TRUE) / 1e8
    }
    }

    
    reg_file <- sub("\\.mosdepth\\.global\\.dist\\.txt$", ".regions.bed.gz", f)

    if (file.exists(reg_file)) {
      med_target <- median_target_from_regions(reg_file)

      tibble(
        SAMPLE = sub("\\.mosdepth\\.global\\.dist\\.txt$", "", basename(f)),
        MEDIAN_COVERAGE = ifelse(is.finite(m), as.integer(m), NA_integer_),
        MEDIAN_TARGET_DEPTH = med_target,
        MEAN_COVERAGE = mean_cov_total
        ) %>%
        bind_cols(as_tibble(frac_cols)) %>%
        mutate(`TOTAL_COVERAGE(1e8)` = tot_cov_100m)

    } else {
      tibble(
        SAMPLE = sub("\\.mosdepth\\.global\\.dist\\.txt$", "", basename(f)),
        MEDIAN_COVERAGE = ifelse(is.finite(m), as.integer(m), NA_integer_),
        MEAN_COVERAGE = mean_cov_total
        ) %>%
        bind_cols(as_tibble(frac_cols)) %>%
        mutate(`TOTAL_COVERAGE(1e8)` = tot_cov_100m)

    }
  })
} else {
  cat("Warning: No mosdepth files found\n")
  medians <- tibble(SAMPLE = character())
}

## duplicate rate from samblaster .err files
err_files <- files[grepl("\\.err$", basename(files))]
cat("Samblaster .err files:", length(err_files), "\n")

# Helper to extract the first "Duplication rate: 16.1217%" hit in a file
extract_dup_rate <- function(f) {
  lines <- readLines(f, warn = FALSE)
  # case-insensitive; tolerates extra spaces
  m <- str_match(lines, regex("Duplication\\s*rate\\s*:\\s*([0-9]+(?:\\.[0-9]+)?)\\s*%", ignore_case = TRUE))
  val <- suppressWarnings(as.numeric(na.omit(m[, 2])))[1]  # first match or NA
  tibble(
    SAMPLE      = tools::file_path_sans_ext(basename(f)),
    DUP_PERCENT = val,           # e.g., 16.1217
    dup_rate_fraction = val / 100     # e.g., 0.161217 (optional)
  )
}

if (length(err_files) > 0) {
  dup_rates <- map_dfr(err_files, extract_dup_rate)
} else {
  cat("Warning: No samblaster .err files found\n")
  dup_rates <- tibble(SAMPLE = character(), DUP_PERCENT = numeric())
}

## alignment metrics
ext <- "\\.alignment_metrics(\\.txt)?$"
align_files <- files[grepl(ext, basename(files))]
cat("Alignment metrics files:", length(align_files), "\n")

if (length(align_files) > 0) {
  alignment_stats <- purrr::map_dfr(align_files, read_picard_metrics, ext_regex = ext, lines = 4) %>%
    filter(CATEGORY == "PAIR") %>%
    select(SAMPLE, PCT_SOFTCLIP, PCT_CHIMERAS) %>%
    mutate(PCT_SOFTCLIP = PCT_SOFTCLIP * 100,
           PCT_CHIMERAS = PCT_CHIMERAS * 100)
} else {
  cat("Warning: No alignment metrics files found\n")
  alignment_stats <- tibble(SAMPLE = character(), PCT_SOFTCLIP = numeric(), PCT_CHIMERAS = numeric())
}

## GC bias metrics
ext <- "\\.gc_summary(\\.txt)?$"
gc_files <- files[grepl(ext, basename(files))]
cat("GC bias files:", length(gc_files), "\n")

if (length(gc_files) > 0) {
  gc_bias <- purrr::map_dfr(gc_files, read_picard_metrics, ext_regex = ext)
  gc_bias <- gc_bias %>%
    select(SAMPLE, GC_DROPOUT, AT_DROPOUT, ALIGNED_READS)
} else {
  cat("Warning: No GC bias files found\n")
  gc_bias <- tibble(SAMPLE = character(), GC_DROPOUT = numeric(), AT_DROPOUT = numeric(), ALIGNED_READS = numeric())
}

## WGS metrics
ext <- "\\.wgs_metrics(\\.txt)?$"
wgs_files <- files[grepl(ext, basename(files))]
cat("WGS metrics files:", length(wgs_files), "\n")

if (length(wgs_files) > 0) {
  wgs_metrics <- purrr::map_dfr(wgs_files, read_picard_metrics, ext_regex = ext) %>%
    select(SAMPLE, FOLD_80_BASE_PENALTY)
} else {
  cat("Warning: No WGS metrics files found\n")
  wgs_metrics <- tibble(SAMPLE = character(), FOLD_80_BASE_PENALTY = numeric())
}

## compile all data and save
cat("\nMerging all metrics...\n")

all_stats <- medians

if (nrow(medians) > 0) {
  all_stats <- all_stats %>% na.omit()
}

if (nrow(dup_rates) > 0) {
  all_stats <- all_stats %>% left_join(dup_rates, by = "SAMPLE")
}

if (nrow(alignment_stats) > 0) {
  all_stats <- all_stats %>% left_join(alignment_stats, by = "SAMPLE")
}

if (nrow(gc_bias) > 0) {
  all_stats <- all_stats %>% left_join(gc_bias, by = "SAMPLE")
}

if (nrow(wgs_metrics) > 0) {
  all_stats <- all_stats %>% left_join(wgs_metrics, by = "SAMPLE")
}

if (nrow(all_stats) == 0) {
  stop("Error: No data found. Check that files exist in the working directory.")
}

all_stats <- all_stats %>%
  mutate(`Sample type` = case_when(
    grepl("WBC|gDNA", SAMPLE, ignore.case = TRUE) ~ "gDNA",
    grepl("utDNA|ctDNA|cfDNA", SAMPLE, ignore.case = TRUE) ~ "cfDNA",
    grepl("Tissue|FiT|FFPE", SAMPLE, ignore.case = TRUE) ~ "Tissue",
    TRUE ~ NA_character_
  )
  )

all_stats <- all_stats %>%
  mutate(across(starts_with("FRAC_GE_"),
                ~ . * 100,
                .names = "{.col}"))

cat("\nProcessed", nrow(all_stats), "samples\n")

## Save outputs in the working directory
output_dir <- working_dir
# summary_qc_bar <- FacetBarQCPlot(seq_stat_table_full = all_stats, target = target_bool)

# plot_file <- file.path(output_dir, "seq_stats_plot.pdf")
table_file <- file.path(output_dir, "seq_stats_all.tsv")

# ggsave(filename = plot_file, plot = summary_qc_bar, width = 16, height = 10)
write_tsv(all_stats, table_file, na='')

cat("\nOutputs saved:\n")
# cat("  Plot:", plot_file, "\n")
cat("  Table:", table_file, "\n")