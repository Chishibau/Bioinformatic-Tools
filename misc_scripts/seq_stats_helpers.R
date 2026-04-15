#!/usr/bin/env Rscript

library(readr)
library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyverse)

theme_set(theme_classic(base_size = 8))

# fxns --------------------------------------------------------------------
BarQCPlot <-  function(df, y, ylabel, percent = F){
  
  hline_data <- df %>%
    group_by(`Sample type`) %>%
    summarise(med = median(!!sym(ylabel)))
  
  plot <- df %>%
    ggplot(aes(x = fct_inorder(SAMPLE), y = {{ y }}, fill = `Sample type`)) +
    geom_bar(stat = "identity") +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    ylab(ylabel) + 
    geom_hline(data = hline_data, aes(yintercept = med, col = `Sample type`), linetype='dotted') + 
    scale_fill_manual(values = c("cfDNA" = "black", "Tissue" = "#67c530", "gDNA" = "red")) + 
    scale_colour_manual(values = c("cfDNA" = "black", "Tissue" = "#67c530", "gDNA" = "red")) 
  
  if (percent){
    plot <- plot + 
      scale_y_continuous(
        limits = c(0, 100),
        expand = expansion(mult = c(0, 0.05))
      )  
  }
  
  return(plot)
}

FacetBarQCPlot <- function(seq_stat_table_full, sample_type = "", target = F, show_legend = T){
  
  if (sample_type != ""){
    seq_stat_table_full <- filter(seq_stat_table_full, `Sample type` == sample_type)
  }

  med_cov <- seq_stat_table_full %>%
    group_by(`Sample type`) %>% 
    arrange((MEDIAN_COVERAGE), .by_group = T) %>% 
    BarQCPlot(MEDIAN_COVERAGE, "MEDIAN_COVERAGE")
  
  total_cov <- seq_stat_table_full %>%
    group_by(`Sample type`) %>% 
    arrange((MEDIAN_COVERAGE), .by_group = T) %>% 
    BarQCPlot(`TOTAL_COVERAGE(1e8)`, "TOTAL_COVERAGE(1e8)", percent = F)
  
  frac_40 <- seq_stat_table_full %>%
    group_by(`Sample type`) %>% 
    arrange((MEDIAN_COVERAGE), .by_group = T) %>% 
    BarQCPlot(FRAC_GE_40, "FRAC_GE_40", percent = T)
  
  frac_90 <- seq_stat_table_full %>%
    group_by(`Sample type`) %>% 
    arrange((MEDIAN_COVERAGE), .by_group = T) %>% 
    BarQCPlot(FRAC_GE_90, "FRAC_GE_90", percent = T)
  
  frac_100 <- seq_stat_table_full %>%
    group_by(`Sample type`) %>% 
    arrange((MEDIAN_COVERAGE), .by_group = T) %>% 
    BarQCPlot(FRAC_GE_100, "FRAC_GE_100", percent = T)
  
  dup_bar <- seq_stat_table_full %>%
    group_by(`Sample type`) %>% 
    arrange((MEDIAN_COVERAGE), .by_group = T) %>% 
    BarQCPlot(DUP_PERCENT, "DUP_PERCENT", percent = T)
  
  at_dropout <- seq_stat_table_full %>%
    group_by(`Sample type`) %>% 
    arrange((MEDIAN_COVERAGE), .by_group = T) %>% 
    BarQCPlot(AT_DROPOUT, "AT_DROPOUT", percent = T)
  
  softclip <- seq_stat_table_full %>%
    group_by(`Sample type`) %>% 
    arrange((MEDIAN_COVERAGE), .by_group = T) %>% 
    BarQCPlot(PCT_SOFTCLIP, "PCT_SOFTCLIP", percent = T)
  
  chimeras <- seq_stat_table_full %>%
    group_by(`Sample type`) %>% 
    arrange((MEDIAN_COVERAGE), .by_group = T) %>% 
    BarQCPlot(PCT_CHIMERAS, "PCT_CHIMERAS", percent = T)
  
  fold80 <- seq_stat_table_full %>%
    group_by(`Sample type`) %>% 
    arrange((MEDIAN_COVERAGE), .by_group = T) %>% 
    BarQCPlot(FOLD_80_BASE_PENALTY, "FOLD_80_BASE_PENALTY")
  
    if (target) {
    med_tar <- seq_stat_table_full %>%
        group_by(`Sample type`) %>% 
        arrange((MEDIAN_COVERAGE), .by_group = T) %>% 
        BarQCPlot(MEDIAN_TARGET_DEPTH, "MEDIAN_TARGET_DEPTH")

    summary_qc_bar <- (
      med_cov + 
      med_tar + 
      total_cov + 
      dup_bar +
      at_dropout +
      softclip +
      chimeras +
      fold80 +
      frac_40 + 
      frac_90 + 
      frac_100
      ) + 
        plot_layout(ncol = 1, guides = "collect") & 
        theme(legend.position = 'bottom',)
  } else {
    summary_qc_bar <- (
      med_cov + 
      total_cov + 
      dup_bar +
      at_dropout +
      softclip +
      chimeras +
      fold80 +
      frac_40 + 
      frac_90 + 
      frac_100
      ) + 
        plot_layout(ncol = 1, guides = "collect") & 
        theme(legend.position = 'bottom',)
  }
}

median_target_from_regions <- function(reg_file) {
  if (!file.exists(reg_file)) return(NA_real_)
  reg <- suppressMessages(read_tsv(reg_file, col_names = FALSE, show_col_types = FALSE))
  if (ncol(reg) < 4) return(NA_real_)              # last column should be the mean depth
  md <- suppressWarnings(as.numeric(reg[[ncol(reg)]]))
  stats::median(md, na.rm = TRUE)
}

read_picard_metrics <- function(f, ext_regex, lines = 2) {
  txt <- readr::read_lines(f)
  
  s <- which(stringr::str_detect(txt, "^## METRICS CLASS"))[1]
  if (is.na(s)) {
    warning("No picard block found in: ", f)
    return(tibble::tibble())
  }
  
  block <- paste(txt[(s + 1):(s + lines)], collapse = "\n")
  
  df <- readr::read_tsv(
    I(block),
    show_col_types = FALSE,
    progress = FALSE
  )
  
  df$SAMPLE <- stringr::str_remove(basename(f), ext_regex)
  df <- dplyr::relocate(df, SAMPLE)
  
  if ("FOLD_90_BASE_PENALTY" %in% names(df)) {
    df$FOLD_90_BASE_PENALTY <- readr::parse_number(as.character(df$FOLD_90_BASE_PENALTY))
  }

  if ("FOLD_80_BASE_PENALTY" %in% names(df)) {
    df$FOLD_80_BASE_PENALTY <- readr::parse_number(as.character(df$FOLD_80_BASE_PENALTY))
  }
  
  df
}