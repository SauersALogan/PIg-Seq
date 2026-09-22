# =============================================================================
# Script:       plotting.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Shared plotting theme, color palette, and size parameters
#               sourced by all analysis notebooks
# =============================================================================

library(ggplot2)

## ---- Size parameters ----
base_sz  <- 7
pt_to_mm <- 2.46
pt_small <- 0.9
pt_mean  <- 2.4
lw_main  <- 0.6
lw_sm    <- 0.3
lw_axis  <- 0.6
arr_len  <- 0.15

## ---- Color palette ----
col_native  <- "#737373"
col_neg     <- "#4393C3"
col_pos     <- "#DD5566"

## ---- Theme ----
theme_for_figures <- function(show_labels = FALSE) {
    base <- theme_classic(base_size = base_sz) +
        theme(
            axis.text    = element_text(colour = "black", size = base_sz),
            legend.text  = element_text(size = base_sz),
            legend.title = element_text(size = base_sz),
            panel.grid   = element_blank(),
            plot.title   = element_text(face = "bold.italic", hjust = 0.5, 
                size = base_sz),
            legend.position = "none",
            axis.line    = element_line(colour = "black", linewidth = lw_axis),
            axis.ticks   = element_line(colour = "black", linewidth = lw_axis),
            plot.margin  = ggplot2::margin(0, 0, 0, 0, "mm"),
            axis.title.y = element_text(margin = ggplot2::margin(
                t = 0, r = 3, b = 0, l = 0)),
            text = element_text(family = "Arial")
        )
    if (show_labels) {
        base
    } else {
        base + theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text    = element_blank()
        )
    }
}
