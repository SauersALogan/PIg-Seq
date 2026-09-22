# =============================================================================
# Script:       04_flow_plots.R
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Flow cytometry plots showing IgA-positive and IgA-negative
#               fraction separation from FlowJo workspace
# Dependencies: plotting.R
# Input:        FlowJo/MRH_Flow_analysis.wsp (https://doi.org/10.48420/33951214)
# Output:       A .svg of flow cytometry separation. Figure 1 panel B
# =============================================================================

## ---- Paths (update for your system) ----
fig_dir  <- "/path/to/figure_output_dir"
flow_dir <- "/path/to/FlowJo"

library(CytoML)
library(ggcyto)
library(flowWorkspace)
library(ggrastr)
library(dplyr)
library(svglite)

source("statistics/functions/plotting.R")

## ---- Load flow data ----
ws <- open_flowjo_xml(file.path(flow_dir, "MRH_Flow_analysis.wsp"))
gs <- flowjo_to_gatingset(ws, name = "All Samples")

## ---- Separation plot ----
gs_sub_pos_neg <- gs[c("09-22_positive.fcs_0000000000010484", 
    "09-22_negative.fcs_0000000000210355")]

pData(gs_sub_pos_neg)$name <- c("Positive", "Negative")

pos_neg_gate_stats <- gs_pop_get_stats(gs_sub_pos_neg, 
    "/SYTO Positive/IgA Bound", type = "percent") %>%
    as.data.frame() %>%
    mutate(label = paste0(round(percent * 100, 1), "%"),
        name = pData(gs_sub_pos_neg)$name)

break_raw <- c(1e3, 1e4, 1e5, 1e6, 1e7)

Pos_Neg <- ggcyto(gs_sub_pos_neg, aes(x = `FL7-H`, y = `FSC-A`), 
        subset = "/SYTO Positive") +
    ggrastr::rasterise(geom_hex(bins = 256), dpi = 900) +
    geom_rect(aes(xmin = gate_rect@min["FL7-H"], xmax = gate_rect@max["FL7-H"],
        ymin = gate_rect@min["FSC-A"], ymax = 2e6),
        inherit.aes = FALSE, fill = NA, colour = "black", linewidth = lw_sm) +
    facet_wrap(~name) +
    ggtitle(NULL) +
    coord_cartesian(ylim = c(0, 2e6)) +
    scale_y_continuous(breaks = c(5e5, 1e6, 1.5e6, 2e6),
        labels = c("500K", "1.0M", "1.5M", "2.0M")) +
    labs(x = "Anti-IgA:PE", y = "Forward Scatter") +
    scale_x_continuous(
        breaks = break_pos_all,
        labels = c("0", "", "", 
            expression(10^3, 10^4, 10^5, 10^6, 10^7))) +
    scale_fill_gradientn(
        colours = c("#0000FF", "#00BFFF", "#44AA99", "#7ad151", "#fde725", "#FF0000"),
        values = c(0, 0.1, 0.3, 0.5, 0.75, 1.0),
        trans = "log10",
        limits = c(1, NA),
        na.value = "transparent",
        name = "Count") +
    theme_for_figures(show_labels = TRUE) +
    theme(strip.background = element_blank(),
        strip.text = element_text(size = base_sz),
        plot.margin = ggplot2::margin(0, 0, 0, 3, "mm"),
        panel.border = element_rect(colour = "black", fill = NA, 
            linewidth = lw_axis),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
    geom_text(data = pos_neg_gate_stats, 
        aes(x = Inf, y = Inf, label = label),
        hjust = 1.2, vjust = 2,
        size = base_sz / .pt,
        inherit.aes = FALSE)
Pos_Neg

# =====================================================================
# Save
# =====================================================================
panel_w <- 84.33
panel_h <- 37

ggsave(file.path(fig_dir, "FlowJo_Plots.svg"), Pos_Neg,
    width = panel_w, height = panel_h, units = "mm",
    device = svglite::svglite, dpi = 600)
