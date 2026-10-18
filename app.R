# Interactive Backward Linkages Network App (Industry x Industry MIP 2018)
# Enhanced R Shiny App: Interactive Sugiyama Layout (Draggable Nodes, Hover Tooltips,
# Gradual Layer Animation, and 3-Part NAICS Manufacturing Selector)

library(shiny)
library(bslib)
library(tidyverse)
library(igraph)
library(ggraph)
library(scales)
library(viridisLite)

# -------------------------------------------------------------------------
# 1. Data Loading & Pre-computation (runs once at app startup)
# -------------------------------------------------------------------------

ruta_mip <- "mip_csv/conjunto de datos"
file_mip_ixi <- "conjunto_de_datos_mip_t_pb_ixi_32018.csv"

mip_ixi_d_m <- read.csv(file.path(ruta_mip, file_mip_ixi), encoding = "UTF-8")

col_id_ixi <- colnames(mip_ixi_d_m)

col_idx_ixi <- tibble(
  index = seq_along(col_id_ixi),
  original_name = col_id_ixi
) %>%
  mutate(
    category = str_extract(original_name, "^[^.]+"),
    sector = str_match(original_name, "\\.(\\d{4})\\.\\.\\.")[, 2]
  )

z_idx_ixi <- col_idx_ixi %>%
  drop_na() %>%
  filter(category == "DI") %>%
  pull(index)

ramas <- col_idx_ixi %>%
  drop_na() %>%
  pull(sector)

# Intermediate demand matrix Z (263 x 263, Total = Domestic + Imported)
Z_ixi <- as.matrix.data.frame(mip_ixi_d_m[z_idx_ixi - 3, z_idx_ixi])
dimnames(Z_ixi) <- list(ramas, ramas)

# Total output X (UTPB = P.1 Production at basic prices)
X_ixi <- mip_ixi_d_m[z_idx_ixi - 3, which(col_idx_ixi$category == "UTPB")]

# Technical requirements matrix A = Z * diag(X)^(-1)
A_ixi <- sweep(Z_ixi, MARGIN = 2, STATS = X_ixi, FUN = "/")

# Leontief inverse L = (I - A)^(-1) and backward linkage multipliers
L_ixi <- solve(diag(nrow(A_ixi)) - A_ixi)
mu_ixi <- tibble(
  rama = ramas,
  mu = colSums(L_ixi)
)

# Extract full & abbreviated descriptions for all 263 sectors
sector_catalog <- tibble(
  raw_desc = as.character(mip_ixi_d_m[z_idx_ixi - 3, 1])
) %>%
  mutate(
    code = str_extract(raw_desc, "^\\d{4}"),
    subsector_2d = str_sub(code, 1, 2),
    full_name = str_replace(raw_desc, "^\\d{4}---", ""),
    short_name = full_name %>%
      str_replace(
        paste0(
          "^Fabricación de |^Elaboración de |^Industria básica del |",
          "^Industria básica de |^Industrias de |^Industria de las |",
          "^Industria del |^Confección de |^Preparación y envasado de |",
          "^Preparación e hilado de "
        ),
        ""
      ) %>%
      str_to_sentence() %>%
      str_trunc(50, ellipsis = "...")
  ) %>%
  left_join(mu_ixi, by = c("code" = "rama"))

# Manufacturing sectors divided into the 3 NAICS groups (31, 32, 33)
mfg_catalog <- sector_catalog %>%
  filter(subsector_2d %in% c("31", "32", "33")) %>%
  mutate(
    naics_group = case_when(
      subsector_2d == "31" ~ "31 - Basic Manufacturing",
      subsector_2d == "32" ~ "32 - Intermediate Processed Goods",
      subsector_2d == "33" ~ "33 - Advanced Manufacturing"
    ),
    label = paste0(code, " - ", short_name)
  )

# Grouped choices list for selectizeInput
get_grouped_mfg_choices <- function(group_filter = "all") {
  df <- mfg_catalog
  if (group_filter %in% c("31", "32", "33")) {
    df <- df %>% filter(subsector_2d == group_filter)
  }
  groups <- unique(df$naics_group)
  res <- lapply(groups, function(g) {
    sub_df <- df %>% filter(naics_group == g)
    setNames(sub_df$code, sub_df$label)
  })
  names(res) <- groups
  res
}

# Edge list representation of A_ixi (excluding self-loops row == col)
ixi_edgelist <- cbind("row" = colnames(A_ixi), as.data.frame(A_ixi)) %>%
  pivot_longer(
    cols = 2:last_col(),
    names_to = "col",
    values_to = "weight"
  ) %>%
  filter(row != col, weight > 0)

# -------------------------------------------------------------------------
# 2. Network Construction, Sugiyama Layout & Static Export Function
# -------------------------------------------------------------------------

build_edges_rama <- function(rama, threshold = 0.01, threshold2 = 0.05) {
  sector <- as.character(rama)

  # Direct inputs (1st layer: suppliers -> target sector)
  edges_direct <- ixi_edgelist %>%
    filter(col == sector, weight >= threshold) %>%
    mutate(round = "Direct (A)", stage_num = 1L)

  suppliers_direct <- edges_direct$row

  # First round of indirect inputs (2nd layer: tier-2 suppliers -> direct suppliers)
  edges_r1 <- ixi_edgelist %>%
    filter(col %in% suppliers_direct, weight >= threshold2) %>%
    mutate(round = "Indirect (A²)", stage_num = 2L)

  # Combine and deduplicate
  edges_rama <- bind_rows(edges_direct, edges_r1) %>%
    distinct(row, col, .keep_all = TRUE)

  list(
    edges = edges_rama,
    n_direct = nrow(edges_direct),
    n_indirect = nrow(edges_r1),
    suppliers_direct = unique(suppliers_direct)
  )
}

# Build interactive JSON payload using exact ggraph Sugiyama coordinates
build_interactive_network <- function(rama, threshold = 0.01, threshold2 = 0.05) {
  sector <- as.character(rama)
  net_data <- build_edges_rama(sector, threshold, threshold2)
  edges_rama <- net_data$edges

  sector_info <- sector_catalog %>% filter(code == sector)
  target_short <- if (nrow(sector_info) > 0) sector_info$short_name[1] else sector

  if (nrow(edges_rama) == 0) {
    return(list(
      empty = TRUE,
      sector = sector,
      target_short = target_short,
      threshold = threshold,
      threshold2 = threshold2
    ))
  }

  # Create igraph and extract exact ggraph Sugiyama layout
  red <- graph_from_data_frame(d = edges_rama, directed = TRUE)
  glay <- create_layout(red, layout = "sugiyama")

  # Map Plasma color scale (direction = -1, matching scale_edge_color_viridis)
  w_range <- range(edges_rama$weight)
  if (diff(w_range) < 1e-9) {
    w_domain <- c(max(0, w_range[1] - 0.01), w_range[2] + 0.01)
  } else {
    w_domain <- w_range
  }
  plasma_pal <- col_numeric(
    palette = viridisLite::plasma(256, direction = -1),
    domain = w_domain
  )

  # Prepare nodes table with coordinates matching coord_flip():
  # In ggraph + coord_flip(), horizontal axis is glay$y (1 = left target, max = right suppliers)
  # and vertical axis is glay$x (larger glay$x = higher up on plot)
  nodes_df <- as_tibble(glay) %>%
    select(code = name, orig_v = x, orig_h = y) %>%
    left_join(sector_catalog, by = "code") %>%
    mutate(
      short_name = coalesce(short_name, code),
      full_name = coalesce(full_name, code),
      mu = round(coalesce(mu, 1), 3),
      role = case_when(
        code == sector ~ "Target Sector",
        code %in% net_data$suppliers_direct ~ "Direct Supplier (A)",
        TRUE ~ "Indirect Supplier (A²)"
      ),
      anim_tier = case_when(
        code == sector ~ 0L,
        code %in% net_data$suppliers_direct ~ 1L,
        TRUE ~ 2L
      )
    )

  # Enrich edges with node names, colors, and animation tier
  edges_df <- edges_rama %>%
    left_join(sector_catalog %>% select(code, source_short = short_name), by = c("row" = "code")) %>%
    left_join(sector_catalog %>% select(code, target_short = short_name), by = c("col" = "code")) %>%
    mutate(
      id = paste0(row, "->", col),
      source = row,
      target = col,
      color = plasma_pal(weight),
      weight_fmt = sprintf("%.4f", weight),
      pct_fmt = sprintf("%.2f%%", weight * 100)
    )

  # Generate legend ticks for the bottom colorbar
  legend_ticks <- pretty(w_range, n = 4)
  legend_ticks <- legend_ticks[legend_ticks >= w_range[1] & legend_ticks <= w_range[2]]
  if (length(legend_ticks) < 2) {
    legend_ticks <- round(seq(w_range[1], w_range[2], length.out = 3), 3)
  }

  list(
    empty = FALSE,
    sector = sector,
    target_short = target_short,
    threshold = threshold,
    threshold2 = threshold2,
    min_weight = w_range[1],
    max_weight = w_range[2],
    legend_ticks = legend_ticks,
    legend_gradient = viridisLite::plasma(32, direction = -1),
    nodes = purrr::transpose(nodes_df %>% select(code, orig_h, orig_v, short_name, full_name, mu, role, anim_tier)),
    edges = purrr::transpose(edges_df %>% select(id, source, target, weight, weight_fmt, pct_fmt, round, stage_num, color, source_short, target_short))
  )
}

# Static ggraph function (for high-res PNG download)
graficar_red_rama <- function(rama, threshold = 0.01, threshold2 = 0.05, panel_bg = "grey80") {
  sector <- as.character(rama)
  net_data <- build_edges_rama(sector, threshold, threshold2)
  edges_rama <- net_data$edges

  if (nrow(edges_rama) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No linkages above thresholds", size = 5) +
        theme_void()
    )
  }

  red <- graph_from_data_frame(d = edges_rama, directed = TRUE)

  ggraph(red, layout = "sugiyama") +
    geom_edge_link(
      aes(color = weight),
      width = 1,
      alpha = 0.75,
      arrow = arrow(length = unit(1.5, "mm"), type = "open"),
      end_cap = circle(6, "mm"),
      start_cap = circle(6, "mm")
    ) +
    geom_node_text(aes(label = name), size = 3.3, fontface = "bold") +
    scale_edge_color_viridis(
      option = "plasma",
      guide = guide_edge_colorbar(
        title = "Coefficient",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(6, "cm"),
        barheight = unit(0.4, "cm")
      ),
      direction = -1
    ) +
    theme(
      panel.background = element_rect(fill = panel_bg, color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 3)),
      plot.subtitle = element_text(face = "italic", size = 11, color = "grey25", margin = margin(b = 8)),
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 14, 10, 10)
    ) +
    labs(
      title = paste0("Backward Linkages: Rama ", sector),
      subtitle = paste0(
        "Threshold: ", threshold,
        "   -   Threshold 2: ", threshold2,
        "   -   Depth: 2"
      )
    ) +
    coord_flip()
}

# -------------------------------------------------------------------------
# 3. User Interface (Light-Themed Minimal GUI + D3 Interactive Canvas)
# -------------------------------------------------------------------------

ui <- page_sidebar(
  title = div(
    style = "display: flex; align-items: baseline; gap: 10px;",
    span("Backward Linkages Explorer", style = "font-weight: 700; font-size: 1.15rem;"),
    span("NAICS 4-Digit Manufacturing (MIP 2018)", style = "font-size: 0.82rem; color: #6c757d;")
  ),
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2b3e50",
    bg = "#ffffff",
    fg = "#212529",
    base_font = font_google("Inter")
  ),
  tags$head(
    tags$script(src = "https://d3js.org/d3.v7.min.js"),
    tags$style(HTML("
      /* NAICS 3-group selector styling */
      .naics-group-pills .btn-group {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 4px;
        width: 100%;
      }
      .naics-group-pills .btn {
        font-size: 0.74rem;
        padding: 4px 6px;
        border-radius: 5px !important;
        text-align: center;
        font-weight: 600;
      }
      /* Interactive Network Wrapper (Google Material 3 / Cloud Visual Style) */
      #net-container {
        position: relative;
        width: 100%;
        height: 660px;
        background: #ffffff;
        display: flex;
        flex-direction: column;
        user-select: none;
        overflow: hidden;
      }
      #net-header {
        padding: 10px 18px 8px 18px;
        display: flex;
        justify-content: space-between;
        align-items: center;
        background: #ffffff;
      }
      #net-title-wrap {
        display: flex;
        flex-direction: column;
        gap: 3px;
      }
      #net-title {
        font-weight: 700;
        font-size: 1.16rem;
        color: #0f172a;
        letter-spacing: -0.015em;
        margin: 0;
        line-height: 1.2;
      }
      #net-subtitle {
        font-size: 0.82rem;
        color: #64748b;
        margin: 0;
        display: flex;
        align-items: center;
        gap: 8px;
      }
      .sub-pill {
        display: inline-flex;
        align-items: center;
        background: #f1f5f9;
        color: #334155;
        font-weight: 600;
        font-size: 0.74rem;
        padding: 2px 8px;
        border-radius: 999px;
        border: 1px solid #e2e8f0;
        font-style: normal;
      }
      #canvas-reset-btn {
        display: inline-flex;
        align-items: center;
        gap: 6px;
        font-size: 0.76rem;
        font-weight: 600;
        color: #1a73e8;
        background: #f0f6ff;
        border: 1px solid #bfdbfe;
        border-radius: 999px;
        padding: 5px 13px;
        transition: all 0.18s ease;
        box-shadow: 0 1px 2px rgba(26, 115, 232, 0.08);
      }
      #canvas-reset-btn:hover {
        background: #e0edff;
        border-color: #93c5fd;
        transform: translateY(-1px);
        box-shadow: 0 3px 8px rgba(26, 115, 232, 0.15);
      }
      #net-canvas-wrap {
        position: relative;
        flex: 1;
        margin: 0 16px;
        background: radial-gradient(circle at 20% 50%, #ffffff 0%, #f8fafc 55%, #f1f5f9 100%);
        border: 1px solid #e2e8f0;
        border-radius: 12px;
        box-shadow: inset 0 2px 6px rgba(15, 23, 42, 0.03), 0 1px 3px rgba(15, 23, 42, 0.04);
        overflow: hidden;
      }
      #net-svg {
        width: 100%;
        height: 100%;
        display: block;
      }
      /* Top-right floating role legend inside canvas */
      .canvas-role-legend {
        position: absolute;
        top: 10px;
        right: 12px;
        display: flex;
        gap: 10px;
        background: rgba(255, 255, 255, 0.88);
        backdrop-filter: blur(6px);
        border: 1px solid #e2e8f0;
        border-radius: 999px;
        padding: 4px 12px;
        font-size: 0.72rem;
        font-weight: 600;
        color: #334155;
        box-shadow: 0 2px 6px rgba(15, 23, 42, 0.05);
        pointer-events: none;
      }
      .role-item {
        display: inline-flex;
        align-items: center;
        gap: 5px;
      }
      .role-dot {
        width: 8px;
        height: 8px;
        border-radius: 50%;
        display: inline-block;
      }
      #net-footer {
        height: 60px;
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: center;
        background: #ffffff;
        padding-top: 4px;
      }
      .colorbar-title {
        font-size: 0.76rem;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.06em;
        color: #475569;
        margin-bottom: 4px;
      }
      .colorbar-bar {
        width: 260px;
        height: 10px;
        border-radius: 999px;
        border: 1px solid rgba(15, 23, 42, 0.12);
        box-shadow: inset 0 1px 2px rgba(0, 0, 0, 0.1);
      }
      .colorbar-ticks {
        width: 260px;
        display: flex;
        justify-content: space-between;
        font-size: 0.72rem;
        font-weight: 600;
        font-variant-numeric: tabular-nums;
        color: #475569;
        margin-top: 3px;
      }
      /* Google Material 3 Floating Tooltip Card */
      #net-tooltip {
        position: absolute;
        pointer-events: none;
        background: rgba(255, 255, 255, 0.98);
        backdrop-filter: blur(8px);
        border: 1px solid #e2e8f0;
        box-shadow: 0 10px 25px -5px rgba(15, 23, 42, 0.14), 0 4px 10px -2px rgba(15, 23, 42, 0.08);
        border-radius: 10px;
        padding: 11px 14px;
        font-size: 0.8rem;
        color: #1e293b;
        max-width: 300px;
        z-index: 100;
        opacity: 0;
        transform: translateY(3px);
        transition: opacity 0.14s ease, transform 0.14s ease;
      }
      #net-tooltip.visible {
        opacity: 1;
        transform: translateY(0);
      }
      /* Node & Edge Interactive Styles */
      .node-group {
        cursor: grab;
      }
      .node-group:active {
        cursor: grabbing;
      }
      .node-inner {
        transition: transform 0.18s cubic-bezier(0.34, 1.56, 0.64, 1);
      }
      .node-group:hover .node-inner,
      .node-group.selected .node-inner {
        transform: scale(1.09);
      }
      .node-capsule {
        transition: stroke 0.18s ease, stroke-width 0.18s ease, filter 0.18s ease;
      }
      .node-label {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
        font-size: 12px;
        font-weight: 700;
        letter-spacing: -0.01em;
        text-anchor: middle;
        dominant-baseline: central;
        pointer-events: none;
      }
      .target-pulse {
        fill: none;
        stroke: #1a73e8;
        stroke-width: 1.8px;
        opacity: 0.35;
        animation: googlePulse 2.6s infinite cubic-bezier(0.4, 0, 0.2, 1);
      }
      @keyframes googlePulse {
        0% { transform: scale(0.96); opacity: 0.45; }
        70% { transform: scale(1.22); opacity: 0; }
        100% { transform: scale(1.22); opacity: 0; }
      }
      .edge-glow {
        fill: none;
        stroke-linecap: round;
        opacity: 0;
        transition: opacity 0.18s ease;
        pointer-events: none;
      }
      .edge-visible {
        fill: none;
        stroke-linecap: round;
        stroke-opacity: 0.82;
        transition: stroke-width 0.18s ease, stroke-opacity 0.18s ease;
      }
      .edge-flow {
        fill: none;
        stroke: #ffffff;
        stroke-width: 1.6px;
        stroke-linecap: round;
        stroke-dasharray: 4 10;
        opacity: 0;
        pointer-events: none;
      }
      @keyframes dashFlow {
        to { stroke-dashoffset: -28; }
      }
      .edge-hit {
        fill: none;
        stroke: transparent;
        stroke-width: 16px;
        cursor: pointer;
      }
      .edge-group:hover .edge-glow,
      .edge-group.selected .edge-glow {
        opacity: 0.32 !important;
      }
      .edge-group:hover .edge-visible,
      .edge-group.selected .edge-visible {
        stroke-opacity: 1 !important;
      }
      .edge-group:hover .edge-flow,
      .edge-group.selected .edge-flow {
        opacity: 0.75;
        animation: dashFlow 0.75s linear infinite;
      }
      .canvas-hint {
        position: absolute;
        bottom: 10px;
        left: 14px;
        font-size: 0.72rem;
        font-weight: 500;
        color: #64748b;
        background: rgba(255, 255, 255, 0.88);
        backdrop-filter: blur(6px);
        border: 1px solid #e2e8f0;
        padding: 4px 11px;
        border-radius: 999px;
        box-shadow: 0 2px 6px rgba(15, 23, 42, 0.04);
        pointer-events: none;
      }
    "))
  ),
  sidebar = sidebar(
    width = 345,
    bg = "#f8f9fa",
    # 1) 3-part NAICS subsector division (31, 32, 33)
    div(
      class = "naics-group-pills",
      radioButtons(
        inputId = "naics_sub",
        label = tags$span(style = "font-weight: 600;", "1. NAICS Manufacturing Group"),
        choices = c(
          "All (31–33)" = "all",
          "31 · Basic Mfg." = "31",
          "32 · Intermediate" = "32",
          "33 · Advanced Mfg." = "33"
        ),
        selected = "all",
        inline = TRUE
      )
    ),
    selectizeInput(
      inputId = "rama",
      label = tags$span(style = "font-weight: 600;", "2. Target Sector (Rama)"),
      choices = get_grouped_mfg_choices("all"),
      selected = "3361",
      options = list(
        placeholder = "Select or type sector...",
        selectOnTab = TRUE
      )
    ),
    hr(style = "margin: 0.65rem 0;"),
    sliderInput(
      inputId = "threshold",
      label = tags$span(style = "font-weight: 600;", "Direct Requirements Threshold (A)"),
      min = 0.005,
      max = 0.15,
      value = 0.01,
      step = 0.005,
      sep = ""
    ),
    sliderInput(
      inputId = "threshold2",
      label = tags$span(style = "font-weight: 600;", "Indirect Requirements Threshold (A²)"),
      min = 0.005,
      max = 0.20,
      value = 0.05,
      step = 0.005,
      sep = ""
    ),
    div(
      style = "display: flex; gap: 6px; flex-wrap: wrap; margin-top: 4px;",
      actionButton(
        inputId = "reset_layout_btn",
        label = "Reset Node Layout",
        class = "btn-sm btn-outline-dark",
        style = "flex: 1;"
      ),
      actionButton(
        inputId = "reset_btn",
        label = "Reset Sliders",
        class = "btn-sm btn-outline-secondary",
        style = "flex: 1;"
      ),
      downloadButton(
        outputId = "download_plot",
        label = "PNG",
        class = "btn-sm btn-outline-primary"
      )
    ),
    hr(style = "margin: 0.75rem 0;"),
    uiOutput("sector_summary")
  ),
  navset_card_underline(
    full_screen = TRUE,
    nav_panel(
      title = "Interactive Network (Sugiyama)",
      div(
        id = "net-container",
        div(
          id = "net-header",
          div(
            id = "net-title-wrap",
            div(id = "net-title", "Backward Linkages: Rama 3361"),
            div(id = "net-subtitle", "Threshold: 0.01   -   Threshold 2: 0.05   -   Depth: 2")
          ),
          tags$button(
            id = "canvas-reset-btn",
            onclick = "if (window.resetSugiyamaPositions) window.resetSugiyamaPositions();",
            HTML('<svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"><path d="M3 12a9 9 0 1 0 9-9 9.75 9.75 0 0 0-6.74 2.74L3 8"/><path d="M3 3v5h5"/></svg> Reset Node Positions')
          )
        ),
        div(
          id = "net-canvas-wrap",
          tags$svg(id = "net-svg"),
          div(
            class = "canvas-role-legend",
            div(class = "role-item", span(class = "role-dot", style = "background: #1a73e8; box-shadow: 0 0 0 2px #dbeafe;"), "Target Industry"),
            div(class = "role-item", span(class = "role-dot", style = "background: #0284c7;"), "Direct Supplier (A)"),
            div(class = "role-item", span(class = "role-dot", style = "background: #7c3aed;"), "Indirect Supplier (A²)")
          ),
          div(id = "net-tooltip"),
          div(
            class = "canvas-hint",
            "Drag any node to rearrange · Hover nodes for sector details · Hover or click edges for coefficient"
          )
        ),
        div(
          id = "net-footer",
          div(class = "colorbar-title", "Coefficient"),
          div(id = "colorbar-gradient", class = "colorbar-bar"),
          div(id = "colorbar-ticks", class = "colorbar-ticks")
        )
      ),
      # D3 Interactive Renderer Script (Google Material 3 / Curved Flow Visuals)
      tags$script(HTML("
        (function() {
          let currentData = null;
          let nodeMap = {};
          const svg = d3.select('#net-svg');
          const tooltip = d3.select('#net-tooltip');
          const wrap = document.getElementById('net-canvas-wrap');

          // Exact intersection with a capsule (rounded rectangle) perimeter so curved arrows
          // start and end cleanly at the pill border regardless of angle
          function clipToCapsule(cx, cy, towardX, towardY, halfW, halfH, extraPad) {
            const dx = towardX - cx;
            const dy = towardY - cy;
            const dist = Math.hypot(dx, dy);
            if (dist < 1e-4) return { x: cx, y: cy };
            const ux = dx / dist;
            const uy = dy / dist;
            const w = halfW + extraPad;
            const h = halfH + extraPad;
            // Ray-ellipse / rounded-box smooth approximation
            const t = 1 / Math.hypot(ux / w, uy / h);
            return {
              x: cx + ux * Math.min(dist * 0.45, t),
              y: cy + uy * Math.min(dist * 0.45, t)
            };
          }

          // Compute smooth cubic Bezier flow path between source and target capsules
          function computeCurvedEdge(sourceNode, targetNode) {
            const x1 = sourceNode.x, y1 = sourceNode.y;
            const x2 = targetNode.x, y2 = targetNode.y;
            const dx = x2 - x1;
            const dy = y2 - y1;

            // Control points for a graceful Sugiyama S-curve
            let c1x, c1y, c2x, c2y;
            if (Math.abs(dx) > 45) {
              c1x = x1 + dx * 0.42;
              c1y = y1 + dy * 0.10;
              c2x = x1 + dx * 0.58;
              c2y = y2 - dy * 0.10;
            } else {
              // Vertical or close neighbors: gentle lateral arc
              const bend = (y2 >= y1 ? 1 : -1) * 18;
              c1x = x1 + dx * 0.35 + bend;
              c1y = y1 + dy * 0.30;
              c2x = x1 + dx * 0.65 + bend;
              c2y = y1 + dy * 0.70;
            }

            const startPt = clipToCapsule(x1, y1, c1x, c1y, sourceNode.hw, sourceNode.hh, 2);
            const endPt = clipToCapsule(x2, y2, c2x, c2y, targetNode.hw, targetNode.hh, 5);

            const sdx = endPt.x - startPt.x;
            const sdy = endPt.y - startPt.y;
            const sc1x = startPt.x + sdx * 0.42;
            const sc1y = startPt.y + sdy * 0.10;
            const sc2x = startPt.x + sdx * 0.58;
            const sc2y = endPt.y - sdy * 0.10;

            return 'M ' + startPt.x.toFixed(2) + ',' + startPt.y.toFixed(2) +
                   ' C ' + sc1x.toFixed(2) + ',' + sc1y.toFixed(2) +
                   ' ' + sc2x.toFixed(2) + ',' + sc2y.toFixed(2) +
                   ' ' + endPt.x.toFixed(2) + ',' + endPt.y.toFixed(2);
          }

          function renderNetwork(data, animateEntrance = true) {
            currentData = data;
            svg.selectAll('*').remove();

            d3.select('#net-title').text('Backward Linkages: Rama ' + data.sector + ' (' + data.target_short + ')');
            d3.select('#net-subtitle').html(
              '<span class=\"sub-pill\">Direct (A) &ge; ' + data.threshold + '</span>' +
              '<span class=\"sub-pill\">Indirect (A&sup2;) &ge; ' + data.threshold2 + '</span>' +
              '<span class=\"sub-pill\">Depth: 2</span>'
            );

            const width = wrap.clientWidth || 860;
            const height = wrap.clientHeight || 530;

            // Define SVG patterns, gradients, and Material 3 elevation shadows
            const defs = svg.append('defs');

            // 1) Architectural dot-matrix background pattern
            const gridPattern = defs.append('pattern')
              .attr('id', 'google-dot-grid')
              .attr('width', 22)
              .attr('height', 22)
              .attr('patternUnits', 'userSpaceOnUse');
            gridPattern.append('circle')
              .attr('cx', 11)
              .attr('cy', 11)
              .attr('r', 1.15)
              .attr('fill', '#cbd5e1')
              .attr('opacity', 0.65);

            // Render dot-grid backdrop
            svg.append('rect')
              .attr('width', width)
              .attr('height', height)
              .attr('fill', 'url(#google-dot-grid)');

            if (data.empty || !data.nodes || data.nodes.length === 0) {
              svg.append('text')
                .attr('x', width / 2)
                .attr('y', height / 2)
                .attr('text-anchor', 'middle')
                .attr('fill', '#64748b')
                .attr('font-weight', '500')
                .attr('font-size', '15px')
                .text('No linkages found above selected thresholds. Try lowering the threshold sliders.');
              return;
            }

            // 2) Material 3 Card Drop Shadows
            const cardShadow = defs.append('filter')
              .attr('id', 'm3-shadow')
              .attr('x', '-25%').attr('y', '-25%')
              .attr('width', '150%').attr('height', '160%');
            cardShadow.append('feDropShadow')
              .attr('dx', 0).attr('dy', 3)
              .attr('stdDeviation', 4)
              .attr('flood-color', '#0f172a')
              .attr('flood-opacity', 0.10);

            const targetShadow = defs.append('filter')
              .attr('id', 'm3-target-shadow')
              .attr('x', '-35%').attr('y', '-35%')
              .attr('width', '170%').attr('height', '180%');
            targetShadow.append('feDropShadow')
              .attr('dx', 0).attr('dy', 4)
              .attr('stdDeviation', 6)
              .attr('flood-color', '#1a73e8')
              .attr('flood-opacity', 0.28);

            // 3) Target Node Google Blue Gradient
            const targetGrad = defs.append('linearGradient')
              .attr('id', 'google-blue-grad')
              .attr('x1', '0%').attr('y1', '0%')
              .attr('x2', '100%').attr('y2', '100%');
            targetGrad.append('stop').attr('offset', '0%').attr('stop-color', '#2563eb');
            targetGrad.append('stop').attr('offset', '100%').attr('stop-color', '#1d4ed8');

            // Update bottom Coefficient colorbar
            if (data.legend_gradient && data.legend_gradient.length > 0) {
              const gradCss = 'linear-gradient(to right, ' + data.legend_gradient.join(', ') + ')';
              document.getElementById('colorbar-gradient').style.background = gradCss;
            }
            const ticksContainer = document.getElementById('colorbar-ticks');
            ticksContainer.innerHTML = '';
            (data.legend_ticks || []).forEach(t => {
              const span = document.createElement('span');
              span.textContent = Number(t).toFixed(2);
              ticksContainer.appendChild(span);
            });

            // Scale Sugiyama coordinates to canvas with generous padding
            const padX = 78;
            const padY = 48;
            const hVals = data.nodes.map(d => d.orig_h);
            const vVals = data.nodes.map(d => d.orig_v);
            const minH = Math.min(...hVals), maxH = Math.max(...hVals);
            const minV = Math.min(...vVals), maxV = Math.max(...vVals);

            const scaleX = d3.scaleLinear()
              .domain([minH === maxH ? minH - 1 : minH, minH === maxH ? maxH + 1 : maxH])
              .range([padX, width - padX]);

            // Larger orig_v in ggraph + coord_flip() is at the top of the plot
            const scaleY = d3.scaleLinear()
              .domain([minV === maxV ? minV - 1 : maxV, minV === maxV ? maxV + 1 : minV])
              .range([padY, height - padY]);

            // Subtle vertical column guide lines for Sugiyama depth layers
            const uniqueCols = Array.from(new Set(hVals)).sort((a, b) => a - b);
            const guideLayer = svg.append('g').attr('class', 'guides-layer');
            uniqueCols.forEach(colVal => {
              const cx = scaleX(colVal);
              guideLayer.append('line')
                .attr('x1', cx).attr('y1', 22)
                .attr('x2', cx).attr('y2', height - 22)
                .attr('stroke', '#e2e8f0')
                .attr('stroke-width', 1)
                .attr('stroke-dasharray', '3 5')
                .attr('opacity', 0.8);
            });

            nodeMap = {};
            const nodes = data.nodes.map(d => {
              const isTarget = d.code === data.sector;
              const nx = scaleX(d.orig_h);
              const ny = scaleY(d.orig_v);
              const obj = Object.assign({}, d, {
                x: nx,
                y: ny,
                origX: nx,
                origY: ny,
                isTarget: isTarget,
                hw: isTarget ? 34 : 30,
                hh: isTarget ? 15.5 : 13.5
              });
              nodeMap[d.code] = obj;
              return obj;
            });

            // Scale edge width smoothly by coefficient weight
            const wScale = d3.scaleLinear()
              .domain([data.min_weight || 0.005, Math.max(data.max_weight || 0.3, 0.05)])
              .range([2.2, 4.5])
              .clamp(true);

            const edges = data.edges.map((e, idx) => Object.assign({}, e, {
              idx: idx,
              strokeW: wScale(e.weight),
              sourceNode: nodeMap[e.source],
              targetNode: nodeMap[e.target]
            }));

            // Define sculpted stealth chevron arrowheads matching each edge's Plasma color
            edges.forEach(e => {
              const marker = defs.append('marker')
                .attr('id', 'arrow-' + e.idx)
                .attr('viewBox', '-2 -5 12 10')
                .attr('refX', 6.8)
                .attr('refY', 0)
                .attr('markerWidth', 5.2)
                .attr('markerHeight', 5.2)
                .attr('orient', 'auto');

              marker.append('path')
                .attr('d', 'M -1,-3.8 L 8.2,0 L -1,3.8 L 1.4,0 Z')
                .attr('fill', e.color);
            });

            const edgeLayer = svg.append('g').attr('class', 'edges-layer');
            const nodeLayer = svg.append('g').attr('class', 'nodes-layer');

            const edgeGroups = edgeLayer.selectAll('.edge-group')
              .data(edges, d => d.id)
              .enter()
              .append('g')
              .attr('class', 'edge-group');

            // 1) Soft luminous halo under-path
            edgeGroups.append('path')
              .attr('class', 'edge-glow')
              .attr('stroke', d => d.color)
              .attr('stroke-width', d => d.strokeW + 5);

            // 2) Main Plasma-colored cubic Bezier path
            edgeGroups.append('path')
              .attr('class', 'edge-visible')
              .attr('stroke', d => d.color)
              .attr('stroke-width', d => d.strokeW)
              .attr('marker-end', d => 'url(#arrow-' + d.idx + ')');

            // 3) Animated flow overlay (visible on hover/selection)
            edgeGroups.append('path')
              .attr('class', 'edge-flow');

            // 4) Wide transparent hit area for effortless hover/click
            edgeGroups.append('path')
              .attr('class', 'edge-hit');

            function updateEdgePaths(selection) {
              selection.each(function(d) {
                const pathStr = computeCurvedEdge(d.sourceNode, d.targetNode);
                const g = d3.select(this);
                g.select('.edge-glow').attr('d', pathStr);
                g.select('.edge-visible').attr('d', pathStr);
                g.select('.edge-flow').attr('d', pathStr);
                g.select('.edge-hit').attr('d', pathStr);
              });
            }

            updateEdgePaths(edgeGroups);

            // Edge hover & click interactions
            edgeGroups
              .on('mouseenter', function(event, d) {
                d3.select(this).raise();
                tooltip.classed('visible', true).html(
                  '<div style=\"display:flex; justify-content:space-between; align-items:center; margin-bottom:5px; border-bottom:1px solid #f1f5f9; padding-bottom:4px;\">' +
                    '<span style=\"font-weight:700; font-size:0.88rem; color:#0f172a;\">' + d.source + ' &rarr; ' + d.target + '</span>' +
                    '<span style=\"font-size:0.68rem; font-weight:600; background:#f1f5f9; color:#475569; padding:2px 7px; border-radius:999px;\">' + d.round + '</span>' +
                  '</div>' +
                  '<div style=\"font-size:0.76rem; color:#475569; margin-bottom:7px; line-height:1.35;\">' +
                    '<div><strong style=\"color:#0f172a;\">Supplier:</strong> ' + d.source_short + '</div>' +
                    '<div><strong style=\"color:#0f172a;\">Buyer:</strong> ' + d.target_short + '</div>' +
                  '</div>' +
                  '<div style=\"display:flex; justify-content:space-between; align-items:center; background:#f8fafc; border:1px solid #e2e8f0; padding:5px 9px; border-radius:6px;\">' +
                    '<span style=\"font-size:0.75rem; font-weight:600; color:#475569;\">Coefficient (a<sub>ij</sub>)</span>' +
                    '<span style=\"font-weight:700; font-size:0.9rem; color:#0f172a;\">' + d.weight_fmt + ' <span style=\"color:#64748b; font-weight:500; font-size:0.78rem;\">(' + d.pct_fmt + ')</span></span>' +
                  '</div>'
                );
              })
              .on('mousemove', function(event) {
                const rect = wrap.getBoundingClientRect();
                let left = event.clientX - rect.left + 15;
                let top = event.clientY - rect.top + 15;
                if (left + 290 > rect.width) left -= 305;
                if (top + 135 > rect.height) top -= 125;
                tooltip.style('left', left + 'px').style('top', top + 'px');
              })
              .on('mouseleave', function() {
                tooltip.classed('visible', false);
              })
              .on('click', function(event, d) {
                event.stopPropagation();
                const isSelected = d3.select(this).classed('selected');
                edgeGroups.classed('selected', false);
                d3.select(this).classed('selected', !isSelected);
              });

            // Build Google Material 3 Capsule Nodes
            const nodeGroups = nodeLayer.selectAll('.node-group')
              .data(nodes, d => d.code)
              .enter()
              .append('g')
              .attr('class', 'node-group')
              .attr('transform', d => 'translate(' + d.x + ',' + d.y + ')');

            const nodeInner = nodeGroups.append('g')
              .attr('class', 'node-inner');

            // Subtle breathing pulse halo on Target Industry node
            nodeInner.filter(d => d.isTarget)
              .append('rect')
              .attr('class', 'target-pulse')
              .attr('x', d => -d.hw - 4)
              .attr('y', d => -d.hh - 4)
              .attr('width', d => (d.hw + 4) * 2)
              .attr('height', d => (d.hh + 4) * 2)
              .attr('rx', d => d.hh + 4)
              .attr('ry', d => d.hh + 4);

            // Elevated Capsule Body
            nodeInner.append('rect')
              .attr('class', 'node-capsule')
              .attr('x', d => -d.hw)
              .attr('y', d => -d.hh)
              .attr('width', d => d.hw * 2)
              .attr('height', d => d.hh * 2)
              .attr('rx', d => d.hh)
              .attr('ry', d => d.hh)
              .attr('fill', d => d.isTarget ? 'url(#google-blue-grad)' : '#ffffff')
              .attr('stroke', d => {
                if (d.isTarget) return '#1d4ed8';
                return d.anim_tier === 1 ? '#93c5fd' : '#cbd5e1';
              })
              .attr('stroke-width', d => d.isTarget ? 2 : 1.5)
              .attr('filter', d => d.isTarget ? 'url(#m3-target-shadow)' : 'url(#m3-shadow)');

            // Status accent indicator dot inside left side of capsule
            nodeInner.append('circle')
              .attr('cx', d => -d.hw + 11)
              .attr('cy', 0)
              .attr('r', d => d.isTarget ? 3.8 : 3.4)
              .attr('fill', d => {
                if (d.isTarget) return '#ffffff';
                return d.anim_tier === 1 ? '#0284c7' : '#7c3aed';
              });

            // 4-Digit NAICS Code Label
            nodeInner.append('text')
              .attr('class', 'node-label')
              .attr('x', 5)
              .attr('y', 0)
              .attr('fill', d => d.isTarget ? '#ffffff' : '#0f172a')
              .text(d => d.code);

            // Node hover interactions (highlight supply paths + show Material tooltip card)
            nodeGroups
              .on('mouseenter', function(event, d) {
                d3.select(this).raise();

                edgeGroups.each(function(e) {
                  const isConnected = (e.source === d.code || e.target === d.code);
                  const eg = d3.select(this);
                  eg.select('.edge-visible')
                    .style('stroke-opacity', isConnected ? 1 : 0.16)
                    .style('stroke-width', isConnected ? (e.strokeW + 1.2) + 'px' : e.strokeW + 'px');
                  eg.select('.edge-glow')
                    .style('opacity', isConnected ? 0.32 : 0);
                  eg.select('.edge-flow')
                    .style('opacity', isConnected ? 0.8 : 0)
                    .style('animation', isConnected ? 'dashFlow 0.75s linear infinite' : 'none');
                });

                const badgeBg = d.isTarget ? '#dbeafe' : (d.anim_tier === 1 ? '#e0f2fe' : '#f3e8ff');
                const badgeFg = d.isTarget ? '#1d4ed8' : (d.anim_tier === 1 ? '#0369a1' : '#6d28d9');

                tooltip.classed('visible', true).html(
                  '<div style=\"display:flex; justify-content:space-between; align-items:center; margin-bottom:5px;\">' +
                    '<span style=\"font-weight:700; font-size:0.92rem; color:#0f172a;\">Rama ' + d.code + '</span>' +
                    '<span style=\"font-size:0.68rem; background:' + badgeBg + '; color:' + badgeFg + '; padding:2px 8px; border-radius:999px; font-weight:700;\">' + d.role + '</span>' +
                  '</div>' +
                  '<div style=\"font-weight:600; color:#1e293b; margin-bottom:4px; font-size:0.83rem;\">' + d.short_name + '</div>' +
                  '<div style=\"font-size:0.73rem; color:#64748b; line-height:1.3; margin-bottom:8px;\">' + d.full_name + '</div>' +
                  '<div style=\"display:flex; justify-content:space-between; align-items:center; background:#f8fafc; border:1px solid #e2e8f0; padding:4px 8px; border-radius:6px; font-size:0.75rem;\">' +
                    '<span style=\"color:#475569; font-weight:600;\">Backward Multiplier (&mu;)</span>' +
                    '<span style=\"font-weight:700; color:#0f172a;\">' + Number(d.mu).toFixed(3) + '</span>' +
                  '</div>'
                );
              })
              .on('mousemove', function(event) {
                const rect = wrap.getBoundingClientRect();
                let left = event.clientX - rect.left + 15;
                let top = event.clientY - rect.top + 15;
                if (left + 295 > rect.width) left -= 310;
                if (top + 150 > rect.height) top -= 140;
                tooltip.style('left', left + 'px').style('top', top + 'px');
              })
              .on('mouseleave', function() {
                edgeGroups.select('.edge-visible')
                  .style('stroke-opacity', null)
                  .style('stroke-width', null);
                edgeGroups.select('.edge-glow').style('opacity', null);
                edgeGroups.select('.edge-flow')
                  .style('opacity', null)
                  .style('animation', null);
                tooltip.classed('visible', false);
              });

            // Smooth 60fps Drag Behavior
            const drag = d3.drag()
              .on('start', function(event, d) {
                d3.select(this).raise().classed('selected', true);
                tooltip.classed('visible', false);
              })
              .on('drag', function(event, d) {
                d.x = Math.max(d.hw + 8, Math.min(width - d.hw - 8, event.x));
                d.y = Math.max(d.hh + 8, Math.min(height - d.hh - 8, event.y));
                d3.select(this).attr('transform', 'translate(' + d.x + ',' + d.y + ')');
                updateEdgePaths(edgeGroups.filter(e => e.source === d.code || e.target === d.code));
              })
              .on('end', function(event, d) {
                d3.select(this).classed('selected', false);
              });

            nodeGroups.call(drag);

            svg.on('click', function() {
              edgeGroups.classed('selected', false);
            });

            // -----------------------------------------------------------------
            // Emphasized 1.0-Second (1000ms) Choreographed Entrance Animation
            // -----------------------------------------------------------------
            if (animateEntrance) {
              // 1) Nodes pop & glide in with elastic overshoot across the 1000ms timeline:
              //    Tier 0 (Target):   0ms -> 380ms
              //    Tier 1 (Direct): 260ms -> 680ms
              //    Tier 2 (Indirect): 540ms -> 1000ms
              nodeGroups
                .style('opacity', 0)
                .attr('transform', d => {
                  const shiftX = d.isTarget ? -18 : 22;
                  return 'translate(' + (d.x + shiftX) + ',' + d.y + ') scale(0.32)';
                })
                .transition()
                .delay(d => d.anim_tier === 0 ? 0 : (d.anim_tier === 1 ? 260 : 540))
                .duration(d => d.anim_tier === 0 ? 380 : (d.anim_tier === 1 ? 420 : 460))
                .ease(d3.easeBackOut.overshoot(1.35))
                .style('opacity', 1)
                .attr('transform', d => 'translate(' + d.x + ',' + d.y + ') scale(1)');

              // 2) Edges progressively draw their curved paths using stroke-dashoffset:
              //    Stage 1 (Direct A):    140ms -> 600ms
              //    Stage 2 (Indirect A²): 480ms -> 1000ms
              edgeGroups.each(function(d) {
                const g = d3.select(this);
                const visPath = g.select('.edge-visible');
                const pathNode = visPath.node();
                const totalLen = pathNode ? pathNode.getTotalLength() : 200;
                const delayMs = d.stage_num === 1 ? 140 : 480;
                const durMs = d.stage_num === 1 ? 460 : 520;

                g.style('opacity', 0);

                visPath
                  .attr('stroke-dasharray', totalLen + ' ' + totalLen)
                  .attr('stroke-dashoffset', totalLen);

                g.transition()
                  .delay(delayMs)
                  .duration(80)
                  .style('opacity', 1);

                visPath.transition()
                  .delay(delayMs)
                  .duration(durMs)
                  .ease(d3.easeCubicInOut)
                  .attr('stroke-dashoffset', 0)
                  .on('end', function() {
                    d3.select(this).attr('stroke-dasharray', null).attr('stroke-dashoffset', null);
                  });
              });
            }

            // Smoothly animate dragged nodes back to their original Sugiyama coordinates
            window.resetSugiyamaPositions = function() {
              const prevCoords = {};
              nodes.forEach(d => {
                prevCoords[d.code] = { x: d.x, y: d.y };
                d.x = d.origX;
                d.y = d.origY;
              });

              nodeGroups.transition()
                .duration(480)
                .ease(d3.easeCubicOut)
                .attr('transform', d => 'translate(' + d.origX + ',' + d.origY + ')');

              edgeGroups.transition()
                .duration(480)
                .ease(d3.easeCubicOut)
                .tween('pathTween', function(d) {
                  const g = d3.select(this);
                  const s0 = prevCoords[d.source], t0 = prevCoords[d.target];
                  const ixS = d3.interpolateNumber(s0.x, d.sourceNode.origX);
                  const iyS = d3.interpolateNumber(s0.y, d.sourceNode.origY);
                  const ixT = d3.interpolateNumber(t0.x, d.targetNode.origX);
                  const iyT = d3.interpolateNumber(t0.y, d.targetNode.origY);

                  return function(t) {
                    const tempS = { x: ixS(t), y: iyS(t), hw: d.sourceNode.hw, hh: d.sourceNode.hh };
                    const tempT = { x: ixT(t), y: iyT(t), hw: d.targetNode.hw, hh: d.targetNode.hh };
                    const pStr = computeCurvedEdge(tempS, tempT);
                    g.select('.edge-glow').attr('d', pStr);
                    g.select('.edge-visible').attr('d', pStr);
                    g.select('.edge-flow').attr('d', pStr);
                    g.select('.edge-hit').attr('d', pStr);
                  };
                });
            };
          }

          Shiny.addCustomMessageHandler('renderSugiyamaNet', function(data) {
            renderNetwork(data, true);
          });

          Shiny.addCustomMessageHandler('resetSugiyamaLayout', function(msg) {
            if (window.resetSugiyamaPositions) window.resetSugiyamaPositions();
          });

          window.addEventListener('resize', function() {
            if (currentData) renderNetwork(currentData, false);
          });
        })();
      "))
    ),
    nav_panel(
      title = "Linkage Table & Sector Directory",
      div(
        style = "padding: 8px; max-height: 630px; overflow-y: auto;",
        tableOutput("edges_table")
      )
    )
  )
)

# -------------------------------------------------------------------------
# 4. Server Logic
# -------------------------------------------------------------------------

server <- function(input, output, session) {

  # Filter dropdown choices into NAICS 31, 32, 33 groups when radio toggle changes
  observeEvent(input$naics_sub, {
    choices <- get_grouped_mfg_choices(input$naics_sub)
    flat_codes <- unlist(choices, use.names = FALSE)
    current_sel <- isolate(input$rama)
    new_sel <- if (!is.null(current_sel) && current_sel %in% flat_codes) {
      current_sel
    } else if ("3361" %in% flat_codes) {
      "3361"
    } else {
      flat_codes[1]
    }
    updateSelectizeInput(session, "rama", choices = choices, selected = new_sel, server = FALSE)
  }, ignoreInit = TRUE)

  # Reset sliders button
  observeEvent(input$reset_btn, {
    updateSliderInput(session, "threshold", value = 0.01)
    updateSliderInput(session, "threshold2", value = 0.05)
  })

  # Reset node positions button
  observeEvent(input$reset_layout_btn, {
    session$sendCustomMessage("resetSugiyamaLayout", list())
  })

  current_net <- reactive({
    req(input$rama)
    build_edges_rama(
      rama = input$rama,
      threshold = input$threshold,
      threshold2 = input$threshold2
    )
  })

  # Send Sugiyama network data to interactive D3 SVG canvas whenever inputs change
  observe({
    req(input$rama)
    payload <- build_interactive_network(
      rama = input$rama,
      threshold = input$threshold,
      threshold2 = input$threshold2
    )
    session$sendCustomMessage("renderSugiyamaNet", payload)
  })

  output$sector_summary <- renderUI({
    req(input$rama)
    info <- sector_catalog %>% filter(code == input$rama)
    net <- current_net()
    nodes_count <- length(unique(c(net$edges$row, net$edges$col)))

    div(
      style = "font-size: 0.83rem; color: #343a40; background: #ffffff; border: 1px solid #dee2e6; border-radius: 6px; padding: 10px 12px;",
      div(style = "font-weight: 700; color: #2b3e50; margin-bottom: 4px;", paste0("Rama ", info$code, " · ", info$short_name)),
      div(style = "color: #495057; margin-bottom: 8px; line-height: 1.3; font-size: 0.77rem;", info$full_name),
      div(
        style = "display: grid; grid-template-columns: 1fr 1fr; gap: 6px; border-top: 1px solid #eee; padding-top: 8px;",
        div(
          span(style = "color: #6c757d; display: block; font-size: 0.74rem;", "Multiplier (μ)"),
          span(style = "font-weight: 600;", sprintf("%.3f", info$mu))
        ),
        div(
          span(style = "color: #6c757d; display: block; font-size: 0.74rem;", "Active Nodes"),
          span(style = "font-weight: 600;", nodes_count)
        ),
        div(
          span(style = "color: #6c757d; display: block; font-size: 0.74rem;", "Direct Links"),
          span(style = "font-weight: 600;", net$n_direct)
        ),
        div(
          span(style = "color: #6c757d; display: block; font-size: 0.74rem;", "Indirect Links"),
          span(style = "font-weight: 600;", net$n_indirect)
        )
      )
    )
  })

  output$edges_table <- renderTable({
    req(input$rama)
    edges <- current_net()$edges
    if (nrow(edges) == 0) return(tibble(Message = "No linkages above selected thresholds."))

    edges %>%
      left_join(sector_catalog %>% select(code, supplier_name = short_name), by = c("row" = "code")) %>%
      left_join(sector_catalog %>% select(code, buyer_name = short_name), by = c("col" = "code")) %>%
      transmute(
        Stage = round,
        Supplier = paste0(row, " - ", supplier_name),
        Buyer = paste0(col, " - ", buyer_name),
        `Coeff (a_ij)` = sprintf("%.4f", weight)
      )
  }, striped = TRUE, hover = TRUE, bordered = FALSE, width = "100%")

  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("red_rama_", input$rama, "_t1_", input$threshold, "_t2_", input$threshold2, ".png")
    },
    content = function(file) {
      p <- graficar_red_rama(
        rama = input$rama,
        threshold = input$threshold,
        threshold2 = input$threshold2,
        panel_bg = "grey80"
      )
      ggsave(file, plot = p, width = 10, height = 7, dpi = 300)
    }
  )
}

shinyApp(ui = ui, server = server)
