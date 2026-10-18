# Interactive Backward Linkages Network App (Industry x Industry MIP 2018)
# Light-themed minimal Shiny GUI using ggraph ("sugiyama" layout)

library(shiny)
library(bslib)
library(tidyverse)
library(igraph)
library(ggraph)

# -------------------------------------------------------------------------
# 1. Data Loading & Pre-computation (runs once at app startup)
# -------------------------------------------------------------------------

ruta_mip <- "mip_csv/conjunto de datos"
file_mip_ixi <- "conjunto_de_datos_mip_t_pb_ixi_32018.csv"

mip_ixi_d_m <- read.csv(file_mip_ixi, encoding = "UTF-8")

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

# Intermediate demand matrix Z (263 x 263)
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

# Extract full descriptions from row descriptors
sector_catalog <- tibble(
  raw_desc = as.character(mip_ixi_d_m[z_idx_ixi - 3, 1])
) %>%
  mutate(
    code = str_extract(raw_desc, "^\\d{4}"),
    full_name = str_replace(raw_desc, "^\\d{4}---", ""),
    # Create clean abbreviated name for UI
    short_name = full_name %>%
      str_replace(
        "^Fabricación de |^Elaboración de |^Industria básica del |^Industria básica de |^Industrias de |^Industria de las |^Industria del |^Confección de |^Preparación y envasado de |^Preparación e hilado de ",
        ""
      ) %>%
      str_to_sentence() %>%
      str_trunc(52, ellipsis = "...")
  ) %>%
  left_join(mu_ixi, by = c("code" = "rama"))

# Dropdown choices: Manufacturing sectors (codes starting with '3')
mfg_catalog <- sector_catalog %>%
  filter(str_starts(code, "3")) %>%
  mutate(label = paste0(code, " - ", short_name))

mfg_choices <- setNames(mfg_catalog$code, mfg_catalog$label)

# Edge list representation of A_ixi (excluding self-loops row == col)
ixi_edgelist <- cbind("row" = colnames(A_ixi), as.data.frame(A_ixi)) %>%
  pivot_longer(
    cols = 2:last_col(),
    names_to = "col",
    values_to = "weight"
  ) %>%
  filter(row != col, weight > 0)

# -------------------------------------------------------------------------
# 2. Network Construction & Plotting Functions
# -------------------------------------------------------------------------

build_edges_rama <- function(rama, threshold = 0.01, threshold2 = 0.05) {
  sector <- as.character(rama)

  # Direct inputs (1st layer: suppliers -> target sector)
  edges_direct <- ixi_edgelist %>%
    filter(col == sector, weight >= threshold) %>%
    mutate(round = "Direct (A)")

  suppliers_direct <- edges_direct$row

  # First round of indirect inputs (2nd layer: tier-2 suppliers -> direct suppliers)
  edges_r1 <- ixi_edgelist %>%
    filter(col %in% suppliers_direct, weight >= threshold2) %>%
    mutate(round = "Indirect (A²)")

  # Combine and deduplicate (if an edge appears in both, keep unique row->col)
  edges_rama <- bind_rows(edges_direct, edges_r1) %>%
    distinct(row, col, .keep_all = TRUE)

  list(
    edges = edges_rama,
    n_direct = nrow(edges_direct),
    n_indirect = nrow(edges_r1),
    n_suppliers_direct = length(unique(suppliers_direct))
  )
}

graficar_red_rama <- function(rama, threshold = 0.01, threshold2 = 0.05, panel_bg = "grey85") {
  sector <- as.character(rama)
  net_data <- build_edges_rama(sector, threshold, threshold2)
  edges_rama <- net_data$edges

  sector_info <- sector_catalog %>% filter(code == sector)
  sector_title <- if (nrow(sector_info) > 0) {
    paste0("Backward Linkages: Rama ", sector, " (", sector_info$short_name[1], ")")
  } else {
    paste0("Backward Linkages: Rama ", sector)
  }

  # Handle empty graph gracefully when thresholds are too high
  if (nrow(edges_rama) == 0) {
    p_empty <- ggplot() +
      annotate(
        "text", x = 0, y = 0,
        label = paste0(
          "No input linkages found for sector ", sector, "\n",
          "at Direct Threshold >= ", threshold,
          " and Indirect Threshold >= ", threshold2, ".\n\n",
          "Try lowering the direct threshold slider."
        ),
        size = 4.5, color = "grey30", fontface = "italic"
      ) +
      labs(
        title = sector_title,
        subtitle = paste0(
          "Threshold (Direct): ", threshold,
          "   |   Threshold (Indirect): ", threshold2,
          "   |   Depth: 2"
        )
      ) +
      theme_void() +
      theme(
        panel.background = element_rect(fill = panel_bg, color = NA),
        plot.title = element_text(face = "bold", size = 14, margin = margin(b = 4)),
        plot.subtitle = element_text(face = "italic", size = 11, margin = margin(b = 10)),
        plot.margin = margin(12, 12, 12, 12)
      )
    return(p_empty)
  }

  # Build directed graph
  red <- graph_from_data_frame(d = edges_rama, directed = TRUE)

  plot_red <- ggraph(red, layout = "sugiyama") +
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

  return(plot_red)
}

# -------------------------------------------------------------------------
# 3. User Interface (Clean, Minimal Light Theme)
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
  sidebar = sidebar(
    width = 330,
    bg = "#f8f9fa",
    selectizeInput(
      inputId = "rama",
      label = tags$span(style = "font-weight: 600;", "Target Manufacturing Sector"),
      choices = mfg_choices,
      selected = "3361",
      options = list(
        placeholder = "Type code or sector name...",
        selectOnTab = TRUE
      )
    ),
    hr(style = "margin: 0.75rem 0;"),
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
      style = "display: flex; justify-content: space-between; align-items: center; margin-top: 4px;",
      actionButton(
        inputId = "reset_btn",
        label = "Reset Defaults",
        class = "btn-sm btn-outline-secondary"
      ),
      downloadButton(
        outputId = "download_plot",
        label = "Save PNG",
        class = "btn-sm btn-outline-primary"
      )
    ),
    hr(style = "margin: 0.85rem 0;"),
    uiOutput("sector_summary")
  ),
  navset_card_underline(
    full_screen = TRUE,
    nav_panel(
      title = "Network Visualization",
      plotOutput("network_plot", height = "650px")
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

  observeEvent(input$reset_btn, {
    updateSliderInput(session, "threshold", value = 0.01)
    updateSliderInput(session, "threshold2", value = 0.05)
  })

  current_net <- reactive({
    req(input$rama)
    build_edges_rama(
      rama = input$rama,
      threshold = input$threshold,
      threshold2 = input$threshold2
    )
  })

  output$sector_summary <- renderUI({
    req(input$rama)
    info <- sector_catalog %>% filter(code == input$rama)
    net <- current_net()
    nodes_count <- length(unique(c(net$edges$row, net$edges$col)))

    div(
      style = "font-size: 0.83rem; color: #343a40; background: #ffffff; border: 1px solid #dee2e6; border-radius: 6px; padding: 10px 12px;",
      div(style = "font-weight: 700; color: #2b3e50; margin-bottom: 4px;", paste0("Rama ", info$code)),
      div(style = "color: #495057; margin-bottom: 8px; line-height: 1.3;", info$full_name),
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

  output$network_plot <- renderPlot({
    req(input$rama)
    graficar_red_rama(
      rama = input$rama,
      threshold = input$threshold,
      threshold2 = input$threshold2,
      panel_bg = "grey80"
    )
  }, res = 115)

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
