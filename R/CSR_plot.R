#' Create a Ternary Plot for CSR Plant Ecological Strategies
#'
#' This function creates a ternary plot to visualize plant ecological strategies
#' based on the CSR (Competitor-Stress tolerator-Ruderal) framework developed by
#' Grime (1974). The plot is built using `ggplot2` and displays the relative
#' proportions of C, S, and R strategies for each species or sample.
#'
#' @param data A data frame containing CSR strategy data. Must include columns:
#'   \describe{
#'     \item{C}{Numeric vector of Competitor strategy values (0-100)}
#'     \item{S}{Numeric vector of Stress-tolerator strategy values (0-100)}
#'     \item{R}{Numeric vector of Ruderal strategy values (0-100)}
#'     \item{type}{Character vector indicating the CSR strategy type/classification}
#'   }
#' @param point_size Numeric value specifying the size of points in the plot.
#'   Default is 3.
#' @param point_shape Numeric value specifying the shape of points in the plot.
#'   Default is 21 (filled circle with border).
#' @param custom_colors Named character vector specifying custom colors for each
#'   CSR strategy type. Default includes 19 predefined colors for all possible
#'   CSR combinations.
#'
#' @return A ggplot object representing a ternary plot with:
#'   \itemize{
#'     \item Points colored by CSR strategy type
#'     \item Ternary coordinate system with C, S, R axes
#'     \item Legend showing strategy types and their colors
#'     \item Grid lines and arrows for better visualization
#'   }
#'
#' @details
#' The CSR strategy framework classifies plants into three primary functional
#' types based on their ecological strategies:
#' \describe{
#'   \item{C (Competitors)}{Species adapted to productive, low-stress environments}
#'   \item{S (Stress-tolerators)}{Species adapted to unproductive, high-stress environments}
#'   \item{R (Ruderals)}{Species adapted to productive, high-disturbance environments}
#' }
#'
#' The ternary plot allows visualization of the relative contribution of each
#' strategy, where each point represents a species positioned according to its
#' C, S, and R values (which sum to 100%).
#'
#' @seealso
#' \code{\link{CSR}} or \code{\link{CSR_hodgson}} for calculating CSR strategies from plant functional traits
#'
#' @importFrom ggplot2 aes arrow coord_fixed element_rect element_text geom_line geom_point geom_polygon geom_segment geom_text guide_legend guides ggplot margin scale_fill_manual theme theme_void unit xlim ylim
#'
#' @references
#' 1. Grime, J.P. (1974). Vegetation classification by reference to strategies. Nature, 250, 26–31.
#' 2. Hodgson, J.G., Wilson, P.J., Hunt, R., Grime, J.P. & Thompson, K. (1999).
#' Allocating CSR plant functional types: a soft approach to a hard problem. Oikos, 85, 282–294.
#' 3. Caccianiga, M., Luzzaro, A., Pierce, S., Ceriani, R.M. & Cerabolini, B. (2006).
#' The functional basis of a primary succession resolved by CSR classification. Oikos, 112, 10–20.
#' 4. Pierce, S., Negreiros, D., Cerabolini, B.E.L., Kattge, J., Díaz, S., et al. (2017).
#' A global method for calculating plant CSR ecological strategies applied across biomes world-wide.
#' Functional Ecology, 31: 444-457.
#'
#' @examples
#' data(PFF)
#' head(PFF)
#' traits <- data.frame(LA=PFF$Leaf_area, LDMC=PFF$LDMC, SLA=PFF$SLA)
#' head(traits)
#' result <- CSR(data = traits)
#' head(result)
#' CSR_plot(data=result)
#'
#' @export
CSR_plot <- function(
    data,
    point_size = 3,
    point_shape = 21,
    custom_colors = c(
      "C" = "#E60D0D", "C/CR" = "#BA0D3B", "C/CS" = "#BA3B0D", "CR" = "#7A0D7A",
      "C/CSR" = "#8A3B3B", "CS" = "#7A7A0D", "CR/CSR" = "#6B2B6B",
      "CS/CSR" = "#6B6B2B", "R/CR" = "#3B0DBA", "CSR" = "#545454",
      "S/CS" = "#3BBA0D", "R/CSR" = "#3B3B8A", "S/CSR" = "#3B8A3B",
      "R" = "#0D0DE6", "SR/CSR" = "#2B6B6B", "S" = "#0DE60D",
      "R/SR" = "#0D3BBA", "S/SR" = "#0DBA3B", "SR" = "#0D7A7A")
){
  # Convert to proportions
  data$R_prop <- data$R / 100
  data$C_prop <- data$C / 100
  data$S_prop <- data$S / 100
  # Coordinate transformation for ternary plot
  height <- sqrt(3) / 2
  data$x <- data$S_prop + 0.5 * data$C_prop
  data$y <- height * data$C_prop
  # Create triangle boundary
  triangle <- data.frame(
    x = c(0, 1, 0.5, 0),
    y = c(0, 0, height, 0)
  )
  # Create grid lines data
  grid_lines <- list()
  line_id <- 1
  # Grid line intervals (20%, 40%, 60%, 80%)
  grid_props <- c(0.2, 0.4, 0.6, 0.8)
  # Lines parallel to R-S edge (constant C value)
  for (c_prop in grid_props) {
    x_left <- 0.5 * c_prop
    y_left <- height * c_prop
    x_right <- 1 - 0.5 * c_prop
    y_right <- height * c_prop
    grid_lines[[length(grid_lines) + 1]] <- data.frame(
      x = c(x_left, x_right),
      y = c(y_left, y_right),
      line_id = line_id,
      type = "C_constant"
    )
    line_id <- line_id + 1
  }
  # Lines parallel to C-S edge (constant R value)
  for (r_prop in grid_props) {
    x_bottom <- r_prop
    y_bottom <- 0
    x_top <- r_prop + 0.5 * (1 - r_prop)
    y_top <- height * (1 - r_prop)

    grid_lines[[length(grid_lines) + 1]] <- data.frame(
      x = c(x_bottom, x_top),
      y = c(y_bottom, y_top),
      line_id = line_id,
      type = "R_constant"
    )
    line_id <- line_id + 1
  }
  # Lines parallel to C-R edge (constant S value)
  for (s_prop in grid_props) {
    x_bottom <- 1 - s_prop
    y_bottom <- 0
    x_top <- 0.5 * (1 - s_prop)
    y_top <- height * (1 - s_prop)

    grid_lines[[length(grid_lines) + 1]] <- data.frame(
      x = c(x_bottom, x_top),
      y = c(y_bottom, y_top),
      line_id = line_id,
      type = "S_constant"
    )
    line_id <- line_id + 1
  }
  # Merge all grid lines
  all_grid_lines <- do.call(rbind, grid_lines)
  # Axis label positions
  axis_labels <- data.frame(
    x = c(0.5, 1.1, -0.1),
    y = c(height + 0.08, -0.08, -0.08),
    label = c("C (%)", "S (%)", "R (%)"),
    hjust = c(0.5, 0.5, 0.5),
    vjust = c(0.5, 1.2, 1.2)
  )
  # Tick labels
  tick_labels <- list()
  tick_values <- seq(0, 100, 20)
  for (val in tick_values) {
    prop <- val / 100
    # C-axis tick
    x_c <- 0.5 * prop - 0.05
    y_c <- height * prop
    tick_labels[[length(tick_labels) + 1]] <- data.frame(
      x = x_c, y = y_c, label = val, axis = "C"
    )
    # S-axis tick
    x_s <- 1 - 0.5 * prop + 0.05
    y_s <- height * prop
    tick_labels[[length(tick_labels) + 1]] <- data.frame(
      x = x_s, y = y_s, label = 100 - val, axis = "S"
    )
    # R-axis tick
    x_r <- prop
    y_r <- -0.05
    tick_labels[[length(tick_labels) + 1]] <- data.frame(
      x = x_r, y = y_r, label = 100 - val, axis = "R"
    )
  }
  all_tick_labels <- do.call(rbind, tick_labels)
  # Arrow data
  arrow_length <- 0.6
  offset <- 0.15
  # R-axis direction arrow
  arrow_data <- data.frame(
    x = c(0.2 + arrow_length),
    y = c(-offset + 0.05),
    xend = c(0.2),
    yend = c(-offset + 0.05),
    axis = c("R_direction")
  )
  # C-axis direction arrow
  arrow_data <- rbind(arrow_data, data.frame(
    x = c(0.15 - offset * cos(pi/6)),
    y = c(height * 0.25),
    xend = c(0.15 - offset * cos(pi/6) + arrow_length * cos(pi/3)),
    yend = c(height * 0.25 + arrow_length * sin(pi/3)),
    axis = c("S_direction")
  ))
  # S-axis direction arrow
  arrow_data <- rbind(arrow_data, data.frame(
    xend = c(0.85 + offset * cos(pi/6)),
    yend = c(height * 0.25),
    x = c(0.85 + offset * cos(pi/6) + arrow_length * cos(2*pi/3)),
    y = c(height * 0.25 + arrow_length * sin(2*pi/3)),
    axis = c("C_direction")
  ))
  # Arrow labels
  arrow_labels <- data.frame(
    x = c(0.2 + arrow_length/2,
          0.15 - offset * cos(pi/6) + arrow_length * cos(pi/3)/2 - 0.03,
          0.85 + offset * cos(pi/6) + arrow_length * cos(2*pi/3)/2 + 0.03),
    y = c(-offset + 0.05 - 0.03,
          height * 0.25 + arrow_length * sin(pi/3)/2 + 0.02,
          height * 0.25 + arrow_length * sin(2*pi/3)/2 + 0.02),
    label = c("R (%)", "C (%)", "S (%)"),
    axis = c("R_label", "C_label", "S_label")
  )
  # Create ggplot
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = triangle,
                          aes(x = .data$x, y = .data$y),
                          fill = NA, color = "black", linewidth = 1) +
    ggplot2::geom_line(data = all_grid_lines,
                       aes(x = .data$x, y = .data$y, group = .data$line_id),
                       color = "gray80", linewidth = 0.3) +
    ggplot2::geom_segment(data = arrow_data[1,],
                          aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm"),
                                                 type = "open"),
                          color = "blue",
                          linewidth = 1) +
    ggplot2::geom_segment(data = arrow_data[2,],
                          aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm"),
                                                 type = "open"),
                          color = "red",
                          linewidth = 1) +
    ggplot2::geom_segment(data = arrow_data[3,],
                          aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm"),
                                                 type = "open"),
                          color = "green",
                          linewidth = 1) +
    # R label
    ggplot2::geom_text(data = arrow_labels[1,],
                       aes(x = .data$x, y = .data$y, label = .data$label),
                       size = 4, color = "blue", fontface = "bold",
                       family = "serif", angle = 0) +
    # C label
    ggplot2::geom_text(data = arrow_labels[2,],
                       aes(x = .data$x, y = .data$y, label = .data$label),
                       size = 4, color = "red", fontface = "bold",
                       family = "serif", angle = 60) +
    # S label
    ggplot2::geom_text(data = arrow_labels[3,],
                       aes(x = .data$x, y = .data$y, label = .data$label),
                       size = 4, color = "green", fontface = "bold",
                       family = "serif", angle = -60) +
    ggplot2::geom_point(data = data,
                        aes(x = .data$x, y = .data$y, fill = .data$type),
                        size = point_size, shape = point_shape, color = "black") +
    ggplot2::geom_text(data = axis_labels,
                       aes(x = .data$x, y = .data$y, label = .data$label),
                       size = 5, family = "serif", fontface = "bold",
                       color = c("red", "green", "blue")) +
    ggplot2::geom_text(data = all_tick_labels,
                       aes(x = .data$x, y = .data$y, label = .data$label),
                       size = 4, family = "serif") +
    ggplot2::scale_fill_manual(values = custom_colors, name = "Types") +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      legend.background = ggplot2::element_rect(color = "black", linewidth = 0.5),
      legend.margin = ggplot2::margin(4, 4, 4, 4),
      legend.key = ggplot2::element_rect(fill = NA, color = NA),
      legend.key.size = ggplot2::unit(0.4, "cm"),
      text = ggplot2::element_text(size = 12, family = "serif"),
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 9, face = "bold"),
      legend.spacing.y = ggplot2::unit(0.2, "cm"),
      legend.box = "vertical",
      plot.margin = ggplot2::margin(20, 20, 20, 20)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      ncol = ifelse(length(unique(data$type)) > 10, 2, 1),
      byrow = TRUE
    )) +
    ggplot2::xlim(-0.1, 1.1) +
    ggplot2::ylim(-0.13, height + 0.1)
  return(p1)
}
