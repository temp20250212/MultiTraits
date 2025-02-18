
library(scatterplot3d)

#' Calculate LHS (Leaf-Height-Seed) Strategy
#'
#' This function calculates the LHS strategy for plant species based on their SLA (Specific Leaf Area),
#' Height, and Seed Mass values.
#'
#' @param data A data frame containing columns for SLA, Height, and SeedMass.
#'
#' @return A data frame with an additional column 'LHS_strategy' containing the calculated LHS strategy.
#'
#' @details The function calculates median values for SLA, Height, and SeedMass from the input data.
#' It then determines the LHS strategy for each species based on whether its trait values are
#' less than or equal to (S) or greater than (L) the median values.
#'
#' The resulting strategy is a combination of three letters (S or L) representing Leaf (SLA),
#' Height, and Seed (SeedMass) respectively, separated by hyphens.
#'
#'
#' @references
#' 1. Westoby, M. (1998). A leaf-height-seed (LHS) plant ecology strategy scheme. Plant and Soil, 199, 213–227. https://doi.org/10.1023/A:1004327224729
#' 2. Yang, J., Wang, Z., Zheng, Y., & Pan, Y. (2022). Shifts in plant ecological strategies in remnant forest patches along urbanization gradients. Forest Ecology and Management, 524, 120540. https://doi.org/10.1016/j.foreco.2022.120540
#'
#' @examples
#' data(PFF)
#' pff <- PFF[, c("SLA", "Height", "SeedMass")]
#' result <- LHS(pff)
#' head(result)
#'
#' @importFrom stats median
#' @export
LHS <- function(data) {
  # Calculating the median SLA, Height and SeedMass
  sla_median <- stats::median(data$SLA, na.rm = TRUE)
  height_median <- stats::median(data$Height, na.rm = TRUE)
  seedmass_median <- stats::median(data$SeedMass, na.rm = TRUE)

  # Define a function to determine the LHS strategy
  get_strategy <- function(sla, height, seedmass) {
    s <- ifelse(sla <= sla_median, "S", "L")
    h <- ifelse(height <= height_median, "S", "L")
    m <- ifelse(seedmass <= seedmass_median, "S", "L")
    paste(s, h, m, sep = "-")
  }

  # Apply functions to data and add LHS strategy column
  data$LHS_strategy <- mapply(get_strategy, data$SLA, data$Height, data$SeedMass)

  return(data)
}


#' Generate a 3D scatterplot of plant traits
#'
#' This function creates a three-dimensional scatterplot of plant traits (Specific Leaf Area, Height, and Seed Mass) based on the Leaf-Height-Seed (LHS) plant ecology strategy scheme.
#'
#' @param data A data frame containing columns: SLA, Height, SeedMass, and LHS_strategy.
#' @param colors A vector of colors for different LHS strategies. Default is a predefined color palette.
#' @param log_transform Logical, indicating whether to log-transform the data. Default is TRUE.
#'
#' @return A 3D scatterplot of SLA, Height, and Seed Mass (optionally log-transformed), with points colored by LHS strategy.
#'
#' @details
#' The function performs the following steps:
#' Checks if the input data contains the required columns.
#' Converts LHS_strategy to a factor.
#' Optionally log-transforms the data based on the log_transform parameter.
#' Creates a 3D scatterplot using (potentially log-transformed) values of SLA, Height, and Seed Mass.
#' Colors points based on LHS strategy.
#' Adds a legend to identify different LHS strategies.
#'
#' @references
#' 1. Westoby, M. (1998). A leaf-height-seed (LHS) plant ecology strategy scheme. Plant and Soil, 199, 213–227. https://doi.org/10.1023/A:1004327224729
#' 2. Yang, J., Wang, Z., Zheng, Y., & Pan, Y. (2022). Shifts in plant ecological strategies in remnant forest patches along urbanization gradients. Forest Ecology and Management, 524, 120540. https://doi.org/10.1016/j.foreco.2022.120540
#'
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom graphics legend
#'
#' @examples
#' data(PFF)
#' pff <- PFF[, c("SLA", "Height", "SeedMass")]
#' result <- LHS(pff)
#' LHS_plot(result)
#' LHS_plot(result, log_transform = FALSE)
#'
#' @export
LHS_plot <- function(data,
                     colors = c("#30123BFF", "#4777EFFF", "#1BD0D5FF", "#62FC6BFF",
                                "#D2E935FF", "#FE9B2DFF", "#DB3A07FF", "#7A0403FF"),
                     log_transform = TRUE) {
  # Check required columns
  required_cols <- c("SLA", "Height", "SeedMass", "LHS_strategy")
  if (!all(required_cols %in% names(data))) {
    missing_cols <- setdiff(required_cols, names(data))
    stop(paste("Data is missing the following columns:", paste(missing_cols, collapse = ", ")))
  }
  # Convert LHS_strategy to factors
  data$LHS_strategy <- as.factor(data$LHS_strategy)
  # Determine if logarithmic conversion is performed
  transform_func <- if (log_transform) log else identity
  # Prepare shaft labels
  x_label <- if(log_transform) expression(log~"(SLA)") else "SLA"
  y_label <- if(log_transform) expression(log~"(Height)") else "Height"
  z_label <- if(log_transform) expression(log~"(Seed mass)") else "Seed mass"
  # Create 3D scatterplots
  s3d <- scatterplot3d::scatterplot3d(transform_func(data$SLA),
                                      transform_func(data$Height),
                                      transform_func(data$SeedMass),
                                      color = colors[as.numeric(data$LHS_strategy)],
                                      bg = colors[as.numeric(data$LHS_strategy)],
                                      xlab = x_label,
                                      ylab = y_label,
                                      zlab = z_label,
                                      type = "h",
                                      pch = 21,
                                      axis = TRUE,
                                      box = FALSE,
                                      grid = TRUE,
                                      angle = 60)
  # Add legend
  graphics::legend("topright",
                   legend = levels(data$LHS_strategy),
                   title = expression(bold("Types")),
                   cex = 0.8,
                   col = "black",
                   pt.bg = colors,
                   pch = 21,
                   box.lty = 1,
                   title.adj = 0.3,
                   y.intersp = 1.05,
                   x.intersp = 1,
                   yjust = 1,
                   bg = "white")
}


#' Create a table of Leaf-Height-Seed (LHS) strategy types
#'
#' This function generates a data frame containing different plant growth strategies
#' based on the Leaf-Height-Seed (LHS) scheme. Each strategy is described by a
#' combination of traits and their corresponding ecological interpretation.
#'
#' @return A data frame with two columns:
#'   \describe{
#'     \item{type}{Character vector of LHS strategy combinations (e.g., "L-L-L", "L-L-S", etc.)}
#'     \item{strategy}{Character vector describing the ecological strategy for each type}
#'   }
#'
#'
#' @references
#' 1. Westoby, M. (1998). A leaf-height-seed (LHS) plant ecology strategy scheme. Plant and Soil, 199, 213–227. https://doi.org/10.1023/A:1004327224729
#' 2. Yang, J., Wang, Z., Zheng, Y., & Pan, Y. (2022). Shifts in plant ecological strategies in remnant forest patches along urbanization gradients. Forest Ecology and Management, 524, 120540. https://doi.org/10.1016/j.foreco.2022.120540
#'
#' @examples
#' LHS_strategy_scheme()
#'
#' @export
LHS_strategy_scheme <- function(){
  LHS_strategy_table <- data.frame(
    type=c( "L-L-L","L-L-S",  "L-S-L","L-S-S", "S-L-L","S-L-S","S-S-L", "S-S-S"),
    strategy=c(
      "Rapid growth, strong survivability and competitiveness",
      "Rapid growth, strong survivability and weak competitiveness",
      "Rapid growth, long-distance dispersal and strong competitiveness",
      "Rapid growth, long-distance dispersal and weak competitiveness",
      "Slow growth, strong survivability and competitiveness",
      "Slow growth, strong survivability and weak competitiveness",
      "Slow growth, long-distance dispersal and strong competitiveness",
      "Slow growth, long-distance dispersal and weak competitiveness")
  )
  return(LHS_strategy_table)
}




