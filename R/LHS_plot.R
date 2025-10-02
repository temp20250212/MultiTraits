#' Create 3D Scatter Plot for LHS Plant Ecological Strategy Data
#'
#' This function creates a three-dimensional scatter plot to visualise the LHS (Leaf-Height-Seed)
#' plant ecological strategy scheme proposed by Westoby (1998). The plot displays species positions
#' in the three-dimensional LHS space defined by log-transformed specific leaf area (SLA),
#' canopy height at maturity, and seed mass.
#'
#' @description
#' The LHS scheme uses three fundamental plant traits that reflect important trade-offs
#' controlling plant ecological strategies :
#' \itemize{
#'   \item \strong{Specific Leaf Area (SLA)}: Light-capturing area deployed per unit dry mass,
#'         reflecting the trade-off between rapid resource acquisition and leaf longevity
#'   \item \strong{Height}: Canopy height at maturity, expressing the amount of growth attempted
#'         between disturbances and competitive ability for light
#'   \item \strong{Seed Mass}: Reflecting the trade-off between seed number and individual seed
#'         provisioning, affecting dispersal capacity and seedling survival
#' }
#'
#' All three axes are log-scaled as they are approximately lognormally distributed between
#' species. The 3D visualisation allows researchers to explore species clustering
#' and relationships within the LHS strategy space, facilitating the identification of
#' functional groups and ecological patterns.
#'
#' @param data A data frame containing LHS analysis results with the following required columns:
#' \describe{
#'   \item{log_SLA}{Natural logarithm of specific leaf area}
#'   \item{log_Height}{Natural logarithm of canopy height at maturity}
#'   \item{log_SeedMass}{Natural logarithm of seed mass}
#' }
#' Typically this would be the output from the \code{LHS()} function. Row names should
#' represent species names or identifiers.
#'
#' @param group Character string specifying the column name to use for grouping points by colour.
#' If \code{NULL} (default), all points will be plotted in the same colour. Common choices
#' include \code{"LHS_strategy"} to colour by the eight LHS strategy types (e.g., "S-S-S",
#' "L-L-L") or any other categorical variable in the dataset.
#'
#' @param show_cube Logical indicating whether to display the 3D cube frame around the plot.
#' Default is \code{TRUE}, which helps with spatial orientation and depth perception.
#'
#' @param colors Character vector of colours to use for different groups. Should contain at
#' least as many colours as there are levels in the grouping variable. Default provides
#' a viridis-inspired colour palette with 8 colours suitable for the 8 LHS strategy types.
#' If only one group is plotted, only the first colour will be used.
#'
#' @param cube_angle Numeric value specifying the viewing angle for the 3D plot in degrees.
#' Default is 60. Values between 40-80 typically provide good visualisation perspectives.
#'
#' @details
#' The function creates a 3D scatter plot using the \code{scatterplot3d} package, with
#' log-transformed trait values on each axis. When grouping is specified, the plot includes
#' an automated legend positioned to the right of the main plot area.
#'
#' The three axes represent the core dimensions of the LHS ecological strategy space:
#' \itemize{
#'   \item X-axis: log(SLA) - reflects the leaf economics spectrum from resource-conservative
#'         to resource-acquisitive strategies
#'   \item Y-axis: log(Height) - represents the plant size spectrum and competitive ability
#'         for light capture
#'   \item Z-axis: log(Seed Mass) - indicates the seed size spectrum affecting dispersal
#'         and establishment success
#' }
#'
#' The function automatically handles layout management when legends are displayed,
#' ensuring optimal use of plotting space.
#'
#' @return Invisibly returns the \code{scatterplot3d} object, which can be used for
#' adding additional graphical elements to the plot if needed.
#'
#'
#' @references
#' 1. Westoby, M. (1998). A leaf-height-seed (LHS) plant ecology strategy scheme.
#' Plant and Soil, 199, 213–227.
#' 2. Yang, J., Wang, Z., Zheng, Y., & Pan, Y. (2022). Shifts in plant ecological strategies
#' in remnant forest patches along urbanization gradients. Forest Ecology and Management, 524, 120540.
#'
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom graphics legend layout par plot.new
#'
#' @examples
#' data(PFF)
#' pff <- PFF[, c("SLA", "Height", "SeedMass")]
#' rownames(pff) <- PFF$species
#' head(pff)
#' result <- LHS(pff)
#' head(result)
#' LHS_plot(result)
#' LHS_plot(result, group = "LHS_strategy", show_cube = TRUE)
#' @export
LHS_plot <- function(data,
                     group = NULL,
                     show_cube = TRUE,
                     colors = c("#30123BFF", "#4777EFFF", "#1BD0D5FF", "#62FC6BFF",
                                "#D2E935FF", "#FE9B2DFF", "#DB3A07FF", "#7A0403FF"),
                     cube_angle = 60
) {
  # Check required columns
  required_cols <- c("log_SLA", "log_Height", "log_SeedMass")
  if (!all(required_cols %in% names(data))) {
    missing_cols <- setdiff(required_cols, names(data))
    stop(paste("Data is missing the following columns:", paste(missing_cols, collapse = ", ")))
  }
  # Check if group column exists
  if (!is.null(group) && !group %in% names(data)) {
    stop(paste("Group column", group, "not found in data"))
  }
  # Prepare data
  x_label <- expression(log~"(SLA)")
  y_label <- expression(log~"(Height)")
  z_label <- expression(log~"(Seed mass)")
  x_data <- data$log_SLA
  y_data <- data$log_Height
  z_data <- data$log_SeedMass
  if (is.null(group)) {
    plot_colors <- rep(colors[1], nrow(data))
    show_legend <- FALSE
  } else {
    group_factor <- as.factor(data[[group]])
    plot_colors <- colors[as.numeric(group_factor)]
    show_legend <- TRUE
  }
  if (show_legend) {
    # Set layout: main plot takes 80% width, legend takes 20% width
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
    # Main plot area
    graphics::par(mar = c(5, 3, 4, 1))
  } else {
    graphics::par(mar = c(5, 3, 4, 3))
  }
  # Create 3D scatter plot
  s3d <- scatterplot3d::scatterplot3d(
    x = x_data,
    y = y_data,
    z = z_data,
    color = plot_colors,
    pch = 21,
    bg = plot_colors,
    xlab = x_label,
    ylab = y_label,
    zlab = z_label,
    type = "h",
    axis = TRUE,
    tick.marks = TRUE,
    label.tick.marks = TRUE,
    grid = TRUE,
    box = show_cube,
    angle = cube_angle,
    cex.symbols = 1
  )
  # Add legend in separate area on the right
  if (show_legend) {
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::legend("left",
                     legend = levels(group_factor),
                     title = "Types",
                     cex = 0.9,
                     col = "black",
                     pt.bg = colors[seq_along(levels(group_factor))],
                     pch = 21,
                     bty = "o",        # Add border ("o" means border is present)
                     title.adj = 0,    # Title left-aligned
                     title.font = 2,   # Title bold (1=normal, 2=bold, 3=italic, 4=bold italic)
                     adj = 0)          # Legend text left-aligned

    # Reset layout
    graphics::layout(1)
  }
  return(invisible(s3d))
}
