#' Calculate LHS (Leaf-Height-Seed) Plant Ecological Strategy Classification
#'
#' This function implements the LHS plant ecology strategy scheme proposed by Westoby (1998),
#' which classifies plant species based on three key functional traits: specific leaf area (SLA),
#' canopy height at maturity, and seed mass. The LHS scheme provides a quantitative framework
#' for comparing plant ecological strategies worldwide.
#'
#' @description
#' The LHS scheme uses three fundamental plant traits that reflect important trade-offs
#' controlling plant strategies:
#' \itemize{
#'   \item \strong{Specific Leaf Area (SLA)}: Light-capturing area deployed per unit dry mass,
#'         reflecting the trade-off between rapid resource acquisition and leaf longevity
#'   \item \strong{Height}: Canopy height at maturity, expressing the amount of growth attempted
#'         between disturbances and competitive ability for light
#'   \item \strong{Seed Mass}: Reflecting the trade-off between seed number and individual seed
#'         provisioning, affecting dispersal capacity and seedling survival
#' }
#'
#' All three axes are log-scaled as they are approximately lognormally distributed between species.
#' Species are classified into eight strategy types based on whether their log-transformed trait
#' values are above (L = Large) or below (S = Small) the median values.
#'
#' @param data A data frame containing plant trait data with the following required columns:
#' \describe{
#'   \item{SLA}{Specific leaf area (area per unit dry mass, typically mm2/mg or m2/kg)}
#'   \item{Height}{Canopy height at maturity (typically in metres)}
#'   \item{SeedMass}{Seed mass (typically in mg or g)}
#' }
#' Row names should represent species names or identifiers.
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item Validates input data for required columns and checks for missing, zero, or negative values
#'   \item Log-transforms all three traits
#'   \item Calculates median values for each log-transformed trait
#'   \item Classifies each species based on whether traits are above (L) or below (S) medians
#'   \item Returns the original data with added log-transformed columns and strategy classification
#' }
#'
#' @return A data frame with the original columns plus:
#' \describe{
#'   \item{log_SLA}{Natural logarithm of SLA}
#'   \item{log_Height}{Natural logarithm of Height}
#'   \item{log_SeedMass}{Natural logarithm of SeedMass}
#'   \item{LHS_strategy}{Character string indicating the LHS strategy type (e.g., "S-L-S")}
#' }
#'
#' @references
#' 1. Westoby, M. (1998). A leaf-height-seed (LHS) plant ecology strategy scheme.
#' Plant and Soil, 199, 213–227.
#' 2. Yang, J., Wang, Z., Zheng, Y., & Pan, Y. (2022). Shifts in plant ecological strategies
#' in remnant forest patches along urbanization gradients. Forest Ecology and Management, 524, 120540.
#'
#' @examples
#' data(PFF)
#' pff <- PFF[, c("SLA", "Height", "SeedMass")]
#' rownames(pff) <- PFF$species
#' head(pff)
#' result <- LHS(pff)
#' head(result)
#'
#' @importFrom stats median
#' @export
LHS <- function(data) {
  # Check if required columns exist
  required_cols <- c("SLA", "Height", "SeedMass")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  # Check for NA values
  for (col in required_cols) {
    if (any(is.na(data[[col]]))) {
      stop(paste("Column", col, "contains NA values, please handle missing data first"))
    }
  }
  # Check for zero or negative values
  for (col in required_cols) {
    if (any(data[[col]] <= 0, na.rm = TRUE)) {
      stop(paste("Column", col, "contains zero or negative values, which are biologically unreasonable"))
    }
  }
  # Log-transform the three traits
  data$log_SLA <- log(data$SLA)
  data$log_Height <- log(data$Height)
  data$log_SeedMass <- log(data$SeedMass)
  # Calculate medians of log-transformed traits
  sla_median <- stats::median(data$log_SLA, na.rm = TRUE)
  height_median <- stats::median(data$log_Height, na.rm = TRUE)
  seedmass_median <- stats::median(data$log_SeedMass, na.rm = TRUE)
  # Define function to determine LHS strategy
  get_strategy <- function(log_sla, log_height, log_seedmass) {
    s <- ifelse(log_sla <= sla_median, "S", "L")
    h <- ifelse(log_height <= height_median, "S", "L")
    m <- ifelse(log_seedmass <= seedmass_median, "S", "L")
    paste(s, h, m, sep = "-")
  }
  # Apply function to data and add LHS strategy column
  data$LHS_strategy <- mapply(get_strategy, data$log_SLA, data$log_Height, data$log_SeedMass)
  return(data)
}
