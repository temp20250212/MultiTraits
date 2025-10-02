#' Classify Plant Strategies using Hodgson et al. (1999) CSR Method
#'
#' Implements the Hodgson et al. (1999) method for allocating plant species
#' into the CSR (Competitor–Stress-tolerator–Ruderal) triangle based on plant
#' functional traits. Also assigns each species to the nearest CSR type.
#'
#' @description
#' This function calculates C, S, and R scores as percentages based on input
#' plant trait data, following the approach of Hodgson et al. (1999) and its
#' application in Caccianiga et al. (2006). Input is a dataframe with specific
#' trait columns, and the output is a new dataframe containing calculated
#' CSR coordinates, percentages, and assigned CSR type.
#'
#' @param data A \code{data.frame} containing the following columns:
#' \describe{
#'   \item{growth_form}{Character vector: plant growth form, "g" for graminoid, "n" for non-graminoid.}
#'   \item{CH}{Numeric: Canopy height (mm).}
#'   \item{LDMC}{Numeric: Leaf dry matter content (percent).}
#'   \item{FP}{Numeric: Flowering period (# of months).}
#'   \item{LS}{Numeric: Lateral spread (six-point classification).}
#'   \item{LDW}{Numeric: Leaf dry weight (mg).}
#'   \item{SLA}{Numeric: Specific leaf area (mm2/mg).}
#'   \item{FS}{Numeric: Flowering start (month).}
#' }
#'
#' @details
#' This implementation:
#' \itemize{
#'   \item Uses different equations for graminoids and non-graminoids to compute raw CSR dimensions.
#'   \item Scales results to coordinate space \eqn{[-2.5, 2.5]}, then shifts to positive and converts to percentages.
#'   \item Assigns the nearest CSR type based on standard reference CSR percentages from Hodgson's scheme.
#' }
#'
#' @return
#' A \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{growth_form, CH, LDMC, FP, LS, LDW, SLA, FS} — copied from input;
#'   \item \code{C, S, R} — calculated CSR percentages;
#'   \item \code{type} — assigned CSR type label (e.g., "C", "CSR", "S/CSR").
#' }
#'
#' @note
#' Input data must not contain \code{NA} values in required columns. If such values
#' are present, the function will stop with an error.
#'
#' @references
#' 1. Hodgson, J.G., Wilson, P.J., Hunt, R., Grime, J.P. & Thompson, K. (1999).
#' Allocating CSR plant functional types: a soft approach to a hard problem. Oikos, 85, 282–294.
#' 2. Caccianiga, M., Luzzaro, A., Pierce, S., Ceriani, R.M. & Cerabolini, B. (2006).
#' The functional basis of a primary succession resolved by CSR classification. Oikos, 112, 10–20.
#'
#'
#' @examples
#' # Example trait dataset
#' traits <- data.frame(
#'   growth_form = c("g", "g", "n", "g", "n"),
#'   CH = c(45.3, 169.7, 13.7, 132.7, 76.0),
#'   LDMC = c(33.0, 37.9, 25.9, 28.0, 15.7),
#'   FP = c(2, 2, 2, 1, 2),
#'   LS = c(3, 5, 4, 2, 5),
#'   LDW = c(1.9, 9.9, 2.3, 7.5, 40.2),
#'   SLA = c(19.0, 20.4, 15.2, 22.6, 21.8),
#'   FS = c(5, 5, 4, 5, 5)
#' )
#'
#' # Run CSR classification
#' result <- CSR_hodgson(traits)
#' print(result)
#'
#' # Plot CSR positions
#' CSR_plot(data = result)
#'
#' @export
CSR_hodgson <- function(data) {
  # Check input data
  if (!is.data.frame(data)) {
    stop("The input must be a dataframe")
  }
  if(!all(c("growth_form","CH","LDMC","FP","LS","LDW","SLA","FS") %in% names(data))) {
    stop("Input data must contain columns: growth_form, CH, LDMC, FP, LS, LDW, SLA, FS")
  }
  # Check for NA values in required columns
  required_columns <- c("growth_form", "CH", "LDMC", "FP", "LS", "LDW", "SLA", "FS")
  na_check <- sapply(data[required_columns], function(x) any(is.na(x)))
  if (any(na_check)) {
    na_cols <- names(na_check)[na_check]
    stop(paste("NA values found in column(s):",
               paste(na_cols, collapse = ", "),
               ". Please remove or handle NA values before processing."))
  }
  results <- data
  # CH classification
  results$CH_class <- with(results, ifelse(CH > 999, 6,
                                           ifelse(CH > 599, 5,
                                                  ifelse(CH > 299, 4,
                                                         ifelse(CH > 99, 3,
                                                                ifelse(CH > 49, 2, 1))))))
  # Variable transformations
  results$LDMC_sqrt <- sqrt(results$LDMC)
  results$LDW_ln <- log(results$LDW) + 3
  results$SLA_sqrt <- sqrt(results$SLA)
  # Initialize Raw dimensions
  results$RawC <- results$RawS <- results$RawR <- NA
  # Graminoid (g)
  g_idx <- results$growth_form == "g"
  if(any(g_idx)){
    results$RawC[g_idx] <- 0.141 * results$CH_class[g_idx]^2 +
      0.09061 * results$LS[g_idx]^2
    results$RawS[g_idx] <- 54.6 -
      1.666 * results$CH_class[g_idx]^2 +
      1.069 * results$LDMC_sqrt[g_idx]^2 -
      2.732 * results$SLA_sqrt[g_idx]^2 +
      1.722 * results$LS[g_idx]^2
    results$RawR[g_idx] <- 2.518 * results$FP[g_idx] -
      2.748 * results$LDW_ln[g_idx] +
      5.37  * results$SLA_sqrt[g_idx]
  }
  # Non-graminoid (n)
  n_idx <- results$growth_form == "n"
  if(any(n_idx)){
    results$RawC[n_idx] <- 0.09245 * results$CH_class[n_idx]^2 +
      0.05631 * results$LS[n_idx]^2 +
      0.01595 * results$LDW_ln[n_idx]^2
    results$RawS[n_idx] <- -39.52 -
      7.581 * results$CH_class[n_idx] +
      2.633 * results$LDMC_sqrt[n_idx]^2 -
      0.351 * results$LDW_ln[n_idx]^2
    results$RawR[n_idx] <- -1.158 * results$LDMC_sqrt[n_idx]^2 +
      3.137 * results$FP[n_idx] +
      3.145 * results$FS[n_idx] -
      0.0849 * results$LDW_ln[n_idx]^2 -
      1.193  * results$SLA_sqrt[n_idx]^2 +
      11.4   * results$SLA_sqrt[n_idx]
  }
  # Map raw values to [-2.5, 2.5]
  results$C_coord <- -2.5 + 0.839  * results$RawC
  results$S_coord <- ifelse(g_idx,
                            -1.103 + 0.0474 * results$RawS,
                            -1.249 + 0.0531 * results$RawS)
  results$R_coord <- -2.5 + 0.119  * results$RawR
  # Limit to [-2.5, 2.5]
  results$C_coord <- pmin(pmax(results$C_coord, -2.5), 2.5)
  results$S_coord <- pmin(pmax(results$S_coord, -2.5), 2.5)
  results$R_coord <- pmin(pmax(results$R_coord, -2.5), 2.5)
  # Shift to positive values
  results$C_pos <- results$C_coord + 2.5
  results$S_pos <- results$S_coord + 2.5
  results$R_pos <- results$R_coord + 2.5
  # Convert to percentages
  sum_vals <- results$C_pos + results$S_pos + results$R_pos
  results$C_percent <- 100 * results$C_pos / sum_vals
  results$S_percent <- 100 * results$S_pos / sum_vals
  results$R_percent <- 100 * results$R_pos / sum_vals
  # Reference CSR strategy types
  C_values <- c(90, 73, 73, 48, 54, 48, 42, 42, 23, 33, 23, 23, 23, 5, 17, 5, 5, 5, 5)
  S_values <- c(5, 5, 23, 5, 23, 48, 17, 42, 5, 33, 73, 23, 54, 5, 42, 90, 23, 73, 48)
  R_values <- c(5, 23, 5, 48, 23, 5, 42, 17, 73, 33, 5, 54, 23, 90, 42, 5, 73, 23, 48)
  csr_types <- c("C", "C/CR", "C/CS", "CR", "C/CSR", "CS", "CR/CSR", "CS/CSR",
                 "R/CR", "CSR", "S/CS", "R/CSR", "S/CSR", "R", "SR/CSR",
                 "S", "R/SR", "S/SR", "SR")
  # Determine closest CSR type
  results$type <- sapply(1:nrow(results), function(i) {
    variances <- (results$C_percent[i] - C_values)^2 +
      (results$S_percent[i] - S_values)^2 +
      (results$R_percent[i] - R_values)^2
    min_variance_index <- which.min(variances)
    return(csr_types[min_variance_index])
  })
  a <- results[, c("growth_form", "CH", "LDMC", "FP", "LS", "LDW", "SLA", "FS",
                   "C_percent", "S_percent", "R_percent", "type")]
  colnames(a) <- c("growth_form", "CH", "LDMC", "FP", "LS", "LDW", "SLA", "FS","C", "S", "R", "type")
  return(a)
}

