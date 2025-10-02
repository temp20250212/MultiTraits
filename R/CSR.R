#' Classify Plant Strategies using Pierce et al. (2017) CSR Method
#'
#' This function calculates CSR (Competitor-Stress tolerator-Ruderal) ecological strategies
#' for plant species based on three key functional traits: Leaf Area (LA), Leaf Dry Matter
#' Content (LDMC), and Specific Leaf Area (SLA). The CSR classification system was originally
#' developed by Grime (1974) and this implementation follows the global method by Pierce et al. (2017).
#'
#' @param data A data frame containing at least three numeric columns:
#'   \code{LA}, \code{LDMC}, and \code{SLA}.
#'   \itemize{
#'     \item \strong{LA}: Leaf area in mm\eqn{^2}.
#'     \item \strong{LDMC}: Leaf dry matter content (%).
#'     \item \strong{SLA}: Specific leaf area in mm\eqn{^2}/mg.
#'   }
#'
#' @return A data frame containing the original input data plus four additional columns:
#'   \itemize{
#'     \item C: Competitor strategy proportion (0-100%)
#'     \item S: Stress-tolerator strategy proportion (0-100%)
#'     \item R: Ruderal strategy proportion (0-100%)
#'     \item type: CSR strategy type classification (character) - one of 19 possible categories:
#'       "C", "C/CR", "C/CS", "CR", "C/CSR", "CS", "CR/CSR", "CS/CSR", "R/CR", "CSR",
#'       "S/CS", "R/CSR", "S/CSR", "R", "SR/CSR", "S", "R/SR", "S/SR", "SR"
#'   }
#'
#' @details
#' The function implements the global CSR classification method which:
#' \enumerate{
#'   \item Transforms the three input traits using species-specific equations
#'   \item Calculates derived traits including leaf dry weight, fresh weight, and succulence index
#'   \item Adjusts LDMC for succulent species (>5 g dm-2)
#'   \item Projects traits onto principal component axes from a global calibration dataset
#'   \item Applies outlier corrections to keep values within calibrated ranges
#'   \item Converts to proportional CSR values that sum to 100%
#'   \item Assigns the closest matching tertiary CSR strategy type
#' }
#'
#' The three strategies represent different ecological approaches:
#' \itemize{
#'   \item \strong{C (Competitor)}: Species adapted to productive, undisturbed environments
#'   \item \strong{S (Stress-tolerator)}: Species adapted to resource-limited environments
#'   \item \strong{R (Ruderal)}: Species adapted to frequently disturbed environments
#' }
#'
#' @note
#' \itemize{
#'   \item All input values must be positive numbers
#'   \item NA values are not permitted and will cause the function to stop with an error
#' }
#'
#' @examples
#' LA <- c(369615.7, 11.8, 55.7, 36061.2, 22391.8, 30068.1, 31059.5, 29895.1)
#' LDMC <- c(25.2, 39.7, 13.3, 35.5, 33.2, 36.1, 35.2, 34.9)
#' SLA <- c(17.4, 6.6, 34.1, 14.5, 8.1, 12.1, 9.4, 10.9)
#' traits <- data.frame(LA, LDMC, SLA)
#' CSR(data = traits)
#'
#' @references
#' 1. Grime, J.P. (1974). Vegetation classification by reference to strategies. Nature, 250, 26–31.
#' 2. Pierce, S., Negreiros, D., Cerabolini, B.E.L., Kattge, J., Díaz, S., et al. (2017). A global method
#' for calculating plant CSR ecological strategies applied across biomes world-wide.
#' Functional Ecology, 31: 444-457.
#'
#' @export
CSR <- function(data) {
  # Check input data
  if (!is.data.frame(data)) {
    stop("The input must be a dataframe")
  }
  required_columns <- c("LA", "LDMC", "SLA")
  if (!all(required_columns %in% colnames(data))) {
    stop("The input data frame must contain LA, LDMC and SLA columns")
  }
  # Check for NA values in required columns
  na_check <- sapply(data[required_columns], function(x) any(is.na(x)))
  if (any(na_check)) {
    na_cols <- names(na_check)[na_check]
    stop(paste("NA values found in column(s):",
               paste(na_cols, collapse = ", "),
               ". Please remove or handle NA values before processing."))
  }
  # This function calculates the CSR strategy
  Strate_CSR <- function(LA, LDMC, SLA) {
    # transformations
    sqrt_LA <- sqrt(LA / 894205) * 100
    log_SLA <- log(SLA)
    # LDW leaf dry weight (mg) [after oven drying to constant weight]
    LDW <- LA/SLA
    # LFW leaf fresh weight (mg) [saturated, turgid fresh weight, after a night wrapped in damp paper and foil at 4 oC]
    LFW <- (100/LDMC)*LDW
    # Succulence index (g water dm-2) [For succulents (>5 g dm2), Leaf Water Content (%) (=100-LDMC) automatically used instead of LDMC for S calculation]
    Succulence_index  <- (LFW-LDW)/(LA/10)
    # LDMC leaf dry matter content adjusted for succulence
    LDMC_adjusted <- ifelse(Succulence_index>5, 100-(LDW*100)/LFW, (LDW*100)/LFW)
    logit_LDMC <- log((LDMC_adjusted / 100) / (1 - (LDMC_adjusted / 100)))

    # Regression against calibration PCA
    pca2_coord_c <- -0.8678 + 1.6464 * sqrt_LA
    pca1_coord_s <- 1.3369 + 0.000010019 * (1 - exp(-0.0000000000022303 * logit_LDMC)) + 4.5835 * (1 - exp(-0.2328 * logit_LDMC))
    pca1_coord_r <- -57.5924 + 62.6802 * exp(-0.0288 * log_SLA)
    min_c <- 0
    min_s <- -0.756451214853076
    min_r <- -11.3467682227961
    Neg_outlier_correction_C <- ifelse(pca2_coord_c < min_c, min_c, pca2_coord_c)
    Neg_outlier_correction_S <- ifelse(pca1_coord_s < min_s, min_s, pca1_coord_s)
    Neg_outlier_correction_R <- ifelse(pca1_coord_r < min_r, min_r, pca1_coord_r)
    max_c <- 57.3756711966087
    max_s <- 5.79158377609218
    max_r <- 1.10795515716546
    Pos_outlier_correction_C <- ifelse(Neg_outlier_correction_C > max_c, max_c, Neg_outlier_correction_C)
    Pos_outlier_correction_S <- ifelse(Neg_outlier_correction_S > max_s, max_s, Neg_outlier_correction_S)
    Pos_outlier_correction_R <- ifelse(Neg_outlier_correction_R > max_r, max_r, Neg_outlier_correction_R)

    # +ve translation (values ranges spanning zero are shifted to all become positive)
    C_positive_translation <- Pos_outlier_correction_C + abs(min_c)
    S_positive_translation <- Pos_outlier_correction_S + abs(min_s)
    R_positive_translation <- Pos_outlier_correction_R + abs(min_r)
    Range_PCA2_C_positive_translation <- max_c + abs(min_c)
    Range_PCA1_S_positive_translation <- max_s + abs(min_s)
    Range_PCA1_R_positive_translation <- max_r + abs(min_r)
    Proportion_of_total_variability_C <- C_positive_translation/Range_PCA2_C_positive_translation*100
    Proportion_of_total_variability_S <- S_positive_translation/Range_PCA1_S_positive_translation*100
    Proportion_of_total_variability_R <- 100-(R_positive_translation/Range_PCA1_R_positive_translation*100)
    Percentage_Conversion_Coefficient <- 100/sum(c(Proportion_of_total_variability_C, Proportion_of_total_variability_S, Proportion_of_total_variability_R))
    Percentage_Conversion_Coefficient

    # Tertiary CSR strategy
    C_values <- c(90, 73, 73, 48, 54, 48, 42, 42, 23, 33, 23, 23, 23, 5, 17, 5, 5, 5, 5)
    S_values <- c(5, 5, 23, 5, 23, 48, 17, 42, 5, 33, 73, 23, 54, 5, 42, 90, 23, 73, 48)
    R_values <- c(5, 23, 5, 48, 23, 5, 42, 17, 73, 33, 5, 54, 23, 90, 42, 5, 73, 23, 48)
    C_output <- Proportion_of_total_variability_C*Percentage_Conversion_Coefficient
    S_output <- Proportion_of_total_variability_S*Percentage_Conversion_Coefficient
    R_output <- Proportion_of_total_variability_R*Percentage_Conversion_Coefficient
    variances <- (C_output - C_values)^2 + (S_output - S_values)^2 + (R_output - R_values)^2
    min_variance_index <- which.min(variances)
    csr_type <- c("C", "C/CR", "C/CS", "CR", "C/CSR", "CS", "CR/CSR", "CS/CSR", "R/CR", "CSR", "S/CS", "R/CSR", "S/CSR", "R", "SR/CSR", "S", "R/SR", "S/SR", "SR")[min_variance_index]

    result <- list(c_proportion = C_output,s_proportion = S_output,r_proportion = R_output,csr_type = csr_type)
    return(result)
  }
  # Apply the Strate_CSR function to each row of data
  results <- t(apply(data, 1, function(row) {
    result <- Strate_CSR(LA = row["LA"], LDMC = row["LDMC"], SLA = row["SLA"])
    unlist(result)
  }))
  # Converting results to dataframes
  results_df <- as.data.frame(results)
  colnames(results_df) <- c("C", "S", "R", "type")
  results_df$C <- as.numeric( results_df$C)
  results_df$S <- as.numeric( results_df$S)
  results_df$R <- as.numeric( results_df$R)
  results <- cbind(data,results_df)
  return(results)
}


