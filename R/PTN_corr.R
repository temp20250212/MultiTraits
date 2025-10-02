#' Calculate and Visualize Plant Trait Correlation Network
#'
#' This function calculates correlation coefficients for given plant traits and generates a correlation network plot.
#' It supports both standard correlation analysis and phylogenetically corrected correlation analysis using
#' phylogenetic independent contrasts.
#'
#' @param traits_matrix A numeric matrix or data frame where each column represents a plant trait and each row
#'   represents a sample (species). Row names should contain species names when using phylogenetic correction.
#' @param rThres Numeric, threshold for correlation coefficient, default is 0.2. Only correlations with
#'   absolute values above this threshold will be displayed in the plot. Must be between 0 and 1.
#' @param pThres Numeric, threshold for p-value, default is 0.05. Only correlations with p-values below
#'   this threshold will be displayed in the plot. Must be between 0 and 1.
#' @param method Character, specifies the correlation method to use: "pearson" (default) or "spearman".
#' @param phylo_correction Logical, whether to apply phylogenetic correction using phylogenetic independent
#'   contrasts. Default is FALSE.
#' @param phylo_tree A phylo object (from the ape package) containing the phylogenetic tree. Required when
#'   phylo_correction = TRUE. Species names in the tree must match row names in traits_matrix.
#'
#' @return Returns a correlation network plot object created by corrplot. The plot displays correlations
#'   as circles, with positive correlations in blue and negative correlations in red. Non-significant
#'   correlations (based on thresholds) are marked with red crosses.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters and data structure
#'   \item If phylogenetic correction is requested, matches species between traits_matrix and phylo_tree
#'   \item Calculates correlation coefficients using either standard correlation or phylogenetic independent contrasts
#'   \item Adjusts p-values using the False Discovery Rate (FDR) method
#'   \item Creates a correlation network plot using hierarchical clustering for trait ordering
#'   \item Marks non-significant correlations with red crosses based on both correlation and p-value thresholds
#' }
#'
#' When phylo_correction = TRUE, the function uses phylogenetic independent contrasts to account for
#' phylogenetic relationships among species, which helps control for the non-independence of species
#' data due to shared evolutionary history.
#'
#' @note
#' \itemize{
#'   \item Missing values (NA) should be handled before using this function
#'   \item When using phylogenetic correction, ensure species names are consistent between traits_matrix
#'     row names and phylo_tree tip labels
#'   \item The function requires the Hmisc, corrplot, and (optionally) ape packages
#' }
#'
#' @references
#' 1. Felsenstein, J. (1985). Phylogenies and the comparative method. The American Naturalist, 125(1), 1-15.
#'
#' 2. He, N., Li, Y., Liu, C., et al. (2020). Plant trait networks: improved resolution of the dimensionality of adaptation.
#' Trends in Ecology & Evolution, 35(10), 908-918.
#'
#' 3. Li, Y., Liu, C., Sack, L., Xu, L., Li, M., Zhang, J., & He, N. (2022). Leaf trait network architecture shifts with
#' species‐richness and climate across forests at continental scale. Ecology Letters, 25(6), 1442-1457.
#'
#'
#' @importFrom Hmisc rcorr
#' @importFrom stats p.adjust
#' @importFrom corrplot corrplot
#' @importFrom ape keep.tip
#'
#' @examples
#' data(PFF)
#' rownames(PFF) <- PFF$species
#' PFF_traits <- PFF[, c("Height", "Leaf_area","LDMC","SLA","SRL","SeedMass","FltDate",
#'                       "FltDur","Leaf_Cmass","Leaf_Nmass","Leaf_CN","Leaf_Pmass",
#'                       "Leaf_NP","Leaf_CP","Root_Cmass","Root_Nmass","Root_CN")]
#' PFF_traits <- na.omit(PFF_traits)
#' head(PFF_traits)
#' PTN_corr(traits_matrix = PFF_traits, rThres = 0.2, pThres = 0.05, method = "pearson")
#'
#' data(PFF_tree)
#' PTN_corr(traits_matrix = PFF_traits, rThres = 0.2, pThres = 0.05, method = "pearson",
#'         phylo_correction = TRUE, phylo_tree = PFF_tree)
#'
#' @export
PTN_corr <- function(traits_matrix, rThres = 0.2, pThres = 0.05, method = "pearson",
                    phylo_correction = FALSE, phylo_tree = NULL) {
  # Validate the method parameter
  method <- match.arg(method, choices = c("pearson", "spearman"))
  # Input validation for traits_matrix
  if (!is.data.frame(traits_matrix) && !is.matrix(traits_matrix)) {
    stop("traits_matrix must be a dataframe or matrix")
  }
  # Phylogenetic correction parameter validation
  if (phylo_correction && is.null(phylo_tree)) {
    stop("phylo_tree must be provided when phylo_correction = TRUE")
  }
  if (phylo_correction) {
    if (!inherits(phylo_tree, "phylo")) {
      stop("phylo_tree must be a phylo object from the ape package")
    }
    # Check if species names in phylogenetic tree match traits_matrix row names
    if (!all(rownames(traits_matrix) %in% phylo_tree$tip.label)) {
      missing_species <- rownames(traits_matrix)[!rownames(traits_matrix) %in% phylo_tree$tip.label]
      stop(paste("The following species are not present in phylo_tree:",
                 paste(missing_species, collapse = ", ")))
    }
    # Ensure data and phylogenetic tree species order are consistent
    common_species <- intersect(rownames(traits_matrix), phylo_tree$tip.label)
    traits_matrix <- traits_matrix[common_species, ]
    phylo_tree <- ape::keep.tip(phylo_tree, common_species)
  }
  # Calculate correlation matrix
  if (phylo_correction) {
    # Calculate phylogenetic independent contrasts correlation
    correlation_matrix <- phylo_correlation(traits_matrix, phylo_tree, method = method)
  } else {
    # Standard correlation analysis
    correlation_matrix <- Hmisc::rcorr(as.matrix(traits_matrix), type = method)
  }
  # Extract correlation and p-value matrices
  r <- correlation_matrix$r
  p <- correlation_matrix$P
  p <- stats::p.adjust(p, method = 'fdr')
  # Create correlation plot
  p1 <- corrplot::corrplot(r, type = "lower", order = "hclust", diag = FALSE,
                           tl.col = "black", tl.srt = 45, method = "circle",
                           insig = "pch", pch = 4, pch.col = "red", pch.cex = 1.5,
                           p.mat = (abs(r) < rThres) | (p >= pThres))
  return(p1)
}
