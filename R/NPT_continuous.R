#' Continuous Niche Classification Based on Periodic Table of Niches
#'
#' This function implements a continuous niche classification scheme based on the
#' periodic table of niches concept. It performs a hierarchical Principal Component
#' Analysis (PCA) approach where separate PCAs are conducted on different niche
#' dimensions, followed by a second-level PCA to integrate results across dimensions.
#'
#' @param data A data frame containing species and their functional trait measurements.
#'   Each row represents a species and columns contain trait values. The data should
#'   include all traits specified in the dimension parameter.
#' @param dimension A named list where each element represents a niche dimension
#'   (e.g., "grow", "survive", "reproductive") and contains a character vector of
#'   column names corresponding to traits associated with that dimension. Each
#'   dimension should contain multiple functionally related traits.
#'
#' @return A list containing three elements:
#'   \item{PCA_first}{A data frame summarizing the first-level PCA results for each
#'     dimension, including variance explained by PC1 and PC2 (as percentages), and
#'     the traits with highest absolute loadings on each principal component axis}
#'   \item{PCA_second}{A matrix containing species scores from the second-level PCA
#'     that integrates all niche dimensions into a unified ordination space}
#'   \item{result}{The complete second-level PCA result object from vegan::rda()
#'     containing detailed ordination results for further analysis}
#'
#' @details
#' The function implements a two-stage hierarchical PCA approach based on the
#' methodology described in Winemiller et al. (2015) for creating continuous niche
#' classification schemes. This approach addresses the challenge that analysis of
#' data sets containing many functionally unrelated measures may fail to detect
#' patterns of covariation that determine species' ecological responses to and
#' effects on their environments.
#'
#' \strong{Stage 1: Dimensional PCA Analysis}
#'
#' Separate Principal Component Analysis is performed on trait data for each niche
#' dimension using the internal \code{pca_first} function. This dimensional approach
#' ensures that all niche dimensions have an equal opportunity to influence the
#' composite niche scheme and species ordinations. For each dimension, the function:
#' \itemize{
#'   \item Performs PCA using \code{vegan::rda()}
#'   \item Calculates variance explained by the first two principal components
#'   \item Identifies traits with highest absolute loadings on PC1 and PC2
#'   \item Extracts species scores on both principal components
#' }
#'
#' \strong{Stage 2: Integration PCA}
#'
#' The species scores from the first two principal components of each dimensional
#' PCA are combined into a new data matrix (with columns named as "pc1.dimension"
#' and "pc2.dimension"). A second PCA is then performed on this matrix to create
#' a two-dimensional continuum integrating patterns (strategies) within each of
#' the niche dimensions. This creates a continuous ordination of species within
#' niche space that can be used for comparative ecological analyses.
#'
#' \strong{Methodological Advantages}
#'
#' This hierarchical PCA method prevents domination by any single type of trait
#' or dimension. The approach allows all niche dimensions to have equal influence
#' on the composite niche scheme, with gradients dominated by those dimensional
#' components having greatest influence on community structure patterns.
#'
#' @note
#' \strong{Important Considerations:}
#' \itemize{
#'   \item Missing values (NA) are automatically removed using \code{na.omit()}
#'   \item A message displays the number of rows removed due to NA values
#'   \item Each niche dimension should contain multiple functionally related traits
#'   \item Users should ensure traits are appropriately scaled/transformed before analysis
#'   \item The function uses scaling = FALSE in PCA, assuming pre-standardized data
#'   \item Column names in the data must exactly match trait names in dimensions list
#' }
#'
#' @references
#' 1. Winemiller, K. O., Fitzgerald, D. B., Bower, L. M., & Pianka, E. R. (2015). Functional traits,
#' convergent evolution, and periodic tables of niches. Ecology letters, 18(8), 737-751.
#' 2. Yu, R., Huang, J., Xu, Y., Ding, Y., & Zang, R. (2020). Plant functional niches in forests across four climatic zones:
#' Exploring the periodic table of niches based on plant functional traits. Frontiers in Plant Science, 11, 841.
#'
#' @examples
#' data(PFF)
#' PFF[,4:21] <- log(PFF[,4:21])
#' traits_dimension <- list(
#'   grow = c("SLA","SRL","Leaf_Nmass","Root_Nmass"),
#'   survive = c("Height","Leaf_CN","Root_CN"),
#'   reproductive = c("SeedMass","FltDate","FltDur")
#' )
#' result <- NPT_continuous(data = PFF, dimension = traits_dimension)
#' result
#'
#' @importFrom vegan rda scores
#' @importFrom stats na.omit
#' @export
NPT_continuous <- function(data, dimension) {
  # Check input data and remove rows with NA values
  original_rows <- nrow(data)
  data <- na.omit(data[, unlist(dimension)])
  removed_rows <- original_rows - nrow(data)
  if (removed_rows > 0) {
    message(paste("Removed", removed_rows, "rows containing NA values."))
  }
  pca_first <- function(data) {
    pca_result <- vegan::rda(data)
    var_explained <- summary(pca_result)$cont$importance[2, ]
    pc1_percent <- var_explained[1] * 100
    pc2_percent <- var_explained[2] * 100
    loadings <- vegan::scores(pca_result, display = "species", scaling = 0)
    pc1_main <- loadings[, 1]
    pc1_max <- sort(abs(pc1_main), decreasing = TRUE)[1]
    pc1_major_eigenvector <- names(pc1_max)
    pc2_main <- loadings[, 2]
    pc2_max <- sort(abs(pc2_main), decreasing = TRUE)[1]
    pc2_major_eigenvector <- names(pc2_max)
    species_scores <- vegan::scores(pca_result, display = "sites", scaling = 0)
    pc1_scores <- species_scores[, 1]
    pc2_scores <- species_scores[, 2]
    result <- list(summary = data.frame(pc1_percent = pc1_percent,
                                        pc1_major_eigenvector = pc1_major_eigenvector,
                                        pc2_percent = pc2_percent,
                                        pc2_major_eigenvector = pc2_major_eigenvector),
                   PC = data.frame(PC1 = pc1_scores, PC2 = pc2_scores)
                   )
    return(result)
  }
  n <- length(dimension)
  S <- NULL
  P <- NULL
  for (i in 1:n) {
    dt <- data[, dimension[[i]]]
    a <- pca_first(data = dt)
    s <- a$summary
    p <- t(a$PC)
    S <- rbind(S, s)
    P <- rbind(P, p)
  }
  rownames(S) <- names(dimension)
  P <- t(P)
  colnames(P) <- paste0(rep(c("pc1", "pc2"), n), ".", rep(names(dimension), each = 2))
  result <- vegan::rda(P)
  sec_scores <- vegan::scores(result, display = "species", scaling = 0)[, 1:2]
  r <- list(PCA_first = S, PCA_second = sec_scores, result = result)
  return(r)
}


