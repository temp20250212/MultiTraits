#' Generate Plant Trait Network
#'
#' This function creates a network graph from a plant trait correlation matrix, applying thresholds for correlation strength and significance. It supports both standard correlations and phylogenetic independent contrasts.
#'
#' @param traits_matrix A numeric matrix or data frame where each column represents a plant trait and each row represents a sample/species. Row names should contain species names when using phylogenetic correction.
#' @param rThres Numeric, threshold for correlation coefficient, default is 0.2. Correlations with absolute values below this threshold are set to zero.
#' @param pThres Numeric, threshold for p-value, default is 0.05. Only correlations with p-values below this threshold are included in the network.
#' @param method Character, specifies the correlation method to use: "pearson" (default) or "spearman".
#' @param phylo_correction Logical, whether to apply phylogenetic correction using phylogenetic independent contrasts, default is FALSE.
#' @param phylo_tree A phylo object from the ape package containing the phylogenetic tree. Required when phylo_correction = TRUE. Species names in the tree must match row names in traits_matrix.
#'
#' @return Returns an igraph object representing the trait network with the following attributes:
#' \itemize{
#'   \item Edge attribute 'correlation': original correlation values (positive or negative)
#'   \item Edge attribute 'weight': absolute correlation values used for network analysis
#'   \item Vertices represent traits that pass the correlation and significance thresholds
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. Validates input parameters and phylogenetic tree compatibility (if applicable).
#' 2. Calculates correlation coefficients and p-values using either standard correlation or phylogenetic independent contrasts.
#' 3. Applies correlation coefficient and p-value thresholds to filter relationships.
#' 4. Adjusts p-values using False Discovery Rate (FDR) correction.
#' 5. Constructs an unweighted undirected graph from the filtered correlation matrix.
#' 6. Removes self-loops and isolated nodes from the graph.
#' 7. Adds correlation coefficients as edge attributes.
#'
#' When phylo_correction = TRUE, the function:
#' - Matches species names between traits_matrix and phylo_tree
#' - Calculates phylogenetic independent contrasts for each trait
#' - Computes correlations between contrasts to control for phylogenetic relatedness
#'
#' @references
#' 1. He, N., Li, Y., Liu, C., et al. (2020). Plant trait networks: improved resolution of the dimensionality of adaptation.
#' Trends in Ecology & Evolution, 35(10), 908-918.
#' 2. Li, Y., Liu, C., Sack, L., Xu, L., Li, M., Zhang, J., & He, N. (2022). Leaf trait network architecture shifts with
#' species‐richness and climate across forests at continental scale. Ecology Letters, 25(6), 1442-1457.
#'
#'
#' @importFrom Hmisc rcorr
#' @importFrom stats p.adjust
#' @importFrom igraph graph_from_adjacency_matrix simplify delete_vertices degree E graph_from_data_frame as_data_frame
#' @importFrom ape keep.tip
#'
#' @examples
#' # Example 1: Standard trait network analysis
#' data(PFF)
#' rownames(PFF) <- PFF$species
#' PFF_traits <- PFF[, c("Height", "Leaf_area","LDMC","SLA","SRL","SeedMass","FltDate",
#'                       "FltDur","Leaf_Cmass","Leaf_Nmass","Leaf_CN","Leaf_Pmass",
#'                       "Leaf_NP","Leaf_CP","Root_Cmass","Root_Nmass","Root_CN")]
#' PFF_traits <- na.omit(PFF_traits)
#' head(PFF_traits)
#'
#' ptn_result <- PTN(traits_matrix = PFF_traits, rThres = 0.2, pThres = 0.05, method = "pearson")
#' ptn_result
#'
#' # Example 2: Phylogenetically corrected trait network analysis
#' data(PFF_tree)
#'
#' # Trait network with phylogenetic correction
#' ptn_phylo_result <- PTN(traits_matrix = PFF_traits,
#'                       rThres = 0.2,
#'                       pThres = 0.05,
#'                       method = "pearson",
#'                       phylo_correction = TRUE,
#'                       phylo_tree = PFF_tree)
#' ptn_phylo_result
#'
#'
#' @export
PTN <- function(traits_matrix, rThres = 0.2, pThres = 0.05, method = "pearson",
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
  # Threshold screening to exclude relationships with correlation coefficients lower than rThres
  r <- correlation_matrix$r
  r[abs(r) < rThres] <- 0
  # Select correlation coefficients with significance p-values less than pThres
  p <- correlation_matrix$P
  p <- stats::p.adjust(p, method = 'fdr')
  p[p >= pThres] <- -1
  p[p < pThres & p >= 0] <- 1
  p[p == -1] <- 0
  # Data retained based on r-values and p-values screened above
  z <- r * p
  diag(z) <- 0
  g <- igraph::graph_from_adjacency_matrix(z, weighted = TRUE, mode = 'undirected')
  # Remove self-loops
  g <- igraph::simplify(g)
  # Delete isolated nodes (nodes with a degree of 0)
  g <- igraph::delete_vertices(g, names(igraph::degree(g)[igraph::degree(g) == 0]))
  # Set edge weights to absolute correlation values and store original correlation values
  igraph::E(g)$correlation <- igraph::E(g)$weight
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  gdf <- igraph::as_data_frame(g, what = "edges")
  gdf <- gdf[,c("from","to","correlation")]
  nodes_name <- unique(c(gdf$from, gdf$to))
  nodes <- data.frame(id = nodes_name, num = 1)
  gg <- igraph::graph_from_data_frame(gdf, vertices = nodes, directed = FALSE)
  return(gg)
}
