#' Discrete Niche Classification Based on Periodic Table of Niches
#'
#' This function implements a discrete niche classification scheme based on the
#' periodic table of niches concept. It performs Principal Component Analysis (PCA)
#' on functional traits grouped by niche dimensions, followed by clustering to
#' create hierarchical niche classifications.
#'
#' @param data A data frame containing species and their functional trait measurements.
#'   Each row represents a species and columns contain trait values.
#' @param dimension A named list where each element represents a niche dimension
#'   (e.g., "Growth", "Survival", "Reproduction") and contains a character vector of
#'   column names corresponding to traits associated with that dimension.
#' @param clustering_method Character string specifying the clustering method to use.
#'   Options are:
#'   \itemize{
#'     \item \code{"CART"} - Classification and Regression Trees (default)
#'     \item \code{"kmeans"} - K-means clustering with automatic k selection
#'   }
#' @param k_max Integer specifying the maximum number of clusters allowed for k-means
#'   clustering. Default is \code{6}. This parameter is only used when
#'   \code{clustering_method = "kmeans"}. The optimal k value will be selected using
#'   the elbow method, constrained by this maximum value.
#'
#' @return A list containing two elements:
#'   \item{niche_classification}{A data frame with species names, cluster assignments
#'     for each dimension, and comprehensive niche codes}
#'   \item{summary}{A summary data frame showing unique niche codes, the number of
#'     species in each niche, and lists of species names}
#'
#' @details
#' The function implements the methodology described in Winemiller et al. (2015)
#' for creating discrete niche classification schemes. The approach follows three
#' main steps:
#'
#' \strong{Step 1: PCA Analysis by Dimension}
#'
#' Separate Principal Component Analysis is performed on trait data for each niche
#' dimension. This dimensional approach prevents functionally unrelated traits from
#' masking important ecological patterns. The first two principal components are
#' retained for each dimension.
#'
#' \strong{Step 2: Clustering Methods}
#'
#' Two clustering approaches are available:
#' \itemize{
#'   \item \strong{CART}: Uses the Euclidean distance from the origin in PCA space
#'     as the response variable and original trait values as predictors to create
#'     regression trees. The tree is pruned using cross-validation.
#'   \item \strong{k-means}: Performs clustering on the two-dimensional PCA space
#'     with automatic optimal k selection using the elbow method. The optimal k
#'     is constrained by the \code{k_max} parameter (default maximum = 6).
#' }
#'
#' \strong{Step 3: Hierarchical Niche Classification}
#'
#' The function combines clustering results from all niche dimensions to create
#' a comprehensive niche classification scheme. Each species receives a niche code
#' representing its cluster membership across all dimensions (e.g., "1,2,1" for
#' cluster 1 in dimension 1, cluster 2 in dimension 2, and cluster 1 in dimension 3).
#'
#' The function also calculates niche occupancy statistics, comparing the number
#' of realized niches to the total number of potential niche combinations.
#'
#' @note
#' \itemize{
#'   \item Missing values should be handled prior to using this function.
#'   \item The function prints diagnostic information during execution, including
#'     variance explained by PCs and clustering results.
#'   \item For k-means clustering, setting a lower \code{k_max} value will force
#'     simpler clustering solutions, while higher values allow for more complex
#'     niche subdivisions.
#'   \item \strong{Randomness warning:} The \code{kmeans} method involves
#'         random initialization of cluster centers, so results may vary
#'         between runs. For reproducibility, set a random seed using
#'         \code{set.seed()} before running this function. CART results may
#'         also differ slightly if predictor splitting involves tie-breaking.
#' }
#'
#'
#' @importFrom vegan rda scores
#' @importFrom rpart rpart prune
#' @importFrom dplyr group_by summarise n
#' @importFrom stats kmeans predict
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @references
#' 1. Winemiller, K. O., Fitzgerald, D. B., Bower, L. M., & Pianka, E. R. (2015).
#'    Functional traits, convergent evolution, and periodic tables of niches.
#'    Ecology letters, 18(8), 737-751.
#' 2. Pianka, E. R., Vitt, L. J., Pelegrin, N., Fitzgerald, D. B., & Winemiller, K. O. (2017).
#'    Toward a periodic table of niches, or exploring the lizard niche hypervolume.
#'    The American Naturalist, 190(5), 601-616.
#'
#' @examples
#' \dontrun{
#' # Load and prepare data
#' data(PFF)
#' rownames(PFF) <- PFF$species
#' PFF_traits <- PFF[, c("SLA", "SRL", "Leaf_Nmass", "Root_Nmass","Height",
#'                       "Leaf_CN", "Root_CN","SeedMass", "FltDate", "FltDur")]
#' # Perform log transformation of data and remove missing values
#' PFF_traits <- log(na.omit(PFF_traits))
#' head(PFF_traits)
#' # Define trait dimensions
#' dimension <- list(Grow = c("SLA", "SRL", "Leaf_Nmass", "Root_Nmass"),
#'                   Survive = c("Height", "Leaf_CN", "Root_CN"),
#'                   Reproductive = c("SeedMass", "FltDate", "FltDur"))
#'
#' set.seed(123)
#' discrete_result <- NPT_discrete(data = PFF_traits, dimension = dimension)
#' head(discrete_result$niche_classification)
#' }
#'
#' @export
NPT_discrete <- function(data,
                         dimension,
                         clustering_method = "CART",
                         k_max = 6) {
  # Data validation
  # Check if data is a dataframe and contains row names
  if (!is.data.frame(data)) {
    stop("Input data must be a dataframe.")
  }
  if (is.null(rownames(data))) {
    stop("Input data must contain row names as species identifiers.")
  }

  # Check if trait columns exist
  all_traits <- unlist(dimension)
  missing_traits <- setdiff(all_traits, names(data))
  if (length(missing_traits) > 0) {
    stop("The following trait columns are missing from data: ", paste(missing_traits, collapse = ", "))
  }

  # Step 1: PCA analysis for each dimension
  pca_results <- list()
  pc_scores <- data.frame(species = rownames(data)) # Use row names as species identifiers
  cat("=== PCA Analysis Results ===\n")
  for(dim_name in names(dimension)) {
    # Extract trait data for current dimension
    traits_data <- data[, dimension[[dim_name]]]
    # Perform PCA analysis using vegan package
    pca <- vegan::rda(traits_data)
    pca_results[[dim_name]] <- pca
    # Get sample scores
    site_scores <- vegan::scores(pca, display = "sites", choices = 1:2)
    # Save first two principal component scores
    pc_scores[[paste0(dim_name, "_PC1")]] <- site_scores[, 1]
    pc_scores[[paste0(dim_name, "_PC2")]] <- site_scores[, 2]
    # Calculate variance explained
    eigenvals <- pca$CA$eig
    total_var <- sum(eigenvals)
    pc1_var <- eigenvals[1] / total_var * 100
    pc2_var <- eigenvals[2] / total_var * 100
    cat("Dimension:", dim_name, "\n")
    cat("PC1 variance explained:", round(pc1_var, 2), "%\n")
    cat("PC2 variance explained:", round(pc2_var, 2), "%\n\n")
  }

  # Step 2: Clustering based on selected method
  dimension_clusters <- list()
  cat("=== ", clustering_method, " Clustering Results ===\n")
  for(dim_name in names(dimension)) {
    pc1_scores <- pc_scores[, paste0(dim_name, "_PC1")]
    pc2_scores <- pc_scores[, paste0(dim_name, "_PC2")]

    if(clustering_method == "CART") {
      # CART method
      predictor_data <- data[, dimension[[dim_name]]]
      pc_distance <- sqrt(pc1_scores^2 + pc2_scores^2)
      cart_distance <- rpart::rpart(pc_distance ~ .,
                                    data = predictor_data,
                                    method = "anova")
      cp_optimal <- cart_distance$cptable[which.min(cart_distance$cptable[,"xerror"]),"CP"]
      cart_distance_pruned <- rpart::prune(cart_distance, cp = cp_optimal)
      distance_predictions <- stats::predict(cart_distance_pruned, predictor_data)
      dimension_clusters[[dim_name]] <- as.numeric(as.factor(distance_predictions))
    } else if(clustering_method == "kmeans") {
      # k-means method
      pc_data <- cbind(pc1_scores, pc2_scores)
      # Use elbow method to determine optimal k value
      wss <- sapply(1:min(10, nrow(pc_data)-1), function(k) {
        if(k == 1) {
          return(sum(scale(pc_data, scale = FALSE)^2))
        }
        stats::kmeans(pc_data, centers = k, nstart = 20)$tot.withinss
      })
      # Select optimal k value using elbow detection
      if(length(wss) >= 3) {
        diffs <- diff(wss)
        k_optimal <- which.min(diff(diffs)) + 2
        k_optimal <- max(2, min(k_optimal, k_max = 6))
      } else {
        k_optimal <- 2
      }
      # Perform k-means clustering
      kmeans_result <- stats::kmeans(pc_data, centers = k_optimal, nstart = 20)
      dimension_clusters[[dim_name]] <- kmeans_result$cluster
    }
    cat("Dimension", dim_name, "cluster count:", length(unique(dimension_clusters[[dim_name]])), "\n")
  }

  # Step 3: Construct hierarchical discrete niche classification
  niche_classification <- data.frame(species = rownames(data))  # Use row names as species identifiers
  # Add clustering results for each dimension
  for(dim_name in names(dimension)) {
    niche_classification[[dim_name]] <- dimension_clusters[[dim_name]]
  }

  # Generate comprehensive niche category codes
  cluster_cols <- names(dimension)
  niche_codes <- apply(niche_classification[, cluster_cols], 1, function(x) {
    paste(x, collapse = ",")
  })
  niche_classification$niche_code <- niche_codes

  # Calculate occupied and potential niche counts
  total_potential_niches <- prod(sapply(dimension_clusters, function(x) length(unique(x))))
  occupied_niches <- length(unique(niche_codes))
  cat("\n=== Niche Classification Results ===\n")
  cat("Total potential niches:", total_potential_niches, "\n")
  cat("Actual occupied niches:", occupied_niches, "\n")
  cat("Niche occupancy rate:", round(occupied_niches/total_potential_niches*100, 2), "%\n")

  # Custom sorting function: sort niche_code numerically from small to large
  sort_niche_code <- function(niche_codes) {
    # Decompose niche_code into numeric vectors for sorting
    niche_matrix <- do.call(rbind, lapply(niche_codes, function(x) {
      as.numeric(strsplit(x, ",")[[1]])
    }))
    # Use do.call and order for multi-column sorting
    order_index <- do.call(order, as.data.frame(niche_matrix))
    return(order_index)
  }

  niche_summary <- niche_classification %>%
    dplyr::group_by(.data$niche_code) %>%
    dplyr::summarise(species_count = n(),
                     species_list = paste(.data$species, collapse = ", "),
                     .groups = 'drop')
  # Sort by niche_code numeric order
  sort_index <- sort_niche_code(niche_summary$niche_code)
  niche_summary <- niche_summary[sort_index, ]

  # Return results list
  result <- list(niche_classification = niche_classification, summary = niche_summary)
  return(result)
}
