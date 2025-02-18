
library(igraph)
library(Hmisc)
library(corrplot)

#' Calculate and Visualize Plant Trait Correlation Network
#'
#' This function calculates correlation coefficients for given plant traits and generates a correlation network plot.
#'
#' @param traits_matrix A numeric matrix where each column represents a plant trait and each row represents a sample.
#' @param rThres Numeric, threshold for correlation coefficient, default is 0.2. Only correlations with absolute values above this threshold will be displayed in the plot.
#' @param pThres Numeric, threshold for p-value, default is 0.05. Only correlations with p-values below this threshold will be displayed in the plot.
#' @param method Character, specifies the correlation method to use: "pearson" (default) or "spearman".
#'
#' @return Returns a correlation network plot object.
#'
#' @details
#' The function first calculates Pearson correlation coefficients between traits, then adjusts p-values using the FDR method.
#' Finally, it plots the correlation network using the corrplot package. The plot displays only correlations that meet both the correlation coefficient and p-value thresholds.
#'
#' @references
#' 1. He, N., Li, Y., Liu, C., et al. (2020). Plant trait networks: improved resolution of the dimensionality of adaptation. Trends in Ecology & Evolution, 35(10), 908-918. https://doi.org/10.1016/j.tree.2020.06.003
#' 2. Li, Y., Liu, C., Sack, L., Xu, L., Li, M., Zhang, J., & He, N. (2022). Leaf trait network architecture shifts with species‐richness and climate across forests at continental scale. Ecology Letters, 25(6), 1442-1457. https://doi.org/10.1111/ele.14009
#'
#' @importFrom Hmisc rcorr
#' @importFrom stats p.adjust
#' @importFrom corrplot corrplot
#'
#' @examples
#' data(PFF)
#' PFF_traits <- PFF[, c("Height", "Leaf_area","LDMC","SLA","SRL","SeedMass","FltDate",
#'                       "FltDur","Leaf_Cmass","Leaf_Nmass","Leaf_CN","Leaf_Pmass",
#'                       "Leaf_NP","Leaf_CP","Root_Cmass","Root_Nmass","Root_CN")]
#' PFF_traits <- na.omit(PFF_traits)
#' head(PFF_traits)
#' TN_corr(traits_matrix = PFF_traits, rThres = 0.3, pThres = 0.01,method = "pearson")
#'
#' @export
TN_corr <- function(traits_matrix, rThres = 0.2, pThres = 0.05,method = "pearson") {
  correlation_matrix <- Hmisc::rcorr(as.matrix(traits_matrix), type = method)
  r <- correlation_matrix$r
  p <- correlation_matrix$P
  p <- stats::p.adjust(p, method = 'fdr')
  p1 <- corrplot::corrplot(r, type = "lower", order = "hclust", diag = FALSE,
                           tl.col = "black", tl.srt = 45, method = "circle",
                           insig = "pch", pch = 4, pch.col = "red", pch.cex = 1.5,
                           p.mat = (abs(r) < rThres) | (p >= pThres))
  return(p1)
}


#' Generate Plant Trait Network
#'
#' This function creates a network graph from a plant trait correlation matrix, applying thresholds for correlation strength and significance.
#'
#' @param traits_matrix A numeric matrix where each column represents a plant trait and each row represents a sample.
#' @param rThres Numeric, threshold for correlation coefficient, default is 0.2. Correlations with absolute values below this threshold are set to zero.
#' @param pThres Numeric, threshold for p-value, default is 0.05. Only correlations with p-values below this threshold are included in the network.
#' @param method Character, specifies the correlation method to use: "pearson" (default) or "spearman".
#'
#' @return Returns an igraph object representing the trait network.
#'
#' @details
#' The function performs the following steps:
#' 1. Calculates Pearson correlation coefficients and p-values for the trait matrix.
#' 2. Applies correlation coefficient and p-value thresholds to filter relationships.
#' 3. Constructs a weighted undirected graph from the filtered correlation matrix.
#' 4. Removes self-loops and isolated nodes from the graph.
#' 5. Adds correlation coefficients as edge attributes.
#'
#' @references
#' 1. He, N., Li, Y., Liu, C., et al. (2020). Plant trait networks: improved resolution of the dimensionality of adaptation. Trends in Ecology & Evolution, 35(10), 908-918. https://doi.org/10.1016/j.tree.2020.06.003
#' 2. Li, Y., Liu, C., Sack, L., Xu, L., Li, M., Zhang, J., & He, N. (2022). Leaf trait network architecture shifts with species‐richness and climate across forests at continental scale. Ecology Letters, 25(6), 1442-1457. https://doi.org/10.1111/ele.14009
#'
#' @importFrom Hmisc rcorr
#' @importFrom stats p.adjust
#' @importFrom igraph graph_from_adjacency_matrix simplify delete_vertices degree E
#'
#' @examples
#' data(PFF)
#' PFF_traits <- PFF[, c("Height", "Leaf_area","LDMC","SLA","SRL","SeedMass","FltDate",
#'                       "FltDur","Leaf_Cmass","Leaf_Nmass","Leaf_CN","Leaf_Pmass",
#'                       "Leaf_NP","Leaf_CP","Root_Cmass","Root_Nmass","Root_CN")]
#' PFF_traits <- na.omit(PFF_traits)
#' head(PFF_traits)
#' Tn_result <- TN(traits_matrix = PFF_traits, rThres = 0.2, pThres = 0.05, method = "pearson")
#' Tn_result
#'
#' @export
TN <- function(traits_matrix, rThres = 0.2, pThres = 0.05, method = "pearson") {
  # Validate the method parameter
  method <- match.arg(method, choices = c("pearson", "spearman"))
  # Calculate the correlation matrix based on the chosen methodology
  correlation_matrix <- Hmisc::rcorr(as.matrix(traits_matrix), type = method)

  # Threshold screening to exclude relationships with Pearson correlation coefficients lower than rThres
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

  # Transform the adjacency matrix into the adjacency list of an igraph network, and construct a weighted undirected network
  g <- igraph::graph_from_adjacency_matrix(z, weighted = TRUE, mode = 'undirected')

  # Remove self-loops
  g <- igraph::simplify(g)

  # Delete isolated nodes (nodes with a degree of 0)
  g <- igraph::delete_vertices(g, names(igraph::degree(g)[igraph::degree(g) == 0]))

  # Set edge weights to absolute correlation values and store original correlation values
  igraph::E(g)$correlation <- igraph::E(g)$weight
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)

  return(g)
}


#' Calculate Node and Global Metrics for Trait Networks
#'
#' This function computes various node and global metrics for a trait network graph.
#'
#' @param graph An igraph object representing the trait network, typically generated by the `TN` function.
#'
#' @return A list containing two data frames:
#'   \item{node}{A data frame with node-level metrics including degree, closeness, betweenness, and local clustering coefficient.}
#'   \item{global}{A data frame with global metrics including edge density, diameter, average path length, average clustering coefficient, and modularity.}
#'
#' @references
#' 1. He, N., Li, Y., Liu, C., et al. (2020). Plant trait networks: improved resolution of the dimensionality of adaptation. Trends in Ecology & Evolution, 35(10), 908-918. https://doi.org/10.1016/j.tree.2020.06.003
#' 2. Li, Y., Liu, C., Sack, L., Xu, L., Li, M., Zhang, J., & He, N. (2022). Leaf trait network architecture shifts with species‐richness and climate across forests at continental scale. Ecology Letters, 25(6), 1442-1457. https://doi.org/10.1111/ele.14009
#'
#' @examples
#' data(PFF)
#' PFF_traits <- PFF[, c("Height", "Leaf_area","LDMC","SLA","SRL","SeedMass","FltDate",
#'                       "FltDur","Leaf_Cmass","Leaf_Nmass","Leaf_CN","Leaf_Pmass",
#'                       "Leaf_NP","Leaf_CP","Root_Cmass","Root_Nmass","Root_CN")]
#' PFF_traits <- na.omit(PFF_traits)
#' head(PFF_traits)
#' Tn_result <- TN(traits_matrix = PFF_traits, rThres = 0.2, pThres = 0.05)
#' TN_metrics(Tn_result)
#'
#' @importFrom igraph degree closeness betweenness transitivity edge_density diameter mean_distance cluster_fast_greedy modularity membership
#' @export
TN_metrics <- function(graph) {
  # Node-level metrics
  degree <- igraph::degree(graph)
  closeness <- igraph::closeness(graph)
  betweenness <- igraph::betweenness(graph)
  local_clustering <- igraph::transitivity(graph, type = "local")
  # Set the clustering coefficient of nodes with degree 0 or 1 to 0
  local_clustering[is.nan(local_clustering) | degree <= 1] <- 0

  # Creating node-level metrics dataframe
  node_metrics <- data.frame(
    degree = degree,
    closeness = closeness,
    betweenness = betweenness,
    clustering_coefficient = local_clustering
  )

  # Global metrics
  edge_density <- igraph::edge_density(graph)
  diameter <- igraph::diameter(graph)
  avg_path_length <- igraph::mean_distance(graph)
  avg_clustering <- igraph::transitivity(graph, type = "global")
  fc <- igraph::cluster_fast_greedy(graph)
  modularity <- igraph::modularity(graph, igraph::membership(fc))

  # Creating global metrics dataframe
  global_metrics <- data.frame(
    edge_density = edge_density,
    diameter = diameter,
    avg_path_length = avg_path_length,
    avg_clustering_coefficient = avg_clustering,
    modularity = modularity
  )

  # Return a list containing both node and global metrics
  return(list(
    node = node_metrics,
    global = global_metrics
  ))
}


#' Plot Trait Network Graph
#'
#' This function visualizes the trait network graph generated by the `TN` function.
#'
#' @param graph An igraph object representing the trait network.
#' @param style A numeric value that determines the plotting style (default is 1).
#' @param vertex.size Numeric value for the size of vertices in the plot (default is 20).
#' @param vertex.label.cex Numeric value for the scaling factor of vertex labels (default is 0.6).
#'
#' @return
#' An object of class `igraph`. This function generates a visualization of the trait network graph.
#' When style = 1, it displays a community structure plot.
#' When style = 2, it displays a circular layout plot where vertex colors represent community membership,
#' edge thickness represents correlation strength, and edge color represents the sign of the correlation (black for positive, red for negative).
#'
#' @details
#' The function uses the `cluster_fast_greedy` algorithm to identify communities
#' in the graph and assigns community membership to vertices. It offers two
#' plotting styles:
#' - Style 1: Plots the community structure.
#' - Style 2: Plots the graph in a circular layout with vertex colors representing communities.
#' The vertex size and label size can be customized using vertex.size and vertex.label.cex parameters respectively.
#'
#' @references
#' 1. He, N., Li, Y., Liu, C., et al. (2020). Plant trait networks: improved resolution of the dimensionality of adaptation. Trends in Ecology & Evolution, 35(10), 908-918. https://doi.org/10.1016/j.tree.2020.06.003
#' 2. Li, Y., Liu, C., Sack, L., Xu, L., Li, M., Zhang, J., & He, N. (2022). Leaf trait network architecture shifts with species‐richness and climate across forests at continental scale. Ecology Letters, 25(6), 1442-1457. https://doi.org/10.1111/ele.14009
#'
#' @examples
#' data(PFF)
#' PFF_traits <- PFF[, c("Height", "Leaf_area","LDMC","SLA","SRL","SeedMass","FltDate",
#'                       "FltDur","Leaf_Cmass","Leaf_Nmass","Leaf_CN","Leaf_Pmass",
#'                       "Leaf_NP","Leaf_CP","Root_Cmass","Root_Nmass","Root_CN")]
#' PFF_traits <- na.omit(PFF_traits)
#' head(PFF_traits)
#' Tn_result <- TN(traits_matrix = PFF_traits, rThres = 0.2, pThres = 0.05)
#' TN_plot(Tn_result, style = 1, vertex.size = 20, vertex.label.cex = 0.6)
#' TN_plot(Tn_result, style = 2, vertex.size = 20, vertex.label.cex = 0.6)
#'
#' @importFrom igraph cluster_fast_greedy membership layout_in_circle V
#' @export
TN_plot <- function(graph, style = 1,vertex.size = 20,vertex.label.cex = 0.6) {
  comm <- igraph::cluster_fast_greedy(graph)
  igraph::V(graph)$community <- igraph::membership(comm)
  layout <- igraph::layout_in_circle(graph, order = order(igraph::V(graph)$community))

  if (style == 1) {
    plot(comm, graph,vertex.size = vertex.size,vertex.label.cex = vertex.label.cex,vertex.label.color = "black")
  } else if (style == 2) {
    # Get the correlation coefficients of the edges
    edge_correlation <- igraph::E(graph)$correlation
    # Set the line width (according to the absolute value of the correlation coefficient)
    edge_width <- abs(edge_correlation) * 5  # The multiplier can be adjusted to change the line thickness range
    # Set the color of the edge (blue for negative correlation, green for positive correlation)
    edge_color <- ifelse(edge_correlation > 0, "green", "blue")
    # Drawing
    plot(graph, layout = layout,
         vertex.color = igraph::V(graph)$community,
         vertex.label.cex = vertex.label.cex,
         vertex.label.dist = 0,
         vertex.size = vertex.size,
         vertex.label.color = "black",
         vertex.label.font = 2,
         edge.width = edge_width,    # Add line width
         edge.color = edge_color)    # Add line colors
  } else {
    stop("Invalid style. Please choose 1 or 2.")
  }
}
