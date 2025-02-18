
library(vegan)
library(ggplot2)
library(ggrepel)


#' Perform nested principal component analysis (PCA of PCAs) for ecological niche periodicity
#'
#' This function conducts a ‘PCA of PCAs’ to analyze
#' ecological niche periodicity based on multiple trait dimensions.
#'
#' @param data A data frame containing species trait data.
#' @param dimension A list of character vectors, each representing a trait dimension
#'        with corresponding column names from the data.
#'
#' @return A list containing three elements:
#'   \item{PCA_first}{A data frame summarizing the first-level PCA results for each dimension}
#'   \item{PCA_second}{A matrix of species scores from the second-level PCA}
#'   \item{result}{The complete rda object from the second-level PCA}
#'
#' @details
#' The function performs the following steps:
#' 1. Checks and cleans input data, removing rows with NA values
#' 2. Conducts first-level PCA for each trait dimension
#' 3. Extracts PC1 and PC2 scores from each first-level PCA
#' 4. Combines all PC scores and performs a second-level PCA
#' 5. Returns summary statistics and results from both levels of analysis
#'
#' @references
#' 1. Winemiller, K. O., Fitzgerald, D. B., Bower, L. M., & Pianka, E. R. (2015). Functional traits, convergent evolution, and periodic tables of niches. Ecology letters, 18(8), 737-751. https://doi.org/10.1111/ele.12462
#' 2. Yu, R., Huang, J., Xu, Y., Ding, Y., & Zang, R. (2020). Plant functional niches in forests across four climatic zones: Exploring the periodic table of niches based on plant functional traits. Frontiers in Plant Science, 11, 841. https://doi.org/10.3389/fpls.2020.00841
#'
#' @examples
#' data(PFF)
#' PFF[,3:20] <- log(PFF[,3:20])
#' traits_dimension <- list(
#'   grow = c("SLA","Leaf_area","LDMC","SRL","Leaf_Nmass","Leaf_Pmass","Root_Nmass"),
#'   survive = c("Height","Leaf_Cmass","Root_Cmass","Leaf_CN","Leaf_NP","Leaf_CP","Root_CN"),
#'   reproductive = c("SeedMass","FltDate","FltDur")
#' )
#' result <- NPT(data = PFF, dimension = traits_dimension)
#' result
#'
#' @importFrom vegan rda scores
#' @importFrom stats na.omit
#' @export
NPT <- function(data, dimension) {
  # Check input data and remove rows with NA values
  original_rows <- nrow(data)
  data <- na.omit(data[, unlist(dimension)])
  removed_rows <- original_rows - nrow(data)

  if (removed_rows > 0) {
    message(paste("Removed", removed_rows, "rows containing NA values."))
  }

  pca_first <- function(data) {
    pca_result <- vegan::rda(data, scale = FALSE)
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
    result <- list(
      summary = data.frame(
        pc1_percent = pc1_percent,
        pc1_major_eigenvector = pc1_major_eigenvector,
        pc2_percent = pc2_percent,
        pc2_major_eigenvector = pc2_major_eigenvector
      ),
      PC = data.frame(
        PC1 = pc1_scores,
        PC2 = pc2_scores
      )
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
  f <- summary(result)
  sec_scores <- vegan::scores(f, display = "species", scaling = 0)[, 1:2]
  r <- list(PCA_first = S, PCA_second = sec_scores, result = result)
  return(r)
}


#' Plot results from nested principal component analysis
#'
#' This function creates a visualization of the results from the nested principal
#' component analysis (NPT function) for ecological niche periodicity.
#'
#' @param pca_obj An rda object returned by the second-level PCA in the NPT function.
#' @param group A vector or factor specifying the grouping of samples.
#'
#' @return A ggplot object representing a biplot of species scores and sample points.
#'
#' @details
#' The function creates a biplot that includes:
#' 1. Sample points colored by group
#' 2. Species scores represented as arrows
#' 3. Species labels positioned using ggrepel to avoid overlapping
#'
#' The plot also includes dashed lines at x=0 and y=0, and displays the percentage
#' of variance explained by each principal component on the axes.
#'
#' @references
#' 1. Winemiller, K. O., Fitzgerald, D. B., Bower, L. M., & Pianka, E. R. (2015). Functional traits, convergent evolution, and periodic tables of niches. Ecology letters, 18(8), 737-751. https://doi.org/10.1111/ele.12462
#' 2. Yu, R., Huang, J., Xu, Y., Ding, Y., & Zang, R. (2020). Plant functional niches in forests across four climatic zones: Exploring the periodic table of niches based on plant functional traits. Frontiers in Plant Science, 11, 841. https://doi.org/10.3389/fpls.2020.00841
#'
#' @examples
#' data(PFF)
#' PFF[,3:20] <- log(PFF[,3:20])
#' PFF <- na.omit(PFF)
#' traits_dimension <- list(
#'   grow = c("SLA","Leaf_area","LDMC","SRL","Leaf_Nmass","Leaf_Pmass","Root_Nmass"),
#'   survive = c("Height","Leaf_Cmass","Root_Cmass","Leaf_CN","Leaf_NP","Leaf_CP","Root_CN"),
#'   reproductive = c("SeedMass","FltDate","FltDur")
#' )
#' npt_result <- NPT(data = PFF, dimension = traits_dimension)
#' dev.new() # A window that is too small will interfere with the drawing.
#'           # Optionally, you can set the drawing window to pop up automatically.
#' NPT_plot(npt_result$result)
#' NPT_plot(npt_result$result, PFF$family)
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom vegan scores
#' @export
NPT_plot <- function(pca_obj, group = NULL) {
  # Extract site scores and add group information
  site_scores <- as.data.frame(vegan::scores(pca_obj, display = "sites"))
  site_scores$sample <- rownames(site_scores)
  site_scores$group <- group

  # Extract species scores
  species_scores <- as.data.frame(vegan::scores(pca_obj, display = "species"))
  species_scores$species <- rownames(species_scores)

  # Calculate variance explained by each PC
  var_explained <- pca_obj$CA$eig / sum(pca_obj$CA$eig) * 100

  # Create the plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = site_scores,
                        ggplot2::aes(x = .data$PC1, y = .data$PC2, fill = group),
                        shape = 21, size = 2, alpha = 0.5) +
    ggplot2::geom_segment(data = species_scores,
                          ggplot2::aes(x = 0, y = 0, xend = .data$PC1, yend = .data$PC2),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
                          color = "black") +
    ggrepel::geom_text_repel(data = species_scores,
                             ggplot2::aes(x = .data$PC1, y = .data$PC2, label = .data$species),
                             color = "black", size = 4, max.overlaps = Inf) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.title = ggplot2::element_text(face = "bold", size = 9),
      legend.text = ggplot2::element_text(size = 8),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      legend.key.size = unit(0.8, "lines"),
      legend.position = "right",
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text = ggplot2::element_text(size = 10, color = "black")
    ) +
    ggplot2::labs(
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      fill = "Group"
    )

  if (is.null(group)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}




