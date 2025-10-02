#' Plot Continuous Niche Classification Results
#'
#' This function creates a biplot visualization of the continuous niche classification
#' results from the hierarchical Principal Component Analysis. It displays species
#' ordination in niche space with optional grouping and shows the contribution of
#' different niche dimensions as arrows.
#'
#' @param pca_obj The PCA result object from \code{NPT_continuous()$result}, which
#'   should be a vegan rda object containing the second-level PCA results that
#'   integrate all niche dimensions into a unified ordination space.
#' @param group Optional vector specifying group membership for each species/sample.
#'   If provided, points will be colored by group. If NULL (default), all points
#'   will have the same color. Length should match the number of rows in the
#'   original data.
#' @param default_fill Character string specifying the default fill color when no
#'   grouping is applied. Default is "#1373D3" (blue).
#'
#' @return A ggplot2 object containing the biplot visualization with:
#'   \itemize{
#'     \item Points representing species/samples in the ordination space
#'     \item Arrows showing the direction and magnitude of niche dimension contributions
#'     \item Labels for niche dimensions (arrows)
#'     \item Variance explained by PC1 and PC2 in axis labels
#'     \item Optional color coding by groups if provided
#'   }
#'
#' @details
#' The function creates a standard PCA biplot where:
#' \itemize{
#'   \item Points represent species positioned in the integrated niche space
#'   \item Red arrows represent the niche dimensions (from first-level PCAs) and
#'     their relative contribution to the ordination axes
#'   \item Arrow length indicates the strength of correlation with the ordination axes
#'   \item Arrow direction shows the gradient direction in niche space
#'   \item Dashed reference lines at x=0 and y=0 help interpret the ordination
#' }
#'
#' The plot helps interpret:
#' \itemize{
#'   \item Species clustering patterns in niche space
#'   \item Which niche dimensions drive the main gradients
#'   \item Relationships between different niche dimensions
#'   \item Group differences in niche occupation (when groups are specified)
#' }
#'
#' @note
#' \itemize{
#'   \item The function requires the result object from \code{NPT_continuous()}
#'   \item A larger plotting window is recommended for better visualization
#'   \item Arrow labels show dimension names (e.g., "pc1.grow", "pc2.survive")
#'   \item The function uses \code{max.overlaps = Inf} to show all labels
#'   \item Group colors are automatically assigned if groups are provided
#' }
#'
#' @references
#' 1. Winemiller, K. O., Fitzgerald, D. B., Bower, L. M., & Pianka, E. R. (2015). Functional traits,
#' convergent evolution, and periodic tables of niches. Ecology letters, 18(8), 737-751.
#' 2. Yu, R., Huang, J., Xu, Y., Ding, Y., & Zang, R. (2020). Plant functional niches in forests across
#' four climatic zones: Exploring the periodic table of niches based on plant functional traits.
#' Frontiers in Plant Science, 11, 841.
#'
#' @examples
#' data(PFF)
#' PFF[,4:21] <- log(PFF[,4:21])
#' PFF <- na.omit(PFF)
#' traits_dimension <- list(
#'   grow = c("SLA","SRL","Leaf_Nmass","Root_Nmass"),
#'   survive = c("Height","Leaf_CN","Root_CN"),
#'   reproductive = c("SeedMass","FltDate","FltDur")
#' )
#' npt_result <- NPT_continuous(data = PFF, dimension = traits_dimension)
#' NPT_continuous_plot(npt_result$result)
#' NPT_continuous_plot(npt_result$result, PFF$family)
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom vegan scores
#' @export
NPT_continuous_plot <- function(pca_obj, group = NULL, default_fill = "#1373D3") {
  # Extract site scores and add group information
  site_scores <- as.data.frame(vegan::scores(pca_obj, display = "sites"))
  site_scores$sample <- rownames(site_scores)
  # Handle group assignment
  if (is.null(group)) {
    site_scores$group <- "All Samples"
  } else {
    site_scores$group <- group
  }
  # Extract species scores
  species_scores <- as.data.frame(vegan::scores(pca_obj, display = "species"))
  species_scores$species <- rownames(species_scores)
  # Calculate variance explained by each PC
  var_explained <- pca_obj$CA$eig / sum(pca_obj$CA$eig) * 100
  # Create the plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = site_scores,
                        ggplot2::aes(x = .data$PC1, y = .data$PC2, fill = .data$group),
                        shape = 21, size = 2, alpha = 0.5) +
    ggplot2::geom_segment(data = species_scores,
                          ggplot2::aes(x = 0, y = 0, xend = .data$PC1, yend = .data$PC2),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.4, "cm")),
                          color = "red", linewidth = 0.8) +
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
      legend.key.size = ggplot2::unit(0.8, "lines"),
      legend.position = "right",
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text = ggplot2::element_text(size = 10, color = "black")
    ) +
    ggplot2::labs(
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      fill = "Group"
    )
  # Handle legend display based on group presence
  if (is.null(group)) {
    p <- p +
      ggplot2::scale_fill_manual(values = c("All Samples" = default_fill)) +
      ggplot2::theme(legend.position = "none")
  }
  return(p)
}
