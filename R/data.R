
#' Plant Functional Traits Dataset from Ponderosa Pine Forests Flora (PFF)
#'
#' @description
#' A dataset containing functional traits for 133 plant species commonly found in
#' southwestern USA Pinus ponderosa var. scopulorum P. & C. Lawson (ponderosa pine) forests.
#'
#' @format A data frame with 137 rows and 20 variables:
#' \describe{
#'   \item{family}{Plant family}
#'   \item{species}{Plant species}
#'   \item{Height}{Canopy height (cm)}
#'   \item{Leaf_area}{Leaf area (mm^2)}
#'   \item{LDMC}{Leaf dry matter content (%)}
#'   \item{SLA}{Specific leaf area (mm^2/mg)}
#'   \item{SRL}{Specific root length (m/g)}
#'   \item{SeedMass}{Seed mass (mg)}
#'   \item{FltDate}{mean Flowering date (Julian day)}
#'   \item{FltDur}{mean Flowering duration (days)}
#'   \item{k_value}{decomposition decay constant, where proportion of original mass remaining = exp(- k-value*0.926)}
#'   \item{Leaf_Cmass}{Leaf carbon content (% dry mass)}
#'   \item{Leaf_Nmass}{Leaf nitrogen content (% dry mass)}
#'   \item{Leaf_CN}{Leaf carbon to nitrogen ratio}
#'   \item{Leaf_Pmass}{Leaf phosphorus content (% dry mass)}
#'   \item{Leaf_NP}{Leaf nitrogen to phosphorus ratio}
#'   \item{Leaf_CP}{Leaf carbon to phosphorus ratio}
#'   \item{Root_Cmass}{Root carbon content (% dry mass)}
#'   \item{Root_Nmass}{Root nitrogen content (% dry mass)}
#'   \item{Root_CN}{Root carbon to nitrogen ratio}
#' }
#'
#' @details
#' This dataset contains measurements of a core set of functional traits that reflect
#' aspects of each species' ability to disperse, establish, acquire water and nutrients,
#' and photosynthesize. Traits include specific leaf area (SLA), height, seed mass,
#' specific root length (SRL), leaf and fine root nitrogen concentration, leaf phosphorus
#' concentration, and leaf dry matter content (LDMC). Julian flowering date and flowering
#' duration were also obtained for each species. Leaf litter decomposition rates were
#' measured on 103 species.
#'
#' @source
#' Laughlin, D. C., Leppert, J. J., Moore, M. M., & Sieg, C. H. (2010). A multi-trait test
#' of the leaf-height-seed plant strategy scheme with 133 species from a pine forest flora.
#' Functional Ecology, 24(3), 485-700.
#'
#' @examples
#' data(PFF)
#' head(PFF)
"PFF"
