## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=5)

## ----class.source = 'fold-show'-----------------------------------------------
# install.packages("MultiTraits")

## ----class.source = 'fold-show'-----------------------------------------------
# if (!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}
# devtools::install_github("biodiversity-monitoring/MultiTraits")

## ----class.source = 'fold-show'-----------------------------------------------
library(MultiTraits)
data(PFF)
# View the structure of the datasets
head(PFF)

## ----class.source = 'fold-show'-----------------------------------------------
# Load packages
library(MultiTraits)
# Load the PFF dataset
data(PFF)
head(PFF)

# Create a data frame with the CSR traits
csr_traits <- data.frame(LA=PFF$Leaf_area, LDMC=PFF$LDMC, SLA=PFF$SLA)
rownames(csr_traits) <- PFF$species
head(csr_traits)

# Perform CSR analysis
csr_result <- CSR(csr_traits)
head(csr_result)

# Create the CSR plot
CSR_plot(csr_result)

## ----class.source = 'fold-show'-----------------------------------------------
# Load packages
library(MultiTraits)
# Load the PFF dataset
data(PFF)
head(PFF)

# Create a data frame with the LHS traits
lhs_traits <- PFF[, c("SLA", "Height", "SeedMass")]
rownames(lhs_traits) <- PFF$species
head(lhs_traits)

# Perform LHS analysis
lhs_result <- LHS(lhs_traits)
head(lhs_result)
table(lhs_result$LHS_strategy)


# Create a visualization plot of the LHS analysis results
LHS_plot(lhs_result, group = "LHS_strategy", show_cube = TRUE)

# Display the LHS strategy scheme diagram
LHS_strategy_scheme()

## ----class.source = 'fold-show'-----------------------------------------------
# Load packages
library(MultiTraits)
# Load and prepare data
data(PFF)
rownames(PFF) <- PFF$species
PFF_traits <- PFF[, c("SLA", "SRL", "Leaf_Nmass", "Root_Nmass","Height",
                      "Leaf_CN", "Root_CN","SeedMass", "FltDate", "FltDur")]
# Perform log transformation of data and remove missing values
PFF_traits <- log(na.omit(PFF_traits))
head(PFF_traits)

# Define trait dimensions
dimension <- list(
  Grow = c("SLA", "SRL", "Leaf_Nmass", "Root_Nmass"),
  Survive = c("Height", "Leaf_CN", "Root_CN"),
  Reproductive = c("SeedMass", "FltDate", "FltDur")
)

# Perform discrete niche scheme
set.seed(123)
discrete_result <- NPT_discrete(data = PFF_traits, dimension = dimension)
discrete_result
NPT_discrete_plot(discrete_result$niche_classification)

# Perform continuous niche scheme
continuous_result <- NPT_continuous(data = PFF_traits, dimension = dimension)
continuous_result 
NPT_continuous_plot(continuous_result$result)

## ----class.source = 'fold-show'-----------------------------------------------
# Load packages
library(MultiTraits)
# Load and prepare data
data(PFF)
rownames(PFF) <- PFF$species
PFF_traits <- PFF[, c("Leaf_area","LDMC","SLA", "Leaf_Cmass","Leaf_Nmass",
                      "Leaf_CN","Leaf_Pmass", "Leaf_NP","Leaf_CP")]
# Perform log transformation of data and remove missing values
PFF_traits <- log(na.omit(PFF_traits))
head(PFF_traits)

data(PFF_tree)
PFF_tree

# Calculate trait correlations using specified thresholds
PTN_corr(traits_matrix=PFF_traits, rThres = 0.2, pThres = 0.05,method = "pearson",
         phylo_correction = TRUE,phylo_tree = PFF_tree)

# Perform Trait Network (TN) analysis
PTN_phylo_result <- PTN(traits_matrix = PFF_traits,
                      rThres = 0.2,
                      pThres = 0.05,
                      method = "pearson",
                      phylo_correction = TRUE,
                      phylo_tree = PFF_tree)

# Calculate network metrics for the trait network
PTN_metrics(PTN_phylo_result)

set.seed(22)
# Create visualization plots of the trait network
PTN_plot(PTN_phylo_result, style = 1, vertex.size = 20, vertex.label.cex = 0.6)

