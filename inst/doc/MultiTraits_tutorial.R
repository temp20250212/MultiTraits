## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=5)

## ----class.source = 'fold-show'-----------------------------------------------
# #####

## ----class.source = 'fold-show'-----------------------------------------------
# if (!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}
# devtools::install_github("#######")

## ----class.source = 'fold-show'-----------------------------------------------
library(MultiTraits)
data(PFF)
# View the structure of the datasets
head(PFF)

## ----class.source = 'fold-show'-----------------------------------------------
# Load the PFF dataset
data(PFF)
head(PFF)

# Select required traits for CSR analysis
traits <- data.frame(LA=PFF$Leaf_area, LDMC=PFF$LDMC, SLA=PFF$SLA)
head(traits)

# Perform CSR analysis
result <- CSR(data = traits)
head(result)

# Visualize CSR strategy results
CSR_plot(data=result, expand_margin = 1)

## ----class.source = 'fold-show'-----------------------------------------------
# Load the PFF dataset
data(PFF)
# Select specific columns (SLA, Height, SeedMass) from the PFF dataset
pff <- PFF[, c("SLA", "Height", "SeedMass")]
head(pff)

# Perform LHS (Leaf-Height-Seed) analysis on the selected data
result <- LHS(pff)
head(result)
table(result$LHS_strategy)

# Create a visualization plot of the LHS analysis results
LHS_plot(result)

# Display the LHS strategy scheme diagram
LHS_strategy_scheme()

## ----class.source = 'fold-show'-----------------------------------------------
# Load the PFF dataset
data(PFF)
# Log-transform columns 3-20 of the dataset
PFF[,3:20] <- log(PFF[,3:20])
# Remove rows with missing values (NA)
PFF <- na.omit(PFF)
head(PFF)

# Define trait dimensions for NPT analysis
traits_dimension <-list(
  grow = c("SLA","Leaf_area","LDMC","SRL","Leaf_Nmass","Leaf_Pmass","Root_Nmass"),
  survive = c("Height","Leaf_Cmass","Root_Cmass","Leaf_CN","Leaf_NP","Leaf_CP","Root_CN"),
  reproductive = c("SeedMass","FltDate","FltDur"))
# Perform NPT analysis using the defined dimensions
npt_result <- NPT(data = PFF, dimension = traits_dimension)
npt_result

# Create a basic visualization of NPT results
NPT_plot(npt_result$result)
# Create a visualization of NPT results colored by plant family
NPT_plot(npt_result$result, PFF$family)

## ----class.source = 'fold-show'-----------------------------------------------
# Load the PFF dataset
data(PFF)

# Select specific trait columns for analysis
PFF_traits <- PFF[, c("Height", "Leaf_area","LDMC","SLA","SRL","SeedMass","FltDate","FltDur","Leaf_Cmass","Leaf_Nmass",
                      "Leaf_CN","Leaf_Pmass","Leaf_NP","Leaf_CP","Root_Cmass","Root_Nmass","Root_CN") ]
# Perform log transformation of data and remove missing values
PFF_traits <- log(na.omit(PFF_traits))
head(PFF_traits)

# Calculate trait correlations using specified thresholds
TN_corr(traits_matrix=PFF_traits, rThres = 0.2, pThres = 0.05,method = "pearson")

# Perform Trait Network (TN) analysis
Tn_result <- TN(traits_matrix = PFF_traits, rThres = 0.2, pThres = 0.05,method = "pearson")
Tn_result

# Calculate network metrics for the trait network
TN_metrics(Tn_result)

# Create visualization plots of the trait network
par(mfrow=c(1,2))
# Style 1
TN_plot(Tn_result, style = 1,vertex.size = 10,vertex.label.cex = 0.6)
# Style 2
TN_plot(Tn_result, style = 2,vertex.size = 20,vertex.label.cex = 0.6)

