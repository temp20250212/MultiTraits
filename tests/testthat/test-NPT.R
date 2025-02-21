# Load required libraries
library(testthat)
library(MultiTraits)

# Test data preparation
test_that("Data preparation works correctly", {
  # Load sample data
  data(PFF)
  PFF[,3:20] <- log(PFF[,3:20])

  # Define trait dimensions based on package examples
  traits_dimension <- list(
    grow = c("SLA","Leaf_area","LDMC","SRL","Leaf_Nmass","Leaf_Pmass","Root_Nmass"),
    survive = c("Height","Leaf_Cmass","Root_Cmass","Leaf_CN","Leaf_NP","Leaf_CP","Root_CN"),
    reproductive = c("SeedMass","FltDate","FltDur")
  )

  # Test input validation
  expect_error(NPT(data = NULL, dimension = traits_dimension))
  expect_error(NPT(data = PFF, dimension = NULL))
})

# Test NPT function core functionality
test_that("NPT function returns correct structure", {
  data(PFF)
  PFF[,3:20] <- log(PFF[,3:20])
  PFF <- na.omit(PFF)

  traits_dimension <- list(
    grow = c("SLA","Leaf_area"),
    survive = c("Height","Leaf_Cmass")
  )

  result <- NPT(data = PFF, dimension = traits_dimension)

  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("PCA_first", "PCA_second", "result"))

  # Check first level PCA results
  expect_s3_class(result$PCA_first, "data.frame")
  expect_true(all(c("pc1_percent", "pc2_percent") %in% colnames(result$PCA_first)))
})

# Test NPT_plot functionality
test_that("NPT_plot function works correctly", {
  data(PFF)
  PFF[,3:20] <- log(PFF[,3:20])
  PFF <- na.omit(PFF)

  traits_dimension <- list(
    grow = c("SLA","Leaf_area"),
    survive = c("Height","Leaf_Cmass")
  )

  result <- NPT(data = PFF, dimension = traits_dimension)

  # Test plot creation without groups
  p1 <- NPT_plot(result$result)
  expect_s3_class(p1, "ggplot")

  # Test plot creation with groups
  p2 <- NPT_plot(result$result, group = PFF$family)
  expect_s3_class(p2, "ggplot")
})

# Test NA handling
test_that("NA handling works correctly", {
  data(PFF)
  PFF[1,3] <- NA

  traits_dimension <- list(
    grow = c("SLA","Leaf_area"),
    survive = c("Height","Leaf_Cmass")
  )

  # Should remove rows with NA values
  expect_message(NPT(data = PFF, dimension = traits_dimension), "Removed")
})
