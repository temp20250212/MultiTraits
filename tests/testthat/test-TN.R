# Load required packages
library(testthat)
library(MultiTraits)

# Test data preparation
test_that("Test data preparation", {
  data(PFF)
  PFF_traits <- PFF[, c("Height", "Leaf_area", "LDMC", "SLA", "SRL",
                        "SeedMass", "FltDate", "FltDur")]

  # Check if data is properly formatted
  expect_s3_class(PFF_traits, "data.frame")
  expect_true(all(sapply(PFF_traits, is.numeric)))
})

# Test TN_corr function
test_that("Test TN_corr function", {
  data(PFF)
  PFF_traits <- PFF[, c("Height", "Leaf_area", "LDMC", "SLA")]
  PFF_traits <- na.omit(PFF_traits)

  # Test function execution
  corr_plot <- TN_corr(PFF_traits, rThres = 0.3, pThres = 0.01)

  # Check output type
  expect_type(corr_plot, "list")
})

# Test TN function
test_that("Test TN function", {
  data(PFF)
  PFF_traits <- PFF[, c("Height", "Leaf_area", "LDMC", "SLA")]
  PFF_traits <- na.omit(PFF_traits)

  # Test function execution with different methods
  g1 <- TN(PFF_traits, rThres = 0.2, pThres = 0.05, method = "pearson")
  g2 <- TN(PFF_traits, rThres = 0.2, pThres = 0.05, method = "spearman")

  # Check output type and structure
  expect_s3_class(g1, "igraph")
  expect_s3_class(g2, "igraph")
})

# Test TN_metrics function
test_that("Test TN_metrics function", {
  data(PFF)
  PFF_traits <- PFF[, c("Height", "Leaf_area", "LDMC", "SLA")]
  PFF_traits <- na.omit(PFF_traits)
  g <- TN(PFF_traits)

  metrics <- TN_metrics(g)

  # Check output structure
  expect_type(metrics, "list")
  expect_true(all(c("node", "global") %in% names(metrics)))
  expect_true(all(c("degree", "closeness", "betweenness", "clustering_coefficient") %in% names(metrics$node)))
})

# Test TN_plot function
test_that("Test TN_plot function", {
  data(PFF)
  PFF_traits <- PFF[, c("Height", "Leaf_area", "LDMC", "SLA")]
  PFF_traits <- na.omit(PFF_traits)
  g <- TN(PFF_traits)

  # Test different plotting styles
  expect_silent(TN_plot(g, style = 1))
  expect_silent(TN_plot(g, style = 2))
})

test_that("TN_plot throws error for invalid style parameter", {
  # Create a simple test graph
  test_graph <- igraph::make_ring(5)

  # Test with invalid style values
  expect_error(TN_plot(test_graph, style = 3),
               "Invalid style. Please choose 1 or 2.")

  expect_error(TN_plot(test_graph, style = 0),
               "Invalid style. Please choose 1 or 2.")

  expect_error(TN_plot(test_graph, style = "invalid"),
               "Invalid style. Please choose 1 or 2.")
})
