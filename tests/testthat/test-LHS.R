# Load required packages
library(testthat)
library(MultiTraits)

# Test data setup
test_that("Test data creation", {
  test_data <- data.frame(
    SLA = c(10, 20, 15, 25),
    Height = c(5, 15, 10, 20),
    SeedMass = c(1, 3, 2, 4)
  )
  expect_s3_class(test_data, "data.frame")
})

# Test LHS function
test_that("Test LHS function behavior", {
  # Test basic functionality
  test_data <- data.frame(
    SLA = c(10, 20, 15, 25),
    Height = c(5, 15, 10, 20),
    SeedMass = c(1, 3, 2, 4)
  )
  result <- LHS(test_data)

  # Check output structure
  expect_true("LHS_strategy" %in% names(result))
  expect_equal(nrow(result), nrow(test_data))

  # Check strategy assignment
  expect_true(all(result$LHS_strategy %in%
                    c("S-S-S", "L-L-L", "S-S-L", "L-L-S", "S-L-S", "L-S-L", "S-L-L", "L-S-S")))
})

# Test LHS_plot function
test_that("Test LHS_plot function", {
  test_data <- data.frame(
    SLA = c(10, 20, 15, 25),
    Height = c(5, 15, 10, 20),
    SeedMass = c(1, 3, 2, 4)
  )
  test_data <- LHS(test_data)

  # Test plot creation
  expect_error(LHS_plot(test_data), NA)
  expect_error(LHS_plot(test_data, log_transform = FALSE), NA)

  # Test error handling for missing columns
  bad_data <- data.frame(SLA = c(1,2))
  expect_error(LHS_plot(bad_data))
})

# Test LHS_strategy_scheme function
test_that("Test LHS_strategy_scheme function", {
  scheme <- LHS_strategy_scheme()

  # Check output structure
  expect_s3_class(scheme, "data.frame")
  expect_equal(ncol(scheme), 2)
  expect_equal(nrow(scheme), 8)
  expect_true(all(c("type", "strategy") %in% names(scheme)))
})
