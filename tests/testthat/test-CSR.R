# Load required packages
library(testthat)
library(MultiTraits)

# Test Strate_CSR function
test_that("Strate_CSR function works correctly", {
  # Test basic functionality
  result <- Strate_CSR(LA = 369615.7, LDMC = 25.2, SLA = 17.4)

  # Check return value type
  expect_type(result, "list")
  expect_named(result, c("c_proportion", "s_proportion", "r_proportion", "csr_type"))

  # Check value ranges
  expect_true(all(result$c_proportion >= 0 & result$c_proportion <= 100))
  expect_true(all(result$s_proportion >= 0 & result$s_proportion <= 100))
  expect_true(all(result$r_proportion >= 0 & result$r_proportion <= 100))
})

# Test CSR function
test_that("CSR function handles data frame input correctly", {
  # Create test data
  test_data <- data.frame(
    LA = c(369615.7, 11.8, 55.7),
    LDMC = c(25.2, 39.7, 13.3),
    SLA = c(17.4, 6.6, 34.1)
  )

  result <- CSR(test_data)

  # Check return structure
  expect_s3_class(result, "data.frame")
  expect_true(all(c("LA", "LDMC", "SLA", "C", "S", "R", "type") %in% colnames(result)))

  # Check error handling
  expect_error(CSR(list()), "The input must be a dataframe")
  expect_error(CSR(data.frame(LA = 1)), "The input data frame must contain LA, LDMC and SLA columns")
})


# Test CSR_plot function
test_that("CSR_plot function creates valid ggplot object", {
  test_data <- data.frame(
    LA = c(369615.7, 11.8),
    LDMC = c(25.2, 39.7),
    SLA = c(17.4, 6.6)
  )
  result <- CSR(test_data)
  plot <- CSR_plot(result)

  expect_s3_class(plot, "ggplot")
})
