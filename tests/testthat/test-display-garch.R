# test-display-garch.R
# Tests for GARCH display functions
# Covers: table.garch.gt.tse, table.garch.cli.tse, table.coef.gt.tse

library(testthat)

# ==============================================================================
# table.garch.gt.tse - Input Validation
# ==============================================================================

test_that("table.garch.gt.tse validates input class", {
  skip_if_not_installed("gt")
  
  expect_error(table.garch.gt.tse(list(a = 1)), "garch.comparison.tse")
  expect_error(table.garch.gt.tse(data.frame(x = 1)), "garch.comparison.tse")
})

test_that("table.garch.gt.tse errors on empty comparison", {
  skip_if_not_installed("gt")
  
  empty_results <- create_empty_garch_comparison()
  expect_error(table.garch.gt.tse(empty_results), "empty")
})


# ==============================================================================
# table.garch.gt.tse - Output
# ==============================================================================

test_that("table.garch.gt.tse returns gt object with mock data", {
  skip_if_not_installed("gt")
  skip_if_not_installed("dplyr")
  
  mock_results <- create_mock_garch_comparison()
  tbl <- table.garch.gt.tse(mock_results)
  expect_s3_class(tbl, "gt_tbl")
})

test_that("table.garch.gt.tse returns gt object with real data", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  y <- generate_garch_data(N_SMALL)
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 0:1)
  )
  
  tbl <- table.garch.gt.tse(results)
  expect_s3_class(tbl, "gt_tbl")
})

test_that("table.garch.gt.tse accepts custom title", {
  skip_if_not_installed("gt")
  skip_if_not_installed("dplyr")
  
  mock_results <- create_mock_garch_comparison()
  tbl <- table.garch.gt.tse(mock_results, title = "Custom Title")
  expect_s3_class(tbl, "gt_tbl")
})

test_that("table.garch.gt.tse accepts custom colors", {
  skip_if_not_installed("gt")
  skip_if_not_installed("dplyr")
  
  mock_results <- create_mock_garch_comparison()
  tbl <- table.garch.gt.tse(
    mock_results,
    color.best = "darkgreen",
    color.fail = "red"
  )
  expect_s3_class(tbl, "gt_tbl")
})

test_that("table.garch.gt.tse handles NA values gracefully", {
  skip_if_not_installed("gt")
  skip_if_not_installed("dplyr")
  
  mock_results <- create_mock_garch_comparison(with_na = TRUE)
  tbl <- table.garch.gt.tse(mock_results)
  expect_s3_class(tbl, "gt_tbl")
})


# ==============================================================================
# table.garch.cli.tse - Input Validation
# ==============================================================================

test_that("table.garch.cli.tse validates input class", {
  skip_if_not_installed("cli")
  
  expect_error(table.garch.cli.tse(list(a = 1)), "garch.comparison.tse")
})

test_that("table.garch.cli.tse handles empty comparison", {
  skip_if_not_installed("cli")
  
  empty_results <- create_empty_garch_comparison()
  
  expect_message(
    result <- table.garch.cli.tse(empty_results),
    "No models to display"
  )
})


# ==============================================================================
# table.garch.cli.tse - Output
# ==============================================================================

test_that("table.garch.cli.tse produces output with mock data", {
  skip_if_not_installed("cli")
  skip_if_not_installed("dplyr")
  
  mock_results <- create_mock_garch_comparison()
  
  expect_output(table.garch.cli.tse(mock_results), "Model")
  expect_output(table.garch.cli.tse(mock_results), "AIC")
  expect_output(table.garch.cli.tse(mock_results), "ARCH|GARCH")
})

test_that("table.garch.cli.tse produces output with real data", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("cli")
  
  y <- generate_garch_data(N_SMALL)
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 0:1)
  )
  
  expect_output(table.garch.cli.tse(results), "Model")
  expect_output(table.garch.cli.tse(results), "AIC")
  expect_output(table.garch.cli.tse(results), "ARCH|GARCH")
})

test_that("table.garch.cli.tse returns comparison invisibly", {
  skip_if_not_installed("cli")
  skip_if_not_installed("dplyr")
  
  mock_results <- create_mock_garch_comparison()
  
  captured <- capture.output(
    returned <- table.garch.cli.tse(mock_results)
  )
  
  expect_s3_class(returned, "tbl_df")
  expect_equal(nrow(returned), 2)
})

test_that("table.garch.cli.tse show.signbias option works", {
  skip_if_not_installed("cli")
  skip_if_not_installed("dplyr")
  
  mock_results <- create_mock_garch_comparison()
  
  output_with <- capture.output(
    table.garch.cli.tse(mock_results, show.signbias = TRUE)
  )
  output_without <- capture.output(
    table.garch.cli.tse(mock_results, show.signbias = FALSE)
  )
  
  expect_true(any(grepl("SignBias", output_with)))
  expect_false(any(grepl("SignBias", output_without)))
})

test_that("table.garch.cli.tse handles NA values gracefully", {
  skip_if_not_installed("cli")
  skip_if_not_installed("dplyr")
  
  mock_results <- create_mock_garch_comparison(with_na = TRUE)
  expect_output(table.garch.cli.tse(mock_results), "GARCH")
})


# ==============================================================================
# table.coef.gt.tse - Input Validation
# ==============================================================================

test_that("table.coef.gt.tse validates input class", {
  skip_if_not_installed("gt")
  
  expect_error(table.coef.gt.tse(list(a = 1)), "uGARCHfit")
})


# ==============================================================================
# table.coef.gt.tse - Output
# ==============================================================================

test_that("table.coef.gt.tse returns gt object", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  y <- generate_garch_data(N_SMALL)
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  fit <- results$fits[["GARCH(1,1)"]]
  tbl <- table.coef.gt.tse(fit)
  
  expect_s3_class(tbl, "gt_tbl")
})

test_that("table.coef.gt.tse accepts custom title", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  y <- generate_garch_data(N_SMALL)
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  fit <- results$fits[["GARCH(1,1)"]]
  tbl <- table.coef.gt.tse(fit, title = "My Custom Title")
  
  expect_s3_class(tbl, "gt_tbl")
})

test_that("table.coef.gt.tse accepts custom colors", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  y <- generate_garch_data(N_SMALL)
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  fit <- results$fits[["GARCH(1,1)"]]
  tbl <- table.coef.gt.tse(fit, color.sig = "blue", color.nonsig = "gray")
  
  expect_s3_class(tbl, "gt_tbl")
})

test_that("table.coef.gt.tse works with different distributions", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  y <- generate_garch_data(N_SMALL)
  
  results_std <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1, distribution = "std")
  )
  
  fit <- results_std$fits[["GARCH(1,1)"]]
  tbl <- table.coef.gt.tse(fit)
  
  expect_s3_class(tbl, "gt_tbl")
})


# ==============================================================================
# Display Functions with All-Significant/Non-Significant Coefficients
# ==============================================================================

test_that("display functions handle all-significant coefficients", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  y <- generate_garch_data(1000)
  
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  tbl <- table.garch.gt.tse(results)
  expect_s3_class(tbl, "gt_tbl")
})

test_that("display functions handle non-significant coefficients", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  set.seed(TEST_SEED)
  y <- rnorm(1000)
  
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  tbl <- table.garch.gt.tse(results)
  expect_s3_class(tbl, "gt_tbl")
})