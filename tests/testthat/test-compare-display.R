# test-garch.R
# Tests for GARCH comparison and display functions

library(testthat)

# ==============================================================================
# Test Parameters
# ==============================================================================

TEST_SEED <- 42
N_SMALL <- 500   # Smaller sample for faster tests
N_LARGE <- 1000  # Larger sample for more reliable fits

# Helper to generate data that GARCH models can actually fit
generate_garch_data <- function(n, seed = TEST_SEED) {
  garch_gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.15, beta = 0.75)
  set.seed(seed)
  garch_gen(n)
}


# ==============================================================================
# Helper Function Tests (No external dependencies)
# ==============================================================================

test_that(".garch_label generates correct labels", {
  # Access internal function
  
  expect_equal(.garch_label(1, 0), "ARCH(1)")
  expect_equal(.garch_label(2, 0), "ARCH(2)")
  expect_equal(.garch_label(1, 1), "GARCH(1,1)")
  expect_equal(.garch_label(2, 1), "GARCH(2,1)")
  expect_equal(.garch_label(1, 2), "GARCH(1,2)")
  expect_equal(.garch_label(0, 1), "GARCH(0,1)")
})

test_that(".fmt_pval formats p-values correctly", {
  expect_equal(.fmt_pval(0.5), "0.500")
  expect_equal(.fmt_pval(0.05), "0.050")
  expect_equal(.fmt_pval(0.001), "0.001")
  expect_equal(.fmt_pval(0.0001), "<0.001")
  expect_equal(.fmt_pval(0.0009), "<0.001")
  expect_equal(.fmt_pval(NA), "NA")
})

test_that(".fmt_pval respects digits argument", {
  expect_equal(.fmt_pval(0.12345, digits = 2), "0.12")
  expect_equal(.fmt_pval(0.12345, digits = 4), "0.1235")
})

test_that(".find_best_indices handles normal data", {
  df <- data.frame(
    AIC = c(100, 90, 95),
    BIC = c(110, 105, 100),
    AICC = c(101, 91, 96),
    WLB1 = c(0.1, 0.5, 0.3),
    WLB2 = c(0.2, 0.4, 0.6),
    WLB3 = c(0.3, 0.3, 0.5),
    Nyblom = c(1.0, 0.5, 0.8),
    SignBias = c(0.1, 0.2, 0.3)
  )
  
  best <- .find_best_indices(df)
  
  expect_equal(best$aic, 2)      # 90 is min
  expect_equal(best$bic, 3)      # 100 is min
  expect_equal(best$aicc, 2)     # 91 is min
  expect_equal(best$wlb1, 2)     # 0.5 is max
  expect_equal(best$wlb2, 3)     # 0.6 is max
  expect_equal(best$wlb3, 3)     # 0.5 is max
  expect_equal(best$nyblom, 2)   # 0.5 is min
  expect_equal(best$signbias, 3) # 0.3 is max
})

test_that(".find_best_indices handles NA values", {
  df <- data.frame(
    AIC = c(100, NA, 95),
    BIC = c(NA, NA, NA),  # All NA
    AICC = c(101, 91, NA),
    WLB1 = c(NA, 0.5, 0.3),
    WLB2 = c(0.2, NA, 0.6),
    WLB3 = c(NA, NA, NA),  # All NA
    Nyblom = c(1.0, 0.5, NA),
    SignBias = c(NA, NA, 0.3)
  )
  
  best <- .find_best_indices(df)
  
  expect_equal(best$aic, 3)           # 95 is min (ignoring NA)
  expect_true(is.na(best$bic))        # All NA returns NA
  expect_equal(best$aicc, 2)          # 91 is min
  expect_equal(best$wlb1, 2)          # 0.5 is max
  expect_equal(best$wlb2, 3)          # 0.6 is max
  expect_true(is.na(best$wlb3))       # All NA returns NA
  expect_equal(best$nyblom, 2)        # 0.5 is min
  expect_equal(best$signbias, 3)      # 0.3 is max (only non-NA)
})


# ==============================================================================
# Tests: compare.garch.tse - Input Validation
# ==============================================================================

test_that("compare.garch.tse validates data input", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  expect_error(compare.garch.tse("not numeric"), "numeric vector")
  expect_error(compare.garch.tse(1:5), "at least 10 observations")
  expect_error(compare.garch.tse(c()), "numeric vector")
})

test_that("compare.garch.tse validates order ranges", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_SMALL)
  
  # Only (0,0) in grid should fail
  expect_error(
    compare.garch.tse(y, arch.range = 0, garch.range = 0),
    "No valid model orders"
  )
})


# ==============================================================================
# Tests: compare.garch.tse - Basic Functionality
# ==============================================================================

test_that("compare.garch.tse returns correct structure", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_SMALL)
  
  # Fit minimal grid for speed
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 0:1)
  )
  
  expect_s3_class(results, "garch.comparison.tse")
  expect_true("fits" %in% names(results))
  expect_true("comparison" %in% names(results))
  expect_true("distribution" %in% names(results))
  expect_true("n" %in% names(results))
  
  expect_type(results$fits, "list")
  expect_s3_class(results$comparison, "tbl_df")
  expect_equal(results$n, N_SMALL)
})

test_that("compare.garch.tse fits correct number of models", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_SMALL)
  
  # arch.range = 0:1, garch.range = 0:1 gives 4 combos minus (0,0) = 3 models
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 0:1, garch.range = 0:1)
  )
  
  expect_equal(nrow(results$comparison), 3)
  expect_equal(length(results$fits), 3)
  
  # Check model names
  expect_true("ARCH(1)" %in% results$comparison$Model)
  expect_true("GARCH(0,1)" %in% results$comparison$Model)
  expect_true("GARCH(1,1)" %in% results$comparison$Model)
})

test_that("compare.garch.tse comparison tibble has correct columns", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_SMALL)
  
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  expected_cols <- c("Model", "ARCH", "GARCH", "AIC", "AICC", "BIC",
                     "WLB1", "WLB2", "WLB3", "Nyblom", "Nyblom_crit",
                     "SignBias", "n_sig", "n_coef")
  
  for (col in expected_cols) {
    expect_true(col %in% names(results$comparison),
                info = paste("Missing column:", col))
  }
})

test_that("compare.garch.tse IC values are finite", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_SMALL)
  
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 0:1)
  )
  
  expect_true(all(is.finite(results$comparison$AIC)))
  expect_true(all(is.finite(results$comparison$BIC)))
  expect_true(all(is.finite(results$comparison$AICC)))
})

test_that("compare.garch.tse AICc > AIC (small sample correction)", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_SMALL)
  
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  # AICc should be greater than AIC due to small sample correction
  expect_true(all(results$comparison$AICC > results$comparison$AIC))
})


# ==============================================================================
# Tests: compare.garch.tse - Distribution Options
# ==============================================================================

test_that("compare.garch.tse respects distribution argument", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_SMALL)
  
  results_norm <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1, distribution = "norm")
  )
  results_std <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1, distribution = "std")
  )
  
  expect_equal(results_norm$distribution, "norm")
  expect_equal(results_std$distribution, "std")
  
  # t-distribution should have different IC (extra parameter)
  expect_false(results_norm$comparison$AIC[1] == results_std$comparison$AIC[1])
})


# ==============================================================================
# Tests: compare.garch.tse - Print and Summary Methods
# ==============================================================================

test_that("print.garch.comparison.tse works", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_SMALL)
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 0:1)
  )
  
  expect_output(print(results), "GARCH Model Comparison")
  expect_output(print(results), "Distribution:")
  expect_output(print(results), "Sample size:")
  expect_output(print(results), "Best by AIC:")
})

test_that("summary.garch.comparison.tse works", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_SMALL)
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  expect_output(summary(results), "Model")
  expect_output(summary(results), "AIC")
})


# ==============================================================================
# Tests: compare.garch.tse - With GARCH Data
# ==============================================================================

test_that("compare.garch.tse works with actual GARCH data", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_LARGE)
  
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 0:1, garch.range = 0:1)
  )
  
  expect_s3_class(results, "garch.comparison.tse")
  expect_gt(nrow(results$comparison), 0)
  
  # GARCH(1,1) should generally fit well - check WLB p-values aren't all tiny
  garch11_row <- results$comparison[results$comparison$Model == "GARCH(1,1)", ]
  expect_gt(garch11_row$WLB1, 0.01)  # Shouldn't strongly reject
})


# ==============================================================================
# Tests: table.garch.gt.tse - Input Validation
# ==============================================================================

test_that("table.garch.gt.tse validates input class", {
  skip_if_not_installed("gt")
  
  expect_error(
    table.garch.gt.tse(list(a = 1)),
    "garch.comparison.tse"
  )
  
  expect_error(
    table.garch.gt.tse(data.frame(x = 1)),
    "garch.comparison.tse"
  )
})

test_that("table.garch.gt.tse errors on empty comparison", {
  skip_if_not_installed("gt")
  
  empty_results <- structure(
    list(
      fits = list(),
      comparison = data.frame(),
      distribution = "norm",
      n = 100
    ),
    class = "garch.comparison.tse"
  )
  
  expect_error(table.garch.gt.tse(empty_results), "empty")
})


# ==============================================================================
# Tests: table.garch.gt.tse - Output
# ==============================================================================

test_that("table.garch.gt.tse returns gt object", {
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
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  y <- generate_garch_data(N_SMALL)
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  tbl <- table.garch.gt.tse(results, title = "Custom Title")
  
  # gt object should be created without error
  expect_s3_class(tbl, "gt_tbl")
})

test_that("table.garch.gt.tse handles NA values gracefully", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  # Create results with some NA values manually
  mock_comparison <- dplyr::tibble(
    Model = c("ARCH(1)", "GARCH(1,1)"),
    ARCH = c(1, 1),
    GARCH = c(0, 1),
    AIC = c(100, 95),
    AICC = c(101, 96),
    BIC = c(105, 100),
    WLB1 = c(0.5, NA),
    WLB2 = c(NA, 0.3),
    WLB3 = c(0.2, 0.4),
    Nyblom = c(0.5, NA),
    Nyblom_crit = c(1.0, 1.0),
    SignBias = c(NA, 0.1),
    n_sig = c(2, 3),
    n_coef = c(2, 3)
  )
  
  mock_results <- structure(
    list(
      fits = list(),
      comparison = mock_comparison,
      distribution = "norm",
      n = 500
    ),
    class = "garch.comparison.tse"
  )
  
  # Should not error with NA values
  tbl <- table.garch.gt.tse(mock_results)
  expect_s3_class(tbl, "gt_tbl")
})


# ==============================================================================
# Tests: table.garch.cli.tse - Input Validation
# ==============================================================================

test_that("table.garch.cli.tse validates input class", {
  skip_if_not_installed("cli")
  
  expect_error(
    table.garch.cli.tse(list(a = 1)),
    "garch.comparison.tse"
  )
})

test_that("table.garch.cli.tse handles empty comparison", {
  skip_if_not_installed("cli")
  
  empty_results <- structure(
    list(
      fits = list(),
      comparison = data.frame(),
      distribution = "norm",
      n = 100
    ),
    class = "garch.comparison.tse"
  )
  
  expect_message(
    result <- table.garch.cli.tse(empty_results),
    "No models to display"
  )
})


# ==============================================================================
# Tests: table.garch.cli.tse - Output
# ==============================================================================
test_that("table.garch.cli.tse produces output", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("cli")
  
  y <- generate_garch_data(N_SMALL)
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 0:1)
  )
  
  # Check for table content (cli headers may not capture in test mode)
  expect_output(table.garch.cli.tse(results), "Model")
  expect_output(table.garch.cli.tse(results), "AIC")
  expect_output(table.garch.cli.tse(results), "ARCH|GARCH")  # Model names
})

test_that("table.garch.cli.tse returns comparison invisibly", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("cli")
  
  y <- generate_garch_data(N_SMALL)
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  captured <- capture.output(
    returned <- table.garch.cli.tse(results)
  )
  
  expect_s3_class(returned, "tbl_df")
  expect_equal(nrow(returned), 1)
})

test_that("table.garch.cli.tse show.signbias option works", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("cli")
  
  y <- generate_garch_data(N_SMALL)
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  output_with <- capture.output(table.garch.cli.tse(results, show.signbias = TRUE))
  output_without <- capture.output(table.garch.cli.tse(results, show.signbias = FALSE))
  
  expect_true(any(grepl("SignBias", output_with)))
  expect_false(any(grepl("SignBias", output_without)))
})

test_that("table.garch.cli.tse handles NA values gracefully", {
  skip_if_not_installed("cli")
  skip_if_not_installed("dplyr")
  
  mock_comparison <- dplyr::tibble(
    Model = c("ARCH(1)", "GARCH(1,1)"),
    ARCH = c(1, 1),
    GARCH = c(0, 1),
    AIC = c(100, 95),
    AICC = c(101, 96),
    BIC = c(105, 100),
    WLB1 = c(0.5, NA),
    WLB2 = c(NA, 0.3),
    WLB3 = c(0.2, 0.4),
    Nyblom = c(0.5, NA),
    Nyblom_crit = c(1.0, 1.0),
    SignBias = c(NA, 0.1),
    n_sig = c(2, 3),
    n_coef = c(2, 3)
  )
  
  mock_results <- structure(
    list(
      fits = list(),
      comparison = mock_comparison,
      distribution = "norm",
      n = 500
    ),
    class = "garch.comparison.tse"
  )
  
  # Should not error with NA values
  expect_output(table.garch.cli.tse(mock_results), "GARCH")
})


# ==============================================================================
# Tests: table.coef.gt.tse - Input Validation
# ==============================================================================

test_that("table.coef.gt.tse validates input class", {
  skip_if_not_installed("gt")
  
  expect_error(
    table.coef.gt.tse(list(a = 1)),
    "uGARCHfit"
  )
})


# ==============================================================================
# Tests: table.coef.gt.tse - Output
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

test_that("table.coef.gt.tse works with different distributions", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  y <- generate_garch_data(N_SMALL)
  
  # Test with t-distribution (has extra shape parameter)
  results_std <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1, distribution = "std")
  )
  
  fit <- results_std$fits[["GARCH(1,1)"]]
  tbl <- table.coef.gt.tse(fit)
  
  expect_s3_class(tbl, "gt_tbl")
})


# ==============================================================================
# Tests: Integration - Full Workflow
# ==============================================================================

test_that("full workflow: generate -> compare -> display (gt)", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  # Generate GARCH data
  y <- generate_garch_data(N_LARGE)
  
  # Compare models
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 0:1, garch.range = 0:1)
  )
  
  # Display comparison table
  tbl_comparison <- table.garch.gt.tse(results)
  expect_s3_class(tbl_comparison, "gt_tbl")
  
  # Display coefficient table for best model
  best_model <- results$comparison$Model[which.min(results$comparison$BIC)]
  best_fit <- results$fits[[best_model]]
  tbl_coef <- table.coef.gt.tse(best_fit)
  expect_s3_class(tbl_coef, "gt_tbl")
})

test_that("full workflow: generate -> compare -> display (cli)", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("cli")
  
  # Generate GARCH data
  y <- generate_garch_data(N_LARGE)
  
  # Compare models
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 0:1, garch.range = 0:1)
  )
  
  # Display in console - check for model names in output
  expect_output(table.garch.cli.tse(results), "GARCH\\(1,1\\)|ARCH\\(1\\)")
})


# ==============================================================================
# Tests: Edge Cases
# ==============================================================================

test_that("compare.garch.tse handles single model in grid", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_SMALL)
  
  # Only fit GARCH(1,1)
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  expect_equal(nrow(results$comparison), 1)
  expect_equal(results$comparison$Model[1], "GARCH(1,1)")
})

test_that("compare.garch.tse handles larger grid", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_SMALL)
  
  # 3x3 grid minus (0,0) = 8 models
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 0:2, garch.range = 0:2)
  )
  
  expect_equal(nrow(results$comparison), 8)
})

test_that("display functions handle all-significant coefficients", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  # Use GARCH data where all params should be significant
  y <- generate_garch_data(N_LARGE)
  
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  # Should not error
  tbl <- table.garch.gt.tse(results)
  expect_s3_class(tbl, "gt_tbl")
})

test_that("display functions handle non-significant coefficients", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  # Use white noise where GARCH params won't be significant
  # But generate enough data that fitting will converge
  set.seed(TEST_SEED)
  y <- rnorm(N_LARGE)
  
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1, garch.range = 1)
  )
  
  # Should not error even with non-significant params
  tbl <- table.garch.gt.tse(results)
  expect_s3_class(tbl, "gt_tbl")
})


# ==============================================================================
# Tests: Consistency Checks
# ==============================================================================

test_that("BIC penalizes complexity more than AIC", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  # Use GARCH data for reliable convergence
  y <- generate_garch_data(N_LARGE)
  
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1:2, garch.range = 0:1)
  )
  
  # For each pair of models, check that more complex model
  # has larger BIC-AIC difference
  comp <- results$comparison
  
  # ARCH(1) vs GARCH(1,1)
  arch1 <- comp[comp$Model == "ARCH(1)", ]
  garch11 <- comp[comp$Model == "GARCH(1,1)", ]
  
  # BIC penalty should be larger for GARCH(1,1) relative to AIC
  # (BIC - AIC) should be larger for more complex model
  bic_aic_diff_arch1 <- arch1$BIC - arch1$AIC
  bic_aic_diff_garch11 <- garch11$BIC - garch11$AIC
  
  expect_gt(bic_aic_diff_garch11, bic_aic_diff_arch1)
})

test_that("coefficient counts are correct", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  y <- generate_garch_data(N_SMALL)
  
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 1:2, garch.range = 0:1)
  )
  
  comp <- results$comparison
  
  # ARCH(1): omega + alpha1 = 2 params
  arch1 <- comp[comp$Model == "ARCH(1)", ]
  expect_equal(arch1$n_coef, 2)
  
  # ARCH(2): omega + alpha1 + alpha2 = 3 params
  arch2 <- comp[comp$Model == "ARCH(2)", ]
  expect_equal(arch2$n_coef, 3)
  
  # GARCH(1,1): omega + alpha1 + beta1 = 3 params
  garch11 <- comp[comp$Model == "GARCH(1,1)", ]
  expect_equal(garch11$n_coef, 3)
  
  # GARCH(2,1): omega + alpha1 + alpha2 + beta1 = 4 params
  garch21 <- comp[comp$Model == "GARCH(2,1)", ]
  expect_equal(garch21$n_coef, 4)
})