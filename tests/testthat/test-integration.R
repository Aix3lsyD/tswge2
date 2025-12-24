# test-integration.R
# Integration tests spanning multiple components
# Tests complete workflows from data generation to display

library(testthat)

# ==============================================================================
# Full Workflow: Generate -> Compare -> Display (gt)
# ==============================================================================

test_that("full workflow: GARCH generation -> comparison -> gt display", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("gt")
  
  # Step 1: Generate GARCH data
  y <- generate_garch_data(1000)
  expect_length(y, 1000)
  expect_true(all(is.finite(y)))
  
  # Step 2: Compare models
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 0:1, garch.range = 0:1)
  )
  expect_s3_class(results, "garch.comparison.tse")
  expect_gt(nrow(results$comparison), 0)
  
  # Step 3: Display comparison table
  tbl_comparison <- table.garch.gt.tse(results)
  expect_s3_class(tbl_comparison, "gt_tbl")
  
  # Step 4: Display coefficient table for best model
  best_model <- results$comparison$Model[which.min(results$comparison$BIC)]
  best_fit <- results$fits[[best_model]]
  tbl_coef <- table.coef.gt.tse(best_fit)
  expect_s3_class(tbl_coef, "gt_tbl")
})


# ==============================================================================
# Full Workflow: Generate -> Compare -> Display (cli)
# ==============================================================================

test_that("full workflow: GARCH generation -> comparison -> cli display", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  skip_if_not_installed("cli")
  
  # Step 1: Generate GARCH data
  y <- generate_garch_data(1000)
  
  # Step 2: Compare models
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 0:1, garch.range = 0:1)
  )
  
  # Step 3: Display in console
  expect_output(table.garch.cli.tse(results), "GARCH\\(1,1\\)|ARCH\\(1\\)")
})


# ==============================================================================
# Full Workflow: ARUMA with GARCH Innovations
# ==============================================================================

test_that("workflow: ARMA-GARCH simulation with custom innovations", {
  skip_if_not_installed("rugarch")
  
  # Create GARCH innovation generator
  garch_gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.15, beta = 0.8)
  
  # Generate ARMA(1,1)-GARCH(1,1) process
  result <- gen.aruma.tse(
    n = 500,
    phi = 0.7,
    theta = 0.3,
    innov_gen = garch_gen,
    sn = 42,
    plot = FALSE
  )
  
  expect_s3_class(result, "aruma")
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
  
  # Verify ARMA structure is present
  acf_vals <- acf(result$y, lag.max = 5, plot = FALSE)$acf[-1]
  expect_true(abs(acf_vals[1]) > 0.3)  # Should show autocorrelation
})


# ==============================================================================
# Full Workflow: Heavy-Tailed Innovations Through ARMA
# ==============================================================================

test_that("workflow: compare tail behavior with different innovations", {
  n <- 2000
  phi <- 0.5
  
  # Normal innovations
  norm_gen <- make.gen.norm.tse()
  result_norm <- gen.aruma.tse(
    n = n, phi = phi, innov_gen = norm_gen, sn = 42, plot = FALSE
  )
  
  # Heavy-tailed innovations
  t_gen <- make.gen.t.tse(df = 4)
  result_t <- gen.aruma.tse(
    n = n, phi = phi, innov_gen = t_gen, sn = 42, plot = FALSE
  )
  
  # Laplace innovations
  lap_gen <- make.gen.laplace.tse()
  result_lap <- gen.aruma.tse(
    n = n, phi = phi, innov_gen = lap_gen, sn = 42, plot = FALSE
  )
  
  # Compare kurtosis
  kurt_norm <- excess_kurtosis(result_norm$y)
  kurt_t <- excess_kurtosis(result_t$y)
  kurt_lap <- excess_kurtosis(result_lap$y)
  
  # t(4) and Laplace should produce heavier tails than normal
  expect_gt(kurt_t, kurt_norm)
  expect_gt(kurt_lap, kurt_norm)
})


# ==============================================================================
# Full Workflow: Seasonal ARUMA
# ==============================================================================

test_that("workflow: seasonal ARUMA with various innovations", {
  n <- 300
  s <- 12
  
  # With normal innovations
  result_norm <- gen.aruma.tse(
    n = n, phi = 0.5, s = s,
    innov_gen = make.gen.norm.tse(),
    sn = 42, plot = FALSE
  )
  expect_length(result_norm$y, n)
  expect_equal(result_norm$s, s)
  
  # With mixture innovations
  mix_gen <- make.gen.mixnorm.tse(sd1 = 1, sd2 = 4, prob1 = 0.9)
  result_mix <- gen.aruma.tse(
    n = n, phi = 0.5, s = s,
    innov_gen = mix_gen,
    sn = 42, plot = FALSE
  )
  expect_length(result_mix$y, n)
})


# ==============================================================================
# Full Workflow: Complex ARIMA Model
# ==============================================================================

test_that("workflow: ARIMA(2,1,1) with scaled t innovations", {
  t_gen <- make.gen.t.tse(df = 6, scale = TRUE)
  
  result <- gen.aruma.tse(
    n = 400,
    phi = c(0.5, -0.2),
    theta = 0.3,
    d = 1,
    innov_gen = t_gen,
    sn = 42,
    plot = FALSE
  )
  
  expect_s3_class(result, "aruma")
  expect_equal(result$p, 2)
  expect_equal(result$q, 1)
  expect_equal(result$d, 1)
  expect_length(result$y, 400)
})


# ==============================================================================
# Full Workflow: Reproducibility Across Components
# ==============================================================================

test_that("workflow: full pipeline is reproducible with seeds", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  # First run
  set.seed(123)
  garch_gen1 <- make.gen.garch.tse(omega = 0.1, alpha = 0.1, beta = 0.8)
  y1 <- garch_gen1(500)
  
  aruma1 <- gen.aruma.tse(
    n = 300, phi = 0.5, innov_gen = make.gen.norm.tse(),
    sn = 456, plot = FALSE
  )
  
  # Second run with same seeds
  
  set.seed(123)
  garch_gen2 <- make.gen.garch.tse(omega = 0.1, alpha = 0.1, beta = 0.8)
  y2 <- garch_gen2(500)
  
  aruma2 <- gen.aruma.tse(
    n = 300, phi = 0.5, innov_gen = make.gen.norm.tse(),
    sn = 456, plot = FALSE
  )
  
  # Verify reproducibility
  expect_equal(y1, y2)
  expect_equal(aruma1$y, aruma2$y)
})


# ==============================================================================
# Full Workflow: Compare Multiple GARCH Specifications
# ==============================================================================

test_that("workflow: GARCH comparison identifies correct model", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("WeightedPortTest")
  
  # Generate data from known GARCH(1,1) process
  garch_gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.15, beta = 0.75)
  set.seed(TEST_SEED)
  y <- garch_gen(1500)
  
  # Fit models
  results <- suppressWarnings(
    compare.garch.tse(y, arch.range = 0:2, garch.range = 0:2)
  )
  
  # GARCH(1,1) should be among the best models by BIC
  best_bic_model <- results$comparison$Model[which.min(results$comparison$BIC)]
  
  # Either GARCH(1,1) or a close variant should win
  # (Due to sampling variability, we check it's in top 3)
  sorted_by_bic <- results$comparison[order(results$comparison$BIC), ]
  top_models <- sorted_by_bic$Model[1:3]
  
  expect_true(
    "GARCH(1,1)" %in% top_models,
    info = paste("Top models by BIC:", paste(top_models, collapse = ", "))
  )
})


# ==============================================================================
# Full Workflow: Extended Generators in ARUMA
# ==============================================================================

test_that("workflow: skew-t innovations through ARMA filter", {
  skip_if_not_installed("sn")
  
  skt_gen <- make.gen.skt.tse(df = 6, alpha = -2, scale = TRUE)
  
  result <- gen.aruma.tse(
    n = 500,
    phi = 0.6,
    theta = 0.2,
    innov_gen = skt_gen,
    sn = 42,
    plot = FALSE
  )
  
  expect_s3_class(result, "aruma")
  expect_length(result$y, 500)
  
  # Output should show some left skewness
  skew <- sample_skewness(result$y)
  expect_lt(skew, 0.5)  # Should be less right-skewed than symmetric
})

test_that("workflow: GED innovations through ARMA filter", {
  skip_if_not_installed("fGarch")
  
  ged_gen <- make.gen.ged.tse(nu = 1.5, sd = 1)
  
  result <- gen.aruma.tse(
    n = 500,
    phi = 0.5,
    innov_gen = ged_gen,
    sn = 42,
    plot = FALSE
  )
  
  expect_s3_class(result, "aruma")
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
  
  # Heavy-tailed GED should produce excess kurtosis
  kurt <- excess_kurtosis(result$y)
  expect_gt(kurt, 0)
})


# ==============================================================================
# Error Handling Across Components
# ==============================================================================

test_that("workflow: graceful handling of convergence issues",
          {
            skip_if_not_installed("rugarch")
            skip_if_not_installed("dplyr")
            skip_if_not_installed("WeightedPortTest")
            
            # Near-degenerate data that may cause fitting issues
            set.seed(999)
            y <- c(rnorm(50), rep(0, 10), rnorm(50))
            
            # Should not error, but may have warnings/messages
            result <- suppressMessages(suppressWarnings(
              compare.garch.tse(y, arch.range = 1, garch.range = 0:1)
            ))
            
            # Should return something usable even if some models failed
            expect_s3_class(result, "garch.comparison.tse")
          })