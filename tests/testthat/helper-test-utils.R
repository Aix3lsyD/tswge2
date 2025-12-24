# helper-test-utils.R
# Shared test utilities, constants, and helper functions
# testthat automatically sources helper-*.R files before running tests

# ==============================================================================
# Test Constants
# ==============================================================================

# Sample sizes
N_SMALL <- 500L    # Smaller sample for faster tests
N_LARGE <- 10000L  # Larger sample for distributional tests

# Tolerances
TOL_MEAN <- 0.05   # Tolerance for mean (should be close to 0
TOL_VAR <- 0.10    # Tolerance for variance (proportion of theoretical)

# Seeds
TEST_SEED <- 42L


# ==============================================================================
# Statistical Helper Functions
# ==============================================================================

#' Check if sample mean is approximately zero
#' @param x Numeric vector
#' @param tol Tolerance (default TOL_MEAN)
#' @return Logical
check_mean_zero <- function(x, tol = TOL_MEAN) {
  
  abs(mean(x)) < tol
}

#' Check if sample variance is close to expected
#' @param x Numeric vector
#' @param expected_var Expected variance
#' @param tol Tolerance as proportion of expected (default TOL_VAR)
#' @return Logical
check_variance <- function(x, expected_var, tol = TOL_VAR) {
  abs(var(x) - expected_var) / expected_var < tol
}

#' Calculate excess kurtosis
#' @param x Numeric vector
#' @return Numeric scalar
excess_kurtosis <- function(x) {
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  sum((x - m)^4) / (n * s^4) - 3
}

#' Calculate sample skewness
#' @param x Numeric vector
#' @return Numeric scalar
sample_skewness <- function(x) {
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  (n / ((n - 1) * (n - 2))) * sum(((x - m) / s)^3)
}

#' QQ correlation for normality check
#' More stable than Shapiro-Wilk for large samples
#' @param x Numeric vector
#' @return Numeric scalar (correlation coefficient)
qq_correlation <- function(x) {
  cor(sort(x), qnorm(ppoints(length(x))))
}


# ==============================================================================
# Data Generation Helpers
# ==============================================================================
#' Generate GARCH(1,1) data for testing
#' Uses fixed parameters known to produce well-behaved data
#' @param n Sample size
#' @param seed Random seed (default TEST_SEED)
#' @return Numeric vector
generate_garch_data <- function(n, seed = TEST_SEED) {
  garch_gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.15, beta = 0.75)
  set.seed(seed)
  garch_gen(n)
}


# ==============================================================================
# Mock Data for Display Tests
# ==============================================================================

#' Create mock garch.comparison.tse object for display testing
#' Avoids fitting actual models when testing display functions
#' @param with_na Logical; include NA values for edge case testing
#' @return Object of class garch.comparison.tse
create_mock_garch_comparison <- function(with_na = FALSE) {
  if (with_na) {
    comparison <- dplyr::tibble(
      Model = c("ARCH(1)", "GARCH(1,1)"),
      ARCH = c(1L, 1L),
      GARCH = c(0L, 1L),
      AIC = c(100, 95),
      AICC = c(101, 96),
      BIC = c(105, 100),
      WLB1 = c(0.5, NA_real_),
      WLB2 = c(NA_real_, 0.3),
      WLB3 = c(0.2, 0.4),
      Nyblom = c(0.5, NA_real_),
      Nyblom_crit = c(1.0, 1.0),
      SignBias = c(NA_real_, 0.1),
      n_sig = c(2L, 3L),
      n_coef = c(2L, 3L)
    )
  } else {
    comparison <- dplyr::tibble(
      Model = c("ARCH(1)", "GARCH(1,1)"),
      ARCH = c(1L, 1L),
      GARCH = c(0L, 1L),
      AIC = c(100, 95),
      AICC = c(101, 96),
      BIC = c(105, 100),
      WLB1 = c(0.5, 0.6),
      WLB2 = c(0.3, 0.4),
      WLB3 = c(0.2, 0.4),
      Nyblom = c(0.5, 0.3),
      Nyblom_crit = c(1.0, 1.0),
      SignBias = c(0.2, 0.1),
      n_sig = c(2L, 3L),
      n_coef = c(2L, 3L)
    )
  }
  
  structure(
    list(
      fits = list(),
      comparison = comparison,
      distribution = "norm",
      n = 500L
    ),
    class = "garch.comparison.tse"
  )
}

#' Create empty garch.comparison.tse for edge case testing
#' @return Object of class garch.comparison.tse with empty comparison
create_empty_garch_comparison <- function() {
  structure(
    list(
      fits = list(),
      comparison = data.frame(),
      distribution = "norm",
      n = 100L
    ),
    class = "garch.comparison.tse"
  )
}