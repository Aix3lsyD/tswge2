# test-generators.R
# Tests for innovation generators and gen.aruma.tse

library(testthat)

# ==============================================================================
# Test Parameters
# ==============================================================================

# Large sample for distributional tests
N_LARGE <- 10000

# Tolerance for mean (should be close to 0)
TOL_MEAN <- 0.05

# Tolerance for variance (proportion of theoretical)
TOL_VAR <- 0.10

# Seed for reproducibility tests
TEST_SEED <- 42


# ==============================================================================
# Helper Functions
# ==============================================================================

#' Check if sample mean is approximately zero
check_mean_zero <- function(x, tol = TOL_MEAN) {
  abs(mean(x)) < tol
}

#' Check if sample variance is close to expected
check_variance <- function(x, expected_var, tol = TOL_VAR) {
  abs(var(x) - expected_var) / expected_var < tol
}

#' Calculate excess kurtosis
excess_kurtosis <- function(x) {
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  sum((x - m)^4) / (n * s^4) - 3
}

#' Calculate sample skewness
sample_skewness <- function(x) {
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  (n / ((n - 1) * (n - 2))) * sum(((x - m) / s)^3)
}

#' QQ correlation for normality check (more stable than Shapiro-Wilk)
qq_correlation <- function(x) {
  cor(sort(x), qnorm(ppoints(length(x))))
}


# ==============================================================================
# Tests: make.gen.norm.tse
# ==============================================================================

test_that("make.gen.norm.tse returns function", {
  gen <- make.gen.norm.tse()
  expect_type(gen, "closure")
})

test_that("make.gen.norm.tse generates correct length", {
  gen <- make.gen.norm.tse()
  expect_length(gen(100), 100)
  expect_length(gen(500), 500)
})

test_that("make.gen.norm.tse has mean zero", {
  gen <- make.gen.norm.tse()
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_mean_zero(x))
})

test_that("make.gen.norm.tse has correct variance", {
  # sd = 1 (default)
  gen1 <- make.gen.norm.tse(sd = 1)
  set.seed(TEST_SEED)
  x1 <- gen1(N_LARGE)
  expect_true(check_variance(x1, expected_var = 1))
  
  # sd = 2 -> variance = 4
  gen2 <- make.gen.norm.tse(sd = 2)
  set.seed(TEST_SEED)
  x2 <- gen2(N_LARGE)
  expect_true(check_variance(x2, expected_var = 4))
})

test_that("make.gen.norm.tse is approximately normal", {
  gen <- make.gen.norm.tse()
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  
  # QQ correlation should be very high for normal data
  expect_gt(qq_correlation(x), 0.999)
})

test_that("make.gen.norm.tse validates sd", {
  expect_error(make.gen.norm.tse(sd = 0), "sd must be positive")
  expect_error(make.gen.norm.tse(sd = -1), "sd must be positive")
})

test_that("make.gen.norm.tse is reproducible with set.seed", {
  gen <- make.gen.norm.tse()
  
  set.seed(TEST_SEED)
  x1 <- gen(100)
  
  set.seed(TEST_SEED)
  x2 <- gen(100)
  
  expect_equal(x1, x2)
})


# ==============================================================================
# Tests: make.gen.t.tse
# ==============================================================================

test_that("make.gen.t.tse returns function", {
  gen <- make.gen.t.tse(df = 5)
  expect_type(gen, "closure")
})

test_that("make.gen.t.tse generates correct length", {
  gen <- make.gen.t.tse(df = 5)
  expect_length(gen(100), 100)
})

test_that("make.gen.t.tse has mean zero", {
  gen <- make.gen.t.tse(df = 5)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_mean_zero(x))
})

test_that("make.gen.t.tse has correct unscaled variance", {
  # t(df) has variance df/(df-2)
  df <- 7
  expected_var <- df / (df - 2)  # 7/5 = 1.4
  
  gen <- make.gen.t.tse(df = df, scale = FALSE)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  
  expect_true(check_variance(x, expected_var = expected_var))
})

test_that("make.gen.t.tse scale=TRUE gives unit variance", {
  gen <- make.gen.t.tse(df = 7, scale = TRUE)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  
  expect_true(check_variance(x, expected_var = 1))
})

test_that("make.gen.t.tse has heavy tails", {
  gen <- make.gen.t.tse(df = 5)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  
  # t(5) has excess kurtosis = 6/(df-4) = 6 for df=5
  # Should be notably positive (heavier than normal)
  kurt <- excess_kurtosis(x)
  expect_gt(kurt, 1)
})

test_that("make.gen.t.tse validates df", {
  expect_error(make.gen.t.tse(df = 0), "df must be positive")
  expect_error(make.gen.t.tse(df = -5), "df must be positive")
})

test_that("make.gen.t.tse warns when df <= 2 with scale", {
  expect_warning(make.gen.t.tse(df = 2, scale = TRUE), "infinite variance")
})


# ==============================================================================
# Tests: make.gen.skt.tse
# ==============================================================================

test_that("make.gen.skt.tse returns function", {
  skip_if_not_installed("sn")
  gen <- make.gen.skt.tse(df = 5)
  expect_type(gen, "closure")
})

test_that("make.gen.skt.tse generates correct length", {
  skip_if_not_installed("sn")
  gen <- make.gen.skt.tse(df = 5)
  expect_length(gen(100), 100)
})

test_that("make.gen.skt.tse has mean approximately zero", {
  skip_if_not_installed("sn")
  
  # Symmetric case (alpha = 0)
  gen_sym <- make.gen.skt.tse(df = 5, alpha = 0)
  set.seed(TEST_SEED)
  x_sym <- gen_sym(N_LARGE)
  expect_true(check_mean_zero(x_sym))
  
  
  # Skewed case (alpha != 0) - should still be centered
  gen_skew <- make.gen.skt.tse(df = 5, alpha = 2)
  set.seed(TEST_SEED)
  x_skew <- gen_skew(N_LARGE)
  expect_true(check_mean_zero(x_skew, tol = 0.1))
})

test_that("make.gen.skt.tse alpha=0 behaves like t", {
  skip_if_not_installed("sn")
  
  gen <- make.gen.skt.tse(df = 7, alpha = 0)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  
  # Should have similar variance to t(7)
  expected_var <- 7 / (7 - 2)
  expect_true(check_variance(x, expected_var = expected_var, tol = 0.15))
})

test_that("make.gen.skt.tse produces skewed output when alpha != 0", {
  skip_if_not_installed("sn")
  
  # Right skew
  gen_right <- make.gen.skt.tse(df = 10, alpha = 5, scale = FALSE)
  set.seed(TEST_SEED)
  x_right <- gen_right(N_LARGE)
  skew_right <- sample_skewness(x_right)
  
  # Left skew
  gen_left <- make.gen.skt.tse(df = 10, alpha = -5, scale = FALSE)
  set.seed(TEST_SEED + 1)
  x_left <- gen_left(N_LARGE)
  skew_left <- sample_skewness(x_left)
  
  # Right skew should be positive, left skew should be negative
  expect_gt(skew_right, 0)
  expect_lt(skew_left, 0)
})

test_that("make.gen.skt.tse validates df", {
  skip_if_not_installed("sn")
  expect_error(make.gen.skt.tse(df = 0), "df must be positive")
})


# ==============================================================================
# Tests: make.gen.ged.tse
# ==============================================================================

test_that("make.gen.ged.tse returns function", {
  skip_if_not_installed("fGarch")
  gen <- make.gen.ged.tse(nu = 2)
  expect_type(gen, "closure")
})

test_that("make.gen.ged.tse generates correct length", {
  skip_if_not_installed("fGarch")
  gen <- make.gen.ged.tse(nu = 2)
  expect_length(gen(100), 100)
})

test_that("make.gen.ged.tse has mean zero", {
  skip_if_not_installed("fGarch")
  gen <- make.gen.ged.tse(nu = 1.5)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_mean_zero(x))
})

test_that("make.gen.ged.tse nu=2 is approximately normal", {
  skip_if_not_installed("fGarch")
  
  gen <- make.gen.ged.tse(nu = 2, sd = 1)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  
  # QQ correlation should be very high for normal data
  expect_gt(qq_correlation(x), 0.999)
})

test_that("make.gen.ged.tse nu<2 has heavy tails", {
  skip_if_not_installed("fGarch")
  
  gen <- make.gen.ged.tse(nu = 1.2)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  
  kurt <- excess_kurtosis(x)
  expect_gt(kurt, 0.5)  # Heavier than normal
})

test_that("make.gen.ged.tse validates parameters", {
  skip_if_not_installed("fGarch")
  expect_error(make.gen.ged.tse(nu = 0), "nu must be positive")
  expect_error(make.gen.ged.tse(nu = 2, sd = -1), "sd must be positive")
})


# ==============================================================================
# Tests: make.gen.laplace.tse
# ==============================================================================

test_that("make.gen.laplace.tse returns function", {
  gen <- make.gen.laplace.tse()
  expect_type(gen, "closure")
})

test_that("make.gen.laplace.tse generates correct length", {
  gen <- make.gen.laplace.tse()
  expect_length(gen(100), 100)
})

test_that("make.gen.laplace.tse has mean zero", {
  gen <- make.gen.laplace.tse()
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_mean_zero(x))
})

test_that("make.gen.laplace.tse default has unit variance", {
  gen <- make.gen.laplace.tse()  # default scale = 1/sqrt(2)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_variance(x, expected_var = 1))
})

test_that("make.gen.laplace.tse has correct variance for custom scale", {
  # Variance = 2 * scale^2
  scale <- 2
  expected_var <- 2 * scale^2  # = 8
  
  gen <- make.gen.laplace.tse(scale = scale)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_variance(x, expected_var = expected_var))
})

test_that("make.gen.laplace.tse has excess kurtosis ~3", {
  gen <- make.gen.laplace.tse()
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  
  kurt <- excess_kurtosis(x)
  # Laplace has theoretical excess kurtosis = 3
  expect_true(abs(kurt - 3) < 0.5)
})

test_that("make.gen.laplace.tse validates scale", {
  expect_error(make.gen.laplace.tse(scale = 0), "scale must be positive")
  expect_error(make.gen.laplace.tse(scale = -1), "scale must be positive")
})

test_that("make.gen.laplace.tse handles edge cases without Inf", {
  # Test that we don't get Inf values from log(0)
  gen <- make.gen.laplace.tse()
  set.seed(TEST_SEED)
  
  # Generate many samples to stress test edge cases
  x <- gen(100000)
  
  expect_true(all(is.finite(x)))
})


# ==============================================================================
# Tests: make.gen.unif.tse
# ==============================================================================

test_that("make.gen.unif.tse returns function", {
  gen <- make.gen.unif.tse()
  expect_type(gen, "closure")
})

test_that("make.gen.unif.tse generates correct length", {
  gen <- make.gen.unif.tse()
  expect_length(gen(100), 100)
})

test_that("make.gen.unif.tse has mean zero", {
  gen <- make.gen.unif.tse()
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_mean_zero(x))
})

test_that("make.gen.unif.tse default has unit variance", {
  gen <- make.gen.unif.tse()  # default half_width = sqrt(3)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_variance(x, expected_var = 1))
})

test_that("make.gen.unif.tse respects bounds", {
  half_width <- 2
  gen <- make.gen.unif.tse(half_width = half_width)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  
  expect_true(all(x >= -half_width))
  expect_true(all(x <= half_width))
})

test_that("make.gen.unif.tse has negative kurtosis", {
  gen <- make.gen.unif.tse()
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  
  kurt <- excess_kurtosis(x)
  # Uniform has theoretical excess kurtosis = -1.2
  expect_lt(kurt, 0)
})

test_that("make.gen.unif.tse validates half_width", {
  expect_error(make.gen.unif.tse(half_width = 0), "half_width must be positive")
})


# ==============================================================================
# Tests: make.gen.mixnorm.tse
# ==============================================================================

test_that("make.gen.mixnorm.tse returns function", {
  gen <- make.gen.mixnorm.tse()
  expect_type(gen, "closure")
})

test_that("make.gen.mixnorm.tse generates correct length", {
  gen <- make.gen.mixnorm.tse()
  expect_length(gen(100), 100)
})

test_that("make.gen.mixnorm.tse has mean zero", {
  gen <- make.gen.mixnorm.tse()
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_mean_zero(x))
})

test_that("make.gen.mixnorm.tse has correct variance", {
  sd1 <- 1
  sd2 <- 3
  prob1 <- 0.9
  # Variance = prob1 * sd1^2 + (1-prob1) * sd2^2
  expected_var <- prob1 * sd1^2 + (1 - prob1) * sd2^2  # 0.9 + 0.9 = 1.8
  
  gen <- make.gen.mixnorm.tse(sd1 = sd1, sd2 = sd2, prob1 = prob1)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_variance(x, expected_var = expected_var))
})

test_that("make.gen.mixnorm.tse has heavy tails from mixture", {
  gen <- make.gen.mixnorm.tse(sd1 = 1, sd2 = 5, prob1 = 0.9)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  
  kurt <- excess_kurtosis(x)
  expect_gt(kurt, 0)  # Heavier than normal
})

test_that("make.gen.mixnorm.tse validates parameters", {
  expect_error(make.gen.mixnorm.tse(sd1 = 0), "sd1 and sd2 must be positive")
  expect_error(make.gen.mixnorm.tse(sd2 = -1), "sd1 and sd2 must be positive")
  expect_error(make.gen.mixnorm.tse(prob1 = 0), "prob1 must be between 0 and 1")
  expect_error(make.gen.mixnorm.tse(prob1 = 1), "prob1 must be between 0 and 1")
})


# ==============================================================================
# Tests: make.gen.garch.tse
# ==============================================================================

test_that("make.gen.garch.tse returns function", {
  skip_if_not_installed("rugarch")
  gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.1, beta = 0.8)
  expect_type(gen, "closure")
})

test_that("make.gen.garch.tse generates correct length", {
  skip_if_not_installed("rugarch")
  gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.1, beta = 0.8)
  expect_length(gen(100), 100)
  expect_length(gen(500), 500)
})

test_that("make.gen.garch.tse has mean approximately zero", {
  skip_if_not_installed("rugarch")
  gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.1, beta = 0.8)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_mean_zero(x, tol = 0.1))
})

test_that("make.gen.garch.tse exhibits volatility clustering", {
  skip_if_not_installed("rugarch")
  
  gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.15, beta = 0.8)
  set.seed(TEST_SEED)
  x <- gen(2000)
  
  # Squared returns should be autocorrelated
  acf_sq <- acf(x^2, lag.max = 5, plot = FALSE)$acf[-1]
  
  # At least first lag should be notably positive
  expect_gt(acf_sq[1], 0.05)
})

test_that("make.gen.garch.tse validates parameters", {
  skip_if_not_installed("rugarch")
  
  expect_error(make.gen.garch.tse(omega = -0.1, alpha = 0.1), 
               "omega must be positive")
  expect_error(make.gen.garch.tse(omega = 0.1, alpha = -0.1), 
               "alpha coefficients must be non-negative")
  expect_error(make.gen.garch.tse(omega = 0.1, alpha = 0.1, beta = -0.5), 
               "beta coefficients must be non-negative")
})

test_that("make.gen.garch.tse warns for non-stationary params", {
  skip_if_not_installed("rugarch")
  
  expect_warning(
    make.gen.garch.tse(omega = 0.1, alpha = 0.5, beta = 0.6),
    "non-stationary"
  )
})


# ==============================================================================
# Tests: Edge Cases for All Generators
# ==============================================================================

test_that("generators handle n = 1", {
  expect_length(make.gen.norm.tse()(1), 1)
  expect_length(make.gen.t.tse(df = 5)(1), 1)
  expect_length(make.gen.laplace.tse()(1), 1)
  expect_length(make.gen.unif.tse()(1), 1)
  expect_length(make.gen.mixnorm.tse()(1), 1)
})

test_that("generators handle large n without error", {
  gen <- make.gen.norm.tse()
  x <- gen(100000)
  expect_length(x, 100000)
  expect_true(all(is.finite(x)))
})

test_that("generator closures capture parameters correctly", {
  # Create two generators with different params
  gen1 <- make.gen.norm.tse(sd = 1)
  gen2 <- make.gen.norm.tse(sd = 10)
  
  set.seed(TEST_SEED)
  x1 <- gen1(N_LARGE)
  
  set.seed(TEST_SEED)
  x2 <- gen2(N_LARGE)
  
  # Same seed, so x2 should be exactly 10 * x1
  expect_equal(x2, x1 * 10)
})

test_that("t generator closure captures df correctly", {
  gen1 <- make.gen.t.tse(df = 5)
  gen2 <- make.gen.t.tse(df = 30)
  
  set.seed(TEST_SEED)
  x1 <- gen1(N_LARGE)
  
  set.seed(TEST_SEED)
  x2 <- gen2(N_LARGE)
  
  # df=5 should have heavier tails (higher kurtosis) than df=30
  expect_gt(excess_kurtosis(x1), excess_kurtosis(x2))
})


# ==============================================================================
# Tests: gen.aruma.tse - Basic Functionality
# ==============================================================================

test_that("gen.aruma.tse returns correct structure", {
  result <- gen.aruma.tse(n = 100, phi = 0.5, plot = FALSE)
  
  expect_s3_class(result, "aruma")
  expect_true("y" %in% names(result))
  expect_true("n" %in% names(result))
  expect_true("p" %in% names(result))
  expect_true("phi" %in% names(result))
  expect_true("plot" %in% names(result))
})

test_that("gen.aruma.tse generates correct length", {
  result <- gen.aruma.tse(n = 200, phi = 0.5, plot = FALSE)
  expect_length(result$y, 200)
  
  result2 <- gen.aruma.tse(n = 500, phi = 0.5, plot = FALSE)
  expect_length(result2$y, 500)
})

test_that("gen.aruma.tse is reproducible with sn parameter", {
  result1 <- gen.aruma.tse(n = 100, phi = 0.7, sn = 123, plot = FALSE)
  result2 <- gen.aruma.tse(n = 100, phi = 0.7, sn = 123, plot = FALSE)
  
  expect_equal(result1$y, result2$y)
})

test_that("gen.aruma.tse different seeds give different results", {
  result1 <- gen.aruma.tse(n = 100, phi = 0.7, sn = 123, plot = FALSE)
  result2 <- gen.aruma.tse(n = 100, phi = 0.7, sn = 456, plot = FALSE)
  
  expect_false(all(result1$y == result2$y))
})

test_that("gen.aruma.tse sn = NULL uses current RNG state", {
  # Two calls with different external seeds should differ
  set.seed(111)
  result1 <- gen.aruma.tse(n = 50, phi = 0.5, sn = NULL, plot = FALSE)
  
  set.seed(222)
  result2 <- gen.aruma.tse(n = 50, phi = 0.5, sn = NULL, plot = FALSE)
  
  expect_false(all(result1$y == result2$y))
})

test_that("gen.aruma.tse AR(1) has correct ACF structure", {
  phi <- 0.8
  result <- gen.aruma.tse(n = 2000, phi = phi, sn = 42, plot = FALSE)
  
  # Theoretical ACF at lag 1 for AR(1) is phi
  acf_vals <- acf(result$y, lag.max = 3, plot = FALSE)$acf
  acf_lag1 <- acf_vals[2]  # Index 1 is lag 0
  
  # Should be close to phi (within sampling variability)
  expect_true(abs(acf_lag1 - phi) < 0.1)
})


# ==============================================================================
# Tests: gen.aruma.tse - Custom Innovation Generators
# ==============================================================================

test_that("gen.aruma.tse uses custom innovation generator", {
  # Use uniform innovations (bounded)
  unif_gen <- make.gen.unif.tse(half_width = 1)
  
  result <- gen.aruma.tse(n = 500, phi = 0.3, innov_gen = unif_gen, 
                          sn = 42, plot = FALSE)
  
  # Series should exist and have finite values
  expect_true(all(is.finite(result$y)))
  expect_length(result$y, 500)
})

test_that("gen.aruma.tse with t-distributed innovations", {
  t_gen <- make.gen.t.tse(df = 5)
  result <- gen.aruma.tse(n = 500, phi = 0.5, innov_gen = t_gen, 
                          sn = 42, plot = FALSE)
  
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse with Laplace innovations", {
  lap_gen <- make.gen.laplace.tse()
  result <- gen.aruma.tse(n = 500, phi = 0.5, innov_gen = lap_gen, 
                          sn = 42, plot = FALSE)
  
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse with mixture normal innovations", {
  mix_gen <- make.gen.mixnorm.tse(sd1 = 1, sd2 = 5, prob1 = 0.9)
  result <- gen.aruma.tse(n = 500, phi = 0.5, innov_gen = mix_gen, 
                          sn = 42, plot = FALSE)
  
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse with GARCH innovations", {
  skip_if_not_installed("rugarch")
  
  garch_gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.1, beta = 0.8)
  result <- gen.aruma.tse(n = 500, phi = 0.5, innov_gen = garch_gen, 
                          sn = 42, plot = FALSE)
  
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse with skew-t innovations", {
  skip_if_not_installed("sn")
  
  skt_gen <- make.gen.skt.tse(df = 5, alpha = -2)
  result <- gen.aruma.tse(n = 500, phi = 0.5, innov_gen = skt_gen, 
                          sn = 42, plot = FALSE)
  
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
})


# ==============================================================================
# Tests: gen.aruma.tse - MA Components
# ==============================================================================

test_that("gen.aruma.tse handles MA(1)", {
  result <- gen.aruma.tse(n = 200, theta = 0.6, sn = 42, plot = FALSE)
  
  expect_equal(result$q, 1)
  expect_length(result$y, 200)
})

test_that("gen.aruma.tse handles MA(2)", {
  result <- gen.aruma.tse(n = 300, theta = c(0.4, 0.2), sn = 42, plot = FALSE)
  
  expect_equal(result$q, 2)
  expect_length(result$y, 300)
})


# ==============================================================================
# Tests: gen.aruma.tse - ARMA Components
# ==============================================================================

test_that("gen.aruma.tse handles ARMA(1,1)", {
  result <- gen.aruma.tse(n = 200, phi = 0.5, theta = 0.3, sn = 42, plot = FALSE)
  
  expect_equal(result$p, 1)
  expect_equal(result$q, 1)
  expect_length(result$y, 200)
})

test_that("gen.aruma.tse handles ARMA(2,1)", {
  result <- gen.aruma.tse(n = 300, phi = c(0.5, -0.3), theta = 0.2, 
                          sn = 42, plot = FALSE)
  
  expect_equal(result$p, 2)
  expect_equal(result$q, 1)
  expect_length(result$y, 300)
})

test_that("gen.aruma.tse handles ARMA(1,2)", {
  result <- gen.aruma.tse(n = 300, phi = 0.5, theta = c(0.3, 0.1), 
                          sn = 42, plot = FALSE)
  
  expect_equal(result$p, 1)
  expect_equal(result$q, 2)
  expect_length(result$y, 300)
})

test_that("gen.aruma.tse handles ARMA(2,2)", {
  result <- gen.aruma.tse(n = 300, phi = c(0.5, -0.2), theta = c(0.3, 0.1),
                          sn = 42, plot = FALSE)
  
  expect_equal(result$p, 2)
  expect_equal(result$q, 2)
  expect_length(result$y, 300)
})


# ==============================================================================
# Tests: gen.aruma.tse - Higher Order AR
# ==============================================================================

test_that("gen.aruma.tse handles AR(2)", {
  result <- gen.aruma.tse(n = 300, phi = c(0.5, -0.3), sn = 42, plot = FALSE)
  
  expect_equal(result$p, 2)
  expect_length(result$y, 300)
})

test_that("gen.aruma.tse handles AR(3)", {
  result <- gen.aruma.tse(n = 300, phi = c(0.4, -0.2, 0.1), sn = 42, plot = FALSE)
  
  expect_equal(result$p, 3)
  expect_length(result$y, 300)
})


# ==============================================================================
# Tests: gen.aruma.tse - Differencing
# ==============================================================================

test_that("gen.aruma.tse handles d = 1", {
  result <- gen.aruma.tse(n = 200, phi = 0.5, d = 1, sn = 42, plot = FALSE)
  
  expect_equal(result$d, 1)
  expect_length(result$y, 200)
})

test_that("gen.aruma.tse handles d = 2", {
  result <- gen.aruma.tse(n = 200, phi = 0.5, d = 2, sn = 42, plot = FALSE)
  
  expect_equal(result$d, 2)
  expect_length(result$y, 200)
})


# ==============================================================================
# Tests: gen.aruma.tse - Seasonal and Nonstationary Components
# ==============================================================================

test_that("gen.aruma.tse handles seasonal differencing (s > 0)", {
  result <- gen.aruma.tse(n = 200, phi = 0.5, s = 12, sn = 42, plot = FALSE)
  
  expect_equal(result$s, 12)
  expect_length(result$y, 200)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse handles s = 4 (quarterly)", {
  result <- gen.aruma.tse(n = 200, phi = 0.5, s = 4, sn = 42, plot = FALSE)
  
  expect_equal(result$s, 4)
  expect_length(result$y, 200)
})

test_that("gen.aruma.tse handles lambda parameter", {
  result <- gen.aruma.tse(n = 200, phi = 0.5, lambda = c(1.2, -0.5), 
                          sn = 42, plot = FALSE)
  
  expect_length(result$y, 200)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse handles lambda + s combined", {
  # This exercises the mult.wge pathway
  
  result <- gen.aruma.tse(n = 200, phi = 0.3, lambda = c(0.8), s = 4, 
                          sn = 42, plot = FALSE)
  
  expect_equal(result$s, 4)
  expect_length(result$y, 200)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse handles d + s combined", {
  result <- gen.aruma.tse(n = 200, phi = 0.3, d = 1, s = 12, 
                          sn = 42, plot = FALSE)
  
  expect_equal(result$d, 1)
  expect_equal(result$s, 12)
  expect_length(result$y, 200)
})


# ==============================================================================
# Tests: gen.aruma.tse - Near Unit Root (Adaptive Burn-in)
# ==============================================================================

test_that("gen.aruma.tse handles near-unit-root AR(1)", {
  # phi = 0.99 needs longer burn-in
  result <- gen.aruma.tse(n = 500, phi = 0.99, sn = 42, plot = FALSE)
  
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
  
  # Series should be stationary (not exploding to extreme values)
  expect_true(var(result$y) < 1e6)
})

test_that("gen.aruma.tse handles near-unit-root AR(2)", {
  # AR(2) with roots near unit circle
  result <- gen.aruma.tse(n = 500, phi = c(0.5, 0.49), sn = 42, plot = FALSE)
  
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
})


# ==============================================================================
# Tests: gen.aruma.tse - Output Methods
# ==============================================================================

test_that("gen.aruma.tse print method works", {
  result <- gen.aruma.tse(n = 100, phi = 0.5, plot = FALSE)
  
  expect_output(print(result), "ARUMA Realization")
  expect_output(print(result), "n = 100")
})

test_that("gen.aruma.tse print method shows phi", {
  result <- gen.aruma.tse(n = 100, phi = 0.5, plot = FALSE)
  expect_output(print(result), "phi")
})

test_that("gen.aruma.tse print method shows theta", {
  result <- gen.aruma.tse(n = 100, theta = 0.3, plot = FALSE)
  expect_output(print(result), "theta")
})

test_that("gen.aruma.tse plot object is valid ggplot", {
  result <- gen.aruma.tse(n = 100, phi = 0.5, plot = FALSE, sn = 42)
  
  expect_s3_class(result$plot, "ggplot")
  
  # Should be able to build without error
  built <- ggplot2::ggplot_build(result$plot)
  expect_true(!is.null(built))
})

test_that("gen.aruma.tse plot method works", {
  result <- gen.aruma.tse(n = 100, phi = 0.5, plot = FALSE, sn = 42)
  
  # plot.aruma should return the ggplot invisibly
  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


# ==============================================================================
# Tests: gen.aruma.tse - Pure Noise (No AR/MA)
# ==============================================================================

test_that("gen.aruma.tse handles pure white noise", {
  result <- gen.aruma.tse(n = 200, sn = 42, plot = FALSE)
  
  expect_equal(result$p, 0)
  expect_equal(result$q, 0)
  expect_length(result$y, 200)
  
  # Should have approximately zero autocorrelation
  acf_vals <- acf(result$y, lag.max = 5, plot = FALSE)$acf[-1]
  expect_true(all(abs(acf_vals) < 0.15))
})

test_that("gen.aruma.tse handles pure noise with custom innovations", {
  t_gen <- make.gen.t.tse(df = 5)
  result <- gen.aruma.tse(n = 500, innov_gen = t_gen, sn = 42, plot = FALSE)
  
  expect_equal(result$p, 0)
  expect_equal(result$q, 0)
  expect_length(result$y, 500)
})


# ==============================================================================
# Tests: gen.aruma.tse - Integration with All Generator Types
# ==============================================================================

test_that("gen.aruma.tse works with all base generators", {
  n <- 200
  phi <- 0.5
  
  # Normal
  result_norm <- gen.aruma.tse(n = n, phi = phi, 
                               innov_gen = make.gen.norm.tse(),
                               sn = 42, plot = FALSE)
  expect_length(result_norm$y, n)
  
  # t
  result_t <- gen.aruma.tse(n = n, phi = phi, 
                            innov_gen = make.gen.t.tse(df = 5),
                            sn = 42, plot = FALSE)
  expect_length(result_t$y, n)
  
  # Laplace
  result_lap <- gen.aruma.tse(n = n, phi = phi, 
                              innov_gen = make.gen.laplace.tse(),
                              sn = 42, plot = FALSE)
  expect_length(result_lap$y, n)
  
  # Uniform
  result_unif <- gen.aruma.tse(n = n, phi = phi, 
                               innov_gen = make.gen.unif.tse(),
                               sn = 42, plot = FALSE)
  expect_length(result_unif$y, n)
  
  # Mixture normal
  result_mix <- gen.aruma.tse(n = n, phi = phi, 
                              innov_gen = make.gen.mixnorm.tse(),
                              sn = 42, plot = FALSE)
  expect_length(result_mix$y, n)
})


# ==============================================================================
# Tests: Distributional Properties Propagate Through ARMA
# ==============================================================================

test_that("heavy-tailed innovations produce heavier-tailed ARMA output", {
  n <- 5000
  phi <- 0.3  # Low persistence to let innovation distribution dominate
  
  # Normal innovations
  set.seed(TEST_SEED)
  result_norm <- gen.aruma.tse(n = n, phi = phi, 
                               innov_gen = make.gen.norm.tse(),
                               plot = FALSE)
  
  # t(3) innovations (very heavy tailed)
  set.seed(TEST_SEED)
  result_t <- gen.aruma.tse(n = n, phi = phi, 
                            innov_gen = make.gen.t.tse(df = 3),
                            plot = FALSE)
  
  # t(3) ARMA should have higher kurtosis
  kurt_norm <- excess_kurtosis(result_norm$y)
  kurt_t <- excess_kurtosis(result_t$y)
  
  expect_gt(kurt_t, kurt_norm)
})