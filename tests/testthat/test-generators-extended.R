# test-generators-extended.R
# Tests for innovation generators requiring external packages
# Covers: make.gen.skt.tse (sn), make.gen.ged.tse (fGarch),
#         make.gen.garch.tse (rugarch)

library(testthat)

# Source test helpers
#source(test_path("helper-test-utils.R"))

# ==============================================================================
# make.gen.skt.tse (requires: sn)
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

test_that("make.gen.skt.tse has mean approximately zero (symmetric)", {
  skip_if_not_installed("sn")
  gen_sym <- make.gen.skt.tse(df = 5, alpha = 0)
  set.seed(TEST_SEED)
  x_sym <- gen_sym(N_LARGE)
  expect_true(check_mean_zero(x_sym))
})

test_that("make.gen.skt.tse has mean approximately zero (skewed)", {
  skip_if_not_installed("sn")
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
  expected_var <- 7 / (7 - 2)
  expect_true(check_variance(x, expected_var = expected_var, tol = 0.15))
})

test_that("make.gen.skt.tse produces skewed output when alpha != 0", {
  skip_if_not_installed("sn")
  
  gen_right <- make.gen.skt.tse(df = 10, alpha = 5, scale = FALSE)
  set.seed(TEST_SEED)
  x_right <- gen_right(N_LARGE)
  skew_right <- sample_skewness(x_right)
  
  gen_left <- make.gen.skt.tse(df = 10, alpha = -5, scale = FALSE)
  set.seed(TEST_SEED + 1)
  x_left <- gen_left(N_LARGE)
  skew_left <- sample_skewness(x_left)
  
  expect_gt(skew_right, 0)
  expect_lt(skew_left, 0)
})

test_that("make.gen.skt.tse validates df", {
  skip_if_not_installed("sn")
  expect_error(make.gen.skt.tse(df = 0), "df.*positive")
})

test_that("make.gen.skt.tse scale=TRUE gives unit variance", {
  skip_if_not_installed("sn")
  gen <- make.gen.skt.tse(df = 7, alpha = 2, scale = TRUE)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_variance(x, expected_var = 1, tol = 0.15))
})


# ==============================================================================
# make.gen.ged.tse (requires: fGarch)
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
  expect_gt(qq_correlation(x), 0.999)
})

test_that("make.gen.ged.tse nu<2 has heavy tails", {
  skip_if_not_installed("fGarch")
  gen <- make.gen.ged.tse(nu = 1.2)
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  kurt <- excess_kurtosis(x)
  expect_gt(kurt, 0.5)
})

test_that("make.gen.ged.tse validates parameters", {
  skip_if_not_installed("fGarch")
  expect_error(make.gen.ged.tse(nu = 0), "nu must be positive")
  expect_error(make.gen.ged.tse(nu = 2, sd = -1), "sd must be positive")
})


# ==============================================================================
# make.gen.garch.tse (requires: rugarch)
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
  acf_sq <- acf(x^2, lag.max = 5, plot = FALSE)$acf[-1]
  expect_gt(acf_sq[1], 0.05)
})

test_that("make.gen.garch.tse validates parameters", {
  skip_if_not_installed("rugarch")
  expect_error(
    make.gen.garch.tse(omega = -0.1, alpha = 0.1),
    "omega must be positive"
  )
  expect_error(
    make.gen.garch.tse(omega = 0.1, alpha = -0.1),
    "alpha coefficients must be non-negative"
  )
  expect_error(
    make.gen.garch.tse(omega = 0.1, alpha = 0.1, beta = -0.5),
    "beta coefficients must be non-negative"
  )
})

test_that("make.gen.garch.tse warns for non-stationary params", {
  skip_if_not_installed("rugarch")
  expect_warning(
    make.gen.garch.tse(omega = 0.1, alpha = 0.5, beta = 0.6),
    "non-stationary"
  )
})

test_that("make.gen.garch.tse works with ARCH only (no beta)", {
  skip_if_not_installed("rugarch")
  gen <- make.gen.garch.tse(omega = 0.2, alpha = c(0.3, 0.2))
  x <- gen(500)
  expect_length(x, 500)
  expect_true(all(is.finite(x)))
})

test_that("make.gen.garch.tse works with t-distribution", {
  skip_if_not_installed("rugarch")
  gen <- make.gen.garch.tse(
    omega = 0.1, alpha = 0.1, beta = 0.85,
    distribution = "std",
    distribution.params = list(shape = 5)
  )
  x <- gen(500)
  expect_length(x, 500)
  expect_true(all(is.finite(x)))
})