# test-generators-base.R
# Tests for innovation generators requiring only base R
# Covers: make.gen.norm.tse, make.gen.t.tse, make.gen.laplace.tse,
#         make.gen.unif.tse, make.gen.mixnorm.tse

library(testthat)

# ==============================================================================
# make.gen.norm.tse
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
  gen1 <- make.gen.norm.tse(sd = 1)
  set.seed(TEST_SEED)
  x1 <- gen1(N_LARGE)
  expect_true(check_variance(x1, expected_var = 1))
  
  gen2 <- make.gen.norm.tse(sd = 2)
  set.seed(TEST_SEED)
  x2 <- gen2(N_LARGE)
  expect_true(check_variance(x2, expected_var = 4))
})

test_that("make.gen.norm.tse is approximately normal", {
  gen <- make.gen.norm.tse()
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
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
# make.gen.t.tse
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
  df <- 7
  expected_var <- df / (df - 2)
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

test_that("t generator closure captures df correctly", {
  gen1 <- make.gen.t.tse(df = 5)
  gen2 <- make.gen.t.tse(df = 30)
  set.seed(TEST_SEED)
  x1 <- gen1(N_LARGE)
  set.seed(TEST_SEED)
  x2 <- gen2(N_LARGE)
  expect_gt(excess_kurtosis(x1), excess_kurtosis(x2))
})


# ==============================================================================
# make.gen.laplace.tse
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
  gen <- make.gen.laplace.tse()
  set.seed(TEST_SEED)
  x <- gen(N_LARGE)
  expect_true(check_variance(x, expected_var = 1))
})

test_that("make.gen.laplace.tse has correct variance for custom scale", {
  scale <- 2
  expected_var <- 2 * scale^2
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
  expect_true(abs(kurt - 3) < 0.5)
})

test_that("make.gen.laplace.tse validates scale", {
  expect_error(make.gen.laplace.tse(scale = 0), "scale must be positive")
  expect_error(make.gen.laplace.tse(scale = -1), "scale must be positive")
})

test_that("make.gen.laplace.tse handles edge cases without Inf", {
  gen <- make.gen.laplace.tse()
  set.seed(TEST_SEED)
  x <- gen(100000)
  expect_true(all(is.finite(x)))
})


# ==============================================================================
# make.gen.unif.tse
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
  gen <- make.gen.unif.tse()
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
  expect_lt(kurt, 0)
})

test_that("make.gen.unif.tse validates half_width", {
  expect_error(make.gen.unif.tse(half_width = 0), "half_width must be positive")
})


# ==============================================================================
# make.gen.mixnorm.tse
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
  expected_var <- prob1 * sd1^2 + (1 - prob1) * sd2^2
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
  expect_gt(kurt, 0)
})

test_that("make.gen.mixnorm.tse validates parameters", {
  expect_error(make.gen.mixnorm.tse(sd1 = 0), "sd1 and sd2 must be positive")
  expect_error(make.gen.mixnorm.tse(sd2 = -1), "sd1 and sd2 must be positive")
  expect_error(make.gen.mixnorm.tse(prob1 = 0), "prob1 must be between 0 and 1")
  expect_error(make.gen.mixnorm.tse(prob1 = 1), "prob1 must be between 0 and 1")
})


# ==============================================================================
# Edge Cases for All Base Generators
# ==============================================================================

test_that("base generators handle n = 1", {
  expect_length(make.gen.norm.tse()(1), 1)
  expect_length(make.gen.t.tse(df = 5)(1), 1)
  expect_length(make.gen.laplace.tse()(1), 1)
  expect_length(make.gen.unif.tse()(1), 1)
  expect_length(make.gen.mixnorm.tse()(1), 1)
})

test_that("base generators handle large n without error", {
  gen <- make.gen.norm.tse()
  x <- gen(100000)
  expect_length(x, 100000)
  expect_true(all(is.finite(x)))
})

test_that("generator closures capture parameters correctly", {
  gen1 <- make.gen.norm.tse(sd = 1)
  gen2 <- make.gen.norm.tse(sd = 10)
  set.seed(TEST_SEED)
  x1 <- gen1(N_LARGE)
  set.seed(TEST_SEED)
  x2 <- gen2(N_LARGE)
  expect_equal(x2, x1 * 10)
})