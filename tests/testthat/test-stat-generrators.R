# Tests for stat generator functions

# --- Test data setup ---
test_x <- c(1.2, 2.5, 2.8, 4.1, 5.3, 5.9, 7.2, 8.1, 9.0, 10.2)
test_x_notrend <- c(0.5, -0.3, 0.8, -0.1, 0.2, -0.4, 0.6, -0.2, 0.3, 0.1)


# --- make.stat.ols.t.tse ---

test_that("make.stat.ols.t.tse returns a function", {
  stat_fn <- make.stat.ols.t.tse()
  expect_type(stat_fn, "closure")
})

test_that("make.stat.ols.t.tse computes correct t-statistic", {
  stat_fn <- make.stat.ols.t.tse()
  result <- stat_fn(test_x)
  
  # Compare with direct calculation
  fit <- lm(test_x ~ seq_along(test_x))
  expected <- summary(fit)$coefficients[2, 3]
  
  expect_equal(result, expected)
})

test_that("make.stat.ols.t.tse returns single numeric", {
  stat_fn <- make.stat.ols.t.tse()
  result <- stat_fn(test_x)
  
  expect_type(result, "double")
  expect_length(result, 1)
})


# --- make.stat.ols.slope.tse ---

test_that("make.stat.ols.slope.tse returns a function", {
  stat_fn <- make.stat.ols.slope.tse()
  expect_type(stat_fn, "closure")
})

test_that("make.stat.ols.slope.tse computes correct slope", {
  stat_fn <- make.stat.ols.slope.tse()
  result <- stat_fn(test_x)
  
  fit <- lm(test_x ~ seq_along(test_x))
  expected <- unname(coef(fit)[2])
  
  expect_equal(unname(result), expected)
})

test_that("make.stat.ols.slope.tse returns positive slope for increasing data", {
  stat_fn <- make.stat.ols.slope.tse()
  result <- stat_fn(test_x)
  
  expect_true(result > 0)
})


# --- make.stat.spearman.tse ---

test_that("make.stat.spearman.tse returns a function", {
  stat_fn <- make.stat.spearman.tse()
  expect_type(stat_fn, "closure")
})

test_that("make.stat.spearman.tse computes correct correlation", {
  stat_fn <- make.stat.spearman.tse()
  result <- stat_fn(test_x)
  
  expected <- cor(test_x, seq_along(test_x), method = "spearman")
  
  expect_equal(result, expected)
})

test_that("make.stat.spearman.tse returns value in [-1, 1]", {
  stat_fn <- make.stat.spearman.tse()
  
  result1 <- stat_fn(test_x)
  expect_true(result1 >= -1 && result1 <= 1)
  
  result2 <- stat_fn(test_x_notrend)
  expect_true(result2 >= -1 && result2 <= 1)
})

test_that("make.stat.spearman.tse returns 1 for perfectly monotonic data", {
  stat_fn <- make.stat.spearman.tse()
  perfect <- 1:10
  result <- stat_fn(perfect)
  
  expect_equal(result, 1)
})


# --- make.stat.sen.tse ---

test_that("make.stat.sen.tse returns a function", {
  stat_fn <- make.stat.sen.tse()
  expect_type(stat_fn, "closure")
})
test_that("make.stat.sen.tse returns single numeric", {
  stat_fn <- make.stat.sen.tse()
  result <- stat_fn(test_x)
  
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("make.stat.sen.tse returns correct slope for linear data", {
  stat_fn <- make.stat.sen.tse()
  linear <- c(2, 4, 6, 8, 10)  # slope = 2
  result <- stat_fn(linear)
  
  expect_equal(result, 2)
})

test_that("make.stat.sen.tse is robust to outliers", {
  stat_fn <- make.stat.sen.tse()
  
  # Linear data with one outlier
  with_outlier <- c(2, 4, 100, 8, 10)
  linear <- c(2, 4, 6, 8, 10)
  
  result_outlier <- stat_fn(with_outlier)
  result_linear <- stat_fn(linear)
  
  # Sen's slope should be more robust than OLS
  ols_fn <- make.stat.ols.slope.tse()
  ols_outlier <- ols_fn(with_outlier)
  
  # The Sen estimate should be closer to 2 than OLS
  expect_true(abs(result_outlier - 2) < abs(ols_outlier - 2))
})


# --- make.stat.bn.tse ---

test_that("make.stat.bn.tse returns a function", {
  stat_fn <- make.stat.bn.tse()
  expect_type(stat_fn, "closure")
})

test_that("make.stat.bn.tse returns single numeric", {
  stat_fn <- make.stat.bn.tse()
  # Need longer series for spectral estimation
  set.seed(123)
  x <- cumsum(rnorm(50)) + 0.1 * (1:50)
  result <- stat_fn(x)
  
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("make.stat.bn.tse respects order.max parameter", {
  stat_fn_default <- make.stat.bn.tse()
  stat_fn_low <- make.stat.bn.tse(order.max = 2)
  
  set.seed(123)
  x <- rnorm(50) + 0.05 * (1:50)
  
  # Both should return valid results
  result_default <- stat_fn_default(x)
  result_low <- stat_fn_low(x)
  
  expect_type(result_default, "double")
  expect_type(result_low, "double")
})


# --- make.stat.lr.tse ---

test_that("make.stat.lr.tse returns a function", {
  stat_fn <- make.stat.lr.tse()
  expect_type(stat_fn, "closure")
})

test_that("make.stat.lr.tse returns non-negative statistic", {
  stat_fn <- make.stat.lr.tse(order = 1)
  set.seed(123)
  x <- rnorm(50) + 0.1 * (1:50)
  result <- stat_fn(x)
  
  # LR statistic should be non-negative
  expect_true(result >= 0)
})

test_that("make.stat.lr.tse respects order parameter", {
  stat_fn_1 <- make.stat.lr.tse(order = 1)
  stat_fn_2 <- make.stat.lr.tse(order = 2)
  
  set.seed(123)
  x <- rnorm(50) + 0.1 * (1:50)
  
  result_1 <- stat_fn_1(x)
  result_2 <- stat_fn_2(x)
  
  # Both should return valid results (may differ due to different AR orders)
  expect_type(result_1, "double")
  expect_type(result_2, "double")
})


# --- make.stat.co.tse ---

test_that("make.stat.co.tse returns a function", {
 
  
  stat_fn <- make.stat.co.tse()
  expect_type(stat_fn, "closure")
})

test_that("make.stat.co.tse returns single numeric", {
 
  
  stat_fn <- make.stat.co.tse()
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  result <- stat_fn(x)
  
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("make.stat.co.tse respects maxp parameter", {
 
  
  stat_fn_3 <- make.stat.co.tse(maxp = 3)
  stat_fn_5 <- make.stat.co.tse(maxp = 5)
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  
  result_3 <- stat_fn_3(x)
  result_5 <- stat_fn_5(x)
  
  expect_type(result_3, "double")
  expect_type(result_5, "double")
})


# --- make.stat.mk.tse ---

test_that("make.stat.mk.tse requires Kendall package", {
  skip_if_not_installed("Kendall")
  
  stat_fn <- make.stat.mk.tse()
  expect_type(stat_fn, "closure")
})

test_that("make.stat.mk.tse returns correct S statistic", {
  skip_if_not_installed("Kendall")
  
  stat_fn <- make.stat.mk.tse()
  result <- stat_fn(test_x)
  
  expected <- Kendall::MannKendall(test_x)$S
  
  expect_equal(result, expected)
})


# --- make.stat.hac.tse ---

test_that("make.stat.hac.tse requires sandwich and lmtest", {
  skip_if_not_installed("sandwich")
  skip_if_not_installed("lmtest")
  
  stat_fn <- make.stat.hac.tse()
  expect_type(stat_fn, "closure")
})

test_that("make.stat.hac.tse returns single numeric", {
  skip_if_not_installed("sandwich")
  skip_if_not_installed("lmtest")
  
  stat_fn <- make.stat.hac.tse()
  result <- stat_fn(test_x)
  
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("make.stat.hac.tse respects lag parameter", {
  skip_if_not_installed("sandwich")
  skip_if_not_installed("lmtest")
  
  stat_fn_auto <- make.stat.hac.tse(lag = NULL)
  stat_fn_fixed <- make.stat.hac.tse(lag = 2)
  
  result_auto <- stat_fn_auto(test_x)
  result_fixed <- stat_fn_fixed(test_x)
  
  expect_type(result_auto, "double")
  expect_type(result_fixed, "double")
})


# --- make.stat.gls.tse ---

test_that("make.stat.gls.tse requires nlme package", {
  skip_if_not_installed("nlme")
  
  stat_fn <- make.stat.gls.tse()
  expect_type(stat_fn, "closure")
})

test_that("make.stat.gls.tse returns single numeric", {
  skip_if_not_installed("nlme")
  
  stat_fn <- make.stat.gls.tse(p = 1)
  set.seed(123)
  x <- rnorm(30) + 0.1 * (1:30)
  result <- stat_fn(x)
  
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("make.stat.gls.tse respects p parameter", {
  skip_if_not_installed("nlme")
  
  stat_fn_1 <- make.stat.gls.tse(p = 1)
  stat_fn_2 <- make.stat.gls.tse(p = 2)
  
  set.seed(123)
  x <- rnorm(30) + 0.1 * (1:30)
  
  result_1 <- stat_fn_1(x)
  result_2 <- stat_fn_2(x)
  
  expect_type(result_1, "double")
  expect_type(result_2, "double")
})


# --- General generator contract tests ---

test_that("all generators return functions with consistent interface", {
 
  
  generators <- list(
    make.stat.co.tse(),
    make.stat.ols.t.tse(),
    make.stat.ols.slope.tse(),
    make.stat.spearman.tse(),
    make.stat.sen.tse(),
    make.stat.bn.tse(),
    make.stat.lr.tse()
  )
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  
  for (stat_fn in generators) {
    result <- stat_fn(x)
    expect_type(result, "double")
    expect_length(result, 1)
    expect_false(is.na(result))
  }
})

test_that("generators with optional packages check dependencies", {
  # These should error informatively if packages missing
  # We just verify the functions exist and are callable when packages present
  
  if (requireNamespace("Kendall", quietly = TRUE)) {
    expect_type(make.stat.mk.tse(), "closure")
  }
  
  if (requireNamespace("sandwich", quietly = TRUE) &&
      requireNamespace("lmtest", quietly = TRUE)) {
    expect_type(make.stat.hac.tse(), "closure")
  }
  
  if (requireNamespace("nlme", quietly = TRUE)) {
    expect_type(make.stat.gls.tse(), "closure")
  }
})