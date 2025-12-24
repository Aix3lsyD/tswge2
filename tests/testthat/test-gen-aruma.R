# test-gen-aruma.R
# Tests for gen.aruma.tse function

library(testthat)

# ==============================================================================
# Basic Functionality
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
  set.seed(111)
  result1 <- gen.aruma.tse(n = 50, phi = 0.5, sn = NULL, plot = FALSE)
  set.seed(222)
  result2 <- gen.aruma.tse(n = 50, phi = 0.5, sn = NULL, plot = FALSE)
  expect_false(all(result1$y == result2$y))
})

test_that("gen.aruma.tse AR(1) has correct ACF structure", {
  phi <- 0.8
  result <- gen.aruma.tse(n = 2000, phi = phi, sn = 42, plot = FALSE)
  acf_vals <- acf(result$y, lag.max = 3, plot = FALSE)$acf
  acf_lag1 <- acf_vals[2]
  expect_true(abs(acf_lag1 - phi) < 0.1)
})


# ==============================================================================
# Custom Innovation Generators
# ==============================================================================

test_that("gen.aruma.tse uses custom innovation generator", {
  unif_gen <- make.gen.unif.tse(half_width = 1)
  result <- gen.aruma.tse(
    n = 500, phi = 0.3, innov_gen = unif_gen, sn = 42, plot = FALSE
  )
  expect_true(all(is.finite(result$y)))
  expect_length(result$y, 500)
})

test_that("gen.aruma.tse with t-distributed innovations", {
  t_gen <- make.gen.t.tse(df = 5)
  result <- gen.aruma.tse(
    n = 500, phi = 0.5, innov_gen = t_gen, sn = 42, plot = FALSE
  )
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse with Laplace innovations", {
  lap_gen <- make.gen.laplace.tse()
  result <- gen.aruma.tse(
    n = 500, phi = 0.5, innov_gen = lap_gen, sn = 42, plot = FALSE
  )
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse with mixture normal innovations", {
  mix_gen <- make.gen.mixnorm.tse(sd1 = 1, sd2 = 5, prob1 = 0.9)
  result <- gen.aruma.tse(
    n = 500, phi = 0.5, innov_gen = mix_gen, sn = 42, plot = FALSE
  )
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse with GARCH innovations", {
  skip_if_not_installed("rugarch")
  garch_gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.1, beta = 0.8)
  result <- gen.aruma.tse(
    n = 500, phi = 0.5, innov_gen = garch_gen, sn = 42, plot = FALSE
  )
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse with skew-t innovations", {
  skip_if_not_installed("sn")
  skt_gen <- make.gen.skt.tse(df = 5, alpha = -2)
  result <- gen.aruma.tse(
    n = 500, phi = 0.5, innov_gen = skt_gen, sn = 42, plot = FALSE
  )
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
})


# ==============================================================================
# MA Components
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
# ARMA Components
# ==============================================================================

test_that("gen.aruma.tse handles ARMA(1,1)", {
  result <- gen.aruma.tse(
    n = 200, phi = 0.5, theta = 0.3, sn = 42, plot = FALSE
  )
  expect_equal(result$p, 1)
  expect_equal(result$q, 1)
  expect_length(result$y, 200)
})

test_that("gen.aruma.tse handles ARMA(2,1)", {
  result <- gen.aruma.tse(
    n = 300, phi = c(0.5, -0.3), theta = 0.2, sn = 42, plot = FALSE
  )
  expect_equal(result$p, 2)
  expect_equal(result$q, 1)
  expect_length(result$y, 300)
})

test_that("gen.aruma.tse handles ARMA(1,2)", {
  result <- gen.aruma.tse(
    n = 300, phi = 0.5, theta = c(0.3, 0.1), sn = 42, plot = FALSE
  )
  expect_equal(result$p, 1)
  expect_equal(result$q, 2)
  expect_length(result$y, 300)
})

test_that("gen.aruma.tse handles ARMA(2,2)", {
  result <- gen.aruma.tse(
    n = 300, phi = c(0.5, -0.2), theta = c(0.3, 0.1), sn = 42, plot = FALSE
  )
  expect_equal(result$p, 2)
  expect_equal(result$q, 2)
  expect_length(result$y, 300)
})


# ==============================================================================
# Higher Order AR
# ==============================================================================

test_that("gen.aruma.tse handles AR(2)", {
  result <- gen.aruma.tse(n = 300, phi = c(0.5, -0.3), sn = 42, plot = FALSE)
  expect_equal(result$p, 2)
  expect_length(result$y, 300)
})

test_that("gen.aruma.tse handles AR(3)", {
  result <- gen.aruma.tse(
    n = 300, phi = c(0.4, -0.2, 0.1), sn = 42, plot = FALSE
  )
  expect_equal(result$p, 3)
  expect_length(result$y, 300)
})


# ==============================================================================
# Differencing
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
# Seasonal and Nonstationary Components
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
  result <- gen.aruma.tse(
    n = 200, phi = 0.5, lambda = c(1.2, -0.5), sn = 42, plot = FALSE
  )
  expect_length(result$y, 200)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse handles lambda + s combined", {
  result <- gen.aruma.tse(
    n = 200, phi = 0.3, lambda = c(0.8), s = 4, sn = 42, plot = FALSE
  )
  expect_equal(result$s, 4)
  expect_length(result$y, 200)
  expect_true(all(is.finite(result$y)))
})

test_that("gen.aruma.tse handles d + s combined", {
  result <- gen.aruma.tse(
    n = 200, phi = 0.3, d = 1, s = 12, sn = 42, plot = FALSE
  )
  expect_equal(result$d, 1)
  expect_equal(result$s, 12)
  expect_length(result$y, 200)
})


# ==============================================================================
# Near Unit Root (Adaptive Burn-in)
# ==============================================================================

test_that("gen.aruma.tse handles near-unit-root AR(1)", {
  result <- gen.aruma.tse(n = 500, phi = 0.99, sn = 42, plot = FALSE)
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
  expect_true(var(result$y) < 1e6)
})

test_that("gen.aruma.tse handles near-unit-root AR(2)", {
  result <- gen.aruma.tse(
    n = 500, phi = c(0.5, 0.49), sn = 42, plot = FALSE
  )
  expect_length(result$y, 500)
  expect_true(all(is.finite(result$y)))
})


# ==============================================================================
# Pure Noise (No AR/MA)
# ==============================================================================

test_that("gen.aruma.tse handles pure white noise", {
  result <- gen.aruma.tse(n = 200, sn = 42, plot = FALSE)
  expect_equal(result$p, 0)
  expect_equal(result$q, 0)
  expect_length(result$y, 200)
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
# Output Methods
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
  built <- ggplot2::ggplot_build(result$plot)
  expect_true(!is.null(built))
})

test_that("gen.aruma.tse plot method works", {
  result <- gen.aruma.tse(n = 100, phi = 0.5, plot = FALSE, sn = 42)
  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


# ==============================================================================
# Works with All Generator Types
# ==============================================================================

test_that("gen.aruma.tse works with all base generators", {
  n <- 200
  phi <- 0.5
  
  result_norm <- gen.aruma.tse(
    n = n, phi = phi, innov_gen = make.gen.norm.tse(), sn = 42, plot = FALSE
  )
  expect_length(result_norm$y, n)
  
  result_t <- gen.aruma.tse(
    n = n, phi = phi, innov_gen = make.gen.t.tse(df = 5), sn = 42, plot = FALSE
  )
  expect_length(result_t$y, n)
  
  result_lap <- gen.aruma.tse(
    n = n, phi = phi, innov_gen = make.gen.laplace.tse(), sn = 42, plot = FALSE
  )
  expect_length(result_lap$y, n)
  
  result_unif <- gen.aruma.tse(
    n = n, phi = phi, innov_gen = make.gen.unif.tse(), sn = 42, plot = FALSE
  )
  expect_length(result_unif$y, n)
  
  result_mix <- gen.aruma.tse(
    n = n, phi = phi, innov_gen = make.gen.mixnorm.tse(), sn = 42, plot = FALSE
  )
  expect_length(result_mix$y, n)
})


# ==============================================================================
# Distributional Properties Propagate
# ==============================================================================

test_that("heavy-tailed innovations produce heavier-tailed ARMA output", {
  n <- 5000
  phi <- 0.3
  
  set.seed(TEST_SEED)
  result_norm <- gen.aruma.tse(
    n = n, phi = phi, innov_gen = make.gen.norm.tse(), plot = FALSE
  )
  
  set.seed(TEST_SEED)
  result_t <- gen.aruma.tse(
    n = n, phi = phi, innov_gen = make.gen.t.tse(df = 3), plot = FALSE
  )
  
  kurt_norm <- excess_kurtosis(result_norm$y)
  kurt_t <- excess_kurtosis(result_t$y)
  
  expect_gt(kurt_t, kurt_norm)
})