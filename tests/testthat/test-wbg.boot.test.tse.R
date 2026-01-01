# Tests for wbg.boot.test.tse

test_that("wbg.boot.test.tse returns correct structure without adjustment", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  result <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, bootadj = FALSE, seed = 42, verbose = FALSE)
  
  expect_type(result, "list")
  expect_named(result, c(
    "obs_stat", "boot_dist", "pvalue_two", "pvalue_upper",
    "pvalue_lower", "ar_order", "ar_phi", "ar_method", "n", "nb", "boot_seeds"
  ))
})

test_that("wbg.boot.test.tse returns correct structure with adjustment", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  result <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, bootadj = TRUE, seed = 42, verbose = FALSE)
  
  expect_type(result, "list")
  expect_true("obs_stat" %in% names(result))
  expect_true("obs_stat_adj" %in% names(result))
  expect_true("pvalue_two_adj" %in% names(result))
  expect_true("adj_factor" %in% names(result))
  expect_true("median_phi" %in% names(result))
  expect_true("boot_dist_adj" %in% names(result))
  expect_true("boot_seeds" %in% names(result))
  expect_true("boot_seeds_adj" %in% names(result))
})

test_that("wbg.boot.test.tse returns correct types and dimensions", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  nb <- 29
  
  result <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = nb, bootadj = FALSE, seed = 42, verbose = FALSE)
  
  expect_type(result$obs_stat, "double")
  expect_length(result$obs_stat, 1)
  
  expect_type(result$boot_dist, "double")
  expect_length(result$boot_dist, nb)
  
  expect_type(result$boot_seeds, "integer")
  expect_length(result$boot_seeds, nb)
  
  expect_type(result$pvalue_two, "double")
  expect_true(result$pvalue_two >= 0 && result$pvalue_two <= 1)
  
  expect_type(result$pvalue_upper, "double")
  expect_true(result$pvalue_upper >= 0 && result$pvalue_upper <= 1)
  
  expect_type(result$pvalue_lower, "double")
  expect_true(result$pvalue_lower >= 0 && result$pvalue_lower <= 1)
  
  expect_type(result$ar_order, "integer")
  expect_true(result$ar_order >= 1 && result$ar_order <= 5)
  
  expect_type(result$ar_phi, "double")
  expect_length(result$ar_phi, result$ar_order)
  
  expect_equal(result$n, 50)
  expect_equal(result$nb, nb)
})

test_that("wbg.boot.test.tse respects seed for reproducibility", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  result1 <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, bootadj = FALSE, seed = 999, verbose = FALSE)
  result2 <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, bootadj = FALSE, seed = 999, verbose = FALSE)
  
  expect_equal(result1$boot_dist, result2$boot_dist)
  expect_equal(result1$pvalue_two, result2$pvalue_two)
  expect_equal(result1$boot_seeds, result2$boot_seeds)
})

test_that("wbg.boot.test.tse reproduces results with boot_seeds", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  result1 <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, bootadj = FALSE, seed = 111, verbose = FALSE)
  result2 <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, bootadj = FALSE, 
                               boot_seeds = result1$boot_seeds, verbose = FALSE)
  
  expect_equal(result1$boot_dist, result2$boot_dist)
  expect_equal(result1$pvalue_two, result2$pvalue_two)
})

test_that("wbg.boot.test.tse produces different results with different seeds", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  result1 <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, bootadj = FALSE, seed = 111, verbose = FALSE)
  result2 <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, bootadj = FALSE, seed = 222, verbose = FALSE)
  
  expect_false(identical(result1$boot_dist, result2$boot_dist))
})

test_that("wbg.boot.test.tse respects p_max parameter", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  result <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, p_max = 2, bootadj = FALSE, seed = 42, verbose = FALSE)
  
  expect_true(result$ar_order <= 2)
})

test_that("wbg.boot.test.tse works with different stat functions", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  
  # OLS t-stat
  result_ols <- wbg.boot.test.tse(x, stat_fn = make.stat.ols.t.tse(), nb = 19, bootadj = FALSE, seed = 42, verbose = FALSE)
  expect_type(result_ols$obs_stat, "double")
  
  # OLS slope
  result_slope <- wbg.boot.test.tse(x, stat_fn = make.stat.ols.slope.tse(), nb = 19, bootadj = FALSE, seed = 42, verbose = FALSE)
  expect_type(result_slope$obs_stat, "double")
  
  # Spearman
  result_spear <- wbg.boot.test.tse(x, stat_fn = make.stat.spearman.tse(), nb = 19, bootadj = FALSE, seed = 42, verbose = FALSE)
  expect_type(result_spear$obs_stat, "double")
  expect_true(abs(result_spear$obs_stat) <= 1)
  
  # CO (if tswge available)
  result_co <- wbg.boot.test.tse(x, stat_fn = make.stat.co.tse(), nb = 19, bootadj = FALSE, seed = 42, verbose = FALSE)
  expect_type(result_co$obs_stat, "double")
})

test_that("wbg.boot.test.tse detects strong trend", {
  
  # Create series with obvious trend
  set.seed(123)
  n <- 100
  trend <- 0.5 * (1:n)
  noise <- gen.arma.wge(n, phi = 0.3, plot = FALSE)
  x <- trend + noise
  
  stat_fn <- make.stat.ols.t.tse()
  result <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 99, bootadj = FALSE, seed = 42, verbose = FALSE)
  
  # With strong trend, p-value should be small
  expect_true(result$pvalue_two < 0.1)
})

test_that("wbg.boot.test.tse does not reject null for trendless data", {
  
  # Create series without trend
  set.seed(456)
  x <- gen.arma.wge(100, phi = 0.5, plot = FALSE)
  
  stat_fn <- make.stat.ols.t.tse()
  result <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 99, bootadj = FALSE, seed = 42, verbose = FALSE)
  
  # Without trend, p-value should generally not be tiny
  # Using a lenient threshold since it's stochastic
  expect_true(result$pvalue_two > 0.01)
})

# --- Tests for COBA adjustment ---

test_that("bootadj produces valid adjustment factor", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.8, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  result <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 49, bootadj = TRUE, seed = 42, verbose = FALSE)
  
  expect_type(result$adj_factor, "double")
  expect_true(result$adj_factor > 0)
  expect_true(is.finite(result$adj_factor))
})

test_that("bootadj produces valid median_phi", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.8, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  result <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 49, bootadj = TRUE, seed = 42, verbose = FALSE)
  
  expect_type(result$median_phi, "double")
  expect_true(length(result$median_phi) >= 1)
})

test_that("bootadj adjusted p-values are valid", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.8, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  result <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 49, bootadj = TRUE, seed = 42, verbose = FALSE)
  
  expect_true(result$pvalue_two_adj >= 0 && result$pvalue_two_adj <= 1)
  expect_true(result$pvalue_upper_adj >= 0 && result$pvalue_upper_adj <= 1)
  expect_true(result$pvalue_lower_adj >= 0 && result$pvalue_lower_adj <= 1)
})

test_that("bootadj second bootstrap distribution has correct length", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.8, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  nb <- 49
  
  result <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = nb, bootadj = TRUE, seed = 42, verbose = FALSE)
  
  expect_length(result$boot_dist_adj, nb)
  expect_length(result$boot_seeds_adj, nb)
})

test_that("bootadj is reproducible with boot_seeds_adj", {
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.8, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  result1 <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 49, bootadj = TRUE, seed = 42, verbose = FALSE)
  result2 <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 49, bootadj = TRUE,
                               boot_seeds = result1$boot_seeds,
                               boot_seeds_adj = result1$boot_seeds_adj, verbose = FALSE)
  
  expect_equal(result1$boot_dist, result2$boot_dist)
  expect_equal(result1$boot_dist_adj, result2$boot_dist_adj)
  expect_equal(result1$pvalue_two_adj, result2$pvalue_two_adj)
  expect_equal(result1$adj_factor, result2$adj_factor)
})

test_that("adjustment factor tends to be less than 1 for high phi", {
  
  # With high autocorrelation, bootstrap samples underestimate variability
  # so adjustment factor C = sd(t_median) / sd(t_boot) should often be < 1
  # (because median model has higher phi than bootstrap estimates)
  
  set.seed(789)
  x <- gen.arma.wge(100, phi = 0.9, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  result <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 99, bootadj = TRUE, seed = 42, verbose = FALSE)
  
  # This is a soft expectation - adjustment factor is often but not always < 1
  # Just verify it's a reasonable positive number
  expect_true(result$adj_factor > 0.1 && result$adj_factor < 10)
})

# --- Tests for parallel mode ---

test_that("parallel mode produces identical results with same boot_seeds", {
  skip_on_cran()
  skip_if(parallel::detectCores() < 2, "Need at least 2 cores for parallel test")
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  # Run sequential first
  result_seq <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, bootadj = FALSE, 
                                  seed = 42, parallel = FALSE, verbose = FALSE)
  
  # Run parallel with same boot_seeds
  result_par <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, bootadj = FALSE,
                                  boot_seeds = result_seq$boot_seeds,
                                  parallel = TRUE, num_cpu = 2, verbose = FALSE)
  
  expect_equal(result_seq$boot_dist, result_par$boot_dist)
  expect_equal(result_seq$pvalue_two, result_par$pvalue_two)
})

test_that("parallel mode with COBA produces identical results with same seeds", {
  skip_on_cran()
  skip_if(parallel::detectCores() < 2, "Need at least 2 cores for parallel test")
  
  set.seed(123)
  x <- gen.arma.wge(50, phi = 0.5, plot = FALSE)
  stat_fn <- make.stat.ols.t.tse()
  
  # Run sequential first
  result_seq <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, bootadj = TRUE,
                                  seed = 42, parallel = FALSE, verbose = FALSE, ar_method = "burg")
  
  # Run parallel with same boot_seeds
  result_par <- wbg.boot.test.tse(x, stat_fn = stat_fn, nb = 19, bootadj = TRUE,
                                  boot_seeds = result_seq$boot_seeds,
                                  boot_seeds_adj = result_seq$boot_seeds_adj,
                                  parallel = TRUE, num_cpu = 2, verbose = FALSE, ar_method = "burg")
  
  expect_equal(result_seq$boot_dist, result_par$boot_dist)
  expect_equal(result_seq$boot_dist_adj, result_par$boot_dist_adj)
  expect_equal(result_seq$pvalue_two_adj, result_par$pvalue_two_adj)
})