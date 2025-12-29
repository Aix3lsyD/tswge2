#' Woodward-Bottone-Gray Bootstrap Test for Trend
#'
#' Performs a bootstrap hypothesis test for trend in a time series while
#' accounting for autocorrelation. Under the null hypothesis, the series
#' is assumed to be a stationary AR process with no trend.
#'
#' @param x Numeric vector. The observed time series.
#' @param stat_fn Function. A function that takes a numeric vector and returns
#'   a single numeric test statistic. Use the \code{make.stat.*.tse} generator
#'   functions to create appropriate statistic functions.
#' @param nb Integer. Number of bootstrap replicates. Default is 399.
#' @param p_max Integer. Maximum AR order to consider when fitting the null
#'   model. Order is selected by AIC. Default is 5.
#' @param bootadj Logical. If TRUE, performs the COBA variance adjustment using
#'   a second bootstrap from the median model. This corrects for bias in AR
#'   parameter estimates that can inflate Type I error when autocorrelation
#'   is strong. Default is TRUE.
#' @param seed Integer. Random seed for reproducibility. If 0 (default), no
#'   seed is set.
#'
#' @return A list containing:
#' \describe{
#'   \item{obs_stat}{The observed test statistic from the original series.}
#'   \item{boot_dist}{Numeric vector of length \code{nb} containing the
#'     bootstrap distribution of the test statistic under the null.}
#'   \item{pvalue_two}{Two-sided p-value: proportion of bootstrap statistics
#'     with absolute value >= |obs_stat|.}
#'   \item{pvalue_upper}{Upper one-sided p-value: proportion of bootstrap
#'     statistics >= obs_stat.}
#'   \item{pvalue_lower}{Lower one-sided p-value: proportion of bootstrap
#'     statistics <= obs_stat.}
#'   \item{ar_order}{Integer. The AR order selected by AIC.}
#'   \item{ar_phi}{Numeric vector. The AR coefficients of the fitted null model.}
#'   \item{n}{Integer. Length of the original series.}
#'   \item{nb}{Integer. Number of bootstrap replicates performed.}
#'   \item{obs_stat_adj}{(If \code{bootadj=TRUE}) The variance-adjusted observed
#'     test statistic.}
#'   \item{boot_dist_adj}{(If \code{bootadj=TRUE}) Bootstrap distribution from
#'     the median model (second bootstrap).}
#'   \item{pvalue_two_adj}{(If \code{bootadj=TRUE}) Adjusted two-sided p-value.}
#'   \item{pvalue_upper_adj}{(If \code{bootadj=TRUE}) Adjusted upper one-sided p-value.}
#'   \item{pvalue_lower_adj}{(If \code{bootadj=TRUE}) Adjusted lower one-sided p-value.}
#'   \item{adj_factor}{(If \code{bootadj=TRUE}) The variance adjustment factor C.}
#'   \item{median_phi}{(If \code{bootadj=TRUE}) AR coefficients of the median model.}
#' }
#'
#' @details
#' The test procedure (COB):
#' \enumerate{
#'   \item Compute the test statistic on the observed series.
#'   \item Fit an AR(p) model to the observed series under H0 (no trend),
#'     with order selected by AIC.
#'   \item Generate \code{nb} bootstrap series from the fitted AR model.
#'   \item Compute the test statistic on each bootstrap series.
#'   \item Calculate p-values by comparing the observed statistic to the
#'     bootstrap distribution.
#' }
#'
#' When \code{bootadj=TRUE}, the COBA adjustment is applied:
#' \enumerate{
#'   \item From the first bootstrap, collect AR coefficient estimates and find
#'     the "median model" (the model whose phi(1) = 1 - sum(phi) is the median).
#'   \item Generate a second set of bootstrap samples from this median model.
#'   \item Compute adjustment factor C = sd(t_median) / sd(t_bootstrap).
#'   \item Scale the observed statistic: t_adj = C * t_obs.
#' }
#'
#' This adjustment corrects for the downward bias in AR parameter estimates
#' that causes bootstrap samples to have less autocorrelation than the true
#' process, which would otherwise inflate Type I error rates.
#'
#' @references
#' Woodward, W. A., Bottone, S., and Gray, H. L. (1997). Improved tests for
#' trend in time series data. \emph{Journal of Agricultural, Biological, and
#' Environmental Statistics}, 2(4), 403-416.
#'
#' @seealso
#' \code{\link{make.stat.co.tse}}, \code{\link{make.stat.ols.t.tse}},
#' \code{\link{make.stat.bn.tse}}, \code{\link{make.stat.mk.tse}}
#'
#' @examples
#' \dontrun{
#' # Generate a series with no trend
#' x <- gen.arma.wge(100, phi = 0.7, plot = FALSE)
#'
#' # Test using Cochrane-Orcutt statistic
#' result <- wbg.boot.test.tse(x, stat_fn = make.stat.co.tse(), nb = 399)
#' result$pvalue_two
#' result$pvalue_two_adj  # variance-adjusted p-value
#'
#' # Visualize the bootstrap distribution
#' hist(result$boot_dist, main = "Bootstrap Distribution")
#' abline(v = result$obs_stat, col = "red", lwd = 2)
#' abline(v = result$obs_stat_adj, col = "blue", lwd = 2, lty = 2)
#'
#' # Compare multiple statistics
#' co_result <- wbg.boot.test.tse(x, stat_fn = make.stat.co.tse())
#' bn_result <- wbg.boot.test.tse(x, stat_fn = make.stat.bn.tse())
#' }
#'
#' @export
wbg.boot.test.tse <- function(x, stat_fn, nb = 399, p_max = 5, bootadj = TRUE, seed = 0) {
  
  if (seed > 0) set.seed(seed)
  
  n <- length(x)
  
  # Observed statistic
  obs_stat <- stat_fn(x)
  
  # Fit AR model under H0: no trend, just stationary AR
  ar_fit <- aic.burg.wge(x, p = 1:p_max)
  
  # First bootstrap: generate distribution and collect AR fits
  boot_stats <- rep(NA, nb)
  boot_phi_list <- vector("list", nb)
  boot_phi1 <- rep(NA, nb)  # phi(1) = 1 - sum(phi) for each bootstrap
  
  for (i in 1:nb) {
    xb <- gen.arma.wge(n, phi = ar_fit$phi, plot = FALSE)
    boot_stats[i] <- stat_fn(xb)
    
    if (bootadj) {
      # Fit AR to bootstrap sample and store coefficients
      ar_boot <- aic.burg.wge(xb, p = 1:p_max)
      boot_phi_list[[i]] <- ar_boot$phi
      boot_phi1[i] <- 1 - sum(ar_boot$phi)
    }
  }
  
  # Compute unadjusted p-values
  pvalue_two <- mean(abs(boot_stats) >= abs(obs_stat))
  pvalue_upper <- mean(boot_stats >= obs_stat)
  pvalue_lower <- mean(boot_stats <= obs_stat)
  
  # Build base result list
  result <- list(
    obs_stat     = obs_stat,
    boot_dist    = boot_stats,
    pvalue_two   = pvalue_two,
    pvalue_upper = pvalue_upper,
    pvalue_lower = pvalue_lower,
    ar_order     = ar_fit$p,
    ar_phi       = ar_fit$phi,
    n            = n,
    nb           = nb
  )
  
  # COBA adjustment: second bootstrap from median model
  if (bootadj) {
    # Find median model: the one whose phi(1) is the median
    median_idx <- which.min(abs(boot_phi1 - median(boot_phi1)))
    median_phi <- boot_phi_list[[median_idx]]
    
    # Second bootstrap from median model
    boot_stats_adj <- rep(NA, nb)
    for (i in 1:nb) {
      xb <- gen.arma.wge(n, phi = median_phi, plot = FALSE)
      boot_stats_adj[i] <- stat_fn(xb)
    }
    
    # Compute adjustment factor: C = sd(t_median) / sd(t_bootstrap)
    adj_factor <- sd(boot_stats_adj) / sd(boot_stats)
    
    # Adjust observed statistic
    obs_stat_adj <- adj_factor * obs_stat
    
    # Compute adjusted p-values (compare adjusted stat to original bootstrap dist)
    pvalue_two_adj <- mean(abs(boot_stats) >= abs(obs_stat_adj))
    pvalue_upper_adj <- mean(boot_stats >= obs_stat_adj)
    pvalue_lower_adj <- mean(boot_stats <= obs_stat_adj)
    
    # Add adjusted results to output
    result$obs_stat_adj     <- obs_stat_adj
    result$boot_dist_adj    <- boot_stats_adj
    result$pvalue_two_adj   <- pvalue_two_adj
    result$pvalue_upper_adj <- pvalue_upper_adj
    result$pvalue_lower_adj <- pvalue_lower_adj
    result$adj_factor       <- adj_factor
    result$median_phi       <- median_phi
  }
  
  result
}