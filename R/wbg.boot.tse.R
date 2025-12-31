#' Bootstrap Test for Trend in Time Series with Correlated Errors
#'
#' Performs a bootstrap-based hypothesis test for linear trend in time series
#' data with autocorrelated residuals. This implements the COB (Cochrane-Orcutt
#' Bootstrap) procedure which provides better control of significance levels
#' than classical methods when residuals are moderately to highly correlated.
#'
#' @param x A numeric vector containing the time series data.
#' @param nb Number of bootstrap replications. Default is 399.
#' @param maxp Maximum AR order for AIC selection in the Cochrane-Orcutt
#'   procedure. Default is 5.
#' @param sn Random seed for reproducibility. If 0 (default), no seed is set.
#' @param maxnull_p Maximum AR order for the null model used to generate
#'   bootstrap samples. If \code{NA} (default), uses \code{maxp}.
#' @param boot_seeds Optional numeric vector of seeds for each bootstrap
#'   replication. If provided, must be of length \code{nb}. If \code{NULL}
#'   (default), seeds are auto-generated (after \code{sn} is set, if provided).
#'   Seeds are always returned for reproducibility. Ignored if
#'   \code{legacy_seed = TRUE}.
#' @param legacy_seed Logical. If \code{TRUE}, uses the legacy seeding approach
#'   (single seed set at start, global RNG state used throughout). This matches
#'   the behavior of \code{wbg.boot.wge()} for verification purposes. Default
#'   is \code{FALSE}.
#'
#' @return A list containing:
#'   \item{p}{AR order selected for the null model.}
#'   \item{phi}{AR coefficients for the null model.}
#'   \item{tco}{Observed t-statistic from the Cochrane-Orcutt procedure.}
#'   \item{pv}{Bootstrap p-value (two-sided).}
#'   \item{boot_dist}{Vector of bootstrap t-statistics.}
#'   \item{boot_seeds}{Vector of seeds used for each bootstrap replication
#'     (returned for reproducibility). \code{NULL} if \code{legacy_seed = TRUE}.}
#'
#' @details
#' The null hypothesis is H0: b = 0 (no linear trend) in the model
#' \deqn{Y_t = a + bt + Z_t}
#' where \eqn{Z_t} is a stationary AR(p) process.
#'
#' Under H0, any trending behavior is attributed to the correlation structure
#' alone. The AR model is fit to the raw series (not detrended) using the Burg
#' (maximum entropy) estimator, which guarantees stationarity. Bootstrap
#' samples are generated from this fitted AR model, and the Cochrane-Orcutt
#' t-statistic is computed for each. The p-value is the proportion of
#' bootstrap |t| values that exceed the observed |t|.
#'
#' @references
#' Woodward, W. A., Bottone, S., and Gray, H. L. (1997). "Improved Tests for
#' Trend in Time Series Data." \emph{Journal of Agricultural, Biological, and
#' Environmental Statistics}, 2(4), 403-416.
#'
#' @seealso \code{\link{co.wge}}, \code{\link{aic.burg.wge}},
#'   \code{\link{gen.arma.wge}}
#'
#' @examples
#' \dontrun{
#' # Test for trend in a time series
#' data(hadley)
#' result = wbg.boot.tse(hadley, nb = 399, maxp = 5, sn = 123)
#' print(result$pv)
#' 
#' # Reproduce exact same results using returned seeds
#' result2 = wbg.boot.tse(hadley, nb = 399, maxp = 5, boot_seeds = result$boot_seeds)
#' identical(result$pv, result2$pv)  # TRUE
#' 
#' # Legacy mode for comparison with wbg.boot.wge()
#' result_legacy = wbg.boot.tse(hadley, nb = 399, maxp = 5, sn = 123, legacy_seed = TRUE)
#' result_old = wbg.boot.wge(hadley, nb = 399, maxp = 5, sn = 123)
#' identical(result_legacy$pv, result_old$pv)  # TRUE
#' }
#'
#' @export
wbg.boot.tse = function(x, nb = 399, maxp = 5, sn = 0, maxnull_p = NA, 
                        boot_seeds = NULL, legacy_seed = FALSE) {
  
  if (sn > 0) set.seed(sn)
  
  if (is.na(maxnull_p)) {
    maxnull_p = maxp
  }
  
  # Handle seeding approach
  if (legacy_seed) {
    # Legacy mode: no per-iteration seeds (matches wbg.boot.wge behavior)
    boot_seeds = NULL
  } else {
    # Modern mode: generate per-iteration seeds if not provided
    if (is.null(boot_seeds)) {
      boot_seeds = sample.int(.Machine$integer.max, nb)
    } else if (length(boot_seeds) != nb) {
      stop("boot_seeds must be of length nb (", nb, "), but got length ", 
           length(boot_seeds))
    }
  }
  
  n = length(x)
  
  
  w = co.wge(x, maxp = maxp)
  obs_t = w$tco
  
  x.aic = aic.burg.wge(x, p = 1:maxnull_p)
  
  boot_t = numeric(nb)
  
  for (ii in 1:nb) {
    iter_seed = if (is.null(boot_seeds)) 0 else boot_seeds[ii]
    xb = gen.arma.wge(n, phi = x.aic$phi, plot = FALSE, sn = iter_seed)
    wb = co.tse(xb, maxp = maxp)
    boot_t[ii] = wb$tco
  }
  
  pv = mean(abs(boot_t) >= abs(obs_t))
  
  list(
    p          = x.aic$p,
    phi        = x.aic$phi,
    tco        = obs_t,
    pv         = pv,
    boot_dist  = boot_t,
    boot_seeds = boot_seeds
  )
}