#' Create Cochrane-Orcutt Statistic Function
#'
#' Creates a function that computes the Cochrane-Orcutt t-statistic for
#' testing trend in a time series with autocorrelated errors.
#'
#' @param maxp Integer. Maximum AR order to consider for the error model.
#'   Default is 5.
#'
#' @return A function that takes a numeric vector \code{x} and returns
#'   the Cochrane-Orcutt t-statistic for the trend coefficient.
#'
#' @seealso \code{\link{wbg.boot.test.tse}}, \code{\link{co.wge}}
#'
#' @examples
#' \dontrun{
#' stat_fn <- make.stat.co.tse(maxp = 3)
#' result <- wbg.boot.test.tse(x, stat_fn = stat_fn)
#' }
#'
#' @export
make.stat.co.tse <- function(maxp = 5) {
  function(x) {
    co.tse(x, maxp = maxp)$tco
  }
}


#' Create OLS t-Statistic Function
#'
#' Creates a function that computes the ordinary least squares t-statistic
#' for the slope in a linear trend regression. This statistic does not
#' account for autocorrelation in the standard error, but can be used
#' within the bootstrap framework which handles the dependence structure.
#'
#' @return A function that takes a numeric vector \code{x} and returns
#'   the OLS t-statistic for the trend coefficient.
#'
#' @seealso \code{\link{wbg.boot.test.tse}}
#'
#' @importFrom stats lm
#'
#' @examples
#' \dontrun{
#' stat_fn <- make.stat.ols.t.tse()
#' result <- wbg.boot.test.tse(x, stat_fn = stat_fn)
#' }
#'
#' @export
make.stat.ols.t.tse <- function() {
  function(x) {
    fit <- lm(x ~ seq_along(x))
    summary(fit)$coefficients[2, 3]
  }
}


#' Create OLS Slope Statistic Function
#'
#' Creates a function that computes the ordinary least squares slope
#' estimate (not t-statistic) for a linear trend.
#'
#' @return A function that takes a numeric vector \code{x} and returns
#'   the OLS slope estimate.
#'
#' @seealso \code{\link{wbg.boot.test.tse}}, \code{\link{make.stat.ols.t.tse}}
#'
#' @examples
#' \dontrun{
#' stat_fn <- make.stat.ols.slope.tse()
#' result <- wbg.boot.test.tse(x, stat_fn = stat_fn)
#' }
#'
#' @export
make.stat.ols.slope.tse <- function() {
  function(x) {
    lm(x ~ seq_along(x))$coefficients[2]
  }
}


#' Create Mann-Kendall Statistic Function
#'
#' Creates a function that computes the Mann-Kendall S statistic for
#' testing monotonic trend. This is a non-parametric, rank-based test
#' that is robust to outliers and non-normality.
#'
#' @return A function that takes a numeric vector \code{x} and returns
#'   the Mann-Kendall S statistic.
#'
#' @note Requires the \pkg{Kendall} package.
#'
#' @seealso \code{\link{wbg.boot.test.tse}}
#'
#' @examples
#' \dontrun{
#' stat_fn <- make.stat.mk.tse()
#' result <- wbg.boot.test.tse(x, stat_fn = stat_fn)
#' }
#'
#' @export
make.stat.mk.tse <- function() {
  if (!requireNamespace("Kendall", quietly = TRUE)) {
    stop("Package 'Kendall' required. Install with install.packages('Kendall')")
  }
  function(x) {
    Kendall::MannKendall(x)$S
  }
}


#' Create Spearman Correlation Statistic Function
#'
#' Creates a function that computes the Spearman rank correlation between
#' the series values and time. This is a non-parametric measure of
#' monotonic trend.
#'
#' @return A function that takes a numeric vector \code{x} and returns
#'   the Spearman correlation with time.
#'
#' @seealso \code{\link{wbg.boot.test.tse}}, \code{\link{make.stat.mk.tse}}
#'
#' @importFrom stats cor
#'
#' @examples
#' \dontrun{
#' stat_fn <- make.stat.spearman.tse()
#' result <- wbg.boot.test.tse(x, stat_fn = stat_fn)
#' }
#'
#' @export
make.stat.spearman.tse <- function() {
  function(x) {
    cor(x, seq_along(x), method = "spearman")
  }
}


#' Create Sen's Slope Statistic Function
#'
#' Creates a function that computes Sen's slope estimator, which is the
#' median of all pairwise slopes. This is a robust, non-parametric
#' estimate of trend magnitude.
#'
#' @return A function that takes a numeric vector \code{x} and returns
#'   Sen's slope estimate.
#'
#' @seealso \code{\link{wbg.boot.test.tse}}, \code{\link{make.stat.mk.tse}}
#'
#' @importFrom stats median
#'
#' @examples
#' \dontrun{
#' stat_fn <- make.stat.sen.tse()
#' result <- wbg.boot.test.tse(x, stat_fn = stat_fn)
#' }
#'
#' @export
make.stat.sen.tse <- function() {
  function(x) {
    n <- length(x)
    slopes <- c()
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        slopes <- c(slopes, (x[j] - x[i]) / (j - i))
      }
    }
    median(slopes)
  }
}


#' Create HAC (Newey-West) t-Statistic Function
#'
#' Creates a function that computes the t-statistic for trend using
#' heteroskedasticity and autocorrelation consistent (HAC) standard errors
#' via the Newey-West estimator.
#'
#' @param lag Integer or NULL. The lag truncation parameter for the
#'   Newey-West estimator. If NULL (default), the bandwidth is selected
#'   automatically.
#'
#' @return A function that takes a numeric vector \code{x} and returns
#'   the HAC-corrected t-statistic for the trend coefficient.
#'
#' @note Requires the \pkg{sandwich} and \pkg{lmtest} packages.
#'
#' @seealso \code{\link{wbg.boot.test.tse}}, \code{\link{make.stat.bn.tse}}
#'
#' @examples
#' \dontrun{
#' stat_fn <- make.stat.hac.tse(lag = 4)
#' result <- wbg.boot.test.tse(x, stat_fn = stat_fn)
#' }
#'
#' @export
make.stat.hac.tse <- function(lag = NULL) {
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("Package 'sandwich' required. Install with install.packages('sandwich')")
  }
  if (!requireNamespace("lmtest", quietly = TRUE)) {
    stop("Package 'lmtest' required. Install with install.packages('lmtest')")
  }
  function(x) {
    fit <- lm(x ~ seq_along(x))
    ct <- lmtest::coeftest(fit, vcov = sandwich::NeweyWest(fit, lag = lag))
    ct[2, 3]
  }
}


#' Create Bloomfield-Nychka t-Statistic Function
#'
#' Creates a function that computes the Bloomfield-Nychka t-statistic,
#' which corrects the standard error for trend using a spectral estimate
#' of the long-run variance (spectral density at frequency zero).
#'
#' @param order.max Integer or NULL. Maximum AR order for spectral
#'   estimation via \code{\link[stats]{spec.ar}}. If NULL (default),
#'   the order is selected automatically.
#'
#' @return A function that takes a numeric vector \code{x} and returns
#'   the Bloomfield-Nychka t-statistic for the trend coefficient.
#'
#' @details
#' The procedure:
#' \enumerate{
#'   \item Fit OLS trend and extract residuals.
#'   \item Estimate spectral density of residuals at frequency zero using
#'     an AR model.
#'   \item Compute corrected variance: Var(slope) = 12 * f(0) / n^3.
#'   \item Return t = slope / sqrt(corrected variance).
#' }
#'
#' @seealso \code{\link{wbg.boot.test.tse}}, \code{\link{make.stat.hac.tse}}
#'
#' @importFrom stats spec.ar coef residuals lm
#'
#' @examples
#' \dontrun{
#' stat_fn <- make.stat.bn.tse(order.max = 10)
#' result <- wbg.boot.test.tse(x, stat_fn = stat_fn)
#' }
#'
#' @export
make.stat.bn.tse <- function(order.max = NULL) {
  function(x) {
    n <- length(x)
    t <- seq_along(x)
    fit <- lm(x ~ t)
    b_hat <- coef(fit)[2]
    resid <- residuals(fit)
    
    sp <- spec.ar(resid, n.freq = 1000, plot = FALSE, order.max = order.max)
    f0 <- sp$spec[1]
    
    var_b <- 12 * f0 / n^3
    b_hat / sqrt(var_b)
  }
}


#' Create Likelihood Ratio Statistic Function
#'
#' Creates a function that computes the likelihood ratio test statistic
#' for trend by comparing AR models with and without a linear trend term.
#'
#' @param order Integer. The AR order for the models. Default is 1.
#'
#' @return A function that takes a numeric vector \code{x} and returns
#'   the likelihood ratio statistic: 2 * (loglik_with_trend - loglik_without).
#'
#' @details
#' Under the null hypothesis of no trend, the likelihood ratio statistic
#' is asymptotically chi-squared with 1 degree of freedom. However, when
#' used within the bootstrap framework, the empirical null distribution
#' is used instead.
#'
#' @seealso \code{\link{wbg.boot.test.tse}}
#'
#' @importFrom stats arima
#'
#' @examples
#' \dontrun{
#' stat_fn <- make.stat.lr.tse(order = 2)
#' result <- wbg.boot.test.tse(x, stat_fn = stat_fn)
#' }
#'
#' @export
make.stat.lr.tse <- function(order = 1) {
  function(x) {
    t <- seq_along(x)
    fit0 <- arima(x, order = c(order, 0, 0))
    fit1 <- arima(x, order = c(order, 0, 0), xreg = t)
    2 * (fit1$loglik - fit0$loglik)
  }
}


#' Create GLS t-Statistic Function
#'
#' Creates a function that computes the generalized least squares t-statistic
#' for trend with AR(p) correlated errors.
#'
#' @param p Integer. The AR order for the error correlation structure.
#'   Default is 1.
#'
#' @return A function that takes a numeric vector \code{x} and returns
#'   the GLS t-statistic for the trend coefficient.
#'
#' @note Requires the \pkg{nlme} package.
#'
#' @seealso \code{\link{wbg.boot.test.tse}}, \code{\link{make.stat.co.tse}}
#'
#' @examples
#' \dontrun{
#' stat_fn <- make.stat.gls.tse(p = 2)
#' result <- wbg.boot.test.tse(x, stat_fn = stat_fn)
#' }
#'
#' @export
make.stat.gls.tse <- function(p = 1) {
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("Package 'nlme' required. Install with install.packages('nlme')")
  }
  function(x) {
    df <- data.frame(x = x, t = seq_along(x))
    fit <- nlme::gls(x ~ t, data = df, correlation = nlme::corARMA(p = p, q = 0))
    summary(fit)$tTable["t", "t-value"]
  }
}