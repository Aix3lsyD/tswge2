#' Cochrane-Orcutt Estimation for Trend with Correlated Errors
#'
#' Estimates linear trend in time series data with autocorrelated residuals
#' using the Cochrane-Orcutt procedure. The residuals are modeled as an AR(p)
#' process with order selected by information criterion.
#'
#' @param x A numeric vector containing the time series data.
#' @param maxp Maximum AR order for model selection. Default is 5.
#' @param ar_method Character. Method for AR estimation: \code{"mle"} or
#'   \code{"burg"}. Default is \code{"mle"}.
#' @param criterion Character. Information criterion for order selection:
#'   \code{"aic"}, \code{"aicc"}, or \code{"bic"}. Default is \code{"aic"}.
#' @param n_best Integer. For MLE, number of top models to check for
#'   stationarity. Default is 3.
#' @param tol Numeric. Tolerance for stationarity check. Default is 1.001.
#' 
#' @importFrom stats resid lm
#'
#' @return A list containing:
#'   \item{x}{The original time series.}
#'   \item{z.x}{Residuals from initial OLS trend fit.}
#'   \item{b0hat}{Cochrane-Orcutt estimate of intercept.}
#'   \item{b1hat}{Cochrane-Orcutt estimate of slope.}
#'   \item{z.order}{AR order selected for the residuals.}
#'   \item{z.phi}{AR coefficients for the residual model.}
#'   \item{pvalue}{P-value for test of H0: slope = 0 (assumes uncorrelated
#'     transformed residuals; use \code{\link{wbg.boot.tse}} for more accurate
#'     p-values with correlated errors).}
#'   \item{tco}{t-statistic for the slope coefficient.}
#'   \item{ar_method}{The AR estimation method used.}
#'
#' @details
#' Fits the model
#' \deqn{Y_t = a + bt + Z_t}
#' where \eqn{Z_t} is a stationary AR(p) process satisfying
#' \eqn{\phi(B)Z_t = a_t}.
#'
#' The procedure:
#' \enumerate{
#'   \item Obtains OLS estimates of \eqn{a} and \eqn{b}
#'   \item Computes residuals and fits AR(p) via the specified method with
#'     information criterion order selection
#'   \item Transforms the data using \eqn{\hat{\phi}(B)} to obtain
#'     \eqn{W_t = \hat{\phi}(B)Y_t}
#'   \item Regresses \eqn{W_t} on the transformed time index
#' }
#'
#' Note: The p-value returned assumes the transformed residuals are
#' uncorrelated, which is only approximate. For small to moderate sample sizes
#' with highly correlated errors, use \code{\link{wbg.boot.tse}} for bootstrap-
#' based inference with better significance level control.
#'
#' @references
#' Cochrane, D. and Orcutt, G. H. (1949). "Application of Least Squares
#' Regression to Relationships Containing Auto-Correlated Error Terms."
#' \emph{Journal of the American Statistical Association}, 44(245), 32-61.
#'
#' Woodward, W. A., Bottone, S., and Gray, H. L. (1997). "Improved Tests for
#' Trend in Time Series Data." \emph{Journal of Agricultural, Biological, and
#' Environmental Statistics}, 2(4), 403-416.
#'
#' @seealso \code{\link{wbg.boot.tse}} for bootstrap-based trend test,
#'   \code{\link{aic.ar.tse}}, \code{\link[tswge]{artrans.wge}}
#'
#' @examples
#' \dontrun{
#' # Estimate trend with Cochrane-Orcutt using MLE
#' data(hadley)
#' result <- co.tse(hadley, maxp = 5, ar_method = "mle")
#' cat("Slope estimate:", result$b1hat, "\n")
#' cat("t-statistic:", result$tco, "\n")
#' 
#' # Compare with Burg estimation
#' result_burg <- co.tse(hadley, maxp = 5, ar_method = "burg")
#' cat("Burg t-statistic:", result_burg$tco, "\n")
#' }
#'
#' @export
co.tse <- function(x, maxp = 5, ar_method = c("mle", "burg"),
                   criterion = c("aic", "aicc", "bic"),
                   n_best = 3, tol = 1.001) {
  
  ar_method <- match.arg(ar_method)
  criterion <- match.arg(criterion)
  
  n <- length(x)
  
  # Step 1: OLS fit and residuals
  t1 <- 1:n
  d <- lm(x ~ t1)
  z.x <- resid(d)
  
  # Step 2: Fit AR model to residuals
  ar_fit <- aic.ar.tse(z.x, p_max = maxp, method = ar_method,
                       criterion = criterion, n_best = n_best, tol = tol)
  pp <- ar_fit$p
  phi <- ar_fit$phi
  
  # Step 3: Transform the data
  x.trans <- artrans.wge(x, phi.tr = phi, plottr = FALSE)
  
  # Step 4: Compute transformed time index
  p1 <- pp + 1
  t.co <- rep(0, n)
  
  for (tt in p1:n) {
    t.co[tt] <- tt - sum(phi * (tt - (1:pp)))
  }
  
  # Step 5: Regress transformed data on transformed time
  d.co <- lm(x.trans ~ t.co[p1:n])
  d.co.sum <- summary(d.co)
  
  list(
    x         = x,
    z.x       = z.x,
    b0hat     = d.co$coefficients[1],
    b1hat     = d.co$coefficients[2],
    z.order   = pp,
    z.phi     = phi,
    pvalue    = d.co.sum$coefficients[2, 4],
    tco       = d.co.sum$coefficients[2, 3],
    ar_method = ar_fit$method_used
  )
}