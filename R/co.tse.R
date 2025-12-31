#' Cochrane-Orcutt Estimation for Trend with Correlated Errors
#'
#' Estimates linear trend in time series data with autocorrelated residuals
#' using the Cochrane-Orcutt procedure. The residuals are modeled as an AR(p)
#' process with order selected by AIC.
#'
#' @param x A numeric vector containing the time series data.
#' @param maxp Maximum AR order for AIC selection. Default is 5.
#'
#' @return A list containing:
#'   \item{x}{The original time series.}
#'   \item{z.x}{Residuals from initial OLS trend fit.}
#'   \item{b0hat}{Cochrane-Orcutt estimate of intercept.}
#'   \item{b1hat}{Cochrane-Orcutt estimate of slope.}
#'   \item{z.order}{AR order selected by AIC for the residuals.}
#'   \item{z.phi}{AR coefficients for the residual model.}
#'   \item{pvalue}{P-value for test of H0: slope = 0 (assumes uncorrelated
#'     transformed residuals; use \code{\link{wbg.boot.tse}} for more accurate
#'     p-values with correlated errors).}
#'   \item{tco}{t-statistic for the slope coefficient.}
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
#'   \item Computes residuals and fits AR(p) via Burg estimation with AIC
#'     order selection
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
#'   \code{\link{aic.burg.wge}}, \code{\link{artrans.wge}}
#'
#' @examples
#' \dontrun{
#' # Estimate trend with Cochrane-Orcutt
#' data(hadley)
#' result = co.tse(hadley, maxp = 5)
#' cat("Slope estimate:", result$b1hat, "\n")
#' cat("t-statistic:", result$tco, "\n")
#' }
#'
#' @export
co.tse = function(x, maxp = 5) {
  
  n = length(x)
  
  
  t1 = 1:n
  d = lm(x ~ t1)
  z.x = resid(d)
  
  aic.z = aic.burg.wge(z.x, p = 1:maxp)
  pp = aic.z$p
  phi = aic.z$phi
  
  x.trans = artrans.wge(x, phi.tr = phi, plottr = FALSE)
  
  p1 = pp + 1
  t.co = rep(0, n)
  
  for (tt in p1:n) {
    t.co[tt] = tt - sum(phi * (tt - (1:pp)))
  }
  
  d.co = lm(x.trans ~ t.co[p1:n])
  d.co.sum = summary(d.co)
  
  list(
    x       = x,
    z.x     = z.x,
    b0hat   = d.co$coefficients[1],
    b1hat   = d.co$coefficients[2],
    z.order = pp,
    z.phi   = phi,
    pvalue  = d.co.sum$coefficients[2, 4],
    tco     = d.co.sum$coefficients[2, 3]
  )
}