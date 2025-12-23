#' @title Innovation and Process Generators
#' @description Factory functions for creating innovation generators and
#'   stochastic process generators for time series simulation.
#' @name generators
#' @keywords internal
"_PACKAGE"


# ==============================================================================
# Normal Innovation Generator
# ==============================================================================

#' Create a Normal Innovation Generator
#'
#' Factory function that returns a normal (Gaussian) innovation generator.
#'
#' @param mean Mean of the distribution (default 0)
#' @param sd Standard deviation (default 1)
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} normal innovations.
#'
#' @export
#'
#' @examples
#' # Default standard normal
#' norm_gen <- make.gen.norm.tse()
#' innovations <- norm_gen(100)
#'
#' # Higher variance innovations
#' norm_gen_hv <- make.gen.norm.tse(sd = 2)
#' innovations_hv <- norm_gen_hv(100)
#'
#' # Use with ARMA simulation
#' # result <- gen.aruma.wge(n = 500, phi = 0.7, innov_gen = make.gen.norm.tse())
make.gen.norm.tse <- function(mean = 0, sd = 1) {
  if (sd <= 0) {
    stop("sd must be positive")
  }
  
  force(mean)
  force(sd)
  
  function(n) {
    rnorm(n, mean = mean, sd = sd)
  }
}


# ==============================================================================
# Student's t Innovation Generator
# ==============================================================================

#' Create a Student's t Innovation Generator
#'
#' Factory function that returns a t-distributed innovation generator.
#' By default, the output is scaled to have unit variance for comparability
#' with normal innovations.
#'
#' @param df Degrees of freedom. Must be > 2 for finite variance when
#'   \code{scale = TRUE}.
#' @param scale Logical. If \code{TRUE} (default), scale to unit variance.
#'   Requires \code{df > 2}.
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} t-distributed innovations.
#'
#' @details
#' The variance of a t-distribution with \code{df} degrees of freedom is
#' \code{df / (df - 2)} for \code{df > 2}. When \code{scale = TRUE}, the
#' output is divided by \code{sqrt(df / (df - 2))} to achieve unit variance.
#'
#' For \code{df <= 2}, the variance is infinite, so scaling is not possible
#' and will be disabled with a warning.
#'
#' @export
#'
#' @examples
#' # Heavy-tailed innovations (df = 5)
#' t_gen <- make.gen.t.tse(df = 5)
#' innovations <- t_gen(100)
#' var(innovations)  # approximately 1
#'
#' # Very heavy tails (df = 3)
#' t_gen_heavy <- make.gen.t.tse(df = 3)
#'
#' # Unscaled t innovations
#' t_gen_raw <- make.gen.t.tse(df = 5, scale = FALSE)
#' var(t_gen_raw(10000))
make.gen.t.tse <- function(df, scale = TRUE) {
  if (df <= 0) {
    stop("df must be positive")
  }
  
  if (df <= 2 && scale) {
    warning("df <= 2 has infinite variance; scale set to FALSE")
    scale <- FALSE
  }
  
  force(df)
  force(scale)
  
  function(n) {
    x <- rt(n, df = df)
    if (scale && df > 2) {
      # Variance of t is df/(df-2), so scale to unit variance
      x <- x * sqrt((df - 2) / df)
    }
    x
  }
}


# ==============================================================================
# GARCH Process Generator
# ==============================================================================

#' Create a GARCH Process Generator
#'
#' Factory function that returns a GARCH process generator using the
#' \pkg{rugarch} package. Supports various GARCH variants and innovation
#' distributions.
#'
#' @param omega The constant term in the variance equation. Must be positive.
#' @param alpha Numeric vector of ARCH coefficients. Length determines the
#'   ARCH order (q).
#' @param beta Numeric vector of GARCH coefficients (default \code{NULL} for
#'   pure ARCH). Length determines the GARCH order (p).
#' @param model Character string specifying the GARCH variant. One of:
#'   \describe{
#'     \item{\code{"sGARCH"}}{Standard GARCH (default)}
#'     \item{\code{"eGARCH"}}{Exponential GARCH (Nelson, 1991)}
#'     \item{\code{"gjrGARCH"}}{GJR-GARCH with leverage (Glosten et al., 1993)}
#'     \item{\code{"apARCH"}}{Asymmetric Power ARCH (Ding et al., 1993)}
#'     \item{\code{"iGARCH"}}{Integrated GARCH}
#'     \item{\code{"csGARCH"}}{Component sGARCH}
#'   }
#' @param distribution Character string specifying the innovation distribution.
#'   One of:
#'   \describe{
#'     \item{\code{"norm"}}{Normal distribution (default)}
#'     \item{\code{"std"}}{Student's t distribution}
#'     \item{\code{"ged"}}{Generalized Error Distribution}
#'     \item{\code{"snorm"}}{Skew normal}
#'     \item{\code{"sstd"}}{Skew Student's t}
#'     \item{\code{"sged"}}{Skew GED}
#'     \item{\code{"nig"}}{Normal Inverse Gaussian}
#'     \item{\code{"jsu"}}{Johnson's SU}
#'   }
#' @param distribution.params Named list of additional distribution parameters
#'   (e.g., \code{list(shape = 5)} for Student's t degrees of freedom).
#'
#' @return A function with signature \code{function(n)} that generates
#'   \code{n} observations from the specified GARCH process.
#'
#' @details
#' The standard GARCH(q, p) model specifies the conditional variance as:
#' \deqn{\sigma_t^2 = \omega + \sum_{i=1}^{q} \alpha_i \epsilon_{t-i}^2 +
#'   \sum_{j=1}^{p} \beta_j \sigma_{t-j}^2}
#'
#' For stationarity of standard GARCH, we require
#' \code{sum(alpha) + sum(beta) < 1}.
#'
#' @seealso \code{\link[rugarch]{ugarchspec}}, \code{\link[rugarch]{ugarchpath}}
#'
#' @export
#'
#' @examples
#' # Standard GARCH(1,1) with normal innovations
#' garch11_gen <- make.gen.garch.tse(
#'   omega = 0.1,
#'   alpha = 0.15,
#'   beta = 0.8
#' )
#' y <- garch11_gen(500)
#' plot(y, type = "l", main = "GARCH(1,1) Simulation")
#'
#' # ARCH(2) process (no GARCH terms)
#' arch2_gen <- make.gen.garch.tse(
#'   omega = 0.2,
#'   alpha = c(0.3, 0.2)
#' )
#'
#' # GARCH(1,1) with Student's t innovations
#' garch_t_gen <- make.gen.garch.tse(
#'   omega = 0.1,
#'   alpha = 0.1,
#'   beta = 0.85,
#'   distribution = "std",
#'   distribution.params = list(shape = 5)
#' )
#'
#' # GJR-GARCH with leverage effect
#' gjr_gen <- make.gen.garch.tse(
#'   omega = 0.1,
#'   alpha = 0.05,
#'   beta = 0.9,
#'   model = "gjrGARCH"
#' )
make.gen.garch.tse <- function(omega,
                               alpha,
                               beta = NULL,
                               model = c("sGARCH", "eGARCH", "gjrGARCH", 
                                         "apARCH", "iGARCH", "csGARCH"),
                               distribution = c("norm", "std", "ged", "snorm",
                                                "sstd", "sged", "nig", "jsu"),
                               distribution.params = list()) {
  
  
  if (!requireNamespace("rugarch", quietly = TRUE)) {
    stop("Package 'rugarch' is required. Install with: install.packages('rugarch')")
    
  }
  
  
  model <- match.arg(model)
  distribution <- match.arg(distribution)
  
  
  if (omega <= 0) {
    stop("omega must be positive")
  }
  if (any(alpha < 0)) {
    stop("alpha coefficients must be non-negative")
  }
  if (!is.null(beta) && any(beta < 0)) {
    stop("beta coefficients must be non-negative")
  }
  
  q <- length(alpha)
  p <- if (is.null(beta)) 0L else length(beta)
  
  if (model == "sGARCH" && (sum(alpha) + sum(beta)) >= 1) {
    warning("sum(alpha) + sum(beta) >= 1: process may be non-stationary")
  }
  
  names(alpha) <- paste0("alpha", seq_len(q))
  if (p > 0) {
    names(beta) <- paste0("beta", seq_len(p))
  }
  
  fixed <- c(omega = omega, alpha)
  if (p > 0) {
    fixed <- c(fixed, beta)
  }
  
  if (length(distribution.params) > 0) {
    fixed <- c(fixed, unlist(distribution.params))
  }
  
  force(model)
  force(distribution)
  force(q)
  force(p)
  force(fixed)
  
  function(n) {
    spec <- rugarch::ugarchspec(
      variance.model = list(model = model, garchOrder = c(q, p)),
      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
      distribution.model = distribution,
      fixed.pars = as.list(fixed)
    )
    
    path <- rugarch::ugarchpath(spec, n.sim = n)
    as.numeric(path@path$seriesSim)
  }
}