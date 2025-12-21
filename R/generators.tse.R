#' Create a GARCH Innovation Generator
#'
#' Factory function that returns a GARCH innovation generator 
#' compatible with gen.aruma.wge's innov_gen parameter.
#'
#' @param alpha0 Constant term (omega) in variance equation
#' @param alpha Vector of ARCH coefficients
#' @param beta Vector of GARCH coefficients
#'
#' @return A function with signature function(n) that generates GARCH innovations
#' @export
#'
#' @examples
#' # Create ARMA(1,1) with GARCH(1,1) errors
#' garch_gen <- make_garch_gen(alpha0 = 0.1, alpha = 0.15, beta = 0.8)
#' result <- gen.aruma.wge(n = 500, phi = 0.7, theta = 0.3, innov_gen = garch_gen)
make_garch_gen <- function(alpha0, alpha, beta) {
  function(n) {
    gen_garch(n, alpha0 = alpha0, alpha = alpha, beta = beta, plot = FALSE)$y
  }
}


#' Create a Student's t Innovation Generator
#'
#' Factory function that returns a t-distributed innovation generator.
#' The output is scaled to have variance 1 for comparability.
#'
#' @param df Degrees of freedom. Must be > 2 for finite variance.
#' @param scale If TRUE (default), scale to unit variance
#'
#' @return A function with signature function(n) that generates t innovations
#' @export
#'
#' @examples
#' # Create ARMA with heavy-tailed errors
#' t_gen <- make_t_gen(df = 5)
#' result <- gen.aruma.wge(n = 500, phi = 0.7, innov_gen = t_gen)
make_t_gen <- function(df, scale = TRUE) {
  if (df <= 2 && scale) {
    warning("df <= 2 has infinite variance; scale set to FALSE")
    scale <- FALSE
  }
  
  function(n) {
    x <- rt(n, df = df)
    if (scale && df > 2) {
      # Variance of t is df/(df-2), so scale to unit variance
      x <- x * sqrt((df - 2) / df)
    }
    x
  }
}


#' Create a Normal Innovation Generator
#'
#' Factory function that returns a normal innovation generator.
#'
#' @param sd Standard deviation (default 1)
#'
#' @return A function with signature function(n) that generates normal innovations
#' @export
#'
#' @examples
#' # Create ARMA with higher variance innovations
#' norm_gen <- make_norm_gen(sd = 2)
#' result <- gen.aruma.wge(n = 500, phi = 0.7, innov_gen = norm_gen)
make_norm_gen <- function(sd = 1) {
  function(n) {
    rnorm(n, mean = 0, sd = sd)
  }
}
