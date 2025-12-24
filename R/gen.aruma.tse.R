#' Generate ARUMA Realizations with Flexible Innovations
#'
#' Generates realizations from ARUMA (ARIMA with possible seasonal and
#' nonstationary components) models with custom innovation distributions.
#'
#' @param n Integer. Length of the realization to generate.
#' @param phi Numeric vector. AR parameters (ATSA sign convention).
#'   Default is 0 (no AR component).
#' @param theta Numeric vector. MA parameters (ATSA sign convention).
#'   Default is 0 (no MA component).
#' @param d Integer. Order of differencing. Default is 0.
#' @param s Integer. Order of seasonal differencing. Default is 0.
#' @param lambda Numeric vector. Parameters in nonstationary operator.
#'   Default is 0.
#' @param innov_gen Function. Innovation generator with signature
#'   \code{function(n)} returning a numeric vector of length n with mean zero.
#'   If \code{NULL} (default), uses standard normal innovations.
#' @param plot Logical. If \code{TRUE} (default), plot the realization.
#' @param sn Integer or NULL. Random seed for reproducibility. 
#'   Default is NULL (no seed set).
#'
#' @return An object of class \code{"aruma"} containing:
#'   \describe{
#'     \item{y}{Numeric vector of the generated realization}
#'     \item{n, p, q, d, s}{Model orders}
#'     \item{phi, theta, lambda}{Model parameters}
#'     \item{plot}{A ggplot object of the realization}
#'   }
#'
#' @seealso \code{\link{make.gen.garch.tse}}, \code{\link{make.gen.t.tse}},
#'   \code{\link{make.gen.norm.tse}}
#'
#' @export
#' @import ggplot2
#' @importFrom stats arima.sim rnorm
#' @importFrom utils head
#'
#' @examples
#' # Simple AR(1) with normal innovations
#' result <- gen.aruma.tse(n = 200, phi = 0.7, sn = 123)
#'
#' # ARMA(1,1) with GARCH(1,1) innovations
#' garch_gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.15, beta = 0.8)
#' result <- gen.aruma.tse(n = 500, phi = 0.7, theta = 0.3,
#'                         innov_gen = garch_gen, sn = 42)
#'
#' # ARIMA(1,1,1) with t-distributed innovations
#' t_gen <- make.gen.t.tse(df = 5)
#' result <- gen.aruma.tse(n = 500, phi = 0.7, theta = 0.3, d = 1,
#'                         innov_gen = t_gen, sn = 42)
gen.aruma.tse <- function(n, phi = 0, theta = 0, d = 0, s = 0, lambda = 0,
                          innov_gen = NULL, plot = TRUE, sn = NULL) {
  # Set seed if provided (before any random generation)
  if (!is.null(sn)) {
    set.seed(sn)
  }
  
  # Default innovation generator: standard normal
  if (is.null(innov_gen)) {
    innov_gen <- function(n) rnorm(n, mean = 0, sd = 1)
  }
  
  # Convert parameters for arima.sim
  # arima.sim uses opposite sign convention for MA
  ar <- phi
  ma <- -theta
  p <- length(ar)
  q <- length(ma)
  dlam <- length(lambda)
  
  if (all(ar == 0)) {
    ar <- NA
    p <- 0
  }
  if (all(ma == 0)) {
    ma <- NA
    q <- 0
  }
  if (all(lambda == 0)) {
    lambda <- NA
    dlam <- 0
  }
  
  # Set up seasonal component
  dlams <- dlam + s
  lambdas <- rep(0, 100)
  seas <- rep(0, 100)
  if (s > 0) seas[s] <- 1
  
  # Calculate total length needed and burn-in
  n.start <- max(100L, 10L * p, 10L * q)
  spin <- 100
  ngen <- n + dlams + spin
  
  # Generate all innovations from custom generator for both burn-in and main series
  total_innov <- innov_gen(ngen + n.start)
  start_innov <- total_innov[1:n.start]
  main_innov <- total_innov[(n.start + 1):(n.start + ngen)]
  
  # Simulate ARMA process with consistent innovations throughout
  if ((p > 0) && (q > 0)) {
    tsdata <- arima.sim(n = ngen, model = list(order = c(p, d, q), ar = ar, ma = ma),
                        innov = main_innov, n.start = n.start, start.innov = start_innov)
  }
  if ((p == 0) && (q > 0)) {
    tsdata <- arima.sim(n = ngen, model = list(order = c(p, d, q), ma = ma),
                        innov = main_innov, n.start = n.start, start.innov = start_innov)
  }
  if ((p > 0) && (q == 0)) {
    tsdata <- arima.sim(n = ngen, model = list(order = c(p, d, q), ar = ar),
                        innov = main_innov, n.start = n.start, start.innov = start_innov)
  }
  if ((p == 0) && (q == 0)) {
    tsdata <- arima.sim(n = ngen, model = list(order = c(0, d, 0)),
                        innov = main_innov, n.start = n.start, start.innov = start_innov)
  }
  
  # Compute the inverse of the nonstationary operator
  y <- as.numeric(tsdata)
  
  if ((dlam > 0) && (s > 0)) {
    temp <- tswge2::mult.wge(fac1 = lambda, fac2 = seas)
    lambdas <- temp$model.coef
  }
  if ((dlam > 0) && (s == 0)) {
    lambdas <- lambda
  }
  if ((dlam == 0) && (s > 0)) {
    lambdas <- seas
  }
  if ((dlam == 0) && (s == 0)) {
    lambdas <- 0
  }
  
  d1 <- d + dlams + 1
  nd <- n + d + dlams + 1
  ndspin <- nd + spin - 1
  xfull <- rep(0, ndspin)
  x <- rep(0, ndspin)
  
  if (dlams == 0) {
    for (i in d1:ndspin) {
      xfull[i] <- y[i]
    }
  }
  
  if (dlams > 0) {
    for (i in d1:ndspin) {
      xfull[i] <- y[i]
      for (j in 1:dlams) {
        xfull[i] <- xfull[i] + lambdas[j] * xfull[i - j]
      }
    }
  }
  
  for (ii in 1:n) {
    x[ii] <- xfull[ii + spin + d1 - 1]
  }
  
  # Final series
  y_final <- as.numeric(x[1:n])
  
  # Build ggplot object
  df <- data.frame(time = 1:n, y = y_final)
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = y)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.5) +
    ggplot2::labs(x = "Time", y = "Value", title = "ARUMA Realization") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Show plot if requested
  if (plot == TRUE) {
    print(gg)
  }
  
  # Build result object
  result <- list(
    y = y_final,
    n = n,
    p = p,
    q = q,
    d = d,
    s = s,
    phi = phi,
    theta = theta,
    lambda = lambda,
    plot = gg
  )
  
  class(result) <- "aruma"
  
  return(result)
}


#' Print Method for aruma Objects
#'
#' @param x An object of class \code{"aruma"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.aruma <- function(x, ...) {
  cat("ARUMA Realization\n")
  cat(sprintf("n = %d, p = %d, q = %d, d = %d, s = %d\n",
              x$n, x$p, x$q, x$d, x$s))
  
  if (x$p > 0) {
    cat("phi =", paste(round(x$phi, 4), collapse = ", "), "\n")
  }
  if (x$q > 0) {
    cat("theta =", paste(round(x$theta, 4), collapse = ", "), "\n")
  }
  if (!all(is.na(x$lambda))) {
    cat("lambda =", paste(round(x$lambda, 4), collapse = ", "), "\n")
  }
  
  cat("\nFirst 6 values of y:\n")
  print(utils::head(x$y))
  
  invisible(x)
}


#' Plot Method for aruma Objects
#'
#' @param x An object of class \code{"aruma"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the ggplot object.
#'
#' @export
plot.aruma <- function(x, ...) {
  print(x$plot)
  invisible(x$plot)
}