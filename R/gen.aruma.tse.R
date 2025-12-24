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
#' @param d Integer. Order of differencing (handled by arima.sim). Default is 0.
#' @param s Integer. Seasonal period for seasonal differencing. Default is 0.
#' @param lambda Numeric vector. Additional nonstationary factor coefficients.
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
#' @details
#' The function generates ARUMA realizations by:
#' \enumerate{
#'   \item Simulating an ARMA(p,q) process with optional differencing via
#'     \code{arima.sim}
#'   \item Applying the inverse of any additional nonstationary operators
#'     (seasonal factor and/or lambda)
#' }
#'
#' Note: \code{arima.sim} with \code{d > 0} returns \code{n + d} values
#' (the first \code{d} are prepended zeros from \code{diffinv}). This
#' function accounts for that offset when extracting the final series.
#'
#' Seasonal differencing \code{s} applies the factor \eqn{(1 - B^s)}.
#' The \code{lambda} parameter specifies additional nonstationary factor
#' coefficients that are multiplied with the seasonal factor.
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
#'
#' # Seasonal model with s = 12
#' result <- gen.aruma.tse(n = 500, phi = 0.5, s = 12, sn = 42)
gen.aruma.tse <- function(n, phi = 0, theta = 0, d = 0, s = 0, lambda = 0,
                          innov_gen = NULL, plot = TRUE, sn = NULL) {
  
  # --------------------------------------------------------------------------
  # Input validation
  # --------------------------------------------------------------------------
  
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("n must be a positive integer")
  }
  n <- as.integer(n)
  
  if (!is.numeric(d) || length(d) != 1 || d < 0) {
    stop("d must be a non-negative integer")
  }
  d <- as.integer(d)
  
  if (!is.numeric(s) || length(s) != 1 || s < 0) {
    stop("s must be a non-negative integer")
  }
  s <- as.integer(s)
  
  # --------------------------------------------------------------------------
  # Setup
  # --------------------------------------------------------------------------
  
  if (!is.null(sn)) set.seed(sn)
  
  if (is.null(innov_gen)) {
    innov_gen <- function(n) rnorm(n)
  }
  
  # --------------------------------------------------------------------------
  # Parse ARMA parameters (convert ATSA -> R sign convention for MA)
  # --------------------------------------------------------------------------
  
  p <- if (all(phi == 0)) 0L else length(phi)
  q <- if (all(theta == 0)) 0L else length(theta)
  
  # Build model list for arima.sim
  model <- list(order = c(p, d, q))
  if (p > 0) model$ar <- phi
  if (q > 0) model$ma <- -theta  # ATSA -> R sign convention
  
  # --------------------------------------------------------------------------
  # Parse nonstationary factors: seasonal (1 - B^s) and lambda
  # These are applied AFTER arima.sim via inverse operator
  # --------------------------------------------------------------------------
  
  has_lambda <- !all(lambda == 0)
  has_seasonal <- s > 0
  
  # Build coefficient vector for inverse nonstationary operator
  # Seasonal factor (1 - B^s) as AR-style coefficients: c(0,...,0,1) length s
  if (has_lambda && has_seasonal) {
    seas_coef <- c(rep(0, s - 1), 1)
    combined <- tswge2::mult.wge(fac1 = lambda, fac2 = seas_coef)
    ns_coef <- combined$model.coef
  } else if (has_lambda) {
    ns_coef <- lambda
  } else if (has_seasonal) {
    ns_coef <- c(rep(0, s - 1), 1)
  } else {
    ns_coef <- NULL
  }
  
  ns_order <- if (is.null(ns_coef)) 0L else length(ns_coef)
  
  # --------------------------------------------------------------------------
  # Calculate lengths and burn-in periods
  # --------------------------------------------------------------------------
  
  # n_start: burn-in for arima.sim (ARMA stabilization)
  # Larger for high-order or near-unit-root processes
  n_start <- max(100L, 10L * p, 10L * q)
  
  # spin: additional burn-in for inverse nonstationary operator
  # Need at least 2 full cycles of the longest nonstationary component
  spin <- max(100L, 2L * ns_order)
  
  # Total length to simulate (this is n for arima.sim, before un-differencing)
  # arima.sim will return n_sim + d values (d prepended by diffinv)
  n_sim <- n + ns_order + spin
  
  # --------------------------------------------------------------------------
  # Generate innovations and simulate ARMA/ARIMA
  # --------------------------------------------------------------------------
  
  # Single call to innov_gen preserves dependence structure (important for GARCH)
  total_innov_needed <- n_sim + n_start
  all_innov <- innov_gen(total_innov_needed)
  
  # arima.sim returns n_sim + d values when d > 0
  # First d values are prepended zeros from diffinv
  tsdata <- arima.sim(
    n = n_sim,
    model = model,
    innov = all_innov[(n_start + 1):total_innov_needed],
    n.start = n_start,
    start.innov = all_innov[1:n_start]
  )
  
  y <- as.numeric(tsdata)
  # y now has length n_sim + d
  
  # --------------------------------------------------------------------------
  # Apply inverse of nonstationary operator (if any)
  # --------------------------------------------------------------------------
  
  if (ns_order > 0) {
    # Apply cumulative recursion: x[t] = y[t] + sum(ns_coef[j] * x[t-j])
    # Pass d so helper can account for diffinv offset
    y_final <- .apply_inverse_ns_operator(y, ns_coef, d, spin, n)
  } else {
    # No nonstationary factors beyond d (which arima.sim handled)
    # Account for d prepended values from diffinv
    start_idx <- d + spin + 1
    y_final <- y[start_idx:(start_idx + n - 1)]
  }
  
  # --------------------------------------------------------------------------
  # Build output
  # --------------------------------------------------------------------------
  
  df_plot <- data.frame(time = seq_len(n), y = y_final)
  
  gg <- ggplot2::ggplot(df_plot, ggplot2::aes(x = time, y = y)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.5) +
    ggplot2::labs(x = "Time", y = "Value", title = "ARUMA Realization") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (isTRUE(plot)) print(gg)
  
  structure(
    list(
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
    ),
    class = "aruma"
  )
}


# ------------------------------------------------------------------------------
# Helper: Apply inverse nonstationary operator
# ------------------------------------------------------------------------------
#' @noRd
.apply_inverse_ns_operator <- function(y, ns_coef, d, spin, n) {
  ns_order <- length(ns_coef)
  len_y <- length(y)
  
  # Expected length: n + ns_order + spin + d (d from diffinv)
  expected_len <- n + ns_order + spin + d
  if (len_y != expected_len) {
    stop(sprintf(
      "Internal error: y has length %d, expected %d (n=%d, ns_order=%d, spin=%d, d=%d)",
      len_y, expected_len, n, ns_order, spin, d
    ))
  }
  
  # Allocate output vector (same length as input)
  x <- numeric(len_y)
  
  # Start recursion after d (diffinv offset) + ns_order (need prior values)
  start_recursion <- d + ns_order + 1
  
  # Apply cumulative sum recursion:
  # x[t] = y[t] + ns_coef[1]*x[t-1] + ns_coef[2]*x[t-2] + ... + ns_coef[k]*x[t-k]
  for (i in start_recursion:len_y) {
    x[i] <- y[i]
    for (j in seq_len(ns_order)) {
      x[i] <- x[i] + ns_coef[j] * x[i - j]
    }
  }
  
  # Extract final n observations after d offset and spin burn-in
  start_extract <- d + spin + ns_order + 1
  end_extract <- start_extract + n - 1
  
  if (end_extract > len_y) {
    stop(sprintf(
      "Internal error: extraction end %d exceeds array length %d",
      end_extract, len_y
    ))
  }
  
  x[start_extract:end_extract]
}


# ------------------------------------------------------------------------------
# S3 Methods
# ------------------------------------------------------------------------------

#' @export
print.aruma <- function(x, ...) {
  cat("ARUMA Realization\n")
  cat(sprintf("n = %d, p = %d, q = %d, d = %d, s = %d\n",
              x$n, x$p, x$q, x$d, x$s))
  
  if (x$p > 0) cat("phi =", paste(round(x$phi, 4), collapse = ", "), "\n")
  if (x$q > 0) cat("theta =", paste(round(x$theta, 4), collapse = ", "), "\n")
  if (!all(x$lambda == 0)) cat("lambda =", paste(round(x$lambda, 4), collapse = ", "), "\n")
  
  cat("\nFirst 6 values of y:\n")
  print(utils::head(x$y))
  
  invisible(x)
}


#' @export
plot.aruma <- function(x, ...) {
  print(x$plot)
  invisible(x$plot)
}