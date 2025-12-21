#' Generate a GARCH(p,q) Realization
#'
#' Generates a realization from a GARCH(p,q) model using fast C++ implementation.
#'
#' @param n Length of the output series
#' @param alpha0 Constant term (omega) in variance equation. Must be positive.
#' @param alpha Vector of ARCH coefficients (length q). These multiply lagged squared observations.
#' @param beta Vector of GARCH coefficients (length p). These multiply lagged conditional variances.
#' @param plot Logical. If TRUE, displays a plot of the realization.
#' @param sn Random seed for reproducibility. If 0, no seed is set.
#' @param burn Length of burn-in period. Default 500.
#'
#' @return An object of class "garch" containing:
#' \itemize{
#'   \item \code{y} - The generated GARCH series
#'   \item \code{sigma2} - Conditional variances
#'   \item \code{sigma} - Conditional standard deviations
#'   \item \code{eps} - The standard normal innovations used
#'   \item \code{alpha0} - The constant term
#'   \item \code{alpha} - ARCH coefficients
#'   \item \code{beta} - GARCH coefficients
#'   \item \code{p} - GARCH order (number of beta terms)
#'   \item \code{q} - ARCH order (number of alpha terms)
#'   \item \code{n} - Length of series
#'   \item \code{persistence} - Sum of alpha and beta (should be < 1 for stationarity)
#'   \item \code{uncond_var} - Unconditional variance
#'   \item \code{plot} - ggplot object
#' }
#'
#' @details
#' The GARCH(p,q) model is defined as:
#' \deqn{a_t = \sigma_t \epsilon_t}
#' \deqn{\sigma_t^2 = \alpha_0 + \sum_{i=1}^{q} \alpha_i a_{t-i}^2 + \sum_{j=1}^{p} \beta_j \sigma_{t-j}^2}
#' where \eqn{\epsilon_t \sim N(0,1)} i.i.d.
#'
#' For stationarity, require \eqn{\sum \alpha_i + \sum \beta_j < 1}.
#'
#' @examples
#' # Simple GARCH(1,1)
#' result <- gen_garch(n = 500, alpha0 = 0.1, alpha = 0.15, beta = 0.8)
#' 
#' # Access components
#' head(result$y)
#' head(result$sigma)
#' 
#' # Customize plot
#' result$plot + ggplot2::labs(title = "My GARCH Series")
#'
#' @export
gen_garch <- function(n, 
                      alpha0, 
                      alpha, 
                      beta, 
                      plot = FALSE, 
                      sn = 0,
                      burn = 500) {
  
  
  # Input validation
  if (n <= 0) stop("n must be a positive integer")
  if (alpha0 <= 0) stop("alpha0 must be positive")
  if (any(alpha < 0)) stop("alpha coefficients must be non-negative")
  if (any(beta < 0)) stop("beta coefficients must be non-negative")
  
  # Ensure alpha and beta are vectors
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  
  # Set seed if provided
  if (sn > 0) set.seed(sn)
  
  # Generate standard normal innovations
  ntot <- n + burn
  eps <- rnorm(ntot, mean = 0, sd = 1)
  
  # Call C++ function
  result <- garch_sim_cpp(n, alpha0, alpha, beta, eps, burn)
  
  # Calculate model properties
  p <- length(beta)
  q <- length(alpha)
  persistence <- sum(alpha) + sum(beta)
  
  uncond_var <- if (persistence < 1) {
    alpha0 / (1 - persistence)
  } else {
    NA  # undefined for IGARCH/explosive
  }
  
  # Build ggplot
  df <- data.frame(
    time = 1:n,
    y = result$a,
    sigma = result$sigma
  )
  
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = y)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.5) +
    ggplot2::labs(
      x = "Time", 
      y = "Value", 
      title = sprintf("GARCH(%d,%d) Realization", p, q)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (plot) print(gg)
  
  # Build output object
  out <- list(
    y = result$a,
    sigma2 = result$sigma2,
    sigma = result$sigma,
    eps = result$eps,
    alpha0 = alpha0,
    alpha = alpha,
    beta = beta,
    p = p,
    q = q,
    n = n,
    persistence = persistence,
    uncond_var = uncond_var,
    plot = gg
  )
  
  class(out) <- "garch"
  return(out)
}


#' @export
print.garch <- function(x, ...) {
  cat(sprintf("GARCH(%d,%d) Realization\n", x$p, x$q))
  cat(sprintf("n = %d\n", x$n))
  cat(sprintf("alpha0 (omega) = %.6f\n", x$alpha0))
  cat("alpha (ARCH) =", paste(round(x$alpha, 6), collapse = ", "), "\n")
  cat("beta (GARCH) =", paste(round(x$beta, 6), collapse = ", "), "\n")
  cat(sprintf("Persistence = %.4f", x$persistence))
  if (x$persistence >= 1) {
    cat(" (IGARCH/Explosive)\n")
  } else {
    cat(" (Stationary)\n")
  }
  if (!is.na(x$uncond_var)) {
    cat(sprintf("Unconditional variance = %.6f\n", x$uncond_var))
  }
  cat("\nFirst 6 values of y:\n")
  print(head(x$y))
  invisible(x)
}


#' @export
plot.garch <- function(x, type = c("series", "volatility", "both"), ...) {
  type <- match.arg(type)
  
  df <- data.frame(
    time = 1:x$n,
    y = x$y,
    sigma = x$sigma
  )
  
  if (type == "series") {
    print(x$plot)
  } else if (type == "volatility") {
    gg <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = sigma)) +
      ggplot2::geom_line(color = "firebrick", linewidth = 0.5) +
      ggplot2::labs(x = "Time", y = expression(sigma[t]), 
                    title = "Conditional Standard Deviation") +
      ggplot2::theme_minimal()
    print(gg)
  } else {
    # Both plots stacked
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = y)) +
      ggplot2::geom_line(color = "steelblue", linewidth = 0.5) +
      ggplot2::labs(x = "", y = "Value", title = sprintf("GARCH(%d,%d)", x$p, x$q)) +
      ggplot2::theme_minimal()
    
    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = sigma)) +
      ggplot2::geom_line(color = "firebrick", linewidth = 0.5) +
      ggplot2::labs(x = "Time", y = expression(sigma[t]), 
                    title = "Conditional Volatility") +
      ggplot2::theme_minimal()
    
    # Requires patchwork or gridExtra
    if (requireNamespace("patchwork", quietly = TRUE)) {
      print(p1 / p2)
    } else {
      print(p1)
      print(p2)
    }
  }
}