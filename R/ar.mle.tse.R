#' Check if AR Model is Stationary
#'
#' Checks whether an AR model with given coefficients is stationary by
#' verifying all roots of the characteristic polynomial lie outside the
#' unit circle.
#'
#' @param phi Numeric vector of AR coefficients.
#' @param tol Numeric. Roots must have modulus > tol to be considered
#'   outside the unit circle. Default is 1.001 to provide a small buffer.
#'
#' @return Logical. TRUE if stationary, FALSE otherwise.
#'
#' @examples
#' \dontrun{
#' check.stationary.tse(0.9)        # TRUE (stationary AR(1))
#' check.stationary.tse(1.01)       # FALSE (non-stationary)
#' check.stationary.tse(c(0.5, 0.3)) # TRUE (stationary AR(2))
#' }
#'
#' @export
check.stationary.tse <- function(phi, tol = 1.001) {
  if (is.null(phi) || length(phi) == 0) return(TRUE)  # white noise
  
  # Characteristic polynomial: 1 - phi[1]*z - phi[2]*z^2 - ...
  poly_coefs <- c(1, -phi)
  roots <- polyroot(poly_coefs)
  
  # All roots must be strictly outside unit circle
  all(Mod(roots) > tol)
}


#' Fit AR Model via MLE with AIC/AICC/BIC Selection
#'
#' Fits AR models of orders 1 through p_max using maximum likelihood
#' estimation, selects the best model according to the specified information
#' criterion. Optionally verifies stationarity and falls back through top 
#' candidates if the best model is non-stationary.
#'
#' @param x Numeric vector. The time series to fit.
#' @param p_max Integer. Maximum AR order to consider. Can also be a vector
#'   of orders to try (e.g., 1:5).
#' @param criterion Character. Information criterion for model selection:
#'   \code{"aic"}, \code{"aicc"} (corrected AIC), or \code{"bic"}. 
#'   Default is \code{"aic"}.
#' @param stationary Logical. If TRUE (default), check that the selected model
#'   is stationary and fall back to next best if not. If FALSE, return the
#'   best model by criterion regardless of stationarity.
#' @param n_best Integer. Number of top models (by criterion) to check for
#'   stationarity before failing. Ignored if \code{stationary = FALSE}. 
#'   Default is 3.
#' @param tol Numeric. Tolerance for stationarity check (roots must have
#'   modulus > tol). Ignored if \code{stationary = FALSE}. Default is 1.001.
#' @param silent Logical. If TRUE, suppress convergence warnings from arima().
#'   Default is TRUE.
#'
#' @return A list containing:
#'   \item{phi}{Numeric vector of AR coefficients.}
#'   \item{p}{Integer. The selected AR order.}
#'   \item{criterion_value}{The value of the selection criterion.}
#'   \item{method_used}{Character. Always "mle" for this function.}
#'   \item{all_criteria}{Named numeric vector of criterion values for all 
#'     successfully fit orders.}
#'
#' @details
#' The function attempts to fit AR(p) models for each order in 1:p_max using
#' \code{arima()} with \code{method = "CSS-ML"} (conditional sum of squares
#' initialization followed by MLE refinement). 
#'
#' Models are ranked by the specified criterion. If \code{stationary = TRUE},
#' the function checks the top \code{n_best} models for stationarity and 
#' returns the best stationary model. If no stationary model is found among 
#' the top candidates, an error is raised.
#'
#' If \code{stationary = FALSE}, the best model by criterion is returned
#' regardless of whether it is stationary.
#'
#' @seealso \code{\link{check.stationary.tse}}, \code{\link{aic.ar.tse}}
#'
#' @examples
#' \dontrun{
#' x <- gen.arma.wge(100, phi = c(0.7, 0.2), sn = 123)
#' fit <- aic.ar.mle.tse(x, p_max = 5, criterion = "bic")
#' print(fit$phi)
#' print(fit$p)
#'
#' # Allow non-stationary estimates (for simulation studies)
#' fit_unrestricted <- aic.ar.mle.tse(x, p_max = 5, stationary = FALSE)
#' }
#'
#' @export
aic.ar.mle.tse <- function(x, p_max, criterion = c("aic", "aicc", "bic"),
                           stationary = TRUE, n_best = 3, tol = 1.001, 
                           silent = TRUE) {
  
  criterion <- match.arg(criterion)
  n <- length(x)
  
  # Handle p_max as either single value or vector
  if (length(p_max) == 1) {
    orders <- 1:p_max
  } else {
    orders <- p_max
  }
  orders <- orders[orders >= 1]  # ensure positive orders
  
  if (length(orders) == 0) {
    stop("No valid AR orders specified")
  }
  
  # Fit models for each order
  results <- list()
  all_criteria <- c()
  
  for (p in orders) {
    fit <- tryCatch({
      if (silent) {
        suppressWarnings(
          arima(x, order = c(p, 0, 0), method = "CSS-ML", include.mean = TRUE)
        )
      } else {
        arima(x, order = c(p, 0, 0), method = "CSS-ML", include.mean = TRUE)
      }
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      # Extract AR coefficients
      coef_names <- names(fit$coef)
      ar_idx <- grep("^ar", coef_names)
      phi <- if (length(ar_idx) > 0) as.numeric(fit$coef[ar_idx]) else numeric(0)
      
      # Calculate information criteria
      k <- p + 2  # AR coefficients + intercept + variance
      loglik <- fit$loglik
      
      aic_val <- -2 * loglik + 2 * k
      aicc_val <- aic_val + (2 * k * (k + 1)) / max(n - k - 1, 1)
      bic_val <- -2 * loglik + k * log(n)
      
      crit_val <- switch(criterion,
                         aic = aic_val,
                         aicc = aicc_val,
                         bic = bic_val)
      
      results[[length(results) + 1]] <- list(
        p = p,
        phi = phi,
        criterion_value = crit_val,
        aic = aic_val,
        aicc = aicc_val,
        bic = bic_val
      )
      
      all_criteria[as.character(p)] <- crit_val
    }
  }
  
  if (length(results) == 0) {
    stop("MLE estimation failed for all AR orders. ",
         "The series may be too short or have numerical issues.")
  }
  
  # Sort by criterion (ascending - lower is better)
  crit_vals <- sapply(results, `[[`, "criterion_value")
  sorted_idx <- order(crit_vals)
  
  # If not enforcing stationarity, return best model by criterion
  if (!stationary) {
    best <- results[[sorted_idx[1]]]
    return(list(
      phi = best$phi,
      p = best$p,
      criterion_value = best$criterion_value,
      method_used = "mle",
      all_criteria = all_criteria
    ))
  }
  
  # Check top N for stationarity
  n_check <- min(n_best, length(results))
  
  for (i in seq_len(n_check)) {
    candidate <- results[[sorted_idx[i]]]
    if (check.stationary.tse(candidate$phi, tol = tol)) {
      return(list(
        phi = candidate$phi,
        p = candidate$p,
        criterion_value = candidate$criterion_value,
        method_used = "mle",
        all_criteria = all_criteria
      ))
    }
  }
  
  stop("No stationary model found among top ", n_best, " candidates by ",
       toupper(criterion), ". Consider increasing n_best or p_max, ",
       "or use method = 'burg' which guarantees stationarity.")
}


#' Unified AR Model Fitting
#'
#' Fits an AR model using either MLE or Burg estimation with automatic
#' order selection via information criteria.
#'
#' @param x Numeric vector. The time series to fit.
#' @param p_max Integer. Maximum AR order to consider.
#' @param method Character. Estimation method: \code{"mle"} or \code{"burg"}.
#'   Default is \code{"mle"}.
#' @param criterion Character. Information criterion for model selection:
#'   \code{"aic"}, \code{"aicc"}, or \code{"bic"}.
#'   Default is \code{"aic"}.
#' @param stationary Logical. If TRUE (default), check that MLE models are
#'   stationary. If FALSE, return best model by criterion regardless.
#'   Ignored for Burg method (always stationary).
#' @param n_best Integer. For MLE, number of top models to check for
#'   stationarity. Ignored for Burg or if \code{stationary = FALSE}. Default is 3.
#' @param tol Numeric. Tolerance for stationarity check. Default is 1.001.
#' @param silent Logical. Suppress warnings. Default is TRUE.
#'
#' @return A list containing:
#'   \item{phi}{Numeric vector of AR coefficients.}
#'   \item{p}{Integer. The selected AR order.}
#'   \item{method_used}{Character. The method actually used.}
#'
#' @details
#' When \code{method = "burg"}, uses \code{aic.burg.wge()} from the tswge
#' package, which guarantees a stationary model but may have bias near
#' the unit circle.
#'
#' When \code{method = "mle"}, uses \code{aic.ar.mle.tse()} which provides
#' less biased estimates but requires stationarity checking (unless
#' \code{stationary = FALSE}).
#'
#' @seealso \code{\link{aic.ar.mle.tse}}, \code{\link[tswge]{aic.burg.wge}}
#'
#' @examples
#' \dontrun{
#' x <- gen.arma.wge(100, phi = 0.95, sn = 123)
#' 
#' # Compare methods
#' fit_burg <- aic.ar.tse(x, p_max = 5, method = "burg")
#' fit_mle <- aic.ar.tse(x, p_max = 5, method = "mle")
#' 
#' cat("Burg phi:", fit_burg$phi, "\n")
#' cat("MLE phi:", fit_mle$phi, "\n")
#' 
#' # Unrestricted MLE for simulation studies
#' fit_mle_free <- aic.ar.tse(x, p_max = 5, method = "mle", stationary = FALSE)
#' }
#'
#' @export
aic.ar.tse <- function(x, p_max, method = c("mle", "burg"),
                       criterion = c("aic", "aicc", "bic"),
                       stationary = TRUE, n_best = 3, tol = 1.001, 
                       silent = TRUE) {
  
  method <- match.arg(method)
  criterion <- match.arg(criterion)
  
  if (method == "burg") {
    # Use tswge's Burg estimation with specified criterion
    fit <- aic.burg.wge(x, p = 1:p_max, type = criterion)
    return(list(
      phi = fit$phi,
      p = fit$p,
      method_used = "burg"
    ))
  }
  
  # MLE method
  aic.ar.mle.tse(x, p_max = p_max, criterion = criterion,
                 stationary = stationary, n_best = n_best, tol = tol, 
                 silent = silent)
}