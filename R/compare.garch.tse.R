#' Compare Multiple GARCH Model Specifications
#'
#' Fits a grid of ARCH/GARCH models to a time series and computes information
#' criteria and diagnostic statistics for model comparison and selection.
#'
#' @param data Numeric vector of returns or residuals to model.
#' @param arch.range Integer vector of ARCH orders to try (default 0:2).
#' @param garch.range Integer vector of GARCH orders to try (default 0:2).
#' @param distribution Character string specifying the innovation distribution.
#'   Default is "norm".
#' @param include.mean Logical. Include a mean term in the model?
#'   Default is FALSE.
#' @param solver Character string for optimization solver. Default is
#'   "hybrid" which tries multiple solvers.
#'
#' @return A list of class "garch.comparison.tse" containing fits (list of
#'   uGARCHfit objects), comparison (tibble with model statistics),
#'   distribution, and n.
#'
#' @details
#' The function fits all combinations of ARCH and GARCH orders (excluding
#' the trivial (0,0) case) and computes AIC, AICc, BIC, Weighted Ljung-Box
#' tests, Nyblom stability test, and Sign Bias test.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' garch_gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.15, beta = 0.8)
#' y <- garch_gen(1000)
#' results <- compare.garch.tse(y)
#' table.garch.gt.tse(results)
#' }
compare.garch.tse <- function(data,
                              arch.range = 0:2,
                              garch.range = 0:2,
                              distribution = "norm",
                              include.mean = FALSE,
                              solver = "hybrid") {
  
  if (!requireNamespace("rugarch", quietly = TRUE)) {
    stop("Package 'rugarch' is required. Install with: install.packages('rugarch')")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Install with: install.packages('dplyr')")
  }
  
  # Input validation
  if (!is.numeric(data) || length(data) < 10) {
    stop("data must be a numeric vector with at least 10 observations")
  }
  
  if (any(is.na(data))) {
    warning("NA values detected in data; these may cause fitting issues")
  }
  
  # Build grid of orders, excluding (0,0)
  orders <- expand.grid(arch = arch.range, garch = garch.range)
  orders <- orders[!(orders$arch == 0 & orders$garch == 0), ]
  
  if (nrow(orders) == 0) {
    stop("No valid model orders to fit (need at least one non-zero order)")
  }
  
  results <- list()
  fits <- list()
  n <- length(data)
  
  for (i in seq_len(nrow(orders))) {
    arch_order <- orders$arch[i]
    garch_order <- orders$garch[i]
    
    model_label <- .garch_label(arch_order, garch_order)
    
    tryCatch({
      spec <- rugarch::ugarchspec(
        variance.model = list(model = "sGARCH", garchOrder = c(arch_order, garch_order)),
        mean.model = list(armaOrder = c(0, 0), include.mean = include.mean),
        distribution.model = distribution
      )
      
      fit <- rugarch::ugarchfit(spec, data, solver = solver)
      
      diag <- .extract_garch_diagnostics(fit, arch_order, garch_order, n)
      
      fits[[model_label]] <- fit
      
      results[[length(results) + 1]] <- dplyr::tibble(
        Model = model_label,
        ARCH = arch_order,
        GARCH = garch_order,
        AIC = diag$aic,
        AICC = diag$aicc,
        BIC = diag$bic,
        WLB1 = diag$wlb[1],
        WLB2 = diag$wlb[2],
        WLB3 = diag$wlb[3],
        Nyblom = diag$nyblom_stat,
        Nyblom_crit = diag$nyblom_crit,
        SignBias = diag$signbias_pval,
        n_sig = diag$n_sig,
        n_coef = diag$n_coef
      )
      
    }, error = function(e) {
      message(paste0(model_label, " failed: ", e$message))
    })
  }
  
  if (length(results) == 0) {
    stop("All models failed to converge. Check your data or try different orders.")
  }
  
  comparison <- dplyr::bind_rows(results)
  
  structure(
    list(
      fits = fits,
      comparison = comparison,
      distribution = distribution,
      n = n
    ),
    class = "garch.comparison.tse"
  )
}


# ==============================================================================
# Helper Functions
# ==============================================================================

#' @noRd
.garch_label <- function(arch_order, garch_order) {
  if (garch_order == 0) {
    paste0("ARCH(", arch_order, ")")
  } else {
    paste0("GARCH(", arch_order, ",", garch_order, ")")
  }
}


#' @noRd
.extract_garch_diagnostics <- function(fit, arch_order, garch_order, n) {
  
  std_resid <- as.numeric(rugarch::residuals(fit, standardize = TRUE))
  k <- arch_order + garch_order + 1
  
  # Information criteria
  ic <- rugarch::infocriteria(fit)
  aic <- ic["Akaike", ]
  bic <- ic["Bayes", ]
  aicc <- aic + (2 * k^2 + 2 * k) / (n - k - 1)
  
  # Weighted Ljung-Box on squared standardized residuals
  df_garch <- arch_order + garch_order
  wlb <- tryCatch({
    box1 <- WeightedPortTest::Weighted.Box.test(std_resid, lag = 1, 
                                                type = "Ljung-Box", fitdf = 0, sqrd.res = TRUE)
    box2 <- WeightedPortTest::Weighted.Box.test(std_resid, 
                                                lag = max(2, 3 * df_garch - 1), 
                                                type = "Ljung-Box", fitdf = df_garch, sqrd.res = TRUE)
    box3 <- WeightedPortTest::Weighted.Box.test(std_resid, 
                                                lag = max(5, 5 * df_garch - 1), 
                                                type = "Ljung-Box", fitdf = df_garch, sqrd.res = TRUE)
    c(box1$p.value, box2$p.value, box3$p.value)
  }, error = function(e) {
    c(NA, NA, NA)
  })
  wlb_pvals <- wlb
  
  # Nyblom stability test
  nyb <- tryCatch({
    rugarch::nyblom(fit)
  }, error = function(e) {
    list(JointStat = NA, JointCritical = c("5%" = NA))
  })
  nyblom_stat <- nyb$JointStat
  nyblom_crit <- nyb$JointCritical["5%"]
  
  # Sign bias test
  sb <- tryCatch({
    rugarch::signbias(fit)
  }, error = function(e) {
    matrix(c(NA, NA, NA, NA, NA, NA, NA, NA), nrow = 4, ncol = 2)
  })
  signbias_pval <- sb[4, 2]
  
  # Coefficient significance
  coef_mat <- fit@fit$matcoef
  n_coef <- nrow(coef_mat)
  n_sig <- sum(coef_mat[, 4] < 0.05)
  
  list(
    aic = aic,
    aicc = aicc,
    bic = bic,
    wlb = wlb_pvals,
    nyblom_stat = nyblom_stat,
    nyblom_crit = nyblom_crit,
    signbias_pval = signbias_pval,
    n_sig = n_sig,
    n_coef = n_coef
  )
}


#' Print method for garch.comparison.tse
#' @param x A garch.comparison.tse object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns x
#' @exportS3Method
print.garch.comparison.tse <- function(x, ...) {
  cat("GARCH Model Comparison\n")
  cat("----------------------\n")
  cat("Distribution:", x$distribution, "\n")
  cat("Sample size:", x$n, "\n")
  cat("Models fitted:", length(x$fits), "\n\n")
  
  # Find best by each criterion
  comp <- x$comparison
  best_aic <- comp$Model[which.min(comp$AIC)]
  best_bic <- comp$Model[which.min(comp$BIC)]
  best_aicc <- comp$Model[which.min(comp$AICC)]
  
  cat("Best by AIC: ", best_aic, "\n")
  cat("Best by BIC: ", best_bic, "\n")
  cat("Best by AICc:", best_aicc, "\n")
  cat("\nUse table.garch.gt.tse() or table.garch.cli.tse() for full comparison.\n")
  
  invisible(x)
}


#' Summary method for garch.comparison.tse
#' @param object A garch.comparison.tse object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns object
#' @exportS3Method
summary.garch.comparison.tse <- function(object, ...) {
  print(object$comparison)
  invisible(object)
}