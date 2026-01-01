#' Woodward-Bottone-Gray Bootstrap Test for Trend
#'
#' Performs a bootstrap hypothesis test for trend in a time series while
#' accounting for autocorrelation. Under the null hypothesis, the series
#' is assumed to be a stationary AR process with no trend.
#' 
#' @importFrom stats sd
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @param x Numeric vector. The observed time series.
#' @param stat_fn Function. A function that takes a numeric vector and returns
#'   a single numeric test statistic. Use \code{\link{make.stat.co.tse}} or 
#'   create your own following the signature \code{function(x) -> numeric(1)}.
#' @param nb Integer. Number of bootstrap replicates. Default is 399.
#' @param p_max Integer. Maximum AR order to consider when fitting the null
#'   model. Order is selected by information criterion. Default is 5.
#' @param ar_method Character. Method for AR estimation of null model:
#'   \code{"mle"} (default) or \code{"burg"}. MLE is less biased near the
#'   unit circle but requires stationarity checking.
#' @param criterion Character. Information criterion for AR order selection:
#'   \code{"aic"} (default), \code{"aicc"}, or \code{"bic"}.
#' @param n_best Integer. For MLE, number of top models to check for
#'   stationarity before failing. Default is 3.
#' @param tol Numeric. Tolerance for stationarity check (roots must have
#'   modulus > tol). Default is 1.001.
#' @param bootadj Logical. If TRUE, performs the COBA variance adjustment.
#'   Default is TRUE. Consider setting to FALSE when using MLE, as the
#'   bias correction may be less necessary.
#' @param seed Integer. Random seed for reproducibility. If 0 (default), no
#'   seed is set. Used to generate \code{boot_seeds} if not provided.
#' @param boot_seeds Optional numeric vector of seeds for each bootstrap
#'   replication. If provided, must be of length \code{nb}. If \code{NULL}
#'   (default), seeds are auto-generated (after \code{seed} is set, if provided).
#'   Enables identical results in sequential and parallel modes.
#' @param boot_seeds_adj Optional numeric vector of seeds for COBA adjustment
#'   bootstrap. If \code{NULL} (default), auto-generated. Ignored if
#'   \code{bootadj = FALSE}.
#' @param parallel Logical. If TRUE, use parallel processing. Default is FALSE.
#' @param num_cpu Integer. Number of CPU cores to use. Default is 1.
#'   Set to 0 to use all available cores.
#' @param verbose Logical. If TRUE, display progress. Default is TRUE.
#'
#' @return A list containing:
#'   \item{obs_stat}{Observed test statistic.}
#'   \item{boot_dist}{Vector of bootstrap test statistics.}
#'   \item{pvalue_two}{Two-sided bootstrap p-value.}
#'   \item{pvalue_upper}{Upper one-sided p-value (H1: stat > 0).}
#'   \item{pvalue_lower}{Lower one-sided p-value (H1: stat < 0).}
#'   \item{ar_order}{AR order selected for null model.}
#'   \item{ar_phi}{AR coefficients for null model.}
#'   \item{ar_method}{AR estimation method used for null model.}
#'   \item{n}{Length of input series.}
#'   \item{nb}{Number of bootstrap replicates.}
#'   \item{boot_seeds}{Seeds used for bootstrap (for reproducibility).}
#'   \item{obs_stat_adj}{(If bootadj) Adjusted observed statistic.}
#'   \item{boot_dist_adj}{(If bootadj) Adjustment bootstrap distribution.}
#'   \item{pvalue_two_adj}{(If bootadj) Adjusted two-sided p-value.}
#'   \item{pvalue_upper_adj}{(If bootadj) Adjusted upper p-value.}
#'   \item{pvalue_lower_adj}{(If bootadj) Adjusted lower p-value.}
#'   \item{adj_factor}{(If bootadj) Variance adjustment factor.}
#'   \item{median_phi}{(If bootadj) AR coefficients of median model.}
#'   \item{boot_seeds_adj}{(If bootadj) Seeds used for adjustment bootstrap.}
#'
#' @references
#' Woodward, W. A., Bottone, S., and Gray, H. L. (1997). "Improved Tests for
#' Trend in Time Series Data." \emph{Journal of Agricultural, Biological, and
#' Environmental Statistics}, 2(4), 403-416.
#'
#' @seealso \code{\link{make.stat.co.tse}}, \code{\link{co.tse}}, 
#'   \code{\link{wbg.boot.tse}}, \code{\link{aic.ar.tse}}
#'
#' @examples
#' \dontrun{
#' # Modern approach: MLE estimation (less biased)
#' stat_fn <- make.stat.co.tse(maxp = 5, ar_method = "mle")
#' result <- wbg.boot.test.tse(x, stat_fn, nb = 399, ar_method = "mle")
#' 
#' # Try without COBA adjustment (MLE may not need it)
#' result_no_adj <- wbg.boot.test.tse(x, stat_fn, nb = 399, 
#'                                     ar_method = "mle", bootadj = FALSE)
#' 
#' # Original WBG paper approach: Burg estimation with COBA
#' stat_fn_burg <- make.stat.co.tse(maxp = 5, ar_method = "burg")
#' result_burg <- wbg.boot.test.tse(x, stat_fn_burg, nb = 399,
#'                                   ar_method = "burg", bootadj = TRUE)
#' 
#' # Reproduce exact results
#' result2 <- wbg.boot.test.tse(x, stat_fn, nb = 399, 
#'                               boot_seeds = result$boot_seeds,
#'                               boot_seeds_adj = result$boot_seeds_adj)
#' identical(result$pvalue_two, result2$pvalue_two)  # TRUE
#' }
#'
#' @export
wbg.boot.test.tse <- function(x, stat_fn, nb = 399, p_max = 5,
                              ar_method = c("mle", "burg"),
                              criterion = c("aic", "aicc", "bic"),
                              n_best = 3, tol = 1.001,
                              bootadj = TRUE,
                              seed = 0, boot_seeds = NULL, boot_seeds_adj = NULL,
                              parallel = FALSE, num_cpu = 1, verbose = TRUE) {
  
  ar_method <- match.arg(ar_method)
  criterion <- match.arg(criterion)
  
  if (parallel) {
    if (num_cpu == 0) {
      num_cpu <- parallel::detectCores()
    }
    num_cpu <- max(1, num_cpu)
  } else {
    num_cpu <- 1
  }
  use_parallel <- parallel && num_cpu > 1
  
  # Set master seed if provided (before generating boot_seeds)
  if (seed > 0) set.seed(seed)
  
  # Generate or validate boot_seeds
  if (is.null(boot_seeds)) {
    boot_seeds <- sample.int(.Machine$integer.max, nb)
  } else if (length(boot_seeds) != nb) {
    stop("boot_seeds must be of length nb (", nb, "), but got length ", 
         length(boot_seeds))
  }
  
  # Generate or validate boot_seeds_adj for COBA
  if (bootadj) {
    if (is.null(boot_seeds_adj)) {
      boot_seeds_adj <- sample.int(.Machine$integer.max, nb)
    } else if (length(boot_seeds_adj) != nb) {
      stop("boot_seeds_adj must be of length nb (", nb, "), but got length ", 
           length(boot_seeds_adj))
    }
  }
  
  n <- length(x)
  obs_stat <- stat_fn(x)
  
  # Fit AR model under H0: no trend, stationary AR
  if (verbose) message("Fitting null model (", ar_method, ")...")
  
  ar_fit <- aic.ar.tse(x, p_max = p_max, method = ar_method,
                       criterion = criterion, n_best = n_best, tol = tol, stationary = TRUE)
  ar_phi <- ar_fit$phi
  ar_p <- ar_fit$p
  ar_method_used <- ar_fit$method_used
  
  if (verbose) {
    message("  Selected AR(", ar_p, ") with phi = [", 
            paste(round(ar_phi, 4), collapse = ", "), "]")
  }
  
  # --- First Bootstrap ---
  if (verbose) message("Running first bootstrap (", nb, " replicates)...")
  
  if (use_parallel) {
    boot_results <- .run_parallel_bootstrap(
      nb = nb, n = n, phi = ar_phi, p_max = p_max,
      ar_method = ar_method, criterion = criterion,
      n_best = n_best, tol = tol,
      stat_fn = stat_fn, bootadj = bootadj, boot_seeds = boot_seeds,
      num_cpu = num_cpu, verbose = verbose
    )
  } else {
    boot_results <- .run_sequential_bootstrap(
      nb = nb, n = n, phi = ar_phi, p_max = p_max,
      ar_method = ar_method, criterion = criterion,
      n_best = n_best, tol = tol,
      stat_fn = stat_fn, bootadj = bootadj, boot_seeds = boot_seeds,
      verbose = verbose
    )
  }
  
  boot_stats <- sapply(boot_results, `[[`, "stat")
  
  # Compute p-values (with +1 correction)
  pvalue_two   <- (sum(abs(boot_stats) >= abs(obs_stat)) + 1) / (nb + 1)
  pvalue_upper <- (sum(boot_stats >= obs_stat) + 1) / (nb + 1)
  pvalue_lower <- (sum(boot_stats <= obs_stat) + 1) / (nb + 1)
  
  result <- list(
    obs_stat     = obs_stat,
    boot_dist    = boot_stats,
    pvalue_two   = pvalue_two,
    pvalue_upper = pvalue_upper,
    pvalue_lower = pvalue_lower,
    ar_order     = ar_p,
    ar_phi       = ar_phi,
    ar_method    = ar_method_used,
    n            = n,
    nb           = nb,
    boot_seeds   = boot_seeds
  )
  
  # --- COBA Adjustment (Second Bootstrap) ---
  if (bootadj) {
    boot_phi1 <- sapply(boot_results, `[[`, "phi1")
    boot_phi_list <- lapply(boot_results, `[[`, "phi")
    
    median_idx <- which.min(abs(boot_phi1 - median(boot_phi1)))
    median_phi <- boot_phi_list[[median_idx]]
    
    if (verbose) message("Running COBA adjustment bootstrap (", nb, " replicates)...")
    
    if (use_parallel) {
      boot_stats_adj <- .run_parallel_bootstrap_simple(
        nb = nb, n = n, phi = median_phi, stat_fn = stat_fn,
        boot_seeds = boot_seeds_adj, num_cpu = num_cpu, verbose = verbose
      )
    } else {
      boot_stats_adj <- .run_sequential_bootstrap_simple(
        nb = nb, n = n, phi = median_phi, stat_fn = stat_fn,
        boot_seeds = boot_seeds_adj, verbose = verbose
      )
    }
    
    adj_factor <- sd(boot_stats_adj) / sd(boot_stats)
    obs_stat_adj <- adj_factor * obs_stat
    
    pvalue_two_adj   <- (sum(abs(boot_stats) >= abs(obs_stat_adj)) + 1) / (nb + 1)
    pvalue_upper_adj <- (sum(boot_stats >= obs_stat_adj) + 1) / (nb + 1)
    pvalue_lower_adj <- (sum(boot_stats <= obs_stat_adj) + 1) / (nb + 1)
    
    result$obs_stat_adj     <- obs_stat_adj
    result$boot_dist_adj    <- boot_stats_adj
    result$pvalue_two_adj   <- pvalue_two_adj
    result$pvalue_upper_adj <- pvalue_upper_adj
    result$pvalue_lower_adj <- pvalue_lower_adj
    result$adj_factor       <- adj_factor
    result$median_phi       <- median_phi
    result$boot_seeds_adj   <- boot_seeds_adj
  }
  
  if (verbose) message("Done.")
  result
}


# --- Internal helper functions ---

.run_sequential_bootstrap <- function(nb, n, phi, p_max, ar_method, criterion,
                                      n_best, tol, stat_fn, bootadj, 
                                      boot_seeds, verbose) {
  if (verbose) pb <- txtProgressBar(min = 0, max = nb, style = 3)
  
  results <- lapply(1:nb, function(i) {
    xb <- gen.arma.wge(n, phi = phi, plot = FALSE, sn = boot_seeds[i])
    stat <- stat_fn(xb)
    
    if (verbose) setTxtProgressBar(pb, i)
    
    if (bootadj) {
      ar_boot <- tryCatch(
        aic.ar.tse(xb, p_max = p_max, method = ar_method,
                   criterion = criterion, n_best = n_best, tol = tol, stationary = TRUE),
        error = function(e) {
          # Fallback to burg if MLE fails
          aic.burg.wge(xb, p = 1:p_max)
        }
      )
      boot_phi <- if (!is.null(ar_boot$phi)) ar_boot$phi else ar_boot$phi
      list(stat = stat, phi = boot_phi, phi1 = 1 - sum(boot_phi))
    } else {
      list(stat = stat)
    }
  })
  
  if (verbose) close(pb)
  results
}

.run_sequential_bootstrap_simple <- function(nb, n, phi, stat_fn, boot_seeds, 
                                             verbose) {
  if (verbose) pb <- txtProgressBar(min = 0, max = nb, style = 3)
  
  results <- sapply(1:nb, function(i) {
    xb <- gen.arma.wge(n, phi = phi, plot = FALSE, sn = boot_seeds[i])
    stat <- stat_fn(xb)
    if (verbose) setTxtProgressBar(pb, i)
    stat
  })
  
  if (verbose) close(pb)
  results
}

.run_parallel_bootstrap <- function(nb, n, phi, p_max, ar_method, criterion,
                                    n_best, tol, stat_fn, bootadj,
                                    boot_seeds, num_cpu, verbose) {
  cl <- parallel::makeCluster(num_cpu)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(cl, c("stat_fn", "phi", "n", "p_max", "bootadj", 
                                "boot_seeds", "ar_method", "criterion",
                                "n_best", "tol"),
                          envir = environment())
  parallel::clusterEvalQ(cl, {
    library(tswge2)
    # Source the MLE functions if they're not in a package
    # This assumes aic.ar.wge and check.stationary are available
  })
  
  if (verbose && requireNamespace("pbapply", quietly = TRUE)) {
    results <- pbapply::pblapply(1:nb, function(i) {
      xb <- gen.arma.wge(n, phi = phi, plot = FALSE, sn = boot_seeds[i])
      stat <- stat_fn(xb)
      
      if (bootadj) {
        ar_boot <- tryCatch(
          aic.ar.tse(xb, p_max = p_max, method = ar_method,
                     criterion = criterion, n_best = n_best, tol = tol, stationary = TRUE),
          error = function(e) aic.burg.wge(xb, p = 1:p_max)
        )
        boot_phi <- ar_boot$phi
        list(stat = stat, phi = boot_phi, phi1 = 1 - sum(boot_phi))
      } else {
        list(stat = stat)
      }
    }, cl = cl)
  } else {
    if (verbose) message("(install 'pbapply' for parallel progress bars)")
    results <- parallel::parLapply(cl, 1:nb, function(i) {
      xb <- gen.arma.wge(n, phi = phi, plot = FALSE, sn = boot_seeds[i])
      stat <- stat_fn(xb)
      
      if (bootadj) {
        ar_boot <- tryCatch(
          aic.ar.tse(xb, p_max = p_max, method = ar_method,
                     criterion = criterion, n_best = n_best, tol = tol, stationary = TRUE),
          error = function(e) aic.burg.wge(xb, p = 1:p_max)
        )
        boot_phi <- ar_boot$phi
        list(stat = stat, phi = boot_phi, phi1 = 1 - sum(boot_phi))
      } else {
        list(stat = stat)
      }
    })
  }
  
  results
}

.run_parallel_bootstrap_simple <- function(nb, n, phi, stat_fn, boot_seeds,
                                           num_cpu, verbose) {
  cl <- parallel::makeCluster(num_cpu)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(cl, c("stat_fn", "phi", "n", "boot_seeds"),
                          envir = environment())
  parallel::clusterEvalQ(cl, library(tswge2))
  
  if (verbose && requireNamespace("pbapply", quietly = TRUE)) {
    results <- unlist(pbapply::pblapply(1:nb, function(i) {
      xb <- gen.arma.wge(n, phi = phi, plot = FALSE, sn = boot_seeds[i])
      stat_fn(xb)
    }, cl = cl))
  } else {
    results <- unlist(parallel::parLapply(cl, 1:nb, function(i) {
      xb <- gen.arma.wge(n, phi = phi, plot = FALSE, sn = boot_seeds[i])
      stat_fn(xb)
    }))
  }
  
  results
}