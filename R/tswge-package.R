#' tswge2: Time Series Analysis
#'
#' Tools for time series generation and analysis, including ARUMA models
#' with flexible innovation distributions (GARCH, t, etc.)
#'
#' @docType package
#' @name tswge2
#' @useDynLib tswge2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Fix ggplot2 NSE bindings (no visible binding for global variable)
utils::globalVariables(c("y", "sigma", "time"))