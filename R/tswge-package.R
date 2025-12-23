#' @keywords internal
"_PACKAGE"

# Export all functions that start with a letter (legacy behavior)
#' @rawNamespace exportPattern("^[[:alpha:]]+")

# ==============================================================================
# Base R Imports
# ==============================================================================

#' @importFrom grDevices gray hcl dev.new
#' @importFrom graphics abline axis hist image layout lines mtext par plot.new points segments title
#' @importFrom stats AIC BIC HoltWinters acf approx ar.burg ar.mle ar.yw arima arima.sim coef fft frequency lm na.pass pacf pchisq predict qnorm rnorm rt sigma spec.pgram time ts var vcov residuals
#' @importFrom utils head

# ==============================================================================
# Package Imports
# ==============================================================================

#' @import magrittr
#' @import zoo
#' @import ggplot2
#' @import dplyr
#' @import tidyverse
#' @import plotrix
#' @import forecast

#' @importFrom nnfor mlp
#' @importFrom astsa ARMAtoAR

# ==============================================================================
# GARCH Packages (rugarch, gt, cli)
# ==============================================================================

#' @importFrom rugarch ugarchspec ugarchfit ugarchpath infocriteria nyblom signbias
#' @importFrom WeightedPortTest Weighted.Box.test
#' @importFrom gt gt tab_header cols_hide cols_label tab_spanner fmt_number cols_align tab_style cell_text cells_body tab_source_note tab_options px pct html
#' @importFrom cli cli_h2 cli_text cli_rule cli_alert_warning col_green col_red col_silver style_bold

# ==============================================================================
# Global Variables (for ggplot2 NSE / R CMD check)
# ==============================================================================

utils::globalVariables(c(
  
  "y", "sigma", "time",
  # table.garch.gt.tse variables
  "ARCH", "GARCH", "AIC", "AICC", "BIC",
  "WLB1", "WLB2", "WLB3", "WLB1_fmt", "WLB2_fmt", "WLB3_fmt",
  "Nyblom", "Nyblom_crit", "Nyblom_fmt", "Nyblom_pass",
  "SignBias", "SignBias_fmt",
  "n_sig", "n_coef", "CoefSig", "CoefSig_pass",
  "Model",
  # table.coef.gt.tse variables
  "Parameter", "Estimate", "Std_Error", "t_value", "p_val", "p_display", "Sig"
))