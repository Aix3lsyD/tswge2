#' @title GARCH Model Display Functions
#' @description Functions for displaying GARCH model comparison results
#'   and coefficient tables in various formats.
#' @name display
#' @keywords internal
"_PACKAGE"


# ==============================================================================
# Shared Utility Functions
# ==============================================================================

#' @noRd
.fmt_pval <- function(p, digits = 3) {
  if (is.na(p)) {
    return("NA")
  }
  if (p < 0.001) {
    "<0.001"
  } else {
    sprintf(paste0("%.", digits, "f"), p)
  }
}


#' @noRd
.find_best_indices <- function(df) {
  list(
    aic = which.min(df$AIC),
    bic = which.min(df$BIC),
    aicc = which.min(df$AICC),
    wlb1 = which.max(df$WLB1),
    wlb2 = which.max(df$WLB2),
    wlb3 = which.max(df$WLB3),
    nyblom = which.min(df$Nyblom),
    signbias = which.max(df$SignBias)
  )
}


# ==============================================================================
# gt Table Display
# ==============================================================================

#' Display GARCH Comparison Table with gt
#'
#' Creates a publication-ready comparison table using the \pkg{gt} package
#' with color-coded highlighting for best values and diagnostic failures.
#'
#' @param results A \code{garch.comparison.tse} object from
#'   \code{\link{compare.garch.tse}}.
#' @param title Character string for table title. Default is
#'   "GARCH Model Comparison".
#' @param color.best Color for best/passing values. Default is "seagreen".
#' @param color.fail Color for failing values. Default is "tomato".
#'
#' @return A \code{gt} table object.
#'
#' @details
#' The table displays:
#' \itemize{
#'   \item \strong{Information Criteria}: AIC, AICc, BIC (lower is better,
#'     minimum highlighted in green)
#'   \item \strong{Weighted Ljung-Box}: P-values at lags 1-3 (higher is better,
#'     red if < 0.05)
#'   \item \strong{Nyblom}: Stability test statistic (lower is better,
#'     red if >= critical value)
#'   \item \strong{Sign Bias}: Joint test p-value (higher is better,
#'     red if < 0.05)
#'   \item \strong{Coef}: Significant/total coefficients (green if all
#'     significant)
#' }
#'
#' @seealso \code{\link{compare.garch.tse}}, \code{\link{table.garch.cli.tse}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' garch_gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.15, beta = 0.8)
#' y <- garch_gen(1000)
#' results <- compare.garch.tse(y)
#'
#' # Display table
#' table.garch.gt.tse(results)
#'
#' # Custom colors
#' table.garch.gt.tse(results, color.best = "darkgreen", color.fail = "red")
#' }
table.garch.gt.tse <- function(results,
                               title = "GARCH Model Comparison",
                               color.best = "seagreen",
                               color.fail = "tomato") {
  
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required. Install with: install.packages('gt')")
  }
  
  if (!inherits(results, "garch.comparison.tse")) {
    stop("results must be a 'garch.comparison.tse' object from compare.garch.tse()")
  }
  
  df <- results$comparison
  
  if (nrow(df) == 0) {
    stop("No models to display (comparison table is empty)")
  }
  
  # Build display dataframe with formatted columns
  display_df <- df
  display_df$WLB1_fmt <- sapply(df$WLB1, .fmt_pval)
  display_df$WLB2_fmt <- sapply(df$WLB2, .fmt_pval)
  display_df$WLB3_fmt <- sapply(df$WLB3, .fmt_pval)
  display_df$Nyblom_fmt <- sprintf("%.2f", df$Nyblom)
  display_df$Nyblom_pass <- df$Nyblom < df$Nyblom_crit
  display_df$SignBias_fmt <- sapply(df$SignBias, .fmt_pval)
  display_df$CoefSig <- paste0(df$n_sig, "/", df$n_coef)
  display_df$CoefSig_pass <- df$n_sig == df$n_coef
  
  tbl <- display_df |>
    gt::gt() |>
    gt::tab_header(
      title = title,
      subtitle = paste0("Distribution: ", results$distribution)
    ) |>
    gt::cols_hide(columns = c(
      ARCH, GARCH, WLB1, WLB2, WLB3, Nyblom, Nyblom_crit,
      Nyblom_pass, SignBias, n_sig, n_coef, CoefSig_pass
    )) |>
    gt::cols_label(
      Model = "Model",
      AIC = "AIC",
      AICC = "AICc",
      BIC = "BIC",
      WLB1_fmt = "Lag 1",
      WLB2_fmt = "Lag 2",
      WLB3_fmt = "Lag 3",
      Nyblom_fmt = "Nyblom",
      SignBias_fmt = "Sign Bias",
      CoefSig = "Coef"
    ) |>
    gt::tab_spanner(label = "Info Criteria", columns = c(AIC, AICC, BIC)) |>
    gt::tab_spanner(label = "Weighted Ljung-Box", columns = c(WLB1_fmt, WLB2_fmt, WLB3_fmt)) |>
    gt::tab_spanner(label = "Diagnostics", columns = c(Nyblom_fmt, SignBias_fmt)) |>
    gt::fmt_number(columns = c(AIC, AICC, BIC), decimals = 2) |>
    gt::cols_align(align = "right", columns = c(AIC, AICC, BIC, WLB1_fmt, WLB2_fmt, WLB3_fmt, Nyblom_fmt, SignBias_fmt)) |>
    gt::cols_align(align = "center", columns = CoefSig) |>
    gt::cols_align(align = "left", columns = Model) |>
    gt::tab_style(style = gt::cell_text(weight = "bold"), locations = gt::cells_body(columns = Model))
  
  # Highlight best IC (green bold)
  tbl <- tbl |>
    gt::tab_style(
      style = gt::cell_text(color = color.best, weight = "bold"),
      locations = gt::cells_body(columns = AIC, rows = AIC == min(AIC))
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = color.best, weight = "bold"),
      locations = gt::cells_body(columns = AICC, rows = AICC == min(AICC))
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = color.best, weight = "bold"),
      locations = gt::cells_body(columns = BIC, rows = BIC == min(BIC))
    )
  
  # Highlight best WLB p-values (green bold) and failures (red)
  tbl <- tbl |>
    gt::tab_style(
      style = gt::cell_text(color = color.best, weight = "bold"),
      locations = gt::cells_body(columns = WLB1_fmt, rows = WLB1 == max(WLB1))
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = color.best, weight = "bold"),
      locations = gt::cells_body(columns = WLB2_fmt, rows = WLB2 == max(WLB2))
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = color.best, weight = "bold"),
      locations = gt::cells_body(columns = WLB3_fmt, rows = WLB3 == max(WLB3))
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = color.fail),
      locations = gt::cells_body(columns = WLB1_fmt, rows = WLB1 < 0.05)
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = color.fail),
      locations = gt::cells_body(columns = WLB2_fmt, rows = WLB2 < 0.05)
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = color.fail),
      locations = gt::cells_body(columns = WLB3_fmt, rows = WLB3 < 0.05)
    )
  
  # Nyblom: lower is better, red if unstable
  tbl <- tbl |>
    gt::tab_style(
      style = gt::cell_text(color = color.best, weight = "bold"),
      locations = gt::cells_body(columns = Nyblom_fmt, rows = Nyblom == min(Nyblom))
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = color.fail),
      locations = gt::cells_body(columns = Nyblom_fmt, rows = Nyblom >= Nyblom_crit)
    )
  
  # Sign Bias: higher p-value is better
  tbl <- tbl |>
    gt::tab_style(
      style = gt::cell_text(color = color.best, weight = "bold"),
      locations = gt::cells_body(columns = SignBias_fmt, rows = SignBias == max(SignBias))
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = color.fail),
      locations = gt::cells_body(columns = SignBias_fmt, rows = SignBias < 0.05)
    )
  
  # Coefficient significance
  tbl <- tbl |>
    gt::tab_style(
      style = gt::cell_text(color = color.best, weight = "bold"),
      locations = gt::cells_body(columns = CoefSig, rows = CoefSig_pass == TRUE)
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = color.fail),
      locations = gt::cells_body(columns = CoefSig, rows = CoefSig_pass == FALSE)
    )
  
  # Footnotes
  tbl <- tbl |>
    gt::tab_source_note("Green = best/pass. Red = fail. Weighted LB lags are model-dependent.") |>
    gt::tab_source_note(gt::html("Nyblom: H<sub>0</sub> = parameters stable. Green = stable (&lt; 5% critical value).")) |>
    gt::tab_source_note(gt::html("Sign Bias: H<sub>0</sub> = symmetric effects. If rejected, consider EGARCH/GJR-GARCH.")) |>
    gt::tab_options(
      table.font.size = gt::px(12),
      heading.title.font.size = gt::px(16),
      heading.subtitle.font.size = gt::px(13),
      column_labels.font.size = gt::px(11),
      column_labels.padding = gt::px(8),
      data_row.padding = gt::px(7),
      source_notes.font.size = gt::px(10),
      table.width = gt::pct(100)
    )
  
  tbl
}


# ==============================================================================
# CLI Table Display
# ==============================================================================

#' Display GARCH Comparison Table in Console
#'
#' Creates a formatted console table using the \pkg{cli} package with
#' color-coded highlighting.
#'
#' @param results A \code{garch.comparison.tse} object from
#'   \code{\link{compare.garch.tse}}.
#' @param show.signbias Logical. Include Sign Bias column? Default is
#'   \code{TRUE}.
#'
#' @return Invisibly returns the comparison tibble.
#'
#' @details
#' Colors in terminal output:
#' \itemize{
#'   \item \strong{Green}: Best value for that criterion, or passing diagnostic
#'   \item \strong{Red}: Failing diagnostic (p < 0.05 for tests, or
#'     Nyblom >= critical value)
#' }
#'
#' @seealso \code{\link{compare.garch.tse}}, \code{\link{table.garch.gt.tse}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' garch_gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.15, beta = 0.8)
#' y <- garch_gen(1000)
#' results <- compare.garch.tse(y)
#'
#' # Display in console
#' table.garch.cli.tse(results)
#' }
table.garch.cli.tse <- function(results, show.signbias = TRUE) {
  
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' is required. Install with: install.packages('cli')")
  }
  
  if (!inherits(results, "garch.comparison.tse")) {
    stop("results must be a 'garch.comparison.tse' object from compare.garch.tse()")
  }
  
  df <- results$comparison
  
  if (nrow(df) == 0) {
    cli::cli_alert_warning("No models to display (comparison table is empty)")
    return(invisible(df))
  }
  
  best <- .find_best_indices(df)
  
  # Column widths
  w_model <- 12
  w_ic <- 8
  w_wlb <- 7
  w_nyb <- 7
  w_sb <- 9
  w_coef <- 5
  sep <- " "
  
  # Formatting helpers (apply styling after sprintf)
  fmt_ic <- function(val, is_best) {
    s <- sprintf("%*.2f", w_ic, val)
    if (is_best) cli::col_green(cli::style_bold(s)) else s
  }
  
  fmt_pval <- function(val, is_best, width = w_wlb) {
    if (is.na(val)) {
      s <- sprintf("%*s", width, "NA")
      return(s)
    }
    if (val < 0.001) {
      s <- sprintf("%*s", width, "<.001")
      cli::col_red(s)
    } else if (val < 0.05) {
      s <- sprintf("%*.3f", width, val)
      cli::col_red(s)
    } else if (is_best) {
      s <- sprintf("%*.3f", width, val)
      cli::col_green(cli::style_bold(s))
    } else {
      sprintf("%*.3f", width, val)
    }
  }
  
  fmt_nyblom <- function(stat, crit, is_best) {
    if (is.na(stat)) {
      s <- sprintf("%*s", w_nyb, "NA")
      return(s)
    }
    s <- sprintf("%*.2f", w_nyb, stat)
    if (is_best) {
      cli::col_green(cli::style_bold(s))
    } else if (stat >= crit) {
      cli::col_red(s)
    } else {
      s
    }
  }
  
  fmt_coef <- function(n_sig, n_coef) {
    s <- sprintf("%*s", w_coef, paste0(n_sig, "/", n_coef))
    if (n_sig == n_coef) cli::col_green(s) else cli::col_red(s)
  }
  
  fmt_model <- function(name) {
    sprintf("%-*s", w_model, name)
  }
  
  # Header
  cli::cli_h2("GARCH Model Comparison")
  cli::cli_text("Distribution: {.val {results$distribution}}")
  cli::cli_text("")
  
  # Build header string
  header <- paste0(
    sprintf("%-*s", w_model, "Model"), sep,
    sprintf("%*s", w_ic, "AIC"), sep,
    sprintf("%*s", w_ic, "BIC"), sep,
    sprintf("%*s", w_ic, "AICc"), sep,
    sprintf("%*s", w_wlb, "WLB1"), sep,
    sprintf("%*s", w_wlb, "WLB2"), sep,
    sprintf("%*s", w_wlb, "WLB3"), sep,
    sprintf("%*s", w_nyb, "Nyblom")
  )
  
  if (show.signbias) {
    header <- paste0(header, sep, sprintf("%*s", w_sb, "SignBias"))
  }
  header <- paste0(header, sep, sprintf("%*s", w_coef, "Coef"))
  
  cat(cli::style_bold(header), "\n")
  cli::cli_rule()
  
  # Data rows
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    
    row <- paste0(
      fmt_model(r$Model), sep,
      fmt_ic(r$AIC, i == best$aic), sep,
      fmt_ic(r$BIC, i == best$bic), sep,
      fmt_ic(r$AICC, i == best$aicc), sep,
      fmt_pval(r$WLB1, i == best$wlb1), sep,
      fmt_pval(r$WLB2, i == best$wlb2), sep,
      fmt_pval(r$WLB3, i == best$wlb3), sep,
      fmt_nyblom(r$Nyblom, r$Nyblom_crit, i == best$nyblom)
    )
    
    if (show.signbias) {
      row <- paste0(row, sep, fmt_pval(r$SignBias, i == best$signbias, width = w_sb))
    }
    row <- paste0(row, sep, fmt_coef(r$n_sig, r$n_coef))
    
    cat(row, "\n")
  }
  
  cli::cli_rule()
  cli::cli_text(cli::col_silver("Green = best/pass | Red = fail | Nyblom: green = stable"))
  
  invisible(df)
}


# ==============================================================================
# Coefficient Table
# ==============================================================================

#' Display GARCH Coefficient Table
#'
#' Creates a formatted table of coefficient estimates, standard errors,
#' t-values, and p-values for a fitted GARCH model.
#'
#' @param fit A \code{uGARCHfit} object from \pkg{rugarch}, typically
#'   extracted from the \code{fits} element of a \code{garch.comparison.tse}
#'   object.
#' @param title Optional custom title. If \code{NULL}, auto-generated from
#'   model specification.
#' @param color.sig Color for significant p-values (< 0.05). Default is
#'   "seagreen".
#' @param color.nonsig Color for non-significant p-values. Default is
#'   "tomato".
#'
#' @return A \code{gt} table object.
#'
#' @details
#' Significance codes in the table:
#' \itemize{
#'   \item \code{***}: p < 0.001
#'   \item \code{**}: p < 0.01
#'   \item \code{*}: p < 0.05
#'   \item \code{.}: p < 0.1
#' }
#'
#' @seealso \code{\link{compare.garch.tse}}, \code{\link{table.garch.gt.tse}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' garch_gen <- make.gen.garch.tse(omega = 0.1, alpha = 0.15, beta = 0.8)
#' y <- garch_gen(1000)
#' results <- compare.garch.tse(y)
#'
#' # Get coefficient table for best model
#' best_fit <- results$fits[["GARCH(1,1)"]]
#' table.coef.gt.tse(best_fit)
#' }
table.coef.gt.tse <- function(fit,
                              title = NULL,
                              color.sig = "seagreen",
                              color.nonsig = "tomato") {
  
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required. Install with: install.packages('gt')")
  }
  
  if (!inherits(fit, "uGARCHfit")) {
    stop("fit must be a 'uGARCHfit' object from rugarch::ugarchfit()")
  }
  
  # Extract model orders for title
  arch_order <- fit@model$modelinc["alpha"]
  garch_order <- fit@model$modelinc["beta"]
  
  if (is.null(title)) {
    title <- paste0(.garch_label(arch_order, garch_order), " Coefficient Estimates")
  }
  
  # Extract coefficient matrix
  coef_mat <- fit@fit$matcoef
  
  df <- as.data.frame(coef_mat)
  colnames(df) <- c("Estimate", "Std_Error", "t_value", "p_val")
  df$Parameter <- rownames(coef_mat)
  
  df <- df[, c("Parameter", "Estimate", "Std_Error", "t_value", "p_val")]
  
  df$p_display <- sapply(df$p_val, function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) "<0.001" else sprintf("%.4f", p)
  })
  
  df$Sig <- sapply(df$p_val, function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else if (p < 0.1) "."
    else ""
  })
  
  tbl <- df |>
    gt::gt() |>
    gt::tab_header(
      title = title,
      subtitle = paste0("Distribution: ", fit@model$modeldesc$distribution)
    ) |>
    gt::cols_hide(columns = c(p_val)) |>
    gt::cols_label(
      Parameter = "Parameter",
      Estimate = "Estimate",
      Std_Error = "Std. Error",
      t_value = "t value",
      p_display = "p-value",
      Sig = ""
    ) |>
    gt::fmt_number(columns = c(Estimate, Std_Error, t_value), decimals = 4) |>
    gt::cols_align(align = "left", columns = Parameter) |>
    gt::cols_align(align = "right", columns = c(Estimate, Std_Error, t_value, p_display)) |>
    gt::cols_align(align = "center", columns = Sig) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(columns = Parameter)
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = color.sig, weight = "bold"),
      locations = gt::cells_body(columns = p_display, rows = p_val < 0.05)
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = color.nonsig),
      locations = gt::cells_body(columns = p_display, rows = p_val >= 0.05)
    ) |>
    gt::tab_source_note("Significance: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1") |>
    gt::tab_options(
      table.font.size = gt::px(12),
      heading.title.font.size = gt::px(15),
      heading.subtitle.font.size = gt::px(12),
      column_labels.padding = gt::px(8),
      data_row.padding = gt::px(6),
      source_notes.font.size = gt::px(10)
    )
  
  tbl
}