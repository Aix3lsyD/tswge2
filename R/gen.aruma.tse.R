library(ggplot2)

gen.aruma.tse <- function(n, phi=0, theta=0, d=0, s=0, lambda=0, 
                          innov_gen=NULL, plot=TRUE, sn=0)
{
  # n: realization length
  # phi: vector of AR parameters (ATSA sign convention)
  # theta: vector of MA parameters (ATSA sign convention)
  # d: order of difference operator
  # s: order of seasonal operator
  # lambda: vector of parameters in nonstationary operator lambda(B)
  # innov_gen: function with signature function(n) returning numeric vector of length n
  #            with mean zero. If NULL, uses rnorm(n, mean=0, sd=1).
  # plot: if TRUE, plot the generated realization
  # sn: seed for reproducibility (0 = no seed)
  
  # Set seed if provided (before any random generation)
  if (sn > 0) {set.seed(sn)}
  
  # Default innovation generator: standard normal
  if (is.null(innov_gen)) {
    innov_gen <- function(n) rnorm(n, mean=0, sd=1)
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
  
  # Calculate total length needed
  spin <- 100
  ngen <- n + dlams + spin
  
  # Generate innovations with buffer for arima.sim burn-in
  arima_burnin <- 100
  innovations <- innov_gen(ngen + arima_burnin)
  
  # Simulate ARMA process with our custom innovations
  if ((p > 0) & (q > 0)) {
    tsdata <- arima.sim(n=ngen, model=list(order=c(p, d, q), ar=ar, ma=ma), innov=innovations)
  }
  if ((p == 0) & (q > 0)) {
    tsdata <- arima.sim(n=ngen, model=list(order=c(p, d, q), ma=ma), innov=innovations)
  }
  if ((p > 0) & (q == 0)) {
    tsdata <- arima.sim(n=ngen, model=list(order=c(p, d, q), ar=ar), innov=innovations)
  }
  if ((p == 0) & (q == 0)) {
    tsdata <- arima.sim(n=ngen, model=list(order=c(0, d, 0)), innov=innovations)
  }
  
  # Compute the inverse of the nonstationary operator (similar to diffinv in R)
  y <- as.numeric(tsdata)
  
  if ((dlam > 0) & (s > 0)) {
    temp <- mult.wge(fac1=lambda, fac2=seas)
    lambdas <- temp$model.coef
  }
  if ((dlam > 0) & (s == 0)) {
    lambdas <- lambda
  }
  if ((dlam == 0) & (s > 0)) {
    lambdas <- seas
  }
  if ((dlam == 0) & (s == 0)) {
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
        xfull[i] <- xfull[i] + lambdas[j] * xfull[i-j]
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
  gg <- ggplot(df, aes(x = time, y = y)) +
    geom_line(color = "steelblue", linewidth = 0.5) +
    labs(x = "Time", y = "Value", title = "ARUMA Realization") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  # Show plot if requested
  if (plot == TRUE) {
    print(gg)
  }
  
  # Build result object
  result <- list(
    y = y_final,
    innovations = innovations[1:n],  # return just the n used for output
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


# Print method for aruma class
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
  print(head(x$y))
  
  invisible(x)
}


# Plot method for aruma class
plot.aruma <- function(x, ...) {
  print(x$plot)
}