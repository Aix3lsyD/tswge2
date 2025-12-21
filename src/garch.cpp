#include <Rcpp.h>
using namespace Rcpp;

//' Generate a GARCH(p,q) realization
 //' 
 //' @param n Length of the output series
 //' @param alpha0 Constant term (omega) in variance equation
 //' @param alpha Vector of ARCH coefficients (length q)
 //' @param beta Vector of GARCH coefficients (length p)
 //' @param eps Vector of standard normal innovations (length n + burn)
 //' @param burn Length of burn-in period to discard
 //' @return List containing the series (a), conditional variances (sigma2), and innovations used
 //' 
 // [[Rcpp::export]]
 List garch_sim_cpp(int n, 
                    double alpha0, 
                    NumericVector alpha, 
                    NumericVector beta,
                    NumericVector eps,
                    int burn = 500) {
   
   int q = alpha.size();  // ARCH order
   int p = beta.size();   // GARCH order
   int ntot = n + burn;
   
   // Check that we have enough innovations
   if (eps.size() < ntot) {
     stop("Not enough innovations provided. Need at least n + burn.");
   }
   
   // Initialize vectors
   NumericVector sigma2(ntot);  // conditional variance
   NumericVector a(ntot);       // the GARCH series (a_t = sigma_t * eps_t)
   
   // Calculate unconditional variance for initialization
   // E[sigma^2] = alpha0 / (1 - sum(alpha) - sum(beta))
   double sum_alpha = sum(alpha);
   double sum_beta = sum(beta);
   double persistence = sum_alpha + sum_beta;
   
   double uncond_var;
   if (persistence < 1.0) {
     uncond_var = alpha0 / (1.0 - persistence);
   } else {
     // For IGARCH or explosive case, just use alpha0 as starting point
     uncond_var = alpha0 / 0.1;  // reasonable default
   }
   
   // Initialize first max(p, q) values with unconditional variance
   int max_pq = std::max(p, q);
   for (int t = 0; t < max_pq; t++) {
     sigma2[t] = uncond_var;
     a[t] = eps[t] * sqrt(sigma2[t]);
   }
   
   // Main GARCH recursion
   for (int t = max_pq; t < ntot; t++) {
     // Start with constant
     sigma2[t] = alpha0;
     
     // Add ARCH terms: alpha_i * a_{t-i}^2
     for (int i = 0; i < q; i++) {
       sigma2[t] += alpha[i] * a[t - i - 1] * a[t - i - 1];
     }
     
     // Add GARCH terms: beta_j * sigma2_{t-j}
     for (int j = 0; j < p; j++) {
       sigma2[t] += beta[j] * sigma2[t - j - 1];
     }
     
     // Generate observation: a_t = eps_t * sigma_t
     a[t] = eps[t] * sqrt(sigma2[t]);
   }
   
   // Extract post-burn-in values
   NumericVector a_out(n);
   NumericVector sigma2_out(n);
   NumericVector eps_out(n);
   
   for (int i = 0; i < n; i++) {
     a_out[i] = a[burn + i];
     sigma2_out[i] = sigma2[burn + i];
     eps_out[i] = eps[burn + i];
   }
   
   return List::create(
     Named("a") = a_out,
     Named("sigma2") = sigma2_out,
     Named("sigma") = sqrt(sigma2_out),
     Named("eps") = eps_out
   );
 }