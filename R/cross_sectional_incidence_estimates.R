#' Infer incidence rates from cross-sectional data
#'
#' @param mu_table 
#' @param mu_estimates 
#' @param mu_CIs_lower 
#' @param mu_CIs_upper 
#' @param X 
#' @param N 
#' @param results 
#' @param results2 
#' @param true_incidence 
#' @param digits 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
estimate_incidence = function(
  mu_table,
  mu_estimates = mu_table$"mu",
  mu_CIs_lower = mu_table$"2.5%",
  mu_CIs_upper = mu_table$"97.5%",
  X,
  N,
  
  results = estimate_incidence_from_cohort(
    mu_estimate = mu_estimates,
    mu_CI_lower = mu_CIs_lower,
    mu_CI_upper = mu_CIs_upper,
    X = X,
    N = N,
    ...),
  
  results2 = round(results, 1),
  
  true_incidence = 1.9,
  
  digits = 1,
  
  ...)

{
  
  {# functions
    
    sum_squares = function(
      alpha, 
      lower_CI, 
      upper_CI, 
      estimate,
      confidence_level = .95,
      
      percentile = 
        (1-confidence_level)/2)
    {
      
      lower_quantile = 
        qgamma(
          shape = alpha, 
          scale = estimate/alpha, 
          p = percentile, 
          lower = TRUE)
      
      upper_quantile = 
        qgamma(
          shape = alpha, 
          scale = estimate/alpha, 
          p = percentile, 
          lower = FALSE)
      
      SS = 
        (lower_quantile - lower_CI)^2/lower_CI^2 + 
        (upper_quantile - upper_CI)^2/upper_CI^2
      
      return(SS)
      
    }
    
    lower = function(x, X, N, shape, scale, alpha, I = x/100/365, p = 1/((scale*N*I) + 1))
    {
      
      pnbinom(
        q = X,
        size = shape,
        prob = p,
        lower = TRUE) - alpha/2
      
    }
    
    upper = function(x, X, N, shape, scale, alpha, I = x/100/365, p = 1/((scale*N*I) + 1))
    {
      
      pnbinom(
        q = X-1,
        size = shape,
        prob = p,
        lower = FALSE) - alpha/2
      
    }
    
    estimate_incidence_from_cohort0 = function(
      mu_estimate,
      mu_CI_lower,
      mu_CI_upper,
      X,
      N,
      incidence_CI_alpha = 0.05,
      gamma_shape_range = c(0,100),
      max_incidence_CI_range = c(0,100)
      # DEM note: it appears we don't have to worry about mis-specifying `max_incidence_CI_range`, 
      # because the `uniroot` function returns an error if the curve does not cross 0.
    )
    {
      
      optimum = optimize(
        sum_squares, 
        interval = gamma_shape_range, 
        estimate = mu_estimate,
        lower_CI = mu_CI_lower,
        upper_CI = mu_CI_upper
      )
      
      shape = optimum$minimum
      scale = mu_estimate/shape
      
      rel_dist_from_edge = (gamma_shape_range[2] - shape)/diff(gamma_shape_range)
      
      if(rel_dist_from_edge < 0.05) stop("Expand range of possible gamma shape parameter values.")
      # DEM: not sure whether this is a reasonable and sufficient check; deserves more thought. 
      # DEM, 2021/07/22: changing this check to only check agains the top edge of the range; the shape parameter can't go any lower than 0.
      
      if(FALSE) # examine the results
      {      
        sum_squares(0, 
                    estimate = mu_estimate,
                    lower_CI = mu_CI_lower,
                    upper_CI = mu_CI_upper
        )
        
        plot(
          function(x) sum_squares(x, 
                                  estimate = mu_estimate,
                                  lower_CI = mu_CI_lower,
                                  upper_CI = mu_CI_upper
          ),
          xlim = gamma_shape_range)
        
        plot(function(x) dgamma(x, shape = shape, scale = scale), xlim = c(mu_CI_lower, mu_CI_upper))
        abline(v = mu_estimate, col = 'red')
        plot(function(x) dgamma(x, shape = 95, scale = mu_estimate/95), xlim = c(mu_CI_lower, mu_CI_upper)/365, add = TRUE, col = "red")
        
      }
      
      
      # DEM note: the distribution of X is negative binomial; could use packaged methods for exact CIs on 'p', and convert to I, if such methods exist.
      upper_CI = uniroot(
        lower, 
        interval = max_incidence_CI_range, 
        X = X, 
        N = N, 
        shape = shape, 
        scale = scale,
        alpha = incidence_CI_alpha)$root
      
      lower_CI = uniroot(
        upper, 
        interval = max_incidence_CI_range, 
        X = X,
        N = N, 
        shape = shape, 
        scale = scale,
        alpha = incidence_CI_alpha)$root
      
      est_I = X/N/mu_estimate * 365 * 100 # yearly incidence, in percent
      
      to_return = data.frame(
        estimated_incidence = est_I,
        CI95_low = lower_CI,
        CI95_high = upper_CI,
        minimum = optimum$minimum,
        objective = optimum$objective)
      
      
      return(
        to_return)
      
    }
    
    estimate_incidence_from_cohort1 = Vectorize(
      estimate_incidence_from_cohort0,
      vectorize.args = c(
        "mu_estimate",
        "mu_CI_lower",
        "mu_CI_upper",
        "X",
        "N"),
      SIMPLIFY = FALSE)
    
    estimate_incidence_from_cohort = function(...)
    {
      
      temp =  estimate_incidence_from_cohort1(...)
      return(do.call(rbind, temp))
      
    }
    
  }
  
  results2$Error = paste(formatC((results$estimated_incidence/true_incidence - 1) * 100, format = "f", digits = digits), "%", sep = "") 
  
  results2$Error_low = paste(formatC((results$CI95_low/true_incidence - 1) * 100, format = "f", digits = digits), "%", sep = "")
  results2$Error_high = paste(formatC((results$CI95_high/true_incidence - 1) * 100, format = "f", digits = digits), "%", sep = "")
  
  
  rownames(results2) = NULL
  
  results2 = cbind(mu_estimates, mu_CIs_lower, mu_CIs_upper, X, N, results2, results)
  
  return(results2)
  
}

