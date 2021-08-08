fit_gamma_distribution = function(
  mean,
  lower,
  upper,
  alpha = .95,
  gamma_shape_range = c(0,100)
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
  
  return(c(shape, scale))
  
}  

  