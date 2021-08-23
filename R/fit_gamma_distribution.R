#' Find gamma distribution parameters to match a point estimate and confidence interval
#'
#' @param mean The point estimate
#' @param lower The lower confidence limit
#' @param upper The upper confidence limit
#' @param confidence_level The coverage probability of the confidence interval
#' @param gamma_shape_range The range of gamma (shape) parameter values to search over (if insufficient, an error will be produced).
#'
#' @return a vector containing the estimated shape and scale parameters
#' @export
#'
#' @importFrom stats optimize
fit_gamma_distribution = function(
  mean,
  lower,
  upper,
  confidence_level = .95,
  gamma_shape_range = c(0,1000)
)
{
  
  optimum = stats::optimize(
    sum_squares, 
    interval = gamma_shape_range, 
    mean = mean,
    lower_CI = lower,
    upper_CI = upper,
    confidence_level = confidence_level)
  
  shape = optimum$minimum
  scale = mean / shape
  
  rel_dist_from_edge = (gamma_shape_range[2] - shape) / diff(gamma_shape_range)
  if(rel_dist_from_edge < 0.05) stop("Expand range of possible gamma shape parameter values.")
  
  return(c(shape = shape, scale = scale))
  
}  

  