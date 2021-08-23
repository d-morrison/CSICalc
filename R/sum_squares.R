
#' Calculate the error in matching a specified confidence interval with a Gamma distribution prior.
#'
#' @param shape The shape ("gamma") parameter of a Gamma distribution
#' @param lower_CI The lower confidence limit that we are trying to match
#' @param upper_CI The upper confidence limit that we are trying to match
#' @param mean The mean of the Gamma distribution (= shape * scale)
#' @param confidence_level The confidence level of the confidence interval (1 - alpha)
#' @param alpha  = 1 - `confidence_level`
#' @return The sum of relative squared errors between the confidence limits and the quantiles of the fitted distribution.

#' @importFrom stats qgamma
sum_squares = function(
  shape, 
  lower_CI, 
  upper_CI, 
  mean,
  confidence_level = .95,
  
  alpha = 
    (1-confidence_level)/2)
{
  
  lower_quantile = 
    stats::qgamma(
      shape = shape, 
      scale = mean/shape, 
      p = alpha, 
      lower = TRUE)
  
  upper_quantile = 
    stats::qgamma(
      shape = shape, 
      scale = mean/shape, 
      p = alpha, 
      lower = FALSE)
  
  SS = 
    (lower_quantile - lower_CI)^2/lower_CI^2 + 
    (upper_quantile - upper_CI)^2/upper_CI^2
  
  return(SS)
  
}