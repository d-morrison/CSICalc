
#' Title
#'
#' @param shape 
#' @param lower_CI 
#' @param upper_CI 
#' @param mean 
#' @param confidence_level 
#' @param percentile 
#'
#' @return
#' @export
#'

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
    qgamma(
      shape = shape, 
      scale = mean/shape, 
      p = alpha, 
      lower = TRUE)
  
  upper_quantile = 
    qgamma(
      shape = shape, 
      scale = mean/shape, 
      p = alpha, 
      lower = FALSE)
  
  SS = 
    (lower_quantile - lower_CI)^2/lower_CI^2 + 
    (upper_quantile - upper_CI)^2/upper_CI^2
  
  return(SS)
  
}