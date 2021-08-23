# note that R models the number of failures until `size` successes, whereas Brookmeyer's formulation gives the number of successes until gamma failures.
#' @importFrom stats pnbinom
lower = function(x, X, N, mu_shape, mu_scale, alpha, I = x, p = 1/((mu_scale*N*I) + 1))
{
 
  stats::pnbinom(
    q = X,
    size = mu_shape,
    prob = p,
    lower = TRUE) - alpha/2
  
}

#' @importFrom stats pnbinom
upper = function(x, X, N, mu_shape, mu_scale, alpha, I = x, p = 1/((mu_scale*N*I) + 1))
{
  
  stats::pnbinom(
    q = X-1,
    size = mu_shape,
    prob = p,
    lower = FALSE) - alpha/2
  
}

#' Estimate incidence from cross-sectional data
#'
#' @param X The number of MAA-positive participants in the cross-sectional sample
#' @param N0 The number of seronegative participatns in the cross-sectional sample
#' @param mu_estimate The estimated mean window period, mu
#' @param mu_CI_low The lower confidence limit for the mean window period, mu (can skip this parameter if providing mu_shape and mu_scale)
#' @param mu_CI_high The upper confidence limit for the mean window period, mu (can skip this parameter if providing mu_shape and mu_scale)
#' @param mu_CI_level The confidence level for the mean window period, mu (can skip this parameter if providing `mu_prior`, or `mu_shape` and `mu_scale`)
#' @param mu_prior The Gamma parameters for the assumed prior distribution on mu (can skip this argument if providing mu_shape and mu_scale directly)
#' @param mu_shape The shape (gamma) parameter for the prior distribution of the mean window period, mu.
#' @param mu_scale The scale (beta) parameter for the prior distribution of the mean window period, mu.
#' @param incidence_CI_level The desired confidence interval coverage probability 
#' @param alpha One minus the desired coverage probability
#' @param max_incidence_CI_range The range of incidence rate values to search over (if not wide enough, an error will be produced).
#' @param tol The desired tolerance for matching the specified confidence level
#' @return a numeric vector containing the estimated incidence rate, lower confidence limit, and upper confidence limit, in units of new infections per 100 person-years.
#' @export
#'

#' @importFrom stats uniroot
est.inc = function(
  X,
  N0,
  mu_estimate,
  mu_CI_low,
  mu_CI_high,
  mu_CI_level = .95,
  mu_prior = fit_gamma_distribution(
    mean = mu_estimate,
    lower = mu_CI_low,
    upper = mu_CI_high,
    confidence_level = mu_CI_level),
  mu_shape = mu_prior[1],
  mu_scale = mu_prior[2],
  incidence_CI_level = .95,
  alpha = 1 - incidence_CI_level,
  tol = 1e-8,
  max_incidence_CI_range = c(0,100)
  # note: we don't have to worry about mis-specifying `max_incidence_CI_range`, 
  # because the `uniroot` function returns an error if the curve does not cross 0.
)
{
  
  # DEM note: the distribution of X is negative binomial; could use packaged methods for exact CIs on 'p', and convert to I, if such methods exist.
  upper_CI = stats::uniroot(
    lower, 
    interval = max_incidence_CI_range, 
    X = X, 
    N = N0, 
    mu_shape = mu_shape, 
    mu_scale = mu_scale,
    alpha = alpha,
    tol = tol)$root
  
  lower_CI = stats::uniroot(
    upper, 
    interval = max_incidence_CI_range, 
    tol = tol,
    X = X,
    N = N0, 
    mu_shape = mu_shape, 
    mu_scale = mu_scale,
    alpha = alpha)$root
  
  to_return = c(
    estimate = X / (N0 * mu_estimate),
    CI_low = lower_CI, 
    CI_high = upper_CI)
  
  
  return(to_return)
  
}
