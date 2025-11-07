# ============================================================
# bin2norm: all-in-one R script for single-threshold vs multiple-threshold
# ============================================================

#' bin2norm: A user-friendly interface to estimate normal distribution parameters
#'    from dichotomized data
#'
#' @description
#' This function handles two data-collection settings for estimating normal parameters
#' from threshold-based (dichotomized) data:
#'
#' \itemize{
#'   \item \strong{Single-threshold per study:}
#'     Each of \eqn{I} studies reports one threshold \eqn{c_i}, a sample size \eqn{n_i},
#'     and the observed proportion \eqn{p_i^{obs}} of samples above that threshold.
#'     We assume one normal distribution \eqn{\mathcal{N}(\mu,\sigma^2)} across all studies.
#'     Methods include \code{"MLE"} and \code{"probit"}.
#'
#'   \item \strong{Multiple-thresholds per study:}
#'     Each study \eqn{i} reports \eqn{K_i} thresholds \eqn{\{c_{ij}\}}, each with an
#'     observed proportion \eqn{p_{ij}^{obs}}. We assume the study-specific mean
#'     \eqn{\mu_i \sim \mathcal{N}(\mu_0,\tau^2)} and within-study variance \eqn{\sigma^2}.
#'     Because each study has multiple cutpoints, one can estimate \eqn{\mu_0, \sigma, \tau}.
#'     Methods include \code{"MLE_integration"}, \code{"GLMM"}, or \code{"Bayesian"} (MCMC).
#' }
#'
#' @param scenario character string, either \code{"single_threshold"} or \code{"multiple_thresholds"}.
#' @param method character string indicating which estimation method to use.
#'   \itemize{
#'     \item For \code{scenario = "single_threshold"}, valid \code{method} are \code{"MLE"} or \code{"probit"}.
#'     \item For \code{scenario = "multiple_thresholds"}, valid \code{method} are
#'           \code{"MLE_integration"}, \code{"GLMM"}, or \code{"Bayesian"}.
#'   }
#' @param n_i,c_i,p_i_obs used \strong{only if \code{scenario="single_threshold"}}. Numeric vectors of
#'   the same length.  \eqn{n_i} is study sample size, \eqn{c_i} is threshold, \eqn{p_i_obs} is observed
#'   proportion above threshold.
#' @param data_list used \strong{only if \code{scenario="multiple_thresholds"}}, a list with:
#'   \itemize{
#'     \item \code{n_i}: numeric vector (length I) of sample sizes
#'     \item \code{c_ij}: list of length I, where \code{c_ij[[i]]} is a numeric vector of thresholds in study i
#'     \item \code{p_ij_obs}: list of length I, where \code{p_ij_obs[[i]]} is a numeric vector of observed proportions above each threshold
#'   }
#' @param ... additional arguments passed to lower-level functions (e.g. \code{use_wols_init},
#'   \code{gh_points}, \code{iter}, \code{chains}, etc.).
#'
#' @return A list of estimated parameters, depending on the data-collection setting
#'   (\code{scenario}) and the chosen method. Typically includes:
#'   \itemize{
#'     \item \code{mu} or \code{mu0}
#'     \item \code{sigma}
#'     \item \code{tau} (only for multiple-threshold methods)
#'   }
#'
#' @examples
#' \dontrun{
#' # Single-threshold example
#' n_i <- c(100, 120, 80)
#' c_i <- c(1.2, 1.0, 1.5)
#' p_i_obs <- c(0.30, 0.25, 0.40)
#' bin2norm(scenario="single_threshold", method="MLE", n_i=n_i, c_i=c_i, p_i_obs=p_i_obs)
#'
#' # Multiple-thresholds example
#' data_list <- list(
#'   n_i = c(100, 120),
#'   c_ij = list(c(1.0,1.2), c(0.8,1.5,2.0)),
#'   p_ij_obs = list(c(0.20,0.30), c(0.15,0.40,0.55))
#' )
#'
#' # MLE with numeric integration
#' bin2norm(scenario="multiple_thresholds", method="MLE_integration",
#'          data_list=data_list, gh_points=5)
#'
#' # GLMM approximation
#' # library(lme4)
#' bin2norm(scenario="multiple_thresholds", method="GLMM",
#'          data_list=data_list, use_lme4=TRUE)
#'
#' # Bayesian MCMC approach
#' # library(rstan)
#' bin2norm(scenario="multiple_thresholds", method="Bayesian",
#'          data_list=data_list, iter=1000, chains=2)
#' }
#' @export
bin2norm <- function(scenario = c("single_threshold","multiple_thresholds"),
                     method  = NULL,
                     n_i     = NULL,
                     c_i     = NULL,
                     p_i_obs = NULL,
                     data_list = NULL,
                     ...)
{
  scenario <- match.arg(scenario)

  if(scenario == "single_threshold"){
    # single threshold per study
    if(is.null(n_i) || is.null(c_i) || is.null(p_i_obs)){
      stop("For single_threshold, must provide n_i, c_i, p_i_obs.")
    }
    if( (length(n_i)!=length(c_i)) || (length(c_i)!=length(p_i_obs)) ){
      stop("single_threshold: n_i, c_i, p_i_obs must have the same length.")
    }
    if(is.null(method)) method <- "MLE"
    if(! method %in% c("MLE","probit")){
      stop("For single_threshold, recognized methods are 'MLE' or 'probit'.")
    }
    if(method=="MLE"){
      return( estimate_singleThresh_MLE(n_i, c_i, p_i_obs, ...) )
    } else {
      return( estimate_singleThresh_probit(n_i, c_i, p_i_obs) )
    }

  } else if(scenario == "multiple_thresholds"){
    # multiple thresholds per study
    if(is.null(data_list)){
      stop("For multiple_thresholds, must provide data_list with n_i, c_ij, p_ij_obs.")
    }
    if(is.null(data_list$n_i) || is.null(data_list$c_ij) || is.null(data_list$p_ij_obs)){
      stop("data_list must contain 'n_i', 'c_ij', and 'p_ij_obs'.")
    }
    I <- length(data_list$n_i)
    if( length(data_list$c_ij) != I || length(data_list$p_ij_obs) != I ){
      stop("In data_list, c_ij and p_ij_obs must each have length = length(n_i).")
    }
    if(is.null(method)) method <- "MLE_integration"
    if(! method %in% c("MLE_integration","GLMM","Bayesian")){
      stop("For multiple_thresholds, recognized methods: 'MLE_integration','GLMM','Bayesian'.")
    }
    if(method=="MLE_integration"){
      return( estimate_multiThresh_MLE(data_list, ...) )
    } else if(method=="GLMM"){
      return( estimate_multiThresh_GLMM(data_list, ...) )
    } else {
      return( estimate_multiThresh_MCMC(data_list, ...) )
    }
  } else {
    stop("scenario must be 'single_threshold' or 'multiple_thresholds'.")
  }
}


# ============================================================
# Single-threshold: MLE and GLM probit
# ============================================================

#' MLE (Single Threshold per Study)
#'
#' @description
#' Treats the count of "above threshold" in study \eqn{i} as binomial with probability
#' \eqn{1 - \Phi((c_i - \mu)/\sigma)}. This uses numerical optimization (\code{optim})
#' to maximize the binomial likelihood. Optionally uses Weighted OLS estimates as starting
#' values to improve convergence.
#'
#' @param n_i numeric vector of sample sizes
#' @param c_i numeric vector of thresholds
#' @param p_i_obs numeric vector of observed proportions above threshold
#' @param use_wols_init logical; if \code{TRUE}, uses Weighted OLS estimates
#'   (\code{\link{estimate_singleThresh_WOLS}}) as initial values in \code{optim}.
#' @return A list with \code{mu}, \code{sigma}, \code{method="MLE"}.
#' @importFrom stats pnorm dbinom optim sd
#' @export
estimate_singleThresh_MLE <- function(n_i, c_i, p_i_obs, use_wols_init = TRUE){
  # 1) filter once on observed proportions and sanity-check inputs
  valid <- is.finite(n_i) & is.finite(c_i) & is.finite(p_i_obs) &
    (p_i_obs > 0) & (p_i_obs < 1)
  n_i     <- as.numeric(n_i[valid])
  c_i     <- as.numeric(c_i[valid])
  p_i_obs <- as.numeric(p_i_obs[valid])
  
  # reconstruct counts consistently and clamp to [0, n]
  k_i <- as.integer(round(n_i * p_i_obs))
  k_i <- pmin(pmax(k_i, 0L), n_i)
  
  big <- 1e100  # large finite penalty instead of Inf
  
  neg_loglik <- function(par){
    mu <- par[1]; sigma <- par[2]
    # always return finite values
    if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) return(big)
    if (any(!is.finite(n_i)) || any(!is.finite(k_i)) || any(!is.finite(c_i))) return(big)
    if (any(k_i < 0 | k_i > n_i)) return(big)
    
    # model p = 1 - Phi((c - mu)/sigma) = Phi((mu - c)/sigma)
    prob_i <- pnorm((mu - c_i) / sigma)
    # 2) clamp to avoid log(0)
    prob_i <- pmin(pmax(prob_i, 1e-12), 1 - 1e-12)
    
    val <- -sum(dbinom(k_i, size = n_i, prob = prob_i, log = TRUE))
    if (!is.finite(val)) big else val
  }
  
  # 3) initializer
  if (use_wols_init){
    wfit <- estimate_singleThresh_WOLS(n_i, c_i, p_i_obs)
    init_par <- c(unname(wfit$mu), abs(unname(wfit$sigma)))
    if (!all(is.finite(init_par)) || init_par[2] <= 0)
      init_par <- c(mean(c_i), max(sd(c_i), 0.1))
  } else {
    init_par <- c(mean(c_i), max(sd(c_i), 0.1))
  }
  
  # 4) optimize with proper lower bound on sigma
  res <- optim(
    par     = init_par,
    fn      = neg_loglik,
    method  = "L-BFGS-B",
    lower   = c(-Inf, 1e-8),
    control = list(maxit = 5000, pgtol = 1e-12)
  )
  
  list(mu = res$par[1],
       sigma = res$par[2],
       
       method = "MLE")
}


#' Weighted OLS (Initial value in Single Threshold per Study MLE)
#'
#' @description
#' Implements the formula \eqn{ c_i = \mu + \sigma * \Phi^{-1}(1 - p_i^{obs}) }
#' in a weighted least-squares sense, with weights = \eqn{n_i}.
#'
#' @param n_i numeric vector
#' @param c_i numeric vector
#' @param p_i_obs numeric vector
#' @return A list with \code{mu}, \code{sigma}.
#' @importFrom stats lm coef qnorm
#' @export
estimate_singleThresh_WOLS <- function(n_i, c_i, p_i_obs){
  # Filter out p_i_obs values that are 0 or 1
  valid_idx <- (p_i_obs > 0) & (p_i_obs < 1)
  n_i <- n_i[valid_idx]
  c_i <- c_i[valid_idx]
  p_i_obs <- p_i_obs[valid_idx]
  
  # Compute inverse normal
  x <- qnorm(1 - p_i_obs)
  fit <- lm(c_i ~ x, weights=n_i)
  mu_hat <- coef(fit)[1]
  sigma_hat <- abs(coef(fit)[2])
  list(mu=mu_hat, sigma=sigma_hat, method="WOLS")
}


#' GLM probit (Single Threshold per Study)
#'
#' @description
#' For each group \eqn{i}, we assume the data follows:
#' \deqn{ \Pr(Y_i = 1) = \Phi\left( \frac{\mu - c_i}{\sigma} \right) }
#' where \eqn{c_i} is a known threshold, and \eqn{\Phi} is the standard normal CDF (the probit link).
#' The function reconstructs individual binary outcomes based on observed probabilities,
#' and estimates the parameters using generalized linear modeling with a probit link.
#'
#' @param n_i numeric vector
#' @param c_i numeric vector
#' @param p_i_obs numeric vector
#'
#' @return A list with \code{mu}, \code{sigma}, \code{method="probit"}.
#' @importFrom stats glm coef vcov binomial
#' @export
estimate_singleThresh_probit <- function(n_i, c_i, p_i_obs) {
  # Keep only valid observations
  valid_idx <- (p_i_obs > 0) & (p_i_obs < 1)
  n_i <- n_i[valid_idx]
  c_i <- c_i[valid_idx]
  p_i_obs <- p_i_obs[valid_idx]
  
  # Reconstruct individual binary outcomes: y = 1 if X > c_i, 0 otherwise
  y <- unlist(mapply(function(n, p) {
    k <- round(n * p)
    c(rep(1, k), rep(0, n - k))
  }, n_i, p_i_obs))
  
  c_rep <- unlist(mapply(function(n, cval) rep(cval, n), n_i, c_i))
  
  
  fit <- glm(y ~ I(-c_rep), family = binomial(link = "probit"))
  
  # Extract coefficients
  alpha <- coef(fit)[1]
  beta  <- coef(fit)[2]
  
  
  mu_hat    <- alpha / beta
  sigma_hat <- abs(1 / beta)
  
  
  vcov_mat <- vcov(fit)
  se_alpha <- sqrt(vcov_mat[1, 1])
  se_beta  <- sqrt(vcov_mat[2, 2])
  cov_ab   <- vcov_mat[1, 2]
  
  se_mu    <- sqrt((1 / beta^2)^2 * se_alpha^2 + (alpha / beta^2)^2 * se_beta^2 + 2 * (alpha / beta^3) * cov_ab)
  se_sigma <- abs(1 / beta^2) * se_beta
  
  list(mu = mu_hat,
       sigma = sigma_hat,
       method = "probit")
}


# ============================================================
# Multiple-thresholds: MLE_integration, GLMM, Bayesian MCMC
# ============================================================

#' MLE with Numeric Integration (Multiple Thresholds per Study)
#'
#' @description
#' Each study \eqn{i} has thresholds \eqn{\{c_{ij}\}}, each with an observed proportion
#' \eqn{p_{ij}^{obs}}. We assume \eqn{\mu_i \sim \mathcal{N}(\mu_0,\tau^2)} and
#' \eqn{X_{ij} \sim \mathcal{N}(\mu_i,\sigma^2)}. The log-likelihood integrates out
#' \eqn{\mu_i} via Gauss-Hermite quadrature.
#'
#' @param data_list A list with:
#'   \itemize{
#'     \item \code{n_i}: numeric vector (length I)
#'     \item \code{c_ij}: list of length I
#'     \item \code{p_ij_obs}: list of length I
#'   }
#' @param gh_points integer; number of Gauss-Hermite points (default 12).
#'
#' @return A list with \code{mu0}, \code{sigma}, \code{tau}, \code{method="MLE_integration"}.
#' @export
estimate_multiThresh_MLE <- function(data_list, gh_points = 20) {
  n_i <- data_list$n_i
  c_ij <- data_list$c_ij
  p_ij_obs <- data_list$p_ij_obs
  I <- length(n_i)
  
  # Recalculate k_ij based on ordered thresholds and p_ij_obs
  k_ij <- vector("list", I)
  for (i in seq_len(I)) {
    p_vec <- p_ij_obs[[i]]
    thresholds <- c_ij[[i]]
    m <- length(p_vec)
    
    # Sort thresholds and reorder p_vec accordingly
    ord <- order(thresholds)
    thresholds <- thresholds[ord]
    p_vec <- p_vec[ord]
    c_ij[[i]] <- thresholds  # Update sorted thresholds in data_list
    
    # Convert cumulative right-tail to interval probabilities
    pi_vec <- numeric(m + 1)
    pi_vec[1] <- 1 - p_vec[1]
    if (m > 1) {
      pi_vec[2:m] <- p_vec[1:(m - 1)] - p_vec[2:m]
    }
    pi_vec[m + 1] <- p_vec[m]
    
    # Normalize and handle numerical edge cases
    pi_vec <- pmax(pmin(pi_vec, 1), 0)
    pi_vec <- pi_vec / sum(pi_vec)
    
    # Compute expected counts
    k_vec <- round(n_i[i] * pi_vec)
    
    # Adjust rounding error
    diff <- n_i[i] - sum(k_vec)
    if (diff != 0) {
      max_idx <- which.max(pi_vec)
      k_vec[max_idx] <- k_vec[max_idx] + diff
    }
    
    k_ij[[i]] <- k_vec
  }
  
  # Gauss-Hermite quadrature
  gh <- gauss.quad(n = gh_points, kind = "hermite")
  xnodes <- gh$nodes
  wnodes <- gh$weights
  
  # Initial parameter estimates
  initial <- estimate_initial_values_from_data(data_list)
  #initial=estimate_multiThresh_GLMM(data_list)
  init_par <- c(mu0 = initial$mu0, sigma = initial$sigma, tau = initial$tau)
  #init_par <- c(mu0 = 1, sigma = 0.1, tau = 0.033)
  # Negative log-likelihood function
  neg_ll <- function(par) {
    mu0 <- par[1]
    sigma <- par[2]
    tau <- par[3]
    if (!is.finite(mu0) || sigma <= 0 || tau < 0) return(1e10)
    
    ll_sum <- 0
    for (i in seq_len(I)) {
      thresholds <- c_ij[[i]]
      counts <- k_ij[[i]]
      m_i <- length(thresholds)
      
      partial_logL <- numeric(gh_points)
      for (g in seq_len(gh_points)) {
        mu_i <- mu0 + sqrt(2) * tau * xnodes[g]
        z <- (thresholds - mu_i) / sigma
        p_vec <- 1 - pnorm(z)
        
        pi_vec <- numeric(m_i + 1)
        pi_vec[1] <- 1 - p_vec[1]
        if (m_i > 1) {
          pi_vec[2:m_i] <- p_vec[1:(m_i - 1)] - p_vec[2:m_i]
        }
        pi_vec[m_i + 1] <- p_vec[m_i]
        pi_vec <- pmin(pmax(pi_vec, 1e-10), 1 - 1e-10)
        pi_vec <- pi_vec / sum(pi_vec)
        
        partial_logL[g] <- log(wnodes[g]) + dmultinom(counts, prob = pi_vec, log = TRUE)
      }
      
      max_logL <- max(partial_logL)
      log_sum <- max_logL + log(sum(exp(partial_logL - max_logL))) - 0.5 * log(pi)
      ll_sum <- ll_sum + log_sum
    }
    
    return(-ll_sum)
  }
  
  # Optimization
  res <- optim(
    init_par, neg_ll,
    method = "L-BFGS-B",
    lower = c(-Inf, 1e-4, 0),
    control = list(maxit = 500)
  )
  # res <- optim(
  #   init_par, neg_ll,
  #   method = "BFGS",        # better on unconstrained problems
  #   control = list(maxit = 1000, reltol = 1e-10, parscale = c(1, 0.1, 0.1))
  # )
  
  return(list(
    mu0 = res$par[1],
    sigma = res$par[2],
    tau = res$par[3],
    method = "MLE_integration",
    init_par = init_par
  ))
}
# ============================================================
# Initial value for MLE_integration, weighted least¨Csquares regression
# ============================================================
#' Get initial values from data
#' @param data_list your inputs
#' @return a named list of initial values
#' @export
estimate_initial_values_from_data<- function(data_list) {
  n_i <- data_list$n_i
  c_ij <- data_list$c_ij
  p_ij_obs <- data_list$p_ij_obs
  I <- length(n_i)
  
  mu_i_est <- numeric(I)
  q_all <- numeric(0)
  c_all <- numeric(0)
  
  for (i in seq_len(I)) {
    p_ij <- p_ij_obs[[i]]
    c_vals <- c_ij[[i]]
    valid <- (p_ij > 0) & (p_ij < 1)
    
    if (sum(valid) < 2) next  # need at least 2 points to estimate
    
    q <- qnorm(1 - p_ij[valid])
    c_use <- c_vals[valid]
    
    fit <- lm(c_use ~ q)
    mu_i_est[i] <- coef(fit)[1]
    q_all <- c(q_all, q)
    c_all <- c(c_all, c_use)
  }
  
  # Refit global sigma from pooled data
  fit_global <- lm(c_all ~ q_all)
  mu0_init <- mean(mu_i_est, na.rm = TRUE)
  sigma_init <- abs(coef(fit_global)[2])
  tau_init <- sd(mu_i_est, na.rm = TRUE)
  
  list(mu0 = mu0_init, sigma = sigma_init, tau = tau_init)
}




#' GLMM (Multiple Thresholds per Study, Probit Link, Random Intercepts)
#'
#' @description
#' Creates a single data frame stacking all thresholds from all studies, then
#' calls \code{lme4::glmer(..., family=binomial(link='probit'))} to fit a random-intercept
#' model:
#' \deqn{k_{ij} \sim \mathrm{Binomial}\bigl(n_i, \Phi(\alpha_i + \beta\, c_{ij})\bigr),}
#' with \eqn{\alpha_i \sim \mathcal{N}(0, \sigma_\alpha^2)}.
#'
#' Interpreting results: \eqn{\sigma = 1/|\,\beta\,|},
#' \eqn{\tau^2 = \sigma^2 \times \sigma_\alpha^2},
#' \eqn{\mu_0 = (\mathrm{Intercept}) \times \sigma} (if not forced to 0).
#'
#' @param data_list same structure: \code{n_i}, \code{c_ij}, \code{p_ij_obs}
#' @param use_lme4 logical; if \code{TRUE}, calls \code{lme4::glmer} with a probit link.
#'
#' @return A list with \code{mu0}, \code{sigma}, \code{tau}, \code{method="GLMM_probit"}.
#' @importFrom lme4 glmer fixef VarCorr
#' @importFrom stats binomial
#' @export
estimate_multiThresh_GLMM <- function(data_list, use_lme4=TRUE){
  n_i <- data_list$n_i
  c_ij <- data_list$c_ij
  p_ij_obs <- data_list$p_ij_obs
  I <- length(n_i)

  # build long data frame
  studyID <- integer(0)
  cVal    <- numeric(0)
  kVal    <- integer(0)
  nVal    <- integer(0)

  for(i in seq_len(I)){
    Ki <- length(c_ij[[i]])
    for(j in seq_len(Ki)){
      studyID <- c(studyID, i)
      cVal    <- c(cVal,  c_ij[[i]][j])
      kVal    <- c(kVal,  round(n_i[i]* p_ij_obs[[i]][j]))
      nVal    <- c(nVal,  n_i[i])
    }
  }
  K <- length(studyID)
  df <- data.frame(study=factor(studyID), cVal=cVal, kVal=kVal, nVal=nVal)

  if(!use_lme4){
    stop("Set use_lme4=TRUE to actually run a probit GLMM with lme4.")
  }
  fit <- lme4::glmer(cbind(kVal, nVal-kVal) ~ cVal + (1|study),
                     data=df, family=binomial(link="probit"))

  fixB <- lme4::fixef(fit)
  intercept <- unname(fixB[1])   # alpha_0
  slope     <- unname(fixB[2])   # beta
  varRandom <- as.numeric(lme4::VarCorr(fit)$study[1,1])  # sigma_alpha^2

  sigma <- 1/abs(slope)
  tau   <- sigma * sqrt(varRandom)
  mu0   <- intercept * sigma

  list(mu0=mu0, sigma=sigma, tau=tau, method="GLMM_probit")
}


#' Bayesian MCMC (Multiple Thresholds per Study) using rstan
#'
#' @description
#' Builds an inline Stan model for multiple thresholds per study. The user must have
#' the \code{rstan} package installed. We place random effects \eqn{\mu_i = \mu_0 + \tau * mu\_raw[i]}
#' and use a binomial likelihood for each threshold. By default, uses simple weakly
#' informative priors.
#'
#' @param data_list same structure as above: \code{n_i}, \code{c_ij}, \code{p_ij_obs}
#' @param iter number of total iterations for each chain (default 2000)
#' @param chains number of MCMC chains (default 2)
#'
#' @return a list containing \code{stan_fit} (the full Stan fit object), plus
#'   \code{mu0_est}, \code{sigma_est}, \code{tau_est} as posterior means, and
#'   \code{method="Bayesian_MCMC"}.
#'
#' @importFrom rstan stan extract
#' @export
estimate_multiThresh_MCMC <- function(data_list, iter=2000, chains=2){
  n_i <- data_list$n_i
  c_ij <- data_list$c_ij
  p_ij_obs <- data_list$p_ij_obs
  I <- length(n_i)

  studyID <- integer(0)
  cVal    <- numeric(0)
  kVal    <- integer(0)
  nVal    <- integer(0)

  for(i in seq_len(I)){
    Ki <- length(c_ij[[i]])
    for(j in seq_len(Ki)){
      studyID <- c(studyID, i)
      cVal    <- c(cVal, c_ij[[i]][j])
      kVal    <- c(kVal, round(n_i[i]* p_ij_obs[[i]][j]))
      nVal    <- c(nVal, n_i[i])
    }
  }
  K <- length(studyID)

  stan_code <- "
data {
  int<lower=1> I;
  int<lower=1> K;
  int<lower=1,upper=I> studyID[K];
  int<lower=0> kVal[K];
  int<lower=0> nVal[K];
  real cVal[K];
}
parameters {
  real mu0;
  real<lower=0> sigma;
  real<lower=0> tau;
  vector[I] mu_raw;
}
transformed parameters {
  vector[I] mu_i;
  for(i in 1:I){
    mu_i[i] = mu0 + tau * mu_raw[i];
  }
}
model {
  // Weakly informative priors
  mu0 ~ normal(0,10);
  sigma ~ cauchy(0,5);
  tau ~ cauchy(0,5);
  mu_raw ~ normal(0,1);

  // Binomial likelihood
  for(k in 1:K){
    real p_ijk = 1 - Phi( (cVal[k] - mu_i[ studyID[k] ]) / sigma );
    kVal[k] ~ binomial(nVal[k], p_ijk);
  }
}
"

stan_data <- list(I=I, K=K, studyID=studyID, kVal=kVal, nVal=nVal, cVal=cVal)
fit <- rstan::stan(model_code=stan_code, data=stan_data, iter=iter, chains=chains)
post <- rstan::extract(fit)

list(
  stan_fit   = fit,
  mu0_est    = mean(post$mu0),
  sigma_est  = mean(post$sigma),
  tau_est    = mean(post$tau),
  method     = "Bayesian_MCMC"
)
}


# ---------------------------------------------------------
# Gauss-Hermite Quadrature (Minimal Demo)
# ---------------------------------------------------------

#' Minimal Gauss-Hermite Quadrature
#'
#' @description
#' Returns \code{(nodes, weights)} for approximating \eqn{\int f(x) e^{-x^2} dx},
#' ignoring any normalizing constant. This is a simple demonstration; for serious
#' applications, more robust libraries or expansions might be used.
#'
#' @param n integer number of quadrature points
#' @return list with \code{nodes} and \code{weights}
#' @export
gaussHermite <- function(n){
  # This is a toy example for n=2 or n=3 only.
  if(n == 2){
    nodes   <- c(-1/sqrt(2), 1/sqrt(2))
    weights <- c(1/sqrt(2),  1/sqrt(2))
  } else if(n == 3){
    nodes   <- c(-sqrt(3/2),  0, sqrt(3/2))
    # approximate weights for demonstration
    weights <- c(0.5, 1, 0.5)
  } else {
    warning("gaussHermite: only minimal demonstration for n=2,3 provided.")
    nodes   <- seq(-1, 1, length.out=n)
    weights <- rep(1, n)
  }
  list(nodes=nodes, weights=weights)
}




