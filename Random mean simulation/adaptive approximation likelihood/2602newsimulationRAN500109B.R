##################
# Load required packages
library(bin2norm)  # Custom package
library(optparse)  # Command-line argument parser
library(lme4)      # Mixed models
library(rstan)     # Bayesian inference
library(statmod)

generate_multiple_thresholds_data <- function(I, sigma, r, sample_size_approach){

  if(sample_size_approach=="A"){
    n_i_vec <- round(runif(I, min=10, max=150))
  } else if(sample_size_approach=="B"){
    n_i_vec <- round(runif(I, min=10, max=15000))
  } else {
    n_i_vec <- round(runif(I, min=10000, max=15000))
  }

  if(r <= 0 || r >= 1){
    if(r==1) stop("r=1 => sigma=0. Not in standard set.")
    if(r==0) stop("Use fixed-mean function for r=0.")
    stop("r must be in (0,1).")
  }
  tau_sq <- (r/(1-r)) * sigma^2
  tau <- sqrt(tau_sq)

  mu0 <- 1
  data_list <- list(
    n_i = n_i_vec,
    c_ij = vector("list", I),
    p_ij_obs = vector("list", I)
  )

  for(i in seq_len(I)){
    mu_i <- rnorm(1, mean=mu0, sd=tau)
    m_i <- sample(5:10, size=1)

    q1 <- qnorm(0.2, mean=mu_i, sd=sigma)
    q2 <- qnorm(0.8, mean=mu_i, sd=sigma)

    c_ij_vec <- numeric(m_i)
    p_ij_obs_vec <- numeric(m_i)

    for(j in seq_len(m_i)){
      c_ij <- runif(1, min=q1, max=q2)
      p_ij <- 1 - pnorm((c_ij - mu_i)/sigma)
      k_ij <- rbinom(1, size=n_i_vec[i], prob=p_ij)
      c_ij_vec[j]     <- c_ij
      p_ij_obs_vec[j] <- k_ij / n_i_vec[i]
    }

    data_list$c_ij[[i]]     <- c_ij_vec
    data_list$p_ij_obs[[i]] <- p_ij_obs_vec
  }

  return(data_list)
}


logLik_multinom_integrate_adaptive_oneStudy <- function(
    ni, thresholds, counts, par,
    K = 10,
    mode_bound = 15,
    abs.tol = 1e-10,
    rel.tol = 1e-10,
    clamp = 1e-12,
    eps_hess = 1e-3,
    min_width = 0.5,
    max_width = 12
) {
  mu0   <- par[1]
  sigma <- par[2]
  tau   <- par[3]
  big <- -1e300  # for log-lik failure
  
  if (!is.finite(mu0) || !is.finite(sigma) || !is.finite(tau) || sigma <= 0 || tau < 0) return(big)
  
  thresholds <- as.numeric(thresholds)
  counts <- as.integer(counts)
  m <- length(thresholds)
  
  if (m <= 0 || length(counts) != (m + 1L) || ni <= 0) return(big)
  
  # Ensure thresholds are sorted and counts aligned
  ord <- order(thresholds)
  thresholds <- thresholds[ord]
  
  # ---- conditional loglik as a function of x (GH x-space) ----
  # mu = mu0 + sqrt(2) * tau * x
  logLi_x <- function(x) {
    mu_i <- mu0 + sqrt(2) * tau * x
    
    # right-tail p_k(mu) = Phi((mu - c_k)/sigma)
    p_vec <- pnorm((mu_i - thresholds) / sigma)
    
    # multinomial cell probs pi_1..pi_{m+1}
    pi_vec <- numeric(m + 1L)
    pi_vec[1] <- 1 - p_vec[1]
    if (m > 1) {
      pi_vec[2:m] <- p_vec[1:(m - 1)] - p_vec[2:m]
    }
    pi_vec[m + 1L] <- p_vec[m]
    
    # numeric safety
    pi_vec <- pmin(pmax(pi_vec, clamp), 1 - clamp)
    pi_vec <- pi_vec / sum(pi_vec)
    
    # multinomial loglik (constant term included; doesn't matter for MLE, but harmless)
    stats::dmultinom(counts, prob = pi_vec, log = TRUE)
  }
  
  # If tau == 0: random effect collapses, no integral needed
  if (tau == 0) {
    ll0 <- logLi_x(0)  # because mu = mu0 when tau=0, corresponds to x=0
    if (!is.finite(ll0)) return(big)
    return(ll0)
  }
  
  h <- function(x) logLi_x(x) - x^2
  
  # ---- 1) find mode x* of h on [-mode_bound, mode_bound] ----
  opt <- optimize(h, interval = c(-mode_bound, mode_bound), maximum = TRUE)
  x_star <- opt$maximum
  h_star <- opt$objective
  if (!is.finite(h_star)) return(big)
  
  # ---- 2) curvature for local width ----
  e <- eps_hess
  hp <- h(x_star + e)
  hm <- h(x_star - e)
  if (!is.finite(hp) || !is.finite(hm)) {
    e <- 5 * eps_hess
    hp <- h(x_star + e)
    hm <- h(x_star - e)
    if (!is.finite(hp) || !is.finite(hm)) return(big)
  }
  h2 <- (hp - 2 * h_star + hm) / (e^2)
  
  if (!is.finite(h2) || (-h2) <= 0) {
    half_width <- 4
  } else {
    s <- 1 / sqrt(-h2)
    half_width <- K * s
    half_width <- max(half_width, min_width)
    half_width <- min(half_width, max_width)
  }
  
  L <- x_star - half_width
  U <- x_star + half_width
  
  # ---- 3) centered integrand exp(h - h*) in [0,1] to avoid under/overflow ----
  integrand <- function(x_vec) {
    hx <- vapply(x_vec, h, numeric(1))
    hx[!is.finite(hx)] <- -Inf
    exp(hx - h_star)
  }
  
  val <- try(stats::integrate(integrand, lower = L, upper = U,
                              abs.tol = abs.tol, rel.tol = rel.tol),
             silent = TRUE)
  if (inherits(val, "try-error") || !is.finite(val$value) || val$value <= 0) return(big)
  
  # log L_i = -0.5 log(pi) + h* + log integral
  logLi_marg <- (-0.5 * log(pi)) + h_star + log(val$value)
  if (!is.finite(logLi_marg)) return(big)
  
  logLi_marg
}



estimate_multiThresh_MLE_adaptive <- function(
    data_list,
    K = 10,
    mode_bound = 15,
    abs.tol = 1e-10,
    rel.tol = 1e-10,
    clamp = 1e-12,
    eps_hess = 1e-3,
    min_width = 0.5,
    max_width = 12
) {
  n_i <- data_list$n_i
  c_ij <- data_list$c_ij
  p_ij_obs <- data_list$p_ij_obs
  I <- length(n_i)
  
  # ---- Recalculate counts K_i from observed right-tail p_ij_obs (same as your GH code) ----
  k_ij <- vector("list", I)
  for (i in seq_len(I)) {
    p_vec <- as.numeric(p_ij_obs[[i]])
    thresholds <- as.numeric(c_ij[[i]])
    m <- length(p_vec)
    
    ord <- order(thresholds)
    thresholds <- thresholds[ord]
    p_vec <- p_vec[ord]
    c_ij[[i]] <- thresholds
    
    # Convert cumulative right-tail to interval probabilities (m+1 cells)
    pi_vec <- numeric(m + 1)
    pi_vec[1] <- 1 - p_vec[1]
    if (m > 1) {
      pi_vec[2:m] <- p_vec[1:(m - 1)] - p_vec[2:m]
    }
    pi_vec[m + 1] <- p_vec[m]
    
    pi_vec <- pmax(pmin(pi_vec, 1), 0)
    pi_vec <- pi_vec / sum(pi_vec)
    
    # counts (same rounding scheme you used)
    k_vec <- round(n_i[i] * pi_vec)
    diff <- n_i[i] - sum(k_vec)
    if (diff != 0) {
      max_idx <- which.max(pi_vec)
      k_vec[max_idx] <- k_vec[max_idx] + diff
    }
    k_ij[[i]] <- as.integer(k_vec)
  }
  
  # ---- Initial parameter estimates ----
  initial <- estimate_initial_values_from_data(data_list)
  init_par <- c(mu0 = initial$mu0, sigma = initial$sigma, tau = initial$tau)
  
  # ---- Adaptive-integrated negative log-likelihood ----
  neg_ll <- function(par) {
    mu0 <- par[1]; sigma <- par[2]; tau <- par[3]
    if (!is.finite(mu0) || sigma <= 0 || tau < 0) return(1e10)
    
    ll_sum <- 0
    for (i in seq_len(I)) {
      thresholds <- c_ij[[i]]
      counts <- k_ij[[i]]
      
      ll_i <- logLik_multinom_integrate_adaptive_oneStudy(
        ni = n_i[i],
        thresholds = thresholds,
        counts = counts,
        par = par,
        K = K,
        mode_bound = mode_bound,
        abs.tol = abs.tol,
        rel.tol = rel.tol,
        clamp = clamp,
        eps_hess = eps_hess,
        min_width = min_width,
        max_width = max_width
      )
      
      if (!is.finite(ll_i) || ll_i < -1e200) return(1e10)
      ll_sum <- ll_sum + ll_i
    }
    -ll_sum
  }
  
  # ---- Optimization ----
  res <- optim(
    init_par, neg_ll,
    method = "L-BFGS-B",
    lower = c(-Inf, 1e-4, 0),
    control = list(maxit = 500)
  )
  
  list(
    mu0 = res$par[1],
    sigma = res$par[2],
    tau = res$par[3],
    method = "MLE_adaptive_integration_multinom",
    init_par = init_par,
    optim = res
  )
}



set.seed(2)

simresult <- matrix(0, ncol = 3, nrow = 1000)
colnames(simresult) <- c(
  "mle_int_mu0","mle_int_sigma","mle_int_tau"
)

for (i in 1:1000) {
  dat <- generate_multiple_thresholds_data(50, 0.1, 0.9, "B")
  MLE  <- estimate_multiThresh_MLE_adaptive(dat)

  simresult[i,] <- c(
    MLE$mu,MLE$sigma,MLE$tau
  )
    cat(i)
}

saveRDS(simresult, file = "2602adaptiveRAN500109B.rds")
