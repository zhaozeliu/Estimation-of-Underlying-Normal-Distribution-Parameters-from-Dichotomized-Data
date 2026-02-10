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
  
  # random-mean model
  # sigma^2 given, then tau^2 so that tau^2/(sigma^2+tau^2)= r
  # => tau^2 = (r/(1-r)) * sigma^2
  if(r <= 0 || r >= 1){
    # if r=0.9 => tau^2= 9*sigma^2, if r=0.1 => tau^2= 0.111... * sigma^2
    # We'll handle boundary conditions if needed, but user sets r in {0.1,0.5,0.9}
    if(r==1) stop("r=1 => sigma=0. Not in standard set.")
    if(r==0) stop("Use fixed-mean function for r=0.")
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
    # sample number of thresholds m_i
    m_i <- sample(5:10, size=1)
    
    # quantiles for threshold sampling
    q1 <- qnorm(0.2, mean=mu_i, sd=sigma)
    q2 <- qnorm(0.8, mean=mu_i, sd=sigma)
    
    c_ij_vec <- numeric(m_i)
    p_ij_obs_vec <- numeric(m_i)
    
    for(j in seq_len(m_i)){
      c_ij <- runif(1, min=q1, max=q2)
      p_ij <- 1 - pnorm((c_ij - mu_i)/sigma)
      k_ij <- rbinom(1, size=n_i_vec[i], prob=p_ij)
      c_ij_vec[j]      <- c_ij
      p_ij_obs_vec[j]  <- k_ij / n_i_vec[i]
    }
    
    data_list$c_ij[[i]]     <- c_ij_vec
    data_list$p_ij_obs[[i]] <- p_ij_obs_vec
  }
  
  return(data_list)
}



set.seed(1)
simresult=matrix(data = c(rep(0,9000)),ncol = 9,nrow = 1000)
colnames(simresult)=c("mle_int_mu0","glmm_mu0","bayes_mu0","mle_int_sigma","glmm_sigma","bayes_sigma","mle_int_tau","glmm_tau","bayes_tau")
for (i in 1:1000) {
  MLEglmmtry=generate_multiple_thresholds_data(50,1.0,0.9,"A")
  MLE=estimate_multiThresh_MLE(MLEglmmtry)
  GLMM=estimate_multiThresh_GLMM(MLEglmmtry)
  MCMC=estimate_multiThresh_MCMC(MLEglmmtry)
  simresult[i,]=c(MLE$mu,GLMM$mu,MCMC$mu0_est,MLE$sigma,GLMM$sigma,MCMC$sigma_est,MLE$tau,GLMM$tau,MCMC$tau_est)
}
saveRDS(simresult, file = "newsimulationRAN501009A.rds")





