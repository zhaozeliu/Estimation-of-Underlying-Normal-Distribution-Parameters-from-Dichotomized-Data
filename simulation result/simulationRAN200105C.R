########################################################################
# R Simulation Script for Compute Canada (SLURM Array Jobs)
########################################################################

##################
# Load required packages
library(bin2norm)  # Custom package
library(optparse)  # Command-line argument parser
library(lme4)      # Mixed models
library(rstan)     # Bayesian inference
library(statmod)
suppressPackageStartupMessages({
  library(optparse)  # or you can just parse commandArgs() manually
})

option_list <- list(
  make_option("--scenario",   type="character", default="fixed_mean_single",
              help="Scenario: 'fixed_mean_single' or 'random_mean_multi'"),
  make_option("--Ival",       type="integer",   default=20,
              help="Number of studies: e.g. 20 or 50"),
  make_option("--sigma",      type="double",    default=0.1,
              help="Within-study standard deviation (0.1 or 1)"),
  make_option("--rval",       type="double",    default=0.1,
              help="Ratio r = tau^2 / (sigma^2 + tau^2). 0 => fixed-mean scenario."),
  make_option("--approach",   type="character", default="A",
              help="Sample-size approach: A,B,C"),
  make_option("--rep_start",  type="integer",   default=1,
              help="Starting replicate index for this job"),
  make_option("--rep_end",    type="integer",   default=100,
              help="Ending replicate index for this job"),
  make_option("--run_glmm",   type="logical",   default=TRUE,
              help="Whether to run GLMM method in multiple-threshold scenario"),
  make_option("--run_bayes",  type="logical",   default=TRUE,
              help="Whether to run Bayesian MCMC in multiple-threshold scenario"),
  make_option("--outdir",     type="character", default=".",
              help="Output directory for results.")
)
opt <- parse_args(OptionParser(option_list=option_list))



scenario   <- opt$scenario
I_         <- opt$Ival
sig_       <- opt$sigma
r_         <- opt$rval
sample_app <- opt$approach
rep_start  <- opt$rep_start
rep_end    <- opt$rep_end
run_glmm   <- opt$run_glmm
run_bayes  <- opt$run_bayes
outdir     <- opt$outdir

################################

cat(sprintf("Starting job with scenario=%s, I=%d, sigma=%g, r=%g, approach=%s, reps=[%d..%d]\n",
            scenario, I_, sig_, r_, sample_app, rep_start, rep_end))

# ------------------------------------------------------------------
# Paste or source your data generation and simulation functions here
# (For brevity, we only include an example snippet).
# ------------------------------------------------------------------

generate_single_threshold_data <- function(I, sigma, sample_size_approach){
  
  # pick the range for n_i according to approach A,B,C
  if(sample_size_approach=="A"){
    # small
    n_i_vec <- round(runif(I, min=20, max=50))
  } else if(sample_size_approach=="B"){
    # mixed
    n_i_vec <- round(runif(I, min=20, max=200))
  } else {
    # large
    n_i_vec <- round(runif(I, min=150, max=200))
  }
  
  mu <- 1   # fixed
  df <- data.frame(
    i = seq_len(I),
    n_i = n_i_vec,
    c_i = NA_real_,
    p_i_obs = NA_real_
  )
  
  # For each study, sample c_i
  # We define q1,q2 as 20%-80% quantiles of Normal(mu, sigma^2)
  # i.e. q1= qnorm(0.2, mean=mu, sd=sigma), q2= qnorm(0.8, mean=mu, sd=sigma)
  q1 <- qnorm(0.2, mean=mu, sd=sigma)
  q2 <- qnorm(0.8, mean=mu, sd=sigma)
  
  for(j in seq_len(I)){
    c_ij <- runif(1, min=q1, max=q2)  # threshold for study j
    p_ij <- 1 - pnorm( (c_ij - mu)/sigma ) 
    k_ij <- rbinom(1, size=df$n_i[j], prob=p_ij)
    df$c_i[j]    = c_ij
    df$p_i_obs[j] = k_ij / df$n_i[j]
  }
  
  return(df)
}

generate_multiple_thresholds_data <- function(I, sigma, r, sample_size_approach){
  
  if(sample_size_approach=="A"){
    n_i_vec <- round(runif(I, min=20, max=50))
  } else if(sample_size_approach=="B"){
    n_i_vec <- round(runif(I, min=20, max=200))
  } else {
    n_i_vec <- round(runif(I, min=150, max=200))
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

run_single_threshold_rep <- function(df){
  # df has columns i, n_i, c_i, p_i_obs
  # We call bin2norm with scenario="single_threshold"
  
  # MLE
  fit_mle <- bin2norm(
    scenario="single_threshold",
    method="MLE",
    n_i = df$n_i,
    c_i = df$c_i,
    p_i_obs = df$p_i_obs
  )
  
  # probit
  fit_probit <- bin2norm(
    scenario="single_threshold",
    method="probit",
    n_i = df$n_i,
    c_i = df$c_i,
    p_i_obs = df$p_i_obs
  )
  
  # Return a named vector or list
  out <- c(
    mu_mle   = fit_mle$mu,
    sig_mle  = fit_mle$sigma,
    mu_probit  = fit_probit$mu,
    sig_probit = fit_probit$sigma
  )
  return(out)
}

run_multiple_thresholds_rep <- function(data_list, run_glmm=FALSE, run_bayes=FALSE){
  
  res <- list()
  
  # MLE_integration
  fit_mle_int <- bin2norm(
    scenario="multiple_thresholds",
    method="MLE_integration",
    data_list=data_list,
    gh_points=3  # for speed
  )
  res$mle_int_mu0   <- fit_mle_int$mu0
  res$mle_int_sigma <- fit_mle_int$sigma
  res$mle_int_tau   <- fit_mle_int$tau
  
  # GLMM
  if(run_glmm){
    # requires lme4
    fit_glmm <- bin2norm(
      scenario="multiple_thresholds",
      method="GLMM",
      data_list=data_list,
      use_lme4=TRUE
    )
    res$glmm_mu0   <- fit_glmm$mu0
    res$glmm_sigma <- fit_glmm$sigma
    res$glmm_tau   <- fit_glmm$tau
  }
  
  # Bayesian
  if(run_bayes){
    # requires rstan
    # note: might be slow if iter is large
    fit_bayes <- bin2norm(
      scenario="multiple_thresholds",
      method="Bayesian",
      data_list=data_list,
      iter=1000,
      chains=2
    )
    res$bayes_mu0   <- fit_bayes$mu0_est
    res$bayes_sigma <- fit_bayes$sigma_est
    res$bayes_tau   <- fit_bayes$tau_est
  }
  
  return(unlist(res))
}




# ------------------------------------------------------------------
# We'll run replicates from rep_start to rep_end in this job
# Store results in a data.frame, then save to .rds
# ------------------------------------------------------------------


set.seed(1)

scenario   <- 'random_mean_multi'
I_         <- 20
sig_       <- 0.1
r_         <- 0.5
sample_app <- "C"
rep_start  <- 1
rep_end    <- 1000
run_glmm   <- TRUE
run_bayes  <- TRUE
outdir     <- "."

cat(sprintf("Starting replicate loop from %d to %d...\n", rep_start, rep_end))
all_res <- list()
idx <- 1
rep_loop_time  <- system.time({
for(rp in seq.int(rep_start, rep_end)){
  if(scenario == "fixed_mean_single"){
    # r should be 0 in principle, we can ignore r_ 
    df_data <- generate_single_threshold_data(I_, sig_, sample_app)
    est_res <- run_single_threshold_rep(df_data)
    # True mu=1, True sigma=sig_
    all_res[[idx]] <- data.frame(
      scenario    = scenario,
      replicate   = rp,
      I           = I_,
      sigma       = sig_,
      approach    = sample_app,
      mu_mle      = est_res[1],
      sig_mle     = est_res[2],
      mu_probit     = est_res[3],
      sig_probit    = est_res[4]
    )
    
  } else if(scenario == "random_mean_multi"){
    # we interpret r_ > 0
    data_list <- generate_multiple_thresholds_data(I_, sig_, r_, sample_app)
    est_res <- run_multiple_thresholds_rep(data_list, run_glmm=run_glmm, run_bayes=run_bayes)
    
    # parse out 
    # e.g. est_res might have "mle_int_mu0","mle_int_sigma","mle_int_tau", etc.
    # be sure to check the names if unlist(...) changes them
    eN <- names(est_res)
    rec <- data.frame(
      scenario   = scenario,
      replicate  = rp,
      I          = I_,
      sigma      = sig_,
      r          = r_,
      approach   = sample_app,
      mle_int_mu0   = est_res[1],
      mle_int_sigma = est_res[2],
      mle_int_tau   = est_res[3]
    )
    # fill in each method if present
    if("mle_int_mu0" %in% eN){
      rec$mle_int_mu0   = est_res[1]
      rec$mle_int_sigma = est_res[2]
      rec$mle_int_tau   = est_res[3]
    }
    if("glmm_mu0" %in% eN){
      rec$glmm_mu0      = est_res["glmm_mu0"]
      rec$glmm_sigma    = est_res["glmm_sigma"]
      rec$glmm_tau      = est_res["glmm_tau"]
    }
    if("bayes_mu0" %in% eN){
      rec$bayes_mu0     = est_res["bayes_mu0"]
      rec$bayes_sigma   = est_res["bayes_sigma"]
      rec$bayes_tau     = est_res["bayes_tau"]
    }
    all_res[[idx]] <- rec
  }
  idx <- idx + 1
}
})
cat(sprintf("Finished replicate loop from %d to %d.\n", rep_start, rep_end))
cat(sprintf("Time taken for replicate loop [%d..%d]: %.2f seconds (elapsed)\n",
            rep_start, rep_end, rep_loop_time["elapsed"]))

final_df <- do.call(rbind, all_res)

# create an output filename 
# e.g.: "results_fixed_mean_single_I20_sigma0.1_approachA_part1.rds"
tag_scenario <- ifelse(scenario=="fixed_mean_single","FixedMean","RandomMean")
tag_glmm     <- ifelse(run_glmm, "GLMM","noGLMM")
tag_bayes    <- ifelse(run_bayes,"Bayes","noBayes")

ofile <- sprintf("results_%s_I%d_sigma%.2f_r%.2f_%s_%s_%s_reps%d-%d.rds",
                 tag_scenario, I_, sig_, r_, sample_app, tag_glmm, tag_bayes, rep_start, rep_end)
ofile_path <- file.path(outdir, ofile)

saveRDS(final_df, file=ofile)
cat(sprintf("Time taken for replicate loop [%d..%d]: %.2f seconds\n",
            rep_start, rep_end, rep_loop_time["elapsed"]))
flush.console()







