
#############################################################
## Bayesian decision for optimal RMST
#############################################################

BayesDecision_RMST <- function(result,
                               prior_expert_list=list(c(5,6),
                                                      c(10,12),
                                                      c(8.5,10)),
                               seed=1234,eta=0.95,model=NULL){
  df_result <- result
  options(warn=-1)
  set.seed(seed)
  # Apply row-wise simulation
  n_sim <- 1000
  tau_grid <- df_result$time
  sims_0 <- mapply(function(mu, sd) {
    rnorm(n = n_sim, mean = mu, sd = sd)
  }, df_result$surv_0, df_result$std_err_0)
  sims_1 <- mapply(function(mu, sd) {
    rnorm(n = n_sim, mean = mu, sd = sd)
  }, df_result$surv_1, df_result$std_err_1)
  delta <- sims_1 - sims_0
  delta_sd <- apply(delta,2,sd, na.rm=TRUE)
  Omega <- delta/delta_sd
  Omega[is.infinite(Omega)] <- 0
  Omega[is.na(Omega)] <- 0
  if(is.null(model)){model = "model"} else{ model = unique(df_result$model)}
  gridOmega <- data.frame(model=model,
                          time=df_result$time,
                          omega_median=apply(Omega,2,median,na.rm=TRUE),
                          omega_q025=apply(Omega,2,quantile,0.025,na.rm=TRUE),
                          omega_q975=apply(Omega,2,quantile,0.975,na.rm=TRUE),
                          omega_mean=apply(Omega,2,mean,na.rm=TRUE),
                          omega_sd=apply(Omega,2,sd,na.rm=TRUE))
  gridOmega <- tibble(gridOmega)
  Omega[is.na(Omega)] <- 0
  # unconstrained - non-informative
  interval_list <- list(range(tau_grid))
  interval_idx <- lapply(interval_list, function(interval) {
    which(tau_grid >= interval[1] & tau_grid <= interval[2])
  })
  tau_max_mat <- sapply(seq_along(interval_idx), function(e) {
    idx <- interval_idx[[e]]
    tau_grid[idx][max.col(Omega[, idx, drop = FALSE])]
  })
  sigma_e <- sapply(interval_idx, function(idx) {
    mean(apply(Omega[, idx, drop = FALSE], 1, sd))
  })
  tau_max_mat <- as.matrix(tau_max_mat[, colSums(!is.na(tau_max_mat)) > 0])
  sigma_e <- sigma_e[complete.cases(sigma_e)]
  w <- exp(-1 / sigma_e)
  w <- as.matrix(w / sum(w,na.rm=TRUE))
  tau_star <- tau_max_mat%*%w
  tau_star <- as.vector(tau_star)
  prob_vec <- numeric(length(tau_star)) # sample size decision rule
  for (j in seq_along(tau_star)) {
    prob_vec[j] <- mean(delta[j, 1:which.min(abs(tau_grid - mean(tau_star, na.rm=TRUE)))] > 0, na.rm = TRUE)
  }
  decision <- mean(as.numeric(prob_vec > eta)) # decision rule
  results <- list()
  results$Unconstrained_Restriction_Time <- list(summary=tibble(x = tau_star) %>%
                                                   summarise(
                                                     median = median(x, na.rm=TRUE),
                                                     q025  = quantile(x, 0.025, na.rm=TRUE),
                                                     q975  = quantile(x, 0.975, na.rm=TRUE),
                                                     mean = mean(x, na.rm=TRUE),
                                                     sd = sd(x, na.rm=TRUE)
                                                   ),
                                                 tau_star = tibble(tau_star, prob_vec),
                                                 decision = decision
  )
  # constrained restriction time - multiple expert knowledge
  interval_list <- c(prior_expert_list, list(range(tau_grid)))
  interval_idx <- lapply(interval_list, function(interval) {
    which(tau_grid >= interval[1] & tau_grid <= interval[2])
  })
  tau_max_mat <- sapply(seq_along(interval_idx), function(e) {
    idx <- interval_idx[[e]]
    tau_grid[idx][max.col(Omega[, idx, drop = FALSE])]
  })
  sigma_e <- sapply(interval_idx, function(idx) {
    mean(apply(Omega[, idx, drop = FALSE], 1, sd))
  })
  tau_max_mat <- as.matrix(tau_max_mat[, colSums(!is.na(tau_max_mat)) > 0])
  sigma_e <- sigma_e[complete.cases(sigma_e)]
  w <- exp(-1 / sigma_e)
  w <- as.matrix(w / sum(w,na.rm=TRUE))
  tau_star <- tau_max_mat%*%w
  tau_star <- as.vector(tau_star)
  prob_vec <- numeric(length(tau_star))
  for (j in seq_along(tau_star)) {
    prob_vec[j] <- mean(delta[j, 1:which.min(abs(tau_grid - mean(tau_star, na.rm=TRUE)))] > 0, na.rm = TRUE)
  }
  decision <- mean(as.numeric(prob_vec > eta))
  results$Constrained_Restriction_Time <- list(summary=tibble(x = tau_star) %>%
                                                 summarise(
                                                   median = median(x, na.rm=TRUE),
                                                   q025  = quantile(x, 0.025, na.rm=TRUE),
                                                   q975  = quantile(x, 0.975, na.rm=TRUE),
                                                   mean = mean(x, na.rm=TRUE),
                                                   sd = sd(x, na.rm=TRUE)
                                                 ),
                                               tau_star = tibble(tau_star, prob_vec),
                                               decision = decision
  )
  #
  return(list(results=results, gridOmega=gridOmega))
  #
}

#############################################################
#############################################################
