

#############################################################
## model pipe-line
#############################################################

survRMST <- function(
    data,
    modelType = c("KM"),
    family = "weibull",
    iter = 2000,
    chains = 1,
    cores = 1,
    control = list(adapt_delta = 0.95),
    seed = 1234,
    prior_expert_list = list(c(5,6), c(10,12), c(8.5,10)),
    eta = 0.95
)
{
  # Fit survival models
  fit <- survModel(
    data = data,
    modelType = modelType,
    family = family,
    iter = iter,
    chains = chains,
    cores = cores,
    control = control,
    seed = seed
  )
  object <- fit[[modelType]]
  # Bayesian decision (tau selection)
  decision <- BayesDecision_RMST(
    result = object,
    prior_expert_list = prior_expert_list,
    seed = seed,
    eta = eta
  )
  # RMST
  rmst <- compute_rmst_table(object=object, decision=decision)
  rmst$model <- modelType
  class(rmst) <- c("rowwise_df","tbl_df","tbl","data.frame","BayesRMST")
  #
  out <- NULL
  out$RMST <- rmst
  out$model_fit <- fit
  out$decision <- decision
  class(out) <- "BayesRMST"
  out
  #
}

#############################################################
#############################################################
