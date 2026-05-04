
#############################################################
## run survival models
#############################################################

survModel <- function(data = df_data,
                      modelType = "BayesPara",
                      family = "weibull",
                      iter = 2000,
                      chains = 1, cores = 1,
                      control = list(adapt_delta = 0.95),
                      seed = 1234){
  #
  # this function is for data with binary treatment/exposure variable
  # model will be fitted separately for each treatment group
  #
  # data = with 3 cols: time, event, trt(binary, 0,1)
  # names(data) <- c("X", "delta", "A") # time, event, trt
  #
  names(data) <- c("X", "delta", "A") # time, event, trt
  trt_names <- names(table(data$A))
  time_grid <- seq(0, max(data$X), length.out = 100)
  df <- setNames(vector("list", length(modelType)), modelType)
  #
  if("KM" %in% modelType){
    #library(survival)
    # Fit survival model
    fit <- survfit(Surv(X, delta) ~ A, data = data)
    # Get summary at specified time points
    s <- summary(fit, times = time_grid)
    df_long <- data.frame(
      time = s$time,
      strata = s$strata,
      surv = s$surv,
      std_err = s$std.err
    )
    # Clean strata names (remove "A=")
    df_long$strata <- sub("A=", "", df_long$strata)
    # Reshape to wide format
    df_wide <- df_long |>
      pivot_wider(
        names_from = strata,
        values_from = c(surv, std_err),
        names_sep = "_"
      )
    df_wide$model = "KM"
    names(df_wide) = c("time","surv_0","surv_1","std_err_0","std_err_1","model")
    df$KM = df_wide
  }
  if("BayesNonPara" %in% modelType){
    #library(survival)
    #library(spBayesSurv)
    # ref: https://onlinelibrary.wiley.com/doi/full/10.1111/j.1541-0420.2008.01166.x
    fitnp0 <- anovaDDP(
      Surv(X, delta) ~ 1,
      data = subset(data, A == trt_names[1]),
      mcmc = list(nburn = 0.5*iter, nsave = iter,
                  nskip = 0, ndisplay = 100)
    )
    fitnp1 <- anovaDDP(
      Surv(X, delta) ~ 1,
      data = subset(data, A == trt_names[2]),
      mcmc = list(nburn = 0.5*iter, nsave = iter,
                  nskip = 0, ndisplay = 100)
    )
    # Posterior draws/summary/surv-prob
    psnp0 <- GetCurves(fitnp0,tgrid = time_grid, PLOT=FALSE)
    psnp0$sd <- (psnp0$Shatup - psnp0$Shatlow) / (2 * 1.96)
    psnp1 <- GetCurves(fitnp1,tgrid = time_grid, PLOT=FALSE)
    psnp1$sd <- (psnp1$Shatup - psnp1$Shatlow) / (2 * 1.96)
    # Combine side-by-side
    df_wide <- data.frame(
      time = time_grid,
      surv_0 = psnp0$Shat,
      surv_1 = psnp1$Shat,
      std_err_0 = psnp0$sd,
      std_err_1 = psnp1$sd
    )
    df_wide$model = "BayesNonPara"
    df$BayesNonPara = tibble(df_wide)
  }
  if("BayesPara" %in% modelType){
    #library(brms)
    # Fit model for A = 0
    fit0 <- brm(
      bf(X | cens(1 - delta) ~ 1),
      data = subset(data, A == trt_names[1]),
      family = family,
      chains = chains, cores = cores,
      iter = iter, seed = seed, control = control
    )
    # Fit model for A = 1
    fit1 <- brm(
      bf(X | cens(1 - delta) ~ 1),
      data = subset(data, A == trt_names[2]),
      family = family,
      chains = chains, cores = cores,
      iter = iter, seed = seed, control = control
    )
    # Posterior draws
    ps0 <- posterior_predict(fit0)
    ps0 <- sapply(time_grid, function(t) colMeans(ps0 > t))
    ps1 <- posterior_predict(fit1)
    ps1 <- sapply(time_grid, function(t) colMeans(ps1 > t))
    # Summaries
    surv0 <- apply(ps0, 2, mean)
    se0   <- apply(ps0, 2, sd)
    surv1 <- apply(ps1, 2, mean)
    se1   <- apply(ps1, 2, sd)
    # Combine side-by-side
    df_wide <- data.frame(
      time = time_grid,
      surv_0 = surv0,
      surv_1 = surv1,
      std_err_0 = se0,
      std_err_1 = se1
    )
    df_wide$model = "BayesPara"
    df$BayesPara = tibble(df_wide)
  }
  if("BayesMspline" %in% modelType){
    # ref: https://doi.org/10.1186/s12874-023-02094-1
    #library(survextrap)
    fitms0 <- survextrap(
      formula = Surv(X, delta) ~ 1,
      data = subset(data, A == trt_names[1]),
      chains = chains, cores = cores,
      iter = iter, seed = seed, control = control
    )
    fitms1 <- survextrap(
      formula = Surv(X, delta) ~ 1,
      data = subset(data, A == trt_names[2]),
      chains = chains, cores = cores,
      iter = iter, seed = seed, control = control
    )
    #
    psms0 <- survival(fitms0, t = time_grid,
                      summ_fns=list(mean=mean, sd=sd)) # ~quantile(.x, probs=c(0.025, 0.975))))
    psms1 <- survival(fitms1, t = time_grid,
                      summ_fns=list(mean=mean, sd=sd))
    # Combine side-by-side
    df_wide <- data.frame(
      time = time_grid,
      surv_0 = psms0$mean,
      surv_1 = psms1$mean,
      std_err_0 = psms0$sd,
      std_err_1 = psms1$sd
    )
    df_wide$model = "BayesMspline"
    df$BayesMspline = tibble(df_wide)
  }
  if(length(df) == 0){
    stop("can take model arguments: KM, BayesNonPara, BayesPara, BayesMspline, ...")
  }
  return(results=df)
}

#############################################################
#############################################################
