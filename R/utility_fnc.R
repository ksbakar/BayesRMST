
#############################################################
#############################################################

print.BayesRMST <- function(x, ...){
  print(x$RMST)
}

#############################################################

plot.BayesRMST <- function(x,
                           survival = FALSE,
                           omega = FALSE,
                           rmst = FALSE,
                           CI = TRUE,
                           ...) {
  #
  # If nothing selected → default to RMST
  if (!survival && !omega && !rmst) {
    rmst <- TRUE
  }
  p1 <- p2 <- p3 <- NULL
  if (survival) {
    p1 <- plot_surv(x$model_fit[[1]], CI = CI)
  }
  if (omega) {
    p2 <- plot_omega(x$decision$gridOmega)
  }
  if (rmst) {
    rmst_tbl <- compute_rmst_table(
      object = x$model_fit[[1]],
      decision = x$decision
    )
    p3 <- plot_rmst(
      object = x$model_fit[[1]],
      decision = x$decision,
      rmst_results = rmst_tbl
    )
  }
  # collect non-null plots
  plots <- list(p1, p2, p3)
  plots <- plots[!sapply(plots, is.null)]
  n <- length(plots)
  # If only one plot: return it
  if (n == 1) {
    return(plots[[1]])
  }
  # If all three: custom layout
  if (n == 3) {
    return((p1 + p2) / p3)
  }
  # If two: side by side
  if (n == 2) {
    return(plots[[1]] + plots[[2]])
  }
}

#############################################################

plot_surv <- function(object, CI=TRUE){
  #
  if(isTRUE(CI)){
    surv_long <- object %>%
      pivot_longer(
        cols = c(surv_0, surv_1, std_err_0, std_err_1),
        names_to = c(".value", "group"),
        names_pattern = "(surv|std_err)_(.)"
      ) %>%
      mutate(
        lower = pmax(0, surv - 1.96 * std_err),
        upper = pmin(1, surv + 1.96 * std_err)
      )
    #
    ggplot(surv_long, aes(x = time, y = surv, color = group, fill = group)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      labs(
        x = "Years",
        y = "Survival Probability",
        color = "",
        fill = ""
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
      )
  }
  else{
    surv_long <- object %>%
      pivot_longer(cols = c(surv_0, surv_1),
                   names_to = "group",
                   values_to = "surv")
    ggplot() +
      geom_line(
        data = surv_long,
        aes(x = time, y = surv, color = group),
        linewidth = 1
      ) +
      labs(
        x = "Years",
        y = "Survival Probability",
        color = ""
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
      )
  }
}

#############################################################

plot_omega <- function(object){
  ggplot(object, aes(x = time, group = model, color = model, fill = model)) +
    # shaded ribbon for 95% interval
    geom_ribbon(aes(ymin = omega_q025, ymax = omega_q975), alpha = 0.2, color = NA) +
    # median line
    geom_line(aes(y = omega_median), size = 1) +
    labs(
      title = "Posterior signal-to-noise",
      x = "Years",
      y = expression(Omega(tau)),
      color = "",
      fill = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
}
#plot_omega(bb$gridOmega) # from decision

#############################################################

compute_rmst_table <- function(object, decision) {
  # reshape to long
  object_long <- object %>%
    pivot_longer(
      cols = c(surv_0, surv_1, std_err_0, std_err_1),
      names_to = c(".value", "group"),
      names_pattern = "(surv|std_err)_(\\d+)"
    ) %>%
    mutate(group = as.integer(group))
  # restriction summaries
  restr_df <- bind_rows(
    decision$results$Unconstrained_Restriction_Time$summary %>%
      mutate(type = "Unconstrained"),
    decision$results$Constrained_Restriction_Time$summary %>%
      mutate(type = "Constrained")
  )
  # trapezoidal RMST helper
  compute_rmst_inner <- function(df, tau) {
    df <- df %>%
      arrange(time) %>%
      filter(time <= tau)
    if (nrow(df) == 0) return(NA_real_)
    if (max(df$time) < tau) {
      last_row <- df[nrow(df), ]
      last_row$time <- tau
      df <- bind_rows(df, last_row)
    }
    sum(diff(df$time) * (head(df$surv, -1) + tail(df$surv, -1)) / 2)
  }
  # compute RMST table
  rmst_results <- restr_df %>%
    rowwise() %>%
    do({
      tau <- .$median
      type <- .$type
      object_long %>%
        group_by(group) %>%
        summarise(
          rmst = compute_rmst_inner(cur_data(), tau),
          .groups = "drop"
        ) %>%
        mutate(type = type, tau = tau)
    }) %>%
    bind_rows()
  return(rmst_results)
}

#############################################################

plot_rmst <- function(object, decision, rmst_results) {
  # reshape to long
  object_long <- object %>%
    pivot_longer(
      cols = c(surv_0, surv_1, std_err_0, std_err_1),
      names_to = c(".value", "group"),
      names_pattern = "(surv|std_err)_(\\d+)"
    ) %>%
    mutate(group = as.integer(group))
  # restriction summaries
  restr_df <- bind_rows(
    decision$results$Unconstrained_Restriction_Time$summary %>%
      mutate(type = "Unconstrained"),
    decision$results$Constrained_Restriction_Time$summary %>%
      mutate(type = "Constrained")
  )
  # prepare plotting data
  surv_long <- object_long %>%
    mutate(group = factor(group))
  restr_plot <- restr_df %>%
    rename(constraint = type, tau = median)
  surv_long <- surv_long %>%
    crossing(restr_plot)
  shade_df <- surv_long %>%
    filter(time <= tau)
  rmst_labels <- rmst_results %>%
    rename(constraint = type) %>%
    mutate(
      group = factor(group),
      label = paste0("RMST (g=", group, ") = ", round(rmst, 2)),
      y_pos = ifelse(group == "0", 0.2, 0.1)
    )
  pp <- ggplot() +
    geom_area(
      data = shade_df,
      aes(x = time, y = surv, fill = group),
      alpha = 0.25,
      position = "identity"
    ) +
    geom_line(
      data = surv_long,
      aes(x = time, y = surv, color = group),
      linewidth = 1
    ) +
    geom_vline(
      data = restr_plot,
      aes(xintercept = tau),
      linetype = "dashed",
      color = "black"
    ) +
    geom_rect(
      data = restr_plot,
      aes(
        xmin = q025,
        xmax = q975,
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = 0.08,
      fill = "grey50",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = rmst_labels,
      aes(
        x = tau * 0.1,
        y = y_pos,
        label = label
      ),
      size = 3,
      hjust = 0
    ) +
    facet_grid(~ constraint) +
    labs(
      x = "Years",
      y = "Survival Probability",
      color = "",
      fill = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(pp)
}

#############################################################

compute_rmst_plot <- function(object, decision, plot = FALSE) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  # reshape to long
  object_long <- object %>%
    pivot_longer(
      cols = c(surv_0, surv_1, std_err_0, std_err_1),
      names_to = c(".value", "group"),
      names_pattern = "(surv|std_err)_(\\d+)"
    ) %>%
    mutate(group = as.integer(group))
  # restriction summaries
  restr_df <- bind_rows(
    decision$results$Unconstrained_Restriction_Time$summary %>%
      mutate(type = "Unconstrained"),
    decision$results$Constrained_Restriction_Time$summary %>%
      mutate(type = "Constrained")
  )
  # trapezoidal RMST
  compute_rmst_inner <- function(df, tau) {
    df <- df %>%
      arrange(time) %>%
      filter(time <= tau)
    if (nrow(df) == 0) return(NA_real_)
    if (max(df$time) < tau) {
      last_row <- df[nrow(df), ]
      last_row$time <- tau
      df <- bind_rows(df, last_row)
    }
    sum(diff(df$time) * (head(df$surv, -1) + tail(df$surv, -1)) / 2)
  }
  # compute RMST
  rmst_results <- restr_df %>%
    rowwise() %>%
    do({
      tau <- .$median
      type <- .$type

      object_long %>%
        group_by(group) %>%
        summarise(
          rmst = compute_rmst_inner(cur_data(), tau),
          .groups = "drop"
        ) %>%
        mutate(type = type, tau = tau)
    }) %>%
    bind_rows()
  # optional plot
  if (plot) {
    # prepare plotting data
    surv_long <- object_long %>%
      mutate(group = factor(group))
    restr_plot <- restr_df %>%
      rename(constraint = type, tau = median)
    surv_long <- surv_long %>%
      crossing(restr_plot)
    shade_df <- surv_long %>%
      filter(time <= tau)
    rmst_labels <- rmst_results %>%
      rename(constraint = type) %>%
      mutate(
        group = factor(group),
        label = paste0("RMST (g=", group, ") = ", round(rmst, 2)),
        y_pos = ifelse(group == "0", 0.2, 0.1)
      )
    pp <- ggplot() +
      geom_area(
        data = shade_df,
        aes(x = time, y = surv, fill = group),
        alpha = 0.25,
        position = "identity"
      ) +
      geom_line(
        data = surv_long,
        aes(x = time, y = surv, color = group),
        linewidth = 1
      ) +
      geom_vline(
        data = restr_plot,
        aes(xintercept = tau),
        linetype = "dashed",
        color = "black"
      ) +
      geom_rect(
        data = restr_plot,
        aes(
          xmin = q025,
          xmax = q975,
          ymin = -Inf,
          ymax = Inf
        ),
        alpha = 0.08,
        fill = "grey50",
        inherit.aes = FALSE
      ) +
      geom_text(
        data = rmst_labels,
        aes(
          x = tau * 0.1,
          y = y_pos,
          label = label
        ),
        size = 3,
        hjust = 0
      ) +
      facet_grid( ~ constraint) +
      labs(
        x = "Years",
        y = "Survival Probability",
        color = "",
        fill = ""
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom"
      )
    print(pp)
  }
  #
  return(rmst_results)
  #
}

#############################################################
#############################################################
