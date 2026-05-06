library(tidyverse)
library(here)
library(janitor)
library(safesink)
library(gganimate)
library(furrr)
library(future)
library(conflicted)

sim_ot_plots <- list()

compute_self_cost <- function(P, C, epsilon) {
  n <- nrow(P)
  log_K <- -C / epsilon
  out <- numeric(n)

  for (i in seq_len(n)) {
    a <- P[i, ]
    sa <- sum(a)

    if (!is.finite(sa) || sa == 0) {
      out[i] <- NA_real_
      next
    }

    a <- a / sa
    P_aa <- sinkhorn_log(a, a, log_K)$plan
    out[i] <- sum(C * P_aa)
  }

  out
}

plan(multisession, workers = parallelly::availableCores(omit = 1))

global_model_fits <- read_rds(here("out", "global_model_fits.rds"))|>
  mutate(
    aa_cost = future_pmap(list(P_obs, C, Temperature), compute_self_cost, .progress = TRUE),
    bb_cost = future_pmap(list(P_hat, C, Temperature), compute_self_cost, .progress = TRUE)
  )

write_rds(global_model_fits, here("out","global_fit_with_self_cost.rds"))


plan(sequential)

