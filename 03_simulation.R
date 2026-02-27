library(reticulate)
ot <- import("ot")

kl_score <- function(P_obs, P_hat, tiny = 1e-15) {

  # Normalize both matrices to ensure they are valid probability distributions.
  # This guards against small numerical drift from earlier computations.
  P_obs <- P_obs / sum(P_obs)
  P_hat <- P_hat / sum(P_hat)

  # Entropic OT should produce strictly positive P_hat,
  # but numerical underflow can produce exact zeros.
  # If P_hat == 0 where P_obs > 0, KL becomes infinite.
  # We therefore bound P_hat away from zero by a tiny constant.
  P_hat <- pmax(P_hat, tiny)

  # Terms where P_obs == 0 contribute 0 to KL by definition:
  #   lim_{x→0} x log(x / y) = 0
  # However, computing 0 * log(0 / y) directly gives NaN.
  # So we only sum over entries where P_obs > 0.
  idx <- P_obs > 0

  # Compute KL divergence:
  #   KL(P_obs || P_hat) == Σ P_obs * log(P_obs / P_hat)
  # This evaluates cross-entropy difference while focusing
  # exclusively on bilateral substitution patterns.
  sum(P_obs[idx] * log(P_obs[idx] / P_hat[idx]))
}

#helper function avoid overflow
log_sum_exp <- function(x) {
  mx <- max(x)
  if (!is.finite(mx)) return(mx)
  mx + log(sum(exp(x - mx)))
}

sinkhorn_log <- function(a, b, C, epsilon, max_iter = 5000, tol = 1e-9, verbose = FALSE) {
  # validate
  stopifnot(is.matrix(C))
  stopifnot(length(a) == nrow(C))
  stopifnot(length(b) == ncol(C))
  # Normalize
  a <- a / sum(a)
  b <- b / sum(b)
  #dimensions
  n <- length(a)
  m <- length(b)
  # K = exp(-C / epsilon) can underflow when C large or epsilon small: use log scale and exponentiate at the very end.
  logK <- -C / epsilon
  # In standard Sinkhorn: P = diag(u) * K * diag(v)
  # In log space: logP = log_u + logK + log_v
  # We initialize log_u and log_v at zero, corresponding to u = v = 1.
  log_u <- rep(0, n)
  log_v <- rep(0, m)
  # Sinkhorn iterations: In log-domain, these become additive updates.
  for (iter in seq_len(max_iter)) {
    log_u_prev <- log_u
    for (i in seq_len(n)) {# Update log_u (row normalization)
      log_u[i] <- log(a[i]) - log_sum_exp(logK[i, ] + log_v)
    }
    for (j in seq_len(m)) { #Update log_v (column normalization)
      log_v[j] <- log(b[j]) - log_sum_exp(logK[, j] + log_u)
    }
    if (max(abs(log_u - log_u_prev)) < tol) {#Convergence in scaling factors implies convergence in P.
      if (verbose) cat("Converged in", iter, "iterations\n")
      break
    }
  }
  # Recover the transport plan: logP = log_u + logK + log_v
  logP <- outer(log_u, log_v, "+") + logK
  P <- exp(logP)
  P <- P / sum(P) # Minor rounding error can accumulate in iterations.
  list(
    plan = P,
    iterations = iter
  )
}


P_pot <- ot$sinkhorn(a, b, C, epsilon)




