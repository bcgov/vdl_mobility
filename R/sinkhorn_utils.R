# sinkhorn_utils.R
# Utility functions for entropy-regularized optimal transport (Sinkhorn)
# Patched version with improved numerical stability and robustness.
# Key improvements:
#  - Safe handling of zero marginals (prevents NaNs in log updates)
#  - KL score renormalization after flooring
#  - Optional matrixStats acceleration when available
#  - Extra finite-value diagnostics


kl_score <- function(P_obs, P_hat, tiny = 1e-15) {

  # Normalize distributions
  P_obs <- P_obs / sum(P_obs)
  P_hat <- P_hat / sum(P_hat)

  # Bound away from zero
  P_hat <- pmax(P_hat, tiny)

  # Renormalize after flooring
  P_hat <- P_hat / sum(P_hat)

  idx <- P_obs > 0

  sum(P_obs[idx] * log(P_obs[idx] / P_hat[idx]))
}

log_sum_exp <- function(x) {
  xmax <- max(x)
  if (is.infinite(xmax)) return(-Inf)
  xmax + log(sum(exp(x - xmax)))
}

# Fast log-sum-exp helpers using matrixStats if available
row_lse <- function(M) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowLogSumExps(M)
  } else {
    apply(M, 1, log_sum_exp)
  }
}

col_lse <- function(M) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colLogSumExps(M)
  } else {
    apply(M, 2, log_sum_exp)
  }
}

sinkhorn_aligned <- function(a, b, C, epsilon, solver, ...) {

  stopifnot(!is.null(names(a)))
  stopifnot(!is.null(names(b)))
  stopifnot(!is.null(rownames(C)))
  stopifnot(!is.null(colnames(C)))
  stopifnot(setequal(names(a), rownames(C)))
  stopifnot(setequal(names(b), colnames(C)))

  C <- C[names(a), names(b), drop = FALSE]

  stopifnot(identical(names(a), rownames(C)))
  stopifnot(identical(names(b), colnames(C)))

  solver(a, b, C, epsilon, ...)
}

sinkhorn_log <- function(a, b, C, epsilon,
                         max_iter = 5000,
                         tol = 1e-9,
                         verbose = FALSE) {

  stopifnot(is.matrix(C))
  stopifnot(length(a) == nrow(C))
  stopifnot(length(b) == ncol(C))
  stopifnot(sum(a) > 0, sum(b) > 0)

  # normalize marginals
  a <- a / sum(a)
  b <- b / sum(b)

  n <- length(a)
  m <- length(b)

  log_a <- ifelse(a > 0, log(a), -Inf)
  log_b <- ifelse(b > 0, log(b), -Inf)

  logK <- -C / epsilon

  log_u <- rep(0, n)
  log_v <- rep(0, m)

  for (iter in seq_len(max_iter)) {

    row_terms <- logK + matrix(log_v, n, m, byrow = TRUE)
    lse_rows <- row_lse(row_terms)

    log_u_new <- log_a - lse_rows
    log_u_new[a == 0] <- -Inf
    log_u <- log_u_new

    col_terms <- logK + matrix(log_u, n, m, byrow = FALSE)
    lse_cols <- col_lse(col_terms)

    log_v_new <- log_b - lse_cols
    log_v_new[b == 0] <- -Inf
    log_v <- log_v_new

    if (iter %% 50 == 0) {

      logP <- outer(log_u, log_v, "+") + logK
      P <- exp(logP)

      err <- max(
        max(abs(rowSums(P) - a)),
        max(abs(colSums(P) - b))
      )

      if (verbose) cat("Iter", iter, "err:", err, "\n")

      if (err < tol) break
    }
  }

  logP <- outer(log_u, log_v, "+") + logK
  P <- exp(logP)
  P <- P / sum(P)

  list(plan = P, log_u = log_u, log_v = log_v)
}

make_safe <- function(x, eps = 1e-12) {
  x <- pmax(x, eps)
  x / sum(x)
}

check_transport <- function(P, a, b, tol = 1e-8) {

  a <- a / sum(a)
  b <- b / sum(b)

  row_err <- max(abs(rowSums(P) - a))
  col_err <- max(abs(colSums(P) - b))
  mass_err <- abs(sum(P) - 1)

  out <- data.frame(
    mass      = sum(P),
    row_err   = row_err,
    col_err   = col_err,
    mass_err  = mass_err,
    min_value = min(P),
    any_neg   = any(P < -tol),
    any_nan   = any(!is.finite(P))
  )

  print(out)

  invisible(out)
}

compare_transport <- function(P1, P2) {

  stopifnot(all(dim(P1) == dim(P2)))

  cor_val <- suppressWarnings(cor(as.vector(P1), as.vector(P2)))

  out <- data.frame(
    max_abs_diff = max(abs(P1 - P2)),
    rmse = sqrt(mean((P1 - P2)^2)),
    correlation = cor_val,
    mass_diff = abs(sum(P1) - sum(P2))
  )

  print(out)

  invisible(out)
}
