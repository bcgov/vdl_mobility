read_data <- function(file_name){
  read_view(file_name)%>%
    clean_names()%>%
    select(o_net_soc_code, element_name, scale_name, data_value)%>%
    pivot_wider(names_from = scale_name, values_from = data_value)%>%
    mutate(score=sqrt(Importance*Level), #geometric mean of importance and level
           #mutate(score=Level,
           category=(str_split(file_name,"\\.")[[1]][1]))%>%
    unite(element_name, category, element_name, sep=": ")%>%
    select(-Importance, -Level)
}

h_dist <- function(tbbl){
  tbbl |>
    mutate(teer_o = as.integer(str_sub(origin, 2, 2)),
           teer_d = as.integer(str_sub(destination, 2, 2)),
           teer_gap = abs(teer_o - teer_d),
           dist = case_when(origin == destination ~ 0,
                            str_sub(origin,1,4) == str_sub(destination,1,4) ~ 1,
                            str_sub(origin,1,3) == str_sub(destination,1,3) ~ 2,
                            str_sub(origin,1,1) == str_sub(destination,1,1) ~ 3 + teer_gap,
                            TRUE ~ 9))|>
    select(origin, destination, distance=dist)
}
offdiag <- function(M) M[row(M) != col(M)]

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

sinkhorn_aligned <- function(a, b, C, epsilon, solver, ...) {
  stopifnot(!is.null(names(a)))
  stopifnot(!is.null(names(b)))
  stopifnot(!is.null(rownames(C)))
  stopifnot(!is.null(colnames(C)))
  stopifnot(setequal(names(a), rownames(C)))
  stopifnot(setequal(names(b), colnames(C)))
  C <- C[names(a), names(b), drop = FALSE]
  stopifnot(identical(names(a), rownames(C))) #this really shouldn't happen now, but just in case.
  stopifnot(identical(names(b), colnames(C))) #this really shouldn't happen now, but just in case.
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

  # normalize
  a <- a / sum(a)
  b <- b / sum(b)

  n <- length(a)
  m <- length(b)

  # log marginals (handle zeros correctly)
  log_a <- ifelse(a > 0, log(a), -Inf)
  log_b <- ifelse(b > 0, log(b), -Inf)

  # log kernel
  logK <- -C / epsilon

  log_u <- rep(0, n)
  log_v <- rep(0, m)

  logsumexp <- function(x) {
    xmax <- max(x)
    xmax + log(sum(exp(x - xmax)))
  }

  for (iter in seq_len(max_iter)) {

    # update log_u
    for (i in seq_len(n)) {
      log_u[i] <- log_a[i] -
        logsumexp(logK[i, ] + log_v)
    }

    # update log_v
    for (j in seq_len(m)) {
      log_v[j] <- log_b[j] -
        logsumexp(logK[, j] + log_u)
    }

    # optional convergence check
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

  list(plan = P, log_u = log_u, log_v = log_v)
}


cumvar_explained <- function(pca_obj, x) {
  if (!inherits(pca_obj, "prcomp")) {
    stop("pca_obj must be a prcomp object.")
  }
  sdev2 <- pca_obj$sdev^2
  prop_var <- sdev2 / sum(sdev2)
  if (x > length(prop_var)) {
    stop("x exceeds number of principal components.")
  }
  sum(prop_var[seq_len(x)]) * 100
}

plot_pca_variance <- function(pca_obj, n_comp = 10, digits = 1) {

  if (!inherits(pca_obj, "prcomp")) {
    stop("Input must be a prcomp object.")
  }

  eigvals  <- pca_obj$sdev^2
  prop_var <- eigvals / sum(eigvals)
  cum_var  <- cumsum(prop_var)

  n <- min(n_comp, length(eigvals))

  percent_vals <- prop_var[seq_len(n)] * 100
  cum_vals     <- cum_var[seq_len(n)] * 100

  df <- data.frame(
    PC = factor(paste0("PC", seq_len(n)),
                levels = paste0("PC", seq_len(n))),
    Percent = percent_vals,
    Cumulative = cum_vals,
    Percent_lab = paste0(round(percent_vals, digits), "%"),
    Cum_lab = paste0(round(cum_vals, digits), "%")
  )

  library(ggplot2)

  ggplot(df, aes(x = PC)) +
    geom_col(aes(y = Percent), width = 0.7) +

    # Bar labels
    geom_text(
      aes(y = Percent, label = Percent_lab),
      vjust = -0.5,
      size = 3.5
    ) +

    # Cumulative line + points
    geom_line(aes(y = Cumulative, group = 1), linewidth = 1) +
    geom_point(aes(y = Cumulative), size = 2) +

    # Cumulative labels (skip PC1)
    geom_text(
      data = df[-1, ],
      aes(y = Cumulative, label = Cum_lab),
      vjust = -0.8,
      size = 3.5
    ) +

    scale_y_continuous(
      name = "Percent Variance Explained",
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.05))
    ) +

    labs(
      title = "PCA Variance Explained",
      x = "Principal Component"
    ) +

    theme_minimal()
}

my_dt <- function(tbbl, round_digits = 3) {
  num_cols <- names(tbbl)[vapply(tbbl, is.numeric, logical(1))]

  DT::datatable(
    tbbl,
    filter = "top",
    extensions = "Buttons",
    rownames = FALSE,
    options = list(
      columnDefs = list(list(className = "dt-center", targets = "_all")),
      paging = TRUE,
      scrollX = TRUE,
      scrollY = TRUE,
      searching = TRUE,
      ordering = TRUE,
      dom = "Btip",
      buttons = list(
        list(extend = "csv", filename = "education_specificity"),
        list(extend = "excel", filename = "education_specificity")
      ),
      pageLength = 25
    )
  ) |>
    DT::formatRound(columns = num_cols, digits = round_digits)
}

extract_margin <- function(tbbl, quoted_age, unquoted_column){
  tbbl|>
    filter(age_broad==quoted_age)|>
    ungroup()|>
    select(noc_plus_title, {{ unquoted_column }})|>
    deframe()
}



