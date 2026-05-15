read_data <- function(file_name){
  read_view(file_name)%>%
    clean_names()%>%
    select(o_net_soc_code, element_name, scale_name, data_value)%>%
    pivot_wider(names_from = scale_name, values_from = data_value)%>%
    mutate(Importance=10*(Importance-1)/4,
           Level=10*Level/7,
           score=sqrt(Importance*Level), #geometric mean of importance and level
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

plot_pca_variance <- function(pca_obj, n_comp = 10, digits = 1, text_size=3.5) {

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
      size = text_size
    ) +

    # Cumulative line + points
    geom_line(aes(y = Cumulative, group = 1), linewidth = 1) +
    geom_point(aes(y = Cumulative), size = 2) +

    # Cumulative labels (skip PC1)
    geom_text(
      data = df[-1, ],
      aes(y = Cumulative, label = Cum_lab),
      vjust = -0.8,
      size = text_size
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

offdiag <- function(M) {
  M[row(M) != col(M)]
}

outer_named <- function(a, b) {
  M <- a %o% b
  rownames(M) <- names(a)
  colnames(M) <- names(b)
  M
}

distance_plot <- function(P_obs, C, subtitle) {

  a <- rowSums(P_obs)
  b <- colSums(P_obs)

  df <- tibble(
    distance = as.vector(C),
    ratio = log(as.vector(P_obs) / as.vector(a %o% b))
  ) |>
    filter(distance > 0)

  ggplot(df, aes(distance, ratio)) +
    geom_jitter(alpha = 0.01) +
    geom_smooth(method="lm", se=FALSE, colour="red")+
#    scale_y_log10() +
    labs(
      x = "Occupational distance (normalized)",
      y = "Log excess transition probability (relative to independence)",
      title = NULL,
      subtitle=subtitle
    )
}

long_to_matrix <- function(tbbl){
  tbbl|>
    pivot_wider(id_cols = origin, names_from = noc, values_from = distance)|>
    column_to_rownames("origin")|>
    matrix()
}

weighted_cost <- function(s, h, h_weight){
   h_weight*h+(1-h_weight)*s
}

noc_prop_plot <- function(tbbl, inital_age, subsequent_age){
  tbbl|>
    filter((syear==2011 & age10==inital_age)| (syear==2021 & age10==subsequent_age))|>
    group_by(syear, age10)|>
    mutate(prop=count/sum(count))|>
    select(syear, noc_plus_title, prop)|>
    pivot_wider(id_cols = noc_plus_title, names_from = syear, values_from = prop)|>
    mutate(diff=`2021`-`2011`,
           abs_diff=abs(diff))|>
    slice_max(abs_diff, n=20)|>
    ggplot(aes(diff, fct_reorder(noc_plus_title, -diff)))+
    geom_col(alpha=.5)+
    scale_x_continuous(limits = c(-.075,.025))+
    labs(title=paste("Change in Occupational Proportions between",inital_age,"and",subsequent_age,"year olds (ten years later)"),
         x="Change in Proportion",
         y=NULL,
         caption="Source: LFS via RTRA")+
#    theme_minimal()+
    theme(
      plot.title.position = "plot"
    )
}


subset_rows <- function(C, idx) {
  C[idx, , drop = FALSE]
}

get_idx <- function(C, sub_regime_vec, label) {
  aligned <- sub_regime_vec[rownames(C)]
  aligned == label
}

subset_score_dist <- function(a){
  keep_rows <- names(a)
  score_dist[keep_rows,]
}


