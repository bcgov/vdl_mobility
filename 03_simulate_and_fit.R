library(tidyverse)
library(here)
library(janitor)
library(vroom)
library(matrixStats)
library(bcgovpond)
library(safesink)
library(conflicted)
conflicts_prefer(vroom::cols)
conflicts_prefer(vroom::col_double)
conflicts_prefer(vroom::col_character)
conflicts_prefer(dplyr::filter)

source(here("R", "other.R"))
true_eps <- 1 #The true temperature of the system
#normalized distance matrices-----------------------------------
skill_dist <- read_rds(here("out","skill_dist.rds"))
hier_dist <- read_rds(here("out","hier_dist.rds"))
binary_dist <- read_rds(here("out","binary_dist.rds"))
#get the noc list....

noc_list <- tibble(noc_plus_title=colnames(skill_dist))|>
  mutate(noc_5=str_sub(noc_plus_title, 1,5))|>
  arrange(noc_5)

lfs_files <- c(resolve_current("age1115p1.csv"),
               resolve_current("age1115p2.csv"),
               resolve_current("age2125p1.csv"),
               resolve_current("age2125p2.csv")
)

lfs_data <- vroom(lfs_files,
                      col_types = cols(
                        SYEAR = col_double(),
                        NOC_5 = col_character(),
                        AGE10 = col_character(),
                        `_COUNT_` = col_double()
                      ))|>
  na.omit()|>
  clean_names()|>
  filter(age10!="nwa",
         noc_5!="no_no",
         syear<max(syear))|>
  mutate(noc_5=if_else(noc_5 %in% paste0("000",11:15), "00018", noc_5))|>
  summarize(count=sum(count)/12, .by=c(syear, age10, noc_5))|>
  inner_join(noc_list)

a <- lfs_data|>
  filter(syear<2015,
         !age10 %in% c("15-24", "55-64"))|> #25-54 in start year
  summarize(count=sum(count), .by = c(noc_plus_title))|>
  mutate(prop=count/sum(count))|>
  select(noc_plus_title, prop)|>
  arrange((noc_plus_title))|>
  with(setNames(prop, noc_plus_title))

b <- lfs_data|>
  filter(syear>2020,
         !age10 %in% c("15-24", "25-34"))|> #35-64 in end year
  summarize(count=sum(count), .by = c(noc_plus_title))|>
  mutate(prop=count/sum(count))|>
  select(noc_plus_title, prop)|>
  arrange(noc_plus_title)|>
  with(setNames(prop, noc_plus_title))


#sanity check before proceeding: we are mixing distance matrices, they need to be in same order.
stopifnot(
  identical(rownames(skill_dist), rownames(hier_dist)),
  identical(rownames(skill_dist), rownames(binary_dist)),
  identical(colnames(skill_dist), colnames(hier_dist)),
  identical(colnames(skill_dist), colnames(binary_dist)),
  setequal(names(a), rownames(skill_dist)), #mismatched vectors ok at this stage align_transport_inputs() before fitting
  setequal(names(b), colnames(skill_dist)) #mismatched vectors ok at this stage align_transport_inputs() before fitting
)

# All objects must share the same ordering: master_ids----------------
if (exists("master_ids", inherits = FALSE)) {
  stop("master_ids already defined — restart R for clean slate")
}

master_ids <- Reduce(intersect, list(
  names(a),
  names(b),
  rownames(skill_dist),
  colnames(skill_dist),
  rownames(hier_dist),
  colnames(hier_dist),
  rownames(binary_dist),
  colnames(binary_dist)
))
lockBinding("master_ids", environment())
#reorder everything-------------------------------

skill_dist <- skill_dist[master_ids, master_ids]
hier_dist  <- hier_dist[master_ids, master_ids]
binary_dist <- binary_dist[master_ids, master_ids]

a <- a[master_ids]
b <- b[master_ids]

global_simulations <- list(
  `DGP: 75% Skill / 25% Hierarchy` = list(
    a = a,
    b = b,
    log_K = -(.75 * skill_dist + .25 * hier_dist)/true_eps
  ),
  `DGP: 50% Skill / 50% Hierarchy` = list(
    a = a,
    b = b,
    log_K = -(.5 * skill_dist + .5 * hier_dist)/true_eps
  ),
  `DGP: 25% Skill / 75% Hierarchy` = list(
    a = a,
    b = b,
    log_K = -(.25 * skill_dist + .75 * hier_dist)/true_eps
  )
)

global_P_fake <- map(global_simulations, \(sim) {
  stopifnot(
    identical(names(sim$a), rownames(sim$log_K)),
    identical(names(sim$b), colnames(sim$log_K))
  )
  sinkhorn_log(sim$a, sim$b, sim$log_K)$plan
})

#have the fake data, ready to fit models

global_model_fits <- crossing(
  dgp = names(global_P_fake),
  `Cost Matrix` = c("Skill","Hierarchy", "Binary"),
  Temperature = 2^seq(-1,1,.1)) |>
  mutate(
    P_obs = map(dgp, ~global_P_fake[[.x]]),
    C = map(`Cost Matrix`, \(cm) {switch(cm, "Skill" = skill_dist,
                                             "Hierarchy" = hier_dist,
                                             "Binary" = binary_dist)}),
    log_K = map2(C, Temperature, \(C, temp) -C / temp),
    a = map(P_obs, rowSums),
    b = map(P_obs, colSums),
    P_hat = pmap(list(a, b, log_K),\(a, b, log_K) {
        stopifnot(
          identical(names(a), rownames(log_K)),
          identical(names(b), colnames(log_K)))
        sinkhorn_log(a, b, log_K)$plan
      }
    ),
    P_ind = map2(a, b, ~ .x %o% .y),
    dgp = factor(dgp,levels = c("DGP: 75% Skill / 25% Hierarchy",
                                "DGP: 50% Skill / 50% Hierarchy",
                                "DGP: 25% Skill / 75% Hierarchy")))|>
  select(dgp, Temperature, C, `Cost Matrix`, P_obs, P_hat, P_ind)


write_rds(global_model_fits, here("out", "global_model_fits.rds"))


# market sub_regimes------------------------------------

noc_specificity <- read_rds(here("out", "noc_specificity.rds")) |>
  arrange(noc) |>
  select(noc_plus_title = noc, sub_regime)

sub_regime_vec <- setNames(noc_specificity$sub_regime,
                           noc_specificity$noc_plus_title)

# enforce canonical ordering
sub_regime_vec <- sub_regime_vec[master_ids]

# ----------------------------
# core validity checks only
# ----------------------------

stopifnot(!anyNA(sub_regime_vec))

stopifnot(all(c(
  "Horizontal (Skill)",
  "Vertical (Hierarchy)",
  "Minimal (Binary)"
) %in% sub_regime_vec))

# ----------------------------
# strict alignment invariants
# (only things that are NOT implied by construction)
# ----------------------------

stopifnot(identical(names(a), master_ids))
stopifnot(identical(names(b), master_ids))
stopifnot(identical(rownames(skill_dist), master_ids))
stopifnot(identical(colnames(skill_dist), master_ids))
stopifnot(identical(rownames(hier_dist), master_ids))
stopifnot(identical(colnames(hier_dist), master_ids))
stopifnot(identical(rownames(binary_dist), master_ids))
stopifnot(identical(colnames(binary_dist), master_ids))

#a, b, and distances all share ordering (as defined by master_ids)

#subset the a vector, and keep appropriate rows of C matrices

a_skill_sub_regime <- a[sub_regime_vec == "Horizontal (Skill)"]
a_hier_sub_regime <- a[sub_regime_vec == "Vertical (Hierarchy)"]
a_binary_sub_regime <- a[sub_regime_vec == "Minimal (Binary)"]

#Horizontal Regime
C_horizontal_skill <- skill_dist[sub_regime_vec == "Horizontal (Skill)", , drop = FALSE]
C_horizontal_hier <- hier_dist[sub_regime_vec == "Horizontal (Skill)", , drop = FALSE]
C_horizontal_binary <- binary_dist[sub_regime_vec == "Horizontal (Skill)", , drop = FALSE]

#Vertical Regime
C_vertical_hier <- hier_dist[sub_regime_vec == "Vertical (Hierarchy)", , drop = FALSE]
C_vertical_skill <- skill_dist[sub_regime_vec == "Vertical (Hierarchy)", , drop = FALSE]
C_vertical_binary <- binary_dist[sub_regime_vec == "Vertical (Hierarchy)", , drop = FALSE]

#Minimal Regime
C_minimal_binary <- binary_dist[sub_regime_vec == "Minimal (Binary)", , drop = FALSE]
C_minimal_skill <- skill_dist[sub_regime_vec == "Minimal (Binary)", , drop = FALSE]
C_minimal_hier <- hier_dist[sub_regime_vec == "Minimal (Binary)", , drop = FALSE]


sub_regime_simulations <- list(`Minimal` = list( a = a_binary_sub_regime, b = b, log_K = -C_minimal_binary/true_eps),
                                `Vertical` = list( a = a_hier_sub_regime, b = b, log_K = -C_vertical_hier/true_eps),
                                `Horizontal` = list( a = a_skill_sub_regime, b = b, log_K = -C_horizontal_skill/true_eps))

sub_regime_P_fake <- map(sub_regime_simulations, \(sim) sinkhorn_log(sim$a, sim$b, sim$log_K)$plan)

#sub regime model fitting-----------------------------------

dgp <- c("Minimal", "Vertical", "Horizontal")
cost <- c("Binary", "Hier", "Skill")
Temperature <- 2^seq(-1, 1, .1)

sub_regime_grid <- tidyr::crossing(
  dgp,
  cost,
  Temperature)|>
  mutate(`Cost Matrix`=paste(dgp, cost, sep="_"))

a_list <- list(
  "Minimal" = a_binary_sub_regime,
  "Vertical" = a_hier_sub_regime,
  "Horizontal" = a_skill_sub_regime
)

C_list <- list(
  "Horizontal_Skill" = C_horizontal_skill,
  "Horizontal_Hier" = C_horizontal_hier,
  "Horizontal_Binary" = C_horizontal_binary,
  "Vertical_Hier" = C_vertical_hier,
  "Vertical_Skill" = C_vertical_skill,
  "Vertical_Binary" = C_vertical_binary,
  "Minimal_Binary" = C_minimal_binary,
  "Minimal_Skill" = C_minimal_skill,
  "Minimal_Hier" = C_minimal_hier
)

sub_regime_grid <- sub_regime_grid |>
  mutate(
    a = map(dgp, ~ a_list[[.x]]),
    b = map(dgp, ~ b),
    C = map(`Cost Matrix`, ~ C_list[[.x]]),
    P_obs =map(dgp, ~sub_regime_P_fake[[.x]])
  )

stopifnot(
  all(map2_lgl(sub_regime_grid$a, sub_regime_grid$C, ~ length(.x) == nrow(.y))),
  all(map2_lgl(sub_regime_grid$b, sub_regime_grid$C, ~ length(.x) == ncol(.y)))
)

sub_regime_model_fits <- sub_regime_grid |>
  mutate(
    P_hat = pmap(
      list(a, b, C, Temperature),
      \(a, b, C, temp){
      log_K <- -C/temp
      stopifnot(
        identical(names(a), rownames(log_K)),
        identical(names(b), colnames(log_K))
      )
      sinkhorn_log(a, b, log_K)$plan
      }
      ),
    P_ind = map2(a, b, ~ .x %o% .y)
  )|>
  select(dgp, Temperature, C, `Cost Matrix`=cost, P_obs, P_hat, P_ind)

write_rds(sub_regime_model_fits, here("out", "sub_regime_model_fits.rds"))

