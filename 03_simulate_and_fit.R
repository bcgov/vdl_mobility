library(tidyverse)
library(here)
library(janitor)
library(vroom)
library(matrixStats)
library(bcgovpond)
library(safesink) #to install pak::pak("bcgov/safesink")
library(conflicted)
conflicts_prefer(vroom::cols)
conflicts_prefer(vroom::col_double)
conflicts_prefer(vroom::col_character)
conflicts_prefer(dplyr::filter)

source(here("R", "other.R"))
true_eps <- 1 #The true epsilon of the system
#normalized distance matrices-----------------------------------
skill_dist <- read_rds(here("out","skill_dist.rds"))
hier_dist <- read_rds(here("out","hier_dist.rds"))
binary_dist <- read_rds(here("out","binary_dist.rds"))
score_dist <- (hier_dist+skill_dist+binary_dist)/3 #average of the 3 distances for scoring
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
score_dist <- score_dist[master_ids, master_ids]

a <- a[master_ids]
b <- b[master_ids]

global_simulations <- list(
  `DGP: 75% Skill / 25% Hierarchy` = list(
    a = a,
    b = b,
    C = .75 * skill_dist + .25 * hier_dist
  ),
  `DGP: 50% Skill / 50% Hierarchy` = list(
    a = a,
    b = b,
    C = .5 * skill_dist + .5 * hier_dist
  ),
  `DGP: 25% Skill / 75% Hierarchy` = list(
    a = a,
    b = b,
    C = .25 * skill_dist + .75 * hier_dist
  )
)

global_P_fake <- map(global_simulations, \(sim) {
  stopifnot(
    identical(names(sim$a), rownames(sim$C)),
    identical(names(sim$b), colnames(sim$C))
  )
  sinkhorn_log(sim$a, sim$b, sim$C, true_eps)$plan
})

#have the fake data, ready to fit models

global_model_fits <- crossing(
  dgp = names(global_P_fake),
  `Cost Matrix` = c("Skill","Hierarchy", "Binary"),
  epsilon = 2^seq(-1,1,.1)) |>
  mutate(
    P_obs = map(dgp, ~global_P_fake[[.x]]),
    C = map(`Cost Matrix`, \(cm) {switch(cm, "Skill" = skill_dist,
                                             "Hierarchy" = hier_dist,
                                             "Binary" = binary_dist)}),
    a = map(P_obs, rowSums),
    b = map(P_obs, colSums),
    P_hat = pmap(list(a, b, C, epsilon),\(a, b, C, e) {
        stopifnot( identical(names(a), rownames(C)), identical(names(b), colnames(C)))
        sinkhorn_log(a, b, C, e)$plan
      }
    ),
    P_ind = map2(a, b, ~ .x %o% .y),
    dgp = factor(dgp, levels = c("DGP: 75% Skill / 25% Hierarchy",
                                "DGP: 50% Skill / 50% Hierarchy",
                                "DGP: 25% Skill / 75% Hierarchy")))|>
  select(dgp, epsilon, C, `Cost Matrix`, P_obs, P_hat, P_ind)|>
  mutate(score_dist = list(score_dist))

write_rds(global_model_fits, here("out", "global_model_fits.rds"))


# MARKET SUB REGIMES------------------------------------

noc_edu_spec <- read_rds(here("out", "noc_edu_spec.rds")) |>
  arrange(noc) |>
  select(noc_plus_title = noc, sub_regime)

sub_regime_vec <- setNames(noc_edu_spec$sub_regime,
                           noc_edu_spec$noc_plus_title)

# enforce canonical ordering
sub_regime_vec <- sub_regime_vec[master_ids]

# ----------------------------
# core validity checks only
# ----------------------------

stopifnot(!anyNA(sub_regime_vec))

stopifnot(all(c(
  "Skill",
  "Hierarchy",
  "Binary"
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

a_skill_sub_regime <- a[sub_regime_vec == "Skill"]
a_hierarchy_sub_regime <- a[sub_regime_vec == "Hierarchy"]
a_binary_sub_regime <- a[sub_regime_vec == "Binary"]

#Skill Regime
C_skill_skill <- skill_dist[sub_regime_vec == "Skill", , drop = FALSE]
C_skill_hierarchy <- hier_dist[sub_regime_vec == "Skill", , drop = FALSE]
C_skill_binary <- binary_dist[sub_regime_vec == "Skill", , drop = FALSE]

#Hierarchy Regime
C_hierarchy_hierarchy <- hier_dist[sub_regime_vec == "Hierarchy", , drop = FALSE]
C_hierarchy_skill <- skill_dist[sub_regime_vec == "Hierarchy", , drop = FALSE]
C_hierarchy_binary <- binary_dist[sub_regime_vec == "Hierarchy", , drop = FALSE]

#Binary Regime
C_binary_binary <- binary_dist[sub_regime_vec == "Binary", , drop = FALSE]
C_binary_skill <- skill_dist[sub_regime_vec == "Binary", , drop = FALSE]
C_binary_hierarchy <- hier_dist[sub_regime_vec == "Binary", , drop = FALSE]


sub_regime_simulations <- list(`Binary` = list( a = a_binary_sub_regime, b = b, C = C_binary_binary),
                               `Hierarchy` = list( a = a_hierarchy_sub_regime, b = b, C = C_hierarchy_hierarchy),
                               `Skill` = list( a = a_skill_sub_regime, b = b, C = C_skill_skill)
                               )

sub_regime_P_fake <- map(sub_regime_simulations, \(sim) sinkhorn_log(sim$a, sim$b, sim$C, true_eps)$plan)

#sub regime model fitting-----------------------------------

dgp <-  c("Binary", "Hierarchy", "Skill")
cost <- dgp
epsilon <- 2^seq(-1, 1, .1)

sub_regime_grid <- tidyr::crossing(
  dgp,
  cost,
  epsilon)|>
  mutate(`Cost Matrix`=paste(dgp, cost, sep="_"))

a_list <- list(
  "Binary" = a_binary_sub_regime,
  "Hierarchy" = a_hierarchy_sub_regime,
  "Skill" = a_skill_sub_regime
)

C_list <- list(
  "Skill_Skill" = C_skill_skill,
  "Skill_Hierarchy" = C_skill_hierarchy,
  "Skill_Binary" = C_skill_binary,
  "Hierarchy_Hierarchy" = C_hierarchy_hierarchy,
  "Hierarchy_Skill" = C_hierarchy_skill,
  "Hierarchy_Binary" = C_hierarchy_binary,
  "Binary_Binary" = C_binary_binary,
  "Binary_Skill" = C_binary_skill,
  "Binary_Hierarchy" = C_binary_hierarchy
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
      list(a, b, C, epsilon),
      \(a, b, C, epsilon){
      stopifnot(
        identical(names(a), rownames(C)),
        identical(names(b), colnames(C))
      )
      sinkhorn_log(a, b, C, epsilon)$plan
      }
      ),
    P_ind = map2(a, b, ~ .x %o% .y),
    score_dist=map(a, subset_score_dist)
  )|>
  select(dgp, a, b, epsilon, C, `Cost Matrix`=cost, P_obs, P_hat, P_ind, score_dist)

write_rds(sub_regime_model_fits, here("out", "sub_regime_model_fits.rds"))

