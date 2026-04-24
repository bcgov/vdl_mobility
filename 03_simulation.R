library(tidyverse)
library(here)
library(janitor)
library(vroom)
library(matrixStats)
library(plotly)
library(bcgovpond)
library(safesink)
library(conflicted)
conflicts_prefer(vroom::cols)
conflicts_prefer(vroom::col_double)
conflicts_prefer(vroom::col_character)
conflicts_prefer(dplyr::filter)

source(here("R", "other.R"))
sim_plots <- list()
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


#sanity check before proceeding
stopifnot(
  identical(rownames(skill_dist), rownames(hier_dist)),
  identical(rownames(skill_dist), rownames(binary_dist)),
  identical(colnames(skill_dist), colnames(hier_dist)),
  identical(colnames(skill_dist), colnames(binary_dist))
)

global_simulations <- list(
  `DGP: 75% Skill / 25% Hierarchy` = list(
    a = a,
    b = b,
    C = .75*skill_dist+.25*hier_dist
  ),
  `DGP: 50% Skill / 50% Hierarchy` = list(
    a = a,
    b = b,
    C = .5*skill_dist+.5*hier_dist
  ),
  `DGP: 25% Skill / 75% Hierarchy` = list(
    a = a,
    b = b,
    C = .25*skill_dist+.75*hier_dist
    )
  )

global_P_fake <- map(global_simulations, \(sim)
              sinkhorn_aligned(sim$a, sim$b, sim$C, 1, sinkhorn_log)$plan
)

#have the fake data, ready to simulate!

global_sim <- crossing(
  dgp = names(global_P_fake),
  `Cost Matrix` = c("Skill","Hierarchy", "Binary"),
  Temperature = 2^seq(-1,1,.1))|>
  mutate(
    P_obs = map(dgp, ~global_P_fake[[.x]]),
    C = case_when(
      `Cost Matrix` == "Skill" ~ list(skill_dist),
      `Cost Matrix` == "Hierarchy"  ~ list(hier_dist),
      `Cost Matrix` == "Binary" ~ list(binary_dist)),
    a = map(P_obs, rowSums),
    b = map(P_obs, colSums),
    P_hat = pmap(
      list(a, b, C, Temperature),
      \(a, b, C, temp)
      sinkhorn_aligned(a, b, C, temp, sinkhorn_log)$plan),
    P_ind = map2(a, b, ~ .x %o% .y),
    `KL-hat` = map2_dbl(P_obs, P_hat, kl_score),
    `KL-ind` = map2_dbl(P_obs, P_ind, kl_score),
    `KL-rel` = (`KL-ind`-`KL-hat`)/`KL-ind`
  )

global_results <- global_sim |>
  select(dgp, `Cost Matrix`, Temperature, `KL-hat`, `KL-ind`,`KL-rel`)|>
  mutate(
    dgp = factor(
      dgp,
      levels = c(
       "DGP: 75% Skill / 25% Hierarchy",
       "DGP: 50% Skill / 50% Hierarchy",
       "DGP: 25% Skill / 75% Hierarchy"
      )
    )
  )

sim_plots$kl_plots <- ggplot(global_results, aes(Temperature, `KL-rel`, colour=`Cost Matrix`))+
  geom_vline(xintercept = 1, colour="white",lwd=2)+
  geom_hline(yintercept = 0, colour="white",lwd=2)+
  geom_line()+
  geom_point()+
  scale_x_continuous(trans="log2")+
  facet_wrap(~dgp, nrow=1)+
  labs(title="Model Performance across DGPs by Cost Matrix and Temperature",
       subtitle="Relative KL improvement: 1 is perfect fit, 0 no improvement over independence, negative worse than independence",
       y="relative KL improvement")



# market sub_regimes------------------------------------

noc_specificity <- read_rds(here("out", "noc_specificity.rds"))|>
  arrange(noc)|>
  select(noc_plus_title=noc, sub_regime)

sub_regime_vec <- setNames(noc_specificity$sub_regime, noc_specificity$noc_plus_title)
master_ids <- names(sub_regime_vec)
sub_regime_vec <- sub_regime_vec[master_ids]
stopifnot(!anyNA(sub_regime_vec))
stopifnot(sum(sub_regime_vec == "Horizontal (Skill)") > 0) #no empty regimes
stopifnot(sum(sub_regime_vec == "Vertical (Hierarchy)") > 0) #no empty regimes
stopifnot(sum(sub_regime_vec == "Minimal (Binary)") > 0) #no empty regimes
stopifnot(!anyDuplicated(master_ids))
stopifnot(!anyNA(master_ids))
stopifnot(all(master_ids %in% names(a)))
stopifnot(all(master_ids %in% names(b)))
stopifnot(all(master_ids %in% rownames(skill_dist)))
stopifnot(all(master_ids %in% colnames(skill_dist)))
stopifnot(all(master_ids %in% rownames(hier_dist)))
stopifnot(all(master_ids %in% colnames(hier_dist)))
stopifnot(all(master_ids %in% rownames(binary_dist)))
stopifnot(all(master_ids %in% colnames(binary_dist)))

a_ordered <- a[master_ids]
b_ordered <- b[master_ids]
skill_dist_ordered <- skill_dist[master_ids, master_ids, drop = FALSE]
hier_dist_ordered <- hier_dist[master_ids, master_ids, drop = FALSE]
binary_dist_ordered <- binary_dist[master_ids, master_ids, drop = FALSE]

stopifnot(identical(names(a_ordered), master_ids))
stopifnot(identical(names(b_ordered), master_ids))
stopifnot(identical(rownames(skill_dist_ordered), master_ids))
stopifnot(identical(colnames(skill_dist_ordered), master_ids))
stopifnot(identical(rownames(hier_dist_ordered), master_ids))
stopifnot(identical(colnames(hier_dist_ordered), master_ids))
stopifnot(identical(rownames(binary_dist_ordered), master_ids))
stopifnot(identical(colnames(binary_dist_ordered), master_ids))

#a, b, and C's all share ordering (as defined by master_ids)

#subset the a vector, and keep appropriate rows of C matrices

a_skill_sub_regime <- a_ordered[sub_regime_vec == "Horizontal (Skill)"]
a_hier_sub_regime <- a_ordered[sub_regime_vec == "Vertical (Hierarchy)"]
a_binary_sub_regime <- a_ordered[sub_regime_vec == "Minimal (Binary)"]

#Horizontal Regime
C_horizontal_skill <- skill_dist_ordered[sub_regime_vec == "Horizontal (Skill)", , drop = FALSE]
C_horizontal_hier <- hier_dist_ordered[sub_regime_vec == "Horizontal (Skill)", , drop = FALSE]
C_horizontal_binary <- binary_dist_ordered[sub_regime_vec == "Horizontal (Skill)", , drop = FALSE]

#Vertical Regime
C_vertical_hier <- hier_dist_ordered[sub_regime_vec == "Vertical (Hierarchy)", , drop = FALSE]
C_vertical_skill <- skill_dist_ordered[sub_regime_vec == "Vertical (Hierarchy)", , drop = FALSE]
C_vertical_binary <- binary_dist_ordered[sub_regime_vec == "Vertical (Hierarchy)", , drop = FALSE]

#Minimal Regime
C_minimal_binary <- binary_dist_ordered[sub_regime_vec == "Minimal (Binary)", , drop = FALSE]
C_minimal_skill <- skill_dist_ordered[sub_regime_vec == "Minimal (Binary)", , drop = FALSE]
C_minimal_hier <- hier_dist_ordered[sub_regime_vec == "Minimal (Binary)", , drop = FALSE]


sub_regime_simulations <- list(`Minimal` = list( a = a_binary_sub_regime, b = b_ordered, C = C_minimal_binary),
                                `Vertical` = list( a = a_hier_sub_regime, b = b_ordered, C = C_vertical_hier),
                                `Horizontal` = list( a = a_skill_sub_regime, b = b_ordered, C = C_horizontal_skill))

sub_regime_P_fake <- map(sub_regime_simulations, \(sim) sinkhorn_aligned(sim$a, sim$b, sim$C, 1, sinkhorn_log)$plan)

#put together the simulation

REGIMES <- c("Minimal", "Vertical", "Horizontal")
COSTS <- c("Binary", "Hier", "Skill")
TEMPS <- 2^seq(-1, 1, .1)

sub_regime_grid <- tidyr::crossing(
  regime = REGIMES,
  cost = COSTS,
  Temperature = TEMPS)|>
  mutate(which_cost=paste(regime,cost,sep="_"))

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
    a = map(regime, ~ a_list[[.x]]),
    b = map(regime, ~ b_ordered),
    C = map(which_cost, ~ C_list[[.x]]),
    P_obs =map(regime, ~sub_regime_P_fake[[.x]])
  )

stopifnot(
  all(map2_lgl(sub_regime_grid$a, sub_regime_grid$C, ~ length(.x) == nrow(.y))),
  all(map2_lgl(sub_regime_grid$b, sub_regime_grid$C, ~ length(.x) == ncol(.y)))
)

sub_regime_grid <- sub_regime_grid |>
  mutate(
    P_hat = pmap(
      list(a, b, C, Temperature),
      \(a, b, C, temp)
      sinkhorn_aligned(a, b, C, temp, sinkhorn_log)$plan
    ),
    P_ind = map2(a, b, ~ .x %o% .y)
  )

sub_regime_grid <- sub_regime_grid |>
  mutate(
    KL_hat = map2_dbl(P_obs, P_hat, kl_score),
    KL_ind = map2_dbl(P_obs, P_ind, kl_score),
    KL_rel = (KL_ind - KL_hat) / KL_ind
  )

sub_regime_results <- sub_regime_grid|>
  select(regime, Temperature, cost, contains("KL"))

sim_plots$dest_gate_kl_plot <- ggplot(sub_regime_grid, aes(Temperature, `KL_rel`, colour=cost))+
  geom_vline(xintercept = 1, colour="white",lwd=2)+
  geom_hline(yintercept = 0, colour="white",lwd=2)+
  geom_hline(yintercept = 1, colour="white",lwd=2)+
  geom_line()+
  geom_point()+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = .1, base = 10), breaks= c(1,0,-50,-100,-150))+
  scale_x_continuous(trans="log2")+
  facet_wrap(~regime, nrow=1)+
  labs(title="Model Performance across DGPs by Cost Matrix and Temperature",
       subtitle="Relative KL improvement: 1 is perfect fit, 0 no improvement over independence, negative worse than independence",
       y="relative KL improvement")


write_rds(sim_plots, here("out", "sim_plots.rds"))


# diagnostic plots (should be linear relationship)

# hierskill <- distance_plot(P_fake[["DGP: 100% Hierarchy"]], C_skill, subtitle="DGP: Hierarchy | Cost: Skill (mis-match)")+
#   theme(axis.title = element_blank())
# hierhier <- distance_plot(P_fake[["DGP: 100% Hierarchy"]], C_hier, subtitle="DGP: Hierarchy | Cost: Hierarchy (matched)")+
#   labs(
#     x = "Occupational distance (normalized)",
#     y = "Log excess transition probability"
#   )
# skillhier <- distance_plot(P_fake[["DGP: 100% Skill"]], C_hier, subtitle="DGP: Skill | Cost: Hierarchy (mis-match)")+
#   theme(axis.title = element_blank())
# skillskill <- distance_plot(P_fake[["DGP: 100% Skill"]], C_skill, subtitle="DGP: Skill | Cost: Skill (matched)")+
#   theme(axis.title = element_blank())
#
# plots$linearity <-((skillskill + skillhier) /(hierhier + hierskill))











