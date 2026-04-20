library(tidyverse)
library(here)
library(janitor)
library(vroom)
library(bcgovpond)
library(patchwork)
library(matrixStats)
library(plotly)
library(safesink)
library(conflicted)
conflicts_prefer(vroom::cols)
conflicts_prefer(vroom::col_double)
conflicts_prefer(vroom::col_character)
conflicts_prefer(dplyr::filter)

source(here("R", "other.R"))
plots <- list()
#read data-------------------
calibration <- read_rds(here("out","calibration.rds"))
skills <- read_rds(here("out","skills.rds"))
hier <- read_rds(here("out","hier.rds"))
binary_mat <- read_rds(here("out","binary_mat.rds"))
#get the noc list....

noc_list <- skills$mapped_nocs|>
  mutate(noc_5=str_sub(noc_plus_title, 1,5))

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
  with(setNames(prop, noc_plus_title))

b <- lfs_data|>
  filter(syear>2020,
         !age10 %in% c("15-24", "25-34"))|> #35-64 in end year
  summarize(count=sum(count), .by = c(noc_plus_title))|>
  mutate(prop=count/sum(count))|>
  select(noc_plus_title, prop)|>
  with(setNames(prop, noc_plus_title))

C_skill <- skills$skills_noc_dist/calibration$s_anchor[calibration$q_cond==.5]
C_hier <- hier$hier_mat/calibration$h_anchor[calibration$q_cond==.5]
C_binary <- binary_mat

global_simulations <- list(
  `DGP: 75% Skill / 25% Hierarchy` = list(
    a = a,
    b = b,
    C = .75*C_skill+.25*C_hier
  ),
  `DGP: 50% Skill / 50% Hierarchy` = list(
    a = a,
    b = b,
    C = .5*C_skill+.5*C_hier
  ),
  `DGP: 25% Skill / 75% Hierarchy` = list(
    a = a,
    b = b,
    C = .25*C_skill+.75*C_hier
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
      `Cost Matrix` == "Skill" ~ list(C_skill),
      `Cost Matrix` == "Hierarchy"  ~ list(C_hier),
      `Cost Matrix` == "Binary" ~ list(C_binary)),
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

# by destination gating------------------------------------

noc_specificity <- read_rds(here("out", "noc_specificity.rds"))|>
  select(noc, gating)

gating_vec <- setNames(noc_specificity$gating, noc_specificity$noc)

b_skill_tertile <- b[gating_vec=="Low"]
b_hier_tertile <- b[gating_vec == "Medium"]
b_binary_tertile <- b[gating_vec == "High"]

C_skill_tertile <- C_skill[, gating_vec == "Low", drop = FALSE]
C_hier_tertile <- C_hier[, gating_vec == "Medium", drop = FALSE]
C_binary_tertile <- C_binary[, gating_vec == "High", drop = FALSE]


tertile_simulations <- list(
  `High gating: DGP=Binary` = list(
    a = a,
    b = b_binary_tertile,
    C = C_binary_tertile
  ),
  `Medium gating: DGP=Hierarchy` = list(
    a = a,
    b = b_hier_tertile,
    C = C_hier_tertile
  ),
  `Low gating: DGP=Skill` = list(
    a = a,
    b = b_skill_tertile,
    C = C_skill_tertile
  )
)

tertile_P_fake <- map(tertile_simulations, \(sim)
                     sinkhorn_aligned(sim$a, sim$b, sim$C, 1, sinkhorn_log)$plan
)

#have the fake data, ready to simulate!

tertile_sim <- crossing(
  dgp = names(tertile_P_fake),
  `Cost Matrix` = c("Skill","Hierarchy", "Binary"),
  Temperature = 2^seq(-1,1,.1)) |>
  mutate(
    C = pmap(
      list(dgp, `Cost Matrix`),
      \(dgp, cm) {
        gate <- case_when(
          dgp == "High gating: DGP=Binary" ~ "High",
          dgp == "Medium gating: DGP=Hierarchy" ~ "Medium",
          dgp == "Low gating: DGP=Skill" ~ "Low"
        )
        C_mat <- switch(cm,
                        "Binary"    = C_binary,
                        "Hierarchy" = C_hier,
                        "Skill"     = C_skill
        )
        C_mat[, gating_vec == gate, drop = FALSE]
      }
    ),
    P_obs = map(dgp, ~tertile_P_fake[[.x]]),
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

tertile_results <- tertile_sim|>
  select(dgp, `Cost Matrix`, Temperature, `KL-hat`, `KL-ind`,`KL-rel`)|>
  mutate(
    dgp = factor(
      dgp,
      levels = c(
        "Low gating: DGP=Skill",
        "Medium gating: DGP=Hierarchy",
        "High gating: DGP=Binary"
      )
    )
  )








#plots---------------------------------------------

#relationship between cost metrics------------------------------------------
c_skill <- offdiag(C_skill)
c_hier  <- offdiag(C_hier)

c_both <- tibble(skill=c_skill, hier=c_hier)|>
  mutate(hier=ordered(hier))

distance_spearman <- round(cor(c_skill, c_hier, method="spearman"),3)

plots$dist_relationship <- c_both|>
  filter(hier<max(hier))|>
  ggplot(aes(hier, skill))+
  geom_jitter(size=.25, alpha=.25)+
  geom_boxplot(fill="red", alpha=.25, outliers=FALSE)+
  labs(title=paste("Hierarchical and skill distances: Spearman correlation", distance_spearman),
       x="Scaled hierarchical distance",
       y="Scaled skill distance",
       caption="The maximal hierarchical distance category is omitted for visual clarity.")+
  theme_minimal()

plots$kl_plots <- ggplot(global_results, aes(Temperature, `KL-rel`, colour=`Cost Matrix`))+
  geom_vline(xintercept = 1, colour="white",lwd=2)+
  geom_hline(yintercept = 0, colour="white",lwd=2)+
  geom_line()+
  geom_point()+
  scale_x_continuous(trans="log2")+
  facet_wrap(~dgp, nrow=1)+
  labs(title="Model Performance across DGPs by Cost Matrix and Temperature",
       subtitle="Relative KL improvement: 1 is perfect fit, 0 no improvement over independence, negative worse than independence",
       y="relative KL improvement")

plots$dest_gate_kl_plot <- ggplot(tertile_results, aes(Temperature, `KL-rel`, colour=`Cost Matrix`))+
  geom_vline(xintercept = 1, colour="white",lwd=2)+
  geom_hline(yintercept = 0, colour="white",lwd=2)+
  geom_hline(yintercept = 1, colour="white",lwd=2)+
  geom_line()+
  geom_point()+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = .1, base = 10), breaks= c(1,0,-50,-100,-150))+
  scale_x_continuous(trans="log2")+
  facet_wrap(~dgp, nrow=1)+
  labs(title="Model Performance across DGPs by Cost Matrix and Temperature",
       subtitle="Relative KL improvement: 1 is perfect fit, 0 no improvement over independence, negative worse than independence",
       y="relative KL improvement")

write_rds(plots, here("out", "plots.rds"))



#distance plot

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








