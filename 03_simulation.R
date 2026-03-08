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
  inner_join(skills$nocs_we_want, by=c("noc_5"="noc_2021"))|> #military with no skill data filtered out.
  select(-noc_5)

a <- lfs_data|>
  filter(syear<2015,
         age10!="55-64")|> #retired in end year
  mutate(age_broad=if_else(age10=="15-24", "young", "old"))|>
  summarize(count=sum(count), .by = c(age_broad, noc_plus_title))|>
  group_by(age_broad)|>
  mutate(prop=count/sum(count))|>
  select(age_broad,
         noc_plus_title,
         a=prop)|>
  nest()|>
  mutate(a=map(data, deframe))|>
  select(-data)

b <- lfs_data|>
  filter(syear>2020,
         age10!="15-24")|> #babies in start year
  mutate(age_broad=if_else(age10=="25-34", "young", "old"))|>
  summarize(count=sum(count), .by = c(age_broad, noc_plus_title))|>
  group_by(age_broad)|>
  mutate(prop=count/sum(count))|>
  select(age_broad,
         noc_plus_title,
         b=prop)|>
  nest()|>
  mutate(b=map(data, deframe))|>
  select(-data)

a_old <- a$a[a$age_broad=="old"][[1]]
b_old <- b$b[b$age_broad=="old"][[1]]
a_young <- a$a[a$age_broad == "young"][[1]]
b_young <- b$b[b$age_broad == "young"][[1]]

C_skill <- skills$skills_noc_dist/calibration$s_anchor[calibration$q_cond==.5]
C_hier <- hier$hier_mat/calibration$h_anchor[calibration$q_cond==.5]

simulations <- list(
    `DGP: 100% Skill` = list(
    a = a_old,
    b = b_old,
    C = C_skill
  ),
  `DGP: 75% Skill / 25% Hierarchy` = list(
    a = a_old,
    b = b_old,
    C = .75*C_skill+.25*C_hier
  ),
  `DGP: 50% Skill / 50% Hierarchy` = list(
    a = a_old,
    b = b_old,
    C = .5*C_skill+.5*C_hier
  ),
  `DGP: 25% Skill / 75% Hierarchy` = list(
    a = a_old,
    b = b_old,
    C = .25*C_skill+.75*C_hier
  ),
  `DGP: 100% Hierarchy` = list(
    a = a_old,
    b = b_old,
    C = C_hier
    )
  )

P_fake <- map(simulations, \(sim)
              sinkhorn_aligned(sim$a, sim$b, sim$C, 1, sinkhorn_log)$plan
)
P_fake$`DGP: Independence` <- outer_named(a_old, b_old)
P_fake$`DGP: Independence` <- P_fake$`DGP: Independence` / sum(P_fake$`DGP: Independence`)

#have the fake data, ready to simulate!

fit_grid <- crossing(
  dgp = names(P_fake),
  `Cost Matrix` = c("Skill","Hierarchy"),
  Temperature = 2^seq(-1,1,.1)[-11]
)

fit_grid <- fit_grid |>
  mutate(
    P_obs = map(dgp, ~P_fake[[.x]]),
    C = case_when(
      `Cost Matrix` == "Skill" ~ list(C_skill),
      `Cost Matrix` == "Hierarchy"  ~ list(C_hier)
    )
  )

fit_grid <- fit_grid |>
  mutate(
    a = map(P_obs, rowSums),
    b = map(P_obs, colSums),

    P_hat = pmap(
      list(a, b, C, Temperature),
      \(a, b, C, temp)
      sinkhorn_aligned(a, b, C, temp, sinkhorn_log)$plan
    )
  )

fit_grid <- fit_grid |>
  mutate(
    `KL-Divergence` = map2_dbl(P_obs, P_hat, kl_score)
  )

results <- fit_grid |>
  select(dgp, `Cost Matrix`, Temperature, `KL-Divergence`)|>
  mutate(
    dgp = factor(
      dgp,
      levels = c(
       "DGP: 100% Skill",
        "DGP: Independence",
        "DGP: 100% Hierarchy",
       "DGP: 75% Skill / 25% Hierarchy",
       "DGP: 50% Skill / 50% Hierarchy",
       "DGP: 25% Skill / 75% Hierarchy"
      )
    )
  )

plots$kl_plots <- ggplot(results, aes(Temperature, `KL-Divergence`, colour=`Cost Matrix`))+
  geom_vline(xintercept = 1, colour="white",lwd=2)+
  geom_line()+
  geom_point()+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log2")+
  facet_wrap(~dgp, nrow=2)+
  labs(title="Sinkhorn Performance across DGPs by Cost Matrix and Temperature",
       subtitle="KL-divergence is the log-likelihood loss: lower values indicate better fit")
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
       caption="The maximal hierarchical distance category pools a large number of heterogeneous occupation pairs and is omitted for visual clarity.")+
  theme_minimal()

#distance plot

hierskill <- distance_plot(P_fake[["DGP: 100% Hierarchy"]], C_skill, subtitle="DGP: Hierarchy | Cost: Skill (mis-match)")+
  theme(axis.title = element_blank())
hierhier <- distance_plot(P_fake[["DGP: 100% Hierarchy"]], C_hier, subtitle="DGP: Hierarchy | Cost: Hierarchy (matched)")+
  labs(
    x = "Occupational distance (normalized)",
    y = "Log excess transition probability (relative to independence)"
  )
skillhier <- distance_plot(P_fake[["DGP: 100% Skill"]], C_hier, subtitle="DGP: Skill | Cost: Hierarchy (mis-match)")+
  theme(axis.title = element_blank())
skillskill <- distance_plot(P_fake[["DGP: 100% Skill"]], C_skill, subtitle="DGP: Skill | Cost: Skill (matched)")+
  theme(axis.title = element_blank())

plots$linearity <-((skillskill + skillhier) /(hierhier + hierskill))

# by destination gating------------------------------------

noc_specificity <- read_rds(here("out", "noc_specificity.rds"))

skill_long <- C_skill|>
  as.data.frame()|>
  rownames_to_column("origin")|>
  pivot_longer(-origin, names_to = "destination", values_to = "distance")

hier_long <- C_hier|>
  as.data.frame()|>
  rownames_to_column("origin")|>
  pivot_longer(-origin, names_to = "destination", values_to = "distance")


base_quartiles <- noc_specificity |>
  select(noc, gating_quartile, hier_weight) |>
  group_by(gating_quartile, hier_weight) |>
  nest() |>
  mutate(skill = map(data, \(d)
                skill_long |>
                  semi_join(d, by = c("destination" = "noc")) |>
                  pivot_wider(names_from = destination, values_from = distance) |>
                  column_to_rownames("origin") |>
                  as.matrix()),
         hier = map(data, \(d)
               hier_long |>
                 semi_join(d, by = c("destination" = "noc")) |>
                 pivot_wider(names_from = destination, values_from = distance) |>
                 column_to_rownames("origin") |>
                 as.matrix()),
         C = pmap(list(skill, hier, hier_weight), weighted_cost),
         a = list(a_old),
         b = map(data, \(d) b_old[d$noc]),
         skill_fit = pmap(list(a, b, skill), \(a, b, C) sinkhorn_aligned(a, b, C, epsilon = 1, solver = sinkhorn_log)$plan),
         hier_fit = pmap(list(a, b, hier), \(a, b, C) sinkhorn_aligned(a, b, C, epsilon = 1, solver = sinkhorn_log)$plan)
         )

by_quartiles <- base_quartiles |>
  crossing(tibble(Temperature = 2^seq(-1,1,.1)[-11])) |>
  mutate(simulated = pmap(list(a, b, C, Temperature),\(a, b, C, T) sinkhorn_aligned(a, b, C, epsilon = T, solver = sinkhorn_log)$plan),
         kl_skill = map2_dbl(simulated, skill_fit, kl_score),
         kl_hier  = map2_dbl(simulated, hier_fit,  kl_score))


kl_by_quartiles <- by_quartiles |>
  select(hier_weight, gating_quartile, Temperature, kl_skill, kl_hier)|>
  pivot_longer(cols= contains("kl"),  names_to = "Cost Matrix", values_to = "KL-Divergence")|>
  mutate(`Cost Matrix`= if_else(`Cost Matrix`=="kl_skill", "Skill", "Hierarchical"))


plots$dest_gate_kl_plot <- ggplot(kl_by_quartiles, aes(Temperature, `KL-Divergence`, colour=`Cost Matrix`))+
  geom_vline(xintercept = 1, colour="white",lwd=2)+
  geom_line()+
  geom_point()+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log2")+
  facet_wrap(~gating_quartile, nrow=2)+
  labs(title="Sinkhorn Performance across DGPs by Cost Matrix and Temperature",
       subtitle="KL-divergence is the log-likelihood loss: lower values indicate better fit")

write_rds(plots, here("out", "plots.rds"))









