library(tidyverse)
library(here)
library(janitor)
library(safesink)
library(plotly)

#initialize storage and load data-------------------------------------------

sim_score_plots <- list()
global_model_fits <- read_rds(here("out", "global_model_fits.rds"))
sub_regime_model_fits <- read_rds(here("out", "sub_regime_model_fits.rds"))

#whole market-------------------------------------

global_measures <- global_model_fits|>
  mutate(KL_global_hat = map2_dbl(P_obs, P_hat, kl_score),
         KL_global_ind = map2_dbl(P_obs, P_ind, kl_score),
         `relative KL improvement` = (KL_global_ind - KL_global_hat) / KL_global_ind,
         CWTVD_global_hat_hier = pmap_dbl(list(P_obs, P_hat, hier_dist), cwtvd),
         CWTVD_global_ind_hier = pmap_dbl(list(P_obs, P_ind, hier_dist), cwtvd),
         `relative CWTVD improvement (Hierarchy)` = (CWTVD_global_ind_hier - CWTVD_global_hat_hier) / CWTVD_global_ind_hier,
         CWTVD_global_hat_skill = pmap_dbl(list(P_obs, P_hat, skill_dist), cwtvd),
         CWTVD_global_ind_skill = pmap_dbl(list(P_obs, P_ind, skill_dist), cwtvd),
         `relative CWTVD improvement (Skill)` = (CWTVD_global_ind_skill - CWTVD_global_hat_skill) / CWTVD_global_ind_skill
         )|>
  select(dgp, epsilon, `Cost Matrix`, contains("relative"))|>
  pivot_longer(cols=contains("relative"),names_to = "Measure")

sim_score_plots$global_measures <- global_measures|>
  ggplot(aes(epsilon, value, colour=`Cost Matrix`))+
  geom_vline(xintercept = 1, colour="white",lwd=2)+
  geom_hline(yintercept = 0, colour="white",lwd=2)+
  geom_hline(yintercept = 1, colour="white",lwd=2)+
  geom_line()+
  geom_point()+
  scale_x_continuous(trans="log10")+
  scale_y_continuous(transform = scales::pseudo_log_trans(sigma = .1, base = 2),
                     labels = scales::label_number())+
  facet_grid(Measure~dgp, scales="free_y")+
  labs(title="Improvement relative to independence by \u03b5 across DGPs and Measure",
       subtitle="1 is perfect fit, 0 no improvement over independence, negative worse than independence",
       x= "\u03b5",
       y="improvement relative to independence")

#sub regimes--------------------------------------------

sub_regime_measures <- sub_regime_model_fits|>
  mutate(KL_sub_hat = map2_dbl(P_obs, P_hat, kl_score),
         KL_sub_ind = map2_dbl(P_obs, P_ind, kl_score),
         `relative KL improvement` = (KL_sub_ind - KL_sub_hat) / KL_sub_ind,
         CWTVD_global_hat = pmap_dbl(list(P_obs, P_hat, score_dist), cwtvd),
         CWTVD_global_ind = pmap_dbl(list(P_obs, P_ind, score_dist), cwtvd),
         `relative CWTVD improvement` = (CWTVD_global_ind - CWTVD_global_hat) / CWTVD_global_ind)|>
  select(dgp, epsilon, `Cost Matrix`, contains("relative"))|>
  pivot_longer(cols=contains("relative"),names_to = "Measure")

sim_score_plots$sub_regime_measures <- sub_regime_measures|>
  ggplot(aes(epsilon, value, colour=`Cost Matrix`))+
  geom_vline(xintercept = 1, colour="white",lwd=2)+
  geom_hline(yintercept = 0, colour="white",lwd=2)+
  geom_hline(yintercept = 1, colour="white",lwd=2)+
  geom_line()+
  geom_point()+
  scale_x_continuous(trans="log10", labels = scales::label_comma())+
  scale_y_continuous(transform = scales::pseudo_log_trans(sigma = .1, base = 10),
                     labels = scales::label_number(),
                     breaks=c(-10,-1,0,1),
                     limits=c(-10,1))+
  facet_grid(Measure~dgp, labeller = labeller(dgp = label_both), scale="free_y")+
  labs(title="Improvement relative to independence by \u03b5 across DGPs and Measure",
       subtitle="1 is perfect fit, 0 no improvement over independence, negative worse than independence",
       x= "\u03b5",
       y="improvement relative to independence")


write_rds(sim_score_plots, here("out", "sim_score_plots.rds"))


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















