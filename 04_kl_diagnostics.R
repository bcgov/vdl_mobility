library(tidyverse)
library(here)
library(janitor)
library(safesink)
library(plotly)
#functions----------------------------

rowwise_kl_plot <- function(dgp, data){
  plt <- ggplot(data, aes(kl_ind, kl_hat, group = origin, text = origin)) +
    geom_point(aes(size = kl_ind, fill = kl_rel),
               shape = 21, colour = "black", stroke = 0.25) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_size_area(
      max_size = 8,
      limits = range(data$kl_ind, na.rm = TRUE)
    ) +
    scale_fill_viridis_c(
      limits = range(data$kl_rel, na.rm = TRUE)
    ) +
    labs(
      x = "KL(truth||independence)\n existing structure",
      y = "KL(truth||model)\nremaining structure",
      title = paste(dgp, "& true temperature=1 \n Model improvement over independence (rowwise KL) Bubble size ∝ row importance; below dashed line = improvement"),
      fill = "Proportion structure explained",
      size ="Row importance"
    ) +
    facet_grid(`Cost Matrix` ~ Temperature, labeller = labeller(Temperature = label_both))+
    theme(
      plot.margin = margin(t = 50, r = 40, b = 70, l = 70),
      axis.title.x = element_text(margin = margin(t = 40)),
      axis.title.y = element_text(margin = margin(r = 10)),
      text=element_text(size=15)
    )
  ggplotly(plt, tooltip = "text")
}

#initialize storage and load data-------------------------------------------

sim_kl_plots <- list()
global_model_fits <- read_rds(here("out", "global_model_fits.rds"))
sub_regime_model_fits <- read_rds(here("out", "sub_regime_model_fits.rds"))

#whole market, global KL-------------------------------------

global_KL_measures <- global_model_fits|>
  mutate(KL_global_hat = map2_dbl(P_obs, P_hat, kl_score),
         KL_global_ind = map2_dbl(P_obs, P_ind, kl_score),
         KL_global_rel = (KL_global_ind - KL_global_hat) / KL_global_ind,
         KL_row_hat = map2(P_obs, P_hat, ~ {
           tibble::enframe(rowwise_kl(.x, .y), name = "origin", value = "kl_hat")
         }),
         KL_row_ind = map2(P_obs, P_ind, ~ {
           tibble::enframe(rowwise_kl(.x, .y), name = "origin", value = "kl_ind")
         }),
         KL_row_rel = map2(KL_row_hat, KL_row_ind, ~ {
           dplyr::left_join(.x, .y, by = "origin") |>
             dplyr::mutate(
               kl_rel = dplyr::if_else(
                 kl_ind > 0,
                 (kl_ind - kl_hat) / kl_ind,
                 NA_real_
               )
             )
         }))


sim_kl_plots$global_kl_rel <- ggplot(global_KL_measures, aes(Temperature, `KL_global_rel`, colour=`Cost Matrix`))+
  geom_vline(xintercept = 1, colour="white",lwd=2)+
  geom_hline(yintercept = 0, colour="white",lwd=2)+
  geom_line()+
  geom_point()+
  scale_x_continuous(trans="log2")+
  facet_wrap(~dgp, nrow=1)+
  labs(title="Relative KL Improvement across DGPs by Cost Matrix and Temperature",
       subtitle="1 is perfect fit, 0 no improvement over independence, negative worse than independence",
       y="relative KL improvement")


sim_kl_plots$global_rowwise_kl <- global_KL_measures |>
  filter(Temperature %in% c(.5,1,2))|>
  select(dgp,
         `Cost Matrix`,
         Temperature,
         KL_row_rel)|>
  unnest(KL_row_rel)|>
  group_by(dgp)|>
  nest()|>
  mutate(plot=map2(dgp,data, rowwise_kl_plot))|>
  pull(plot)


#sub regimes--------------------------------------------

sub_regime_kl_measures <- sub_regime_model_fits|>
  mutate(KL_global_hat = map2_dbl(P_obs, P_hat, kl_score),
         KL_global_ind = map2_dbl(P_obs, P_ind, kl_score),
         KL_global_rel = (KL_global_ind - KL_global_hat) / KL_global_ind,
         KL_row_hat = map2(P_obs, P_hat, ~ {
           tibble::enframe(rowwise_kl(.x, .y), name = "origin", value = "kl_hat")
         }),
         KL_row_ind = map2(P_obs, P_ind, ~ {
           tibble::enframe(rowwise_kl(.x, .y), name = "origin", value = "kl_ind")
         }),
         KL_row_rel = map2(KL_row_hat, KL_row_ind, ~ {
           dplyr::left_join(.x, .y, by = "origin") |>
             dplyr::mutate(
               kl_rel = dplyr::if_else(
                 kl_ind > 0,
                 (kl_ind - kl_hat) / kl_ind,
                 NA_real_
               )
             )
         }))


sim_kl_plots$sub_regime_kl_rel <- ggplot(sub_regime_kl_measures, aes(Temperature, KL_global_rel, colour=`Cost Matrix`))+
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


#sub regime rowwise KL-----------------------------------------------

sim_kl_plots$sub_regime_rowwise_kl <- sub_regime_kl_measures |>
  filter(Temperature %in% c(.5,1,2))|>
  select(dgp,
         `Cost Matrix`,
         Temperature,
         KL_row_rel)|>
  unnest(KL_row_rel)|>
  group_by(dgp)|>
  nest()|>
  mutate(plot=map2(dgp,data, rowwise_kl_plot))|>
  pull(plot)




write_rds(sim_kl_plots, here("out", "sim_kl_plots.rds"))


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















