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

C_skill <- skills$skills_noc_dist/calibration$s_anchor[calibration$q_cond==.5]
C_hier <- hier$hier_mat/calibration$h_anchor[calibration$q_cond==.5]
C_binary <- binary_mat

#sanity check before proceeding
# stopifnot(
#   identical(rownames(C_skill), rownames(C_hier)),
#   identical(rownames(C_skill), rownames(C_binary)),
#   identical(colnames(C_skill), colnames(C_hier)),
#   identical(colnames(C_skill), colnames(C_binary))
# )
#
# global_simulations <- list(
#   `DGP: 75% Skill / 25% Hierarchy` = list(
#     a = a,
#     b = b,
#     C = .75*C_skill+.25*C_hier
#   ),
#   `DGP: 50% Skill / 50% Hierarchy` = list(
#     a = a,
#     b = b,
#     C = .5*C_skill+.5*C_hier
#   ),
#   `DGP: 25% Skill / 75% Hierarchy` = list(
#     a = a,
#     b = b,
#     C = .25*C_skill+.75*C_hier
#     )
#   )
#
# global_P_fake <- map(global_simulations, \(sim)
#               sinkhorn_aligned(sim$a, sim$b, sim$C, 1, sinkhorn_log)$plan
# )

#have the fake data, ready to simulate!

# global_sim <- crossing(
#   dgp = names(global_P_fake),
#   `Cost Matrix` = c("Skill","Hierarchy", "Binary"),
#   Temperature = 2^seq(-1,1,.1))|>
#   mutate(
#     P_obs = map(dgp, ~global_P_fake[[.x]]),
#     C = case_when(
#       `Cost Matrix` == "Skill" ~ list(C_skill),
#       `Cost Matrix` == "Hierarchy"  ~ list(C_hier),
#       `Cost Matrix` == "Binary" ~ list(C_binary)),
#     a = map(P_obs, rowSums),
#     b = map(P_obs, colSums),
#     P_hat = pmap(
#       list(a, b, C, Temperature),
#       \(a, b, C, temp)
#       sinkhorn_aligned(a, b, C, temp, sinkhorn_log)$plan),
#     P_ind = map2(a, b, ~ .x %o% .y),
#     `KL-hat` = map2_dbl(P_obs, P_hat, kl_score),
#     `KL-ind` = map2_dbl(P_obs, P_ind, kl_score),
#     `KL-rel` = (`KL-ind`-`KL-hat`)/`KL-ind`
#   )
#
# global_results <- global_sim |>
#   select(dgp, `Cost Matrix`, Temperature, `KL-hat`, `KL-ind`,`KL-rel`)|>
#   mutate(
#     dgp = factor(
#       dgp,
#       levels = c(
#        "DGP: 75% Skill / 25% Hierarchy",
#        "DGP: 50% Skill / 50% Hierarchy",
#        "DGP: 25% Skill / 75% Hierarchy"
#       )
#     )
#   )
#
# plots$kl_plots <- ggplot(global_results, aes(Temperature, `KL-rel`, colour=`Cost Matrix`))+
#   geom_vline(xintercept = 1, colour="white",lwd=2)+
#   geom_hline(yintercept = 0, colour="white",lwd=2)+
#   geom_line()+
#   geom_point()+
#   scale_x_continuous(trans="log2")+
#   facet_wrap(~dgp, nrow=1)+
#   labs(title="Model Performance across DGPs by Cost Matrix and Temperature",
#        subtitle="Relative KL improvement: 1 is perfect fit, 0 no improvement over independence, negative worse than independence",
#        y="relative KL improvement")
#


# market sub_regimes------------------------------------

noc_specificity <- read_rds(here("out", "noc_specificity.rds"))|>
  arrange(noc)|>
  select(noc_plus_title=noc, sub_regime)

sub_regime_vec <- setNames(noc_specificity$sub_regime, noc_specificity$noc_plus_title)
master_ids <- names(sub_regime_vec)

stopifnot(sum(sub_regime_vec == "Horizontal (Skill)") > 0) #no empty regimes
stopifnot(sum(sub_regime_vec == "Vertical (Hierarchy)") > 0) #no empty regimes
stopifnot(sum(sub_regime_vec == "Minimal (Binary)") > 0) #no empty regimes
stopifnot(!anyDuplicated(master_ids))
stopifnot(!anyNA(master_ids))
stopifnot(all(master_ids %in% names(a)))
stopifnot(all(master_ids %in% names(b)))
stopifnot(all(master_ids %in% rownames(C_skill)))
stopifnot(all(master_ids %in% colnames(C_skill)))
stopifnot(all(master_ids %in% rownames(C_hier)))
stopifnot(all(master_ids %in% colnames(C_hier)))
stopifnot(all(master_ids %in% rownames(C_binary)))
stopifnot(all(master_ids %in% colnames(C_binary)))

a_ordered <- a[master_ids]
b_ordered <- b[master_ids]
C_skill_ordered <- C_skill[master_ids, master_ids, drop = FALSE]
C_hier_ordered <- C_hier[master_ids, master_ids, drop = FALSE]
C_binary_ordered <- C_binary[master_ids, master_ids, drop = FALSE]

stopifnot(identical(names(a_ordered), master_ids))
stopifnot(identical(names(b_ordered), master_ids))
stopifnot(identical(rownames(C_skill_ordered), master_ids))
stopifnot(identical(colnames(C_skill_ordered), master_ids))
stopifnot(identical(rownames(C_hier_ordered), master_ids))
stopifnot(identical(colnames(C_hier_ordered), master_ids))
stopifnot(identical(rownames(C_binary_ordered), master_ids))
stopifnot(identical(colnames(C_binary_ordered), master_ids))

#a, b, and C's all share ordering (as defined by master_ids)

#subset the a vector, and keep appropriate rows of C matrices

a_skill_sub_regime <- a_ordered[sub_regime_vec == "Horizontal (Skill)"]
a_hier_sub_regime <- a_ordered[sub_regime_vec == "Vertical (Hierarchy)"]
a_binary_sub_regime <- a_ordered[sub_regime_vec == "Minimal (Binary)"]
C_skill_skill <- C_skill_ordered[sub_regime_vec == "Horizontal (Skill)", , drop = FALSE]
C_hier_hier <- C_hier_ordered[sub_regime_vec == "Vertical (Hierarchy)", , drop = FALSE]
C_binary_binary <- C_binary_ordered[sub_regime_vec == "Minimal (Binary)", , drop = FALSE]

sub_regime_simulations <- list(`Minimal (Binary)` = list( a = a_binary_sub_regime, b = b_ordered, C = C_binary_binary),
                                `Vertical (Hierarchy)` = list( a = a_hier_sub_regime, b = b_ordered, C = C_hier_hier),
                                `Horizontal (Skill)` = list( a = a_skill_sub_regime, b = b_ordered, C = C_skill_skill))

sub_regime_P_fake <- map(sub_regime_simulations, \(sim) sinkhorn_aligned(sim$a, sim$b, sim$C, 1, sinkhorn_log)$plan)






# plots$dest_gate_kl_plot <- ggplot(sub_regime_results, aes(Temperature, `KL-rel`, colour=cost_name))+
#   geom_vline(xintercept = 1, colour="white",lwd=2)+
#   geom_hline(yintercept = 0, colour="white",lwd=2)+
#   geom_hline(yintercept = 1, colour="white",lwd=2)+
#   geom_line()+
#   geom_point()+
#   scale_y_continuous(trans = scales::pseudo_log_trans(sigma = .1, base = 10), breaks= c(1,0,-50,-100,-150))+
#   scale_x_continuous(trans="log2")+
#   facet_wrap(~dgp_name, nrow=1)+
#   labs(title="Model Performance across DGPs by Cost Matrix and Temperature",
#        subtitle="Relative KL improvement: 1 is perfect fit, 0 no improvement over independence, negative worse than independence",
#        y="relative KL improvement")


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

write_rds(plots, here("out", "plots.rds"))









