 # library(reticulate)
 # reticulate::py_require("POT")
 # ot <- import("ot")
library(tidyverse)
library(here)
library(janitor)
library(vroom)
library(bcgovpond)
library(patchwork)
library(matrixStats)
library(conflicted)
conflicts_prefer(vroom::cols)
conflicts_prefer(vroom::col_double)
conflicts_prefer(vroom::col_character)
conflicts_prefer(dplyr::filter)

source(here("R","sinkhorn_utils.R"))
source(here("R", "other.R"))
lfs <- list()
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
  semi_join(skills$nocs_we_want, by=c("noc_5"="noc_2021")) #military with no skill data filtered out.


from <- lfs_data|>
  filter(syear<2015,
         age10!="55-64")|> #retired in end year
  mutate(age_broad=if_else(age10=="15-24", "young", "old"))|>
  summarize(count=sum(count), .by = c(age_broad, noc_5))|>
  group_by(age_broad)|>
  mutate(prop=count/sum(count))|>
  select(age_broad,
         noc_5,
         from_prop=prop)

to <- lfs_data|>
  filter(syear>2020,
         age10!="15-24")|> #babies in start year
  mutate(age_broad=if_else(age10=="25-34", "young", "old"))|>
  summarize(count=sum(count), .by = c(age_broad, noc_5))|>
  group_by(age_broad)|>
  mutate(prop=count/sum(count))|>
  select(age_broad,
         noc_5,
         to_prop=prop)

margins <- full_join(from, to, by = join_by(age_broad, noc_5))|>
  full_join(skills$nocs_we_want, by=c("noc_5"="noc_2021"))|>
  mutate(diff=to_prop-from_prop,
         abs_diff=abs(diff))

#visualizations

# young_change <- margins|>
#   filter(age_broad=="young")|>
#   slice_max(abs_diff, n=20)|>
#   ggplot(aes(diff, fct_reorder(noc_plus_title, diff)))+
#   geom_col()+
#   coord_cartesian(xlim = c(-.075, .025))+
#   labs(title="Occupational changes for young",
#        x="Change in Occupational proportion",
#        y=NULL
#        )
#
# old_change <- margins|>
#   filter(age_broad=="old")|>
#   slice_max(abs_diff, n=20)|>
#   ggplot(aes(diff, fct_reorder(noc_plus_title, diff)))+
#   geom_col()+
#   coord_cartesian(xlim = c(-.075, .025))+
#   labs(title="Occupational changes for old",
#        x="Change in Occupational proportion",
#        y=NULL
#   )
#
# young_change+old_change

C_skill <- skills$skills_noc_dist
C_hier <- hier[["hier_mat"]]

epsilon_skill <- calibration$s_anchor[calibration$q_cond==.5]
epsilon_hier <- calibration$h_anchor[calibration$q_cond==.5]

a_young <- extract_margin(margins, "young", from_prop)
b_young <- extract_margin(margins, "young", to_prop)
a_old <- extract_margin(margins, "old", from_prop)
b_old <- extract_margin(margins, "old", to_prop)

 # a_safe <- make_safe(a_young)
 # b_safe <- make_safe(b_young)

P <- sinkhorn_aligned(a_young, b_young, C_skill, epsilon=epsilon_skill, solver=sinkhorn_log)
#P_p <- sinkhorn_aligned(a_safe, b_safe, C_skill, epsilon=epsilon_skill, solver=ot$sinkhorn)

check_transport(P$plan, a_safe, b_safe)
 # check_transport(P_p, a_safe, b_safe)
 # compare_transport(P$plan, P_p)







