library(tidyverse)
library(here)
library(bcgovpond)
library(janitor)
library(conflicted)
conflicts_prefer(dplyr::filter)
#functions------------------------------------
source(here("R", "other.R"))
#constants----------------------------------------------
skills <- list()
hier <- list()

#skill distance stuff------------------------------------------
skills$onet_2019_to_soc_2018 <- read_view("2019_to_SOC_Crosswalk.xlsx", skip=3)|>
  select(o_net_soc_code=`O*NET-SOC 2019 Code`, soc_2018=`2018 SOC Code`)

skills$soc_2018_to_noc_2016 <- read_view("noc2016v1_3-soc2018us-eng.csv")|>
  select(soc_2018=`SOC 2018 (US) Code`, noc_2016=`NOC 2016  Version 1.3 Code`)

skills$noc_2016_to_noc_2021 <- read_view("noc2016v1_3-noc2021v1_0-eng.csv")|>
  select(noc_2016=`NOC 2016 V1.3 Code`, noc_2021=`NOC 2021 V1.0 Code`, noc2021_title=`NOC 2021 V1.0 Title`)|>
  mutate(noc_2016=str_pad(noc_2016, width=4, pad="0"),
         noc_2021=str_pad(noc_2021, width=5, pad="0"),
         noc2021_title=if_else(noc_2021 %in% c("00011", "00012", "00013", "00014", "00015"), "Senior managers - public and private sector", noc2021_title),
         noc_2021=if_else(noc_2021 %in% c("00011", "00012", "00013", "00014", "00015"), "00018", noc_2021)
        )

skills$mapping <- left_join(skills$onet_2019_to_soc_2018, skills$soc_2018_to_noc_2016)|>
  left_join(skills$noc_2016_to_noc_2021)%>%
  select(noc_2021, noc2021_title, o_net_soc_code)|>
  arrange(o_net_soc_code, noc_2021)|>
  distinct()

skills$nocs_we_want <- skills$mapping |>
  select(noc_2021, noc2021_title) |>
  distinct()|>
  mutate(noc_plus_title=paste(noc_2021, noc2021_title, sep=": "))|>
  select(noc_2021, noc_plus_title)|>
  filter(noc_2021!="44200")

skills$onet_raw <- tibble(file=c("Skills.xlsx", "Abilities.xlsx", "Knowledge.xlsx", "Work Activities.xlsx"))%>%
  mutate(data=map(file, read_data))%>%
  select(-file)%>%
  unnest(data)%>%
  pivot_wider(id_cols = o_net_soc_code, names_from = element_name, values_from = score)%>%
  inner_join(skills$mapping)%>%
  ungroup()%>%
  select(-o_net_soc_code, -noc2021_title)%>%
  select(noc_2021, everything())

skills$onet_mapped <- skills$onet_raw|>
  group_by(noc_2021)%>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

skills$four_digit <- skills$onet_raw|>
  mutate(noc_four=str_sub(noc_2021,1, 4))|>
  group_by(noc_four)|>
  summarise(across(contains(":"), ~mean(.x, na.rm = TRUE)))

skills$missing_four <- anti_join(skills$nocs_we_want, skills$onet_raw|>select(noc_2021))|>
  mutate(noc_four=str_sub(noc_2021, 1, 4),
         .after=noc_2021)|>
  inner_join(skills$four_digit)|>
  select(-noc_four, -noc_plus_title)

skills$onet_full <- bind_rows(skills$onet_mapped, skills$missing_four)|>
  inner_join(skills$nocs_we_want)|>
  ungroup()|>
  arrange(noc_2021)|>
  select(-noc_2021)|>
  column_to_rownames("noc_plus_title")

skills$onet_pca <- prcomp(skills$onet_full, center=TRUE, scale=TRUE)
skills$noc_coords <- skills$onet_pca$x[, 1:10]#keep first 10 components
skills$skills_noc_dist<- dist(skills$noc_coords, method = "euclidean")|>
  as.matrix()
skills$mean_dist <- skills$skills_noc_dist|>
  mean()

skills$mds2 <- cmdscale(skills$skills_noc_dist, k = 2)|>
  as.data.frame()|>
  rownames_to_column("noc_2021")|>
  mutate(noc_2021=str_sub(noc_2021, 1, 5))


#hierarchical distance--------------------------------------------
hier$hierarchy_distance <- crossing(origin=skills$nocs_we_want$noc_plus_title,
                               destination=skills$nocs_we_want$noc_plus_title)|>
  h_dist()|>
  select(origin, destination, distance)

hier$hier_counts <- hier$hierarchy_distance|>
  group_by(distance)|>
  count()|>
  mutate(distance=factor(distance))

hier$hier_wide <- hier$hierarchy_distance|>
  pivot_wider(names_from = destination, values_from = distance)|>
  column_to_rownames("origin")

hier$hier_mat <- hier$hier_wide|>as.matrix()
#sanity checks
stopifnot(hier$hier_mat==t(hier$hier_mat))
stopifnot(all(diag(hier$hier_mat) == 0))
stopifnot(min(hier$hier_mat) == 0, max(hier$hier_mat) == 9)

#destination gating

cip_noc <- bcgovpond::read_view("9810040401.csv") |>
  select(highest = starts_with("Highest"),
         cip     = starts_with("Major"),
         noc     = starts_with("Occupation"),
         value   = VALUE) |>
  filter(highest != "Total - Highest certificate, diploma or degree")

#education proportions (aggregation across nocs)
dest_p0 <- cip_noc |>
  summarise(total = sum(value), .by = c(highest, cip)) |>
  mutate(p0 = total / sum(total)) |>
  select(cip, highest, p0) |>
  unite(education, highest, cip, sep=": ")|>
  tibble::deframe()

noc_specificity <- cip_noc |>
  filter(value > 0) |> #no zeros anyways, but just incase
  add_count(noc, wt = value, name = "T")|> # noc total across educations
  mutate(
    p  = value / T, #probability of education conditional on noc
    education = paste(highest, cip, sep = ": "), #create variable education out of level and cip
    dest_p0 = dest_p0[education], # the unconditional probability of education
    kl_contribution = p * log(p / dest_p0) #the education's contribution to the distinctness of this occupation's education mix.
  )|>
  summarise(
    KL = sum(kl_contribution),
    T  = first(T),
    .by = noc
  )

dest_fit_kl <- lm(log(KL) ~ log(T), data = noc_specificity)

noc_specificity <- noc_specificity |>
  mutate(
    specificity = log(KL) - predict(dest_fit_kl),
    TEER=str_sub(noc, 2,2))

#education specificity

spec_p0 <- cip_noc |>
  summarise(total = sum(value), .by = noc) |>
  mutate(p0 = total / sum(total)) |>
  select(noc, p0) |>
  tibble::deframe()

educ_specificity <- cip_noc |>
  filter(value > 0) |>
  add_count(highest, cip, wt = value, name = "T") |>
  mutate(
    p  = value / T,
    spec_p0 = spec_p0[noc],
    kl_part = p * log(p / spec_p0)
  ) |>
  summarise(
    KL = sum(kl_part),
    T  = first(T),
    .by = c(highest, cip)
  )

spec_fit_kl <- lm(log(KL) ~ log(T), data = educ_specificity)

educ_specificity <- educ_specificity |>
  mutate(
    specificity = log(KL) - predict(spec_fit_kl)
  ) |>
  mutate(attain_bin = case_when(
    highest %in% c(
      "Postsecondary certificate, diploma or degree",
      "Postsecondary certificate or diploma below bachelor level",
      "Apprenticeship or trades certificate or diploma",
      "Non-apprenticeship trades certificate or diploma",
      "Apprenticeship certificate",
      "College, CEGEP or other non-university certificate or diploma",
      "University certificate or diploma below bachelor level"
    ) ~ "Non-university postsecondary",

    highest %in% c(
      "Bachelor's degree",
      "Bachelor’s degree or higher"
    ) ~ "Bachelor’s",

    highest == "University certificate or diploma above bachelor level" ~
      "Post-bachelor certificate",

    highest == "Master's degree" ~ "Master’s",

    highest %in% c(
      "Degree in medicine, dentistry, veterinary medicine or optometry",
      "Earned doctorate"
    ) ~ "Professional / Doctorate"
  ),
  attain_bin = factor(
    attain_bin,
    levels = c(
      "Non-university postsecondary",
      "Bachelor’s",
      "Post-bachelor certificate",
      "Master’s",
      "Professional / Doctorate"
    ),
    ordered = TRUE
  ))

#anchoring epsilon-------------------------

max_h <- max(hier$hier_mat, na.rm = TRUE)

h_all <- offdiag(hier$hier_mat)
s_all <- offdiag(skills$skills_noc_dist)

# Define the "non-maximal" hierarchical subset (exclude max distance and optionally zeros)
h_nonmax <- h_all[h_all < max_h & h_all > 0]

# Conditional quantile anchors within the informative (non-maximal) region
q_cond <- c(0.25, 0.50, 0.75)
h_anchor <- as.numeric(quantile(h_nonmax, probs = q_cond, na.rm = TRUE, type = 7))

# Map each anchor value to its unconditional percentile in the FULL hierarchical distribution
# (including max-distance pairs). This is the empirical CDF value F_H(h_anchor).
p_uncond <- sapply(h_anchor, function(x) mean(h_all <= x, na.rm = TRUE))

# Apply those percentiles to the skill distance distribution
s_anchor <- as.numeric(quantile(s_all, probs = p_uncond, na.rm = TRUE, type = 7))

# Summarize
calibration <- tibble::tibble(
  q_cond = q_cond,
  h_anchor = h_anchor,
  p_uncond = p_uncond,
  s_anchor = s_anchor
)

#write objects to disk--------------------------------
write_rds(skills, here("out", "skills.rds"))
write_rds(hier, here("out", "hier.rds"))
write_rds(noc_specificity, here("out", "noc_specificity.rds"))
write_rds(educ_specificity, here("out", "educ_specificity.rds"))
write_rds(calibration, here("out", "calibration.rds"))





