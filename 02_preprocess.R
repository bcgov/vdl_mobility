library(tidyverse)
library(here)
library(bcgovpond) #to install pak::pak("bcgov/bcgovpond")
library(janitor)
library(conflicted)
library(stringr)
conflicts_prefer(dplyr::count)
conflicts_prefer(dplyr::filter)
#functions------------------------------------
source(here("R", "other.R"))
#constants----------------------------------------------

skills <- list()
hier <- list()

#skill distance stuff------------------------------------------

mapping <- read_csv("https://raw.githubusercontent.com/bcgov/onet-noc2021-crosswalk/main/output/3.0.0/onet_to_noc2021_mapping.csv")

mapped_nocs <- mapping|>
  select(noc_plus_title)|>
  distinct()

onet_raw <- tibble(file=c("Skills.xlsx", "Abilities.xlsx", "Knowledge.xlsx", "Work Activities.xlsx"))%>%
  mutate(data=map(file, read_data))%>%
  select(-file)%>%
  unnest(data)%>%
  pivot_wider(id_cols = o_net_soc_code, names_from = element_name, values_from = score)

skill_cols <- setdiff(names(onet_raw), "o_net_soc_code")

noc_skills <- mapping|>
  ungroup()|>
  left_join(onet_raw, by = c("onet_soc_code"="o_net_soc_code"))|>
  mutate(across(all_of(skill_cols), ~ .x * down_weight))|>
  group_by(noc_plus_title)|>
  summarise(
    across(all_of(skill_cols), \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  )|>
  column_to_rownames("noc_plus_title")|>
  as.matrix()

onet_pca <- prcomp(noc_skills, center=TRUE, scale=TRUE)

D_full_vec <- scale(noc_skills, center = onet_pca$center, scale = onet_pca$scale)|>
  dist()|>
  as.vector()
max_k <- ncol(onet_pca$x)

k_vs_spearman <- map_dfr(1:max_k, function(k) {
  X_k <- onet_pca$x[, 1:k, drop = FALSE] #keep as a matrix even if k=1
  D_k <- dist(X_k)
  D_k_vec <- as.vector(D_k)

  tibble(
    k = k,
    spearman = cor(D_full_vec, D_k_vec, method = "spearman")
  )
}
)

retain_num_pca <- k_vs_spearman|>
  filter(spearman>.99)|>
  filter(spearman==min(spearman))|>
  pull(k)

noc_coords <- onet_pca$x[, 1:retain_num_pca]#trim redundant
skills_noc_dist<- dist(noc_coords, method = "euclidean")|>
  as.matrix()

#hierarchical distance--------------------------------------------
hier$hierarchy_distance <- crossing(origin=mapped_nocs$noc_plus_title,
                               destination=mapped_nocs$noc_plus_title)|>
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

# binary distance------------------------------------------

binary_mat <- crossing(origin=mapped_nocs$noc_plus_title,
                                    destination=mapped_nocs$noc_plus_title)|>
  mutate(distance=if_else(origin==destination, 0,1))|>
  pivot_wider(names_from = destination, values_from = distance)|>
  column_to_rownames("origin")|>
  as.matrix()

#destination gating----------------------------------

cip_noc <- bcgovpond::read_view("9810040401.csv") |>
  select(highest = starts_with("Highest"),
         cip     = starts_with("Major"),
         noc     = starts_with("Occupation"),
         value   = VALUE) |>
  filter(highest != "Total - Highest certificate, diploma or degree")|>
  mutate(noc = str_replace(noc, "^(.{5})", "\\1:"),
         noc = str_replace(noc, "Seniors", "Senior")
         )|>
  inner_join(mapped_nocs, by = c("noc"="noc_plus_title"))


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
    TEER = str_sub(noc, 2, 2),
    tertile = ntile(specificity, 3),
    sub_regime = case_when(TEER %in% c(0,5) ~ 1, #TEER 0 and 5 governed by skills
                       tertile==3 & TEER %in% c(2,3,4) ~ 2, #even if highly specific education, TEERS 2,3,4 can move within hierarchy
                       TRUE ~ tertile #otherwise regime determined by education specificity
                      ),
    sub_regime = factor(
      sub_regime,
      labels = c("Horizontal (Skill)", "Vertical (Hierarchy)", "Minimal (Binary)")
    )
  )

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
s_all <- offdiag(skills_noc_dist)

# Define the "non-maximal" hierarchical subset (exclude max distance and zeros)
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

skills$skills_noc_dist <- skills_noc_dist
skills$onet_pca <- onet_pca
skills$mapped_nocs <- mapped_nocs

write_rds(skills, here("out", "skills.rds"))
write_rds(hier, here("out", "hier.rds"))
write_rds(binary_mat, here("out", "binary_mat.rds"))
write_rds(noc_specificity, here("out", "noc_specificity.rds"))
write_rds(educ_specificity, here("out", "educ_specificity.rds"))
write_rds(calibration, here("out", "calibration.rds"))
write_rds(k_vs_spearman, here("out", "k_vs_spearman.rds"))




