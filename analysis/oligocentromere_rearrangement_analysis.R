##### Packages

library(tidyverse)
library(ggplot2)
library(ape)
library(brms)
library(loo)
library(dplyr)
library(tidyr)


######################### Plotting species-level data for Figure 3

desired_order <- c(
  "Carex_elata",
  "Carex_nigra",
  "Carex_acutiformis",
  "Carex_hirta",
  "Carex_riparia",
  "Carex_rostrata",
  "Carex_distans",
  "Carex_extensa",
  "Carex_laevigata",
  "Carex_pendula",
  "Carex_sylvatica",
  "Carex_depauperata",
  "Carex_caryophyllea",
  "Carex_divulsa",
  "Carex_spicata",
  "Carex_arenaria",
  "Carex_echinata",
  "Carex_littledalei",
  "Carex_myosuroides",
  "Eriophorum_angustifolium",
  "Eriophorum_vaginatum",
  "Scirpus_sylvaticus",
  "Trichophorum_cespitosum",
  "Cyperus_fuscus",
  "Cyperus_rotundus",
  "x_Bolboschoenoplectus_mariqueter",
  "Bolboschoenus_planiculmis",
  "Rhynchospora_breviuscula",
  "Rhynchospora_tenuis",
  "Schoenus_nigricans",
  "Juncus_effusus",
  "Juncus_inflexus",
  "Juncus_squarrosus",
  "Luzula_pallescens",
  "Luzula_sylvatica"
)


#Read in combined dataset of array lengths & gap lengths:
array_gap_lengths <- read.csv("array_gap_lengths.csv")

array_gap_lengths$species <- factor(
  array_gap_lengths$species,
  levels = desired_order
)

remove_species <- c("Juncus_effusus", "Juncus_inflexus", "Juncus_squarrosus", "Cyperus_rotundus"
)

array_gap_lengths <- array_gap_lengths |>
  dplyr::filter(!species %in% remove_species)

#boxplot of just array lengths:

array_lengths <- array_gap_lengths %>%
  filter(type == "array") %>%
  mutate(species = factor(species, levels = desired_order))

ggplot(array_lengths, aes(x = species, y = length_kb, fill = species)) +
  geom_boxplot(outlier.alpha = 0.4) +
  scale_y_log10() +
  scale_fill_viridis_d(
    option = "plasma",
    direction = -1,
    guide = "none"  
  ) +
  labs(
    x = "Species",
    y = "Array length (kb)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

#boxplot of just inter-array gap length

gap_lengths <- array_gap_lengths %>%
  filter(type == "gap") %>%
  mutate(species = factor(species, levels = desired_order))

ggplot(gap_lengths, aes(x = species, y = length_kb, fill = species)) +
  geom_boxplot(outlier.alpha = 0.4) +
  scale_y_log10() +
  scale_fill_viridis_d(
    option = "plasma",
    direction = -1,
    guide = "none"
  ) +
  labs(
    x = "Species",
    y = "Inter-array gap length (kb)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

###### plotting columns of rearrangement rates

chunks_rearrangement_data_new <- read.csv("chunk_stats_plus_rearrangements.csv")

rearrangement_plot_df <- chunks_rearrangement_data_new %>%
  select(
    species,
    Proximal.fusion.rate,
    Proximal.fission.rate
  ) %>%
  pivot_longer(
    cols = c(Proximal.fusion.rate, Proximal.fission.rate),
    names_to = "rate_type",
    values_to = "rate"
  ) %>%
  mutate(
    species = factor(species, levels = rev(desired_order)),
    rate_type = recode(
      rate_type,
      Proximal.fusion.rate  = "Fusion",
      Proximal.fission.rate = "Fission"
    )
  )

ggplot(rearrangement_plot_df, aes(x = rate_type, y = species, fill = log(rate))) +
  geom_tile(color = "grey80", width = 0.9, height = 0.9) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "Rate"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11)
  )


########################### Species-level analyses


to_remove <- c("Juncus_effusus", "Juncus_squarrosus", "Juncus_inflexus", "Cyperus_rotundus")
chunks_rearrangement_data_new <- 
  chunks_rearrangement_data_new[ 
    !chunks_rearrangement_data_new$species %in% to_remove, 
  ]

#Read the tree
tree <- read.tree("cyperid33.newick.txt")

#Change Species column to match tip labels
chunks_rearrangement_data_new$species <-
  gsub(" ", "_", chunks_rearrangement_data_new$species)

#Drop tips not in dataset
tips_to_drop <- setdiff(tree$tip.label, chunks_rearrangement_data_new$species)

tree_pruned <- drop.tip(tree, tips_to_drop)

#re-order data
chunks_rearrangement_data_new <- chunks_rearrangement_data_new[
  match(tree_pruned$tip.label, chunks_rearrangement_data_new$species),
]

#make Species factor
chunks_rearrangement_data_new$species <-
  factor(chunks_rearrangement_data_new$species,
         levels = tree_pruned$tip.label)

#phylo covariance matrix
A <- ape::vcv(
  tree_pruned,
  corr = TRUE
)

#scaling the predictors


chunks_rearrangement_data_new$mean_gap_sc <- scale(chunks_rearrangement_data_new$mean_gap_kb)
chunks_rearrangement_data_new$median_len_sc <- scale(chunks_rearrangement_data_new$median_array_length_kb)
chunks_rearrangement_data_new$oligocentromere_share_sc <- scale(chunks_rearrangement_data_new$Oligocentromere.proportion)
chunks_rearrangement_data_new$array_density_sc <- scale(chunks_rearrangement_data_new$Array.density.5kb)

#ensure there is no zero data so as no problems with lognormal error distribution

chunks_rearrangement_data_new$rate_eps <-
  chunks_rearrangement_data_new$Proximal.rearrangement.rate + 1e-6
chunks_rearrangement_data_new$fusion_rate_eps <-
  chunks_rearrangement_data_new$Proximal.fusion.rate + 1e-6
chunks_rearrangement_data_new$fission_rate_eps <-
  chunks_rearrangement_data_new$Proximal.fission.rate + 1e-6

#### model comparison - fusion

prior_null <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)

prior_full <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 1), class = "Intercept"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)

# Null model - phylogeny alone explaining rates.
fit_null_fusion <- brm(
  fusion_rate_eps ~ 1 + (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior=prior_null,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99)
)

#Individual model - oligocentromere share
fit_share_fusion <- brm(
  fusion_rate_eps ~ oligocentromere_share_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior=prior_full,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

#Individual model - array density
fit_density_fusion <- brm(
  fusion_rate_eps ~ array_density_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior=prior_full,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

#Individual model - gap length
fit_gaps_fusion <- brm(
  fusion_rate_eps ~ mean_gap_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior=prior_full,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

#Individual model - array length
fit_length_fusion <- brm(
  fusion_rate_eps ~ median_len_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior=prior_full,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

#Oligocentromere density & share of the genome. 
fit_amount_fusion <- brm(
  fusion_rate_eps ~ oligocentromere_share_sc + array_density_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior=prior_full,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99)
)

#Oligocentromere organisation - gap between them, and their length. 
fit_geometry_fusion <- brm(
  fusion_rate_eps ~ mean_gap_sc + median_len_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior = prior_full,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99)
)

#Full model - amount of oligocentromere AND structure of arrays (gaps & lengths)
fit_full_fusion <- brm(
  fusion_rate_eps ~ mean_gap_sc + median_len_sc +
    oligocentromere_share_sc + array_density_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior = prior_full,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995)
)

loo_null_fusion <- loo(fit_null_fusion, reloo = TRUE)
loo_share_fusion <- loo(fit_share_fusion, reloo = TRUE)
loo_density_fusion <- loo(fit_density_fusion, reloo = TRUE)
loo_gaps_fusion <- loo(fit_gaps_fusion, reloo = TRUE)
loo_length_fusion <- loo(fit_length_fusion, reloo = TRUE)
loo_amount_fusion <- loo(fit_amount_fusion, reloo = TRUE)
loo_geometry_fusion <- loo(fit_geometry_fusion, reloo = TRUE)
loo_full_fusion <- loo(fit_full_fusion, reloo = TRUE)

loo_compare(
  loo_null_fusion,
  loo_share_fusion,
  loo_density_fusion,
  loo_gaps_fusion,
  loo_length_fusion,
  loo_amount_fusion,
  loo_geometry_fusion,
  loo_full_fusion
)

stack_fusion <- loo_model_weights(
  list(
    null = loo_null_fusion,
    share = loo_share_fusion,
    density = loo_density_fusion,
    gaps = loo_gaps_fusion,
    length = loo_length_fusion,
    amount = loo_amount_fusion,
    geometry = loo_geometry_fusion,
    full = loo_full_fusion
  ),
  method = "stacking"
)

stack_fusion

#Evaluating predictive performance of the array density model with posterior predictive intervals for each species
post_preds_density <- posterior_predict(fit_density_fusion, ndraws = 1000)
pred_lower_density <- apply(post_preds_density, 2, quantile, probs = 0.025)
pred_upper_density <- apply(post_preds_density, 2, quantile, probs = 0.975)
interval_width_density <- pred_upper_density - pred_lower_density
median(interval_width_density)
mean(interval_width_density)

posterior_summary(fit_amount_fusion, pars = "^b_")

#### model comparison - fission.

# Null model - phylogeny alone explaining rates.
fit_null_fission <- brm(
  fission_rate_eps ~ 1 + (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99)
)

#Individual model - oligocentromere share
fit_share_fission <- brm(
  fission_rate_eps ~ oligocentromere_share_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99)
)

#Individual model - array density
fit_density_fission <- brm(
  fission_rate_eps ~ array_density_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99)
)

#Individual model - gap length
fit_gaps_fission <- brm(
  fission_rate_eps ~ mean_gap_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99)
)

#Individual model - array length
fit_length_fission <- brm(
  fission_rate_eps ~ median_len_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99)
)

#Oligocentromere number & share of the genome. More sites available for fissions if they occur on oligocentromeres
fit_amount_fission <- brm(
  fission_rate_eps ~ oligocentromere_share_sc + array_density_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99)
)

#Oligocentromere organisation - gap between them, and their length. 
fit_geometry_fission <- brm(
  fission_rate_eps ~ mean_gap_sc + median_len_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99)
)

#Full model - amount of oligocentromere AND structure of arrays (gaps & lengths)
fit_full_fission <- brm(
  fission_rate_eps ~ mean_gap_sc + median_len_sc +
    oligocentromere_share_sc + array_density_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99)
)

loo_null_fission <- loo(fit_null_fission, reloo = TRUE)
loo_share_fission <- loo(fit_share_fission, reloo = TRUE)
loo_density_fission <- loo(fit_density_fission, reloo = TRUE)
loo_gaps_fission <- loo(fit_gaps_fission, reloo = TRUE)
loo_length_fission <- loo(fit_length_fission, reloo = TRUE)
loo_amount_fission <- loo(fit_amount_fission, reloo = TRUE)
loo_geometry_fission <- loo(fit_geometry_fission, reloo = TRUE)
loo_full_fission <- loo(fit_full_fission, reloo = TRUE)

loo_compare(
  loo_null_fission,
  loo_share_fission,
  loo_density_fission,
  loo_gaps_fission,
  loo_length_fission,
  loo_amount_fission,
  loo_geometry_fission,
  loo_full_fission
)

stack_fission <- loo_model_weights(
  list(
    null = loo_null_fission,
    share = loo_share_fission,
    density = loo_density_fission,
    gaps = loo_gaps_fission,
    length = loo_length_fission,
    amount = loo_amount_fission,
    geometry = loo_geometry_fission,
    full = loo_full_fission
  ),
  method = "stacking"
)

stack_fission

posterior_summary(fit_gaps_fission, pars = "^b_")

#### model comparison - overall rearrangement

rearrangement_priors <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(2), class = "sd", group = "species"),
  prior(exponential(1), class = "sigma")
)

# Null model - phylogeny alone explaining rates.
fit_null_rearrangement <- brm(
  rate_eps ~ 1 + (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior = rearrangement_priors,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

rearrangement_priors_share <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(2), class = "sd", group = "species"),
  prior(exponential(1), class = "sigma"),
  prior(normal(0, 0.5), class = "b", coef = "oligocentromere_share_sc")
)

#Individual model - oligocentromere share
fit_share_rearrangement <- brm(
  rate_eps ~ oligocentromere_share_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  prior = rearrangement_priors_share,
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

#Individual model - array density
fit_density_rearrangement <- brm(
  rate_eps ~ array_density_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior = rearrangement_priors,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

#Individual model - gap length
fit_gaps_rearrangement <- brm(
  rate_eps ~ mean_gap_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior = rearrangement_priors,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

#Individual model - array length
fit_length_rearrangement <- brm(
  rate_eps ~ median_len_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior = rearrangement_priors,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

#Oligocentromere number & share of the genome. More sites available for fissions if they occur on oligocentromeres
rearrangement_priors_amount <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(2), class = "sd", group = "species"),
  prior(exponential(1), class = "sigma"),
  prior(normal(0, 0.5), class = "b", coef = "oligocentromere_share_sc"),
  prior(normal(0, 0.5), class = "b", coef = "array_density_sc")
)

fit_amount_rearrangement <- brm(
  rate_eps ~ oligocentromere_share_sc + array_density_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior = rearrangement_priors_amount,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

#Oligocentromere organisation - gap between them, and their length. 
fit_geometry_rearrangement <- brm(
  rate_eps ~ mean_gap_sc + median_len_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior = rearrangement_priors,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

#Full model - amount of oligocentromere AND structure of arrays (gaps & lengths)
rearrangement_priors_full <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(2), class = "sd", group = "Species"),
  prior(exponential(1), class = "sigma"),
  prior(normal(0, 0.5), class = "b", coef = "oligocentromere_share_sc"),
  prior(normal(0, 0.5), class = "b", coef = "array_density_sc"),
  prior(normal(0, 0.5), class = "b", coef = "mean_gap_sc"),
  prior(normal(0, 0.5), class = "b", coef = "median_len_sc")
)

fit_full_rearrangement <- brm(
  rate_eps ~ mean_gap_sc + median_len_sc +
    oligocentromere_share_sc + array_density_sc +
    (1 | gr(species, cov = A)),
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior = rearrangement_priors,
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.9995, max_treedepth = 20)
)

loo_null_rearrangement <- loo(fit_null_rearrangement, reloo = TRUE)
loo_share_rearrangement <- loo(fit_share_rearrangement, reloo = TRUE)
loo_density_rearrangement <- loo(fit_density_rearrangement, reloo = TRUE)
loo_gaps_rearrangement <- loo(fit_gaps_rearrangement, reloo = TRUE)
loo_length_rearrangement <- loo(fit_length_rearrangement, reloo = TRUE)
loo_amount_rearrangement <- loo(fit_amount_rearrangement, reloo = TRUE)
loo_geometry_rearrangement <- loo(fit_geometry_rearrangement, reloo = TRUE)
loo_full_rearrangement <- loo(fit_full_rearrangement, reloo = TRUE)


loo_compare(
  loo_null_rearrangement,
  loo_share_rearrangement,
  loo_density_rearrangement,
  loo_gaps_rearrangement,
  loo_length_rearrangement,
  loo_amount_rearrangement,
  loo_geometry_rearrangement,
  loo_full_rearrangement
)

stack_rearrangement <- loo_model_weights(
  list(
    null = loo_null_rearrangement,
    share = loo_share_rearrangement,
    density = loo_density_rearrangement,
    gaps = loo_gaps_rearrangement,
    length = loo_length_rearrangement,
    amount = loo_amount_rearrangement,
    geometry = loo_geometry_rearrangement,
    full = loo_full_rearrangement
  ),
  method = "stacking"
)

stack_rearrangement

########## plotting across-species significant results from these models

effects_df <- tibble::tibble(
  relationship = c(
    "Array density ↔ fusion rate",
    "Oligocentromere share ↔ fusion rate",
    "Inter-array gap length ↔ fission rate"
  ),
  estimate = c(-0.5296936, -0.2303389, -2.794599),
  lower    = c(-1.432656, -1.146347, -6.170118),
  upper    = c( 0.3859139,  0.7110758,  0.5306196)
)

ggplot(effects_df,
       aes(x = estimate, y = relationship)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.2
  ) +
  geom_point(size = 3) +
  labs(
    x = "Standardised effect size",
    y = ""
  ) +
  theme_classic()


######################### Plotting chromosome-level trends for Figure 3
chromosome_data_extended <- read.csv("per_chromosome_array_number_extended.csv")

chromosome_data_extended <- chromosome_data_extended %>%
  filter(!species %in% c("Cyperus_rotundus", "Juncus_effusus", "Juncus_inflexus", "Juncus_squarrosus"))

chromosome_data_extended <- chromosome_data_extended %>%
  mutate(non_oligo_chunks = chromosome_size_chunks - num_oligocentromere_chunks)

chromosome_data_extended <- chromosome_data_extended %>%
  mutate(oligo_share = num_oligocentromere_chunks / non_oligo_chunks)

chromosome_data_extended <- chromosome_data_extended %>%
  mutate(array_density = (num_oligocentromere_arrays / non_oligo_chunks) * 2000)

chromosome_data_extended$species <- factor(
  chromosome_data_extended$species,
  levels = desired_order
)

ggplot(chromosome_data_extended, aes(x = species, y = oligo_share, fill = species)) +
  geom_boxplot(outlier.alpha = 0.4) +
  scale_y_log10() +
  scale_fill_viridis_d(
    option = "plasma",
    direction = -1,
    guide = "none" 
  ) +
  labs(
    x = "Species",
    y = "Oligocentromere share of the chromosome"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )


ggplot(chromosome_data_extended, aes(x = species, y = array_density, fill = species)) +
  geom_boxplot(outlier.alpha = 0.4) +
  scale_y_log10() +
  scale_fill_viridis_d(
    option = "plasma",
    direction = -1,
    guide = "none"
  ) +
  labs(
    x = "Species",
    y = "Array density (arrays per 10 Mbp)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

######################### Model comparison with the within-species chromosome-level data

chromosome_data_extended <- read.csv("per_chromosome_array_number_extended.csv")

chromosome_data_extended <- chromosome_data_extended %>%
  mutate(chromosome_size_kb = chromosome_size_chunks * 5)

head(chromosome_data_extended)

chromosome_data_extended <- chromosome_data_extended %>%
  mutate(non_oligo_chromosome_size_kb = non_oligo_chunks * 5)

chromosome_data_extended <- chromosome_data_extended %>%
  filter(!species %in% c("Cyperus_rotundus", "Juncus_effusus", "Juncus_inflexus", "Juncus_squarrosus"))

fusion_fission_status <- read.csv("chromosome_fusion_fission_status.csv")

chromosome_data_extended <- chromosome_data_extended %>%
  left_join(
    fusion_fission_status,
    by = c(
      "species"    = "species",
      "chromosome" = "extant_chromosome"
    )
  )

table(is.na(chromosome_data_extended$fusion))
chromosome_data_extended[!complete.cases(chromosome_data_extended), ]

##need to remove NAs here.
library(ape)

tree <- read.tree("cyperid33.newick.txt")

tips_to_drop <- setdiff(tree$tip.label, chromosome_data_extended$species)

tree_pruned_new <- drop.tip(tree, tips_to_drop)

#phylo covariance matrix
A <- ape::vcv(
  tree_pruned_new,
  corr = TRUE
)

#scaling the predictors

chromosome_data_extended <- chromosome_data_extended %>%
  mutate(oligo_share = num_oligocentromere_chunks / non_oligo_chunks)
chromosome_data_extended <- chromosome_data_extended %>%
  mutate(array_density = num_oligocentromere_arrays / non_oligo_chunks)


chromosome_data_extended$mean_gap_chunks_sc <- scale(chromosome_data_extended$mean_gap_chunks)
chromosome_data_extended$median_array_length_chunks_sc <- scale(chromosome_data_extended$median_array_length_chunks)
chromosome_data_extended$oligo_share_sc <- scale(chromosome_data_extended$oligo_share)
chromosome_data_extended$array_density_sc <- scale(chromosome_data_extended$array_density)

#scaling and logging chromosome size

chromosome_data_extended$log_nonoligo_chrsize_sc <-
  scale(log(chromosome_data_extended$non_oligo_chromosome_size_kb))[,1]

#make sure species match tips of tree

chromosome_data_extended$species <-
  factor(chromosome_data_extended$species,
         levels = tree_pruned_new$tip.label)

#### model comparison - chromosome size

fit_null_chrsize <- brm(
  log_nonoligo_chrsize_sc ~ 1 +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  family = gaussian(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.999, max_treedepth=20)
)

fit_amount_chrsize <- brm(
  log_nonoligo_chrsize_sc ~ array_density_sc + oligo_share_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  family = gaussian(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.999, max_treedepth=20)
)

fit_geometry_chrsize <- brm(
  log_nonoligo_chrsize_sc ~ mean_gap_chunks_sc + median_array_length_chunks_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  family = gaussian(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.999, max_treedepth=20)
)

fit_gaponly_chrsize <- brm(
  log_nonoligo_chrsize_sc ~ mean_gap_chunks_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  family = gaussian(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.999, max_treedepth=20)
)

fit_lengthonly_chrsize <- brm(
  log_nonoligo_chrsize_sc ~ median_array_length_chunks_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  family = gaussian(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.999, max_treedepth=20)
)

fit_densityonly_chrsize <- brm(
  log_nonoligo_chrsize_sc ~ array_density_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  family = gaussian(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.999, max_treedepth=20)
)

fit_shareonly_chrsize <- brm(
  log_nonoligo_chrsize_sc ~ oligo_share_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  family = gaussian(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.999, max_treedepth=20)
)

fit_full_chrsize <- brm(
  log_nonoligo_chrsize_sc ~ array_density_sc + oligo_share_sc + mean_gap_chunks_sc + median_array_length_chunks_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  family = gaussian(),
  data2 = list(A = A),
  save_pars = save_pars(all = TRUE),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.999, max_treedepth=20)
)

loo_null_chrsize <- loo(fit_null_chrsize, reloo = TRUE)
loo_amount_chrsize <- loo(fit_amount_chrsize, reloo = TRUE)
loo_geometry_chrsize <- loo(fit_geometry_chrsize, reloo = TRUE)
loo_gaponly_chrsize <- loo(fit_gaponly_chrsize, reloo = TRUE)
loo_lengthonly_chrsize <- loo(fit_lengthonly_chrsize, reloo = TRUE)
loo_densityonly_chrsize <- loo(fit_densityonly_chrsize, reloo = TRUE)
loo_shareonly_chrsize <- loo(fit_shareonly_chrsize, reloo = TRUE)
loo_full_chrsize <- loo(fit_full_chrsize, reloo = TRUE)

loo_compare(
  loo_null_chrsize,
  loo_amount_chrsize,
  loo_geometry_chrsize,
  loo_gaponly_chrsize,
  loo_lengthonly_chrsize,
  loo_densityonly_chrsize,
  loo_shareonly_chrsize,
  loo_full_chrsize
)

stack_chrsize <- loo_model_weights(
  list(
    null = loo_null_chrsize,
    amount = loo_amount_chrsize,
    geometry = loo_geometry_chrsize,
    gaponly = loo_gaponly_chrsize,
    lengthonly = loo_lengthonly_chrsize,
    densityonly = loo_densityonly_chrsize,
    shareonly = loo_shareonly_chrsize,
    full = loo_full_chrsize
  ),
  method = "stacking"
)

stack_chrsize

#posterior predictive intervals for the share model of chr. size.
post_preds_share <- posterior_predict(fit_shareonly_chrsize, ndraws = 1000)
pred_lower_share <- apply(post_preds_share, 2, quantile, probs = 0.025)
pred_upper_share <- apply(post_preds_share, 2, quantile, probs = 0.975)
interval_width_share <- pred_upper_share - pred_lower_share
median(interval_width_share)
mean(interval_width_share)
#relative to the observed range...
max(chromosome_data_extended$log_nonoligo_chrsize_sc)
min(chromosome_data_extended$log_nonoligo_chrsize_sc)

posterior_summary(fit_geometry_chrsize, pars = "^b_")

#### plotting of the above.

ggplot(chromosome_data_extended,
       aes(x = log10(non_oligo_chromosome_size_kb),
           y = log10(oligo_share),
           colour = species,
           group = species)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d(option = "plasma") +
  theme_classic() +
  labs(
    x = "Log of non-oligocentromeric chromosome size (kb)",
    y = "Log of oligocentromeric share"
  )

ggplot(chromosome_data_extended,
       aes(x = non_oligo_chromosome_size_kb,
           y = oligo_share,
           colour = species,
           group = species)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_viridis_d(option = "plasma") +
  theme_classic() +
  labs(
    x = "Non-oligocentromeric chromosome size (kb, log scale)",
    y = "Oligocentromeric share (log scale)"
  )


ggplot(chromosome_data_extended,
       aes(x = (non_oligo_chromosome_size_kb),
           y = array_density*2000,
           colour = species,
           group = species)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_viridis_d(option = "plasma") +
  theme_classic() +
  labs(
    x = "Log of non-oligocentromeric chromosome size (kb) (log scale)",
    y = "Array density (arrays per 10 Mbp) (log scale)"
  )

ggplot(chromosome_data_extended,
       aes(x = log(non_oligo_chromosome_size_kb),
           y = log(mean_gap_chunks*5),
           colour = species,
           group = species)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d(option = "plasma") +
  theme_classic() +
  labs(
    x = "Log of non-oligocentromeric chromosome size (kb)",
    y = "Log of mean inter-array gap length (kb)"
  )

ggplot(chromosome_data_extended,
       aes(x = log(non_oligo_chromosome_size_kb),
           y = log(median_array_length_chunks*5),
           colour = species,
           group = species)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d(option = "plasma") +
  theme_classic() +
  labs(
    x = "Log of non-oligocentromeric chromosome size (kb)",
    y = "Log of median array length (kb)"
  )

#### fitting models of organisation variables with fission and fusion history as predictors (no model comparison approach here)
#(and controlling for chromosome size)

chromosome_data_extended[!complete.cases(chromosome_data_extended), ]


fit_share_history <- brm(
  oligo_share_sc ~ fusion + fission + log_nonoligo_chrsize_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  data2 = list(A = A),
  family = gaussian(),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

fit_density_history <- brm(
  array_density_sc ~ fusion + fission + log_nonoligo_chrsize_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  data2 = list(A = A),
  family = gaussian(),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

fit_gaps_history <- brm(
  mean_gap_chunks_sc ~ fusion + fission + log_nonoligo_chrsize_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  data2 = list(A = A),
  family = gaussian(),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

fit_length_history <- brm(
  median_array_length_chunks_sc ~ fusion + fission + log_nonoligo_chrsize_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  data2 = list(A = A),
  family = gaussian(),
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

posterior_summary(fit_share_history, pars = "^b_")
posterior_summary(fit_density_history, pars = "^b_")
posterior_summary(fit_gaps_history, pars = "^b_")
posterior_summary(fit_length_history, pars = "^b_")

############ plotting array length & fission relationship

species_to_remove <- c("Juncus_squarrosus", "Juncus_inflexus", "Juncus_effusus", "Cyperus_rotundus")

plot_data <- chromosome_data_extended %>%
  mutate(species = as.character(species)) %>%   # temporarily character
  filter(!species %in% species_to_remove) %>%
  mutate(species = factor(species))

remaining_order <- setdiff(desired_order, species_to_remove)


ggplot(
  plot_data,
  aes(
    x = median_array_length_chunks_sc,
    y = fission,
    colour = species,
    group = species
  )
) +
  geom_point(
    size = 2,
    alpha = 0.4,
    position = position_jitter(height = 0.04)
  ) +
  stat_smooth(
    method = "glm",
    method.args = list(family = binomial),
    se = FALSE
  ) +
  scale_colour_viridis_d(
    option = "plasma",
    limits = rev(remaining_order),
    drop = TRUE
  ) +
  xlab("Median array length (scaled)") +
  ylab("Fission history") +
  theme_classic()

ggplot(
  plot_data,
  aes(
    x = median_array_length_chunks,
    y = fission,
    colour = species,
    group = species
  )
) +
  geom_point(
    size = 2,
    alpha = 0.4,
    position = position_jitter(height = 0.04)
  ) +
  stat_smooth(
    method = "glm",
    method.args = list(family = binomial),
    se = FALSE
  ) +
  scale_colour_viridis_d(
    option = "plasma",
    limits = rev(remaining_order),
    drop = TRUE
  ) +
  xlab("Median array length") +
  ylab("Fission history") +
  theme_classic()


######################## looking at rearrangement history - combined fissions & fusions

chromosome_data_extended <- chromosome_data_extended %>%
  mutate(rearrangement = as.integer(fusion == 1 | fission == 1))

prior_rearr_history <- c(
  prior(normal(0, 1), class = "b"),    
  prior(normal(0, 1.5), class = "Intercept"),
  prior(exponential(2), class = "sd")
)

fit_share_rearr_history <- brm(
  rearrangement ~ oligo_share_sc + log_nonoligo_chrsize_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  family = bernoulli(link = "logit"),
  data2 = list(A = A),
  prior = prior_rearr_history,
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

fit_density_rearr_history <- brm(
  rearrangement ~ array_density_sc + log_nonoligo_chrsize_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  family = bernoulli(link = "logit"),
  data2 = list(A = A),
  prior = prior_rearr_history,
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

fit_gaps_rearr_history <- brm(
  rearrangement ~ mean_gap_chunks_sc + log_nonoligo_chrsize_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  family = bernoulli(link = "logit"),
  data2 = list(A = A),
  prior = prior_rearr_history,
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

fit_length_rearr_history <- brm(
  rearrangement ~ median_array_length_chunks_sc + log_nonoligo_chrsize_sc +
    (1 | gr(species, cov = A)),
  data = chromosome_data_extended,
  family = bernoulli(link = "logit"),
  data2 = list(A = A),
  prior = prior_rearr_history,
  chains = 4, cores = 4, seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

posterior_summary(fit_share_rearr_history, pars = "^b_")
posterior_summary(fit_density_rearr_history, pars = "^b_")
posterior_summary(fit_gaps_rearr_history, pars = "^b_")
posterior_summary(fit_length_rearr_history, pars = "^b_")
