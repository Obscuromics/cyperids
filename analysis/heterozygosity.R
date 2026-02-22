# packages
library(ggplot2)
library(brms)
library(ape)
library(dplyr)

heterozygosity <- read.csv("heterozygosity.csv")

het_hist <- ggplot(data = heterozygosity, aes(x = Heterozygosity)) +
  geom_histogram(bins = 100)

het_hist

######### Phylogenetic regression models

chunks_rearrangement_data_new <- read.csv("chunk_stats_plus_rearrangements.csv")

head(chunks_rearrangement_data_new)

chunks_rearrangement_data_new$Species <-
  gsub(" ", "_", chunks_rearrangement_data_new$Species)

#need to remove species without heterozygosity data

to_remove <- c("Juncus effusus", "Juncus squarrosus", "Juncus inflexus", "Carex_littledalei")
chunks_rearrangement_data_new <- 
  chunks_rearrangement_data_new[ 
    !chunks_rearrangement_data_new$Species %in% to_remove, 
  ]


tree <- read.tree("cyperid33.newick.txt")

tips_to_drop <- setdiff(tree$tip.label, chunks_rearrangement_data_new$Species)

tree_pruned_het <- drop.tip(tree, tips_to_drop)

#re-order data
chunks_rearrangement_data_new <- chunks_rearrangement_data_new[
  match(tree_pruned$tip.label, chunks_rearrangement_data_new$Species),
]

chunks_rearrangement_data_new$Species <-
  factor(chunks_rearrangement_data_new$Species,
         levels = tree_pruned$tip.label)


#phylo covariance matrix
A <- ape::vcv(
  tree_pruned,
  corr = TRUE
)

#adding and scaling the predictor

chunks_rearrangement_data_new <- chunks_rearrangement_data_new %>%
  filter(!is.na(Species)) #removing stray NA rows

chunks_rearrangement_data_new <- chunks_rearrangement_data_new %>%
  left_join(
    heterozygosity %>% select(Species, Heterozygosity),
    by = "Species"
  )

chunks_rearrangement_data_new$heterozygosity_sc <- scale(chunks_rearrangement_data_new$Heterozygosity)

chunks_rearrangement_data_new$rate_eps <-
  chunks_rearrangement_data_new$Proximal.rearrangement.rate + 1e-6
chunks_rearrangement_data_new$fusion_rate_eps <-
  chunks_rearrangement_data_new$Proximal.fusion.rate + 1e-6
chunks_rearrangement_data_new$fission_rate_eps <-
  chunks_rearrangement_data_new$Proximal.fission.rate + 1e-6

library(brms)

het_prior <- c(
  prior(normal(-3, 1.5), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(1), class = "sd", group = "Species"),
  prior(exponential(1), class = "sigma")
)

fit_rearrangements_heterozygosity <- brm(
  rate_eps ~
    heterozygosity_sc +
    (1 | gr(Species, cov = A)),
  
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  prior = het_prior,
  chains = 4,
  cores = 4,
  seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995,
                 max_treedepth = 20)
)

fit_fusions_heterozygosity <- brm(
  fusion_rate_eps ~
    heterozygosity_sc +
    (1 | gr(Species, cov = A)),
  
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  
  prior = het_prior,
  
  chains = 4,
  cores = 4,
  seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995,
                 max_treedepth = 20)
)

fit_fissions_heterozygosity <- brm(
  fission_rate_eps ~
    heterozygosity_sc +
    (1 | gr(Species, cov = A)),
  
  data = chunks_rearrangement_data_new,
  family = lognormal(),
  data2 = list(A = A),
  
  prior = het_prior,
  
  chains = 4,
  cores = 4,
  seed = 12345,
  iter = 6000,
  control = list(adapt_delta = 0.995,
                 max_treedepth = 20)
)

summary(fit_rearrangements_heterozygosity)
summary(fit_fissions_heterozygosity)
summary(fit_fusions_heterozygosity)
