library(ape)
library(brms)

data <- read.csv("Cyperaceae_Juncaceae_rearrangements.csv")

tree <- read.tree("cyperid33.newick.txt")

# keep only species in the dataset
tree <- drop.tip(tree, setdiff(tree$tip.label, data$species))

data <- data[match(tree$tip.label, data$species), ]
rownames(data) <- data$species

# phylogenetic correlation matrix
A <- ape::vcv.phylo(tree, corr = TRUE)

### firstly, for Cyperaceae vs Juncaceae, rearrangements generally

prior <- c(
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd"),      # phylogenetic SD
  prior(exponential(1), class = "shape")    # negbin shape
)

fit <- brm(
  rearrangements ~ family + offset(log(time)) + (1 | gr(species, cov = A)),
  data = data,
  family = negbinomial(),
  data2 = list(A = A),
  prior = prior,
  chains = 4, cores = 4, iter = 6000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(fit)

### now for Cyperaceae vs Juncaceae, fusions only

fit_fusions <- brm(
  fusions ~ family + offset(log(time)) + (1 | gr(species, cov = A)),
  data = data,
  family = negbinomial(),
  data2 = list(A = A),
  prior = prior,
  chains = 4, cores = 4, iter = 6000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(fit_fusions)

### Cyperaceae vs Juncaceae, fissions only

fit_fissions <- brm(
  fissions ~ family + offset(log(time)) + (1 | gr(species, cov = A)),
  data = data,
  family = negbinomial(),
  data2 = list(A = A),
  prior = prior,
  chains = 4, cores = 4, iter = 6000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)


summary(fit_fissions)

###### Carex vs non-Carex model.

data_Carex <- read.csv("Carex_nonCarex_rearrangements.csv")

# keep only species in the dataset
tree_Cyperaceae <- drop.tip(tree, setdiff(tree$tip.label, data_Carex$species))

data_Carex <- data_Carex[match(tree_Cyperaceae$tip.label, data_Carex$species), ]
rownames(data_Carex) <- data_Carex$species

# phylogenetic correlation matrix
A <- ape::vcv.phylo(tree_Cyperaceae, corr = TRUE)

fit_Carex <- brm(
  rearrangements ~ genus + (1 | gr(species, cov = A)),
  data = data_Carex,
  family = negbinomial(),
  data2 = list(A = A),
  prior = prior,
  chains = 4, cores = 4, iter = 6000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(fit_Carex)

fit_Carex_fusions <- brm(
  fusions ~ genus + (1 | gr(species, cov = A)),
  data = data_Carex,
  family = negbinomial(),
  data2 = list(A = A),
  prior = prior,
  chains = 4, cores = 4, iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 15) 
)

summary(fit_Carex_fusions)

fit_Carex_fissions <- brm(
  fissions ~ genus + (1 | gr(species, cov = A)),
  data = data_Carex,
  family = negbinomial(),
  data2 = list(A = A),
  prior = prior,
  chains = 4, cores = 4, iter = 6000,
  control = list(adapt_delta = 0.995, max_treedepth = 15) 
)

summary(fit_Carex_fissions)
