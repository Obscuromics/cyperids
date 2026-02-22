library(ape)
library(dplyr)
library(ggplot2)

tree <- read.tree("Carex.txt")

df <- read.csv("fissions_fusions_plus_Carex.csv", stringsAsFactors=F)

#combining inferred number of rearrangements with the Carex phylogeny

node_labels <- c(tree$tip.label, tree$node.label)

edges_df <- data.frame(
  parent = node_labels[tree$edge[,1]],
  child  = node_labels[tree$edge[,2]],
  edge_index = seq_len(nrow(tree$edge)),
  stringsAsFactors = FALSE
)

edges_joined <- edges_df %>%
  left_join(df, by = "child")

rearr_total <- edges_joined$rearrangements

root_chr <- 26 #number of Carex ALGs, as a starting point for the random walks.

#### simulating tip karyotypes with fixed probability of fission and fusion calculated from observation

simulate_tip_karyotypes_unequalprob <- function(tree, rearr_total, root_chr) {
  
  n_nodes <- length(tree$tip.label) + tree$Nnode
  chr <- numeric(n_nodes)
  
  # identify root
  root_node <- setdiff(tree$edge[,1], tree$edge[,2])
  if (length(root_node) != 1) {
    stop("Tree does not appear to be singly rooted")
  }
  
  chr[root_node] <- root_chr
  
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i,1]
    child  <- tree$edge[i,2]
    
    n_events <- rearr_total[i]
    
    delta <- if (n_events > 0) {
      sum(sample(c(-1, 1), n_events, replace = TRUE, prob = c(0.54079254079, 0.4592074592))) #fusions more common than fissions
    } else {
      0
    }
    
    chr[child] <- chr[parent] + delta
  }
  
  chr[seq_along(tree$tip.label)]
}

n_sim <- 10000
null_variance_unequalprob <- numeric(n_sim)

for (i in seq_len(n_sim)) {
  tips <- simulate_tip_karyotypes_unequalprob(tree, rearr_total, root_chr)
  null_variance_unequalprob[i] <- var(tips)
}

df_unequalprob <- data.frame(variance = null_variance_unequalprob)

obs_variance <- 55.813

mean(df_unequalprob$variance <= obs_variance, na.rm = TRUE) #empirical quantile

#visualisation
ggplot(df_unequalprob, aes(x = variance)) +
  geom_density(fill = "grey85", colour = "grey30") +
  geom_vline(xintercept = obs_variance,
             colour = "red", linewidth = 1) +
  labs(
    x = "Variance in tip chromosome number",
    y = "Density"
  ) +
  theme_classic()

#### simulating tip karyotypes with fusion probability drawn from a Beta distribution parameterised using Beta-Binomial likelihood

Carex_branch_rearrangements <- read.csv("Carex_branch_rearrangements_df.csv")

#calculating proportion of fusions out of total rearrangements along each branch
pfusions_df <- Carex_branch_rearrangements %>%
  mutate(
    N = rearrangements,
    p_obs = ifelse(N > 0, fusions / N, NA_real_)
  ) %>%
  filter(!is.na(p_obs))

#making sure no values are exactly 0 or 1 as Beta distribution defined on (0,1)
eps <- 1e-6
p_obs_trim <- pmin(pmax(pfusions_df$p_obs, eps), 1 - eps)

#estimating alpha and beta hyperparameters using Beta-Binomial likelihood
loglik_beta_binom <- function(par, F, N) {
  alpha <- par[1]
  beta  <- par[2]
  
  if (alpha <= 0 || beta <= 0) return(Inf)
  
  -sum(
    lbeta(F + alpha, N - F + beta) -
      lbeta(alpha, beta)
  )
}

fit <- optim(
  par = c(2, 2),
  fn = loglik_beta_binom,
  F = Carex_branch_rearrangements$fusions,
  N = Carex_branch_rearrangements$rearrangements,
  method = "L-BFGS-B",
  lower = c(1e-3, 1e-3)
)

alpha_hat <- fit$par[1]
beta_hat  <- fit$par[2]

#plotting the shape of the Beta distribution with these hyperparameters
curve(dbeta(x, alpha_hat, beta_hat),
      from = 0, to = 1,
      lwd = 2,
      xlab = "Fusion probability",
      ylab = "Density",
      main = "Beta prior estimated via Beta–Binomial model")

#running the simulation

simulate_tip_karyotypes_branch_hetero_beta <- function(tree, rearr_total, p_i, root_chr) {
  
  n_nodes <- length(tree$tip.label) + tree$Nnode
  chr <- numeric(n_nodes)
  
  root_node <- setdiff(tree$edge[,1], tree$edge[,2])
  chr[root_node] <- root_chr
  
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i,1]
    child  <- tree$edge[i,2]
    
    n_events <- rearr_total[i]
    
    delta <- if (n_events > 0) {
      sum(sample(c(-1, 1),
                 n_events,
                 replace = TRUE,
                 prob = c(p_i[i], 1 - p_i[i])))
    } else {
      0
    }
    
    chr[child] <- chr[parent] + delta
  }
  
  chr[seq_along(tree$tip.label)]
}

n_sim <- 10000
null_variance_EB_resampled_BetaBinom <- numeric(n_sim)

for (k in seq_len(n_sim)) {
  
  p_i_sim <- rbeta(length(rearr_total),
                   shape1 = alpha_hat,
                   shape2 = beta_hat)
  
  tips <- simulate_tip_karyotypes_branch_hetero_beta(
    tree = tree,
    rearr_total = rearr_total,
    p_i = p_i_sim,
    root_chr = root_chr
  )
  
  null_variance_EB_resampled_BetaBinom[k] <- var(tips)
}

df_beta <- data.frame(variance = null_variance_EB_resampled_BetaBinom)

mean(df_beta$variance <= obs_variance, na.rm = TRUE) #empirical quantile

ggplot(df_beta, aes(x = variance)) +
  geom_density(fill = "grey85", colour = "grey30") +
  geom_vline(xintercept = obs_variance,
             colour = "red", linewidth = 1) +
  labs(
    x = "Variance in tip chromosome number",
    y = "Density"
  ) +
  theme_classic()


#### making combined plot

df_both <- bind_rows(
  data.frame(
    variance = null_variance_unequalprob,
    model = "Fixed fusion probability"
  ),
  data.frame(
    variance = null_variance_EB_resampled_BetaBinom,
    model = "Branch-heterogeneous (Beta–Binomial)"
  )
)

ggplot(df_both, aes(x = variance, fill = model, colour = model)) +
  geom_density(alpha = 0.35, linewidth = 1) +
  geom_vline(
    xintercept = obs_variance,
    colour = "red",
    linewidth = 1.2
  ) +
  labs(
    x = "Variance in tip chromosome number",
    y = "Density",
    title = "Null expectations for variance in tip chromosome number",
    fill = "Model",
    colour = "Model"
  ) +
  theme_classic()

ggplot(df_both, aes(x = variance, fill = model, colour = model)) +
  geom_density(alpha = 0.35, linewidth = 1) +
  geom_vline(
    xintercept = obs_variance,
    colour = "black",
    linewidth = 1.2
  ) +
  scale_fill_manual(values = c("#0D0887FF", "#FCA636FF")) +
  scale_colour_manual(values = c("#0D0887FF", "#FCA636FF")) +
  labs(
    x = "Variance in tip chromosome number",
    y = "Density",
    title = "Null expectations for variance in tip chromosome number",
    fill = "Model",
    colour = "Model"
  ) +
  theme_classic()

#### a third simulation - what if fusion-heavy branches are followed by fission-heavy branches, usually? What would a simulation like this look like?

logit <- function(x) log(x / (1 - x))
inv_logit <- function(x) 1 / (1 + exp(-x))

simulate_tip_karyotypes_stabilising <- function(
    tree,
    rearr_total,
    root_chr,
    p0,
    C_opt,
    gamma
) {
  
  n_nodes <- length(tree$tip.label) + tree$Nnode
  chr <- numeric(n_nodes)
  
  root_node <- setdiff(tree$edge[,1], tree$edge[,2])
  chr[root_node] <- root_chr
  
  for (i in seq_len(nrow(tree$edge))) {
    
    parent <- tree$edge[i,1]
    child  <- tree$edge[i,2]
    n_events <- rearr_total[i]
    
    C_parent <- chr[parent]
    
    # stabilising fusion probability
    p_fusion <- inv_logit(
      logit(p0) + gamma * (C_opt - C_parent)
    )
    
    delta <- if (n_events > 0) {
      sum(sample(
        c(-1, 1),
        n_events,
        replace = TRUE,
        prob = c(p_fusion, 1 - p_fusion)
      ))
    } else {
      0
    }
    
    chr[child] <- C_parent + delta
  }
  
  chr[seq_along(tree$tip.label)]
}

## choosing the gamma value by branch-based likelihood

write.csv(Carex_branch_rearrangements, "Carex_branch_rearrangements_df.csv")

Carex_branch_rearrangements <- read.csv("Carex_branch_rearrangements_df.csv") #manually added internal node karyotypes

F <- Carex_branch_rearrangements$fusions
N <- Carex_branch_rearrangements$rearrangements
C <- Carex_branch_rearrangements$C_parent   # chromosome number at branch start
p0    <- 0.54    # baseline fusion probability at optimum
C_opt <- 26      # stabilising optimum

p_branch <- function(gamma, C, p0, C_opt) {
  eta0 <- logit(p0)
  eta  <- eta0 - gamma * (C - C_opt)
  inv_logit(eta)
}

loglik_gamma <- function(gamma, F, N, C, p0, C_opt) {
  
  # enforce gamma >= 0
  if (gamma < 0) return(-Inf)
  
  p <- p_branch(gamma, C, p0, C_opt)
  
  # avoid log(0)
  eps <- 1e-10
  p <- pmin(pmax(p, eps), 1 - eps)
  
  sum(dbinom(F, size = N, prob = p, log = TRUE))
}

fit_gamma <- optim(
  par = 0.01,
  fn  = function(g) -loglik_gamma(g, F, N, C, p0, C_opt),
  method = "L-BFGS-B",
  lower = 0,
  upper = 5
)

gamma_hat <- fit_gamma$par
gamma_hat

#gamma_hat is 0! 

## choosing the gamma value by simulation-based calibration
gamma_grid <- seq(0, 0.2, by = 0.01)
loglik_gamma <- numeric(length(gamma_grid))

for (g in seq_along(gamma_grid)) {
  
  gamma <- gamma_grid[g]
  
  sim_var <- replicate(
    3000,
    var(simulate_tip_karyotypes_stabilising(
      tree, rearr_total, root_chr,
      p0 = 0.54,
      C_opt = 26,
      gamma = gamma
    ))
  )
  
  dens <- density(sim_var)
  loglik_gamma[g] <- log(
    approx(dens$x, dens$y, xout = obs_variance)$y
  )
}

gamma_hat <- gamma_grid[which.max(loglik_gamma)]

## running the simulation with chosen gamma
n_sim <- 10000
null_variance_stabilising <- numeric(n_sim)

for (k in seq_len(n_sim)) {
  
  tips <- simulate_tip_karyotypes_stabilising(
    tree = tree,
    rearr_total = rearr_total,
    root_chr = root_chr,
    p0 = 0.54,          # baseline
    C_opt = 26,         # optimum
    gamma = 0.03        # tuning parameter
  )
  
  null_variance_stabilising[k] <- var(tips)
}

df_stabilising <- data.frame(variance = null_variance_stabilising)

ggplot(df_stabilising, aes(x = variance)) +
  geom_density(fill = "grey85", colour = "grey30") +
  geom_vline(xintercept = obs_variance,
             colour = "red", linewidth = 1) +
  labs(
    x = "Variance in tip chromosome number",
    y = "Density"
  ) +
  theme_classic()

#### three-dist. combined plot

df_all <- bind_rows(
  data.frame(
    variance = null_variance_unequalprob,
    model = "Fixed fusion probability"
  ),
  data.frame(
    variance = null_variance_EB_resampled_BetaBinom,
    model = "Branch-heterogeneous (Beta–Binomial)"
  ),
  data.frame(
    variance = null_variance_stabilising,
    model = "Stabilising bias"
  )
)

ggplot(df_all, aes(x = variance, fill = model, colour = model)) +
  geom_density(alpha = 0.35, linewidth = 1) +
  geom_vline(
    xintercept = obs_variance,
    colour = "red",
    linewidth = 1.2
  ) +
  labs(
    x = "Variance in tip chromosome number",
    y = "Density",
    title = "Null expectations for variance in tip chromosome number",
    fill = "Model",
    colour = "Model"
  ) +
  theme_classic()


obs_variance

######### fitting the "finite-memory heterogeneity" model. 

p_hat <- 0.54079254079
eta0 <- qlogis(0.54079254079)

df <- Carex_branch_rearrangements
df$branch_id <- seq_len(nrow(df))
child_to_branch <- setNames(df$branch_id, df$child)
parent_branch <- child_to_branch[ df$parent ]
keep <- !is.na(parent_branch)
parent <- parent_branch[keep]

eps <- 1e-6
p_i_hat <- pmin(pmax(df$pi_hat, eps), 1 - eps)
eta_hat_all <- qlogis(p_i_hat)
eta_hat <- eta_hat_all[keep]

anyNA(eta_hat)
anyNA(parent)
range(parent)
length(eta_hat)

loglik_rho_sigma <- function(par, eta, parent, eta0) {
  
  rho   <- par[1]
  sigma <- par[2]
  
  if (!is.finite(rho) || !is.finite(sigma) ||
      rho < 0 || rho > 1 || sigma <= 0)
    return(-1e12)
  
  mu <- (1 - rho) * eta0 + rho * eta[parent]
  
  if (any(!is.finite(mu)))
    return(-1e12)
  
  ll <- dnorm(eta, mean = mu, sd = sigma, log = TRUE)
  
  if (any(!is.finite(ll)))
    return(-1e12)
  
  sum(ll)
}

fit <- optim(
  par = c(rho = 0.3, sigma = 0.2),
  fn  = function(p) -loglik_rho_sigma(p, eta_hat, parent, eta0),
  method = "L-BFGS-B",
  lower = c(0, 1e-3),
  upper = c(1, 2)
)

rho_hat   <- fit$par[1]
sigma_hat <- fit$par[2]

######### supplementary figure - effect of accounting for rearrangement underestimation

simulate_tip_karyotypes_1.5rate <- function(tree, rearr_total, root_chr) {
  
  n_nodes <- length(tree$tip.label) + tree$Nnode
  chr <- numeric(n_nodes)
  
  # identify root
  root_node <- setdiff(tree$edge[,1], tree$edge[,2])
  if (length(root_node) != 1) {
    stop("Tree does not appear to be singly rooted")
  }
  
  chr[root_node] <- root_chr
  
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i,1]
    child  <- tree$edge[i,2]
    
    n_events <- rearr_total[i]*1.5
    
    delta <- if (n_events > 0) {
      sum(sample(c(-1, 1), n_events, replace = TRUE, prob = c(0.54079254079, 0.4592074592))) #fusions more common than fissions
    } else {
      0
    }
    
    chr[child] <- chr[parent] + delta
  }
  
  chr[seq_along(tree$tip.label)]
}

n_sim <- 10000
null_variance_1.5rate <- numeric(n_sim)

for (i in seq_len(n_sim)) {
  tips <- simulate_tip_karyotypes_1.5rate(tree, rearr_total, root_chr)
  null_variance_1.5rate[i] <- var(tips)
}

df_1.5rate <- data.frame(variance = null_variance_1.5rate)

## rate doubled

simulate_tip_karyotypes_doublerate <- function(tree, rearr_total, root_chr) {
  
  n_nodes <- length(tree$tip.label) + tree$Nnode
  chr <- numeric(n_nodes)
  
  # identify root
  root_node <- setdiff(tree$edge[,1], tree$edge[,2])
  if (length(root_node) != 1) {
    stop("Tree does not appear to be singly rooted")
  }
  
  chr[root_node] <- root_chr
  
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i,1]
    child  <- tree$edge[i,2]
    
    n_events <- rearr_total[i]*2
    
    delta <- if (n_events > 0) {
      sum(sample(c(-1, 1), n_events, replace = TRUE, prob = c(0.54079254079, 0.4592074592))) #fusions more common than fissions
    } else {
      0
    }
    
    chr[child] <- chr[parent] + delta
  }
  
  chr[seq_along(tree$tip.label)]
}

n_sim <- 10000
null_variance_doublerate <- numeric(n_sim)

for (i in seq_len(n_sim)) {
  tips <- simulate_tip_karyotypes_doublerate(tree, rearr_total, root_chr)
  null_variance_doublerate[i] <- var(tips)
}

df_doublerate <- data.frame(variance = null_variance_doublerate)

## combined plot

df_all_rates <- bind_rows(
  data.frame(
    variance = null_variance_unequalprob,
    model = "Inferred rearrangement rate"
  ),
  data.frame(
    variance = null_variance_1.5rate,
    model = "Inferred rearrangement rate * 1.5"
  ),
  data.frame(
    variance = null_variance_doublerate,
    model = "Inferred rearrangement rate * 2"
  )
)

ggplot(df_all_rates, aes(x = variance, fill = model, colour = model)) +
  geom_density(alpha = 0.35, linewidth = 1) +
  geom_vline(
    xintercept = obs_variance,
    colour = "red",
    linewidth = 1.2
  ) + scale_fill_manual(values = c("#0D0887FF", "#B12A90FF", "#FCA636FF")) +
  scale_colour_manual(values = c("#0D0887FF", "#B12A90FF", "#FCA636FF")) +
  labs(
    x = "Variance in tip chromosome number",
    y = "Density",
    title = "Null expectations for variance in tip chromosome number",
    fill = "Rearrangement rate",
    colour = "Rearrangement rate"
  ) +
  theme_classic()

#quantiles

mean(df_1.5rate$variance <= obs_variance, na.rm = TRUE)
mean(df_doublerate$variance <= obs_variance, na.rm = TRUE)

#maximum likelihood estimation of the rate-multiplying parameter

simulate_tip_karyotypes <- function(tree, rearr_total, root_chr, mult) {
  
  n_nodes <- length(tree$tip.label) + tree$Nnode
  chr <- numeric(n_nodes)
  
  # identify root
  root_node <- setdiff(tree$edge[,1], tree$edge[,2])
  if (length(root_node) != 1) {
    stop("Tree does not appear to be singly rooted")
  }
  
  chr[root_node] <- root_chr
  
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i,1]
    child  <- tree$edge[i,2]
    
    n_events <- round(rearr_total[i] * mult)
    
    delta <- if (n_events > 0) {
      sum(sample(c(-1, 1), n_events, replace = TRUE,
                 prob = c(0.54079254079, 0.4592074592)))
    } else {
      0
    }
    
    chr[child] <- chr[parent] + delta
  }
  
  chr[seq_along(tree$tip.label)]
}

sim_var_dist <- function(mult, n_sim = 2000) {
  replicate(n_sim, {
    tips <- simulate_tip_karyotypes(tree, rearr_total, root_chr, mult)
    var(tips)
  })
}

loglik_mult <- function(mult) {
  sims <- sim_var_dist(mult, n_sim = 2000)
  
  log_sims <- log(sims)
  mu  <- mean(log_sims)
  sdv <- sd(log_sims)
  
  dnorm(log(obs_variance), mean = mu, sd = sdv, log = TRUE)
}

fit <- optimize(function(m) -loglik_mult(m),
                interval = c(0.1, 10))

fit$minimum

## combined plot of original rearrangement rate, *1.611058 rate, and double rate

# first, simulation run for *1.611058 rate

n_sim <- 10000
null_variance_MLErate <- numeric(n_sim)

for (i in seq_len(n_sim)) {
  tips <- simulate_tip_karyotypes(tree, rearr_total, root_chr, mult = 1.611058)
  null_variance_MLErate[i] <- var(tips)
}

df_MLErate <- data.frame(variance = null_variance_MLErate)

df_all_rates_alt <- bind_rows(
  data.frame(
    variance = null_variance_unequalprob,
    model = "Inferred rearrangement rate"
  ),
  data.frame(
    variance = null_variance_MLErate,
    model = "Inferred rearrangement rate * 1.611058"
  ),
  data.frame(
    variance = null_variance_doublerate,
    model = "Inferred rearrangement rate * 2"
  )
)

ggplot(df_all_rates_alt, aes(x = variance, fill = model, colour = model)) +
  geom_density(alpha = 0.35, linewidth = 1) +
  geom_vline(
    xintercept = obs_variance,
    colour = "red",
    linewidth = 1.2
  ) + scale_fill_manual(values = c("#0D0887FF", "#B12A90FF", "#FCA636FF")) +
  scale_colour_manual(values = c("#0D0887FF", "#B12A90FF", "#FCA636FF")) +
  labs(
    x = "Variance in tip chromosome number",
    y = "Density",
    title = "Null expectations for variance in tip chromosome number",
    fill = "Rearrangement rate",
    colour = "Rearrangement rate"
  ) +
  theme_classic()

######### supplementary figure - shape of the Beta distribution overlaid on the histogram of observed p(fusion) values

beta_binom_df <- data.frame(
  p = p_grid,
  density = dbeta(p_grid, alpha_hat, beta_hat)
)

ggplot() +
  geom_histogram(aes(x = p_obs_trim, y = after_stat(density)),
                 bins = 30, fill = "grey80", colour = "black") +
  geom_line(data = beta_binom_df,
            aes(x = p, y = density),
            linewidth = 1.2) +
  labs(
    x = "Fusion probability",
    y = "Density",
    title = "Empirical distribution of branch fusion probabilities\nand fitted Beta prior"
  ) +
  theme_classic()
