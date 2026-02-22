############### Rearrangement heat-phylogeny

library(ape)
library(ggtree)
library(dplyr)
library(readr)

rearrangement_data <- read.csv("fissions_fusions_plus.csv")

tree <- read.tree("cyperid33.newick.txt")

node_labels <- read_tsv(
  "cyperid33.clusters.tsv",
  col_names = FALSE
)[, 2:3] %>%
  rename(node = X2, karyotype = X3)

internal_nodes <- node_labels %>%
  filter(grepl("^n[0-9]+$", node))

internal_df <- tibble(
  node = tree$node.label
) %>%
  left_join(internal_nodes, by = "node") %>%
  filter(!is.na(karyotype)) %>%
  mutate(type = "internal")

View(internal_df)

tip_labels <- read.csv("tip_labels_cyperids.csv")

tip_df <- tibble(
  node = tree$tip.label
) %>%
  left_join(
    tip_labels,
    by = c("node" = "species")
  ) %>%
  transmute(
    node = node,
    karyotype = haploid_number,
    type = "tip"
  )

node_karyotype_table <- bind_rows(internal_df, tip_df)

rearrangement_data <- rearrangement_data %>%
  left_join(
    node_karyotype_table %>% select(node, karyotype),
    by = c("child" = "node")
  )

head (rearrangement_data)

##

rearrangement_data <- rearrangement_data %>%
  rename(parent_node = parent) %>% #ggtree doesn't like the word parent
  mutate(node = child) #change name of "child" column to "node" as ggtree expects node-based mapping.

setdiff(rearrangement_data$node, c(tree$tip.label, tree$node.label))

p <- ggtree(tree)
p$data <- p$data %>%
  left_join(rearrangement_data, by = c("label" = "node"))

library(stringr)

p$data <- p$data %>%
  mutate(
    label_clean = ifelse(
      isTip,
      str_replace_all(label, "_", " "),
      label
    )
  ) #remove underscores

##making certain branches black, identified by their child node.

black_nodes <- c(
  "Ananas_comosus",
  "Oryza_sativa",
  "n2",
  "n4",
  "n5",
  "n6"
)

p$data <- p$data %>%
  mutate(
    force_black = label %in% black_nodes
  )

library(ggplot2)

p +
  ## base layer: all branches coloured by rearrangement rate
  geom_tree(
    aes(color = log(rearrangement_rate)),
    size = 1.2
  ) +
  
  ## overlay layer: forced-black branches only
  geom_tree(
    data = function(x) x[x$force_black, ],
    color = "black",
    size = 1.2
  ) +
  
  ## tip labels
  geom_tiplab(
    aes(label = label_clean),
    fontface = "italic",
    size = 3
  ) +
  
  scale_color_viridis_c(
    option = "plasma",
    na.value = "grey80",
    name = "Log of rearrangement rate"
  ) +
  
  ## karyotype circles (tips + internal nodes)
  geom_point(
    data = function(x) x[!is.na(x$karyotype), ],
    aes(x = x, y = y),
    shape = 21,
    fill = "white",
    colour = "black",
    size = 3.8,
    stroke = 0.4
  ) +
  
  ## karyotype numbers inside circles
  geom_text(
    data = function(x) x[!is.na(x$karyotype), ],
    aes(x = x, y = y, label = karyotype),
    colour = "black",
    size = 2
  ) +
  
  theme_tree2()



#### fusion rate heat-phylogeny

p +
  ## base layer: all branches coloured by rearrangement rate
  geom_tree(
    aes(color = log(fusion_rate)),
    size = 1.2
  ) +
  
  ## overlay layer: forced-black branches only
  geom_tree(
    data = function(x) x[x$force_black, ],
    color = "black",
    size = 1.2
  ) +
  
  ## tip labels
  geom_tiplab(
    aes(label = label_clean),
    fontface = "italic",
    size = 3
  ) +
  
  scale_color_viridis_c(
    option = "plasma",
    na.value = "grey80",
    name = "Log of fusion rate"
  ) +
  
  theme_tree2()


#### fission rate heat-phylogeny

p +
  ## base layer: all branches coloured by rearrangement rate
  geom_tree(
    aes(color = log(fission_rate)),
    size = 1.2
  ) +
  
  ## overlay layer: forced-black branches only
  geom_tree(
    data = function(x) x[x$force_black, ],
    color = "black",
    size = 1.2
  ) +
  
  ## tip labels
  geom_tiplab(
    aes(label = label_clean),
    fontface = "italic",
    size = 3
  ) +
  
  scale_color_viridis_c(
    option = "plasma",
    na.value = "grey80",
    name = "Log of fission rate"
  ) +
  
  theme_tree2()


#### Rearrangement violin

violin_data <- read.csv("rearrangements_for_violin.csv")

violin_data <- violin_data %>%
  mutate(group = case_when(
    Family == "Juncaceae" ~ "Juncaceae",
    Family == "Cyperaceae" & Genus == "Carex" ~ "Cyperaceae (Carex)",
    Family == "Cyperaceae" & Genus == "non-Carex" ~ "Cyperaceae (non-Carex)"
  ))

#for plotting order:

violin_data$group <- factor(violin_data$group,
                   levels = c("Juncaceae",
                              "Cyperaceae (non-Carex)",
                              "Cyperaceae (Carex)"))

violin_long <- violin_data %>%
  pivot_longer(
    cols = c(rearr_per_my, fusions_per_my, fissions_per_my),
    names_to = "rate_type",
    values_to = "rate"
  ) %>%
  mutate(rate_type = factor(rate_type,
                            levels = c("rearr_per_my", "fusions_per_my", "fissions_per_my"),
                            labels = c("Rearrangements", "Fusions", "Fissions")
  ))

ggplot(violin_long, aes(x = group, y = rate)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.6, size = 2) +
  facet_wrap(~ rate_type, ncol = 3, scales = "fixed") +
  labs(
    x = "",
    y = "Rate (per Myr)",
    title = "Comparison of chromosomal rearrangement rates"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 13)
  )

