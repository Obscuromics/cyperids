library(tidyverse)
library(dplyr)

files <- list.files(
  pattern = "_chunks_augmented_CDS.csv$",
  full.names = TRUE
)

chunks <- files %>%
  map_dfr(~ {
    read_csv(.x, show_col_types = FALSE) %>%
      mutate(
        species = str_remove(basename(.x), "_chunks_augmented_CDS.csv")
      )
  })

chunks <- chunks %>%
  mutate(
    array_chunk = if_else(
      !is.na(repeat_classes) & repeat_classes != "",
      1L, 0L
    )
  )

#keeping only cleanest TE data in dedicated column

chunks <- chunks %>%
  mutate(
    TE_class_clean = case_when(
      str_detect(TE_class, "^(LTR|LINE|SINE|DNA|RC)") ~ TE_class,
      TRUE ~ NA_character_
    ),
    TE_overlap_clean = if_else(!is.na(TE_class_clean), 1L, 0L)
  )

#compute distance to nearest array

chunks <- chunks %>%
  group_by(species, chromosome) %>%
  arrange(chunk_id) %>%
  mutate(
    dist_to_array = {
      array_idx <- which(array_chunk == 1)
      if (length(array_idx) == 0) {
        rep(NA_integer_, n())
      } else {
        sapply(seq_along(chunk_id), function(i) {
          min(abs(i - array_idx))
        })
      }
    }
  ) %>%
  ungroup()

#summarise CDS presence by distance

decay_data <- chunks %>%
  filter(!is.na(dist_to_array)) %>%
  group_by(species, dist_to_array) %>%
  summarise(
    prop_CDS = mean(CDS_overlap),
    n = n(),
    .groups = "drop"
  )


#plot with smoothers
ggplot(decay_data, aes(dist_to_array * 5, prop_CDS)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(
    method = "loess",
    span = 0.3,
    se = FALSE
  ) +
  facet_wrap(~ species, scales = "free_y") +
  coord_cartesian(xlim = c(0, 5000)) +
  labs(
    x = "Distance to nearest oligocentromeric array (kb)",
    y = "Proportion of chunks with CDS"
  ) +
  theme_classic()

#adding mean inter-array gap length

species_chunk_data <- read.csv("chunk_stats_plus_rearrangements.csv")

decay_data2 <- decay_data %>%
  left_join(
    species_chunk_data %>%
      select(species, mean_gap_kb),
    by = "species"
  )

decay_data2 <- decay_data2 %>%
  filter(dist_to_array * 5 <= mean_gap_kb)

#removing the first 5kb by an array

decay_data3 <- decay_data2 %>%
  filter(dist_to_array * 5 > 5,
         dist_to_array * 5 <= mean_gap_kb)

ggplot(decay_data3, aes(dist_to_array * 5, prop_CDS)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(
    method = "loess",
    span = 0.3,
    se = FALSE
  ) +
  geom_vline(
    data = species_chunk_data,
    aes(xintercept = mean_gap_kb / 2),
    colour = "red",
    linetype = "dotted"
  ) +
  facet_wrap(~ species, scales = "free_x") +
  labs(
    x = "Distance to nearest oligocentromeric array (kb) (w/o first 5 kb)",
    y = "Proportion of chunks with CDS"
  ) +
  theme_classic()


######### Now for TEs - TE enrichment rather than gene decay

enrichment_data <- chunks %>%
  filter(!is.na(dist_to_array)) %>%
  group_by(species, dist_to_array) %>%
  summarise(
    prop_TE = mean(TE_overlap_clean),
    n = n(),
    .groups = "drop"
  )

enrichment_data2 <- enrichment_data %>%
  left_join(
    species_chunk_data %>%
      select(species, mean_gap_kb),
    by = "species"
  )

enrichment_data2 <- enrichment_data2 %>%
  filter(dist_to_array * 5 <= mean_gap_kb)

#### without the first 5 kb

enrichment_data3 <- enrichment_data2 %>%
  filter(dist_to_array * 5 > 5,
         dist_to_array * 5 <= mean_gap_kb)

ggplot(enrichment_data3, aes(dist_to_array * 5, prop_TE)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(
    method = "loess",
    span = 0.3,
    se = FALSE
  ) +
  geom_vline(
    data = species_chunk_data,
    aes(xintercept = mean_gap_kb / 2),
    colour = "red",
    linetype = "dotted"
  ) +
  facet_wrap(~ species, scales = "free_x") +
  labs(
    x = "Distance to nearest oligocentromeric array (kb) (w/o first 5 kb)",
    y = "Proportion of chunks with TEs"
  ) +
  theme_classic()


#diagnostic plot - number of chunks contributing at each distance
ggplot(decay_data, aes(dist_to_array * 5, n)) +
  geom_line() +
  facet_wrap(~ species, scales = "free_y") +
  labs(
    x = "Distance to array (kb)",
    y = "Number of chunks"
  ) +
  theme_bw()

################### Statistical tests

##### Per-species GLMs

library(broom)

CDS_slopes <- chunks %>%
  filter(!is.na(dist_to_array)) %>%
  mutate(
    dist_kb = dist_to_array * 5
  ) %>%
  left_join(
    species_chunk_data %>%
      select(species, mean_gap_kb),
    by = "species"
  ) %>%
  filter(
    dist_kb > 5,                       # remove first 5 kb
    dist_kb <= mean_gap_kb / 2          # restrict to half inter-array gap
  ) %>%
  mutate(
    log_dist = log(dist_kb + 1)
  ) %>%
  group_by(species) %>%
  do(
    tidy(
      glm(
        CDS_overlap ~ log_dist,
        family = binomial,
        data = .
      ),
      conf.int = TRUE
    )
  ) %>%
  ungroup() %>%
  filter(term == "log_dist")

#summary
CDS_slopes %>%
  mutate(significant = p.value < 0.05) %>%
  arrange(p.value)

#plotting effect sizes across species

bold_species <- c(
  "Juncus effusus",
  "Juncus inflexus",
  "Juncus squarrosus"
)

library(stringr)

species_labels <- function(x) {
  sapply(x, function(sp) {
    sp_clean <- str_replace_all(sp, "_", " ")
    if (sp %in% bold_species) {
      paste0("bold(italic('", sp_clean, "'))")
    } else {
      paste0("italic('", sp_clean, "')")
    }
  })
}

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

CDS_slopes <- CDS_slopes %>%
  mutate(
    species = factor(species, levels = rev(desired_order))
  )

ggplot(CDS_slopes,
       aes(x = species, y = estimate)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high, group = species),
    width = 0.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_x_discrete(
    labels = function(x) parse(text = species_labels(x))
  ) +
  labs(
    x = "Species",
    y = "Effect of log-distance from the oligocentromere on CDS presence"
  ) +
  theme_classic()


######## TEs

te_slopes <- chunks %>%
  filter(!is.na(dist_to_array)) %>%
  mutate(
    dist_kb = dist_to_array * 5
  ) %>%
  left_join(
    species_chunk_data %>%
      select(species, mean_gap_kb),
    by = "species"
  ) %>%
  filter(
    dist_kb > 5,                       # remove first 5 kb
    dist_kb <= mean_gap_kb / 2          # restrict to half inter-array gap
  ) %>%
  mutate(
    log_dist = log(dist_kb + 1)
  ) %>%
  group_by(species) %>%
  do(
    tidy(
      glm(
        TE_overlap_clean ~ log_dist,
        family = binomial,
        data = .
      ),
      conf.int = TRUE
    )
  ) %>%
  ungroup() %>%
  filter(term == "log_dist")

#summary
te_slopes %>%
  mutate(significant = p.value < 0.05) %>%
  arrange(p.value)

#plotting effect sizes across species
ggplot(te_slopes,
       aes(x = reorder(species, estimate),
           y = estimate)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    x = "Species",
    y = "Effect of log-distance from the oligocentromere on TE presence"
  ) +
  theme_classic()

####### plot

te_slopes <- te_slopes %>%
  mutate(
    species = factor(species, levels = rev(desired_order))
  )

ggplot(te_slopes,
       aes(x = species, y = estimate)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high, group = species),
    width = 0.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_x_discrete(
    labels = function(x) parse(text = species_labels(x))
  ) +
  labs(
    x = "Species",
    y = "Effect of log-distance from the oligocentromere on TE presence"
  ) +
  theme_classic()


################################# for Carex myosuroides, only Kelly repeats


chunks2 <- chunks %>%
  mutate(
    is_177_array = grepl(
      "^177_",
      ifelse(is.na(repeat_classes), "", repeat_classes)
    )
  )

chunks2 <- chunks2 %>%
  filter(species == "Carex_myosuroides") %>%
  group_by(chromosome) %>%
  filter(any(is_177_array)) %>%
  ungroup()

chunks2 <- chunks2 %>%
  filter(chromosome %in% c(
    "CM051465.1",
    "CM051472.1",
    "CM051474.1"
  ))

chunks2 <- chunks2 %>%
  group_by(chromosome) %>%
  mutate(
    dist_to_177_chunks = sapply(
      chunk_id,
      function(x) min(abs(x - chunk_id[is_177_array]))
    ),
    dist_to_177_kb = dist_to_177_chunks * 5
  ) %>%
  ungroup()

cds_decay_177 <- chunks2 %>%
  group_by(dist_to_177_kb) %>%
  summarise(
    cds_density = mean(CDS_overlap),
    n_chunks = n(),
    .groups = "drop"
  )


ggplot(cds_decay_177, aes(x = dist_to_177_kb, y = cds_density)) +
  geom_point(shape = 16, size = 1.5) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 8000)) +
  labs(
    x = "Distance to nearest 177 array (kb)",
    y = "CDS density"
  ) +
  theme_classic()

te_enrichment_177 <- chunks2 %>%
  group_by(dist_to_177_kb) %>%
  summarise(
    te_density = mean(TE_overlap_clean),
    n_chunks = n(),
    .groups = "drop"
  )

ggplot(te_enrichment_177, aes(x = dist_to_177_kb, y = te_density)) +
  geom_point(shape = 16, size = 1.5) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 8000)) +
  labs(
    x = "Distance to nearest 177 array (kb)",
    y = "TE density"
  ) +
  theme_classic()

####### new individual plots for the Juncus species

juncus_species <- c(
  "Juncus_effusus",
  "Juncus_squarrosus",
  "Juncus_inflexus"
)

chunks_juncus <- chunks %>%
  filter(species %in% juncus_species) %>%
  mutate(
    is_repeat = !is.na(repeat_classes)
  )

chunks_juncus <- chunks_juncus %>%
  group_by(species, chromosome) %>%
  filter(any(is_repeat)) %>%
  ungroup()

chunks_juncus <- chunks_juncus %>%
  group_by(species, chromosome) %>%
  mutate(
    dist_to_repeat_chunks = sapply(
      chunk_id,
      function(x) min(abs(x - chunk_id[is_repeat]))
    ),
    dist_to_repeat_kb = dist_to_repeat_chunks * 5
  ) %>%
  ungroup()

cds_decay_juncus <- chunks_juncus %>%
  group_by(species, dist_to_repeat_kb) %>%
  summarise(
    cds_density = mean(CDS_overlap),
    n_chunks = n(),
    .groups = "drop"
  )

ggplot(cds_decay_juncus,
       aes(x = dist_to_repeat_kb, y = cds_density)) +
  geom_point(shape = 16, size = 1.5) +
  coord_cartesian(xlim = c(0, 8000)) +
  facet_wrap(~ species, ncol = 1) +
  labs(
    x = "Distance to nearest repeat (kb)",
    y = "CDS density"
  ) +
  theme_classic()

te_enrichment_juncus <- chunks_juncus %>%
  group_by(species, dist_to_repeat_kb) %>%
  summarise(
    te_density = mean(TE_overlap_clean),
    n_chunks = n(),
    .groups = "drop"
  )

ggplot(te_enrichment_juncus,
       aes(x = dist_to_repeat_kb, y = te_density)) +
  geom_point(shape = 16, size = 1.5) +
  coord_cartesian(xlim = c(0, 8000)) +
  facet_wrap(~ species, ncol = 1) +
  labs(
    x = "Distance to nearest repeat (kb)",
    y = "TE density"
  ) +
  theme_classic()

######### now combined plots

cds_decay_myosuroides <- cds_decay_177

colnames(cds_decay_myosuroides)[1] <- "dist_to_repeat_kb"

cds_decay_all <- bind_rows(
  cds_decay_juncus,
  cds_decay_myosuroides
)


ggplot(cds_decay_all,
       aes(x = dist_to_repeat_kb, y = cds_density)) +
  geom_point(shape = 16, size = 1.3, alpha = 0.7) +
  geom_smooth(se = FALSE, method = "loess", span = 0.4) +
  coord_cartesian(xlim = c(0, 4000)) +
  facet_wrap(~ species, ncol = 1) +
  labs(
    x = "Distance to nearest centromere (kb)",
    y = "CDS density"
  ) +
  theme_classic()

te_enrichment_myosuroides <- te_enrichment_177

colnames(te_enrichment_myosuroides)[1] <- "dist_to_repeat_kb"

te_enrichment_all <- bind_rows(
  te_enrichment_juncus,
  te_enrichment_myosuroides
)


ggplot(te_enrichment_all,
       aes(x = dist_to_repeat_kb, y = te_density)) +
  geom_point(shape = 16, size = 1.3, alpha = 0.7) +
  geom_smooth(se = FALSE, method = "loess", span = 0.4) +
  coord_cartesian(xlim = c(0, 4000)) +
  facet_wrap(~ species, ncol = 1) +
  labs(
    x = "Distance to nearest centromere (kb)",
    y = "TE density"
  ) +
  theme_classic()


########## statistically testing CDS decay and TE enrichment by these neo-monocentromeres (within 4000 kb)

carex_test <- chunks2 %>%
  filter(dist_to_177_kb <= 4000)
View(chunks2)

cds_model <- glm(
  CDS_overlap ~ dist_to_177_kb,
  data = carex_test,
  family = binomial
)

summary(cds_model)

te_model <- glm(
  TE_overlap_clean ~ dist_to_177_kb,
  data = carex_test,
  family = binomial
)

summary(te_model)
