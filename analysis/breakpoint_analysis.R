########################################################################### Synteny plotting function

library(tidyverse)

#' Create a synteny plot between two species
#'
#' @param file1 Path to first TSV file (format: Genus_species_m33_*.tsv)
#' @param file2 Path to second TSV file (format: Genus_species_m33_*.tsv)
#' @param output_prefix Optional custom prefix for output files (default: auto-generated from species names)
#' @param width Plot width in inches (default: 16)
#' @param height Plot height in inches (default: 8)
#' @param save_plot Whether to save the plot to files (default: TRUE)
#' @return A ggplot object
#' @export
plot_synteny <- function(file1, file2, output_prefix = NULL, 
                         width = 16, height = 8, save_plot = TRUE) {
  
  # Extract species names from filenames (part before _m33)
  extract_species <- function(filename) {
    # Get basename without directory
    basename_file <- basename(filename)
    # Extract part between "Genus_" and "_m33" (any genus)
    species <- str_match(basename_file, "^[^_]+_(.+?)_m33")[,2]
    if (is.na(species)) {
      stop("Filename must be in format: Genus_speciesname_m33_*.tsv")
    }
    return(species)
  }
  
  sp1 <- extract_species(file1)
  sp2 <- extract_species(file2)
  
  # Also extract genus name for labels
  extract_genus <- function(filename) {
    basename_file <- basename(filename)
    genus <- str_match(basename_file, "^([^_]+)_")[,2]
    return(genus)
  }
  
  genus1 <- extract_genus(file1)
  genus2 <- extract_genus(file2)
  
  cat("Species 1:", genus1, sp1, "\n")
  cat("Species 2:", genus2, sp2, "\n")
  
  # Read the data files
  buscos_sp1 <- read_tsv(file1, show_col_types = FALSE)
  buscos_sp2 <- read_tsv(file2, show_col_types = FALSE)
  
  # Join the datasets on buscoID to get matched BUSCOs
  matched_buscos <- buscos_sp1 %>%
    inner_join(buscos_sp2, by = "buscoID", suffix = c("_sp1", "_sp2"))
  
  cat("Total BUSCOs in", sp1, ":", nrow(buscos_sp1), "\n")
  cat("Total BUSCOs in", sp2, ":", nrow(buscos_sp2), "\n")
  cat("Matched BUSCOs:", nrow(matched_buscos), "\n")
  
  # Function to order chromosomes by ALG grouping
  order_chromosomes_by_ALG <- function(busco_data, chr_col, alg_col) {
    # For each chromosome, find its dominant ALG
    chr_alg_summary <- busco_data %>%
      group_by(!!sym(chr_col), !!sym(alg_col)) %>%
      summarise(n_buscos = n(), .groups = "drop") %>%
      group_by(!!sym(chr_col)) %>%
      slice_max(n_buscos, n = 1) %>%
      ungroup() %>%
      select(!!sym(chr_col), dominant_alg = !!sym(alg_col))
    
    # Order chromosomes first by ALG, then alphabetically within ALG
    chr_order <- chr_alg_summary %>%
      arrange(dominant_alg, !!sym(chr_col)) %>%
      pull(!!sym(chr_col))
    
    return(chr_order)
  }
  
  # Get chromosome orders
  sp1_chr_order <- order_chromosomes_by_ALG(buscos_sp1, "query_chr", "assigned_chr")
  sp2_chr_order <- order_chromosomes_by_ALG(buscos_sp2, "query_chr", "assigned_chr")
  
  cat("\n", sp1, "chromosome order:", paste(sp1_chr_order, collapse = ", "), "\n")
  cat(sp2, "chromosome order:", paste(sp2_chr_order, collapse = ", "), "\n")
  
  # Create chromosome length information
  sp1_chr_lengths <- buscos_sp1 %>%
    group_by(query_chr) %>%
    summarise(max_pos = max(position)) %>%
    mutate(chr_order = match(query_chr, sp1_chr_order)) %>%
    arrange(chr_order) %>%
    mutate(chr_start = cumsum(c(0, head(max_pos, -1))) + (chr_order - 1) * 1e6,
           chr_end = chr_start + max_pos,
           chr_mid = (chr_start + chr_end) / 2)
  
  sp2_chr_lengths <- buscos_sp2 %>%
    group_by(query_chr) %>%
    summarise(max_pos = max(position)) %>%
    mutate(chr_order = match(query_chr, sp2_chr_order)) %>%
    arrange(chr_order) %>%
    mutate(chr_start = cumsum(c(0, head(max_pos, -1))) + (chr_order - 1) * 1e6,
           chr_end = chr_start + max_pos,
           chr_mid = (chr_start + chr_end) / 2)
  
  # Add genome-wide positions to busco data
  sp1_with_genomepos <- buscos_sp1 %>%
    left_join(sp1_chr_lengths %>% select(query_chr, chr_start), by = "query_chr") %>%
    mutate(genome_pos = chr_start + position)
  
  sp2_with_genomepos <- buscos_sp2 %>%
    left_join(sp2_chr_lengths %>% select(query_chr, chr_start), by = "query_chr") %>%
    mutate(genome_pos = chr_start + position)
  
  # Update matched buscos with genome positions
  matched_with_pos <- matched_buscos %>%
    left_join(sp1_with_genomepos %>% select(buscoID, genome_pos_sp1 = genome_pos), 
              by = "buscoID") %>%
    left_join(sp2_with_genomepos %>% select(buscoID, genome_pos_sp2 = genome_pos), 
              by = "buscoID")
  
  # Create an alternating color palette for ALGs
  algs <- sort(unique(c(buscos_sp1$assigned_chr, buscos_sp2$assigned_chr)))
  n_algs <- length(algs)
  
  # Create two distinct color palettes and interleave them
  palette1 <- rainbow(ceiling(n_algs/2), s = 0.8, v = 0.8, start = 0, end = 0.5)
  palette2 <- rainbow(floor(n_algs/2), s = 0.8, v = 0.7, start = 0.5, end = 1)
  
  # Interleave the palettes for maximum contrast between adjacent ALGs
  alternating_colors <- character(n_algs)
  alternating_colors[seq(1, n_algs, 2)] <- palette1[1:length(seq(1, n_algs, 2))]
  alternating_colors[seq(2, n_algs, 2)] <- palette2[1:length(seq(2, n_algs, 2))]
  
  alg_colors <- setNames(alternating_colors, algs)
  
  # Create the synteny plot
  p <- ggplot() +
    # Ribbons connecting matched BUSCOs
    geom_segment(data = matched_with_pos,
                 aes(x = genome_pos_sp1, xend = genome_pos_sp2,
                     y = 1, yend = 0,
                     color = assigned_chr_sp1),
                 alpha = 0.3, linewidth = 0.3) +
    
    # Points for sp1 BUSCOs
    geom_point(data = sp1_with_genomepos,
               aes(x = genome_pos, y = 1, color = assigned_chr),
               size = 1, alpha = 0.6) +
    
    # Points for sp2 BUSCOs
    geom_point(data = sp2_with_genomepos,
               aes(x = genome_pos, y = 0, color = assigned_chr),
               size = 1, alpha = 0.6) +
    
    # Grey rectangles outlining chromosomes for sp1
    geom_rect(data = sp1_chr_lengths,
              aes(xmin = chr_start, xmax = chr_end,
                  ymin = 0.98, ymax = 1.02),
              fill = NA, color = "grey50", linewidth = 0.5) +
    
    # Grey rectangles outlining chromosomes for sp2
    geom_rect(data = sp2_chr_lengths,
              aes(xmin = chr_start, xmax = chr_end,
                  ymin = -0.02, ymax = 0.02),
              fill = NA, color = "grey50", linewidth = 0.5) +
    
    # Chromosome labels
    geom_text(data = sp1_chr_lengths,
              aes(x = chr_mid, y = 1.05, label = query_chr),
              angle = 45, hjust = 0, size = 2.5) +
    geom_text(data = sp2_chr_lengths,
              aes(x = chr_mid, y = -0.05, label = query_chr),
              angle = 45, hjust = 1, size = 2.5) +
    
    scale_color_manual(values = alg_colors, name = "Ancestral Linkage Group") +
    scale_y_continuous(breaks = c(0, 1), 
                       labels = c(paste0(substr(genus2, 1, 1), ". ", sp2), 
                                  paste0(substr(genus1, 1, 1), ". ", sp1))) +
    labs(title = paste("Synteny between", genus1, sp1, "and", genus2, sp2),
         subtitle = "BUSCOs colored by ancestral linkage group (ALG)",
         x = "Genomic Position",
         y = "") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    coord_cartesian(ylim = c(-0.15, 1.15))
  
  # Save the plot if requested
  if (save_plot) {
    if (is.null(output_prefix)) {
      output_prefix <- paste0(sp1, "_", sp2, "_synteny_plot")
    }
    ggsave(paste0(output_prefix, ".png"), p, width = width, height = height, dpi = 300)
    ggsave(paste0(output_prefix, ".pdf"), p, width = width, height = height)
    cat("\nPlots saved as '", output_prefix, ".png' and '", output_prefix, ".pdf'\n", sep = "")
  }
  
  # Print summary statistics
  cat("\n=== Summary Statistics ===\n")
  cat("Number of unique ALGs:", n_algs, "\n")
  
  crossing_analysis <- matched_with_pos %>%
    arrange(genome_pos_sp1) %>%
    mutate(crossing = genome_pos_sp2 < lag(genome_pos_sp2, default = -Inf)) %>%
    summarise(total_crossings = sum(crossing, na.rm = TRUE))
  
  cat("Approximate number of line crossings:", crossing_analysis$total_crossings, "\n")
  
  # Show ALG distribution across chromosomes
  cat("\n===", substr(genus1, 1, 1), ".", sp1, ": Chromosomes by dominant ALG ===\n", sep = " ")
  buscos_sp1 %>%
    group_by(query_chr, assigned_chr) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(query_chr) %>%
    slice_max(n, n = 1) %>%
    arrange(assigned_chr, query_chr) %>%
    print(n = Inf)
  
  cat("\n===", substr(genus2, 1, 1), ".", sp2, ": Chromosomes by dominant ALG ===\n", sep = " ")
  buscos_sp2 %>%
    group_by(query_chr, assigned_chr) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(query_chr) %>%
    slice_max(n, n = 1) %>%
    arrange(assigned_chr, query_chr) %>%
    print(n = Inf)
  
  cat("\nDone!\n")
  
  # Return the plot object
  return(invisible(p))
}

plot_synteny("Carex_elata_m33_n68_complete_location.tsv", "Carex_nigra_m33_n68_complete_location.tsv")

plot_synteny("Carex_pendula_m33_n57_complete_location.tsv", "Carex_sylvatica_m33_n57_complete_location.tsv")

plot_synteny("Carex_divulsa_m33_n52_complete_location.tsv", "Carex_spicata_m33_n52_complete_location.tsv")

plot_synteny("Carex_littledalei_m33_n39_complete_location.tsv", "Carex_myosuroides_m33_n39_complete_location.tsv")

plot_synteny("Rhynchospora_breviuscula_m33_n15_complete_location.tsv", "Rhynchospora_tenuis_m33_n15_complete_location.tsv")

plot_synteny("Eriophorum_angustifolium_m33_n38_complete_location.tsv", "Eriophorum_vaginatum_m33_n38_complete_location.tsv")

####################################################### Generalisable functions for overlaying breakpoints & oligocentromeres

### Helper functions for the below

# Helper function to normalise species names
normalise_species_name <- function(species_name) {
  # Remove spaces and replace with underscores
  species_name <- gsub(" ", "_", species_name)
  # Ensure format is Genus_species
  if (!grepl("_", species_name)) {
    stop("Species name must be in format 'Genus_species' or 'Genus species'")
  }
  return(species_name)
}

# Helper function to match species name to what's in the data
match_species_name <- function(species_input, chunks_data) {
  # Normalize the input
  normalised_input <- normalise_species_name(species_input)
  
  # Get unique species names from chunks data
  available_species <- unique(chunks_data$species)
  
  # Try exact match first
  if (normalised_input %in% available_species) {
    return(normalised_input)
  }
  
  # If not found, try to match by removing/adding underscores
  # Convert input to both formats
  with_underscore <- gsub(" ", "_", species_input)
  without_underscore <- gsub("_", " ", species_input)
  
  # Check both versions
  if (with_underscore %in% available_species) {
    return(with_underscore)
  }
  if (without_underscore %in% available_species) {
    return(without_underscore)
  }
  
  # If still not found, provide helpful error
  stop(paste0("Species '", species_input, "' not found in chunks data.\n",
              "Available species: ", paste(available_species, collapse = ", ")))
}


#' Create a single chromosome cartoon plot for one species
#'
#' @param chunks_file Either a file path to chunks data (TSV) OR a pre-loaded data frame with columns: species, chromosome, chunk_id, repeat_classes
#' @param species_name Name of the species (e.g., "Carex_distans", "Carex distans", or "distans")
#' @param breakpoints_file Path to breakpoints file (TSV with: chr, break_start, break_end, class)
#' @param other_species_name Name of comparison species (for subtitle) - can also be with/without underscore
#' @param chunk_size Size of chunks in bp (default: 5000)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 10)
#' @param save_plot Whether to save plot to file (default: TRUE)
#' @return A ggplot object
#' @export
plot_single_chromosome_cartoon <- function(chunks_file,
                                           species_name,
                                           breakpoints_file,
                                           other_species_name = NULL,
                                           chunk_size = 5000,
                                           width = 12,
                                           height = 10,
                                           save_plot = TRUE) {
  
  # Read the chunks data (or use provided data frame)
  if (is.character(chunks_file)) {
    # It's a file path, read it
    chunks <- read_tsv(chunks_file, show_col_types = FALSE)
  } else if (is.data.frame(chunks_file)) {
    # It's already a data frame, use it directly
    chunks <- chunks_file
  } else {
    stop("chunks_file must be either a file path (character) or a data frame")
  }
  
  # Match species name to what's actually in the chunks data
  species_matched <- match_species_name(species_name, chunks)
  
  cat("Matched species name in data:\n")
  cat("  Input:", species_name, "-> Matched:", species_matched, "\n")
  
  # Read breakpoints
  breakpoints <- read_tsv(breakpoints_file, show_col_types = FALSE)
  
  # Extract species names
  extract_epithet <- function(sp_name) {
    if (grepl("_", sp_name)) {
      return(str_match(sp_name, "_(.+)$")[,2])
    } else if (grepl(" ", sp_name)) {
      return(str_match(sp_name, " (.+)$")[,2])
    } else {
      return(sp_name)
    }
  }
  
  extract_genus <- function(sp_name) {
    if (grepl("_", sp_name)) {
      return(str_match(sp_name, "^([^_]+)_")[,2])
    } else if (grepl(" ", sp_name)) {
      return(str_match(sp_name, "^([^ ]+) ")[,2])
    } else {
      return("")
    }
  }
  
  current_epithet <- extract_epithet(species_matched)
  genus <- extract_genus(species_matched)
  
  if (!is.null(other_species_name)) {
    other_epithet <- extract_epithet(other_species_name)
  } else {
    other_epithet <- "[other species]"
  }
  
  cat("\nCreating plot for", species_matched, "\n")
  
  # Filter chunks for this species (using matched name)
  chunks_sp <- chunks %>% 
    filter(species == species_matched)
  
  # Calculate chunk start and end positions
  chunks_sp <- chunks_sp %>%
    mutate(chunk_start = (chunk_id - 1) * chunk_size,
           chunk_end = chunk_id * chunk_size)
  
  # Initialize the breakpoint class column
  chunks_sp$breakpoint_class <- NA
  
  # For each breakpoint region, find overlapping chunks
  for(i in 1:nrow(breakpoints)){
    chr <- breakpoints$chr[i]
    break_start <- breakpoints$break_start[i]
    break_end <- breakpoints$break_end[i]
    class <- breakpoints$class[i]
    
    # Find chunks that overlap this breakpoint
    overlapping <- which(
      chunks_sp$chromosome == chr & 
        chunks_sp$chunk_end >= break_start & 
        chunks_sp$chunk_start <= break_end
    )
    
    # Assign the class to overlapping chunks
    chunks_sp$breakpoint_class[overlapping] <- class
  }
  
  # Print summary
  cat("\nSummary:\n")
  cat("Total chunks:", nrow(chunks_sp), "\n")
  cat("Chunks with repeats:", sum(!is.na(chunks_sp$repeat_classes)), "\n")
  cat("Chunks with fusions:", sum(grepl("fusion", chunks_sp$breakpoint_class, ignore.case = TRUE), na.rm = TRUE), "\n")
  cat("Chunks with fissions:", sum(grepl("fission", chunks_sp$breakpoint_class, ignore.case = TRUE), na.rm = TRUE), "\n")
  
  # Prepare the data with colors
  chunks_plot <- chunks_sp %>%
    mutate(
      color = case_when(
        !is.na(repeat_classes) ~ "black",
        grepl("fission", breakpoint_class, ignore.case = TRUE) ~ "red",
        grepl("fusion", breakpoint_class, ignore.case = TRUE) ~ "blue",
        TRUE ~ "lightgray"
      ),
      position_mb = chunk_start / 1e6
    )
  
  # Get chromosome order and sizes
  chr_info <- chunks_plot %>%
    group_by(chromosome) %>%
    summarise(max_pos = max(chunk_end) / 1e6) %>%
    arrange(desc(max_pos))
  
  # Set factor levels for plotting order
  chunks_plot$chromosome <- factor(chunks_plot$chromosome, 
                                   levels = rev(chr_info$chromosome))
  
  # Create the plot
  p <- ggplot(chunks_plot, aes(x = position_mb, y = chromosome)) +
    geom_segment(aes(x = position_mb, xend = position_mb + (chunk_size/1e6),
                     y = chromosome, yend = chromosome, color = color),
                 linewidth = 2) +
    geom_segment(data = filter(chunks_plot, color == "black"),
                 aes(x = position_mb, xend = position_mb + (chunk_size/1e6),
                     y = chromosome, yend = chromosome),
                 color = "black",
                 linewidth = 2) +
    scale_color_identity() +
    labs(x = "Position (Mb)", 
         y = "Chromosome",
         title = paste(genus, current_epithet, "chromosomes"),
         subtitle = paste0("Black: Repeat arrays | Blue: fusions in ", genus, " ", current_epithet, 
                           " | Red: fissions in ", genus, " ", other_epithet)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 8))
  
  # Save plot if requested
  if (save_plot) {
    output_prefix <- paste0(tolower(current_epithet), "_chromosome_cartoons")
    ggsave(paste0(output_prefix, ".pdf"), p, width = width, height = height)
    ggsave(paste0(output_prefix, ".png"), p, width = width, height = height, dpi = 300, bg = "white")
    cat("\nPlot saved as:", output_prefix, ".pdf/.png\n")
  }
  
  return(invisible(p))
}

############################################################################### Function for statistical testing enrichment

# Statistical Test for Oligocentromere Enrichment in Breakpoint Regions

library(tidyverse)

# Helper function to normalize species names (reuse from cartoon function)
normalise_species_name <- function(species_name) {
  species_name <- gsub(" ", "_", species_name)
  if (!grepl("_", species_name)) {
    stop("Species name must be in format 'Genus_species' or 'Genus species'")
  }
  return(species_name)
}

# Helper function to match species name to data
match_species_name <- function(species_input, chunks_data) {
  normalised_input <- normalise_species_name(species_input)
  available_species <- unique(chunks_data$species)
  
  if (normalised_input %in% available_species) {
    return(normalised_input)
  }
  
  with_underscore <- gsub(" ", "_", species_input)
  without_underscore <- gsub("_", " ", species_input)
  
  if (with_underscore %in% available_species) {
    return(with_underscore)
  }
  if (without_underscore %in% available_species) {
    return(without_underscore)
  }
  
  stop(paste0("Species '", species_input, "' not found in chunks data.\n",
              "Available species: ", paste(available_species, collapse = ", ")))
}

#' Test for oligocentromere enrichment in breakpoint regions
#'
#' @param breakpoints_file Path to breakpoints file (TSV with: chr, break_start, break_end, class)
#' @param chunks_file Either a file path to chunks data (TSV) OR a pre-loaded data frame
#' @param species_name Name of the species (e.g., "Carex_distans", "Carex distans", or "distans")
#' @param breakpoint_type Which breakpoints to test: "fusion", "fission", or "both" (default: "both")
#' @param chunk_size Size of chunks in bp (default: 5000)
#' @param verbose Print detailed output (default: TRUE)
#' @return A list with test results and summary statistics
#' @export
test_oligocentromere_enrichment <- function(breakpoints_file,
                                            chunks_file,
                                            species_name,
                                            breakpoint_type = "both",
                                            chunk_size = 5000,
                                            verbose = TRUE) {
  
  # Validate breakpoint_type
  if (!breakpoint_type %in% c("fusion", "fission", "both")) {
    stop("breakpoint_type must be one of: 'fusion', 'fission', or 'both'")
  }
  
  # Read the chunks data (or use provided data frame)
  if (is.character(chunks_file)) {
    chunks <- read_tsv(chunks_file, show_col_types = FALSE)
  } else if (is.data.frame(chunks_file)) {
    chunks <- chunks_file
  } else {
    stop("chunks_file must be either a file path (character) or a data frame")
  }
  
  # Match species name
  species_matched <- match_species_name(species_name, chunks)
  
  if (verbose) {
    cat("===========================================\n")
    cat("Oligocentromere Enrichment Test\n")
    cat("===========================================\n")
    cat("Species:", species_matched, "\n")
    cat("Testing breakpoint type:", breakpoint_type, "\n\n")
  }
  
  # Read breakpoints
  breakpoints <- read_tsv(breakpoints_file, show_col_types = FALSE)
  
  # Filter breakpoints by type if needed
  if (breakpoint_type == "fusion") {
    breakpoints <- breakpoints %>%
      filter(grepl("fusion", class, ignore.case = TRUE))
  } else if (breakpoint_type == "fission") {
    breakpoints <- breakpoints %>%
      filter(grepl("fission", class, ignore.case = TRUE))
  }
  # If "both", use all breakpoints
  
  if (nrow(breakpoints) == 0) {
    stop("No breakpoints found for the specified type")
  }
  
  if (verbose) {
    cat("Number of breakpoint regions:", nrow(breakpoints), "\n")
  }
  
  # Filter chunks for this species
  chunks_sp <- chunks %>%
    filter(species == species_matched)
  
  # Calculate chunk start and end positions
  chunks_sp <- chunks_sp %>%
    mutate(chunk_start = (chunk_id - 1) * chunk_size,
           chunk_end = chunk_id * chunk_size)
  
  # Classify each chunk as in breakpoint region or not
  chunks_sp$in_breakpoint <- FALSE
  
  for (i in 1:nrow(breakpoints)) {
    chr <- breakpoints$chr[i]
    break_start <- breakpoints$break_start[i]
    break_end <- breakpoints$break_end[i]
    
    # Find chunks that overlap this breakpoint
    overlapping <- which(
      chunks_sp$chromosome == chr & 
        chunks_sp$chunk_end >= break_start & 
        chunks_sp$chunk_start <= break_end
    )
    
    chunks_sp$in_breakpoint[overlapping] <- TRUE
  }
  
  # Identify chunks with oligocentromeres (non-NA repeat_classes)
  chunks_sp <- chunks_sp %>%
    mutate(has_oligocentromere = !is.na(repeat_classes))
  
  # Create contingency table
  # Rows: in/out of breakpoint region
  # Cols: has/doesn't have oligocentromere
  contingency_table <- table(
    Breakpoint = chunks_sp$in_breakpoint,
    Oligocentromere = chunks_sp$has_oligocentromere
  )
  
  if (verbose) {
    cat("\n--- Contingency Table ---\n")
    print(contingency_table)
    cat("\n")
  }
  
  # Calculate proportions
  in_breakpoint <- chunks_sp %>% filter(in_breakpoint)
  out_breakpoint <- chunks_sp %>% filter(!in_breakpoint)
  
  prop_in <- mean(in_breakpoint$has_oligocentromere)
  prop_out <- mean(out_breakpoint$has_oligocentromere)
  
  if (verbose) {
    cat("--- Proportions ---\n")
    cat("Chunks in breakpoint regions with oligocentromeres:", 
        sum(in_breakpoint$has_oligocentromere), "/", nrow(in_breakpoint),
        sprintf(" (%.2f%%)\n", prop_in * 100))
    cat("Chunks outside breakpoint regions with oligocentromeres:", 
        sum(out_breakpoint$has_oligocentromere), "/", nrow(out_breakpoint),
        sprintf(" (%.2f%%)\n", prop_out * 100))
    cat("Fold enrichment:", sprintf("%.2f\n", prop_in / prop_out))
    cat("\n")
  }
  
  # Perform Fisher's Exact Test
  fisher_test <- fisher.test(contingency_table)
  
  # Perform Chi-square test (if cell counts are sufficient)
  chisq_test <- NULL
  chisq_warning <- NULL
  tryCatch({
    chisq_test <- chisq.test(contingency_table)
  }, warning = function(w) {
    chisq_warning <<- w$message
  })
  
  if (verbose) {
    cat("--- Fisher's Exact Test ---\n")
    cat("Odds Ratio:", sprintf("%.3f\n", fisher_test$estimate))
    cat("95% CI:", sprintf("[%.3f, %.3f]\n", 
                           fisher_test$conf.int[1], fisher_test$conf.int[2]))
    cat("P-value:", format.pval(fisher_test$p.value, digits = 3), "\n")
    
    if (fisher_test$p.value < 0.001) {
      cat("*** Highly significant enrichment\n")
    } else if (fisher_test$p.value < 0.01) {
      cat("** Significant enrichment\n")
    } else if (fisher_test$p.value < 0.05) {
      cat("* Significant enrichment\n")
    } else {
      cat("Not significant\n")
    }
    cat("\n")
    
    if (!is.null(chisq_test)) {
      cat("--- Chi-square Test ---\n")
      cat("Chi-square statistic:", sprintf("%.3f\n", chisq_test$statistic))
      cat("P-value:", format.pval(chisq_test$p.value, digits = 3), "\n\n")
    } else if (!is.null(chisq_warning)) {
      cat("--- Chi-square Test ---\n")
      cat("Warning:", chisq_warning, "\n")
      cat("(Use Fisher's test for small sample sizes)\n\n")
    }
  }
  
  # Permutation test for additional validation
  if (verbose) {
    cat("--- Permutation Test (10,000 iterations) ---\n")
    cat("Running permutations...")
  }
  
  n_permutations <- 10000
  observed_diff <- prop_in - prop_out
  
  permuted_diffs <- replicate(n_permutations, {
    # Randomly shuffle the breakpoint labels
    shuffled <- sample(chunks_sp$in_breakpoint)
    prop_in_perm <- mean(chunks_sp$has_oligocentromere[shuffled])
    prop_out_perm <- mean(chunks_sp$has_oligocentromere[!shuffled])
    prop_in_perm - prop_out_perm
  })
  
  # Two-tailed p-value
  perm_pvalue <- mean(abs(permuted_diffs) >= abs(observed_diff))
  
  if (verbose) {
    cat(" Done!\n")
    cat("Observed difference in proportions:", sprintf("%.4f\n", observed_diff))
    cat("Permutation p-value:", format.pval(perm_pvalue, digits = 3), "\n")
    cat("\n")
  }
  
  # Summary interpretation
  if (verbose) {
    cat("===========================================\n")
    cat("SUMMARY\n")
    cat("===========================================\n")
    
    if (fisher_test$p.value < 0.05 && prop_in > prop_out) {
      cat("CONCLUSION: Oligocentromeres are SIGNIFICANTLY ENRICHED\n")
      cat("in", breakpoint_type, "breakpoint regions.\n")
    } else if (fisher_test$p.value < 0.05 && prop_in < prop_out) {
      cat("CONCLUSION: Oligocentromeres are SIGNIFICANTLY DEPLETED\n")
      cat("in", breakpoint_type, "breakpoint regions.\n")
    } else {
      cat("CONCLUSION: No significant enrichment or depletion of\n")
      cat("oligocentromeres in", breakpoint_type, "breakpoint regions.\n")
    }
    cat("\n")
  }
  
  # Return results
  results <- list(
    species = species_matched,
    breakpoint_type = breakpoint_type,
    n_breakpoint_regions = nrow(breakpoints),
    n_chunks_total = nrow(chunks_sp),
    n_chunks_in_breakpoint = nrow(in_breakpoint),
    n_chunks_out_breakpoint = nrow(out_breakpoint),
    
    # Oligocentromere counts
    n_oligo_in_breakpoint = sum(in_breakpoint$has_oligocentromere),
    n_oligo_out_breakpoint = sum(out_breakpoint$has_oligocentromere),
    
    # Proportions
    prop_oligo_in_breakpoint = prop_in,
    prop_oligo_out_breakpoint = prop_out,
    fold_enrichment = prop_in / prop_out,
    
    # Contingency table
    contingency_table = contingency_table,
    
    # Statistical tests
    fisher_test = fisher_test,
    fisher_pvalue = fisher_test$p.value,
    fisher_odds_ratio = as.numeric(fisher_test$estimate),
    fisher_ci_lower = fisher_test$conf.int[1],
    fisher_ci_upper = fisher_test$conf.int[2],
    
    chisq_test = chisq_test,
    chisq_pvalue = if(!is.null(chisq_test)) chisq_test$p.value else NA,
    
    permutation_pvalue = perm_pvalue,
    permutation_diffs = permuted_diffs,
    observed_diff = observed_diff,
    
    # Detailed chunk data for further analysis
    chunks_data = chunks_sp
  )
  
  class(results) <- c("oligocentromere_enrichment_test", "list")
  return(results)
}


#' Print method for enrichment test results
#' @export
print.oligocentromere_enrichment_test <- function(x, ...) {
  cat("Oligocentromere Enrichment Test Results\n")
  cat("========================================\n")
  cat("Species:", x$species, "\n")
  cat("Breakpoint type:", x$breakpoint_type, "\n")
  cat("Fisher's exact test p-value:", format.pval(x$fisher_pvalue, digits = 3), "\n")
  cat("Odds ratio:", sprintf("%.2f", x$fisher_odds_ratio), "\n")
  cat("Fold enrichment:", sprintf("%.2f", x$fold_enrichment), "\n")
  
  if (x$fisher_pvalue < 0.05) {
    if (x$prop_oligo_in_breakpoint > x$prop_oligo_out_breakpoint) {
      cat("\n*** Significant ENRICHMENT of oligocentromeres in breakpoint regions ***\n")
    } else {
      cat("\n*** Significant DEPLETION of oligocentromeres in breakpoint regions ***\n")
    }
  } else {
    cat("\nNo significant enrichment or depletion\n")
  }
}


#' Plot permutation test results
#' @export
plot_permutation_test <- function(test_results) {
  if (!inherits(test_results, "oligocentromere_enrichment_test")) {
    stop("Input must be results from test_oligocentromere_enrichment()")
  }
  
  df <- data.frame(diff = test_results$permutation_diffs)
  
  p <- ggplot(df, aes(x = diff)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = test_results$observed_diff, 
               color = "red", linewidth = 1.5, linetype = "dashed") +
    annotate("text", 
             x = test_results$observed_diff, 
             y = Inf, 
             label = paste("Observed\ndifference\n", 
                           sprintf("%.4f", test_results$observed_diff)),
             vjust = 1.5, hjust = -0.1, color = "red", size = 4) +
    labs(title = "Permutation Test Results",
         subtitle = paste("Species:", test_results$species, 
                          "| Breakpoint type:", test_results$breakpoint_type),
         x = "Difference in proportions (permuted)",
         y = "Frequency",
         caption = paste("P-value:", format.pval(test_results$permutation_pvalue, digits = 3))) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(p)
}


#' Compare enrichment across multiple breakpoint types
#' @export
compare_enrichment_types <- function(breakpoints_file,
                                     chunks_file,
                                     species_name,
                                     chunk_size = 5000) {
  
  cat("Testing enrichment for different breakpoint types...\n\n")
  
  # Test each type
  results_fusion <- test_oligocentromere_enrichment(
    breakpoints_file, chunks_file, species_name, 
    breakpoint_type = "fusion", chunk_size = chunk_size, verbose = FALSE
  )
  
  results_fission <- test_oligocentromere_enrichment(
    breakpoints_file, chunks_file, species_name, 
    breakpoint_type = "fission", chunk_size = chunk_size, verbose = FALSE
  )
  
  results_both <- test_oligocentromere_enrichment(
    breakpoints_file, chunks_file, species_name, 
    breakpoint_type = "both", chunk_size = chunk_size, verbose = FALSE
  )
  
  # Create summary table
  summary_df <- data.frame(
    Breakpoint_Type = c("Fusion", "Fission", "Both"),
    N_Regions = c(results_fusion$n_breakpoint_regions,
                  results_fission$n_breakpoint_regions,
                  results_both$n_breakpoint_regions),
    Prop_with_Oligo = c(results_fusion$prop_oligo_in_breakpoint,
                        results_fission$prop_oligo_in_breakpoint,
                        results_both$prop_oligo_in_breakpoint),
    Fold_Enrichment = c(results_fusion$fold_enrichment,
                        results_fission$fold_enrichment,
                        results_both$fold_enrichment),
    Odds_Ratio = c(results_fusion$fisher_odds_ratio,
                   results_fission$fisher_odds_ratio,
                   results_both$fisher_odds_ratio),
    Fisher_Pvalue = c(results_fusion$fisher_pvalue,
                      results_fission$fisher_pvalue,
                      results_both$fisher_pvalue),
    Significant = c(results_fusion$fisher_pvalue < 0.05,
                    results_fission$fisher_pvalue < 0.05,
                    results_both$fisher_pvalue < 0.05)
  )
  
  cat("===========================================\n")
  cat("Comparison of Enrichment Across Types\n")
  cat("===========================================\n")
  cat("Species:", results_both$species, "\n")
  cat("Background proportion (outside breakpoints):", 
      sprintf("%.2f%%\n\n", results_both$prop_oligo_out_breakpoint * 100))
  
  print(summary_df, row.names = FALSE)
  
  cat("\n")
  
  return(list(
    fusion = results_fusion,
    fission = results_fission,
    both = results_both,
    summary = summary_df
  ))
}


############################################################################ performing for the species pairs

elata_plot <- plot_single_chromosome_cartoon(chunks, "Carex elata", "break_intervals_classified_elataquery.txt", other_species = "Carex nigra")

nigra_plot <- plot_single_chromosome_cartoon(chunks, "Carex nigra", "break_intervals_classified_nigraquery.txt", other_species = "Carex elata")

elata_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_elataquery.txt",
  chunks_file = chunks,
  species_name = "Carex elata", 
  breakpoint_type = "both"
)

nigra_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_elataquery.txt",
  chunks_file = chunks,
  species_name = "Carex elata", 
  breakpoint_type = "fission"
)

elata_fusions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_elataquery.txt",
  chunks_file = chunks,
  species_name = "Carex elata", 
  breakpoint_type = "fusion"
)

nigra_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_nigraquery.txt",
  chunks_file = chunks,
  species_name = "Carex nigra", 
  breakpoint_type = "both"
)

elata_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_nigraquery.txt",
  chunks_file = chunks,
  species_name = "Carex nigra", 
  breakpoint_type = "fission"
)

nigra_fusions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_nigraquery.txt",
  chunks_file = chunks,
  species_name = "Carex nigra", 
  breakpoint_type = "fusion"
)


## hirta and riparia

hirta_plot <- plot_single_chromosome_cartoon(chunks, "Carex hirta", "break_intervals_classified_hirtaquery.txt", other_species = "Carex riparia")

riparia_plot <- plot_single_chromosome_cartoon(chunks, "Carex riparia", "break_intervals_classified_ripariaquery.txt", other_species = "Carex hirta")

hirta_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_hirtaquery.txt",
  chunks_file = chunks,
  species_name = "Carex hirta",  
  breakpoint_type = "both"
)

riparia_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_hirtaquery.txt",
  chunks_file = chunks,
  species_name = "Carex hirta",
  breakpoint_type = "fission" 
)

hirta_fusions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_hirtaquery.txt",
  chunks_file = chunks,
  species_name = "Carex hirta",  
  breakpoint_type = "fusion"
)

riparia_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_ripariaquery.txt",
  chunks_file = chunks,
  species_name = "Carex riparia", 
  breakpoint_type = "both"
)

hirta_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_ripariaquery.txt",
  chunks_file = chunks,
  species_name = "Carex riparia", 
  breakpoint_type = "fission"
)

riparia_fusions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_ripariaquery.txt",
  chunks_file = chunks,
  species_name = "Carex riparia", 
  breakpoint_type = "fusion"
)

## distans and extensa

distans_plot <- plot_single_chromosome_cartoon(chunks, "Carex distans", "break_intervals_classified_distansquery.txt", other_species = "Carex extensa")

extensa_plot <- plot_single_chromosome_cartoon(chunks, "Carex extensa", "break_intervals_classified_extensaquery.txt", other_species = "Carex distans")

distans_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_distansquery.txt",
  chunks_file = chunks,
  species_name = "Carex distans", 
  breakpoint_type = "both"
)

extensa_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_distansquery.txt",
  chunks_file = chunks,
  species_name = "Carex distans", 
  breakpoint_type = "fission"
)

distans_fusions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_distansquery.txt",
  chunks_file = chunks,
  species_name = "Carex distans", 
  breakpoint_type = "fusion"
)

extensa_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_extensaquery.txt",
  chunks_file = chunks,
  species_name = "Carex extensa", 
  breakpoint_type = "both"
)

distans_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_extensaquery.txt",
  chunks_file = chunks,
  species_name = "Carex extensa", 
  breakpoint_type = "fission"
)

extensa_fusions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_extensaquery.txt",
  chunks_file = chunks,
  species_name = "Carex extensa", 
  breakpoint_type = "fusion"
)

## Carex pendula & Carex sylvatica

pendula_plot <- plot_single_chromosome_cartoon(chunks, "Carex pendula", "break_intervals_classified_pendulaquery.txt", other_species = "Carex sylvatica")

sylvatica_plot <- plot_single_chromosome_cartoon(chunks, "Carex sylvatica", "break_intervals_classified_sylvaticaquery.txt", other_species = "Carex pendula")

pendula_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_pendulaquery.txt",
  chunks_file = chunks,
  species_name = "Carex pendula", 
  breakpoint_type = "both"
)

sylvatica_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_pendulaquery.txt",
  chunks_file = chunks,
  species_name = "Carex pendula", 
  breakpoint_type = "fission"
)

pendula_fusion_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_pendulaquery.txt",
  chunks_file = chunks,
  species_name = "Carex pendula", 
  breakpoint_type = "fusion"
)

sylvatica_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_sylvaticaquery.txt",
  chunks_file = chunks,
  species_name = "Carex sylvatica", 
  breakpoint_type = "both"
)

pendula_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_sylvaticaquery.txt",
  chunks_file = chunks,
  species_name = "Carex sylvatica", 
  breakpoint_type = "fission"
)

sylvatica_fusions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_sylvaticaquery.txt",
  chunks_file = chunks,
  species_name = "Carex sylvatica", 
  breakpoint_type = "fusion"
)

# Carex divulsa and Carex spicata

divulsa_plot <- plot_single_chromosome_cartoon(chunks, "Carex divulsa", "break_intervals_classified_divulsaquery.txt", other_species = "Carex spicata")

spicata_plot <- plot_single_chromosome_cartoon(chunks, "Carex spicata", "break_intervals_classified_spicataquery.txt", other_species = "Carex divulsa")

divulsa_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_divulsaquery.txt",
  chunks_file = chunks,
  species_name = "Carex divulsa", 
  breakpoint_type = "both"
)

spicata_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_divulsaquery.txt",
  chunks_file = chunks,
  species_name = "Carex divulsa", 
  breakpoint_type = "fission"
)

divulsa_fusion_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_divulsaquery.txt",
  chunks_file = chunks,
  species_name = "Carex divulsa", 
  breakpoint_type = "fusion"
)

spicata_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_spicataquery.txt",
  chunks_file = chunks,
  species_name = "Carex spicata", 
  breakpoint_type = "both"
)

divulsa_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_spicataquery.txt",
  chunks_file = chunks,
  species_name = "Carex spicata", 
  breakpoint_type = "fission"
)

spicata_fusions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_spicataquery.txt",
  chunks_file = chunks,
  species_name = "Carex spicata", 
  breakpoint_type = "fusion"
)

## Carex littledalei and Carex myosuroides

littledalei_plot <- plot_single_chromosome_cartoon(chunks, "Carex littledalei", "break_intervals_classified_littledaleiquery.txt", other_species = "Carex myosuroides")

myosuroides_plot <- plot_single_chromosome_cartoon(chunks, "Carex myosuroides", "break_intervals_classified_myosuroidesquery.txt", other_species = "Carex littledalei")

littledalei_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_littledaleiquery.txt",
  chunks_file = chunks,
  species_name = "Carex littledalei", 
  breakpoint_type = "both"
)

myosuroides_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_littledaleiquery.txt",
  chunks_file = chunks,
  species_name = "Carex littledalei", 
  breakpoint_type = "fission"
)

littledalei_fusions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_littledaleiquery.txt",
  chunks_file = chunks,
  species_name = "Carex littledalei", 
  breakpoint_type = "fusion"
)

myosuroides_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_myosuroidesquery.txt",
  chunks_file = chunks,
  species_name = "Carex myosuroides", 
  breakpoint_type = "both"
)

littledalei_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_myosuroidesquery.txt",
  chunks_file = chunks,
  species_name = "Carex myosuroides", 
  breakpoint_type = "fission"
)

myosuroides_fusions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_myosuroidesquery.txt",
  chunks_file = chunks,
  species_name = "Carex myosuroides", 
  breakpoint_type = "fusion"
)

## Eriophorum angustifolium and Eriophorum vaginatum

angustifolium_plot <- plot_single_chromosome_cartoon(chunks, "Eriophorum angustifolium", "break_intervals_classified_angustifoliumquery.txt", other_species = "Eriophorum angustifolium")

vaginatum_plot <- plot_single_chromosome_cartoon(chunks, "Eriophorum vaginatum", "break_intervals_classified_vaginatumquery.txt", other_species = "Eriophorum vaginatum")

angustifolium_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_angustifoliumquery.txt",
  chunks_file = chunks,
  species_name = "Eriophorum angustifolium", 
  breakpoint_type = "both"
)

vaginatum_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_angustifoliumquery.txt",
  chunks_file = chunks,
  species_name = "Eriophorum angustifolium", 
  breakpoint_type = "fission"
)

angustifolium_fusions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_angustifoliumquery.txt",
  chunks_file = chunks,
  species_name = "Eriophorum angustifolium", 
  breakpoint_type = "fusion"
)

vaginatum_regions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_vaginatumquery.txt",
  chunks_file = chunks,
  species_name = "Eriophorum vaginatum", 
  breakpoint_type = "both"
)

angustifolium_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_vaginatumquery.txt",
  chunks_file = chunks,
  species_name = "Eriophorum vaginatum", 
  breakpoint_type = "fission"
)

vaginatum_fusions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_vaginatumquery.txt",
  chunks_file = chunks,
  species_name = "Eriophorum vaginatum", 
  breakpoint_type = "fusion"
)

## Rhynchospora breviuscula and Rhynchospora tenuis

breviuscula_plot <- plot_single_chromosome_cartoon(chunks, "Rhynchospora breviuscula", "break_intervals_classified_breviusculaquery.txt", other_species = "Rhynchospora tenuis")

tenuis_plot <- plot_single_chromosome_cartoon(chunks, "Rhynchospora tenuis", "break_intervals_classified_tenuisquery.txt", other_species = "Rhynchospora breviuscula")

breviuscula_fissions_test <- test_oligocentromere_enrichment(
  breakpoints_file = "break_intervals_classified_tenuisquery.txt",
  chunks_file = chunks,
  species_name = "Rhynchospora tenuis", 
  breakpoint_type = "fission"
)




########################################################################################### Forest plot for Carex

desired_order_breakpoints <- c("Carex elata",
                               "Carex nigra",
                               "Carex hirta",
                               "Carex riparia",
                               "Carex distans",
                               "Carex extensa",
                               "Carex pendula",
                               "Carex sylvatica",
                               "Carex divulsa",
                               "Carex spicata",
                               "Carex littledalei",
                               "Carex myosuroides")

# Read the data
breakpoint_data <- read_tsv("Carex_breakpoints_data.tsv")

# Clean column names
colnames(breakpoint_data) <- gsub(" ", "_", colnames(breakpoint_data))

# Prepare fission data
fission_data <- breakpoint_data %>%
  select(Species,
         n_regions = Fission_regions,
         share_in = Share_within_fissions,
         share_out = Share_outside_fissions,
         OR = Odds_ratio...12,  # Second OR column (for fissions)
         CI_lower = Lower_bound...13,
         CI_upper = Upper_bound...14,
         P = `P-value...15`) %>%
  # Remove species with no fission data
  filter(!is.na(OR)) %>%
  # Parse p-values (handle "<2E-16" format)
  mutate(
    P_numeric = case_when(
      grepl("<", P) ~ 0,
      TRUE ~ as.numeric(P)
    ),
    # Create significance categories
    sig_level = case_when(
      P_numeric < 0.001 ~ "P < 0.001",
      P_numeric < 0.01 ~ "P < 0.01",
      P_numeric < 0.05 ~ "P < 0.05",
      TRUE ~ "P ≥ 0.05"
    ),
    sig_level = factor(sig_level, levels = c("P < 0.001", "P < 0.01", "P < 0.05", "P ≥ 0.05")),
    # Flag studies with few regions
    low_power = n_regions < 3
  ) %>%
  # Order species by OR (smallest to largest)
  arrange(OR) %>%
  mutate(Species = factor(Species, levels = Species))

fission_data$Species <- factor(fission_data$Species, levels = rev(desired_order_breakpoints))
fission_data$Species_label <- gsub("_", " ", fission_data$Species)
fission_data$Species_label <- paste0("italic('", fission_data$Species_label, "')")


# Create the forest plot
p <- ggplot(fission_data, aes(y = Species)) +
  
  scale_y_discrete(labels = function(x) {
    x <- gsub("_", " ", x)
    parse(text = paste0("italic('", x, "')"))
  }) +
  
  # Reference line at OR = 1 (no enrichment)
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  
  # Confidence interval error bars
  geom_linerange(aes(xmin = CI_lower, xmax = CI_upper, 
                     #color = sig_level
                     ),
                 linewidth = 0.8, alpha = 0.7) +
  
  # Point estimates (odds ratios)
  geom_point(aes(x = OR, 
                 #color = sig_level
                 ), size = 3, alpha = 0.9) +
  
  # Color scheme for significance levels
  #scale_color_manual(
    #values = c("P < 0.001" = "#E41A1C",   # Red
               #"P < 0.01" = "#FF7F00",     # Orange
               #"P < 0.05" = "#FFD700",     # Yellow/Gold
               #"P ≥ 0.05" = "black"),     # Gray
    #name = "Significance"
  #) +
  
  # Log scale for x-axis
  scale_x_log10(
    breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16),
    labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8", "16")
  ) +
  
  # Labels
  labs(
    title = "Oligocentromere enrichment in fission breakpoint regions",
    x = "Odds ratio (95% CI, log scale)",
    y = NULL
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray95", linewidth = 0.3),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 9),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30"),
    plot.caption = element_text(size = 8, color = "gray50")
  )

# View the plot
print(p)

# Save the plot
ggsave("forest_plot_fissions.png", p, width = 8, height = 8, dpi = 300, bg = "white")
ggsave("forest_plot_fissions.pdf", p, width = 8, height = 8)

#### Fusions plot

# Prepare fusion data
fusion_data <- breakpoint_data %>%
  select(Species,
         n_regions = Fusion_regions,
         share_in = Share_within_fusions,
         share_out = Share_outside_fusions,
         OR = Odds_ratio...5,  # First OR column (for fusions)
         CI_lower = Lower_bound...6,
         CI_upper = Upper_bound...7,
         P = `P-value...8`) %>%
  # Remove species with no fusion data
  filter(!is.na(OR)) %>%
  # Parse p-values (handle "<2E-16" format)
  mutate(
    P_numeric = case_when(
      grepl("<", P) ~ 0,
      TRUE ~ as.numeric(P)
    ),
    # Create significance categories
    sig_level = case_when(
      P_numeric < 0.001 ~ "P < 0.001",
      P_numeric < 0.01 ~ "P < 0.01",
      P_numeric < 0.05 ~ "P < 0.05",
      TRUE ~ "P ≥ 0.05"
    ),
    sig_level = factor(sig_level, levels = c("P < 0.001", "P < 0.01", "P < 0.05", "P ≥ 0.05")),
    # Flag studies with few regions
    low_power = n_regions < 3
  ) %>%
  # Order species by OR (smallest to largest)
  arrange(OR) %>%
  mutate(Species = factor(Species, levels = Species))

fusion_data$Species <- factor(fusion_data$Species, levels = rev(desired_order_breakpoints))
fusion_data$Species_label <- gsub("_", " ", fusion_data$Species)
fusion_data$Species_label <- paste0("italic('", fusion_data$Species_label, "')")


# Create the forest plot
p <- ggplot(fusion_data, aes(y = Species)) +
  
  scale_y_discrete(labels = function(x) {
    x <- gsub("_", " ", x)
    parse(text = paste0("italic('", x, "')"))
  }) +
  
  # Reference line at OR = 1 (no enrichment)
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  
  # Confidence interval error bars
  geom_linerange(aes(xmin = CI_lower, xmax = CI_upper, 
                     #color = sig_level
  ),
  linewidth = 0.8, alpha = 0.7) +
  
  # Point estimates (odds ratios)
  geom_point(aes(x = OR, 
                 #color = sig_level
  ), size = 3, alpha = 0.9) +
  
  # Color scheme for significance levels
  #scale_color_manual(
  #values = c("P < 0.001" = "#E41A1C",   # Red
  #"P < 0.01" = "#FF7F00",     # Orange
  #"P < 0.05" = "#FFD700",     # Yellow/Gold
  #"P ≥ 0.05" = "black"),     # Gray
  #name = "Significance"
  #) +
  
  #Log scale for x-axis
  scale_x_log10(
  breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16),
  labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8", "16")
  ) +
  
  # Labels
  labs(
    title = "Oligocentromere enrichment in fusion breakpoint regions",
    x = "Odds ratio (95% CI, log scale)",
    y = NULL
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray95", linewidth = 0.3),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 9),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30"),
    plot.caption = element_text(size = 8, color = "gray50")
  )

# View the plot
print(p)

# Save the plot
ggsave("forest_plot_fusions.png", p, width = 8, height = 8, dpi = 300, bg = "white")
ggsave("forest_plot_fusions.pdf", p, width = 8, height = 8)


################################################################################## Chromosome cartoon plots


### Helper functions for the below

# Helper function to normalise species names
normalise_species_name <- function(species_name) {
  # Remove spaces and replace with underscores
  species_name <- gsub(" ", "_", species_name)
  # Ensure format is Genus_species
  if (!grepl("_", species_name)) {
    stop("Species name must be in format 'Genus_species' or 'Genus species'")
  }
  return(species_name)
}

# Helper function to match species name to what's in the data
match_species_name <- function(species_input, chunks_data) {
  # Normalize the input
  normalised_input <- normalise_species_name(species_input)
  
  # Get unique species names from chunks data
  available_species <- unique(chunks_data$species)
  
  # Try exact match first
  if (normalised_input %in% available_species) {
    return(normalised_input)
  }
  
  # If not found, try to match by removing/adding underscores
  # Convert input to both formats
  with_underscore <- gsub(" ", "_", species_input)
  without_underscore <- gsub("_", " ", species_input)
  
  # Check both versions
  if (with_underscore %in% available_species) {
    return(with_underscore)
  }
  if (without_underscore %in% available_species) {
    return(without_underscore)
  }
  
  # If still not found, provide helpful error
  stop(paste0("Species '", species_input, "' not found in chunks data.\n",
              "Available species: ", paste(available_species, collapse = ", ")))
}


#' Merge consecutive chunks into contiguous regions
#'
#' Takes chunks with features and merges consecutive chunks of the same type
#' into single regions to minimize the number of rectangles to plot
#'
#' @param chunks_df Data frame with chunk information
#' @return Data frame with merged regions
merge_consecutive_chunks <- function(chunks_df) {
  if (nrow(chunks_df) == 0) {
    return(data.frame(
      chromosome = character(),
      region_start = numeric(),
      region_end = numeric(),
      color = character()
    ))
  }
  
  # Sort by chromosome and position
  chunks_sorted <- chunks_df %>%
    arrange(chromosome, chunk_start)
  
  # Initialize results list
  regions <- list()
  
  # Process each chromosome separately
  for (chr in unique(chunks_sorted$chromosome)) {
    chr_chunks <- chunks_sorted %>% filter(chromosome == chr)
    
    if (nrow(chr_chunks) == 0) next
    
    # Start first region
    current_start <- chr_chunks$chunk_start[1]
    current_end <- chr_chunks$chunk_end[1]
    current_color <- chr_chunks$color[1]
    current_chr <- chr
    
    for (i in 2:nrow(chr_chunks)) {
      # Check if this chunk is consecutive and same color
      if (chr_chunks$chunk_start[i] == current_end && 
          chr_chunks$color[i] == current_color) {
        # Extend current region
        current_end <- chr_chunks$chunk_end[i]
      } else {
        # Save current region and start new one
        regions[[length(regions) + 1]] <- data.frame(
          chromosome = current_chr,
          region_start = current_start,
          region_end = current_end,
          color = current_color
        )
        current_start <- chr_chunks$chunk_start[i]
        current_end <- chr_chunks$chunk_end[i]
        current_color <- chr_chunks$color[i]
      }
    }
    
    # Don't forget the last region
    regions[[length(regions) + 1]] <- data.frame(
      chromosome = current_chr,
      region_start = current_start,
      region_end = current_end,
      color = current_color
    )
  }
  
  # Combine all regions
  if (length(regions) > 0) {
    merged <- bind_rows(regions)
    return(merged)
  } else {
    return(data.frame(
      chromosome = character(),
      region_start = numeric(),
      region_end = numeric(),
      color = character()
    ))
  }
}


#' Create chromosome ideogram with only feature regions plotted
#'
#' This function creates simple chromosome ideograms showing only repeat arrays
#' and breakpoint regions. The chromosomes are drawn as solid bars with features
#' overlaid - NO individual chunks are plotted.
#'
#' @param chunks_file Either a file path to chunks data (TSV) OR a pre-loaded data frame
#' @param species_name Name of the species
#' @param breakpoints_file Path to breakpoints file (TSV with: chr, break_start, break_end, class)
#' @param other_species_name Name of comparison species (for subtitle)
#' @param chunk_size Size of chunks in bp (default: 5000)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 10)
#' @param chr_height Height of chromosome bars (default: 0.6)
#' @param fusion_color Color for fusion regions (default: "#6BAED6")
#' @param fission_color Color for fission regions (default: "#FC8D62")
#' @param repeat_color Color for repeat arrays (default: "black")
#' @param background_color Color for chromosome background (default: "#E5E5E5")
#' @param save_plot Whether to save plot to file (default: TRUE)
#' @return A ggplot object
#' @export
plot_chromosome_ideogram <- function(chunks_file,
                                     species_name,
                                     breakpoints_file,
                                     other_species_name = NULL,
                                     chunk_size = 5000,
                                     width = 12,
                                     height = 10,
                                     chr_height = 0.6,
                                     fusion_color = "#6BAED6",
                                     fission_color = "#FC8D62",
                                     repeat_color = "black",
                                     background_color = "#E5E5E5",
                                     save_plot = TRUE) {
  
  # Read the chunks data (or use provided data frame)
  if (is.character(chunks_file)) {
    chunks <- read_tsv(chunks_file, show_col_types = FALSE)
  } else if (is.data.frame(chunks_file)) {
    chunks <- chunks_file
  } else {
    stop("chunks_file must be either a file path (character) or a data frame")
  }
  
  # Match species name
  species_matched <- match_species_name(species_name, chunks)
  cat("Matched species name:", species_matched, "\n")
  
  # Read breakpoints
  breakpoints <- read_tsv(breakpoints_file, show_col_types = FALSE)
  
  # Extract species names
  extract_epithet <- function(sp_name) {
    if (grepl("_", sp_name)) {
      return(str_match(sp_name, "_(.+)$")[,2])
    } else if (grepl(" ", sp_name)) {
      return(str_match(sp_name, " (.+)$")[,2])
    } else {
      return(sp_name)
    }
  }
  
  extract_genus <- function(sp_name) {
    if (grepl("_", sp_name)) {
      return(str_match(sp_name, "^([^_]+)_")[,2])
    } else if (grepl(" ", sp_name)) {
      return(str_match(sp_name, "^([^ ]+) ")[,2])
    } else {
      return("")
    }
  }
  
  current_epithet <- extract_epithet(species_matched)
  genus <- extract_genus(species_matched)
  
  if (!is.null(other_species_name)) {
    other_epithet <- extract_epithet(other_species_name)
  } else {
    other_epithet <- "[other species]"
  }
  
  # Filter chunks for this species
  chunks_sp <- chunks %>% 
    filter(species == species_matched)
  
  # Calculate chunk positions
  chunks_sp <- chunks_sp %>%
    mutate(chunk_start = (chunk_id - 1) * chunk_size,
           chunk_end = chunk_id * chunk_size)
  
  # Assign breakpoint classes
  chunks_sp$breakpoint_class <- NA
  
  for(i in 1:nrow(breakpoints)){
    chr <- breakpoints$chr[i]
    break_start <- breakpoints$break_start[i]
    break_end <- breakpoints$break_end[i]
    class <- breakpoints$class[i]
    
    overlapping <- which(
      chunks_sp$chromosome == chr & 
        chunks_sp$chunk_end >= break_start & 
        chunks_sp$chunk_start <= break_end
    )
    
    chunks_sp$breakpoint_class[overlapping] <- class
  }
  
  # Print summary
  cat("\nSummary:\n")
  cat("Total chunks:", nrow(chunks_sp), "\n")
  cat("Chunks with repeats:", sum(!is.na(chunks_sp$repeat_classes)), "\n")
  cat("Chunks with fusions:", sum(grepl("fusion", chunks_sp$breakpoint_class, ignore.case = TRUE), na.rm = TRUE), "\n")
  cat("Chunks with fissions:", sum(grepl("fission", chunks_sp$breakpoint_class, ignore.case = TRUE), na.rm = TRUE), "\n")
  
  # Get chromosome lengths (from first to last chunk position)
  chr_sizes <- chunks_sp %>%
    group_by(chromosome) %>%
    summarise(
      chr_start = 0,
      chr_end = max(chunk_end)
    ) %>%
    arrange(desc(chr_end)) %>%
    mutate(chr_num = n():1)
  
  cat("\nChromosomes:", nrow(chr_sizes), "\n")
  
  # Filter to only chunks with features
  chunks_features <- chunks_sp %>%
    filter(!is.na(repeat_classes) | !is.na(breakpoint_class)) %>%
    mutate(
      color = case_when(
        !is.na(repeat_classes) ~ repeat_color,
        grepl("fission", breakpoint_class, ignore.case = TRUE) ~ fission_color,
        grepl("fusion", breakpoint_class, ignore.case = TRUE) ~ fusion_color,
        TRUE ~ background_color
      )
    )
  
  cat("Chunks with features:", nrow(chunks_features), "\n")
  
  # Merge consecutive chunks into regions
  cat("Merging consecutive chunks into regions...\n")
  regions <- merge_consecutive_chunks(chunks_features)
  
  cat("Created", nrow(regions), "feature regions\n")
  
  # Add chromosome numbers and convert to Mb
  regions <- regions %>%
    left_join(chr_sizes %>% select(chromosome, chr_num), by = "chromosome") %>%
    mutate(
      region_start_mb = region_start / 1e6,
      region_end_mb = region_end / 1e6
    )
  
  # Convert chromosome sizes to Mb
  chr_sizes <- chr_sizes %>%
    mutate(
      chr_start_mb = chr_start / 1e6,
      chr_end_mb = chr_end / 1e6
    )
  
  # Separate regions by type
  regions_fusion <- regions %>% filter(color == fusion_color)
  regions_fission <- regions %>% filter(color == fission_color)
  regions_repeat <- regions %>% filter(color == repeat_color)
  
  cat("\nFinal counts:\n")
  cat("  Chromosomes:", nrow(chr_sizes), "\n")
  cat("  Fusion regions:", nrow(regions_fusion), "\n")
  cat("  Fission regions:", nrow(regions_fission), "\n")
  cat("  Repeat regions:", nrow(regions_repeat), "\n")
  cat("  TOTAL OBJECTS TO PLOT:", nrow(chr_sizes) + nrow(regions), "\n")
  
  # Create the plot
  p <- ggplot() +
    # Chromosome backgrounds (one rectangle per chromosome)
    geom_rect(data = chr_sizes,
              aes(xmin = chr_start_mb, xmax = chr_end_mb,
                  ymin = chr_num - chr_height/2, ymax = chr_num + chr_height/2),
              fill = background_color,
              color = NA) +
    # Fusion regions
    geom_rect(data = regions_fusion,
              aes(xmin = region_start_mb, xmax = region_end_mb,
                  ymin = chr_num - chr_height/2, ymax = chr_num + chr_height/2),
              fill = fusion_color,
              color = NA) +
    # Fission regions
    geom_rect(data = regions_fission,
              aes(xmin = region_start_mb, xmax = region_end_mb,
                  ymin = chr_num - chr_height/2, ymax = chr_num + chr_height/2),
              fill = fission_color,
              color = NA) +
    # Repeat regions (on top)
    geom_rect(data = regions_repeat,
              aes(xmin = region_start_mb, xmax = region_end_mb,
                  ymin = chr_num - chr_height/2, ymax = chr_num + chr_height/2),
              fill = repeat_color,
              color = NA) +
    scale_y_continuous(breaks = chr_sizes$chr_num,
                       labels = chr_sizes$chromosome) +
    labs(x = "Position (Mb)", 
         y = "Chromosome",
         title = bquote(italic(.(genus))~.(current_epithet)~"chromosomes"),
         subtitle = paste0("Black: Repeat arrays | Blue: fusions in ", genus, " ", current_epithet, 
                           " | Red: fissions in ", genus, " ", other_epithet)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10))
  
  # Save plot if requested
  if (save_plot) {
    output_prefix <- paste0(tolower(current_epithet), "_chromosome_ideogram")
    ggsave(paste0(output_prefix, ".pdf"), p, width = width, height = height)
    ggsave(paste0(output_prefix, ".png"), p, width = width, height = height, dpi = 300, bg = "white")
    ggsave(paste0(output_prefix, ".svg"), p, width = width, height = height, bg = "white")
    cat("\nPlot saved as:", output_prefix, ".pdf/.png/.svg\n")
  }
  
  return(invisible(p))
}

riparia_cartoon <- plot_chromosome_ideogram(chunks, "Carex riparia", "break_intervals_classified_ripariaquery.txt",
                                            "Carex hirta")

print(riparia_cartoon)


#################################################################################
##################### Enrichment of genes #######################################
#################################################################################

# Helper function to normalize species names (reuse from cartoon function)
normalise_species_name <- function(species_name) {
  species_name <- gsub(" ", "_", species_name)
  if (!grepl("_", species_name)) {
    stop("Species name must be in format 'Genus_species' or 'Genus species'")
  }
  return(species_name)
}

# Helper function to match species name to data
match_species_name <- function(species_input, chunks_data) {
  normalised_input <- normalise_species_name(species_input)
  available_species <- unique(chunks_data$species)
  
  if (normalised_input %in% available_species) {
    return(normalised_input)
  }
  
  with_underscore <- gsub(" ", "_", species_input)
  without_underscore <- gsub("_", " ", species_input)
  
  if (with_underscore %in% available_species) {
    return(with_underscore)
  }
  if (without_underscore %in% available_species) {
    return(without_underscore)
  }
  
  stop(paste0("Species '", species_input, "' not found in chunks data.\n",
              "Available species: ", paste(available_species, collapse = ", ")))
}

#' Test for gene enrichment in breakpoint regions
#'
#' @param breakpoints_file Path to breakpoints file (TSV with: chr, break_start, break_end, class)
#' @param chunks_file Either a file path to chunks data (TSV) OR a pre-loaded data frame
#' @param species_name Name of the species (e.g., "Carex_distans", "Carex distans", or "distans")
#' @param breakpoint_type Which breakpoints to test: "fusion", "fission", or "both" (default: "both")
#' @param chunk_size Size of chunks in bp (default: 5000)
#' @param verbose Print detailed output (default: TRUE)
#' @return A list with test results and summary statistics
#' @export
test_CDS_enrichment <- function(breakpoints_file,
                                            chunks_file,
                                            species_name,
                                            breakpoint_type = "both",
                                            chunk_size = 5000,
                                            verbose = TRUE) {
  
  # Validate breakpoint_type
  if (!breakpoint_type %in% c("fusion", "fission", "both")) {
    stop("breakpoint_type must be one of: 'fusion', 'fission', or 'both'")
  }
  
  # Read the chunks data (or use provided data frame)
  if (is.character(chunks_file)) {
    chunks <- read_tsv(chunks_file, show_col_types = FALSE)
  } else if (is.data.frame(chunks_file)) {
    chunks <- chunks_file
  } else {
    stop("chunks_file must be either a file path (character) or a data frame")
  }
  
  # Match species name
  species_matched <- match_species_name(species_name, chunks)
  
  if (verbose) {
    cat("===========================================\n")
    cat("Oligocentromere Enrichment Test\n")
    cat("===========================================\n")
    cat("Species:", species_matched, "\n")
    cat("Testing breakpoint type:", breakpoint_type, "\n\n")
  }
  
  # Read breakpoints
  breakpoints <- read_tsv(breakpoints_file, show_col_types = FALSE)
  
  # Filter breakpoints by type if needed
  if (breakpoint_type == "fusion") {
    breakpoints <- breakpoints %>%
      filter(grepl("fusion", class, ignore.case = TRUE))
  } else if (breakpoint_type == "fission") {
    breakpoints <- breakpoints %>%
      filter(grepl("fission", class, ignore.case = TRUE))
  }
  # If "both", use all breakpoints
  
  if (nrow(breakpoints) == 0) {
    stop("No breakpoints found for the specified type")
  }
  
  if (verbose) {
    cat("Number of breakpoint regions:", nrow(breakpoints), "\n")
  }
  
  # Filter chunks for this species
  chunks_sp <- chunks %>%
    filter(species == species_matched)
  
  # Calculate chunk start and end positions
  chunks_sp <- chunks_sp %>%
    mutate(chunk_start = (chunk_id - 1) * chunk_size,
           chunk_end = chunk_id * chunk_size)
  
  # Classify each chunk as in breakpoint region or not
  chunks_sp$in_breakpoint <- FALSE
  
  for (i in 1:nrow(breakpoints)) {
    chr <- breakpoints$chr[i]
    break_start <- breakpoints$break_start[i]
    break_end <- breakpoints$break_end[i]
    
    # Find chunks that overlap this breakpoint
    overlapping <- which(
      chunks_sp$chromosome == chr & 
        chunks_sp$chunk_end >= break_start & 
        chunks_sp$chunk_start <= break_end
    )
    
    chunks_sp$in_breakpoint[overlapping] <- TRUE
  }
  
  # Identify chunks with genes (non-NA repeat_classes)
  chunks_sp <- chunks_sp %>%
    mutate(has_CDS = CDS_overlap == 1)
  
  # Create contingency table
  # Rows: in/out of breakpoint region
  # Cols: has/doesn't have oligocentromere
  contingency_table <- table(
    Breakpoint = chunks_sp$in_breakpoint,
    CDS = chunks_sp$has_CDS
  )
  
  if (verbose) {
    cat("\n--- Contingency Table ---\n")
    print(contingency_table)
    cat("\n")
  }
  
  # Calculate proportions
  in_breakpoint <- chunks_sp %>% filter(in_breakpoint)
  out_breakpoint <- chunks_sp %>% filter(!in_breakpoint)
  
  prop_in <- mean(in_breakpoint$has_CDS)
  prop_out <- mean(out_breakpoint$has_CDS)
  
  if (verbose) {
    cat("--- Proportions ---\n")
    cat("Chunks in breakpoint regions with CDS:", 
        sum(in_breakpoint$has_CDS), "/", nrow(in_breakpoint),
        sprintf(" (%.2f%%)\n", prop_in * 100))
    cat("Chunks outside breakpoint regions with CDS:", 
        sum(out_breakpoint$has_CDS), "/", nrow(out_breakpoint),
        sprintf(" (%.2f%%)\n", prop_out * 100))
    cat("Fold enrichment:", sprintf("%.2f\n", prop_in / prop_out))
    cat("\n")
  }
  
  # Perform Fisher's Exact Test
  fisher_test <- fisher.test(contingency_table)
  
  # Perform Chi-square test (if cell counts are sufficient)
  chisq_test <- NULL
  chisq_warning <- NULL
  tryCatch({
    chisq_test <- chisq.test(contingency_table)
  }, warning = function(w) {
    chisq_warning <<- w$message
  })
  
  if (verbose) {
    cat("--- Fisher's Exact Test ---\n")
    cat("Odds Ratio:", sprintf("%.3f\n", fisher_test$estimate))
    cat("95% CI:", sprintf("[%.3f, %.3f]\n", 
                           fisher_test$conf.int[1], fisher_test$conf.int[2]))
    cat("P-value:", format.pval(fisher_test$p.value, digits = 3), "\n")
    
    if (fisher_test$p.value < 0.001) {
      cat("*** Highly significant enrichment\n")
    } else if (fisher_test$p.value < 0.01) {
      cat("** Significant enrichment\n")
    } else if (fisher_test$p.value < 0.05) {
      cat("* Significant enrichment\n")
    } else {
      cat("Not significant\n")
    }
    cat("\n")
    
    if (!is.null(chisq_test)) {
      cat("--- Chi-square Test ---\n")
      cat("Chi-square statistic:", sprintf("%.3f\n", chisq_test$statistic))
      cat("P-value:", format.pval(chisq_test$p.value, digits = 3), "\n\n")
    } else if (!is.null(chisq_warning)) {
      cat("--- Chi-square Test ---\n")
      cat("Warning:", chisq_warning, "\n")
      cat("(Use Fisher's test for small sample sizes)\n\n")
    }
  }
  
  # Permutation test for additional validation
  if (verbose) {
    cat("--- Permutation Test (10,000 iterations) ---\n")
    cat("Running permutations...")
  }
  
  n_permutations <- 10000
  observed_diff <- prop_in - prop_out
  
  permuted_diffs <- replicate(n_permutations, {
    # Randomly shuffle the breakpoint labels
    shuffled <- sample(chunks_sp$in_breakpoint)
    prop_in_perm <- mean(chunks_sp$has_CDS[shuffled])
    prop_out_perm <- mean(chunks_sp$has_CDS[!shuffled])
    prop_in_perm - prop_out_perm
  })
  
  # Two-tailed p-value
  perm_pvalue <- mean(abs(permuted_diffs) >= abs(observed_diff))
  
  if (verbose) {
    cat(" Done!\n")
    cat("Observed difference in proportions:", sprintf("%.4f\n", observed_diff))
    cat("Permutation p-value:", format.pval(perm_pvalue, digits = 3), "\n")
    cat("\n")
  }
  
  # Summary interpretation
  if (verbose) {
    cat("===========================================\n")
    cat("SUMMARY\n")
    cat("===========================================\n")
    
    if (fisher_test$p.value < 0.05 && prop_in > prop_out) {
      cat("CONCLUSION: CDSs are SIGNIFICANTLY ENRICHED\n")
      cat("in", breakpoint_type, "breakpoint regions.\n")
    } else if (fisher_test$p.value < 0.05 && prop_in < prop_out) {
      cat("CONCLUSION: CDSs are SIGNIFICANTLY DEPLETED\n")
      cat("in", breakpoint_type, "breakpoint regions.\n")
    } else {
      cat("CONCLUSION: No significant enrichment or depletion of\n")
      cat("CDSs in", breakpoint_type, "breakpoint regions.\n")
    }
    cat("\n")
  }
  
  # Return results
  results <- list(
    species = species_matched,
    breakpoint_type = breakpoint_type,
    n_breakpoint_regions = nrow(breakpoints),
    n_chunks_total = nrow(chunks_sp),
    n_chunks_in_breakpoint = nrow(in_breakpoint),
    n_chunks_out_breakpoint = nrow(out_breakpoint),
    
    # Oligocentromere counts
    n_CDS_in_breakpoint = sum(in_breakpoint$has_CDS),
    n_CDS_out_breakpoint = sum(out_breakpoint$has_CDS),
    
    # Proportions
    prop_CDS_in_breakpoint = prop_in,
    prop_CDS_out_breakpoint = prop_out,
    fold_enrichment = prop_in / prop_out,
    
    # Contingency table
    contingency_table = contingency_table,
    
    # Statistical tests
    fisher_test = fisher_test,
    fisher_pvalue = fisher_test$p.value,
    fisher_odds_ratio = as.numeric(fisher_test$estimate),
    fisher_ci_lower = fisher_test$conf.int[1],
    fisher_ci_upper = fisher_test$conf.int[2],
    
    chisq_test = chisq_test,
    chisq_pvalue = if(!is.null(chisq_test)) chisq_test$p.value else NA,
    
    permutation_pvalue = perm_pvalue,
    permutation_diffs = permuted_diffs,
    observed_diff = observed_diff,
    
    # Detailed chunk data for further analysis
    chunks_data = chunks_sp
  )
  
  class(results) <- c("CDS_enrichment_test", "list")
  return(results)
}

## testing

elata_fusions_gene_test <- test_CDS_enrichment("break_intervals_classified_elataquery.txt", chunks, 
                                               "Carex elata", breakpoint_type = "fusion")

nigra_fissions_gene_test <- test_CDS_enrichment("break_intervals_classified_elataquery.txt", chunks, 
                                                "Carex elata", breakpoint_type = "fission")

nigra_fusions_gene_test <- test_CDS_enrichment("break_intervals_classified_nigraquery.txt", chunks, 
                                               "Carex nigra", breakpoint_type = "fusion")

elata_fissions_gene_test <- test_CDS_enrichment("break_intervals_classified_nigraquery.txt", chunks, 
                                                "Carex nigra", breakpoint_type = "fission")
##

hirta_fusions_gene_test <- test_CDS_enrichment("break_intervals_classified_hirtaquery.txt", chunks, 
                                               "Carex hirta", breakpoint_type = "fusion")

riparia_fissions_gene_test <- test_CDS_enrichment("break_intervals_classified_hirtaquery.txt", chunks, 
                                                  "Carex hirta", breakpoint_type = "fission")

riparia_fusions_gene_test <- test_CDS_enrichment("break_intervals_classified_ripariaquery.txt", chunks, 
                                                 "Carex riparia", breakpoint_type = "fusion")

hirta_fissions_gene_test <- test_CDS_enrichment("break_intervals_classified_ripariaquery.txt", chunks, 
                                                "Carex riparia", breakpoint_type = "fission")

##

distans_fusions_gene_test <- test_CDS_enrichment("break_intervals_classified_distansquery.txt", chunks, 
                                                 "Carex distans", breakpoint_type = "fusion")

extensa_fissions_gene_test <- test_CDS_enrichment("break_intervals_classified_distansquery.txt", chunks, 
                                                  "Carex distans", breakpoint_type = "fission")

extensa_fusions_gene_test <- test_CDS_enrichment("break_intervals_classified_extensaquery.txt", chunks, 
                                                 "Carex extensa", breakpoint_type = "fusion")

distans_fissions_gene_test <- test_CDS_enrichment("break_intervals_classified_extensaquery.txt", chunks, 
                                                  "Carex extensa", breakpoint_type = "fission")

##

pendula_fusions_gene_test <- test_CDS_enrichment("break_intervals_classified_pendulaquery.txt", chunks, 
                                                 "Carex pendula", breakpoint_type = "fusion")

sylvatica_fissions_gene_test <- test_CDS_enrichment("break_intervals_classified_pendulaquery.txt", chunks, 
                                                    "Carex pendula", breakpoint_type = "fission")

sylvatica_fusions_gene_test <- test_CDS_enrichment("break_intervals_classified_sylvaticaquery.txt", chunks, 
                                                   "Carex sylvatica", breakpoint_type = "fusion")

pendula_fissions_gene_test <- test_CDS_enrichment("break_intervals_classified_sylvaticaquery.txt", chunks, 
                                                  "Carex sylvatica", breakpoint_type = "fission")

##

divulsa_fusions_gene_test <- test_CDS_enrichment("break_intervals_classified_divulsaquery.txt", chunks, 
                                                 "Carex divulsa", breakpoint_type = "fusion")

spicata_fissions_gene_test <- test_CDS_enrichment("break_intervals_classified_divulsaquery.txt", chunks, 
                                                  "Carex divulsa", breakpoint_type = "fission")

spicata_fusions_gene_test <- test_CDS_enrichment("break_intervals_classified_spicataquery.txt", chunks, 
                                                 "Carex spicata", breakpoint_type = "fusion")

divulsa_fissions_gene_test <- test_CDS_enrichment("break_intervals_classified_spicataquery.txt", chunks, 
                                                  "Carex spicata", breakpoint_type = "fission")

##

littledalei_fusions_gene_test <- test_CDS_enrichment("break_intervals_classified_littledaleiquery.txt", chunks, 
                                                     "Carex littledalei", breakpoint_type = "fusion")

myosuroides_fissions_gene_test <- test_CDS_enrichment("break_intervals_classified_littledaleiquery.txt", chunks, 
                                                     "Carex littledalei", breakpoint_type = "fission")

myosuroides_fusions_gene_test <- test_CDS_enrichment("break_intervals_classified_myosuroidesquery.txt", chunks, 
                                                     "Carex myosuroides", breakpoint_type = "fusion")

littledalei_fissions_gene_test <- test_CDS_enrichment("break_intervals_classified_myosuroidesquery.txt", chunks, 
                                                      "Carex myosuroides", breakpoint_type = "fission")

#################################################################################
##################### Enrichment of TEs #########################################
#################################################################################

# Helper function to normalise species names (reuse from cartoon function)
normalise_species_name <- function(species_name) {
  species_name <- gsub(" ", "_", species_name)
  if (!grepl("_", species_name)) {
    stop("Species name must be in format 'Genus_species' or 'Genus species'")
  }
  return(species_name)
}

# Helper function to match species name to data
match_species_name <- function(species_input, chunks_data) {
  normalised_input <- normalise_species_name(species_input)
  available_species <- unique(chunks_data$species)
  
  if (normalised_input %in% available_species) {
    return(normalised_input)
  }
  
  with_underscore <- gsub(" ", "_", species_input)
  without_underscore <- gsub("_", " ", species_input)
  
  if (with_underscore %in% available_species) {
    return(with_underscore)
  }
  if (without_underscore %in% available_species) {
    return(without_underscore)
  }
  
  stop(paste0("Species '", species_input, "' not found in chunks data.\n",
              "Available species: ", paste(available_species, collapse = ", ")))
}

#' Test for TE enrichment in breakpoint regions
#'
#' @param breakpoints_file Path to breakpoints file (TSV with: chr, break_start, break_end, class)
#' @param chunks_file Either a file path to chunks data (TSV) OR a pre-loaded data frame
#' @param species_name Name of the species (e.g., "Carex_distans", "Carex distans", or "distans")
#' @param breakpoint_type Which breakpoints to test: "fusion", "fission", or "both" (default: "both")
#' @param chunk_size Size of chunks in bp (default: 5000)
#' @param verbose Print detailed output (default: TRUE)
#' @return A list with test results and summary statistics
#' @export
test_TE_enrichment <- function(breakpoints_file,
                                chunks_file,
                                species_name,
                                breakpoint_type = "both",
                                chunk_size = 5000,
                                verbose = TRUE) {
  
  # Validate breakpoint_type
  if (!breakpoint_type %in% c("fusion", "fission", "both")) {
    stop("breakpoint_type must be one of: 'fusion', 'fission', or 'both'")
  }
  
  # Read the chunks data (or use provided data frame)
  if (is.character(chunks_file)) {
    chunks <- read_tsv(chunks_file, show_col_types = FALSE)
  } else if (is.data.frame(chunks_file)) {
    chunks <- chunks_file
  } else {
    stop("chunks_file must be either a file path (character) or a data frame")
  }
  
  # Match species name
  species_matched <- match_species_name(species_name, chunks)
  
  if (verbose) {
    cat("===========================================\n")
    cat("TE Enrichment Test\n")
    cat("===========================================\n")
    cat("Species:", species_matched, "\n")
    cat("Testing breakpoint type:", breakpoint_type, "\n\n")
  }
  
  # Read breakpoints
  breakpoints <- read_tsv(breakpoints_file, show_col_types = FALSE)
  
  # Filter breakpoints by type if needed
  if (breakpoint_type == "fusion") {
    breakpoints <- breakpoints %>%
      filter(grepl("fusion", class, ignore.case = TRUE))
  } else if (breakpoint_type == "fission") {
    breakpoints <- breakpoints %>%
      filter(grepl("fission", class, ignore.case = TRUE))
  }
  # If "both", use all breakpoints
  
  if (nrow(breakpoints) == 0) {
    stop("No breakpoints found for the specified type")
  }
  
  if (verbose) {
    cat("Number of breakpoint regions:", nrow(breakpoints), "\n")
  }
  
  # Filter chunks for this species
  chunks_sp <- chunks %>%
    filter(species == species_matched)
  
  # Calculate chunk start and end positions
  chunks_sp <- chunks_sp %>%
    mutate(chunk_start = (chunk_id - 1) * chunk_size,
           chunk_end = chunk_id * chunk_size)
  
  # Classify each chunk as in breakpoint region or not
  chunks_sp$in_breakpoint <- FALSE
  
  for (i in 1:nrow(breakpoints)) {
    chr <- breakpoints$chr[i]
    break_start <- breakpoints$break_start[i]
    break_end <- breakpoints$break_end[i]
    
    # Find chunks that overlap this breakpoint
    overlapping <- which(
      chunks_sp$chromosome == chr & 
        chunks_sp$chunk_end >= break_start & 
        chunks_sp$chunk_start <= break_end
    )
    
    chunks_sp$in_breakpoint[overlapping] <- TRUE
  }
  
  # Identify chunks with genes (non-NA repeat_classes)
  chunks_sp <- chunks_sp %>%
    mutate(has_TE = !is.na(TE_class_clean))
  
  # Create contingency table
  # Rows: in/out of breakpoint region
  # Cols: has/doesn't have oligocentromere
  contingency_table <- table(
    Breakpoint = chunks_sp$in_breakpoint,
    TE = chunks_sp$has_TE
  )
  
  if (verbose) {
    cat("\n--- Contingency Table ---\n")
    print(contingency_table)
    cat("\n")
  }
  
  # Calculate proportions
  in_breakpoint <- chunks_sp %>% filter(in_breakpoint)
  out_breakpoint <- chunks_sp %>% filter(!in_breakpoint)
  
  prop_in <- mean(in_breakpoint$has_TE)
  prop_out <- mean(out_breakpoint$has_TE)
  
  if (verbose) {
    cat("--- Proportions ---\n")
    cat("Chunks in breakpoint regions with TE:", 
        sum(in_breakpoint$has_TE), "/", nrow(in_breakpoint),
        sprintf(" (%.2f%%)\n", prop_in * 100))
    cat("Chunks outside breakpoint regions with TE:", 
        sum(out_breakpoint$has_TE), "/", nrow(out_breakpoint),
        sprintf(" (%.2f%%)\n", prop_out * 100))
    cat("Fold enrichment:", sprintf("%.2f\n", prop_in / prop_out))
    cat("\n")
  }
  
  # Perform Fisher's Exact Test
  fisher_test <- fisher.test(contingency_table)
  
  # Perform Chi-square test (if cell counts are sufficient)
  chisq_test <- NULL
  chisq_warning <- NULL
  tryCatch({
    chisq_test <- chisq.test(contingency_table)
  }, warning = function(w) {
    chisq_warning <<- w$message
  })
  
  if (verbose) {
    cat("--- Fisher's Exact Test ---\n")
    cat("Odds Ratio:", sprintf("%.3f\n", fisher_test$estimate))
    cat("95% CI:", sprintf("[%.3f, %.3f]\n", 
                           fisher_test$conf.int[1], fisher_test$conf.int[2]))
    cat("P-value:", format.pval(fisher_test$p.value, digits = 3), "\n")
    
    if (fisher_test$p.value < 0.001) {
      cat("*** Highly significant enrichment\n")
    } else if (fisher_test$p.value < 0.01) {
      cat("** Significant enrichment\n")
    } else if (fisher_test$p.value < 0.05) {
      cat("* Significant enrichment\n")
    } else {
      cat("Not significant\n")
    }
    cat("\n")
    
    if (!is.null(chisq_test)) {
      cat("--- Chi-square Test ---\n")
      cat("Chi-square statistic:", sprintf("%.3f\n", chisq_test$statistic))
      cat("P-value:", format.pval(chisq_test$p.value, digits = 3), "\n\n")
    } else if (!is.null(chisq_warning)) {
      cat("--- Chi-square Test ---\n")
      cat("Warning:", chisq_warning, "\n")
      cat("(Use Fisher's test for small sample sizes)\n\n")
    }
  }
  
  # Permutation test for additional validation
  if (verbose) {
    cat("--- Permutation Test (10,000 iterations) ---\n")
    cat("Running permutations...")
  }
  
  n_permutations <- 10000
  observed_diff <- prop_in - prop_out
  
  permuted_diffs <- replicate(n_permutations, {
    # Randomly shuffle the breakpoint labels
    shuffled <- sample(chunks_sp$in_breakpoint)
    prop_in_perm <- mean(chunks_sp$has_TE[shuffled])
    prop_out_perm <- mean(chunks_sp$has_TE[!shuffled])
    prop_in_perm - prop_out_perm
  })
  
  # Two-tailed p-value
  perm_pvalue <- mean(abs(permuted_diffs) >= abs(observed_diff))
  
  if (verbose) {
    cat(" Done!\n")
    cat("Observed difference in proportions:", sprintf("%.4f\n", observed_diff))
    cat("Permutation p-value:", format.pval(perm_pvalue, digits = 3), "\n")
    cat("\n")
  }
  
  # Summary interpretation
  if (verbose) {
    cat("===========================================\n")
    cat("SUMMARY\n")
    cat("===========================================\n")
    
    if (fisher_test$p.value < 0.05 && prop_in > prop_out) {
      cat("CONCLUSION: TEs are SIGNIFICANTLY ENRICHED\n")
      cat("in", breakpoint_type, "breakpoint regions.\n")
    } else if (fisher_test$p.value < 0.05 && prop_in < prop_out) {
      cat("CONCLUSION: TEs are SIGNIFICANTLY DEPLETED\n")
      cat("in", breakpoint_type, "breakpoint regions.\n")
    } else {
      cat("CONCLUSION: No significant enrichment or depletion of\n")
      cat("TEs in", breakpoint_type, "breakpoint regions.\n")
    }
    cat("\n")
  }
  
  # Return results
  results <- list(
    species = species_matched,
    breakpoint_type = breakpoint_type,
    n_breakpoint_regions = nrow(breakpoints),
    n_chunks_total = nrow(chunks_sp),
    n_chunks_in_breakpoint = nrow(in_breakpoint),
    n_chunks_out_breakpoint = nrow(out_breakpoint),
    
    # Oligocentromere counts
    n_TE_in_breakpoint = sum(in_breakpoint$has_TE),
    n_TE_out_breakpoint = sum(out_breakpoint$has_TE),
    
    # Proportions
    prop_TE_in_breakpoint = prop_in,
    prop_TE_out_breakpoint = prop_out,
    fold_enrichment = prop_in / prop_out,
    
    # Contingency table
    contingency_table = contingency_table,
    
    # Statistical tests
    fisher_test = fisher_test,
    fisher_pvalue = fisher_test$p.value,
    fisher_odds_ratio = as.numeric(fisher_test$estimate),
    fisher_ci_lower = fisher_test$conf.int[1],
    fisher_ci_upper = fisher_test$conf.int[2],
    
    chisq_test = chisq_test,
    chisq_pvalue = if(!is.null(chisq_test)) chisq_test$p.value else NA,
    
    permutation_pvalue = perm_pvalue,
    permutation_diffs = permuted_diffs,
    observed_diff = observed_diff,
    
    # Detailed chunk data for further analysis
    chunks_data = chunks_sp
  )
  
  class(results) <- c("TE_enrichment_test", "list")
  return(results)
}

###

elata_fusions_TE_test <- test_TE_enrichment("break_intervals_classified_elataquery.txt", chunks, 
                                            "Carex_elata", breakpoint_type = "fusion")

nigra_fissions_TE_test <- test_TE_enrichment("break_intervals_classified_elataquery.txt", chunks, 
                                             "Carex_elata", breakpoint_type = "fission")

nigra_fusions_TE_test <- test_TE_enrichment("break_intervals_classified_nigraquery.txt", chunks, 
                                            "Carex_nigra", breakpoint_type = "fusion")

elata_fissions_TE_test <- test_TE_enrichment("break_intervals_classified_nigraquery.txt", chunks, 
                                             "Carex_nigra", breakpoint_type = "fission")

##

hirta_fusions_TE_test <- test_TE_enrichment("break_intervals_classified_hirtaquery.txt", chunks, 
                                            "Carex_hirta", breakpoint_type = "fusion")

riparia_fissions_TE_test <- test_TE_enrichment("break_intervals_classified_hirtaquery.txt", chunks, 
                                               "Carex_hirta", breakpoint_type = "fission")

riparia_fusions_TE_test <- test_TE_enrichment("break_intervals_classified_ripariaquery.txt", chunks, 
                                              "Carex_riparia", breakpoint_type = "fusion")

hirta_fissions_TE_test <- test_TE_enrichment("break_intervals_classified_ripariaquery.txt", chunks, 
                                             "Carex_riparia", breakpoint_type = "fission")

##

distans_fusions_TE_test <- test_TE_enrichment("break_intervals_classified_distansquery.txt", chunks, 
                                            "Carex_distans", breakpoint_type = "fusion")

extensa_fissions_TE_test <- test_TE_enrichment("break_intervals_classified_distansquery.txt", chunks, 
                                               "Carex_distans", breakpoint_type = "fission")

extensa_fusions_TE_test <- test_TE_enrichment("break_intervals_classified_extensaquery.txt", chunks, 
                                              "Carex_extensa", breakpoint_type = "fusion")

distans_fissions_TE_test <- test_TE_enrichment("break_intervals_classified_extensaquery.txt", chunks, 
                                               "Carex_extensa", breakpoint_type = "fission")

##

pendula_fusions_TE_test <- test_TE_enrichment("break_intervals_classified_pendulaquery.txt", chunks, 
                                              "Carex_pendula", breakpoint_type = "fusion")

sylvatica_fissions_TE_test <- test_TE_enrichment("break_intervals_classified_pendulaquery.txt", chunks, 
                                                 "Carex_pendula", breakpoint_type = "fission")

sylvatica_fusions_TE_test <- test_TE_enrichment("break_intervals_classified_sylvaticaquery.txt", chunks, 
                                                "Carex_sylvatica", breakpoint_type = "fusion")

pendula_fissions_TE_test <- test_TE_enrichment("break_intervals_classified_sylvaticaquery.txt", chunks, 
                                               "Carex_sylvatica", breakpoint_type = "fission")

##

divulsa_fusions_TE_test <- test_TE_enrichment("break_intervals_classified_divulsaquery.txt", chunks, 
                                              "Carex_divulsa", breakpoint_type = "fusion")

spicata_fissions_TE_test <- test_TE_enrichment("break_intervals_classified_divulsaquery.txt", chunks, 
                                               "Carex_divulsa", breakpoint_type = "fission")

spicata_fusions_TE_test <- test_TE_enrichment("break_intervals_classified_spicataquery.txt", chunks, 
                                              "Carex_spicata", breakpoint_type = "fusion")

divulsa_fissions_TE_test <- test_TE_enrichment("break_intervals_classified_spicataquery.txt", chunks, 
                                               "Carex_spicata", breakpoint_type = "fission")

##

littledalei_fusions_TE_test <- test_TE_enrichment("break_intervals_classified_littledaleiquery.txt", chunks, 
                                                  "Carex_littledalei", breakpoint_type = "fusion")

myosuroides_fissions_TE_test <- test_TE_enrichment("break_intervals_classified_littledaleiquery.txt", chunks, 
                                                   "Carex_littledalei", breakpoint_type = "fission")

myosuroides_fusions_TE_test <- test_TE_enrichment("break_intervals_classified_myosuroidesquery.txt", chunks, 
                                                  "Carex_myosuroides", breakpoint_type = "fusion")

littledalei_fissions_TE_test <- test_TE_enrichment("break_intervals_classified_myosuroidesquery.txt", chunks, 
                                                   "Carex_myosuroides", breakpoint_type = "fission")
