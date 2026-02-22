species_consensi <- read.csv("species_consensi.csv")

# function incorporating reverse strand

revcomp <- function(seq) {
  chartr("ACGT", "TGCA", paste(rev(strsplit(seq, "")[[1]]), collapse = ""))
}

get_kmers_canonical <- function(seq, k) {
  seq <- toupper(gsub("\\s+", "", seq))
  if (nchar(seq) < k) return(character(0))
  
  kmers <- sapply(1:(nchar(seq) - k + 1),
                  function(i) substr(seq, i, i + k - 1))
  
  rc <- sapply(kmers, revcomp)
  
  canonical <- ifelse(kmers < rc, kmers, rc)
  unique(canonical)
}

#applying function
k <- 8
kmers_list <- lapply(species_consensi$Consensus, get_kmers_canonical, k = k)
kmers_set  <- lapply(kmers_list, unique)

#now I have k-mer set for each species-level consensus

### create Jaccard matrix

#function
jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

#apply
n <- length(kmers_set)
jac_mat <- matrix(NA, n, n)

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    jac_mat[i, j] <- jaccard(kmers_set[[i]], kmers_set[[j]])
  }
}
labels <- paste(species_consensi$Species,
                species_consensi$Family,
                sep = "_")

rownames(jac_mat) <- labels
colnames(jac_mat) <- labels

## printing similarity matrix as csv

write.csv(
  jac_mat,
  file = "jaccard_similarity_k8.csv",
  quote = FALSE
)

#average k-mer similarity between all non-Kobresia Carol consensus sequences

carol_idx <- grepl("_Carol$", colnames(jac_mat))
carol_mat <- jac_mat[carol_idx, carol_idx]
carol_vals <- carol_mat[upper.tri(carol_mat)]  # take each pair once
mean(carol_vals, na.rm = TRUE)

exclude <- c("Carex littledalei_Carol", "Carex myosuroides_Carol")

carol_keep <- grepl("_Carol$", colnames(jac_mat)) &
  !colnames(jac_mat) %in% exclude

carol_mat_filt <- jac_mat[carol_keep, carol_keep]

carol_vals_filt <- carol_mat_filt[upper.tri(carol_mat_filt)]
mean(carol_vals_filt, na.rm = TRUE)

## making a histogram
# Convert similarity matrix to vector (excluding diagonal)
mat_vec <- as.vector(jac_mat[upper.tri(jac_mat)])

# Exploratory histogram
hist(mat_vec,
     breaks = 100,                  # number of bins
     main = "Histogram of pairwise Jaccard distances",
     xlab = "Jaccard similarity",
     col = "steelblue",
     border = "white")


############ plotting 


library(igraph)

threshold <- 0.1

adj <- jac_mat >= threshold
diag(adj) <- FALSE

g <- graph_from_adjacency_matrix(adj, mode = "undirected")
comp <- components(g)

ord <- names(comp$membership)[order(comp$membership)]

jac_family <- jac_masked[ord, ord]

family_sizes <- table(comp$membership)
gaps <- cumsum(family_sizes)


species_labels <- sub("_.*$", "", rownames(jac_family))
species_labels_italic <- parse(
  text = paste0("italic('", species_labels, "')")
)
pheatmap(
  jac_family,
  color = colorRampPalette(c("white", "red"))(100),
  na_col = "white",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gaps,
  gaps_col = gaps,
  border_color = NA,
  labels_row = species_labels_italic,
  fontsize_row = 7
)
