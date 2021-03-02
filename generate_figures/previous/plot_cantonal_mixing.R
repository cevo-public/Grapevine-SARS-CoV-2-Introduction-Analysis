# Plot network of which cantons cluster with which other cantons.
# Highlight significant edges which occur more often than 90th percentile value
# in permutation analysis (see script permute_swiss_clusters.R). 90th percentile
# b/c it's a one-tailed test.

require(ggplot2)
require(dplyr)
require(igraph)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
PREFIX_DATA <- "rep_3_n_sim_1000_n_imports_padded_0"
PREFIX <- "rep_3_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_F_viollier_only"
METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")

OUTDIR <- paste(WORKDIR, "figures/cantonal_mixing", sep = "/")
CANTONAL_MIXING_RESULTS <- paste(WORKDIR, "/cantonal_mixing/", PREFIX, "_permutation_results.txt", sep = "")

source(paste(WORKDIR, "figures/scripts/plotting_utility_functions.R", sep = "/"))
system(command = paste("mkdir -p", OUTDIR))

# Load data
metadata <- read.delim(file = METADATA, quote = "")
permutation_results <- read.delim(file = CANTONAL_MIXING_RESULTS, sep = "\t")

# Complete permutation results with cantons that didn't mix but do have seqs
cantons <- as.character(unique(unlist(metadata %>% filter(originating_lab == "Viollier AG") %>% select(division))))
permutation_results$canton_A <- factor(
  x = permutation_results$canton_A, 
  levels = cantons)
permutation_results$canton_B <- factor(
  x = permutation_results$canton_B, 
  levels = cantons)
permutation_results <- permutation_results %>%
  complete(canton_A, canton_B)

true_mixing <- pivot_wider(
  data = permutation_results_2[c("canton_A", "canton_B", "true_mixing_tally")],
  names_from = "canton_B",
  values_from = "true_mixing_tally")
rownames(true_mixing) <- true_mixing$canton_A
true_mixing <- as.matrix(true_mixing)

true_mixing[is.na(true_mixing)] <- 0

# Create network
net <- igraph::graph_from_adjacency_matrix(
  true_mixing[, 2:ncol(true_mixing)],  # eliminate 1st column with canton names 
  weighted = T,
  mode = "undirected")

# Assign significance to edges
sig_edges <- c()
for(edge in E(net)){
  verts = ends(net, edge)
  if (verts[2] != verts[1]) {
    if (permutation_results[
      permutation_results$canton_A == verts[2] & 
      permutation_results$canton_B == verts[1], "more_mixing_significance_alpha_1"] == "TRUE") {
      sig_edges <- c(sig_edges, "#E31A1C")
    } else if (permutation_results[
      permutation_results$canton_A == verts[2] & 
      permutation_results$canton_B == verts[1], "more_mixing_significance_alpha_2"] == "TRUE") {
      sig_edges <- c(sig_edges, "#FB9A99")
    } else if (permutation_results[
      permutation_results$canton_A == verts[2] & 
      permutation_results$canton_B == verts[1], "less_mixing_significance_alpha_1"] == "TRUE") {
      sig_edges <- c(sig_edges, "#1F78B4")
    } else if (permutation_results[
      permutation_results$canton_A == verts[2] & 
      permutation_results$canton_B == verts[1], "less_mixing_significance_alpha_2"] == "TRUE") {
      sig_edges <- c(sig_edges, "#A6CEE3")
    } else {
      sig_edges <- c(sig_edges, adjustcolor("grey", alpha.f = .5))
    }
  }
}

seqs_per_canton <- metadata %>% 
  filter(originating_lab == "Viollier AG") %>%
  group_by(division) %>%
  summarize(n_seqs = n())

NODE_SF <- 10
node.size <- setNames(seqs_per_canton$n_seqs / NODE_SF, seqs_per_canton$division)
node.size <- node.size[match(names(V(net)), names(node.size))]

canton_loc_data <- canton_loc_data[match(names(V(net)), canton_loc_data$cantons), ]
capitals_layout <- canton_loc_data[c("x_new", "y_new")]

png(
  file = paste(OUTDIR, paste(PREFIX, "_viollier_only_cantonal_mixing.png", sep = ""), sep = "/"),
  width = 6.5, height = 7.5, units = "in", res = 300)

par(mar=c(0,0,0,0)+.1)
plot.igraph(
  net,
  layout = as.matrix(capitals_layout),
  edge.width = E(net)$weight,
  edge.color = sig_edges,
  vertex.color = adjustcolor("grey", alpha.f = .5),
  vertex.label.family="Helvetica",
  vertex.label.font=2,
  vertex.label.color = "black",
  vertex.size = node.size)
dev.off()


