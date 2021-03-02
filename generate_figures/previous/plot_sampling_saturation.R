require(ggplot2)
require(dplyr)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
PREFIX_DATA <- "rep_1_n_sim_1000_n_imports_padded_0"
PREFIX <- "rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1"
PREFIX_1 <- "rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_T"
PREFIX_2 <- "rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_F"

CLUSTERS_1 <- paste(WORKDIR, "/clusters/", PREFIX_1, "_clusters.txt", sep = "")
CLUSTERS_2 <- paste(WORKDIR, "/clusters/", PREFIX_2, "_clusters.txt", sep = "")
METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
OUTDIR <- paste(WORKDIR, "/figures/sampling_saturation", sep = "")

source(paste(WORKDIR, "figures/scripts/plotting_utility_functions.R", sep = "/"))
system(command = paste("mkdir -p", OUTDIR))

clusters_1 <- read.delim(file = CLUSTERS_1)
clusters_2 <- read.delim(file = CLUSTERS_2)
metadata <- read.delim(file = METADATA, stringsAsFactors = F, sep = "\t", quote = "")

downsample_data <- function(cluster_data_by_tip, n_seqs) {
  selected_idxs <- sample(
    x = 1:nrow(cluster_data_by_tip), replace = F, size = n_seqs)
  return(cluster_data_by_tip[selected_idxs, ])
}

get_clusters_in_downsampled_data <- function(downsampled_data) {
  clusters <- downsampled_data %>%
    group_by(cluster_idx) %>%
    summarise(size = n())
  return(clusters)
}

get_n_singletons_in_downsampled_data <- function(downsampled_clusters) {
  return(sum(downsampled_clusters$size == 1))
}

get_n_chains_in_downsampled_data <- function(downsampled_clusters) {
  return(sum(downsampled_clusters$size > 1))
}

# Plot # clusters with increasing # swiss sequences
make_plot <- function(clusters, metadata) {
  cluster_data_by_tip <- get_cluster_data_by_tip(
    cluster_data = clusters, metadata = metadata)
  
  n_swiss_seqs <- floor(seq(from = 0.1, to = 1, by = 0.1) * nrow(cluster_data_by_tip))
  
  downsamples <- lapply(
    X = n_swiss_seqs, 
    FUN = downsample_data, 
    cluster_data_by_tip = cluster_data_by_tip)
  
  downsamples_clusters <- lapply(
    X = downsamples,
    FUN = get_clusters_in_downsampled_data)
  
  n_singletons <- unlist(lapply(
    X = downsamples_clusters,
    FUN = get_n_singletons_in_downsampled_data))
  
  n_chains <- unlist(lapply(
    X = downsamples_clusters,
    FUN = get_n_chains_in_downsampled_data))
  
  downsample_summary <- data.frame(
    n_swiss_seqs = n_swiss_seqs,
    n_singletons = n_singletons,
    n_chains = n_chains)
  
  p <- ggplot(
    data = downsample_summary,
    aes(x = n_swiss_seqs)) + 
    geom_point(aes(y = n_singletons, shape = "No. singletons")) +
    geom_point(aes(y = n_chains, shape = "No. transmission chains")) + 
    scale_shape_manual(
      name = element_blank(),
      values = list("No. singletons" = 8, 
                    "No. transmission chains" = 3)) + 
    labs(x = "No. Swiss sequences considered", y = element_blank()) + 
    theme_bw()
  
  return(p)
}

p1 <- make_plot(clusters = clusters_1, metadata = metadata)
p2 <- make_plot(clusters = clusters_2, metadata = metadata)

library(gtable)
library(grid)
g1 <- ggplotGrob(p1  + ggtitle("a)"))
g2 <- ggplotGrob(p2 + ggtitle("b)"))
g <- rbind(g1, g2, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths)

png(
  file = paste(OUTDIR, "/", PREFIX, "_sampling_saturation.png", sep = ""),
  width = 6.5, height = 4, units = "in", res = 300)
grid.newpage()
grid.draw(g)
dev.off()

# png(
#   file = paste(OUTDIR, paste(PREFIX, "_sampling_saturation.png", sep = ""), sep = "/"),
#   width = 6.5, height = 2, units = "in", res = 300)
# show(p)
# dev.off()
