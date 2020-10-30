# This script takes clusters, metadata about which canton each tip is from & 
# tallies the # times samples from each canton cluster together. Results are 
# compared to results from randomly permuting samples amongst clusters,  
# conditioned on clusters having samples from at least 2 different cantons. Two 
# types of permutation are done, one permuting samples from qualifying clusters
# and one collapsing samples from the same canton in the same cluster, then 
# permuting. The first assumes the max. # exchanges bewteen cantons has occurred
# and the second assumes the minimum. 
# Results are written out as a table listing the significance of the true # 
# co-occurances compared to the # co-occurances in permuted clusters for each 
# pair of cantons. 

# Hardcoded things 
WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
PREFIX_DATA <- "rep_2_n_sim_1000_n_imports_padded_0"
PREFIX <- "rep_2_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_F"
OUTDIR <- paste(WORKDIR, "cantonal_mixing", sep = "/")

CLUSTERS <- paste(WORKDIR, "/clusters/", PREFIX, "_clusters.txt", sep = "")
METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
N_PERMUTATIONS <- 100
ONE_TAIL_ALPHA_1 <- 0.05  
ONE_TAIL_ALPHA_2 <- 0.1
N_1 <- N_PERMUTATIONS - ONE_TAIL_ALPHA_1 * N_PERMUTATIONS # true value must be greater/smaller than Nth largest/smallest value from permutation analysis for edge to be significant
N_2 <- N_PERMUTATIONS - ONE_TAIL_ALPHA_2 * N_PERMUTATIONS # true value must be greater/smaller than Nth largest/smallest value from permutation analysis for edge to be significant

ASSUME_MIN_N_EXCHANGES <- T

source(paste(WORKDIR, "scripts/analyze_clusters/cluster_permutation_utility_functions.R", sep = "/"))
source(paste(WORKDIR, "figures/scripts/plotting_utility_functions.R", sep = "/"))
system(command = paste("mkdir -p", OUTDIR))

# Load data 
cluster_data <- read.delim(file = CLUSTERS)
metadata <- read.table(file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
cluster_data_by_tip <- get_cluster_data_by_tip(
  cluster_data = cluster_data, 
  metadata = metadata) %>% 
  filter(originating_lab == "Viollier AG") %>%  # take only viollier samples
  group_by(cluster_idx) %>%
  filter(length(unique(division)) > 1)  # only consider clusters with at least 2 cantons

get_permutation_results <- function(cluster_data_by_tip, ASSUME_MIN_N_EXCHANGES, N_PERMUTATIONS) {
  clusters <- split(cluster_data_by_tip$division, cluster_data_by_tip$cluster_idx)  # a list, names are cluster idxs and values are each seq's canton
  cantons <- unique(cluster_data_by_tip$division)
  
  if (ASSUME_MIN_N_EXCHANGES) {
    clusters_w_cantons_collapsed <- lapply(X = clusters, FUN = unique)  # collapse seqs from same canton in each cluster to a single seq
    clusters <- clusters_w_cantons_collapsed
  }
  # Run permutations 
  for (i in 1:N_PERMUTATIONS) {
    if (i %% 10 == 0) {
      print(paste("Conducting permutation", i))
    }
    permutation_i_results <- run_permutation(
      clusters = clusters, cantons = cantons)
    permutation_i_results$permutation_idx <- i
    if (i == 1) {
      permutation_results <- permutation_i_results
    } else {
      permutation_results <- rbind(permutation_results, permutation_i_results)
    }
  }
  return(permutation_results)
}

get_significance_bounds <-function(permutation_results, ASSUME_MIN_N_EXCHANGES, N_1, N_2, ONE_TAIL_ALPHA_1, ONE_TAIL_ALPHA_2) {
  # Summarize permutations
  if (ASSUME_MIN_N_EXCHANGES) {
    # one-tailed test for less mixing than expected
    significance_bounds <- permutation_results %>%
      group_by(canton_A, canton_B) %>%
      summarize(
        n_permutations = n(),
        permuted_values = paste0(x = sort(x = mixing_tally, decreasing = T), collapse = ", "),
        n_1_th_smallest_value = sort(x = mixing_tally, decreasing = T)[N_1],
        n_2_th_smallest_value = sort(x = mixing_tally, decreasing = T)[N_2],
        n_1 = N_1, n_2 = N_2,
        alpha_1 = ONE_TAIL_ALPHA_1, alpha_2 = ONE_TAIL_ALPHA_2)
  } else {
    # one-tailed test for less mixing than expected
    significance_bounds <- permutation_results %>%
      group_by(canton_A, canton_B) %>%
      summarize(
        n_permutations = n(),
        permuted_values = paste0(x = sort(x = mixing_tally, decreasing = T), collapse = ", "),
        n_1_th_largest_value = sort(x = mixing_tally)[N_1],
        n_2_th_largest_value = sort(x = mixing_tally)[N_2],
        n_1 = N_1, n_2 = N_2,
        alpha_1 = ONE_TAIL_ALPHA_1, alpha_2 = ONE_TAIL_ALPHA_2)
  }
  return(significance_bounds)
}

perm_results_less_mixing <- get_permutation_results(
  cluster_data_by_tip = cluster_data_by_tip, 
  N_PERMUTATIONS = N_PERMUTATIONS, ASSUME_MIN_N_EXCHANGES = T)
perm_results_more_mixing <- get_permutation_results(
  cluster_data_by_tip = cluster_data_by_tip, 
  N_PERMUTATIONS = N_PERMUTATIONS, ASSUME_MIN_N_EXCHANGES = F)

sig_bounds_less_mixing <- get_significance_bounds(
  ASSUME_MIN_N_EXCHANGES = T,
  permutation_results = perm_results_less_mixing, 
  N_1, N_2, ONE_TAIL_ALPHA_1, ONE_TAIL_ALPHA_2)

sig_bounds_more_mixing <- get_significance_bounds(
  ASSUME_MIN_N_EXCHANGES = F,
  permutation_results = perm_results_more_mixing, 
  N_1, N_2, ONE_TAIL_ALPHA_1, ONE_TAIL_ALPHA_2)

# Compare true mixing to permutations
cantons <- unique(cluster_data_by_tip$division)
clusters <- split(cluster_data_by_tip$division, cluster_data_by_tip$cluster_idx)  # a list, names are cluster idxs and values are each seq's canton

true_mixing <- tally_mixing(clusters = clusters, cantons = cantons)
true_mixing_df <- as.data.frame(true_mixing)
true_mixing_df$canton_A <- rownames(true_mixing)

true_mixing_results <- true_mixing_df %>% 
  pivot_longer(
    cols = 1:length(cantons), 
    names_to = "canton_B", 
    values_to = "true_mixing_tally")

results <- merge(
  x = sig_bounds_less_mixing %>% select(-c(permuted_values)), 
  y = sig_bounds_more_mixing %>% select(-c(permuted_values)), all = T,
  by = c("canton_A", "canton_B", "n_permutations", "n_1", "n_2", "alpha_1","alpha_2"))
results <- merge(
  x = results, y = true_mixing_results, all = T)

results$more_mixing_significance_alpha_1 <- results$true_mixing_tally > results$n_1_th_largest_value
results$more_mixing_significance_alpha_2 <- results$true_mixing_tally > results$n_2_th_largest_value
results$less_mixing_significance_alpha_1 <- results$true_mixing_tally < results$n_1_th_smallest_value
results$less_mixing_significance_alpha_2 <- results$true_mixing_tally < results$n_2_th_smallest_value

# Write out results
OUTFILE <- paste(OUTDIR, "/", PREFIX, "_viollier_only_permutation_results.txt", sep = "")

write.table(
  x = results,
  file = OUTFILE,
  quote = F, sep = "\t", row.names = F, col.names = T)

