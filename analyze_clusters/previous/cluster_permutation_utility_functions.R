run_permutation <- function(clusters, cantons) {
  # Permute samples within clusters. Return dataframe with cols: canton_A, 
  # canton_B, mixing_tally (# times cantons co-occur in a cluster)
  permuted_clusters <- permute_clusters(clusters)
  permuted_mixing <- tally_mixing(clusters = permuted_clusters, cantons = cantons)
  permuted_mixing_df <- as.data.frame(permuted_mixing)
  permuted_mixing_df$canton_A <- rownames(permuted_mixing)
  
  permutation_results <- permuted_mixing_df %>% 
    pivot_longer(
      cols = 1:length(cantons), 
      names_to = "canton_B", 
      values_to = "mixing_tally")
  
  return(permutation_results)
}

permute_clusters <- function(clusters, timout = Inf, verbose = F) {
  # Permute samples in clusters until a valid configuration is found; return it
  is_valid_draw <- F
  i <- 0
  while (!is_valid_draw & i < timout) {
    i <- i + 1
    bucket <- unlist(clusters)
    scrambled_bucket <- sample(bucket)
    cluster_sizes <- unlist(lapply(X = clusters, FUN = length))
    draw <- split_at(vec = scrambled_bucket, idxs = cumsum(cluster_sizes) + 1)
    is_valid_draw <- get_is_valid_draw_samples(draw)
  }
  if (is_valid_draw) {
    if (verbose) {
      print(paste(i, "shuffles required to find a valid configuration."))
    }
    return(draw)
  }
  stop(paste("Timeout after", i, "shuffles without a valid configuration."))
}

split_at <- function(vec, idxs) {
  # Returns a list, which is the vec split at idxs 
  # First entry always starts at 1, subsequent entries start at idxs
  return(unname(split(vec, cumsum(seq_along(vec) %in% idxs)))) 
}

get_is_valid_draw_samples <- function(draw) {
  # Returns whether all clusters are valid
  return(all(unlist(lapply(X = draw, FUN = get_is_valid_cluster_samples))))
}

get_is_valid_cluster_samples <- function(cluster) {
  # A samples cluster is valid if the whole cluster isn't from same canton
  return(length(unique(cluster)) > 1)
}

tally_mixing <- function(clusters, cantons) {
  # Return symmetric matrix of # times cantons co-occur in clusters
  mixing <- matrix(data = 0, ncol = length(cantons), nrow = length(cantons))
  rownames(mixing) <- colnames(mixing) <- cantons
  for (cluster in clusters) {
    unique_cantons <- unique(cluster)
    for (canton in unique_cantons) {
      mixing[canton, unique_cantons] <- mixing[canton, unique_cantons] + 1
    }
  }
  diag(mixing) <- NA
  return(mixing)
}

get_n_mixes <- function(canton_1, canton_2, clusters) {
  # For testing purposes, return how often 2 cantons co-occur in a cluster
  n_mixes <- 0
  for (cluster in clusters) {
    if (canton_1 %in% cluster & canton_2 %in% cluster) {
      n_mixes <- n_mixes + 1
    }
  }
  return(n_mixes)
}