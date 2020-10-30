# This script takes a tree and tip location information and returns 
# cluster assignments based on a user-defined max number of monophyletic foreign
# sub-clades (exports) and a max number of foreign clades budding from any internal
# branch. Note: the tree is NOT allowed to be one big cluster (at the highest 
# level of aggregation, clades descending from the root will be clusters). Also
# note that the number of budding (pendant) subclades is calculated as the maximum 
# # pendant (budding subclades) to any swiss child plus the number of non-swiss
# siblings. This is conservative because it assumes polytomies would be resolved
# such that the maximum # non-swiss clades bud off in a row.

require(treeio)
require(tidytree)

IS_TEST_CASE <- F

# --- start test case code ---

# WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-22_ch_cluster_analysis"
# TEST_CASE_NAME <- 1
# MAX_NONFOCAL_SUBCLADES <- 3
# MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES <- 1
# POLYTOMIES_ARE_SWISS <- T
#
# TREE <- paste(WORKDIR, "/test_cases/test_case_", TEST_CASE_NAME, "/tree.nex", sep = "")
# METADATA <- paste(WORKDIR, "/test_cases/test_case_", TEST_CASE_NAME, "/metadata.txt", sep = "")
# OUTDIR <- paste(WORKDIR, "/test_cases/test_case_", TEST_CASE_NAME, "/results", sep = "")
# PREFIX <- paste("m", MAX_NONFOCAL_SUBCLADES, "p", MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES, "l", POLYTOMIES_ARE_SWISS, sep = "_")
# IS_VERBOSE <- T
# FOCAL_LOC <- "Switzerland"
# IS_TEST_CASE <- T

# WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-06_ch_cluster_analysis_final"
# PREFIX_DATA <- "rep_1_n_sim_1000_n_imports_padded_0"
# MAX_NONFOCAL_SUBCLADES <- 3
# MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES <- 1
# POLYTOMIES_ARE_SWISS <- T
# 
# TREE <- paste(WORKDIR, "/dated_trees/", PREFIX_DATA, ".timetree.nex", sep = "")
# METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
# OUTDIR <- paste(WORKDIR, "/clusters", sep = "")
# PREFIX <- paste("rep_1_n_sim_1000_n_imports_padded_0_m_", MAX_NONFOCAL_SUBCLADES, "_p_", MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES, "_swiss_polytomies_", POLYTOMIES_ARE_SWISS, sep = "")
# IS_VERBOSE <- T
# FOCAL_LOC <- "Switzerland"
# IS_TEST_CASE <- T
# UTILITY_FUNCTIONS <- paste(WORKDIR, "scripts/utility_functions.R", sep = "/")

# --- end test case code ---

if (!IS_TEST_CASE) {
  parser <- argparse::ArgumentParser()
  parser$add_argument("--tree", type="character")
  parser$add_argument("--metadata", type="character")
  parser$add_argument("--outdir", type="character")
  parser$add_argument("--maxtotalsubclades", type="double")
  parser$add_argument("--maxconsecutivesubclades", type="double")
  parser$add_argument("--prefix", type="character")
  parser$add_argument("--polytomiesareswiss", action = "store_true")
  parser$add_argument("--utilityfunctions", type = "character")
  
  args <- parser$parse_args()
  
  TREE <- args$tree
  METADATA <- args$metadata
  OUTDIR <- args$outdir
  MAX_NONFOCAL_SUBCLADES <- args$maxtotalsubclades
  MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES <- args$maxconsecutivesubclades
  PREFIX <- args$prefix
  POLYTOMIES_ARE_SWISS <- args$polytomiesareswiss
  UTILITY_FUNCTIONS <- args$utilityfunctions
  
  FOCAL_LOC <- "Switzerland"
  IS_VERBOSE <- F
}

print(paste("MAX_NONFOCAL_SUBCLADES:", MAX_NONFOCAL_SUBCLADES))
print(paste("MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES:", MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES))
print(paste("POLYTOMIES_ARE_SWISS:", POLYTOMIES_ARE_SWISS))

system(command = paste("mkdir -p", OUTDIR))

source(UTILITY_FUNCTIONS)

# Load data
tree <- treeio::read.beast(file = TREE)
metadata <- read.delim(file = METADATA, stringsAsFactors = F, sep = "\t", quote = "") %>% unique.data.frame()  # remove duplicates (to account for a concatenate mistake in metadata)
tree_data <- tidytree::as_tibble(tree)
tree_data <- merge(
  x = tree_data, y = metadata[c("strain", "country_recoded")],
  by.x = "label", by.y = "strain", all.x = T)
tree_ids <- tree_data$label[!is.na(tree_data$label)]  # tips only
n_tips <- length(tree_ids)

# Functions
get_dummy_cluster_for_nonfocal_node <- function(node_data, verbose = F) {
  # Returns a placeholder cluster representing a non-swiss tip or clade 
  if (verbose) {
    print(paste("...node", node_data$node,  "is a non-Swiss clade"))
  }
  return(data.frame(
    "parent_node" = NA, 
    "size" = NA, 
    "foreign_mrca" = NA,
    "foreign_tmrca" = NA,
    "foreign_tmrca_CI" = NA, 
    "foreign_neighbors" = NA,
    "ch_mrca" = NA,
    "ch_tmrca" = NA, 
    "ch_tmrca_CI" = NA, 
    "tips" = NA,
    "tip_nodes" = as.character(node_data$node),
    "n_basal_nonswiss_clades" = 0,
    "n_nonfocal_subclades" = 1))
}

initiate_swiss_cluster <- function(node_data, verbose = F) {
  # Takes node data in tidytree treedata format for a single Swiss tip
  # Returns a cluster object for the single tip
  if (verbose) {
    print(paste("...node", node_data$node, "is Swiss tip"))
  }
  return(data.frame(
    "size" = 1, 
    "ch_mrca" = node_data$node,
    "ch_tmrca" = node_data$date, 
    "ch_tmrca_CI" = node_data$CI_date,
    "tips" = as.character(node_data$label),
    "tip_nodes" = as.character(node_data$node),
    "n_basal_nonswiss_clades" = 0,
    "n_nonfocal_subclades" = 0))
}

propogate_swiss_cluster <- function(child_clusters, node_data, n_nonfocal_subclades, n_basal_nonswiss_clades, verbose = F) {
  # Given  there's only one child with a swiss cluster, update information on # 
  # non-focal subclades and basal non-swiss sequences.
  if (verbose) {
    print(paste("...only one child of",
                node_data$node, "is swiss but criteria not exceeded by siblings, propogating the cluster up with", 
                n_nonfocal_subclades, "nonfocal_subclades and",
                n_basal_nonswiss_clades, "basal_nonswiss_clades"))
  }
  child_clusters[[1]][["n_basal_nonswiss_clades"]] <- n_basal_nonswiss_clades
  child_clusters[[1]][["n_nonfocal_subclades"]] <- n_nonfocal_subclades
  if (n_basal_nonswiss_clades == 1) {
    # If connecting swiss cluster to first non-swiss node, also update foreign mrca
    if (verbose) {
      print(paste("updating foreign mrca to", node_data$node))
    }
    child_clusters[[1]][["foreign_mrca"]] <- node_data$node
    child_clusters[[1]][["foreign_tmrca"]] <- node_data$date
    child_clusters[[1]][["foreign_tmrca_CI"]] <- I(list(node_data$CI_date))
  }
  return(child_clusters[[1]])
}

merge_child_swiss_clusters <- function(child_clusters, node_data, n_nonfocal_subclades, n_basal_nonswiss_clades, verbose = F) {
  # Given cluster data for clusters descending from the parent node, and the 
  # parent node data in treedata format, merge the child clusters. If there's 
  # only one child with a swiss cluster, update information on # non-focal 
  # subclades and basal non-swiss sequences.
  if (verbose) {
    print(paste("...some children of",
                node_data$node, "are swiss & merge criteria satisfied, returning merged cluster with", 
                n_nonfocal_subclades, "nonfocal_subclades and",
                n_basal_nonswiss_clades, "basal_nonswiss_clades"))
  }
  cluster_sizes <- unlist(lapply(X = child_clusters, FUN = get_cluster_size))
  cluster_tips <- unlist(lapply(X = child_clusters, FUN = get_cluster_tips))
  cluster_tip_nodes <- unlist(lapply(X = child_clusters, FUN = get_cluster_tip_nodes))
  return(data.frame(
    "size" = sum(cluster_sizes),
    "ch_mrca" = node_data$node,
    "ch_tmrca" = node_data$date,
    "ch_tmrca_CI" = I(list(node_data$CI_date)),
    "tips" = paste0(cluster_tips, collapse = ", "),
    "tip_nodes" = paste0(cluster_tip_nodes, collapse = ", "),
    "n_nonfocal_subclades" = n_nonfocal_subclades,
    "n_basal_nonswiss_clades" = n_basal_nonswiss_clades))
}

finish_cluster <- function(child_cluster, child_node_data, node_data, verbose = F) {
  # Given a parent node and a child swiss cluster, return cluster information for the child.
  if (verbose) {
    print(paste("...finishing child", child_node_data$node, "of node", node_data$node)) 
  }
  if (child_cluster$size > 1)  {
    cluster_finished <- data.frame(
      "size" = child_cluster$size,
      "ch_mrca" = child_cluster$ch_mrca,
      "ch_tmrca" = child_cluster$ch_tmrca,
      "ch_tmrca_CI" = I(child_cluster$ch_tmrca_CI),
      "tips" = child_cluster$tips,
      "tip_nodes" = child_cluster$tip_nodes,
      "n_nonfocal_subclades" = child_cluster$n_nonfocal_subclades,
      "n_basal_nonswiss_clades" = 0,
      "foreign_tmrca" = node_data$date,
      "foreign_tmrca_CI" = I(list(c(node_data$CI_date))),
      "foreign_mrca" = node_data$node)
    return(cluster_finished)
  }
  cluster_finished <- data.frame(
    "size" = child_cluster$size, 
    "foreign_tmrca" = node_data$date,
    "foreign_tmrca_CI" =  I(list(c(node_data$CI_date))),
    "foreign_mrca" = node_data$node,
    "ch_tmrca" = child_cluster$ch_tmrca,
    "ch_tmrca_CI" = NA,
    "ch_mrca" = child_cluster$ch_mrca,
    "tips" = child_cluster$tips,
    "tip_nodes" = child_cluster$tip_nodes,
    "n_nonfocal_subclades" = child_cluster$n_nonfocal_subclades,
    "n_basal_nonswiss_clades" = 0)
  return(cluster_finished)
}

characterize_ch_tmrcas <- function(node_data, n_tips, tree_data, max_nonfocal_subclades, max_consecutive_budding_nonfocal_subclades, verbose = F) {
  # Traverse tree recursively from node described by node_data, find and summarize swiss clusters
  node_data <- as.list(node_data)
  if (verbose) {
    print(paste("considering node", node_data$node))
  }
  if (is_tip(node_data$node, n_tips) & !is.na(node_data$country_recoded) & node_data$country_recoded == FOCAL_LOC) {
    # Base case 1: if tip is swiss return cluster of size 1
    return(initiate_swiss_cluster(node_data, verbose = verbose))
  } else if (is_tip(node_data$node, n_tips)) {
    # Base case 2: if tip is not swiss return nonfocal clade
    return(get_dummy_cluster_for_nonfocal_node(node_data = node_data, verbose = verbose))
  } else {
    # Recursively call function on all children of node
    child_node_data <- get_child_node_data(node_data$node, tree_data)
    child_clusters <- apply(
      X = child_node_data, FUN = characterize_ch_tmrcas, MARGIN = 1, n_tips, 
      tree_data, max_nonfocal_subclades, max_consecutive_budding_nonfocal_subclades, verbose = verbose)  # a list
  
    child_swissness <- unlist(lapply(X = child_clusters, FUN = is_swiss_cluster))
    n_nonfocal_subclades_in_children <- unlist(lapply(X = child_clusters, FUN = get_n_nonfocal_subclades))
    n_pendant_subclades_to_swiss_children <- unlist(lapply(X = child_clusters, FUN = get_n_basal_nonswiss_clades))
    n_nonswiss_children <- sum(!child_swissness)  
    
    # Get number of pendant non-swiss subclades to consider for cluster merge
    n_pendant_subclades <- max(n_pendant_subclades_to_swiss_children) + n_nonswiss_children
    if (verbose) {
      print(paste("Checking to merge clusters at node", 
                  node_data$node, "which has", 
                  max(n_pendant_subclades_to_swiss_children), 
                  "from max to a swiss child and", n_nonswiss_children, 
                  "from non-swiss children"))
    }
    
    # Recursive case 1: no children are swiss, return 1 non-focal subclade
    if (!any(child_swissness)) {
      return(get_dummy_cluster_for_nonfocal_node(node_data = node_data, verbose = verbose))
    }
    # Recursive cases 2: merging swiss clusters wouldn't violate cluster criteria, return merged swiss cluster
    # Except if the node is the root, which I don't allow to be Swiss
    if (sum(n_nonfocal_subclades_in_children) <= max_nonfocal_subclades & node_data$node != n_tips + 1) {  # merging children doesn't violate MAX_NONFOCAL_SUBCLADES
      if (n_pendant_subclades <= max_consecutive_budding_nonfocal_subclades) {  # merging children doesn't violate MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES
        if (sum(child_swissness) > 1) {
          merged_cluster <- merge_child_swiss_clusters(
            child_clusters[which(child_swissness)], 
            node_data, 
            n_nonfocal_subclades = sum(n_nonfocal_subclades_in_children), 
            n_basal_nonswiss_clades = n_pendant_subclades, 
            verbose = verbose)
        } else {
          merged_cluster <- propogate_swiss_cluster(
            child_clusters[which(child_swissness)], 
            node_data, 
            n_nonfocal_subclades = sum(n_nonfocal_subclades_in_children), 
            n_basal_nonswiss_clades = n_pendant_subclades, 
            verbose = verbose)
        }
        return(merged_cluster)
      }
    }
    # Recursive cases 3a: if MAX_NONFOCAL_SUBCLADES unviolated, merge swiss descendents of node to one cluster
    if (POLYTOMIES_ARE_SWISS) {
      if (sum(child_swissness) > 1) {
        merged_cluster <- merge_child_swiss_clusters(
          child_clusters[which(child_swissness)], 
          node_data, 
          n_nonfocal_subclades = sum(n_nonfocal_subclades_in_children), 
          n_basal_nonswiss_clades = n_pendant_subclades, 
          verbose = verbose)
      } else {
        merged_cluster <- propogate_swiss_cluster(
          child_clusters[which(child_swissness)], 
          node_data, 
          n_nonfocal_subclades = sum(n_nonfocal_subclades_in_children), 
          n_basal_nonswiss_clades = n_pendant_subclades, 
          verbose = verbose)
      }
      cluster_info <- finish_cluster(
        merged_cluster,
        node_data,
        node_data,
        verbose = verbose)
      clusters <<- rbind(clusters, cluster_info, stringsAsFactors = F)
    } else {
      # Recursive cases 3b: finish daughter clusters individually
      for (child_idx in 1:length(child_clusters)) {
        if (child_swissness[child_idx]) {
          child_cluster_info <- finish_cluster(
            child_clusters[[child_idx]], 
            child_node_data[child_idx, ], 
            node_data, 
            verbose = verbose)
          clusters <<- rbind(clusters, child_cluster_info, stringsAsFactors = F)
        }
      }
    }
    # Return node as 1 non-focal subclade
    return(get_dummy_cluster_for_nonfocal_node(node_data = node_data, verbose = verbose))  # collapse all descendents to one export
  }
}
  
# Delete zero-length internal branches 
root_node <- unname(unlist(tree_data %>% filter(parent == node) %>% select(node)))
tree_data <- delete_internal_zero_branches(tree_data = tree_data, root_node = root_node)

# Initialize results data frames
CLUSTER_INFO <- c(
  "foreign_mrca", 
  "size", 
  "foreign_tmrca", 
  "foreign_tmrca_CI", 
  "ch_tmrca", 
  "ch_tmrca_CI", 
  "ch_mrca",
  "tips",
  "tip_nodes",
  "n_nonfocal_subclades",
  "n_basal_nonswiss_clades")
clusters <- as.data.frame(matrix(nrow = 0, ncol = length(CLUSTER_INFO)))
names(clusters) <- CLUSTER_INFO

root_data <- tree_data[tree_data$node == root_node, ]
root_cluster <- characterize_ch_tmrcas(
  node_data = root_data, 
  n_tips = n_tips,
  tree_data = tree_data, 
  max_nonfocal_subclades = MAX_NONFOCAL_SUBCLADES,
  max_consecutive_budding_nonfocal_subclades = MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES,
  verbose = IS_VERBOSE)

# Add information on how clusters were generated
clusters$max_nonfocal_subclades <- MAX_NONFOCAL_SUBCLADES
clusters$max_consecutive_budding_nonfocal_subclades <- MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES

# Date ranges at root are list of depth 2, which is problematic, so correct
# TODO: not sure why this happens
depth <- function(this,thisdepth=0){
  if(!is.list(this)){
    return(thisdepth)
  }else{
    return(max(unlist(lapply(this,depth,thisdepth=thisdepth+1))))    
  }
}

for (i in 1:nrow(clusters)) {
  if(depth(clusters[i, "foreign_tmrca_CI"]) == 2) {
    clusters[[i, "foreign_tmrca_CI"]] <- I(unlist(clusters[i, "foreign_tmrca_CI"]))
  }
  if(depth(clusters[i, "ch_tmrca_CI"]) == 2) {
    clusters[[i, "ch_tmrca_CI"]] <- I(unlist(clusters[i, "ch_tmrca_CI"]))
  }
}

clusters$foreign_tmrca_CI <- as.character(clusters$foreign_tmrca_CI)
clusters$ch_tmrca_CI <- as.character(clusters$ch_tmrca_CI)

# Write out cluster data 
write.table(
  x = clusters,
  file = paste(OUTDIR, paste(PREFIX, "_clusters.txt", sep = ""), sep = "/"),
  quote = F, row.names = F, col.names = T, sep = "\t")

