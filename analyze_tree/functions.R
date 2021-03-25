#' @return Summarize characteristics of swiss transmission chains below input node.
#' @param tree_data Tree data in tidytree treedata format.
#' @param m maximum number of international subclades per transmission chain.
#' @param p maximum number of international subclaes budding from an internal branch
#' of a transmission chain.
#' @param s boolean indicating whether to assume swiss descendents of a 
#' polytomy are from the same introduction or not.
#' @return chains Dataframe with fields describing each transmission chain.
pick_chains <- function(tree_data, m, p, s, verbose) {
  chains <<- initialize_chains()
  
  root_node <- unname(unlist(
    tree_data %>% filter(parent == node) %>% select(node)))
  tree_data <- delete_internal_zero_branches(
    tree_data = tree_data, root_node = root_node)
  root_data <- tree_data[tree_data$node == root_node, ]
  
  n_tips <- length(tree_data$label[!is.na(tree_data$label)])
  
  root_cluster <- get_cluster_from_node(
    node_data = root_data, 
    n_tips = n_tips,
    tree_data = tree_data, 
    m = m, p = p, s = s,
    verbose = verbose)  # this recursive function updates the chains object
  
  return(chains)
}

#' Initialize Dataframe with fields to describe each transmission chain.
initialize_chains <- function() {
  chain_info <- c(
    "foreign_mrca", 
    "size", 
    "foreign_tmrca", 
    "foreign_tmrca_CI", 
    "ch_tmrca", 
    "ch_tmrca_CI", 
    "ch_mrca",
    "tips",
    "tip_nodes",
    "n_intl_subclades",
    "n_basal_intl_clades")
  chains <- as.data.frame(matrix(nrow = 0, ncol = length(chain_info)))
  names(chains) <- chain_info
  return(chains)
}


#' Summarize cluster information for the node; a recursive function. 
#' @param node_data Node data in tidytree treedata format for a single node.
#' @return cluster Dataframe with fields describing the cluster at the node.
get_cluster_from_node <- function(
  node_data, n_tips, tree_data, m, p, s, verbose = F
) {
  # Traverse tree recursively from node described by node_data, find and summarize swiss clusters
  node_data <- as.list(node_data)
  node <- node_data$node
  iso_country <- node_data$iso_country
  
  if (verbose) {
    print(paste("considering node", node))
  }
  
  # Base case 1: if tip is swiss return cluster of size 1
  if (is_tip(node, n_tips) & !is.na(iso_country) & iso_country == "CHE") {
    swiss_cluster <- initiate_swiss_cluster(node_data, verbose = verbose)
    return(swiss_cluster)
    
    # Base case 2: if tip is not swiss return intl clade
  } else if (is_tip(node, n_tips)) { 
    intl_cluster <- get_intl_cluster(node_data = node_data, verbose)
    return(intl_cluster)
    
    # Recursive case: get cluster information from all child nodes
  } else {
    child_node_data <- get_child_node_data(node, tree_data)
    child_clusters <- apply(
      X = child_node_data, 
      FUN = get_cluster_from_node, 
      MARGIN = 1, 
      n_tips, tree_data, m, p, s, verbose = verbose)  # a list
    child_swissness <- unlist(lapply(
      X = child_clusters, 
      FUN = is_swiss_cluster))
    n_intl_subclades_in_children <- unlist(lapply(
      X = child_clusters, 
      FUN = get_n_intl_subclades))
    n_pendant_subclades_to_swiss_children <- unlist(lapply(
      X = child_clusters, 
      FUN = get_n_basal_intl_clades))
    child_sizes <- unlist(lapply(
      X = child_clusters,
      FUN = get_cluster_size))
    n_intl_children <- sum(!child_swissness)  
    n_pendant_subclades <- case_when(
      sum(child_swissness) == 1 ~ max(n_pendant_subclades_to_swiss_children) + 1,
      T ~ max(n_pendant_subclades_to_swiss_children)) # "pendant" means branching off in a row. The spine of the clade is proposed to be swiss, so we only add a pendant subclade if all the siblings are non-swiss
    
    if (verbose) {
      print(paste("Checking whether to merge clusters at node", 
                  node, "which would introduce a Swiss cluser with", 
                  n_pendant_subclades, "pendant subclades."))
    }
    
    # Recursive case 1: no children are swiss, return 1 international subclade
    if (!any(child_swissness)) {
      intl_cluster <- get_intl_cluster(node_data = node_data, verbose = verbose)
      return(intl_cluster)
    }
    
    m_condition <- sum(n_intl_subclades_in_children) <= m
    p_condition <- n_pendant_subclades <= p
    root_condition <- node != n_tips + 1  # root not allowed to be swiss
    
    # Recursive case 2: merging all child clusters to same transmission chain doesn't violate transmission chain criteria
    if (m_condition & p_condition & root_condition) {  
      if (sum(child_swissness) > 1) {
        merged_cluster <- merge_child_swiss_clusters(
          child_clusters[which(child_swissness)], 
          node_data, 
          n_intl_subclades = sum(n_intl_subclades_in_children), 
          n_basal_intl_clades = n_pendant_subclades, 
          verbose = verbose)
      } else {
        merged_cluster <- propogate_swiss_cluster(
          child_clusters[which(child_swissness)], 
          node_data, 
          n_intl_subclades = sum(n_intl_subclades_in_children), 
          n_basal_intl_clades = n_pendant_subclades, 
          verbose = verbose)
      }
      return(merged_cluster)
      
    # Recursive cases 3a: merging all child clusters to same transmission chain violates transmission chain criteria,
    # but we assume polytomies is resolved such that as many swiss clusters as possible are merged
    # Don't need to consider p condition here because not merging in any intl siblings, only combining swiss clades
    } else if (s & root_condition & sum(child_swissness) > 1) {
      child_merge_list <- get_child_idxs_to_merge(
        ms = n_intl_subclades_in_children,
        sizes = child_sizes,
        swissness = child_swissness,
        m = m, verbose = verbose)
      for (i in 1:length(child_merge_list)) {
        child_idxs_to_merge <- child_merge_list[[i]]
        if (verbose) {
          print("Merging children")
          print(child_idxs_to_merge)
        }
        if (length(child_idxs_to_merge) > 1) {
          merged_cluster <- merge_child_swiss_clusters(
            child_clusters[child_idxs_to_merge],
            node_data,
            n_intl_subclades = sum(n_intl_subclades_in_children[child_idxs_to_merge]),
            n_basal_intl_clades = n_pendant_subclades,
            verbose = verbose)
          cluster <- finish_cluster(
            merged_cluster,
            child_node_data = node_data,
            node_data = node_data,
            verbose = verbose)
          chains <<- rbind(chains, cluster, stringsAsFactors = F)
        } else {
          child_cluster <- finish_cluster(
            child_clusters[[child_idxs_to_merge]],
            child_node_data = child_node_data[child_idxs_to_merge, ],
            node_data = node_data,
            verbose = verbose)
          chains <<- rbind(chains, child_cluster, stringsAsFactors = F)
        }
      }

    # Recursive cases 3b: merging all child clusters to same transmission chain violates transmission chain criteria,
    # and we assume polytomy is not in Switzerland: finish each child of polytomy seperately.
    } else {
      for (child_idx in (1:length(child_clusters))[child_swissness]) {
        child_cluster <- finish_cluster(
          child_clusters[[child_idx]],
          child_node_data = child_node_data[child_idx, ],
          node_data = node_data,
          verbose = verbose)
        chains <<- rbind(chains, child_cluster, stringsAsFactors = F)
      }
    }
    intl_cluster <- get_intl_cluster(node_data = node_data, verbose)  # collapse all descendents to one export
    return(intl_cluster)  
  }
}

#' At a polytomy, generate the largest Swiss cluster possible by combining 
#' daughter clades in order of size (# Swiss descendents) as long as transmission chain
#' criteria m is satisfied. When no more daughter clades can be added, 
#' finish the set to merge and start another merge, working from largest to smallest size.
#' @param ms A vector of the number intl subclades in each child.
#' @param sizes Sizes (# Swiss descendents) in each child.
#' @param swissness Boolean vector whether each child is Swiss or not.
#' @param m Maximum # non-Swiss subclades in a transmission chain.
#' @return A vector of indexes in list child_clusters that can be merged.
get_child_idxs_to_merge <- function(
  ms, sizes, swissness, m, verbose
) {
  nodes <- 1:length(ms)
  merge_data <- data.frame(
    node = nodes[swissness],
    size = sizes[swissness],
    m = ms[swissness]) %>%
    arrange(desc(size))
  
  m_tally <- 0
  merge_idx <- "A"
  merge_list <- list()
  for (i in 1:nrow(merge_data)) {
    m_tally_proposed <- m_tally + merge_data[i, "m"]
    if (m_tally_proposed <= m) {
      m_tally <- m_tally_proposed
    } else {
      m_tally <-  merge_data[i, "m"]
      merge_idx <- intToUtf8(utf8ToInt(merge_idx) + 1)
      merge_list[[merge_idx]] <- c()
    }
    merge_list[[merge_idx]] <- c(merge_list[[merge_idx]], merge_data[i, "node"])
  }
  if (verbose) {
    print(merge_data)
    print("Results in the following merges:")
    print(merge_list)
  }
  return(merge_list)
}

#' Create a placeholder cluster representing an international tip or clade
#' @return Dataframe with fields describing the cluster
get_intl_cluster <- function(node_data, verbose = F) {
  if (verbose) {
    print(paste("...node", node_data$node,  "is a international clade"))
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
    "n_basal_intl_clades" = 0,
    "n_intl_subclades" = 1))
}

#' Create a cluster representing a swiss tip
#' @param node_data Node data in tidytree treedata format for a single swiss tip
#' @return Dataframe with fields describing the cluster
initiate_swiss_cluster <- function(node_data, verbose = F) {
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
    "n_basal_intl_clades" = 0,
    "n_intl_subclades" = 0))
}

#' Update swiss cluster that has only international siblings: n_intl_subclades,
#' n_basal_intl_clades, and if relevant foreign mrca fields are updated.
#' @return Dataframe with fields describing the cluster. 
propogate_swiss_cluster <- function(
  child_clusters, node_data, n_intl_subclades, n_basal_intl_clades, 
  verbose = F
) {
  if (verbose) {
    print(paste(
      "...only one child of",
      node_data$node, "is swiss but criteria not exceeded by siblings,",
      "propogating the cluster up with", 
      n_intl_subclades, "intl_subclades and",
      n_basal_intl_clades, "basal_intl_clades"))
  }
  child_clusters[[1]][["n_basal_intl_clades"]] <- n_basal_intl_clades
  child_clusters[[1]][["n_intl_subclades"]] <- n_intl_subclades
  if (n_basal_intl_clades == 1) {
    # If connecting swiss cluster to first international node, update foreign mrca
    if (verbose) {
      print(paste("updating foreign mrca to", node_data$node))
    }
    child_clusters[[1]][["foreign_mrca"]] <- node_data$node
    child_clusters[[1]][["foreign_tmrca"]] <- node_data$date
    child_clusters[[1]][["foreign_tmrca_CI"]] <- I(list(node_data$CI_date))
  }
  return(child_clusters[[1]])
}

#' Merge sibling swiss clusters. 
#' @return  Dataframe with fields describing the merged cluster.
merge_child_swiss_clusters <- function(
  child_clusters, node_data, n_intl_subclades, n_basal_intl_clades, 
  verbose = F
) {
  if (verbose) {
    print(paste(
      "...some children of",
      node_data$node, "are swiss & merge criteria satisfied,",
      "returning merged cluster with", 
      n_intl_subclades, "intl_subclades and",
      n_basal_intl_clades, "basal_intl_clades"))
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
    "n_intl_subclades" = n_intl_subclades,
    "n_basal_intl_clades" = n_basal_intl_clades))
}

#' Re-format cluster information once the cluster is finished (now termed a 
#' transmission chain).
#' @return Dataframe with fields describing the transmission chain.
finish_cluster <- function(
  child_cluster, child_node_data, node_data, verbose = F
) {
  if (verbose) {
    print(paste(
      "...finishing child", child_node_data$node, "of node", node_data$node)) 
  }
  if (child_cluster$n_basal_intl_clades > 0) {
    foreign_mrca_node_data <- child_node_data
  } else {
    foreign_mrca_node_data <- node_data
  }
  if (child_cluster$size > 1)  {
    cluster_finished <- data.frame(
      "size" = child_cluster$size,
      "ch_mrca" = child_cluster$ch_mrca,
      "ch_tmrca" = child_cluster$ch_tmrca,
      "ch_tmrca_CI" = I(child_cluster$ch_tmrca_CI),
      "tips" = child_cluster$tips,
      "tip_nodes" = child_cluster$tip_nodes,
      "n_intl_subclades" = child_cluster$n_intl_subclades,
      "n_basal_intl_clades" = 0,
      "foreign_tmrca" = foreign_mrca_node_data$date,
      "foreign_tmrca_CI" = I(list(c(foreign_mrca_node_data$CI_date))),
      "foreign_mrca" = foreign_mrca_node_data$node)
    return(cluster_finished)
  }
  cluster_finished <- data.frame(
    "size" = child_cluster$size, 
    "foreign_tmrca" = foreign_mrca_node_data$date,
    "foreign_tmrca_CI" =  I(list(c(foreign_mrca_node_data$CI_date))),
    "foreign_mrca" = foreign_mrca_node_data$node,
    "ch_tmrca" = child_cluster$ch_tmrca,
    "ch_tmrca_CI" = NA,
    "ch_mrca" = child_cluster$ch_mrca,
    "tips" = child_cluster$tips,
    "tip_nodes" = child_cluster$tip_nodes,
    "n_intl_subclades" = child_cluster$n_intl_subclades,
    "n_basal_intl_clades" = 0)
  return(cluster_finished)
}

get_n_intl_subclades <- function(cluster) {
  return(cluster$n_intl_subclades)
}

get_n_basal_intl_clades <- function(cluster) {
  return(cluster$n_basal_intl_clades)
}

is_swiss_cluster <- function(cluster) {
  return(!(is.na(cluster$size)))  # If cluster object has null size, it's a dummy cluster 
}

get_cluster_size <- function(cluster) {
  return(cluster$size)
}

get_cluster_tips <- function(cluster) {
  return(cluster$tips)
}

get_cluster_tip_nodes <- function(cluster) {
  return(cluster$tip_nodes)
}

#' Correct chains because date ranges at root are list of depth 2, which is problematic.
#' Not sure why this happens.
clean_chains <- function(chains) {
  for (i in 1:nrow(chains)) {
    if(get_depth(chains[i, "foreign_tmrca_CI"]) == 2) {
      chains[[i, "foreign_tmrca_CI"]] <- I(unlist(chains[i, "foreign_tmrca_CI"]))
    }
    if(get_depth(chains[i, "ch_tmrca_CI"]) == 2) {
      chains[[i, "ch_tmrca_CI"]] <- I(unlist(chains[i, "ch_tmrca_CI"]))
    }
  }
  
  chains$foreign_tmrca_CI <- as.character(chains$foreign_tmrca_CI)
  chains$ch_tmrca_CI <- as.character(chains$ch_tmrca_CI)
  return(chains)
}

#' Return depth of list.
get_depth <- function(my_list, thisdepth = 0){
  if (!is.list(my_list)) {
    return(thisdepth)
  } else {
    return(max(unlist(lapply(my_list, get_depth, thisdepth = thisdepth + 1))))    
  }
}

plot_chains_on_tree <- function(chains, max_chains_to_plot = 4) {
  chain_data_by_tip <- chains %>%
    mutate(chain_idx = 1:n()) %>%
    tidyr::separate(
      col = tips,
      sep = ", ",
      into = paste("tip", 1:max(chains$size), sep = "_"),
      remove = F,
      fill = "right") %>% 
    tidyr::pivot_longer(
      cols = paste("tip", 1:max(chains$size), sep = "_"),
      values_to = "tip",
      names_to = "tip_idx",
      names_prefix = "tip_") %>%
    filter(!(is.na(tip)))
  
  chain_idxs <- unique(chain_data_by_tip$chain_idx)
  if (length(chains) > max_chains_to_plot) {
    chain_data_by_tip <- chain_data_by_tip %>%
      mutate(
        chain_idx = as.character(chain_idx),
        chain_idx = case_when(
        chain_idx %in% as.character(1:max_chains_to_plot) ~ chain_idx,
        T ~ "other (max colors exceeded)"
      ))
    chain_idxs <- unique(chain_data_by_tip$chain_idx)
  }
  tree_data_2 <- merge(x = tree_data, y = chain_data_by_tip, 
                       by.x = "label", by.y = "tip", all.x = T)
  my.cols <- RColorBrewer::brewer.pal(length(chain_idxs), "Dark2")
  names(my.cols) <- chain_idxs
  tree_plot <- ggtree::ggtree(tr = tree) %<+% tree_data_2 + 
    geom_tippoint(aes(
      color = case_when(
        !is.na(chain_idx) ~ as.character(chain_idx),
        T ~ "foreign sequence"))) + 
    geom_nodelab(aes(label = node)) +
    geom_tiplab(aes(label = node)) + 
    scale_color_manual(name = "transmission chain",
                       values = c("foreign sequence" = "grey", my.cols))
  return(tree_plot)
}

# tree <- "/Users/nadeaus/Repos/grapevine/dont_commit/grapevine_for_testing_scripts/tmp/lsd/B.1.1.74.timetree.nex"
# metadata <- "/Users/nadeaus/Repos/grapevine/dont_commit/grapevine_for_testing_scripts/tmp/alignments/B.1.1.74_metadata.csv"
# 
# # Load data
# tree <- treeio::read.beast(file = tree)
# metadata <- read.csv(file = metadata, stringsAsFactors = F)
# tree_data <- tidytree::as_tibble(tree)
# tree_data <- merge(
#   x = tree_data, y = metadata %>% select(-c(date)),
#   by.x = "label", by.y = "tree_label", all.x = T)
# tree_data[tree_data$node %in% c(27, 13, 14), "iso_country"] <- "CHE"
# 
# ggtree(tr = tree) %<+% tree_data + 
#   geom_tippoint(aes(color = iso_country == "CHE"))
# 
# source("utility_functions.R")
# 
# # Get transmission chains
# chains <- pick_chains(
#   tree_data = tree_data,
#   m = 3,
#   p = 1,
#   s = T,
#   verbose = T
# )
# 
# plot_chains_on_tree(chains = chains, max_chains_to_plot = 10)

