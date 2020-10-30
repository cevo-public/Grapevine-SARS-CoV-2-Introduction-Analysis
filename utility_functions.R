# Utility functions for project, mostly having to do with manipulating data and 
# trees.

is_tip <- function(node, n_tips) {
  # Return T if the node is a tip based on the default node #ing and # tips
  return(as.numeric(node) <= n_tips)
}

get_child_node_data <- function(node, tree_data) {
  # Take node # and tree data (tidytree treedata format, must have cols node, parent)
  # Return data frame in treedata format for all nodes w/ node as parent 
  child_node_data <- tree_data[tree_data$parent == node, ]
  # If node is root, do not return root as child
  child_node_data <- child_node_data[child_node_data$parent != child_node_data$node, ]  
  return(child_node_data)
}

get_n_descendent_swiss_mrcas <- function(node, tree_data, cluster_data, n_tips) {
  # Returns 1 if the node is a Swiss MRCA, otherwise returns the # Swiss MRCAs
  # descending from the node. Counts only the uppermost MRCA if there additional
  # Swiss MRCAs are nested under a Swiss MRCA higher in the tree.
  if (node %in% cluster_data$ch_mrca) {
    return(1)
  } else if (node <= n_tips) {
    return(0)
  } else {
    child_node_data <- get_child_node_data(node = node, tree_data = tree_data)
  }
  n_descendent_swiss_mrcas <- sum(unlist(lapply(
    X = child_node_data$node,
    FUN = get_n_descendent_swiss_mrcas,
    tree_data = tree_data, 
    cluster_data = cluster_data, 
    n_tips = n_tips)))
  return(n_descendent_swiss_mrcas)
}

collapse_list <- function(list) {
  return(paste0(list, collapse = ", "))
}

get_n_context_children <- function(node, tree_data, context_node_set, n_tips) {
  child_node_data <- get_child_node_data(node, tree_data)
  child_has_context_seqs <- unlist(lapply(
    X = child_node_data$node, 
    FUN = get_node_has_context_seqs,
    tree_data = tree_data,
    context_node_set = context_node_set,
    n_tips = n_tips))
  return(sum(child_has_context_seqs))
}

recode_colnames <- function(col_name) {
  # Recode ASR locations (format "First.Second_loc_weight") to match country_recoded names
  if (grepl(x = col_name, pattern = "_loc_weight")) {
    col_name_new <- gsub(x = col_name, pattern = "_loc_weight", replacement = "")
    col_name_new <- gsub(x = col_name_new, pattern = "\\.", replacement = " ")
    return(col_name_new)
  }
  return(col_name)
}

get_node_has_context_seqs <- function(tree_data, node, context_node_set, n_tips) {
  # Return whether the node or any of its descendents are in the context node set
  is_tip <- node <= n_tips
  has_context <- node %in% context_node_set 
  while (!(is_tip | has_context)) {
    child_node_data <- get_child_node_data(node = node, tree_data = tree_data)
    has_context <- any(unlist(lapply(
      X = child_node_data$node, FUN = get_node_has_context_seqs, 
      tree_data = tree_data, context_node_set = context_node_set, n_tips = n_tips)))
  }
  return(has_context)
}

is_swiss_cluster <- function(cluster) {
  # If cluster object has null size, it's a dummy cluster 
  return(!(is.na(cluster$size)))
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

get_n_nonfocal_subclades <- function(cluster) {
  return(cluster$n_nonfocal_subclades)
}

get_n_basal_nonswiss_clades <- function(cluster) {
  return(cluster$n_basal_nonswiss_clades)
}

check_all_tips_in_metadata <- function(metadata_ids, tree_ids) {
  # Check all samples in tree present in metadata
  missing_in_metadata <- tree_ids[!(tree_ids %in% metadata_ids)]
  if (length(missing_in_metadata) > 0) {
    message <- paste(length(missing_in_metadata), "seqs missing from metadata:")
    stop(paste(message, paste0(head(missing_in_metadata, 5), collapse = ", ")))
  }
}

get_CI_date_min <- function(CI_date) {
  # tidytree::get_data returns node date ranges as list. Get the lower bound.
  if (is.null(CI_date)) {
    return(NA)
  }
  return(CI_date[1])
}

get_CI_date_max <- function(CI_date) {
  # tidytree::get_data returns node date ranges as list. Get the upper bound.
  if (is.null(CI_date)) {
    return(NA)
  }
  return(CI_date[2])
}

get_tips_under_node <- function(node, tree_data, n_tips) {
  # Base case: node is tip, return tip
  if (is_tip(node = node, n_tips = n_tips)) {
    return(node)
  }
  # Recursive case: paste together tips under children
  child_node_data <- get_child_node_data(node = node, tree_data = tree_data)
  return(unlist(lapply(FUN = get_tips_under_node, X = child_node_data$node, tree_data = tree_data, n_tips = n_tips)))
}

delete_internal_zero_branches <- function(tree_data, root_node) {
  # Delete zero-length internal branches 
  nodes_to_delete <- c()
  for (i in 1:nrow(tree_data)) {
    bl <- tree_data$branch.length[i]
    node_idx <- i
    while (bl == 0 & tree_data$node[node_idx] > root_node) {  # if a zero length branch b/w 2 internal nodes
      parent_idx <- which(tree_data$node == tree_data$parent[node_idx])
      bl <- tree_data[parent_idx, "branch.length"]
      node_idx <- parent_idx
      node <- tree_data[i, "node"]
      parent_node <- tree_data[parent_idx, "node"]
      print(paste("Internal node", node, "is a daughter of a zero length branch"))
      print(paste("it will be deleted and its children will have parent node", parent_node))
      tree_data[tree_data$parent == node, "parent"] <- parent_node
      nodes_to_delete <- c(nodes_to_delete, node)
    }
  }
  tree_data <- tree_data[!(tree_data$node %in% nodes_to_delete), ]
  return(tree_data)
}

get_locs_of_descendent_tips <- function(node, tree_data) {
  # Return unique locations amongst descendent tips of the node
  child_data <- get_child_node_data(node = node, tree_data = tree_data)
  return(unique(child_data$country_recoded[!is.na(child_data$country_recoded)]))
}

# Mapping from unambiguous canton codes to Nextstrain division names
kanton_code_dict <- list(
  "AG" = "Aargau",
  "BE" = "Bern",
  "BL" = "Basel-Land",
  "BS" = "Basel-Stadt",
  "FR" = "Fribourg",
  "GE" = "Geneva",
  "VD" = "Vaud",
  "GL" = "Glarus",
  "GR" = "Graubünden",
  "JU" = "Jura",
  "LU" = "Lucerne",
  "SG" = "Sankt Gallen",
  "SH" = "Schaffhausen",
  "SO" = "Solothurn",
  "SZ" = "Schwyz",
  "TG" = "Thurgau",
  "TI" = "Ticino",
  "UR" = "Uri",
  "VS" = "Valais",
  "ZG" = "Zug",
  "ZH" = "Zürich")

get_ns_division_from_canton_code <- function(canton_code) {
  if (!(canton_code %in% names(kanton_code_dict))) {
    stop(paste("Canton", canton_code, "not known."))
  }
  return(kanton_code_dict[canton_code])
}

get_CI_min <- function(CI) {
  CI_str <- gsub(
    x = gsub(x = CI, pattern = "c\\(", replacement = ""),
    pattern = "\\)", replacement = "")
  return(strsplit(CI_str, split = ", ")[[1]][1])
}
get_CI_max <- function(CI) {
  CI_str <- gsub(
    x = gsub(x = CI, pattern = "c\\(", replacement = ""),
    pattern = "\\)", replacement = "")
  return(strsplit(CI_str, split = ", ")[[1]][2])
}

get_days_from_190101 <- function(date) {
  # A hack because scatterpie doesn't seem to play nicely with dates
  # Translate dates into a numberic value (days since Jan. 1st 2019)
  return(as.numeric(as.Date(date) - as.Date("2019-01-01")))
}

