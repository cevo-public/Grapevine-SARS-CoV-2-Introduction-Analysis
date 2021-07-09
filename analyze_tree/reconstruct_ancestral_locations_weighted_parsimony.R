# This script infers ancestral locations on the tree. For this, we use only tips 
# in the context dataset and ignore tips in the similarity dataset. 
# Otherwise, over-sequenced locations would show up too often as sources. 
# Our procedure is as follows: we begin by labelling the MRCAs of focal
# transmission chains and all subtending nodes in the tree (excepting exported 
# clades) as focal. Focal country is excluded as a possible location for all
# other nodes. Given these constraints, we calculate the parsimony score for 
# each possible location at each remaining node. Among the possible locations, 
# we also include a “dummy” location to record the maximum parsimony score. So, 
# in a first step we calculate parsimony scores up the tree. In a second step, 
# we covert the parsimony scores at each node to location weights by taking the 
# difference to the “dummy” score. Thus, locations with no support in the 
# subtending tree are given a weight 0 and supported locations are weighted by 
# the number of location changes they prevent.

# This script takes a tree (possibly with polytomies) and calculates parsimony
# scores (number of downstream changes necessary) for all possible loctions at ancestral nodes. 
# It's a modification of the Sankoff algorithm described by Joe Felsenstein:
# https://evolution.gs.washington.edu/gs541/2010/lecture1.pdf
# The modification is that all the nodes that are ancestors of tips classified 
# as in a focal cluster (as defined in the chains produced by
# pick_focal_transmission_chains.R) are treated like tips in that their locations are fixed.
# Ancestral state reconstruction then proceeds from nodes that are not ancestral 
# to tips classified as in a focal cluster without regard to to the rest of the
# tree. Additionally, at binary splits we keep only the weights for minimally-
# scoring (most parsimonious if this were the root) locations and set all other
# weights to Inf. At polytomies, we keep all weights because we want to keep some
# uncertainty reflecting that we don't know how to resolve the polytomy. 
# 
# To then assign locations based on these weights, we calculate a score as the 
# number of changes avoided. This is normalized across all locations. 

suppressMessages(suppressWarnings(require(ape)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(treeio)))
suppressMessages(suppressWarnings(require(ggplot2)))
suppressMessages(suppressWarnings(require(ggtree)))
suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(magrittr)))

# tree_file <- "/Users/nadeaus/Repos/cov-swiss-phylogenetics/results_all/validation/jan-dec_-005_max_sampling_-25_travel_-5_sim_context-sf_111_travel-wt/tmp/lsd/B.1.timetree.nex"
# metadata_file <- "/Users/nadeaus/Repos/cov-swiss-phylogenetics/results_all/validation/jan-dec_-005_max_sampling_-25_travel_-5_sim_context-sf_111_travel-wt/tmp/alignments/B.1_metadata.tsv"
# chains_file <- "/Users/nadeaus/Repos/cov-swiss-phylogenetics/results_all/validation/jan-dec_-005_max_sampling_-25_travel_-5_sim_context-sf_111_travel-wt/tmp/chains/B.1_m_3_p_1_s_F_chains.txt"
# s <- F
# outdir <- "~/Downloads"
# verbose <- T
# plot_tree <- T
# write_scores <- F
# prefix <- paste("test_s_", s, sep = "")
# 
parser <- argparse::ArgumentParser()
parser$add_argument("--tree", type="character")
parser$add_argument("--metadata", type="character")
parser$add_argument("--contextmetadata", type = "character")
parser$add_argument("--chains", type="character")
parser$add_argument("--outdir", type="character")
parser$add_argument("--prefix", type="character")
parser$add_argument("--plottree", action="store_true")
parser$add_argument("--verbose", action="store_true")
parser$add_argument("--focalcountry", type="character")
parser$add_argument("--writescores", action="store_true")

args <- parser$parse_args()

tree_file <- args$tree
metadata_file <- args$metadata
chains_file <- args$chains
outdir <- args$outdir
prefix <- args$prefix
verbose <- args$verbose
plot_tree <- args$plottree
focal_country <- args$focalcountry
write_scores <- args$writescores

source("utility_functions.R")
system(command = paste("mkdir -p", outdir))

# Load data
tree <- treeio::read.beast(file = tree_file)
tree_data <- tidytree::as_tibble(tree)
metadata <- read.table(file = metadata_file, stringsAsFactors = F, sep = "\t", header = T, quote = "")
chains <- read.delim(file = chains_file, stringsAsFactors = F)

# Make sure all tips present in metadata
tree_ids <- tree_data$label[!is.na(tree_data$label)]  # tips only
check_all_tips_in_metadata(
  metadata_ids = metadata$tree_label, 
  tree_ids = tree_ids)

# Get list of nodes in travel set (ignore similarity set, outgroup)
tree_data <- merge(x = tree_data, y = metadata %>% select(-c(date)), 
                   by.x = "label", by.y = "tree_label", all.x = T)
context_node_set <- unlist(tree_data %>% 
  filter(travel_context) %>%
  select(node))
print(paste(tree_file, "has", length(context_node_set), "travel context sequences"))
focal_node_set <- unlist(tree_data %>%
  filter(focal_sequence) %>%
  select(node))
print(paste(tree_file, "has", length(focal_node_set), "focal sequences"))

# Initialize tree data with ancestral state information and cluster status
tree_data$asr_loc <- tree_data$iso_country
tree_data[tree_data$node %in% chains$ch_mrca, "asr_loc"] <- focal_country

# Delete internal zero-length branches as they screw up parsimony ASR
n_tips <- length(tree_ids)
root_node <- unname(unlist(tree_data %>% filter(parent == node) %>% select(node)))
tree_data <- delete_internal_zero_branches(tree_data = tree_data, root_node = root_node)

set_tip_score <- function(node, loc_to_idx_mapping, verbose) {
  # If tip is in context seq set, initialize score for tip loc to be 0, all others Inf;
  # Otherwise initialize score for all tips to be 0 (no information provided); 
  # Mark node as visited
  node_data <- tree_data[tree_data$node == node, ]
  if (node %in% c(context_node_set, focal_node_set)) {
    if (verbose) {
      print(paste("visiting context tip node", node))
    }
    tip_loc <- node_data$asr_loc
    node_data$candidate_asr_scores[[1]][loc_to_idx_mapping[tip_loc]] <- 0
  } else {
    if (verbose) {
      print(paste("visiting priority tip node", node, "which will be given 0 score for all locs"))
    }
    node_data$candidate_asr_scores[[1]] <- rep(0, length(loc_to_idx_mapping))
  }
  node_data$is_visited <- T
  tree_data[tree_data$node == node, ] <<- node_data
}

set_mrca_score <- function(node, loc_to_idx_mapping, verbose) {
  # Initialize score for assigned ASR loc to be 0, all others Inf; mark node as visited
  if (verbose) {
    print(paste("visiting focal MRCA node", node))
  }
  node_data <- tree_data[tree_data$node == node, ]
  tip_loc <- node_data$asr_loc
  node_data$candidate_asr_scores[[1]][loc_to_idx_mapping[tip_loc]] <- 0
  node_data$is_visited <- T
  tree_data[tree_data$node == node, ] <<- node_data
}

is_ch_mrca <- function(node, n_tips, focal_country) {
  # Return if ancestral location has been fixed to focal country in the tree data
  if (is.na(tree_data[tree_data$node == node, "asr_loc"]) | is_tip(node = node, n_tips = n_tips)) {
    return(F)
  }
  return(tree_data[tree_data$node == node, "asr_loc"] == focal_country)
}

visit_children <- function(node, verbose) {
  # Visit all children, calculating their scores w/ additional recursive calls if necessary
  child_node_data <- get_child_node_data(node = node, tree_data)
  if (verbose) {
    print(paste("node", node, "has children", paste0(child_node_data$node, collapse = ", ")))
  }
  # Visit unvisited children
  if (!(all(child_node_data$is_visited))) {
    for (i in which(!(child_node_data$is_visited))) {
      if (verbose) {
        print(paste("node", child_node_data[i, "node"], "is an unvisited child of", node, "visiting it"))
      }
      get_asr_scores(node = child_node_data[i, "node"], n_tips = n_tips, loc_to_idx_mapping = loc_to_idx_mapping, verbose = verbose)
    }
  }
  child_node_data <- get_child_node_data(node = node, tree_data)  # fetch again now that child tree_data has been updated
  return(child_node_data)
}

visit_nonfocal_node <- function(node, candidate_asr_locs, child_node_data, verbose) {
  if (verbose) {
    print(paste("calculating scores for node", node))
  }
  # initialize data structures to hold scores & pointers for each possible ancestral state
  S_a <- rep(NA, length(candidate_asr_locs))
  S_a_pointers <- rep(NA, length(candidate_asr_locs))
  # for every possible ancestral state
  for (i in 1:length(candidate_asr_locs)) {
    S_a_i <- 0  # initialize a score for this ancestral state
    S_a_i_pointers <- rep(NA, nrow(child_node_data))  # initialize a list for pointers to which child state(s) this score came from
    # for every child node, add that child's contribution to the score
    for (c in 1:nrow(child_node_data)) {
      child_data_c <- child_node_data[c, ]
      contribution_c <- Inf  # initialize a minimum contribution 
      pointers_c <- c()
      # update minimum, considering every possible state at the child node
      for (j in 1:length(candidate_asr_locs)) {
        c_ij <- ifelse(test = i == j, yes = 0, no = 1)  # cost for any transition is 1
        if (i == loc_to_idx_mapping[focal_country]) {
          c_ij <- Inf  # cost for transition from focal country is prohibitily large (focal country not allowed as location at non-focal nodes)
        }
        S_c_j <- child_data_c$candidate_asr_scores[[1]][j]
        if (c_ij + S_c_j < contribution_c) {
          contribution_c <- c_ij + S_c_j
          pointers_c <- j
        } else if (c_ij + S_c_j == contribution_c) {
          pointers_c <- c(pointers_c, j)
        }
      }
      # add minimum of child's possible contributions to ancestral state score
      S_a_i <- S_a_i + contribution_c
      # track which child state(s) this score came from
      S_a_i_pointers[c] <- I(list(pointers_c))
    }
    if (verbose) {
      print(paste("score for", candidate_asr_locs[i], "at node", node, "is", S_a_i, "which comes from pointers to"))
      print(S_a_i_pointers)
    }
    S_a[i] <- S_a_i
    S_a_pointers[i] <- I(list(S_a_i_pointers))
  }
  
  # Record this in the data structure
  tree_data[tree_data$node == node, "candidate_asr_scores"][[1]] <<- I(list(S_a))
  tree_data[tree_data$node == node, "candidate_asr_pointers"][[1]] <<- I(list(S_a_pointers))
  tree_data[tree_data$node == node, "is_visited"] <<- T
}

assign_node_as_focal <- function(node, cluster_tips, loc_to_idx_mapping, verbose, focal_country) {
  if (verbose) {
    print(paste("node", node, "is to be assigned as focal. Checking whether its children contain any of", paste0(cluster_tips, collapse = ", ")))
  }
  tree_data[tree_data$node == node, "asr_loc"] <<- focal_country
  child_node_data <- get_child_node_data(node = node, tree_data = tree_data)
  for (i in 1:nrow(child_node_data)) {
    child_node <- child_node_data[i, "node"]
    if (is_tip(node = child_node, n_tips = n_tips)) {
      # Base case: node is tip, initialize score list 
      set_tip_score(node = child_node, loc_to_idx_mapping = loc_to_idx_mapping, verbose = verbose)
    } else {
      if (verbose) {
        print(paste("checking child", child_node, "for descendent cluster tips"))
      }
      tips_under_child <- get_tips_under_node(node = child_node, tree_data = tree_data, n_tips)
      if (any(tips_under_child %in% cluster_tips)) {
        if (verbose) {
          print(paste("one of the tips under", child_node, "is found to be in the cluster"))
        }
        assign_node_as_focal(node = child_node, cluster_tips, loc_to_idx_mapping, verbose, focal_country)
      } else {
        if (verbose) {
          print("no cluster tips descendent, treating as new tree root")
        }
        get_asr_scores(node = child_node, n_tips = n_tips, loc_to_idx_mapping = loc_to_idx_mapping, verbose = verbose)
      }
    }
  }
  set_mrca_score(node = node, loc_to_idx_mapping = loc_to_idx_mapping, verbose = verbose)  # set MRCA location to focal country
}

get_cluster_tips <- function(ch_mrca, chains) {
  tips_str <- as.character(chains[chains$ch_mrca == ch_mrca, "tip_nodes"])
  tips <- strsplit(x = tips_str, split = ", ")[[1]]
  return(as.numeric(tips))
}

# Recursively traverse tree, assigning scores and pointers as you go
get_asr_scores <- function(node, n_tips, loc_to_idx_mapping, verbose) {
  if (verbose) {
    print(paste("considering node", node, "for score"))
  }
  if (is_tip(node = node, n_tips = n_tips)) {
    # Base case: node is tip, initialize score list 
    set_tip_score(node = node, loc_to_idx_mapping = loc_to_idx_mapping, verbose = verbose)
    return()
  } 
  if (is_ch_mrca(node = node, n_tips = n_tips, focal_country = focal_country)) {
    # Recursive case: node is a focal mrca, continue on to score children as if they are the root of a tree
    cluster_tips = get_cluster_tips(ch_mrca = node, chains = chains)
    assign_node_as_focal(node, cluster_tips = cluster_tips, loc_to_idx_mapping = loc_to_idx_mapping, verbose, focal_country)
    return()
  }
  # Recursive case, visit unvisited children
  if (verbose) {
    print(paste("visiting children of node", node))
  }
  child_node_data <- visit_children(node = node, verbose)
  if (verbose) {
    print(paste("visiting node", node))
  }
  visit_nonfocal_node(node = node, candidate_asr_locs = names(loc_to_idx_mapping), child_node_data = child_node_data, verbose)
  return()
}

# Make treedata structure with locations, scores, and pointers to update during DP algorithm
n_tips <- length(tree_ids)
candidate_asr_locs <- c(unique(metadata$iso_country), "dummy_loc")  # include dummy location to keep track of the max. # changes required from any node
loc_to_idx_mapping <- 1:length(candidate_asr_locs)
names(loc_to_idx_mapping) <- candidate_asr_locs

if (!(focal_country %in% names(loc_to_idx_mapping))) {  # if no focal tips in tree, add focal country to candidate locations anyways so code doesn't break
  candidate_asr_locs <- c(candidate_asr_locs, focal_country)
  loc_to_idx_mapping[[focal_country]] <- length(loc_to_idx_mapping) + 1
}

tree_data$candidate_asr_locs <- I(list(candidate_asr_locs))
tree_data$candidate_asr_scores <- I(list(rep(Inf, length(candidate_asr_locs))))
tree_data$candidate_asr_pointers <- NA

# Set all nodes to unvisited
tree_data$is_visited <- F

# Bottom-up tree traversal to assign parsimony scores and collect pointers to which child assignments generated the scores
print("Getting potential ancestral state reconstruction scores")
get_asr_scores(
  node = n_tips + 1, n_tips = n_tips, 
  loc_to_idx_mapping = loc_to_idx_mapping, verbose = verbose)

# Write out scores
if (write_scores) {
  tree_data_to_print <- tree_data 
  tree_data_to_print$candidate_asr_scores <- unlist(lapply(
    FUN = collapse_list, X = tree_data_to_print$candidate_asr_scores))
  tree_data_to_print$candidate_asr_locs <- unlist(lapply(
    FUN = collapse_list, X = tree_data_to_print$candidate_asr_locs))
  tree_data_to_print$candidate_asr_pointers <- unlist(lapply(
    FUN = collapse_list, X = tree_data_to_print$candidate_asr_pointers))
  tree_data_to_print <- tree_data_to_print %>%
    select(-c(CI_date, CI_height))
  
  write.table(
    x = tree_data_to_print,
    file = paste(outdir, paste(prefix, "_tree_data_with_scores.txt", sep = ""), sep = "/"),
    quote = F, row.names = F, col.names = T, sep = "\t")
}

# Generate data frame with one row per node, 1 column per candidate ancestral location,
# where the entries are the ASR weights normalized to 1 for each node
temp <- tree_data[c("node", "candidate_asr_scores")]

temp$candidate_asr_scores <- unlist(lapply(
  X = temp$candidate_asr_scores, FUN = paste0, collapse = ", "))
asr_pie_data <- temp %>%
  tidyr::separate(
    col = candidate_asr_scores,
    sep = ", ",
    into = candidate_asr_locs,
    remove = T)
asr_pie_data <- apply(X = asr_pie_data, FUN = as.numeric, MARGIN = 2)

get_max_non_inf_asr_score <- function(asr_scores_for_node) {
  return(asr_scores_for_node[loc_to_idx_mapping[["dummy_loc"]]])
}

# Get max # changes (= parsimony score for a dummy location that has no tips in the tree)
max_changes <- apply(
  X = asr_pie_data[, 2:ncol(asr_pie_data)],
  FUN = get_max_non_inf_asr_score,
  MARGIN = 1)

# Assign ASR score for each valid location to be: max_changes - n_changes
# Where min_changes is the # descendent focal MRCAs
min_changes <- unlist(lapply(
  X = asr_pie_data[, 1],
  FUN = get_n_descendent_focal_mrcas,
  tree_data = tree_data, chains = chains, n_tips = n_tips))

calc_score <- function(max_changes, min_changes, num_changes) {
  # If the location is not allowed, return score 0
  if (num_changes == Inf) {
    return(0)
  # If all neighbor sequences are priority, then the max # changes = min # changes and return 0 for all locations
  } else if (max_changes <= min_changes) {
    return(0)
  # If the location is perfect, return score 1 (encompasses the case where focal country is enforced)
  } else if (num_changes == 0) {
    return(1)
  # If the location the best possible, return score 1 (encompasses the case in which the loc is the only valid option)
  } else if (num_changes == min_changes) {
    return(max_changes - num_changes)
  # Otherwise, return score as # changes avoided
  } else {
    return(max_changes - num_changes)
  }
}

for (i in 1:nrow(asr_pie_data)) {
  scores <- unlist(lapply(
    X = asr_pie_data[i, 2:ncol(asr_pie_data)],
    FUN = calc_score,
    max_changes = max_changes[i],
    min_changes = min_changes[i]))
  asr_pie_data[i, 2:ncol(asr_pie_data)] <- scores
}

# Disregard dummy location, normalize weights to one
# If no locs have positive score, e.g. a node with only priority nodes descending,
# all locs will get an NA weight (the node's location is unassignable)
asr_pie_data[, ncol(asr_pie_data)] <- 0  # set score for dummy loc to be always 0
asr_pie_data[, 2:ncol(asr_pie_data)] <- sweep(
  x = asr_pie_data[, 2:ncol(asr_pie_data)], 
  MARGIN = 1, 
  FUN = "/",
  rowSums(asr_pie_data[, 2:ncol(asr_pie_data)]))
asr_pie_data <- as.data.frame(asr_pie_data)

# Write out tree data with ancestral node location weights
asr_pie_data_to_merge <- asr_pie_data
colnames(asr_pie_data_to_merge)[2:ncol(asr_pie_data)] <- 
  paste(colnames(asr_pie_data)[2:ncol(asr_pie_data)], "loc_weight", sep = "_")

tree_data_to_print <- merge(
  x = tree_data %>% select(-c(
    "asr_loc", "candidate_asr_locs", "candidate_asr_scores", 
    "candidate_asr_pointers")),
  y = asr_pie_data_to_merge %>% select(-c("dummy_loc_loc_weight")), 
  by = "node")

# Have to coerce list elements to character for printing
tree_data_to_print[c("CI_date", "CI_height")] <- apply(
  X = tree_data_to_print[c("CI_date", "CI_height")],
  FUN = as.character,
  MARGIN = 2)

write.table(
  x = tree_data_to_print,
  file = paste(outdir, paste(prefix, "_tree_data_with_asr.txt", sep = ""), sep = "/"),
  row.names = F, col.names = T, quote = F, sep = "\t")

# Plot the tree with ASR pie charts at nodes
if (plot_tree) {
  tree_data <- tree_data %>% mutate(
    node_annotation = case_when(
      node %in% chains$ch_mrca ~ "(Focal)",
      T ~ ""),
    is_context  = case_when(
      node %in% context_node_set ~ "Context",
      T ~ "Priority"))
  tree_2 <- full_join(tree, tree_data)
  pies <- ggtree::nodepie(asr_pie_data, cols = 2:ncol(asr_pie_data))
  p <- ggtree(tr = tree_2) +
    geom_inset(insets = pies) +
    geom_tiplab(aes(
      label = ifelse(
        test = is.na(node_annotation),
        yes = "",
        no = paste(strain, ": ", country, sep = "")),
      color = country,
      alpha = is_context),
      hjust = -0.01,
      size = 1) +
    scale_alpha_manual(
      values = c("Priority" = 0.5, "Context" = 1),
      name = "Sample type",
      na.translate = F) +
    scale_color_discrete(guide = F) +
    theme(legend.position = "bottom") +
    geom_nodelab(aes(
      label = ifelse(
        test = is.na(node_annotation),
        yes = "",
        no = paste(node, node_annotation))),
      hjust = 0.01)
  
  ggsave(
    file = paste(outdir, paste(prefix, "tree_with_asr.pdf", sep = "_"), sep = "/"), 
    plot = p,
    dpi = 400
    )
}

