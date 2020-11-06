# Get at how many different lineages seeded the second wave
# Merge swiss seqs until DATE_THRESHOLD, then count how many different lineages there 
# still are and how many samples were produced by each one.

# DATE_THRESHOLD <- "2020-07-01"
# TREE <- "/Users/nadeaus/Downloads//tmp/lsd/rep_1_n_sim_1000_n_imports_padded_0.timetree.nex"
# CLADE_DATA <- "/Users/nadeaus/Downloads//clades/swiss_alignment_filtered2_masked_oneline_clades.tsv"
# METADATA <- "/Users/nadeaus/Downloads//tmp/alignments/rep_1_n_sim_1000_n_imports_padded_0/rep_1_n_sim_1000_n_imports_padded_0_tree_metadata.txt"
# UTILITY_FUNCTIONS <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/automated//utility_functions.R"

require(dplyr)
require(ape)

parser <- argparse::ArgumentParser()
parser$add_argument("--datethreshold", type="character")
parser$add_argument("--tree", type="character")
parser$add_argument("--clades", type="character")
parser$add_argument("--metadata", type="character")
parser$add_argument("--prefixdata", type="character")
parser$add_argument("--outdir", type="character")
parser$add_argument("--utilityfunctions", type="character")

args <- parser$parse_args()

DATE_THRESHOLD <- args$datethreshold
TREE <- args$tree
CLADE_DATA <- args$clades
METADATA <- args$metadata
OUTDIR <- args$outdir
UTILITY_FUNCTIONS <- args$utilityfunctions
PREFIX_DATA <- args$prefixdata

source(UTILITY_FUNCTIONS)
system(command = paste("mkdir -p", OUTDIR))

tree <- treeio::read.beast(file = TREE)
tree_data <- tidytree::as_tibble(tree)
metadata <- read.table(file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
clades <- read.delim(file = CLADE_DATA, sep = "\t")

metadata <- merge(
  x = metadata, y = clades[c("seqName", "clade")],
  by.x = "old_strain", by.y = "seqName", all.x = T)

metadata <- metadata %>% 
  mutate(
    is_viollier = case_when(
      originating_lab == "Viollier AG" ~ "viollier",
      T ~ "other")) # count only viollier samples
    
tree_data <- merge(
  x = tree_data, y = metadata[c("strain", "country_recoded", "is_viollier", "clade")],
  by.x = "label", by.y = "strain", all.x = T)

n_tips <- sum(!is.na(tree_data$label))
tree_data <- delete_internal_zero_branches(tree_data = tree_data, root_node = n_tips + 1, verbose = F)

# Get tips clustered into lineages present on July 1

# Recursively traverse tree from the root, as soon as you hit a node 
# with estimated date after July 1st, if some tips under the node are swiss
# record the node and tips under it
get_tips_by_threshold_lineage <- function(node, tree_data, n_tips, tree, date_threshold) {
  # Base case 1: node is tip and is it's own July 1st lineage
  node_date <- as.Date(tree_data[tree_data$node == node, "date"])
  tips_below <- get_tips_under_node(node = node, tree_data = tree_data, n_tips = n_tips)
  if (is_tip(node, n_tips = n_tips) & node_date > as.Date(date_threshold)) {
    if (tree_data[tree_data$node == node, "is_viollier"] == "viollier") {
      return(data.frame(
        threshold_node = node, 
        swiss_tip = node, 
        all_tips_below_threshold_node = 1))
    }
  } 
  # Base case 2: node is after date threshold, return tips under it as a July 1st lineage
  if (node_date > as.Date(date_threshold)) {
    viollier_tips_below <- tree_data[
      tree_data$node %in% tips_below & tree_data$is_viollier == "viollier", "node"] 
    if (length(viollier_tips_below) > 1) {
      return(data.frame(
        threshold_node = node, 
        swiss_tip = viollier_tips_below, 
        all_tips_below_threshold_node = length(tips_below)))
    }
  }
  # Recursive case: node before date threshold, check daughters
  child_nodes <- get_child_node_data(node = node, tree_data = tree_data)$node
  for (child in child_nodes) {
    child_threshold_lineages <- get_tips_by_threshold_lineage(
      node = child, tree_data, n_tips, tree, date_threshold)
    threshold_lineages <<- rbind(threshold_lineages, child_threshold_lineages)
  }
}

threshold_lineages <- data.frame(
  threshold_node = c(), 
  swiss_tip = c(),
  all_tips_below_threshold_node = c())
get_tips_by_threshold_lineage(
  node = n_tips + 1, 
  tree_data = tree_data, tree = tree, n_tips = n_tips,
  date_threshold = DATE_THRESHOLD)

threshold_lineages <- merge(
  x = threshold_lineages, y = tree_data[c("node", "clade")], 
  by.x = "swiss_tip", by.y = "node", all.x = T)

threshold_lineage_summary <-threshold_lineages %>%
  group_by(threshold_node, all_tips_below_threshold_node) %>%
  summarize(n_viollier_tips_below = n(),
            n_20A = sum(clade == "20A"),
            n_20B = sum(clade == "20B"),
            n_20C = sum(clade == "20C"),
            n_19A = sum(clade == "19A"),
            n_other = sum(!(clade %in% c("20A", "20B", "20C", "19A")))) %>%
  arrange(desc(n_viollier_tips_below))

if (any(threshold_lineage_summary$n_other > 0)) {
  stop("some unknown clades")
}

assign_clade_to_lineage <- function(n_20A, n_20B, n_20C, n_19A) {
  if (n_20A != 0 & n_20B == 0 & n_20C == 0) {
    return("20A")
  } else if (n_20A == 0 & n_20B != 0 & n_20C == 0) {
    return("20B")
  } else if (n_20A == 0 & n_20B == 0 & n_20C != 0) {
    return("20C")
  } else if (n_20A == 0 & n_20B == 0 & n_20C == 0 & n_19A != 0) {
    return("19A")
  } else {
    print(paste("n_20A:", n_20A, "n_20B:", n_20B, "n_20C:", n_20C, "n_19A:", n_19A))
    clade <- c("20A", "20B", "20C", "n_19A")[which(c(n_20A, n_20B, n_20C, n_19A) == max(c(n_20A, n_20B, n_20C, n_19A)))][1]
    print(paste("is assigned clade", clade))
    return(clade)
  }
}

threshold_lineage_summary$clade <- mapply(
  FUN = assign_clade_to_lineage,
  n_20A = threshold_lineage_summary$n_20A,
  n_20B = threshold_lineage_summary$n_20B,
  n_20C = threshold_lineage_summary$n_20C,
  n_19A = threshold_lineage_summary$n_19A)

write.table(
  x = threshold_lineage_summary,
  file = paste(OUTDIR, "/", PREFIX_DATA, "_lineages_crossing_", DATE_THRESHOLD, ".txt", sep = ""),
  row.names = F, col.names = T, sep = "\t", quote = F)
