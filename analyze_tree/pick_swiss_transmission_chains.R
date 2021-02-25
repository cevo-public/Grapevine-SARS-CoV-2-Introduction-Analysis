# This script takes a tree and tip location information and returns 
# chain assignments based on a user-defined max number of monophyletic foreign
# sub-clades (exports) and a max number of foreign clades budding from any internal
# branch. Note: the tree is NOT allowed to be one big chain (at the highest 
# level of aggregation, clades descending from the root will be chains). Also
# note that the number of budding (pendant) subclades is calculated as the maximum 
# # pendant (budding subclades) to any swiss child plus the number of non-swiss
# siblings. This is conservative because it assumes polytomies would be resolved
# such that the maximum # non-swiss clades bud off in a row.

require(treeio)
require(tidytree)
require(ggplot2)
require(ggtree)
require(dplyr)

# tree <- "/Users/nadeaus/Repos/grapevine/dont_commit/test/tmp/lsd/B.1.1.277.timetree.nex"
# metadata <- "/Users/nadeaus/Repos/grapevine/dont_commit/test/tmp/alignments/B.1.1.277_metadata.csv"
# outdir <- "~/Downloads"
# verbose <- T
# m <- 3
# p <- 1
# s <- T
# prefix <- paste("test_s_", s, sep = "")

parser <- argparse::ArgumentParser()
parser$add_argument("--tree", type="character")
parser$add_argument("--metadata", type="character")
parser$add_argument("--outdir", type="character")
parser$add_argument("--maxtotalsubclades", type="double")
parser$add_argument("--maxconsecutivesubclades", type="double")
parser$add_argument("--prefix", type="character")
parser$add_argument("--polytomiesareswiss", action = "store_true")

args <- parser$parse_args()

tree <- args$tree
metadata <- args$metadata
outdir <- args$outdir
m <- args$maxtotalsubclades
p <- args$maxconsecutivesubclades
prefix <- args$prefix
s <- args$polytomiesareswiss
verbose <- F

system(command = paste("mkdir -p", outdir))
source("utility_functions.R")
source("analyze_tree/functions.R")

# Load data
tree <- treeio::read.beast(file = tree)
metadata <- read.csv(file = metadata, stringsAsFactors = F)
tree_data <- tidytree::as_tibble(tree)
tree_data <- merge(
  x = tree_data, y = metadata %>% select(-c(date)),
  by.x = "label", by.y = "tree_label", all.x = T)

# Get transmission chains
chains <- pick_chains(
  tree_data = tree_data, 
  m = m, 
  p = p, 
  s = s, 
  verbose = verbose
)

# Clean transmission chains (messy handling of some list depth, type issues)
chains <- clean_chains(
  chains = chains
)

# Plot transmission chains on tree
tree_plot <- plot_chains_on_tree(
  chains = chains
)

# Write out chain data 
write.table(
  x = chains,
  file = paste(outdir, paste(prefix, "_chains.txt", sep = ""), sep = "/"),
  quote = F, row.names = F, col.names = T, sep = "\t")
ggsave(
  file = paste(outdir, paste(prefix, "_chains.png", sep = ""), sep = "/"),
  plot = tree_plot)
