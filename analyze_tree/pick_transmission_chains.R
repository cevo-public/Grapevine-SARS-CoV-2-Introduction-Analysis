# This script takes a tree and tip location information and returns 
# chain assignments based on a user-defined max number of monophyletic foreign
# sub-clades (exports) and a max number of foreign clades budding from the spine
# of a focal country transmission chain in a row.
# Note: the tree is NOT allowed to be one big chain (at the highest 
# level of aggregation, clades descending from the root will be separate 
# transmission chains). 
# Note: When calculating the number of foreign clades that bud off from a focal
# transmission chain in a row, foreign descendents of polytomies are assumed to 
# be one clade. The rationale is that a single export from a focal chain seems
# more likely than a series of many exports in a row.
# Note: When calculating the largest plausible chains (polytomies are
# assumed to be focal, l = T), at large polytomies focal descendents are
# aggregated into transmission chains in size order. This doesn't guarantee the
# largest possible transmission chains are formed but should promote the 
# formation of large transmission chains.  

suppressMessages(suppressWarnings(require(treeio)))
suppressMessages(suppressWarnings(require(tidytree)))
suppressMessages(suppressWarnings(require(ggplot2)))
suppressMessages(suppressWarnings(require(ggtree)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(purrr)))
suppressMessages(suppressWarnings(require(readr)))

# tree <- "/Users/nadeaus/Repos/cov-swiss-phylogenetics/results_all/2021-07-08_nzl/tmp/lsd/A.2.timetree.nex"
# metadata <- "~/Downloads/metadata_quote.tsv"
# outdir <- "~/Downloads"
# verbose <- T
# m <- 3
# p <- 1
# l <- T
# prefix <- paste("test_l_", l, sep = "")
# focalcountry <- "NZL"
# dont_plot_tree <- F

parser <- argparse::ArgumentParser()
parser$add_argument("--tree", type="character")
parser$add_argument("--metadata", type="character")
parser$add_argument("--outdir", type="character")
parser$add_argument("--maxtotalsubclades", type="double")
parser$add_argument("--maxconsecutivesubclades", type="double")
parser$add_argument("--prefix", type="character")
parser$add_argument("--largestchains", action = "store_true")
parser$add_argument("--focalcountry", type="character")
parser$add_argument("--dontplot", action = "store_true")

args <- parser$parse_args()

tree <- args$tree
metadata_file <- args$metadata
outdir <- args$outdir
m <- args$maxtotalsubclades
p <- args$maxconsecutivesubclades
prefix <- args$prefix
l <- args$largestchains
dont_plot_tree <- args$dontplot
focal_country <- args$focalcountry
verbose <- F

system(command = paste("mkdir -p", outdir))
source("utility_functions.R")
source("analyze_tree/functions.R")

# Load data
tree <- read.beast(file = tree)
metadata <- readr::read_delim(file = metadata_file, delim = "\t", escape_backslash=TRUE, escape_double = FALSE, col_types = cols())
tree_data <- as_tibble(tree)
tree_data <- merge(
  x = tree_data, y = metadata %>% select(-c(date)),
  by.x = "label", by.y = "tree_label", all.x = T)

# Get transmission chains
chains <- pick_chains(
  tree_data = tree_data, 
  m = m, 
  p = p, 
  l = l,
  focal_country = focal_country,
  verbose = verbose
)

# Clean transmission chains (messy handling of some list depth, type issues)
chains <- clean_chains(
  chains = chains
)

# Write out chain data 
write.table(
  x = chains,
  file = paste(outdir, paste(prefix, "_chains.txt", sep = ""), sep = "/"),
  quote = F, row.names = F, col.names = T, sep = "\t")

# Plot transmission chains on tree
if (!dont_plot_tree) {
  tree_plot <- plot_chains_on_tree(
    chains = chains
  )
  ggsave(
    file = paste(outdir, paste(prefix, "_chains.png", sep = ""), sep = "/"),
    plot = tree_plot)
}