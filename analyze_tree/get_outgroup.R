# This script outputs the outgroup file for LSD program (provide via -g option).

require(argparse)
require(ape)
require(ggplot2)
require(gridExtra)
require(ggtree)
require(dplyr)

TREEFILE <- "/Users/nadeaus/Repos/grapevine/dont_commit/aug_1/iqtree/B.1.1.1.treefile"
OUTGROUP_STR <- "EPI_ISL_402125, EPI_ISL_406798"
METADATA <- "/Users/nadeaus/Repos/grapevine/dont_commit/aug_1/alignments/B.1.1.1_metadata.csv"
OUTDIR <- "~/Downloads"
PLOT <- F

# parser <- argparse::ArgumentParser()
# parser$add_argument("--treefile", type="character")
# parser$add_argument("--metadata", type="character")
# parser$add_argument("--outgroup", type="character")
# parser$add_argument("--outdir", type="character")
# parser$add_argument("--plottree", action="store_true", default=F)
# 
# args <- parser$parse_args()
# 
# TREEFILE <- args$treefile
# METADATA <- args$metadata
# OUTGROUP_STR <- args$outgroup
# OUTDIR <- args$outdir
# PLOT <- args$plottree

TREENAME <- strsplit(x = TREEFILE, split = "/")[[1]]
TREENAME <- TREENAME[length(TREENAME)]
TREENAME <- gsub(x = TREENAME, pattern = ".treefile", replacement = "")

print(paste("TREENAME:", TREENAME))
print(paste("OUTGROUP_STR:", OUTGROUP_STR))
print(paste("OUTDIR:", OUTDIR))

# Load data
print("loading data")
tree <- ape::read.tree(file = TREEFILE)
metadata <- read.csv(file = METADATA)
outgroup_defining_tips <- strsplit(x = OUTGROUP_STR, split = ", ")[[1]]
tree_data <- tidytree::as_tibble(tree)

# If outgroup tips are monophyletic, root on branch between the 2 outgroup tips
outgroup_defining_tips_parents <- unlist(tree_data %>%
  mutate(outgroup_defining_tip = label %in% !! outgroup_defining_tips) %>%
  filter(outgroup_defining_tip) %>%
  select(parent))
if (length(unique(outgroup_defining_tips_parents)) == 1) {
  print(paste("Outgroup-defining tips", 
              paste0(outgroup_defining_tips, collapse = ", "),
              "are monophyletic. Rooting on the branch between them",
              "by specifying one as the outgroup."))
  outgroup_tips <- outgroup_defining_tips[1]
  rooted_tree <- ape::root(phy = tree, outgroup = outgroup_tips)
} else {
  stop("Tree", TREEFILE, "does not have outgroup-defining tips as monophyletic. Don't know how to root the tree.")
}

# Plot tree rooted with defined outgroup
if (PLOT) {
  tree_figname <- paste(OUTDIR, paste(TREENAME, "_rooted.png", sep = ""), sep = "/")
  unrooted_tree_figname <- paste(OUTDIR, paste(TREENAME, "_unrooted.png", sep = ""), sep = "/")
  png(file = unrooted_tree_figname, height = 15, width = 15, units = "in", res = 250)
  ggtree(tr = tree, layout="fan", open.angle = 120) + 
    geom_tippoint(aes(color = case_when(
      label %in% outgroup_defining_tips ~ "outgroup_defining_seq",
      T ~ "context or swiss"))) + 
    scale_color_discrete(name = "tip type")
  dev.off()
  png(file = tree_figname, height = 15, width = 15, units = "in", res = 250)
  ggtree(tr = rooted_tree, layout="fan", open.angle = 120) + 
    geom_tippoint(aes(color = case_when(
      label %in% outgroup_defining_tips ~ "outgroup_defining_seq",
      label %in% outgroup_tips ~ "in outgroup",
      T ~ "in ingroup"))) + 
    scale_color_discrete(name = "tip type")
  dev.off()
}

# Write out outgroup seqs for IQ-TREE
outgroup_str <- paste("'", paste0(outgroup_tips, collapse = ","), "'", sep = "")
filename <- paste(OUTDIR, paste(TREENAME, ".outgroup.txt", sep = ""), sep = "/")
writeLines(
  text = outgroup_str,
  con = file(filename))
