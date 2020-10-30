# This script outputs the outgroup file for LSD program (provide via -g option).

require(argparse)
require(ape)
require(ggplot2)
require(gridExtra)
require(ggtree)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
PREFIX <- "rep_1_n_sim_1000_n_imports_padded_0"
TREEFILE <- paste(WORKDIR, "/ml_trees/", PREFIX, ".treefile", sep = "")
OUTGROUP_STR <- "Wuhan_Hu-1_2019|EPI_ISL_402125|2019-12-26, Wuhan_WH01_2019|EPI_ISL_406798|2019-12-26"
METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX, "/", PREFIX, "_tree_metadata.txt", sep = "")
OUTDIR <- "~/Downloads"

parser <- argparse::ArgumentParser()
parser$add_argument("--treefile", type="character")
parser$add_argument("--metadata", type="character")
parser$add_argument("--outgroup", type="character")
parser$add_argument("--outdir", type="character")

args <- parser$parse_args()

TREEFILE <- args$treefile
METADATA <- args$metadata
OUTGROUP_STR <- args$outgroup
OUTDIR <- args$outdir

TREENAME <- strsplit(x = TREEFILE, split = "/")[[1]]
TREENAME <- TREENAME[length(TREENAME)]

print(paste("TREENAME:", TREENAME))
print(paste("OUTGROUP_STR:", OUTGROUP_STR))
print(paste("OUTDIR:", OUTDIR))

# Load data
print("loading data")
tree <- ape::read.tree(file = TREEFILE)
metadata <- read.table(file = METADATA, sep = "\t", quote = "", fill = T, header = T)
outgroup_defining_tips <- strsplit(x = OUTGROUP_STR, split = ", ")[[1]]

# Work my way up from the first outgroup-defining tip to find the node defining the 
# largest possible ingroup. All tips not in this ingroup are in the outgroup.
# This is a messy work-around b/c I couldn't find a way to get the MRCA of 
# 2 nodes on an unrooted tree (e.g. ape would just return the tree root node to me,
# which puts all tips in the outgroup).
print("finding largest possible ingroup")
tree_data <- tidytree::as_tibble(tree)
node <- unlist(tree_data[!is.na(tree_data$label) & 
                    tree_data$label == outgroup_defining_tips[1], "node"][[1]])
outgroup_seqs_found <- F
biggest_sib_node_size <- 0
while(!outgroup_seqs_found) {
  parent_node <- unname(unlist(tree_data[tree_data$node == node, "parent"]))
  siblings <- tree_data[tree_data$parent == parent_node, "node"]$node
  for (sibling in siblings[siblings != node & siblings > length(tree$tip.label)]) {
    sibling_tips <- ape::extract.clade(phy = tree, node = sibling)$tip.label
    if (outgroup_defining_tips[2] %in% sibling_tips) {
      outgroup_seqs_found <- T
    } else {
      if (length(sibling_tips) > biggest_sib_node_size) {
        biggest_sib_node_size <- length(sibling_tips)
        biggest_sib_node <- sibling
        biggest_sib_node_tips <- sibling_tips
        print("Assigning outgroup tips")
        outgroup_tips <- tree$tip.label[!(tree$tip.label %in% sibling_tips)]
      } 
    }
  }
  node <- parent_node
}

# if searching up from one outgroup-defining tip hit the root before finding non-outgroup siblings, start from the other
if (!exists("outgroup_tips")) {
  node <- unlist(tree_data[!is.na(tree_data$label) & 
                             tree_data$label == outgroup_defining_tips[2], "node"][[1]])
  outgroup_seqs_found <- F
  biggest_sib_node_size <- 0
  while(!outgroup_seqs_found) {
    parent_node <- unname(unlist(tree_data[tree_data$node == node, "parent"]))
    siblings <- tree_data[tree_data$parent == parent_node, "node"]$node
    for (sibling in siblings[siblings != node & siblings > length(tree$tip.label)]) {
      sibling_tips <- ape::extract.clade(phy = tree, node = sibling)$tip.label
      if (outgroup_defining_tips[1] %in% sibling_tips) {
        outgroup_seqs_found <- T
      } else {
        if (length(sibling_tips) > biggest_sib_node_size) {
          biggest_sib_node_size <- length(sibling_tips)
          biggest_sib_node <- sibling
          biggest_sib_node_tips <- sibling_tips
          outgroup_tips <- tree$tip.label[!(tree$tip.label %in% sibling_tips)]
        } 
      }
    }
    node <- parent_node
  }
}

# If the ingroup is smaller than the outgroup, make the ingroup the outgroup
if (length(outgroup_tips) > biggest_sib_node_size) {
  outgroup_tips <- biggest_sib_node_tips
}

print(paste(length(outgroup_tips), "outgroup tips"))

rooted_tree <- ape::root(phy = tree, node = biggest_sib_node)

# Plot tree rooted with defined outgroup
tree_figname <- paste(OUTDIR, paste(TREENAME, "_rooted.png", sep = ""), sep = "/")
png(file = tree_figname, height = 15, width = 15, units = "in", res = 250)
ggtree(tr = tree, layout="fan", open.angle = 120) + 
  geom_tippoint(aes(color = label %in% outgroup_tips))
dev.off()

# Plot clade assignment of outgroup seqs vs. ingroup seqs
metadata[metadata$strain %in% outgroup_tips, "group"] <- "outgroup"
p1 <- ggplot(data = metadata, aes(x = GISAID_clade)) + 
  ggtitle(paste(TREENAME, "GISAID clade assignments")) + 
  geom_bar(aes(fill = group))

p2 <- ggplot(data = metadata, aes(x = substr(x = pangolin_lineage, start = 1, stop = 3))) + 
  ggtitle(paste(TREENAME, "pangolin lineage assignments (to 1st digit)")) + 
  geom_bar(aes(fill = group))

figname <- paste(OUTDIR, paste(TREENAME, "_outgroup_clades.png", sep = ""), sep = "/")
png(file = figname, height = 5, width = 6.5, res = 300, units = "in")
gridExtra::grid.arrange(p1, p2, nrow = 2)
dev.off()

# Write out outgroup seqs for IQ-TREE
outgroup_str <- paste("'", paste0(outgroup_tips, collapse = ","), "'", sep = "")
filename <- paste(OUTDIR, paste(TREENAME, ".outgroup.txt", sep = ""), sep = "/")
writeLines(
  text = outgroup_str,
  con = file(filename))
