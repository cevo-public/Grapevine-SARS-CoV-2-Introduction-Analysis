require(ape)
require(dplyr)
require(treeio)
require(ggplot2)
require(ggtree)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-29_ch_cluster_analysis"
PREFIX_DATA <- "rep_1_n_sim_1000_n_imports_padded_0"
PREFIX <- "rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_F"
METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
CONTEXT_METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_context_metadata.txt", sep = "")
TREE_DATA_WITH_ASR <- paste(WORKDIR, "/asr/", PREFIX, "_tree_data_with_asr.txt", sep = "")
CLUSTER_DATA <- CLUSTER_DATA <- paste(WORKDIR, "/clusters/", PREFIX, "_clusters.txt", sep = "")
TREE <- paste(WORKDIR, "/dated_trees/", PREFIX_DATA, ".timetree.nex", sep = "")
OUTDIR <- paste(WORKDIR, "figures/trainsmission_chain_clades", PREFIX, sep = "/")
FONT_SIZE <- 3

system(command = paste("mkdir -p", OUTDIR))

source(paste(WORKDIR, "scripts/utility_functions.R", sep = "/"))
source(paste(WORKDIR, "figures/scripts/plotting_utility_functions.R", sep = "/"))

# Load data
tree <- treeio::read.beast(file = TREE)
tree_data_with_asr <- read.delim(file = TREE_DATA_WITH_ASR, stringsAsFactors = F)
metadata <- read.table(file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
context_metadata <- read.table(file = CONTEXT_METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
cluster_data <- read.delim(file = CLUSTER_DATA, stringsAsFactors = F)

# Manipulate date format for foreign tmrca
if (all(grepl(x = cluster_data$foreign_tmrca, pattern = "\\."))) {
  print("Date format with '.'")
  cluster_data$foreign_tmrca <- as.Date(cluster_data$foreign_tmrca, format = "%d.%m.%Y")
} else {
  cluster_data$foreign_tmrca <- as.Date(cluster_data$foreign_tmrca)
}

# Make sure all tips present in metadata
tree_ids <- tree_data_with_asr$label[!is.na(tree_data_with_asr$label)]  # tips only
check_all_tips_in_metadata(
  metadata_ids = metadata$strain, tree_ids = tree_ids)

tree_data_with_asr <- merge(
  x = tree_data_with_asr, y = metadata %>% select(c("strain", "division")),
  by.x = "label", by.y = "strain", all.x = T)

# Add old_node number and division (canton) to tree
tree_data_with_asr$old_node <- tree_data_with_asr$node
tree_with_data <- full_join(
  tree, 
  tree_data_with_asr %>% select(c("node", "old_node", "division")), by = "node")

# Get list of nodes in context set (ignore all others as they are only in priority set)
context_data <- metadata %>% 
  filter(old_strain %in% context_metadata$strain | country_recoded == "Switzerland")

# Recode ASR locations (format "First.Second_loc_weight") to match country_recoded names
loc_colnames <- colnames(tree_data_with_asr)[grep(
  x = colnames(tree_data_with_asr),
  pattern = "_loc_weight")]

colnames(tree_data_with_asr) <- unlist(lapply(
  X =  colnames(tree_data_with_asr), 
  FUN = recode_colnames))
loc_colnames <- unlist(lapply(X = loc_colnames, FUN = recode_colnames))

countries_to_plot <- c("Switzerland", "Germany", "Austria", "France", "United Kingdom", "United States", "Serbia", "Belgium")
other_countries <- c(loc_colnames[!(loc_colnames %in% countries_to_plot)], "priority")
color_scale <- c(
  RColorBrewer::brewer.pal(n = length(countries_to_plot), name = "Dark2"), 
  rep("black", length(other_countries) - 1),
  "grey")
names(color_scale) <- c(countries_to_plot, other_countries)

# Subset tree & metadata to cluster
unique_parents <- unique(cluster_data$foreign_mrca)

# Specify a set of parent nodes to plot
# cluster_data_by_tip <- get_cluster_data_by_tip(cluster_data = cluster_data, metadata = metadata)
# clusters_of_interest <- cluster_data_by_tip %>%
#   group_by(cluster_idx) %>%
#   filter(all(c("Aargau", "Bern") %in% division))
# foreign_parents_of_interest <- unique(clusters_of_interest$foreign_mrca)

# clusters_of_interest <- cluster_data %>%
#   filter(foreign_tmrca >= as.Date("2020-02-17"), foreign_tmrca <= as.Date("2020-03-03"), size == 1)
# foreign_parents_of_interest <- unique(clusters_of_interest$foreign_mrca)

clusters_of_interest <- cluster_data %>%
  filter(size > 1)
foreign_parents_of_interest <- unique(clusters_of_interest$foreign_mrca)

print(paste("Making plots for", length(foreign_parents_of_interest), "parent nodes of Swiss introductions"))

for (foreign_parent in foreign_parents_of_interest) {
  subset_tree <- treeio::tree_subset(tree = tree_with_data, node = foreign_parent, levels_back = 0)
  subset_tree_data <- tidytree::as_tibble(subset_tree)
  subset_tree_data$date <- as.Date(subset_tree_data$date)
  
  # Merge in ASR information to tree data for plotting
  subset_tree_data_2 <- merge(
    x = subset_tree_data, 
    y = tree_data_with_asr %>% select(c("node", "country_recoded", loc_colnames)), 
    by.x = "old_node", 
    by.y = "node")

  # Add annotations for plotting to tree data
  n_tips <- sum(!is.na(subset_tree_data_2$label))
  subset_tree_data_2 <- subset_tree_data_2 %>% mutate(
    node_annotation = case_when(
      old_node %in% cluster_data$ch_mrca ~ "Swiss"),
    is_context  = case_when(
      label %in% context_data$strain ~ "Context",
      country_recoded == "Switzerland" ~ "Context",
      node <= n_tips ~ "Priority"),
    tip_label = case_when(
      country_recoded == "Switzerland" ~ division,
      !(label %in% context_data$strain) ~ paste("(p)", country_recoded, sep = ""),
      T ~ country_recoded),
    color_by = case_when(
      !(label %in% context_data$strain) ~ "priority",
      T ~ country_recoded))
  
  # Plot cluster
  pies <- ggtree::nodepie(subset_tree_data_2, cols = loc_colnames)
  internal_pies <- pies[(n_tips + 1):length(pies)]
  internal_pies <- lapply(
    X = internal_pies, FUN = function(g) g + scale_fill_manual(values = color_scale))
  
  # Get reasonable date range for plot
  min_date <- cluster_data[cluster_data$foreign_mrca == foreign_parent, "foreign_tmrca"][1]
  max_date <- max(subset_tree_data_2$date)
  
  p <- ggtree(
    tr = subset_tree, 
    size = 0.2,
    mrsd = max(subset_tree_data_2$date), 
    as.Date = T) %<+% subset_tree_data_2 +
    theme_tree2() +
    geom_inset(insets = internal_pies) +
    geom_tiplab(aes(
      label = paste(old_node, ": ", tip_label, sep = ""),
      color = color_by),
      # color = country_recoded),
      # alpha = is_context),
      size = FONT_SIZE, hjust = -0.01) +
    # scale_alpha_manual(
    #   values = c("Priority" = 0.4, "Context" = 1),
    #   labels = c("Priority" = "(iii) genetically similar", "Context" = "(ii) for source\nlocation inference"),
    #   name = "Sample type",
    #   na.translate = F) +
    scale_color_manual(guide = F, values = color_scale) +
    theme(legend.position = "right") +
    geom_nodelab(aes(
      label = ifelse(
        test = is.na(node_annotation),
        yes = old_node,
        no = paste(old_node, "(", node_annotation, ")", sep = ""))),
      size = FONT_SIZE, hjust = 0.5) + 
    scale_x_date(
      date_labels = "%b-%d", 
      limits = c(
        min_date - floor((max_date - min_date) * 0.3),
        max_date + floor((max_date - min_date) * 0.6)))  # add some days at end so labels hopefully don't get cut off
  
  pdf(
    file = paste(OUTDIR, paste("node_", foreign_parent, ".pdf", sep = ""), sep = "/"), 
    height = (n_tips / 9) + 1, width = 6.5)
  show(p)
  dev.off()
}

