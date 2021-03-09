require(ggplot2)
require(argparse)
require(dplyr)
require(tidyr)
require(scatterpie)
require(lubridate)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-29_ch_cluster_analysis"
PREFIX_DATA <- "rep_1_n_sim_1000_n_imports_padded_1"
PREFIX <- "rep_1_n_sim_1000_n_imports_padded_1_m_3_p_1_swiss_polytomies_F"
MIN_CLUSTER_SIZE <- NULL  # NULL or integer (2 for transmission chains)
PLOT_SINGLETONS <- T # NULL or T

TREE_DATA_WITH_ASR <- paste(WORKDIR, "/asr/", PREFIX, "_tree_data_with_asr.txt", sep = "")
CLUSTERS <- paste(WORKDIR, "/clusters/", PREFIX, "_clusters.txt", sep = "")
METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
OUTDIR <- paste(WORKDIR, "/figures/figure_2/", sep = "")

source(paste(WORKDIR, "scripts/utility_functions.R", sep = "/"))
source(paste(WORKDIR, "figures/scripts/plotting_utility_functions.R", sep = "/"))
system(command = paste("mkdir -p", OUTDIR))

metadata <- read.table(
  file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
cluster_data <- read.delim(file = CLUSTERS)

cluster_data_by_tip <- get_cluster_data_by_tip(
  cluster_data = cluster_data, metadata = metadata)

cluster_summary <- cluster_data_by_tip %>%
  group_by(cluster_idx, foreign_mrca, foreign_tmrca, ch_mrca, ch_tmrca, foreign_tmrca_CI, ch_tmrca_CI) %>%
  summarize(first_sample = min(date),
            last_sample = max(date))

cluster_summary$foreign_tmrca_CI_min <- as.Date(unlist(lapply(FUN = get_CI_min, X = cluster_summary$foreign_tmrca_CI)))
cluster_summary$foreign_tmrca_CI_max <- as.Date(unlist(lapply(FUN = get_CI_max, X = cluster_summary$foreign_tmrca_CI)))
cluster_summary$ch_tmrca_CI_min <- as.Date(unlist(lapply(FUN = get_CI_min, X = cluster_summary$ch_tmrca_CI)))
cluster_summary$ch_tmrca_CI_max <- as.Date(unlist(lapply(FUN = get_CI_max, X = cluster_summary$ch_tmrca_CI)))

# unique foreign tMRCA dates?
foreign_tmrca_summary <- cluster_summary %>% 
  group_by(foreign_mrca, foreign_tmrca) %>% 
  summarize(n_clusters = n()) %>%
  mutate(week = lubridate::week(foreign_tmrca))
ggplot(
  data = foreign_tmrca_summary,
  aes(x = week)) + 
  geom_bar()

# CH MRCA dates?
ch_tmrca_summary <- cluster_summary %>% 
  group_by(ch_mrca, ch_tmrca) %>% 
  mutate(week = lubridate::week(ch_tmrca))
ggplot(
  data = ch_tmrca_summary,
  aes(x = week)) + 
  geom_bar()


week(as.Date("2020-03-01"))
week(as.Date("2020-03-10"))

week(as.Date("2020-03-11"))
