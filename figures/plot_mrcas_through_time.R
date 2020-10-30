# Plot MRCAs through time

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
PREFIX_DATA <- "rep_1_n_sim_1000_n_imports_padded_0"
PREFIX_PLOT <- "rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1"
PREFIX_T <- "rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_T"
PREFIX_F <- "rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_F"

TREE_DATA_WITH_ASR_T <- paste(WORKDIR, "/asr/", PREFIX_T, "_tree_data_with_asr.txt", sep = "")
CLUSTERS_T <- paste(WORKDIR, "/clusters/", PREFIX_T, "_clusters.txt", sep = "")
TREE_DATA_WITH_ASR_F <- paste(WORKDIR, "/asr/", PREFIX_F, "_tree_data_with_asr.txt", sep = "")
CLUSTERS_F <- paste(WORKDIR, "/clusters/", PREFIX_F, "_clusters.txt", sep = "")
METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
OUTDIR <- paste(WORKDIR, "/figures/mrca_dates/", sep = "")

source(paste(WORKDIR, "scripts/utility_functions.R", sep = "/"))
source(paste(WORKDIR, "figures/scripts/plotting_utility_functions.R", sep = "/"))
system(command = paste("mkdir -p", OUTDIR))

# Preliminaries ----------------------------------------------------------------
week_to_date_ranges <- get_week_to_date_range_mapping(
  min_date = "2019-12-01", 
  max_date = "2020-08-31")
get_week_start <- function(date_range) {
  return(strsplit(x = date_range, split = " - ")[[1]][1])
}
week_to_date_ranges$week_start <- unlist(lapply(
  X = week_to_date_ranges$date_range, 
  FUN = get_week_start))

color_scale <- scale_fill_manual(
  values = c("#FDC086", "#BEAED4"), 
  labels = c("FALSE" = "Singletons", "TRUE" = "Transmission chains"),
  name = element_blank())

# Load data --------------------------------------------------------------------
metadata <- read.table(
  file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
cluster_data_T <- read.delim(file = CLUSTERS_T)
cluster_data_F <- read.delim(file = CLUSTERS_F)
tree_data_with_asr_T <- read.delim(file = TREE_DATA_WITH_ASR_T, stringsAsFactors = F)
tree_data_with_asr_F <- read.delim(file = TREE_DATA_WITH_ASR_F, stringsAsFactors = F)

# Get cluster information ------------------------------------------------------
get_cluster_information <- function(cluster_data, metadata) {
  cluster_data_by_tip <- get_cluster_data_by_tip(
    cluster_data = cluster_data, metadata = metadata)
  
  cluster_summary <- cluster_data_by_tip %>%
    group_by(cluster_idx, foreign_tmrca, ch_tmrca, size) %>%
    summarize(first_sample = min(date),
              last_sample = max(date))
  
  cluster_summary$ch_tmrca_week <- lubridate::week(cluster_summary$ch_tmrca)
  cluster_summary$ch_tmrca_week <- factor(
    x = cluster_summary$ch_tmrca_week,
    levels = week_to_date_ranges$week,
    labels = week_to_date_ranges$date_range)
  
  cluster_summary$foreign_tmrca_week <- lubridate::week(cluster_summary$foreign_tmrca)
  cluster_summary$foreign_tmrca_week <- factor(
    x = cluster_summary$foreign_tmrca_week,
    levels = week_to_date_ranges$week,
    labels = week_to_date_ranges$date_range)
  
  return(cluster_summary)
}

cluster_summary_T <- get_cluster_information(cluster_data_T, metadata)
cluster_summary_F <- get_cluster_information(cluster_data_F, metadata)

cluster_summary_T$polyswiss <- "minimum introductions,\nmaximum local transmission"
cluster_summary_F$polyswiss <- "maximum introductions,\nminimum local transmission"

cluster_summary <- rbind(cluster_summary_T, cluster_summary_F)

ch_mrcas_through_time <- ggplot(
  data = cluster_summary,
  aes(x = ch_tmrca_week, fill = size > 1)) + 
  geom_bar(position = "stack") + 
  geom_text(stat='count', aes(label=..count.., group = size > 1), size = 2, position = position_stack(vjust = 0.5)) + 
  theme_bw() + 
  labs(x = element_blank(), y = "Count") + 
  color_scale + 
  scale_x_discrete(labels = week_to_date_ranges$week_start, limits = week_to_date_ranges$date_range) + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "bottom") + 
  facet_wrap(. ~ polyswiss, nrow = 1, scales = "free_y")
# show(ch_mrcas_through_time)

foreign_mrcas_through_time <- ggplot(
  data = cluster_summary,
  aes(x = foreign_tmrca_week, fill = size > 1)) + 
  geom_bar(position = "stack") + 
  geom_text(stat='count', aes(label=..count.., group = size > 1), size = 2, position = position_stack(vjust = 0.5)) + 
  theme_bw() + 
  labs(x = element_blank(), y = "Count") + 
  scale_x_discrete(labels = week_to_date_ranges$week_start, limits = week_to_date_ranges$date_range) + 
  color_scale +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "bottom") + 
  facet_wrap(. ~ polyswiss, nrow = 1, scales = "free_y")
# show(foreign_mrcas_through_time)

library(gtable)
library(grid)
g1 <- ggplotGrob(foreign_mrcas_through_time  + ggtitle("a)") + theme(legend.position = "none", axis.text.x = element_blank()))
g2 <- ggplotGrob(ch_mrcas_through_time + ggtitle("b)") + labs(x = "Beginning of sampling week"))
g <- rbind(g1, g2, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths)

png(
  file = paste(OUTDIR, "/", PREFIX_PLOT, "mrcas_through_time.png", sep = ""),
  width = 6.5, height = 6.5, units = "in", res = 300)
grid.newpage()
grid.draw(g)
dev.off()

longest_lived_chains <- cluster_summary %>% 
  mutate(chain_longevity_days = as.Date(last_sample) - as.Date(first_sample)) %>%
  arrange(desc(chain_longevity_days)) %>%
  filter(as.Date(first_sample) < as.Date("2020-06-01"), as.Date(last_sample) > as.Date("2020-08-01"))

write.table(
  x = longest_lived_chains,
  file = paste(OUTDIR, "/", PREFIX_PLOT, "_longest_chain_info.txt", sep = ""),
  row.names = F, col.names = T, quote = F, sep = "\t")

# oldest_mrca <- min(as.Date(cluster_summary$ch_tmrca))
# oldest_mrca_uncertainty <- cluster_data[as.Date(cluster_data$ch_tmrca) == min(as.Date(cluster_data$ch_tmrca)), "ch_tmrca_CI"]
