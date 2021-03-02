# Plot MRCAs through time

WORKDIR <- "dont_commit/2021_01_18_all"
PREFIX_DATA <- "rep_1_n_sim_3000_n_imports_padded_1"
PREFIX_PLOT <- "rep_1_n_sim_3000_n_imports_padded_1_m_3_p_1"
PREFIX_T <- "rep_1_n_sim_3000_n_imports_padded_1_m_3_p_1_swiss_polytomies_T"
PREFIX_F <- "rep_1_n_sim_3000_n_imports_padded_1_m_3_p_1_swiss_polytomies_F"

CLUSTERS_T <- paste(WORKDIR, "/clusters/", PREFIX_T, "_clusters.txt", sep = "")
CLUSTERS_F <- paste(WORKDIR, "/clusters/", PREFIX_F, "_clusters.txt", sep = "")
METADATA <- paste(WORKDIR, "/alignment_metadata/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
OUTDIR <- paste(WORKDIR, "/figures/mrca_dates/", sep = "")

require(ggplot2)
require(lubridate)
source("utility_functions.R")
source("figures/plotting_utility_functions.R")
system(command = paste("mkdir -p", OUTDIR))

# Load data --------------------------------------------------------------------
metadata <- read.table(
  file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
cluster_data_T <- read.delim(file = CLUSTERS_T, stringsAsFactors = F)
cluster_data_F <- read.delim(file = CLUSTERS_F, stringsAsFactors = F)

# Preliminaries ----------------------------------------------------------------

color_scale <- scale_fill_manual(
  values = c("#FDC086", "#BEAED4"), 
  labels = c("FALSE" = "Singletons", "TRUE" = "Transmission chains"),
  name = element_blank())

dates <- seq.Date(from = as.Date("2020-01-01"), to = as.Date("2021-01-18"), by = 7)
weeks_to_date_range_mapping <- get_weeks_since_to_date_range_mapping(
  weeks_since = get_week_since_epidemic_start(dates))

# Get cluster information ------------------------------------------------------
get_cluster_information <- function(cluster_data, metadata) {
  cluster_data_by_tip <- get_cluster_data_by_tip(
    cluster_data = cluster_data, metadata = metadata)
  
  cluster_summary <- cluster_data_by_tip %>%
    group_by(cluster_idx, foreign_tmrca, ch_tmrca, size) %>%
    summarize(first_sample = min(date),
              last_sample = max(date))
  
  cluster_summary$ch_tmrca[cluster_summary$ch_tmrca == "2021"] <- "2021-01-01"
  # cluster_summary$ch_tmrca_week <- lubridate::week(cluster_summary$ch_tmrca)
  cluster_summary$ch_tmrca_week <- get_week_since_epidemic_start(cluster_summary$ch_tmrca)
  cluster_summary$ch_tmrca_week <- factor(
    x = cluster_summary$ch_tmrca_week,
    levels = weeks_to_date_range_mapping$weeks_since,
    labels = weeks_to_date_range_mapping$date_range)
  
  cluster_summary$foreign_tmrca_week <- get_week_since_epidemic_start(cluster_summary$ch_tmrca)
  cluster_summary$foreign_tmrca_week <- factor(
    x = cluster_summary$foreign_tmrca_week,
    levels = weeks_to_date_range_mapping$weeks_since,
    labels = weeks_to_date_range_mapping$date_range)
  
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
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "bottom") + 
  facet_wrap(. ~ polyswiss, nrow = 1, scales = "free_y")
# show(ch_mrcas_through_time)

presentation_fig <- ggplot(
  data = cluster_summary %>% 
    filter(polyswiss == "maximum introductions,\nminimum local transmission"),
  aes(x = ch_tmrca_week, fill = size > 1)) + 
  geom_bar(position = "stack") + 
  # geom_text(stat='count', aes(label=..count.., group = size > 1), size = 2, position = position_stack(vjust = 0.5)) + 
  theme_bw() + 
  labs(x = element_blank(), y = "Count") + 
  color_scale + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "none") #+ 
  # ggtitle("Date and number of Swiss TMRCAs")
show(presentation_fig)
ggsave(
  presentation_fig, 
  path = "/Users/nadeaus/Documents/Presentations/21_01_20_PHB_infectious_disease_cluster",
  filename = "ch_mrcas_from_start.png",
  width = 6, height = 2, units = "in")

foreign_mrcas_through_time <- ggplot(
  data = cluster_summary,
  aes(x = foreign_tmrca_week, fill = size > 1)) + 
  geom_bar(position = "stack") + 
  geom_text(stat='count', aes(label=..count.., group = size > 1), size = 2, position = position_stack(vjust = 0.5)) + 
  theme_bw() + 
  labs(x = element_blank(), y = "Count") + 
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
