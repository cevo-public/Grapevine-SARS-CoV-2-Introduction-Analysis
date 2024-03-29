require(ggplot2)
require(dplyr)
require(lubridate)
require(tidyr)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
PREFIX_DATA <- "rep_1_n_sim_1000_n_imports_padded_0"
PREFIX_BASE <- "rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1"
PREFIX_1 <- paste(PREFIX_BASE, "_swiss_polytomies_F", sep = "")
PREFIX_2 <- paste(PREFIX_BASE, "_swiss_polytomies_T", sep = "")

NOTIFICATION_DELAY <- 0 
TRANSMISSION_TO_TEST_DELAY <- 10 
LEGEND_TEXT_SIZE <- 6

CLUSTERS_1 <- paste(WORKDIR, "/clusters/", PREFIX_1, "_clusters.txt", sep = "")
CLUSTERS_2 <- paste(WORKDIR, "/clusters/", PREFIX_2, "_clusters.txt", sep = "")
METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
PREFIX_OUT <- "figure_3_presentation"
OUTDIR <- paste(WORKDIR, "/figures", sep = "")

source(paste(WORKDIR, "figures/scripts/plotting_utility_functions.R", sep = "/"))  # has list of UZH_SAMPLES to exclude
system(command = paste("mkdir -p", OUTDIR))

# Load data --------------------------------------------------------------------
clusters_1 <- read.delim(file = CLUSTERS_1, stringsAsFactors = F)
clusters_2 <- read.delim(file = CLUSTERS_2, stringsAsFactors = F)
metadata <- read.delim(file = METADATA, stringsAsFactors = F, sep = "\t", quote = "")

cluster_1_data_by_tip <- get_cluster_data_by_tip(
  cluster_data = clusters_1, metadata = metadata) %>%
  mutate(week = week(date), swiss_polytomies = F)

cluster_2_data_by_tip <- get_cluster_data_by_tip(
  cluster_data = clusters_2, metadata = metadata) %>%
  mutate(week = week(date), swiss_polytomies = T)

cluster_data_by_tip <- rbind(
  cluster_1_data_by_tip,
  cluster_2_data_by_tip)

cluster_data_by_tip <- cluster_data_by_tip %>% 
  filter(
    country_recoded == "Switzerland", 
    authors == "Christian Beisel et al", 
    !(gisaid_epi_isl %in% UHZ_SAMPLES))  # take only viollier samples

# Preliminaries ----------------------------------------------------------------
week_to_date_ranges <- get_week_to_date_range_mapping(
  min_date = "2020-01-01", 
  max_date = max(cluster_data_by_tip$date)) 

get_week_start <- function(date_range) {
  return(strsplit(x = date_range, split = " - ")[[1]][1])
}
week_to_date_ranges$week_start <- unlist(lapply(
  X = week_to_date_ranges$date_range, 
  FUN = get_week_start))

cluster_discovery_weeks <- cluster_data_by_tip %>% 
  group_by(cluster_idx, swiss_polytomies) %>%
  summarize(week_first_sampled = week(min(date))) 

# Panel 1: # introductions, % introductions through time  ----------------------

introduction_data <- cluster_data_by_tip %>%
  group_by(week, cluster_idx, swiss_polytomies) %>%
  mutate(ordering = 1:n()) %>%
  top_n(n = 1, wt = ordering) %>%  # keep one row for each introduction sampled in each week
  select(-c("ordering"))

get_is_introduction_newly_discovered <- function(week, cluster_idx, cluster_discovery_weeks, swiss_polytomies) {
  introduction_discovery_week <- cluster_discovery_weeks[
    cluster_discovery_weeks$cluster_idx == cluster_idx & cluster_discovery_weeks$swiss_polytomies == swiss_polytomies, 
    "week_first_sampled"]
  return(week == introduction_discovery_week)
}

introduction_data$is_new_introduction <- unlist(mapply(
  week = week(introduction_data$date),
  cluster_idx = introduction_data$cluster_idx,
  swiss_polytomies = introduction_data$swiss_polytomies,
  MoreArgs = list("cluster_discovery_weeks" = cluster_discovery_weeks),
  FUN = get_is_introduction_newly_discovered))

introduction_info <- introduction_data %>%
  group_by(week, swiss_polytomies) %>%
  summarize(n_new_introductions = sum(is_new_introduction))
introduction_info <- introduction_info %>% 
  pivot_wider(
    names_from = swiss_polytomies, 
    names_prefix = "n_new_introductions_swiss_polytomies_",
    values_from = n_new_introductions)

get_is_sample_from_newly_discovered_introduction <- function(cluster_idx, date, cluster_discovery_weeks, swiss_polytomies) {
  introduction_discovery_week <- cluster_discovery_weeks[
    cluster_discovery_weeks$cluster_idx == cluster_idx & cluster_discovery_weeks$swiss_polytomies == swiss_polytomies, 
    "week_first_sampled"]
  return(week(date) == introduction_discovery_week)
}

cluster_data_by_tip$is_from_newly_discovered_introduction <- unlist(mapply(
  date = cluster_data_by_tip$date,
  cluster_idx = cluster_data_by_tip$cluster_idx,
  swiss_polytomies = cluster_data_by_tip$swiss_polytomies,
  MoreArgs = list("cluster_discovery_weeks" = cluster_discovery_weeks),
  FUN = get_is_sample_from_newly_discovered_introduction))

# Keep only one introduction index case, assign all others to local transmission
cluster_data_by_tip$week <- week(cluster_data_by_tip$date)
cluster_data_by_tip <- cluster_data_by_tip %>%
  group_by(cluster_idx, week, swiss_polytomies) %>%
  arrange(date) %>%
  mutate("sample_idx_in_cluster_week" = 1:n()) %>%
  mutate(
    is_first_sample_in_introduction = case_when(
      is_from_newly_discovered_introduction & sample_idx_in_cluster_week == 1 ~ T,
      T ~ F))

sample_info <- cluster_data_by_tip %>% 
  group_by(week, swiss_polytomies) %>%
  summarize(
    per_local_transmission = sum(!is_first_sample_in_introduction) * 100 / n())
sample_info <- sample_info %>% 
  pivot_wider(
    names_from = swiss_polytomies, 
    names_prefix = "per_local_transmission_swiss_polytomies_",
    values_from = per_local_transmission)

# Plot upper/lower bounds on percent samples each week that are introduction index cases 
# vs. local transmission cases

# sample_focused_plot <- introduction_plot + 
sample_focused_plot <- ggplot() + 
  geom_errorbar(
    data = sample_info,
    aes(x = week, 
        ymin = per_local_transmission_swiss_polytomies_FALSE,
        ymax = per_local_transmission_swiss_polytomies_TRUE)) + 
  scale_x_continuous(breaks = week_to_date_ranges$week, labels = week_to_date_ranges$week_start) +
  scale_y_continuous(
    name = "% of samples from\nlocal transmission",
    limits = c(0, 100)) +
  theme_bw()
show(sample_focused_plot)

# Panel 2: # samples, % preventable samples through time  ----------------------

cluster_contact_tracing_cutoffs <- cluster_data_by_tip %>%
  group_by(cluster_idx, swiss_polytomies) %>%
  arrange(date) %>%
  mutate(sample_idx_in_cluster = 1:n()) %>%
  filter(sample_idx_in_cluster == 1) %>%
  summarize(first_sample_date = date) %>%
  mutate(expected_end_sample_date = as.Date(first_sample_date) + 
           NOTIFICATION_DELAY + TRANSMISSION_TO_TEST_DELAY)

get_is_sample_past_expected_end_sample_date <- function(cluster_idx, date, swiss_polytomies, cluster_contact_tracing_cutoffs) {
  cluster_contact_tracing_cutoff <- cluster_contact_tracing_cutoffs[
    cluster_contact_tracing_cutoffs$cluster_idx == cluster_idx & cluster_contact_tracing_cutoffs$swiss_polytomies == swiss_polytomies, "expected_end_sample_date"]
  return(date > cluster_contact_tracing_cutoff)
}

# Make sample-focused data frame about contact tracing failures: for each week, 
# No. & % of samples that should not have been infected under perfect test, trace, isolate implementation  
cluster_data_by_tip$is_preventable <- unlist(mapply(
  date = cluster_data_by_tip$date,
  swiss_polytomies = cluster_data_by_tip$swiss_polytomies,
  cluster_idx = cluster_data_by_tip$cluster_idx,
  MoreArgs = list("cluster_contact_tracing_cutoffs" = cluster_contact_tracing_cutoffs),
  FUN = get_is_sample_past_expected_end_sample_date))

prevention_info <- cluster_data_by_tip %>% 
  group_by(week, swiss_polytomies) %>%
  summarize(
    n_samples = n(),
    per_preventable_samples = 100 * sum(is_preventable) / n(),
    n_preventable_samples = sum(is_preventable))

prevention_info_wide <- prevention_info %>% select(-n_preventable_samples)
  pivot_wider(
    names_from = swiss_polytomies, 
    names_prefix = "per_is_preventable",
    values_from = per_preventable_samples)

prevention_focused_plot <- ggplot() +
  geom_errorbar(
    data = prevention_info_wide,
    aes(x = week, 
        ymin = per_is_preventableFALSE,
        ymax = per_is_preventableTRUE)) + 
  scale_x_continuous(breaks = week_to_date_ranges$week, labels = week_to_date_ranges$week_start) +
  scale_y_continuous(
    name = "% of samples\npotentially preventable",
    limits = c(0, 100)) +
  theme_bw()
show(prevention_focused_plot)

label_1 <- "sample within 10 days of\nintroduction index case"
label_2 <- "sample after 10 days of\nintroduction index case"
prevention_focused_plot_2 <- ggplot() +
  geom_bar(
    data = cluster_data_by_tip %>% filter(!swiss_polytomies),
    aes(x = week, fill = is_preventable)) + 
  scale_x_continuous(breaks = week_to_date_ranges$week, labels = week_to_date_ranges$week_start, name = element_blank()) +
  scale_y_continuous(
    name = "Number of\nsequenced cases") +
  scale_fill_manual(
    name = element_blank(),
    values = c("TRUE" = "#E31A1C", "FALSE" = "#FB9A99"),
    labels = c("TRUE" = label_2,
               "FALSE" = label_1)) + 
  theme_bw()
  
show(prevention_focused_plot_2)

# Arrange panels ---------------------------------------------------------------

sample_focused_plot_legend <- cowplot::get_legend(sample_focused_plot + guides(shape = "none"))
prevention_plot_legend <- cowplot::get_legend(prevention_focused_plot + guides(shape = "none"))

bottom_plot_theme <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
shape_legend <- cowplot::get_legend(prevention_focused_plot + guides(linetype = "none", fill = "none") + theme(legend.position = "bottom"))

png(
  file = paste(OUTDIR, paste(PREFIX_OUT, ".png", sep = ""), sep = "/"),
  width = 6.5, height = 4, units = "in", res = 300)

egg::ggarrange(
  sample_focused_plot + 
    bottom_plot_theme + 
    guides(shape = "none") +  
    labs(x = element_blank()),
  prevention_focused_plot_2 + 
    bottom_plot_theme +
    theme(legend.spacing.y = unit(0, 'cm')),
  nrow = 2, 
  heights = c(1, 1))

dev.off()

