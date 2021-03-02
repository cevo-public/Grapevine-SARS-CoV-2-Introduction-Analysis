require(ggplot2)
require(dplyr)
require(tidyr)

# MAX_CLUSTERS <- "/Users/nadeaus/Downloads//tmp/clusters/rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_F_clusters.txt"
# MIN_CLUSTERS <- "/Users/nadeaus/Downloads//tmp/clusters/rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_T_clusters.txt"
# METADATA <- "/Users/nadeaus/Downloads//tmp/alignments/rep_1_n_sim_1000_n_imports_padded_0/rep_1_n_sim_1000_n_imports_padded_0_tree_metadata.txt"
# UTILITY_FUNCTIONS <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/automated//utility_functions.R"
# TRANSMISSION_TO_TEST_DELAY <- 10

parser <- argparse::ArgumentParser()
parser$add_argument("--maxclusters", type="character")
parser$add_argument("--minclusters", type="character")
parser$add_argument("--metadata", type="character")
parser$add_argument("--outdir", type="character")
parser$add_argument("--utilityfunctions", type="character")
parser$add_argument("--transmissiontestdelay", type="double")
parser$add_argument("--prefix", type="character")

args <- parser$parse_args()

MAX_CLUSTERS <- args$maxclusters
MIN_CLUSTERS <- args$minclusters
METADATA <- args$metadata
OUTDIR <- args$outdir
UTILITY_FUNCTIONS <- args$utilityfunctions
TRANSMISSION_TO_TEST_DELAY <- args$transmissiontestdelay
PREFIX <- args$prefix

source(UTILITY_FUNCTIONS)
system(command = paste("mkdir -p", OUTDIR))

# Load data --------------------------------------------------------------------
max_clusters <- read.delim(file = MAX_CLUSTERS, stringsAsFactors = F)
min_clusters <- read.delim(file = MIN_CLUSTERS, stringsAsFactors = F)
metadata <- read.delim(file = METADATA, stringsAsFactors = F, sep = "\t", quote = "")

max_cluster_data_by_tip <- get_cluster_data_by_tip(
  cluster_data = max_clusters, metadata = metadata) %>%
  mutate(week = get_week_since_epidemic_start(date), swiss_polytomies = F)

min_cluster_data_by_tip <- get_cluster_data_by_tip(
  cluster_data = min_clusters, metadata = metadata) %>%
  mutate(week = get_week_since_epidemic_start(date), swiss_polytomies = T)

cluster_data_by_tip <- rbind(
  max_cluster_data_by_tip,
  min_cluster_data_by_tip)

cluster_data_by_tip <- cluster_data_by_tip %>% 
  filter(originating_lab == "Viollier AG")

# Preliminaries ----------------------------------------------------------------
cluster_discovery_weeks <- cluster_data_by_tip %>% 
  group_by(cluster_idx, swiss_polytomies) %>%
  summarize(week_first_sampled = get_week_since_epidemic_start(min(date))) 

weeks_to_date_range_mapping <- get_weeks_since_to_date_range_mapping(
  weeks_since = get_week_since_epidemic_start(cluster_data_by_tip$date))

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
  week = get_week_since_epidemic_start(introduction_data$date),
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

introduction_plot <- ggplot() + 
  geom_errorbar(
    data = introduction_info,
    aes(x = week, ymin = n_new_introductions_swiss_polytomies_TRUE,
        ymax = n_new_introductions_swiss_polytomies_FALSE,
        color = "No. newly sampled introductions")) + 
  scale_color_manual(
    name = element_blank(),
    values = c("No. newly sampled introductions" = "#FDC086"))
# show(introduction_plot)

get_is_sample_from_newly_discovered_introduction <- function(cluster_idx, date, cluster_discovery_weeks, swiss_polytomies) {
  introduction_discovery_week <- cluster_discovery_weeks[
    cluster_discovery_weeks$cluster_idx == cluster_idx & cluster_discovery_weeks$swiss_polytomies == swiss_polytomies, 
    "week_first_sampled"]
  return(get_week_since_epidemic_start(date) == introduction_discovery_week)
}

cluster_data_by_tip$is_from_newly_discovered_introduction <- unlist(mapply(
  date = cluster_data_by_tip$date,
  cluster_idx = cluster_data_by_tip$cluster_idx,
  swiss_polytomies = cluster_data_by_tip$swiss_polytomies,
  MoreArgs = list("cluster_discovery_weeks" = cluster_discovery_weeks),
  FUN = get_is_sample_from_newly_discovered_introduction))

# Keep only one introduction index case, assign all others to local transmission
cluster_data_by_tip$week <- get_week_since_epidemic_start(cluster_data_by_tip$date)
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
legend_title <- "Assumption:"
line_label <- "% samples from local transmission:"
label_1a <- "minimum introductions"
label_2a <- "maximum introductions"
max_introductions <- max(introduction_info$n_new_introductions_swiss_polytomies_FALSE)
SAMPLE_PER_MULT <- max_introductions / 100
sample_focused_plot <- introduction_plot + 
  geom_point(
    data = sample_info,
    aes(x = week, y = SAMPLE_PER_MULT * per_local_transmission_swiss_polytomies_FALSE, 
        shape = label_2a)) +
  geom_line(
    data = sample_info,
    aes(x = week, y = SAMPLE_PER_MULT * per_local_transmission_swiss_polytomies_FALSE, group = 1, 
        linetype = line_label)) + 
  geom_point(
    data = sample_info,
    aes(x = week, y = SAMPLE_PER_MULT * per_local_transmission_swiss_polytomies_TRUE, 
        shape = label_1a)) +
  geom_line(
    data = sample_info,
    aes(x = week, y = SAMPLE_PER_MULT * per_local_transmission_swiss_polytomies_TRUE, group = 1, 
        linetype = line_label)) +
  scale_shape_manual(name = legend_title, values = c(4, 6)) + 
  scale_linetype_manual(name = element_blank(), values = c(1)) + 
  scale_x_continuous(
    breaks = weeks_to_date_range_mapping$weeks_since, 
    labels = weeks_to_date_range_mapping$date_range) + 
  scale_y_continuous(
    name = "No. introductions",
    limits = c(0, max_introductions + 1),
    sec.axis = sec_axis(trans = ~./SAMPLE_PER_MULT, name = "% of samples")) +
  theme_bw() + 
  guides(shape = "none",
         linetype = guide_legend(order=1))

# Panel 2: # samples, % preventable samples through time  ----------------------

cluster_contact_tracing_cutoffs <- cluster_data_by_tip %>%
  group_by(cluster_idx, swiss_polytomies) %>%
  arrange(date) %>%
  mutate(sample_idx_in_cluster = 1:n()) %>%
  filter(sample_idx_in_cluster == 1) %>%
  summarize(first_sample_date = date) %>%
  mutate(expected_end_sample_date = as.Date(first_sample_date) + 
           TRANSMISSION_TO_TEST_DELAY)

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
    n_preventable_samples = sum(is_preventable),
    n_non_preventable_samples = sum(!is_preventable))

max_samples <- max(prevention_info$n_samples)
PREVENTABLE_PER_MULT <- max_samples / 100
legend_title_1 <- "Assuming minimum introductions:"
legend_title_2 <- expression("% samples" <= "10 days of index sample:")
line_label_2 <- "% un-preventable samples"
label_1 <- "un-preventable sample"
label_2 <- "possibly preventable sample"

cluster_data_by_tip$is_preventable <- factor(
  x = cluster_data_by_tip$is_preventable,
  levels = c("TRUE", "FALSE"))

prevention_focused_plot <- ggplot() +
  geom_bar(
    data = cluster_data_by_tip %>% filter(swiss_polytomies),
    aes(x = week, fill = is_preventable)) +
  geom_point(
    data = prevention_info %>% filter(swiss_polytomies),
    aes(x = week, y = PREVENTABLE_PER_MULT * n_non_preventable_samples * 100 / n_samples,
        shape = label_1a)) +
  geom_line(
    data = prevention_info %>% filter(swiss_polytomies),
    aes(x = week, y = PREVENTABLE_PER_MULT * n_non_preventable_samples * 100 / n_samples,
        linetype = line_label_2, group = 1)) +
  geom_point(
    data = prevention_info %>% filter(!swiss_polytomies),
    aes(x = week, y = PREVENTABLE_PER_MULT * n_non_preventable_samples * 100 / n_samples,
        shape = label_2a)) +
  geom_line(
    data = prevention_info %>% filter(!swiss_polytomies),
    aes(x = week, y = PREVENTABLE_PER_MULT * n_non_preventable_samples * 100 / n_samples,
        linetype = line_label_2, group = 1)) +
  scale_linetype_manual(name = element_blank(), values = c(1)) +
  scale_shape_manual(name = legend_title, values = c(4, 6)) +
  scale_fill_manual(
    name = legend_title_1,
    values = c("FALSE" = "#E31A1C", "TRUE" = "#FB9A99"),
    labels = c("TRUE" = label_2,
               "FALSE" = label_1)) +
  guides(fill = guide_legend(order=3),
         shape = guide_legend(order=1),
         linetype = guide_legend(order=2)) +
  scale_x_continuous(
    breaks = weeks_to_date_range_mapping$weeks_since, 
    labels = weeks_to_date_range_mapping$date_range) +
  scale_y_continuous(
    name = "No. samples",
    limits = c(0, max_samples + 1),
    sec.axis = sec_axis(trans = ~./PREVENTABLE_PER_MULT, name = "% of samples")) +
  labs(x = "Beginning of sampling week", y = element_blank()) +
  theme_bw()

# Arrange panels ---------------------------------------------------------------

sample_focused_plot_legend <- cowplot::get_legend(sample_focused_plot + guides(shape = "none"))
prevention_plot_legend <- cowplot::get_legend(prevention_focused_plot + guides(shape = "none"))

bottom_plot_theme <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
shape_legend <- cowplot::get_legend(prevention_focused_plot + guides(linetype = "none", fill = "none") + theme(legend.position = "bottom"))

png(
  file = paste(OUTDIR, paste(PREFIX, "_viollier_only_unpreventable_figure_3.png", sep = ""), sep = "/"),
  width = 6.5, height = 3.5, units = "in", res = 300)
egg::ggarrange(
  sample_focused_plot + 
    theme(axis.text.x = element_blank(), legend.spacing.y = unit(0, 'cm')) + 
    guides(shape = "none") +  
    labs(x = element_blank()) + 
    ggtitle("a)"),
  prevention_focused_plot + 
    bottom_plot_theme +
    theme(legend.spacing.y = unit(0, 'cm')) + 
    ggtitle("b)"),
  nrow = 2, 
  heights = c(1, 1))
dev.off()

