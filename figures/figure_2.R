require(ggplot2)
require(argparse)
require(dplyr)
require(tidyr)
require(scatterpie)
require(lubridate)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
PREFIX_DATA <- "rep_1_n_sim_1000_n_imports_padded_0"
PREFIX <- "rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_F"
MIN_CLUSTER_SIZE <- NULL  # NULL or integer (2 for transmission chains)
PLOT_SINGLETONS <- T # NULL or T

TREE_DATA_WITH_ASR <- paste(WORKDIR, "/asr/", PREFIX, "_tree_data_with_asr.txt", sep = "")
CLUSTERS <- paste(WORKDIR, "/clusters/", PREFIX, "_clusters.txt", sep = "")
METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
OUTDIR <- paste(WORKDIR, "/figures/figure_2/", sep = "")

source(paste(WORKDIR, "scripts/utility_functions.R", sep = "/"))
source(paste(WORKDIR, "figures/scripts/plotting_utility_functions.R", sep = "/"))
system(command = paste("mkdir -p", OUTDIR))

# Hardcoded things -------------------------------------------------------------
countries_to_plot <- c("Italy", "Germany", "France", "United Kingdom", "United States", "Belgium", "Austria", "Netherlands")
color_scale <- c(
  RColorBrewer::brewer.pal(n = length(countries_to_plot), name = "Paired"), "brown", "black")
names(color_scale) <- c(countries_to_plot, "other", "Switzerland")

ch_mrca_color <- "light grey"
foreign_mrca_color <- "dark grey"

# Preliminaries ----------------------------------------------------------------

week_to_date_ranges <- get_week_to_date_range_mapping(
  min_date = "2020-01-01", 
  max_date = "2020-10-15")

# Load data --------------------------------------------------------------------
metadata <- read.table(
  file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
cluster_data <- read.delim(file = CLUSTERS)
tree_data_with_asr <- read.delim(file = TREE_DATA_WITH_ASR, stringsAsFactors = F)

# Get cluster information ------------------------------------------------------
cluster_data_by_tip <- get_cluster_data_by_tip(
  cluster_data = cluster_data, metadata = metadata)

cluster_summary <- cluster_data_by_tip %>%
  group_by(cluster_idx) %>%
  summarize(first_sample = min(date),
            last_sample = max(date))

# Add tmrca information to cluster summary
cluster_data <- cluster_data %>% mutate(cluster_idx = 1:n())  # same ordering as get_cluster_data_by_tip
cluster_summary <- merge(
  x = cluster_summary, y = cluster_data, 
  by = "cluster_idx", all = T) 

if (!is.null(MIN_CLUSTER_SIZE)) {
  cluster_summary_for_plot <- cluster_summary %>%
    filter(size >= MIN_CLUSTER_SIZE) %>%
    group_by(foreign_mrca) %>%
    arrange(first_sample) %>%
    mutate(cluster_order = 1:n())
  cluster_data_by_tip_for_plot <- merge(
    x = cluster_data_by_tip %>% filter(size >= MIN_CLUSTER_SIZE), 
    y = cluster_summary_for_plot[c("cluster_idx", "cluster_order")], 
    by = c("cluster_idx"), all.x = T) 
} else if (!is.null(PLOT_SINGLETONS)) {
  cluster_summary_for_plot <- cluster_summary %>%
    filter(size == 1) %>%
    group_by(foreign_mrca) %>%
    arrange(first_sample) %>%
    mutate(cluster_order = 1:n())
  cluster_data_by_tip_for_plot <- merge(
    x = cluster_data_by_tip %>% filter(size == 1), 
    y = cluster_summary_for_plot[c("cluster_idx", "cluster_order")], 
    by = c("cluster_idx"), all.x = T) 
} else {
  stop("Must specify either minclustersize or plotsingletons.")
}

cluster_summary_for_plot$foreign_tmrca_CI_min <- as.Date(unlist(lapply(FUN = get_CI_min, X = cluster_summary_for_plot$foreign_tmrca_CI)))
cluster_summary_for_plot$foreign_tmrca_CI_max <- as.Date(unlist(lapply(FUN = get_CI_max, X = cluster_summary_for_plot$foreign_tmrca_CI)))
cluster_summary_for_plot$ch_tmrca_CI_min <- as.Date(unlist(lapply(FUN = get_CI_min, X = cluster_summary_for_plot$ch_tmrca_CI)))
cluster_summary_for_plot$ch_tmrca_CI_max <- as.Date(unlist(lapply(FUN = get_CI_max, X = cluster_summary_for_plot$ch_tmrca_CI)))

# Get MRCA location information for MRCA pie charts ----------------------------
mrca_asr_data <- cluster_summary_for_plot %>%
  group_by(foreign_mrca, foreign_tmrca, foreign_tmrca_CI_min) %>%
  summarise(cluster_order_midpoint = max(cluster_order) / 2 + 0.5)

mrca_asr_data <- merge(
  x = mrca_asr_data, y = tree_data_with_asr,
  by.x = "foreign_mrca", by.y = "node", all.x = T)

# Recode ASR locations (format "First.Second_loc_weight") to match country_recoded names
loc_colnames <- colnames(mrca_asr_data)[grep(
  x = colnames(mrca_asr_data),
  pattern = "_loc_weight")]
loc_colnames <- unlist(lapply(X = loc_colnames, FUN = recode_colnames))

colnames(mrca_asr_data) <- unlist(lapply(
  X =  colnames(mrca_asr_data), 
  FUN = recode_colnames))

other_countries <- loc_colnames[!(loc_colnames %in% countries_to_plot)]

mrca_asr_data$other <- rowSums(mrca_asr_data[other_countries])
mrca_asr_data_2 <- mrca_asr_data %>% 
  select(c(foreign_mrca, foreign_tmrca, foreign_tmrca_CI_min, cluster_order_midpoint, "other", all_of(countries_to_plot)))

# Construct transmission chain plot --------------------------------------------
custom_theme_elements <- theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.border = element_blank(), 
  panel.background = element_blank(), 
  panel.grid = element_blank(), 
  panel.spacing.x = unit(0,"line"),
  strip.text.y = element_blank(),
  legend.position = c(0.07, 0.45),
  legend.background = element_rect(fill = "transparent"),
  panel.spacing = unit(0, "lines"))
axis_labs <- labs(x = element_blank(), y = element_blank())

dates_to_plot <- c("2019-12-01", "2020-01-01", "2020-02-01", "2020-03-01", "2020-04-01", "2020-05-01", "2020-06-01", "2020-07-01", "2020-08-01", "2020-09-01", "2020-10-01")
x_scale <- scale_x_continuous(
  breaks = get_days_from_190101(as.Date(dates_to_plot)),
  labels = format(as.Date(dates_to_plot), "%b."),
  limits = c(get_days_from_190101("2019-11-01"), get_days_from_190101("2020-09-15")),
  expand = c(0, 14))

cluster_data_by_tip_for_plot$foreign_mrca <- factor(
  x = cluster_data_by_tip_for_plot$foreign_mrca,
  levels = unique(unlist(cluster_data_by_tip_for_plot %>% arrange(foreign_tmrca) %>% select(foreign_mrca))))
cluster_summary_for_plot$foreign_mrca <- factor(
  x = cluster_summary_for_plot$foreign_mrca,
  levels = unique(unlist(cluster_summary_for_plot %>% arrange(foreign_tmrca) %>% select(foreign_mrca))))
mrca_asr_data_2$foreign_mrca <- factor(
  x = mrca_asr_data_2$foreign_mrca,
  levels = unique(unlist(mrca_asr_data_2 %>% arrange(foreign_tmrca) %>% select(foreign_mrca))))

mrca_asr_data_2$x_axis_days <- get_days_from_190101(mrca_asr_data_2$foreign_tmrca_CI_min) - 7

cluster_data_by_tip_for_plot <- cluster_data_by_tip_for_plot %>% mutate(
  exposure_country_recoded = case_when(
    country_exposure %in% c(countries_to_plot, "Switzerland") ~ country_exposure,
    T ~ "other"))

CROSS_SIZE <- 0.75
POINT_SIZE <- 0.25
LINE_SIZE <- 0.25
if (!is.null(MIN_CLUSTER_SIZE)) {
  transmission_chain_plot <- ggplot() + 
    geom_segment(
      data = cluster_summary_for_plot,
      aes(x = get_days_from_190101(first_sample), xend = get_days_from_190101(last_sample),
          y = cluster_order, yend = cluster_order),
      linetype = "longdash", size = LINE_SIZE) +
    geom_point(
      data = cluster_summary_for_plot,
      aes(x = get_days_from_190101(foreign_tmrca), y = cluster_order),
      shape = 4, color = foreign_mrca_color, size = CROSS_SIZE) +
    geom_point(
      data = cluster_summary_for_plot,
      aes(x = get_days_from_190101(ch_tmrca), y = cluster_order), 
      shape = 4, color = ch_mrca_color, size = CROSS_SIZE) +
    geom_errorbarh(
      data = cluster_summary_for_plot,
      aes(xmin = get_days_from_190101(ch_tmrca_CI_min), 
          xmax = get_days_from_190101(ch_tmrca_CI_max), y = cluster_order),
      color = ch_mrca_color, size = LINE_SIZE) + 
    geom_errorbarh(
      data = cluster_summary_for_plot,
      aes(xmin = get_days_from_190101(foreign_tmrca_CI_min), 
          xmax = get_days_from_190101(foreign_tmrca_CI_max), y = cluster_order),
      color = foreign_mrca_color, size = LINE_SIZE) +
    geom_point(
      data = cluster_data_by_tip_for_plot,
      aes(x = get_days_from_190101(date), y = cluster_order, 
          color = exposure_country_recoded), size = POINT_SIZE) +
    # Plot colored points over the top so they stand out more
    geom_point(
      data = cluster_data_by_tip_for_plot %>% filter(country_exposure != country_recoded),
      aes(x = get_days_from_190101(date), y = cluster_order, 
          color = exposure_country_recoded), size = POINT_SIZE) +
    geom_scatterpie(
      data = mrca_asr_data_2,
      aes(x = x_axis_days, y = cluster_order_midpoint, r = 4),
      cols = c(countries_to_plot, "other"),
      color = NA) +
    theme_bw() + 
    custom_theme_elements +
    scale_fill_manual(
      values = color_scale, 
      name = element_blank(), 
      aesthetics = c("fill", "color"),
      breaks = c(countries_to_plot[countries_to_plot != "Switzerland"], "other")) + 
    x_scale +
    labs(y = "\n\n", x = element_blank()) + 
    facet_grid(foreign_mrca ~ ., scales = "free_y", space = "free_y")
} else {
  transmission_chain_plot <- ggplot() + 
    geom_point(
      data = cluster_data_by_tip_for_plot,
      aes(x = get_days_from_190101(date), y = 1, 
          color = exposure_country_recoded), size = POINT_SIZE) +
    geom_point(
      data = cluster_summary_for_plot,
      aes(x = get_days_from_190101(foreign_tmrca), y = 1),
      shape = 4, color = foreign_mrca_color, size = POINT_SIZE) +
    geom_errorbarh(
      data = cluster_summary_for_plot,
      aes(xmin = get_days_from_190101(foreign_tmrca_CI_min), 
          xmax = get_days_from_190101(foreign_tmrca_CI_max), y = 1),
      color = foreign_mrca_color, size = LINE_SIZE) +
    geom_scatterpie(
      data = mrca_asr_data_2,
      aes(x = x_axis_days, y = 1, r = 4),
      cols = c(countries_to_plot, "other"),
      color = NA) +
    # Show 1st detected sample
    # geom_point(
    #   data = cluster_data_by_tip_for_plot %>% filter(date == as.Date("2020-02-24")),
    #   aes(x = get_days_from_190101(date), y = cluster_order), color = "red", size = 3) +
    theme_bw() + 
    custom_theme_elements +
    scale_fill_manual(
      values = color_scale, 
      name = element_blank(), 
      aesthetics = c("fill", "color"),
      breaks = c(countries_to_plot[countries_to_plot != "Switzerland"], "other")) + 
    x_scale +
    labs(y = "\n\n", x = element_blank()) + 
    facet_grid(foreign_mrca ~ .)
}

# Construct ASR through time plot ----------------------------------------------

cluster_summary$foreign_tmrca_CI_min <- as.Date(unlist(lapply(FUN = get_CI_min, X = cluster_summary$foreign_tmrca_CI)))
cluster_summary$foreign_tmrca_CI_max <- as.Date(unlist(lapply(FUN = get_CI_max, X = cluster_summary$foreign_tmrca_CI)))

# Get MRCA locations including singletons for bottom plot
mrca_asr_data_all <- cluster_summary %>%
  group_by(foreign_mrca, foreign_tmrca, foreign_tmrca_CI_min, foreign_tmrca_CI_max) %>%
  summarise(n_swiss_descendents = sum(size))

mrca_asr_data_all <- merge(
  x = mrca_asr_data_all, y = tree_data_with_asr,
  by.x = "foreign_mrca", by.y = "node", all.x = T)

# Recode ASR locations (format "First.Second_loc_weight") to match country_recoded names
loc_colnames_2 <- colnames(mrca_asr_data_all)[grep(
  x = colnames(mrca_asr_data_all),
  pattern = "_loc_weight")]
loc_colnames_2 <- unlist(lapply(X = loc_colnames_2, FUN = recode_colnames))

colnames(mrca_asr_data_all) <- unlist(lapply(
  X =  colnames(mrca_asr_data_all), 
  FUN = recode_colnames))

other_countries_2 <- loc_colnames_2[!(loc_colnames_2 %in% countries_to_plot)]

mrca_asr_data_all$other <- rowSums(mrca_asr_data_all[other_countries_2])
mrca_asr_data_all_2 <- mrca_asr_data_all %>% 
  select(c(foreign_mrca, foreign_tmrca, foreign_tmrca_CI_min, foreign_tmrca_CI_min, 
           n_swiss_descendents, "other", all_of(countries_to_plot))) %>%
  mutate(foreign_tmrca_week = lubridate::week(foreign_tmrca))

week_to_days_since_190101_mapping <- data.frame(
  week = c(50, 51, 52, 1:43),
  days_since_190101 = seq(from = 346, to = 662, by = 7))  # This will break once the year wraps around!!

get_days_since_190101_for_bin <- function(week) {
  return(week_to_days_since_190101_mapping[
    week_to_days_since_190101_mapping$week == week, "days_since_190101"])
}

mrca_asr_data_all_2$foreign_tmrca_days_since_190101 <- unlist(lapply(
  X = mrca_asr_data_all_2$foreign_tmrca,
  FUN = get_days_from_190101))

mrca_asr_data_all_2$foreign_tmrca_days_since_190101_bin <- unlist(lapply(
  X = mrca_asr_data_all_2$foreign_tmrca_week,
  FUN = get_days_since_190101_for_bin))

mrca_asr_data_all_2_long <- mrca_asr_data_all_2 %>% tidyr::pivot_longer(
  cols = c("other", all_of(countries_to_plot)),
  names_to = "country",
  values_to = "asr_contribution")

country_order <- mrca_asr_data_all_2_long %>% 
  group_by(country) %>% 
  summarise(total_asr_contribution = sum(asr_contribution, na.rm = T)) %>%
  arrange(total_asr_contribution)

mrca_asr_data_all_2_long$country <- factor(
  x = mrca_asr_data_all_2_long$country,
  levels = country_order$country)

# Get number of nodes in time period and number of nodes with ASR information
n_nodes_data <- mrca_asr_data_all_2 %>% 
  group_by(foreign_tmrca_days_since_190101_bin) %>%
  summarize(n_nodes = n(), n_nodes_with_asr = sum(!is.na(other)))  # nodes without ASR have NA for all locations, including "other"
n_nodes_data$n_nodes_label <- paste(
  n_nodes_data$n_nodes_with_asr, "(", n_nodes_data$n_nodes, ")", sep = "")

asr_contributions_through_time_plot <- ggplot(
  data = mrca_asr_data_all_2_long,
  aes(x = foreign_tmrca_days_since_190101_bin)) + 
  geom_col(aes(y = asr_contribution, fill = country), position = "fill") +
  geom_text(
    data = n_nodes_data,
    aes(x = foreign_tmrca_days_since_190101_bin, 
        y = 1.07, 
        label = n_nodes_label),
    angle = 90, hjust = 0, size = 2) + 
  x_scale + 
  scale_y_continuous(
    limits = c(0, 1.75), 
    expand = c(0, 0), 
    breaks = c(0, 0.5, 1)) + 
  scale_fill_manual(values = color_scale) + 
  theme_classic() + 
  labs(x = element_blank(), y = "Frequency") + 
  theme(legend.position = "none", 
        plot.margin=unit(c(-0.25,0.18,0,0), "cm"))

# Assemble figure --------------------------------------------------------------

if (!is.null(MIN_CLUSTER_SIZE)) {
  SUFFIX <- paste("_transmission_chains_minsize_", MIN_CLUSTER_SIZE, ".png", sep = "")
  png(
    file = paste(OUTDIR, paste(PREFIX, SUFFIX, sep = ""), sep = "/"),
    width = 5.5, height = 8, res = 300, units = "in")
  gridExtra::grid.arrange(
    transmission_chain_plot + ggtitle("a)"),
    asr_contributions_through_time_plot + ggtitle("b)"),
    layout_matrix = cbind(c(1, NA), c(1, 2)),
    nrow = 2, ncol = 2,
    heights = c(6.75, 1),
    widths = c(0.03, 1))
  # heights = c(9, 1),
  # widths = c(0.015, 1))
  dev.off()
} else {
  SUFFIX <- "_singletons.png"
  png(
    file = paste(OUTDIR, paste(PREFIX, SUFFIX, sep = ""), sep = "/"),
    width = 5.5, height = 8, res = 300, units = "in")
  show(transmission_chain_plot)
  dev.off()
}


# Things to note for text ------------------------------------------------------

# Earliest local transmission
cluster_data %>% filter(size > 1) %>% arrange(ch_tmrca) %>% head(1)

# Longest lived cluster
cluster_data_by_tip %>% 
  group_by(cluster_idx) %>%
  summarise(
    timespan = as.Date(max(date)) - as.Date(min(date)),
    min_date = min(date),
    max_date = max(date)) %>%
  arrange(desc(timespan)) %>% 
  head(5)

