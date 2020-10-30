require(ggplot2)
require(dplyr)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
DATE_THRESHOLD <- "2020-06-01"
PREFIX_DATA <- "rep_1_n_sim_1000_n_imports_padded_0"
OUTDIR <- paste(WORKDIR, "figures/figure_1", sep = "/")

FOUNDING_LINEAGE_DATA <- paste(WORKDIR, "/second_wave_lineages/", PREFIX_DATA, "_lineages_crossing_", DATE_THRESHOLD, ".txt", sep = "")
METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
CLADE_DATA <- paste(WORKDIR, "/data/clades/swiss_alignment_filtered2_masked_oneline_clades.tsv", sep = "")

source(paste(WORKDIR, "figures/scripts/plotting_utility_functions.R", sep = "/"))  # has list of UZH_SAMPLES to exclude
system(command = paste("mkdir -p", OUTDIR))

# Load data --------------------------------------------------------------------
founding_lineage_data <- read.delim(file = FOUNDING_LINEAGE_DATA, sep = "\t")
metadata <- read.table(file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
clades <- read.delim(file = CLADE_DATA, sep = ";")
metadata <- merge(
  x = metadata, y = clades[c("seqName", "clade")],
  by.x = "old_strain", by.y = "seqName", all.x = T)
viollier_only_metadata <- metadata %>% 
  filter(
    country_recoded == "Switzerland", 
    authors == "Christian Beisel et al", 
    !(gisaid_epi_isl %in% UHZ_SAMPLES)) # take only viollier samples

# Hardcoded things -------------------------------------------------------------
clade_color_scale <- RColorBrewer::brewer.pal(
  n = length(unique(viollier_only_metadata$clade)), 
  name = "Accent")
names(clade_color_scale) <- c("20A", "20B", "20C", "19A", "19B")

# Panel 1: clades over time ----------------------------------------------------
week_to_date_ranges <- get_week_to_date_range_mapping(
  min_date = "2020-01-01", max_date = "2020-08-31")
get_week_start <- function(date_range) {
  return(strsplit(x = date_range, split = " - ")[[1]][1])
}
week_to_date_ranges$week_start <- unlist(lapply(
  X = week_to_date_ranges$date_range, 
  FUN = get_week_start))

clade_order <- viollier_only_metadata %>% 
  group_by(clade) %>% 
  summarise(n_samples = n()) %>%
  arrange(desc(n_samples))

viollier_only_metadata$clade <- factor(
  x = viollier_only_metadata$clade, levels = clade_order$clade)

clades_over_time_plot <- ggplot(
  data = viollier_only_metadata,
  aes(x = lubridate::week(date), fill = clade)) + 
  geom_bar(position = "fill", color='black') + 
  scale_fill_manual(values = clade_color_scale, name = "Clade") + 
  theme_classic() + 
  scale_x_continuous(
    breaks = as.numeric(week_to_date_ranges$week),
    labels = week_to_date_ranges$week_start) +
  scale_y_continuous(expand = c(0, 0)) + 
  labs(x = "Beginning of sampling week", y = "Frequency") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90, hjust = 1))
show(clades_over_time_plot)

# Panel 2: 2nd wave founding lineages ------------------------------------------
# Get pie chart of # samples caused by each founding lineage after July 1st
founding_lineage_data <- founding_lineage_data %>% mutate(
  per_second_wave_samples = n_viollier_tips_below / sum(n_viollier_tips_below),
  slice_label = case_when(
    per_second_wave_samples > .03 ~ as.character(n_viollier_tips_below),
    T ~ ""))

founding_clade_data <- founding_lineage_data %>%
  group_by(clade) %>%
  summarize(per_second_wave_samples = sum(per_second_wave_samples))
founding_clade_data$clade <- factor(
  x = founding_clade_data$clade, levels = clade_order$clade)
founding_clade_data <- founding_clade_data %>% 
  arrange(desc(clade)) %>%
  mutate(y_pos_max = cumsum(per_second_wave_samples))
founding_clade_data$y_pos_min <- c(
  0, lag(founding_clade_data$y_pos_max)[2:nrow(founding_clade_data)])

second_wave_founding_lineage_plot <- ggplot(
  data = founding_lineage_data,
  aes(x = "", y = per_second_wave_samples, fill = clade)) + 
  geom_bar(stat = "identity", color='black') + 
  coord_polar("y", start = 0) + 
  theme_void() + 
  geom_text(
    aes(label = slice_label, x = 1.8), 
    position = position_stack(vjust = 0.5), 
    size = 1.8) +
  labs(x = element_blank(), y = element_blank(), fill = "Clade") + 
  geom_rect(
    data = founding_clade_data,
    aes(ymin = y_pos_min, ymax = y_pos_max, fill = clade),
    xmin = 1.53, xmax = 1.6, color = "black") + 
  scale_fill_manual(values = clade_color_scale)
show(second_wave_founding_lineage_plot)

# Arrange panels ---------------------------------------------------------------
clades_over_time_plot_legend <- cowplot::get_legend(
  clades_over_time_plot + 
    guides(fill = guide_legend(nrow = 3)))

png(
  file = paste(OUTDIR, "/", PREFIX_DATA, "_second_wave_start_", DATE_THRESHOLD, "_figure_1.png", sep = ""),
  width = 6.5, height = 2.5, units = "in", res = 300)
gridExtra::grid.arrange(
  clades_over_time_plot + 
    theme(legend.position = "none") + 
    ggtitle("a)"),
  second_wave_founding_lineage_plot + 
    theme(legend.position = "none") + 
    ggtitle("b)"),
  clades_over_time_plot_legend,
  nrow = 2, ncol = 2,
  widths = c(2, 1), heights = c(2, 1.75),
  layout_matrix = cbind(c(1, 1), c(2, 3)))
dev.off()

# Write out information for figure legend --------------------------------------
legend_info_outfile_con <- file(
  open = "a",
  paste(OUTDIR, "/", PREFIX_DATA, "_second_wave_start_", DATE_THRESHOLD, "_figure_1_legend_info.txt", sep = ""))
writeLines(
  con = legend_info_outfile_con,
  text = print(paste(
    sum(founding_lineage_data$n_viollier_tips_below),
    "Viollier tips sampled after", DATE_THRESHOLD)))
writeLines(
  con = legend_info_outfile_con,
  text = print(paste(
    sum(as.numeric(founding_lineage_data$slice_label), na.rm = T),
    "samples generated by the labelled lineages.")))
writeLines(
  con = legend_info_outfile_con,
  text = print(paste(
    nrow(founding_lineage_data %>% filter(clade == "20A")), 
    "20A lineages with Viollier descendents cross", 
    DATE_THRESHOLD, 
    "in the tree.")))
writeLines(
  con = legend_info_outfile_con,
  text = print(paste(
    nrow(founding_lineage_data), 
    "total lineages with Viollier descendents cross", 
    DATE_THRESHOLD, 
    "in the tree.")))
close(legend_info_outfile_con)

