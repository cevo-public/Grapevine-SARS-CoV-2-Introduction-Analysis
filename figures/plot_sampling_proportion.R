# The weeks of Mar. 18 - 24, Apr. 8 - 14, and May 20 - 26 are particularly under-sampled.
# Do we have any more unsequenced samples from these weeks?

require(lubridate)
require(tidyr)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
VIOLLIER_METADATA <- "/Users/nadeaus/Documents/2019-ncov-data/ViollierConfidential/Metadata/merged_metadata/merged_metadata_14.txt"
GISAID_METADATA <- "/Volumes/nadeaus/2019-ncov-data/gisaid_data/2020-10-07_ch_cluster_analysis/nextmeta_with_unreleased.tsv"
SEQUENCED_SAMPLE_DIR <- "/Volumes/shared/covid19-pangolin/pangolin/consensus_data/batch/samples"
CASE_DATA_PATH <- paste(WORKDIR, "data", "covid19_cases_switzerland_openzh.csv", sep = "/")
CASE_DATA_LINK <- "https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_cases_switzerland_openzh.csv"
OUTDIR <- paste(WORKDIR, "figures/ch_undersampling", sep = "/")

source("/Volumes/shared/covid19-pangolin/pangolin/consensus_qc_scripts/utility_functions.R")
source(paste(WORKDIR, "/figures/scripts/plotting_utility_functions.R", sep = "/"))
system(command = paste("mkdir -p", OUTDIR))

# Load data
viollier_metadata <- read.delim(file = VIOLLIER_METADATA)
seq_sample_dirnames <- list.dirs(
  path = SEQUENCED_SAMPLE_DIR, recursive = F, full.names = F)
sequenced_eth_ids <- unlist(lapply(
  X = seq_sample_dirnames, FUN = get_eth_id_from_sample_name))
sample_metadata <- read.table(file = GISAID_METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)

download.file(url = CASE_DATA_LINK, destfile = CASE_DATA_PATH)
case_data <- read.delim(file = CASE_DATA_PATH, sep = ",")

# Summarize total # seqs
ch_sequences_by_week <- sample_metadata %>%
  filter(country == "Switzerland") %>%
  mutate(week = week(as.Date(date))) %>%
  group_by(week) %>%
  summarize(n_ch_sequences = n())

# Summarize # seqs from viollier
viollier_samples_by_week <- sample_metadata %>%
  filter(authors == "Christian Beisel et al", !gisaid_epi_isl %in% UHZ_SAMPLES) %>%
  mutate(week = week(as.Date(date))) %>%
  group_by(week) %>%
  summarize(n_viollier_sequences = n())

# Summarize confirmed case data
conf_cases_by_week <- case_data %>% 
  select(CH, Date) %>%
  mutate(n_new_cases = CH - lag(CH, n = 1)) %>%
  mutate(week = week(Date)) %>%
  group_by(week) %>%
  summarise(n_conf_cases = sum(n_new_cases))

# Summarize unsequenced positives from viollier
viollier_unseqed_by_week <- viollier_metadata %>%
  mutate(
    is_sequenced = ETH.ID %in% sequenced_eth_ids,
    week = week(as.Date(Order.date, format = "%d.%m.%Y"))) %>% 
  group_by(week) %>%
  summarize(n_viollier_unseqed = n() - sum(is_sequenced))

# Merge all data by week
sampling_by_week <- merge(x = ch_sequences_by_week, y = viollier_samples_by_week, all = T)
sampling_by_week <- merge(x = sampling_by_week, y = conf_cases_by_week, all = T) 
sampling_by_week <- merge(x = sampling_by_week, y = viollier_unseqed_by_week, all = T) 
sampling_by_week[is.na(sampling_by_week)] <- 0

# Make plots
MIN_DATE <- as.Date("2020-02-26")
MAX_DATE <- as.Date("2020-09-15")

week_to_date_ranges <- get_week_to_date_range_mapping(
  min_date = MIN_DATE, max_date = MAX_DATE,
  date_format = "%b. %d") 

sampling_by_week$week <- factor(
  x = sampling_by_week$week,
  levels = week_to_date_ranges$week,
  labels = week_to_date_ranges$date_range)

sampling_by_week <- sampling_by_week %>% filter(!is.na(week))
sampling_by_week <- sampling_by_week %>%
  mutate(
    n_unseq_other = n_conf_cases - n_ch_sequences - n_viollier_unseqed,
    n_seq_other = n_ch_sequences - n_viollier_sequences)
sampling_by_week <- sampling_by_week %>% pivot_longer(
  cols = c("n_unseq_other", "n_viollier_unseqed", "n_viollier_sequences", "n_seq_other"),
  names_to = "sample_type",
  values_to = "count")

sampling_by_week$sample_type <- factor(
  x = sampling_by_week$sample_type,
  levels = c("n_unseq_other", "n_seq_other", "n_viollier_unseqed", "n_viollier_sequences"))

villlier_seq_color <- RColorBrewer::brewer.pal(name = "Dark2", n = 3)[1]
viollier_unseq_color <- RColorBrewer::brewer.pal(name = "Dark2", n = 3)[2]
other_seq_color <- RColorBrewer::brewer.pal(name = "Dark2", n = 3)[3]
unseq_other_color <- "grey"

p1 <- ggplot(
  data = sampling_by_week, 
  aes(x = week, y = count, fill = sample_type)) + 
  geom_col(position = position_fill()) +
  labs(x = element_blank(), y = "Frequency") + 
  theme_bw() + 
  scale_x_discrete(
    breaks = week_to_date_ranges$date_range,
    limits = week_to_date_ranges$date_range) + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "right") + 
  scale_fill_manual(
    values = c("n_viollier_sequences" = villlier_seq_color, 
               "n_viollier_unseqed" = viollier_unseq_color,
               "n_seq_other" = other_seq_color,
               "n_unseq_other" = unseq_other_color),
    labels = c("n_viollier_sequences" = "Sequenced sample\nfrom Viollier", 
               "n_viollier_unseqed" = "Unsequenced sample\nfrom Viollier",
               "n_seq_other" = "Sequenced by others",
               "n_unseq_other" = "Other confirmed case"),
    name = element_blank())
show(p1)

png(
  file = paste(OUTDIR, "sampling_vs_conf_cases.png", sep = "/"),
  width = 6.5, height = 3, units = "in", res = 300)
show(p1)
dev.off()

# What are the actual % values?
sampling_by_week <- sampling_by_week %>% mutate(
  percent_of_conf_cases = 100 * count / n_conf_cases)

write.table(
  x = sampling_by_week,
  file = paste(OUTDIR, "sampling_by_week.txt", sep = "/"),
  row.names = F, col.names = T, quote = F, sep = "\t")

viollier_as_per_conf_cases_by_week <- unlist(sampling_by_week %>% 
  filter(sample_type == "n_viollier_sequences") %>% 
  select(percent_of_conf_cases))

other_as_per_conf_cases_by_week <- unlist(
  sampling_by_week %>%
  filter(sample_type == "n_seq_other") %>%
    select(percent_of_conf_cases))

range(viollier_as_per_conf_cases_by_week)
mean(viollier_as_per_conf_cases_by_week)
quantile(viollier_as_per_conf_cases_by_week, c(0.25, 0.75))

