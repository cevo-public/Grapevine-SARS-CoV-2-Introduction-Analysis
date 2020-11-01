# This script is to determine how many seqs to take from each country in each month

require(ggplot2)
require(ggrepel)
require(argparse)
require(dplyr)
require(lubridate)

# WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-29_ch_cluster_analysis_better_downsampling"
# IMPORTS_PER_COUNTRY_MONTH <- paste(WORKDIR, "data/est_imports/estimated_imports_per_country_month.txt", sep = "/")
# METADATA <- paste(WORKDIR, "data/est_imports/metadata_all.txt", sep = "/")
# OUTDIR_DATA <- paste(WORKDIR, "data", sep = "/")
# OUTDIR <- paste(WORKDIR, "figures/n_samples", sep = "/")
# PADDING <- 1
# N_CONTEXT_SEQS <- 1000 # try out different values until approx. correct # seqs sampled (not exact b/c of rounding # imports)
# MAX_DATE <- 2020-08

parser <- argparse::ArgumentParser()
parser$add_argument("--importsdata", type="character")
parser$add_argument("--metadata", type="character")
parser$add_argument("--padding", type="integer")
parser$add_argument("--approxncontextseqs", type="integer")
parser$add_argument("--outdirdata", type="character")
parser$add_argument("--outdirfigs", type="character")
parser$add_argument("--maxdate", type="character")

args <- parser$parse_args()

IMPORTS_PER_COUNTRY_MONTH <- args$importsdata
METADATA <- args$metadata
PADDING <- args$padding
N_CONTEXT_SEQS <- args$approxncontextseqs
OUTDIR_DATA <- args$outdirdata
OUTDIR <- args$outdirfigs
MAX_DATE <- parse_date_time(args$maxdate, "%y-%m-%d")

system(command = paste("mkdir -p", OUTDIR))

# Load data
metadata <- read.table(file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
imports <- read.delim(file = IMPORTS_PER_COUNTRY_MONTH, stringsAsFactors = F, sep = "\t")
imports$year_month_date <- parse_date_time(imports$year_month, "%y-%m")

# Divide total # context sequences between countries according to # cases per country month
imports <- imports %>% filter(year_month_date <= MAX_DATE)

imports$est_n_imports_padded <- imports$est_n_imports + PADDING
imports$n_seqs_ideal <- round(
  imports$est_n_imports_padded * N_CONTEXT_SEQS / sum(imports$est_n_imports_padded, na.rm = T))

seqs_by_country_year_month <- metadata %>% 
  group_by(country_recoded, year_month) %>%
  summarize(n_seqs_available = n())

sampling_info <- merge(x = imports, y = seqs_by_country_year_month, all.x = T)
sampling_info$n_seqs_available[is.na(sampling_info$n_seqs_available)] <- 0
sampling_info <- sampling_info %>% filter(!is.na(est_n_imports))  # don't sample from countries with no import data

# When possible, take sequences from appropriate month
sampling_info$n_seqs_actual <- ifelse(
  test = sampling_info$n_seqs_ideal < sampling_info$n_seqs_available,
  yes = sampling_info$n_seqs_ideal,
  no = 0)
sampling_info <- sampling_info %>% 
  mutate(n_seqs_remaining = n_seqs_available - n_seqs_actual,
         n_seqs_accounted_for = n_seqs_actual)

# Otherwise take from next month
is_first <- T
for (country_c in unique(sampling_info$country_recoded)) {
  sampling_info_i <- sampling_info %>% 
    filter(country_recoded == country_c) %>%
    arrange(year_month_date)
  for (j in 1:(nrow(sampling_info_i) - 1)) {
    n_missing <- sampling_info_i[j, "n_seqs_ideal"] - sampling_info_i[j, "n_seqs_accounted_for"]
    k <- j + 1
    n_extra_taken <- min(n_missing, sampling_info_i[k, "n_seqs_remaining"])
    sampling_info_i[k, "n_seqs_actual"] <- sampling_info_i[k, "n_seqs_actual"] + n_extra_taken
    sampling_info_i[j, "n_seqs_accounted_for"] <- sampling_info_i[j, "n_seqs_accounted_for"] + n_extra_taken
    sampling_info_i[k, "n_seqs_remaining"] <- sampling_info_i[k, "n_seqs_remaining"] - n_extra_taken
    n_missing <- n_missing - n_extra_taken
  }
  if (is_first) {
    sampling_info_fill_forward <- sampling_info_i
    is_first <- F
  } else {
    sampling_info_fill_forward <- rbind(sampling_info_fill_forward, sampling_info_i)
  }
}

# Otherwise take from previous month
is_first <- T
for (country_c in unique(sampling_info_fill_forward$country_recoded)) {
  sampling_info_i <- sampling_info_fill_forward %>% 
    filter(country_recoded == country_c) %>%
    arrange(year_month_date)
  for (j in nrow(sampling_info_i):2) {
    n_missing <- sampling_info_i[j, "n_seqs_ideal"] - sampling_info_i[j, "n_seqs_accounted_for"]
    k <- j - 1
    n_extra_taken <- min(n_missing, sampling_info_i[k, "n_seqs_remaining"])
    sampling_info_i[k, "n_seqs_actual"] <- sampling_info_i[k, "n_seqs_actual"] + n_extra_taken
    sampling_info_i[j, "n_seqs_accounted_for"] <- sampling_info_i[j, "n_seqs_accounted_for"] + n_extra_taken
    sampling_info_i[k, "n_seqs_remaining"] <- sampling_info_i[k, "n_seqs_remaining"] - n_extra_taken
    n_missing <- n_missing - n_extra_taken
  }
  if (is_first) {
    sampling_info_fill_backward <- sampling_info_i
    is_first <- F
  } else {
    sampling_info_fill_backward <- rbind(sampling_info_fill_backward, sampling_info_i)
  }
}

sampling_info_summary <- sampling_info_fill_backward %>%
  group_by(country_recoded) %>%
  summarise(actual_total_seqs = sum(n_seqs_actual),
            ideal_total_seqs = sum(n_seqs_ideal))

undersampled_countries <- sampling_info_summary %>% 
  filter(ideal_total_seqs > actual_total_seqs)

if (nrow(undersampled_countries) > 0) {
  print("Some countries are undersampled:")
  print(undersampled_countries)
} 

n_context_seqs <- sum(sampling_info_fill_backward$n_seqs_actual)
print(paste(n_context_seqs, "seqs to be included as context."))

write.table(
  x = sampling_info_fill_backward,
  file = paste(
    OUTDIR_DATA,
    paste("samples_per_country_w_import_padding_", PADDING, ".txt", sep = ""),
    sep = "/"),
  row.names = F, col.names = T, quote = F, sep = "\t")

sampling_info_fill_backward <- sampling_info_fill_backward %>% mutate(
  label = case_when(
    n_seqs_actual > 5 ~ country_recoded,
    est_n_imports > 5 ~ country_recoded,
    est_n_imports == 0 & n_seqs_actual > 2 ~ country_recoded,
    T ~ ""))  # label only some points so that plot not too busy

p1 <- ggplot(
  data = sampling_info_fill_backward,
  aes(x = est_n_imports, y = n_seqs_actual)) +
  scale_x_continuous(
    trans=scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 10, 20, 40, 80, 160, 320, 640)) +
  scale_y_continuous(
    trans=scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 10, 20, 40, 80, 160, 320, 640)) +
  geom_point(aes(color = year_month)) +
  geom_text_repel(aes(label = label), size = 2.5) +
  labs(
    x = "Estimated # infected individuals arrived", 
    y = "No. sequences included in analysis") +
  scale_color_discrete(name = element_blank()) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "bottom",
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5)) + 
  ggtitle(paste(n_context_seqs, "total context sequences"))

png(
  file = paste(
    OUTDIR,
    paste("samples_per_country_month_w_import_padding_", PADDING, ".png", sep = ""),
    sep = "/"),
  width = 5, height = 4, units = "in", res = 300)
show(p1)
dev.off()

sampling_info_fill_backward_country <- sampling_info_fill_backward %>%
  group_by(country_recoded) %>%
  summarize(est_n_imports = sum(est_n_imports),
            n_seqs_actual = sum(n_seqs_actual)) %>%
  mutate(label = case_when(
    n_seqs_actual > 0 ~ country_recoded,
    T ~ ""))

p2 <- ggplot(
  data = sampling_info_fill_backward_country,
  aes(x = est_n_imports, y = n_seqs_actual)) +
  scale_x_continuous(
    trans=scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 10, 20, 40, 80, 160, 320, 640)) +
  scale_y_continuous(
    trans=scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 10, 20, 40, 80, 160, 320, 640)) +
  geom_point() +
  geom_text_repel(aes(label = label), size = 2.5) +
  labs(
    x = "Estimated # infected individuals arrived", 
    y = "No. sequences included in analysis") +
  scale_color_discrete(name = element_blank()) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "bottom",
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5)) + 
  ggtitle(paste(n_context_seqs, "total context sequences"))

png(
  file = paste(
    OUTDIR,
    paste("samples_per_country_w_import_padding_", PADDING, ".png", sep = ""),
    sep = "/"),
  width = 5, height = 4, units = "in", res = 300)
show(p2)
dev.off()




