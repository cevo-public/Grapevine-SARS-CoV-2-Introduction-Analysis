# Plot sample numbers vs. estimated imports from each country in each month

require(dplyr)
require(ggplot2)
require(ggrepel)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
PREFIX <- "rep_1_n_sim_1000_n_imports_padded_0"
N_CONTEXT_SEQ_DATA <- paste(WORKDIR, "data/est_imports/samples_per_country_w_import_padding_0.txt", sep = "/")

METADATA_CONTEXT <- paste(WORKDIR, "/data/alignments/", PREFIX, "/", PREFIX, "_context_metadata.txt", sep = "")
METADATA_PRIORITY <- paste(WORKDIR, "/data/alignments/", PREFIX, "/", PREFIX, "_priority_metadata.txt", sep = "")
METADATA_ALL <- paste(WORKDIR, "/data/alignments/", PREFIX, "/", PREFIX, "_tree_metadata.txt", sep = "")
OUTDIR <- paste(WORKDIR, "figures/n_samples", sep = "/")

system(command = paste("mkdir -p", OUTDIR))

# Load data
metadata <- read.delim(file = METADATA_ALL, stringsAsFactors = F, sep = "\t", quote = "")
priority_metadata <- read.delim(file = METADATA_PRIORITY, stringsAsFactors = F, sep = "\t", quote = "")
context_metadata <- read.delim(file = METADATA_CONTEXT, stringsAsFactors = F, sep = "\t", quote = "")
samples_per_country_month <- read.delim(file = N_CONTEXT_SEQ_DATA, stringsAsFactors = F, sep = "\t")

sim_dataset_name <- "Similarity data set"
context_dataset_name <- "Context data set"

context_metadata$analysis <- "context"
priority_metadata$analysis <- "priority"
metadata_w_analysis <- merge(context_metadata, priority_metadata, all = T)

samples_per_country_month_for_both_analyses <- rbind(
  samples_per_country_month %>% mutate(analysis = "context"),
  samples_per_country_month %>% mutate(analysis = "priority"))

metadata_country_year_month <- metadata_w_analysis %>%
  group_by(analysis, country_recoded, year_month) %>% 
  summarize(n_samples = n())
metadata_country_year_month <- merge(
  x = metadata_country_year_month, 
  y = samples_per_country_month_for_both_analyses,
  all = T)

# After merge, year_months without samples included should be changed from NA > 0
metadata_country_year_month$n_samples[is.na(metadata_country_year_month$n_samples)] <- 0

# year_months with samples but no import data also get 0 for estimated imports
metadata_country_year_month$est_n_imports[is.na(metadata_country_year_month$est_n_imports)] <- 0

metadata_country_year_month$year_month <- factor(
  x = metadata_country_year_month$year_month,
  levels = c("2019-12", "2020-1", "2020-2", "2020-3", "2020-4", "2020-5", "2020-6", "2020-7", "2020-8", "2020-9"),
  labels = c("Dec. 2019", "Jan. 2020", "Feb. 2020", "Mar. 2020", "Apr. 2020", "May 2020", "Jun. 2020", "Jul. 2020", "Aug. 2020", "Sept. 2020"))

metadata_country_year_month$analysis <- factor(
  x = metadata_country_year_month$analysis,
  levels = c("priority", "context"),
  labels = c(sim_dataset_name, context_dataset_name))

get_month <- function(year_month) {
  return(strsplit(x = year_month, split = " ")[[1]][1])
}

metadata_country_year_month$month <- 
  unlist(lapply(
    X = as.character(metadata_country_year_month$year_month),
    FUN = get_month))

metadata_country_year_month <- metadata_country_year_month %>% mutate(
  label = case_when(
    n_samples > 5 ~ paste(country_recoded, " (", month, ")", sep = ""),
    est_n_imports > 1 ~ paste(country_recoded, " (", month, ")", sep = ""),
    n_samples > 1 & analysis == context_dataset_name ~ paste(country_recoded, " (", month, ")", sep = ""),
    T ~ ""))  # label only some points so that plot not too busy

p1 <- ggplot(
  data = metadata_country_year_month,
  aes(x = est_n_imports, y = n_samples)) +
  scale_x_continuous(
    trans=scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 10, 20, 40, 80, 160, 320, 640)) +
  scale_y_continuous(
    trans=scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 10, 20, 40, 80, 160, 320, 640)) +
  geom_point() + 
  geom_text_repel(aes(label = label), size = 2.5) +
  facet_wrap(. ~ analysis, scales = "free_x") +
  labs(x = "Estimated # infected individuals arrived", y = "No. sequences included in analysis") + 
  # scale_color_discrete(name = element_blank()) + 
  theme_bw() + 
  theme(
    text = element_text(size = 12), 
    legend.position = "bottom",
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5))

metadata_country <- metadata_country_year_month %>%
  group_by(country_recoded, analysis) %>% 
  summarize(
    est_n_imports = sum(est_n_imports), 
    n_samples = sum(n_samples))

metadata_country <- metadata_country %>% mutate(
  label = case_when(
    n_samples > 5 ~ country_recoded,
    est_n_imports > 5 ~ country_recoded,
    est_n_imports == 0 & analysis == sim_dataset_name & n_samples > 2 ~ country_recoded,
    n_samples > 0 & analysis == context_dataset_name ~ country_recoded,
    T ~ ""))  # label only some points so that plot not too busy

p2 <- ggplot(
  data = metadata_country,
  aes(x = est_n_imports, y = n_samples)) +
  scale_x_continuous(
    trans=scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 10, 20, 40, 80, 160, 320, 640)) +
  scale_y_continuous(
    trans=scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 5, 10, 20, 40, 80, 160, 320, 640)) +
  geom_point() + 
  geom_text_repel(aes(label = label), size = 2.5) +
  facet_wrap(. ~ analysis, scales = "free_x") +
  labs(x = "Estimated # infected individuals arrived", y = "No. sequences included in analysis") + 
  theme_bw() + 
  theme(
    text = element_text(size = 12), 
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5))

metadata_country_year_month <- metadata_country_year_month %>%
  mutate(country_group = case_when(
    country_recoded %in% c("Belgium", "Netherlands") ~ country_recoded,
    T ~ "other"))

# Was thinking about Emma's co-author comments email:
# metadata_country_year_month$country_group <- factor(
#   x = metadata_country_year_month$country_group,
#   levels = c("other", "Belgium", "Netherlands"))
# 
# ggplot(
#   data = metadata_country_year_month,
#   aes(x = year_month, y = n_samples, fill = country_group)) + 
#   facet_wrap(. ~ analysis, scales = "free_y") + 
#   geom_col(position = position_fill())

png(
  file = paste(OUTDIR, paste(PREFIX, "n_samples_per_country_year_month.png", sep = "_"), sep = "/"), 
  width = 6.5, height = 5, res = 300, units = "in")
print(p1)
dev.off()

png(
  file = paste(OUTDIR, paste(PREFIX, "n_samples_per_country.png", sep = "_"), sep = "/"), 
  width = 6.5, height = 3.5, res = 300, units = "in")
print(p2)
dev.off()
