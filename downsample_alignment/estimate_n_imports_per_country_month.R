require(dplyr)
require(lubridate)
require(ape)
require(ggplot2)
require(argparse)
require(tidyr)

# WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-29_ch_cluster_analysis"
# PER_INFECTIOUS_POP_DATA <- paste(WORKDIR, "/data/est_imports/infectious_pop_by_country_month.txt", sep = "/")
# TRAVEL_DATA <- paste(WORKDIR, "/data/est_imports/travel_per_country_month.txt", sep = "/")
# PRIORITY_INFO <- paste(WORKDIR, "/data/qc_master_alignment/priorities.txt", sep = "/")
# SWISS_SEQS <- paste(WORKDIR, "/data/qc_master_alignment/swiss_alignment_filtered2_masked_oneline.fasta", sep = "/")
# METADATA_RAW <- "/Volumes/groups/2019-ncov-data/gisaid_data/2020-09-29_ch_cluster_analysis/with_not_yet_released_ch_data_metadata_2020-09-29_07-11.tsv"
# OUTDIR_DATA <- paste(WORKDIR, "data/est_imports", sep = "/")
# OUTDIR <- paste(WORKDIR, "figures/est_imports", sep = "/")

parser <- argparse::ArgumentParser()
parser$add_argument("--infectiouspopdata", type="character")
parser$add_argument("--arrivaldata", type="character")
parser$add_argument("--prioritydata", type="character")
parser$add_argument("--metadata", type="character")
parser$add_argument("--swissseqs", type="character")
parser$add_argument("--outdirdata", type="character")
parser$add_argument("--outdirfigs", type="character")

args <- parser$parse_args()

PER_INFECTIOUS_POP_DATA <- args$infectiouspopdata
TRAVEL_DATA <- args$arrivaldata
PRIORITY_INFO <- args$prioritydata
METADATA_RAW <- args$metadata
SWISS_SEQS <- args$swissseqs
OUTDIR_DATA <- args$outdirdata
OUTDIR <- args$outdirfigs

per_infectious_pop <- read.delim(file = PER_INFECTIOUS_POP_DATA, stringsAsFactors = F)
travel <- read.delim(file = TRAVEL_DATA, stringsAsFactors = F)
priorities <- read.delim(file = PRIORITY_INFO, stringsAsFactors = F, header = F)
metadata_raw <- read.delim(file = METADATA_RAW, stringsAsFactors = F) 
swiss_seqs <- ape::read.FASTA(file = SWISS_SEQS)

countries_to_plot <- c("Italy", "Germany", "France", "United Kingdom", "United States", "Belgium", "Austria", "Netherlands")
color_scale <- c(
  RColorBrewer::brewer.pal(n = length(countries_to_plot), name = "Paired"), "brown")
names(color_scale) <- c(countries_to_plot, "other")

# Filter metadata to context sequence set and double-check no bad dates (should be filtered out during QC anyways)
colnames(priorities) <- c("strain", "priority")
metadata_2 <- metadata_raw %>% 
  filter(strain %in% priorities$strain | strain %in% names(swiss_seqs)) %>%
  mutate(date = as.Date(date))%>%
  filter(!is.na(date))

# Merge in priority information
metadata_3 <- merge(x = metadata_2, y = priorities, all.x = T)
metadata_final <- metadata_3 %>%
  mutate(
    year_month = paste(year(date), month(date), sep = "-"),
    country_recoded = recode(
      .x = country,
      `Hong Kong` = "China",
      USA = "United States"))

write.table(
  x = metadata_final, 
  file = paste(OUTDIR_DATA, "metadata_all.txt", sep = "/"),
  row.names = F, col.names = T, quote = F, sep = "\t")

travel_final <- travel %>% 
  mutate(
    year_month = paste(year, month, sep = "-"),
    country_recoded = recode(
      .x = source_country,
      Czechia = "Czech Republic",
      `Taiwan  Chinese Taipei ` = "Taiwan",
      Irland = "Ireland",
      `New Zealand  Oceania` = "New Zealand",
      `Hong Kong` = "China")) %>%
  group_by(country_recoded, year_month) %>%
  summarise(n_arrivals = sum(n_arrivals))

per_infectious_pop_2 <- per_infectious_pop %>% 
  mutate(
    country_recoded = gsub(x = country, pattern = "_", replacement = " ")) %>%
  mutate(
    country_recoded = recode(
      .x = country_recoded,
      Burma = "Myanmar",
      Czechia = "Czech Republic")) %>%
  select(-c(country))

per_infectious_pop_final <- per_infectious_pop_2 %>%
  complete(expand(per_infectious_pop_2, country_recoded, year_month)) %>%
  mutate(
    avg_per_pop_infectious = replace_na(
      data = avg_per_pop_infectious,
      replace = 0))

# Check data completeness
print("Locations with arrivals into CH but no sequences:")
travel_final %>%
  filter(
    n_arrivals > 0, 
    !(country_recoded %in% metadata_final$country_recoded)) %>%
  group_by(country_recoded) %>%
  summarize(n_arrivals = sum(n_arrivals))

print("Locations with arrivals into CH but no information of source country infection load:")
travel_final %>%
  filter(
    n_arrivals > 0,
    !(country_recoded %in% per_infectious_pop_final$country_recoded)) %>%
  group_by(country_recoded) %>%
  summarize(n_arrivals = sum(n_arrivals))

# Get estimated # imports from each country in each month
import_data <- merge(
  x = travel_final, 
  y = per_infectious_pop_final, 
  all.x = T)
import_data_final <- import_data %>%
  mutate(est_n_imports = n_arrivals * avg_per_pop_infectious)

write.table(
  x = import_data_final,
  file = paste(OUTDIR_DATA, "estimated_imports_per_country_month.txt", sep = "/"),
  sep = "\t", row.names = F, col.names = T, quote = F)

# Plot estimated import cases over time for sanity-check & supplemental
early_timeframe <- "Dec. 2019 -\n Feb. 2020"
late_timeframe <- "Mar. 2020 - Sept. 2020"
import_data_to_plot <- import_data_final %>%
  mutate(
    source_location = case_when(
      country_recoded %in% countries_to_plot ~ country_recoded,
      T ~ "other"),
    timeframe = case_when(
      year_month %in% c("2019-12", "2020-1", "2020-2") ~ early_timeframe,
      T ~ late_timeframe))

print(paste("estimated", sum(import_data_to_plot$est_n_imports, na.rm = T), "total imports."))

import_data_to_plot$source_location <- factor(
  x = import_data_to_plot$source_location,
  levels = c(countries_to_plot, "other"))

import_data_to_plot$timeframe <- factor(
  x = import_data_to_plot$timeframe,
  levels = c(early_timeframe, late_timeframe))

p <- ggplot(
  data = import_data_to_plot,
  aes(x = year_month, y = est_n_imports)) +
  geom_col(aes(fill = source_location)) + 
  theme_bw() + 
  labs(x = element_blank(), y = "Estimated No. individuals") + 
  facet_wrap(. ~ timeframe, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  scale_fill_manual(values = color_scale)

# From: https://stackoverflow.com/questions/52341385/how-to-automatically-adjust-the-width-of-each-facet-for-facet-wrap/52422707
gp <- ggplotGrob(p)

# get gtable columns corresponding to the facets (5 & 9, in this case)
facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]

# get the number of unique x-axis values per facet (1 & 3, in this case)
x.var <- sapply(ggplot_build(p)$layout$panel_scales_x,
                function(l) length(l$range$range))

# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var

png(
  file = paste(OUTDIR, "estimated_imports_per_country_month.png", sep = "/"), 
  width = 5, height = 3, units = "in", res = 300)
grid::grid.draw(gp)
dev.off()
