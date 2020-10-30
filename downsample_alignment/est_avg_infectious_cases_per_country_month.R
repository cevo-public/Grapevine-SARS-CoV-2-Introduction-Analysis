require(dplyr)
require(tidyr)
require(ggplot2)
require(argparse)
require(lubridate)

# WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-29_ch_cluster_analysis"
# CASE_DATA_LINK <- "https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"
# OUTDIR_DATA <- paste(WORKDIR, "data/est_imports", sep = "/")
# OUTDIR <- paste(WORKDIR, "figures/est_imports", sep = "/")

parser <- argparse::ArgumentParser()
parser$add_argument("--casedatalink", type="character")
parser$add_argument("--outdirdata", type="character")
parser$add_argument("--outdirfigs", type="character")

args <- parser$parse_args()

CASE_DATA_LINK <- args$casedatalink
OUTDIR_DATA <- args$outdirdata
OUTDIR <- args$outdirfigs

TODAYS_DATE <- format(Sys.time(), "%Y-%m-%d")
CASE_DATA_FN <- paste("ecdc_global_cases", TODAYS_DATE, ".csv", sep = "")
CASE_DATA_FP <- paste(OUTDIR_DATA, CASE_DATA_FN, sep = "/")

system(command = paste("mkdir -p", OUTDIR_DATA))
system(command = paste("mkdir -p", OUTDIR))
download.file(url = CASE_DATA_LINK, destfile = CASE_DATA_FP)

# load case data from ECDC
case_data <- read.delim(file = CASE_DATA_FP, sep = ",")
case_data$date <- as.Date(case_data$dateRep, format = "%d/%m/%Y")

case_data <- case_data %>% 
  mutate(
    country = recode(
      .x = countriesAndTerritories,
      `United_States_of_America` = "United States",
      `United_Kingdom` = "United Kingdom"))

# Get number of infectious individuals on each day, assuming individuals are infectious beginning 
# 10 days before the diagnostic test and ending on the day of the diagnostic test
# Because of data corrections, some days have negative #s infectious individuals
# but should hopefully cancel out to roughly correct monthly averages
case_data_2 <- case_data %>%
  group_by(country) %>%
  arrange(date) %>%
  mutate(n_cumul_cases = cumsum(cases)) %>%
  mutate(n_infectious_cases = lead(n_cumul_cases, n = 10) - n_cumul_cases) %>%
  mutate(per_infectious_individuals = n_infectious_cases / popData2019)

case_data_final <- case_data_2 %>%
  mutate(year_month = paste(year(date), month(date), sep = "-")) %>%
  group_by(country, year_month) %>%
  summarise(
    avg_daily_infectious_individuals = mean(n_infectious_cases, na.rm = T),
    avg_population = mean(popData2019, na.rm = T),
    n_cases = sum(cases)) %>%
  mutate(
    avg_per_pop_infectious = avg_daily_infectious_individuals / avg_population)

# Write out results
write.table(
  x = case_data_final,
  file = paste(OUTDIR_DATA, "infectious_pop_by_country_month.txt", sep = "/"),
  row.names = F, col.names = T, quote = F, sep = "\t")

# Plot percent of population infectious over time for sanity-check & supplemental
countries_to_plot <- c("Italy", "Germany", "France", "Austria", "United States", "China", "United Kingdom", "Switzerland")
case_data_to_plot <- case_data_final %>%
  filter(country %in% countries_to_plot)

case_data_to_plot$country <- factor(
  x = case_data_to_plot$country,
  levels = countries_to_plot)

p <- ggplot(
  data = case_data_to_plot,
  aes(x = year_month, y = avg_per_pop_infectious)) +
  geom_col(aes(fill = country)) + 
  theme_bw() + 
  labs(x = element_blank(), y = "Estimated % of population infectious") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none") + 
  facet_wrap(country ~ ., scales = "free_y")

png(
  file = paste(OUTDIR, "infectious_pop_by_country_month.png", sep = "/"), 
  width = 6.5, height = 4, units = "in", res = 300)
show(p)
dev.off()

p2 <- ggplot(
  data = case_data_2 %>% filter(country %in% countries_to_plot),
  aes(x = date, y = cases)) +
  geom_col(aes(fill = country)) + 
  theme_bw() + 
  labs(x = element_blank(), y = "Daily cases") + 
  scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none") + 
  facet_wrap(country ~ ., scales = "free_y")

png(
  file = paste(OUTDIR, "daily_cases_by_country.png", sep = "/"), 
  width = 6.5, height = 4, units = "in", res = 300)
show(p2)
dev.off()
