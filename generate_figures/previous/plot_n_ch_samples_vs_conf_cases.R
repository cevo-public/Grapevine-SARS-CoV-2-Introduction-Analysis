# Plot genetic samples vs. confirmed case counts 

require(ggplot2)
require(lubridate)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
PREFIX_DATA <- "rep_1_n_sim_1000_n_imports_padded_0"

METADATA <- paste(WORKDIR, "/data/alignments/", PREFIX_DATA, "/", PREFIX_DATA, "_tree_metadata.txt", sep = "")
CASE_DATA_LINK <- "https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_cases_switzerland_openzh.csv"
CASE_DATA_PATH <- paste(WORKDIR, "data/covid19_cases_switzerland_dprobst.csv", sep = "/")
OUTDIR <- paste(WORKDIR, "figures/n_samples", sep = "/")

system(command = paste("mkdir -p", OUTDIR))
source(paste(WORKDIR, "/figures/scripts/plotting_utility_functions.R", sep = "/"))

# Load data
download.file(url = CASE_DATA_LINK, destfile = CASE_DATA_PATH)
case_data <- read.delim(file = CASE_DATA_PATH, sep = ",")
metadata <- read.table(file = METADATA, sep = "\t", quote = "", fill = T, header = T, stringsAsFactors = F)
metadata$date <- as.Date(metadata$date)

case_data_long <- tidyr::pivot_longer(
  data = case_data, 
  cols = colnames(case_data)[colnames(case_data) != "Date"],
  names_to = "Kanton",
  values_to = "n_conf_cases") %>% 
  filter(!is.na(n_conf_cases))
case_data_long$Date <- as.Date(case_data_long$Date)

Aug_31_confirmed_cases <- case_data_long %>%
  group_by(Kanton) %>%
  filter(Date <= as.Date("2020-08-31")) %>%
  top_n(n = 1, wt = Date)

most_recent_confirmed_cases <- case_data_long %>%
  group_by(Kanton) %>%
  top_n(n = 1, wt = Date)

# What % of confirmed cases is the Swiss sequencing set?
n_swiss_seqs <- nrow(metadata %>% filter(country_recoded == "Switzerland"))
n_conf_cases <- unlist(most_recent_confirmed_cases[most_recent_confirmed_cases$Kanton == "CH", "n_conf_cases"])

per_conf_sequenced <- n_swiss_seqs / n_conf_cases
print(paste(per_conf_sequenced * 100, "% of all confirmed cases are in our sequence set."))

n_conf_cases_Aug_31 <- unlist(Aug_31_confirmed_cases[Aug_31_confirmed_cases$Kanton == "CH", "n_conf_cases"]) 
per_conf_sequenced_Aug_31 <- n_swiss_seqs / n_conf_cases_Aug_31
print(paste(per_conf_sequenced_Aug_31 * 100, "% of confirmed cases until Aug. 31 are in our sequence set."))

IS_VIOLLIER_ONLY <- T

if (IS_VIOLLIER_ONLY) {
  metadata <- metadata %>% 
    filter(
      country_recoded == "Switzerland", 
      authors == "Christian Beisel et al", 
      !(gisaid_epi_isl %in% UHZ_SAMPLES))  # take only viollier samples
}

# Are we sampling proportional to confirmed cases by canton?
sampling_by_location_data <- metadata %>% 
  filter(country_recoded == "Switzerland", division != "Switzerland") %>%  # one sample has non-specific division
  group_by(division) %>%
  summarise(n_samples = n())

sampling_by_location_data$Kanton <- unlist(lapply(
  X = sampling_by_location_data$division,
  FUN = get_kanton_code_from_fullnames))

sampling_by_location_data_2 <- merge(
  x = sampling_by_location_data, 
  y = Aug_31_confirmed_cases[c("Kanton", "n_conf_cases")])

p <- ggplot(
  data = sampling_by_location_data_2,
  aes(x = n_conf_cases, y = n_samples)) + 
  geom_text(aes(label = Kanton), size = 2) + 
  theme_bw() + 
  labs(x = "No. confirmed cases as of Aug. 31", y = "No. genome samples")

png(
  file = paste(OUTDIR, "/viollier_only_", IS_VIOLLIER_ONLY, "_n_ch_samples_vs_conf_cases_aug_31.png", sep = ""),
  width = 4, height = 3, units = "in", res = 300)
show(p)
dev.off()
