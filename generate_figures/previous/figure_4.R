require(tidyverse)
require(lubridate)
require(ggfortify)

# ANALYSIS <- "Samp2"
# DATASET <- "viollier_polyF"
# WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
# WORKDIR_BDSKY <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-22_ch_cluster_analysis/bdsky"
# OUTDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-22_ch_cluster_analysis"
# LINE_LIST_RE_EST_LINK <- "https://raw.githubusercontent.com/covid-19-Re/dailyRe-Data/master/CHE-estimates.csv"
# CASE_DATA <- paste(WORKDIR, "data/covid19_cases_switzerland_dprobst.csv", sep = "/")
# BDSKY_DATE_FILE <- paste(WORKDIR_BDSKY, "/sequences/", DATASET, ".dates.txt", sep = "")
# BDSKY_LOG_FILES <- list.files(
#   path = paste(WORKDIR_BDSKY, "results/logfiles_Re_metaprior", sep = "/"),
#   pattern = paste(ANALYSIS, sep = "."),
#   full.names = T)
# BDSKY_SEQS <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-22_ch_cluster_analysis/data/cluster_alignments/rep_1_n_sim_1000_n_imports_padded_0_m_3_p_1_swiss_polytomies_F_viollier_only_swiss_clusters.fasta"
# OUTDIR <- paste(OUTDIR, "figures/figure_4", sep = "/")
# COMP_TITLE <- "without additional Mar. seqs"

ANALYSIS <- "Samp2"
DATASET <- "swiss2rep2_polyF"
COMP_TITLE <- "with additional Mar. seqs, inferred skyride variance"
WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
WORKDIR_BDSKY <- "/Users/nadeaus/Documents/2019-ncov-data/our_papers/swiss_sars_cov_2/analyses/bdsky"
OUTDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-10-07_ch_cluster_analysis"
LINE_LIST_RE_EST_LINK <- "https://raw.githubusercontent.com/covid-19-Re/dailyRe-Data/master/CHE-estimates.csv"
CASE_DATA <- paste(WORKDIR, "data/covid19_cases_switzerland_dprobst.csv", sep = "/")
BDSKY_DATE_FILE <- paste(WORKDIR_BDSKY, "/sequences/", DATASET, ".dates.txt", sep = "")
BDSKY_LOG_FILES <- list.files(
  path = paste(WORKDIR_BDSKY, "results", sep = "/"),
  pattern = paste(ANALYSIS, DATASET, sep = "."),
  full.names = T)
BDSKY_SEQS <- paste(WORKDIR_BDSKY, "/sequences/", DATASET, ".fasta", sep = "")
OUTDIR <- paste(OUTDIR, "figures/figure_4", sep = "/")
MIN_DATE <- as.Date("2020-05-01")
MAX_DATE <- as.Date("2020-08-31")

SAMP_PROP_PRIOR_BETA_SHAPE_1 <- 1
SAMP_PROP_PRIOR_BETA_SHAPE_2 <- 99
SAMP_PROP_Y_LIM <- 0.2

RE_PRIOR_LOGNORM_MEAN <- 0.8
RE_PRIOR_LOGNORM_SD <- 0.5
RE_Y_LIM <- 3

# MIN_DATE <- NA
# MAX_DATE <- NA

source(paste(WORKDIR, "figures/scripts/plotting_utility_functions.R", sep = "/"))
system(command = paste("mkdir -p", OUTDIR))

conf_case_color <- "#E7298A"
bdsky_case_color <- "#E6AB02"
line_list_re_color <- "#E7298A"
bdsky_re_color <- "#E6AB02"
underreporting_color <- "#D95F02"
  
# Load data --------------------------------------------------------------------
download.file(
  url = LINE_LIST_RE_EST_LINK,
  destfile = paste(WORKDIR, "data/line_list_re_estimates.csv", sep = "/"))
line_list_re_data <- read.csv(
  file = paste(WORKDIR, "data/line_list_re_estimates.csv", sep = "/"))
line_list_re_data <- line_list_re_data %>%
  filter(region == "CHE", data_type == "Confirmed cases", estimate_type == "Cori_slidingWindow") %>%
  select(c(date, median_R_mean, median_R_highHPD, median_R_lowHPD)) %>%
  mutate(date = as.Date(date)) %>%
  rename(
    line_list_median = median_R_mean,
    line_list_high = median_R_highHPD,
    line_list_low = median_R_lowHPD)

case_data <- read.delim(file = CASE_DATA, sep = ",")
case_data <- case_data %>%
  select("Date", "CH") %>%
  mutate("new_cases" = CH - lag(CH, n = 1))

bdsky_data <- loadLogFiles(
  logFiles = BDSKY_LOG_FILES,
  datesFile = BDSKY_DATE_FILE,
  sampling = "Skyline",
  sequences = "Viollier",
  burninFrac = 0.1)
bdsky_re_data <- bdsky_data$R0Data %>%
  group_by(Date) %>%
  summarize(
    bdsky_median = median(R0),
    bdsky_low = quantile(R0, 0.025),
    bdsky_high = quantile(R0, 0.975))
bdsky_sampling_prop_data <- bdsky_data$samplingPropData %>%
  group_by(Date) %>%
  summarize(
    bdsky_median = median(sampProp),
    bdsky_low = quantile(sampProp, 0.025),
    bdsky_high = quantile(sampProp, 0.975)) %>%
  mutate(interval_start_date = Date - 6)

bdsky_seqs <- grep(
  x = readLines(BDSKY_SEQS),
  pattern = "^>",
  value = T)
bdsky_sample_data <- data.frame(header = bdsky_seqs) %>%
  tidyr::separate(
    col = header,
    into = c("strain", "accession_id", "date", "cluster"),
    sep = "\\|") %>%
  mutate(date = as.Date(date))

get_samples_per_interval <- function(interval_start, interval_end) {
  n_samples <- nrow(
    bdsky_sample_data %>%
    filter(date <= as.Date(interval_end), date > as.Date(interval_start)))
  return(n_samples)
}

get_confirmed_cases_per_interval <- function(interval_start, interval_end) {
  interval_case_data <- case_data %>%
    mutate(Date = as.Date(Date)) %>%
    filter(Date <= as.Date(interval_end), Date > as.Date(interval_start))
  return(sum(interval_case_data$new_cases))
}

bdsky_sampling_prop_data$n_seqs <- mapply(
  FUN = get_samples_per_interval,
  interval_start = as.Date(bdsky_sampling_prop_data$interval_start_date),
  interval_end = as.Date(bdsky_sampling_prop_data$Date))

bdsky_sampling_prop_data$n_cases <- mapply(
  FUN = get_confirmed_cases_per_interval,
  interval_start = as.Date(bdsky_sampling_prop_data$interval_start_date),
  interval_end = as.Date(bdsky_sampling_prop_data$Date))

bdsky_sampling_prop_data <- bdsky_sampling_prop_data %>%
  mutate(
    n_cases_median = n_seqs / bdsky_median,
    n_cases_high = n_seqs / bdsky_low,
    n_cases_low = n_seqs / bdsky_high)

# 1st facet: Re through time overlayed -----------------------------------------

if (!is.na(MIN_DATE)) {
  ANALYSIS <- paste("lim_date_", ANALYSIS, sep = "")
  x_scale <- scale_x_date(
    date_breaks = "1 week", date_labels = "%b. %d",
    limits = c(MIN_DATE, MAX_DATE),
    expand = c(0, 0))
} else {
  x_scale <- scale_x_date(
    date_breaks = "1 week", date_labels = "%b. %d",
    limits = c(as.Date("2020-02-24"), MAX_DATE),
    expand = c(0, 0))
}

re_plot <- ggplot() + 
  geom_line(
    data = bdsky_re_data, 
    aes(x = Date - 7, y = bdsky_median), color = bdsky_case_color) +  # skyline intervals are logged by the week's end date
  geom_ribbon(
    data = bdsky_re_data, 
    aes(x = Date - 7, ymin = bdsky_low, ymax = bdsky_high), 
    alpha = 0.5, fill = bdsky_case_color) + 
  geom_line(
    data = line_list_re_data, 
    aes(x = date, y = line_list_median), color = line_list_re_color) + 
  geom_ribbon(
    data = line_list_re_data, 
    aes(x = date, ymin = line_list_low, ymax = line_list_high), 
    alpha = 0.5, fill = line_list_re_color) + 
  x_scale + 
  geom_hline(yintercept = 1, linetype = 3) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  labs(x = element_blank(), y = expression(R[e])) + 
  scale_y_continuous(limits = c(0, RE_Y_LIM))
# show(re_plot)

# 2nd facet: Confirmed vs. estimated cases through time ------------------------

if (!is.na(MIN_DATE)) {
  plot_data <- bdsky_sampling_prop_data %>% filter(Date >= MIN_DATE)
} else {
  plot_data <- bdsky_sampling_prop_data
}

case_plot <- ggplot(data = plot_data) + 
  geom_line(
    aes(x = Date - 3.5, y = n_cases), color = conf_case_color) + 
  geom_line(
    aes(x = Date - 3.5, y = n_cases_median), color = bdsky_case_color) + 
  geom_ribbon(
    aes(x = Date - 3.5, ymin = n_cases_low, ymax = n_cases_high), 
    fill = bdsky_case_color, alpha = 0.5) + 
  scale_y_continuous(limits = c(0, 2500)) + 
  x_scale + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  labs(x = element_blank(), y = "No. infected\nindividuals")
# show(case_plot)

# For our own information: Undersampling through time --------------------------

undersampling_data <- bdsky_sampling_prop_data %>% mutate(
  per_underreporting_median = 100 * (n_cases_median - n_cases) / n_cases_median,
  per_underreporting_high = 100 * (n_cases_high - n_cases) / n_cases_high,
  per_underreporting_low = 100 * (n_cases_low - n_cases) / n_cases_low)

if (!is.na(MIN_DATE)) {
  underreporting_plot_data <- undersampling_data %>% filter(Date >= MIN_DATE)
} else {
  underreporting_plot_data <- undersampling_data
}

underreporting_plot <- ggplot(data = underreporting_plot_data) + 
  geom_line(
    aes(x = Date - 3.5, y = per_underreporting_median), color = underreporting_color) + 
  geom_ribbon(
    aes(x = Date - 3.5, ymin = per_underreporting_low, ymax = per_underreporting_high), 
    fill = underreporting_color, alpha = 0.5) + 
  geom_hline(yintercept = 0) + 
  x_scale + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  geom_vline(xintercept = as.Date("2020-04-27"), linetype = 3) +  # testing of all symptomatic patients resumed
  geom_vline(xintercept = as.Date("2020-06-25"), linetype = 3) +  # tests became free
  labs(x = element_blank(), y = "% Underreporting")
# show(underreporting_plot)

# 2nd facet: Sampling proportion through time ----------------------------------

if (!is.na(MIN_DATE)) {
  bdsky_sampling_prop_plot_data <- bdsky_sampling_prop_data %>% filter(Date >= MIN_DATE)
} else {
  bdsky_sampling_prop_plot_data <- bdsky_sampling_prop_data
}

sampling_prop_plot <- ggplot(data = bdsky_sampling_prop_plot_data) + 
  geom_line(
    aes(x = Date - 3.5, y = bdsky_median), color = bdsky_re_color) + 
  geom_line(
    aes(x = Date, y = n_seqs / n_cases), color = conf_case_color) + 
  geom_ribbon(
    aes(x = Date - 3.5, ymin = bdsky_low, ymax = bdsky_high), 
    fill = bdsky_re_color, alpha = 0.5) + 
  x_scale + 
  scale_y_continuous(limits = c(0, SAMP_PROP_Y_LIM)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  geom_vline(xintercept = as.Date("2020-04-27"), linetype = 3) +  # testing of all symptomatic patients resumed
  geom_vline(xintercept = as.Date("2020-06-25"), linetype = 3) +  # tests became free
  labs(x = element_blank(), y = "Sampling proportion")
show(sampling_prop_plot)

# For comparision of results: with sampling proportion through time ------------

library(gtable)
library(grid)
g1 <- ggplotGrob(re_plot + theme(axis.text.x = element_blank()) + ggtitle(COMP_TITLE))
g2 <- ggplotGrob(sampling_prop_plot + theme(axis.text.x = element_blank()))
g3 <- ggplotGrob(underreporting_plot)
g <- rbind(g1, g2, g3, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths)

png(
  file = paste(OUTDIR, "/", ANALYSIS, ".", DATASET, "_figure_4_comparison.png", sep = ""),
  width = 4, height = 6, units = "in", res = 300)
grid.newpage()
grid.draw(g)
dev.off()

# Sampling proportion prior ----------------------------------------------------

samp_prop_metaprior <- ggplot() +
  geom_area(
    stat = "function",
    fun = dbeta, 
    fill = "grey", 
    aes(x = c(0, SAMP_PROP_Y_LIM)),
    args = list(shape1 = SAMP_PROP_PRIOR_BETA_SHAPE_1, 
                shape2 = SAMP_PROP_PRIOR_BETA_SHAPE_2),
    colour = "grey") + 
  coord_flip() + 
  labs(x = element_blank(), y = element_blank()) + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + 
  scale_x_continuous(limits = c(0, SAMP_PROP_Y_LIM))

show(samp_prop_metaprior)

re_metaprior <- ggplot() +
  geom_area(
    stat = "function",
    fun = dlnorm, 
    fill = "grey",
    aes(x = c(0, RE_Y_LIM)), 
    args = list(meanlog = RE_PRIOR_LOGNORM_MEAN, 
                sdlog = RE_PRIOR_LOGNORM_SD),
    colour = "grey") + 
  coord_flip() + 
  labs(x = element_blank(), y = element_blank()) + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank())

show(re_metaprior)

library(gtable)
library(grid)
g1 <- ggplotGrob(sampling_prop_plot + ggtitle("a)"))
g2 <- ggplotGrob(samp_prop_metaprior + ggtitle("b)"))
g <- cbind(g1, g2, size = "last")
g$heights <- unit.pmax(g1$heights, g2$heights)

png(
  file = paste(OUTDIR, "/", ANALYSIS, ".", DATASET, "_sampling_proportion.png", sep = ""),
  width = 6.5, height = 2.5, units = "in", res = 300)
grid.draw(g)
dev.off()

# For main text figure: with sampling proportion, priors -----------------------

png(
  file = paste(OUTDIR, "/", ANALYSIS, ".", DATASET, "_figure_4.png", sep = ""),
  width = 5, height = 3.5, units = "in", res = 300)
egg::ggarrange(
  re_plot + theme(axis.text.x = element_blank()) + ggtitle("a)") , 
  re_metaprior + ggtitle("c)"), 
  sampling_prop_plot + ggtitle("b)") + labs(x = "Beginning of sampling week"), 
  samp_prop_metaprior + ggtitle("d)"),
  nrow =2,
  widths = c(4, 1))
dev.off()



