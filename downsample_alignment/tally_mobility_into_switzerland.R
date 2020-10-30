# Get data on arrivals in Switzerland for various countries over time.
# Result: data frame with # border crossings into CH from various countries each month
# Method: I have data on the number of tourist arrivals by month from various 
# countries and the number of cross-border work permits by quarter from the countries
# surrounding Switzerland. I tally arrivals per country per month and add 1 for every existing
# work permit for Switzerland for that country (i.e. each work permit counts the 
# same as 1 tourist arrival). Since the data only goes until July, I assume constant mobility 
# from July through present.

require(tidyr)
require(ggplot2)
require(lubridate)
require(dplyr)
require(argparse)

# WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/analyses/2020-09-29_ch_cluster_analysis"
# TOURIST_ARRIVAL_DATA <- paste(WORKDIR, "/data/est_imports/FSO_tourist_arrival_statistics_clean.csv", sep = "/")
# CROSS_BORDER_COMMUTER_DATA <- paste(WORKDIR, "data/est_imports/FSO_grenzgaenger_statistics_clean.csv", sep = "/")
# OUTDIR_DATA <- paste(WORKDIR, "data/est_imports", sep = "/")
# OUTDIR <- paste(WORKDIR, "figures/est_imports", sep = "/")

parser <- argparse::ArgumentParser()
parser$add_argument("--tourists", type="character")
parser$add_argument("--commuters", type="character")
parser$add_argument("--outdirdata", type="character")
parser$add_argument("--outdirfigs", type="character")

args <- parser$parse_args()

TOURIST_ARRIVAL_DATA <- args$tourists
CROSS_BORDER_COMMUTER_DATA <- args$commuters
OUTDIR_DATA <- args$outdirdata
OUTDIR <- args$outdirfigs

system(command = paste("mkdir -p", OUTDIR_DATA))
system(command = paste("mkdir -p", OUTDIR))

tourists <- read.delim(file = TOURIST_ARRIVAL_DATA, sep = ",")
commuters <- read.delim(file = CROSS_BORDER_COMMUTER_DATA, sep = ",")

# Massage data
timeframe_row_filter <- !(tourists$year == 2019 & tourists$month < 12) & 
  !(tourists$year == 2020 & tourists$month > 7)
tourists_2 <- tourists[timeframe_row_filter, ]

# Copy-paste July to all months until the present
july_data <- tourists_2[tourists_2$year == 2020 & tourists_2$month == 7, ]
todays_month <- strsplit(x = as.character(Sys.Date()), split = "-")[[1]][2]
july_data_index <- which(tourists_2$year == 2020 & tourists_2$month == 7)
fill_data <- tourists_2[rep(july_data_index, each = as.numeric(todays_month) - 7), ]
fill_data$month <- (7 + 1):todays_month
month_names <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
fill_data$month_name <- month_names[fill_data$month]
tourists_3 <- rbind(tourists_2, fill_data)

tourists_final <- tourists_3 %>%
  pivot_longer(
    cols = 4:ncol(tourists_3), 
    names_to = "source_loc", 
    values_to = "n_tourist_arrivals") %>%
  select(-c("month_name")) %>%
  mutate(
    quarter = case_when(
      month %in% 1:3 ~ 1,
      month %in% 4:6 ~ 2,
      month %in% 7:9 ~ 3,
      month %in% 10:12 ~ 4),
    n_tourist_arrivals = as.numeric(as.character(n_tourist_arrivals)))

tourists_final$source_country <- gsub(
  tourists_final$source_loc, pattern = "\\.", replacement = " ")

# Copy-paste 2020 Q2 to 2020 Q3 and Q4
commuters_q3 <- commuters[, 4]
commuters_q4 <- commuters[, 4]
commuters_2 <- cbind(commuters, commuters_q3, commuters_q4)
colnames(commuters_2) <- c(colnames(commuters), "X2020Q3", "X2020Q4")

commuters_3 <- commuters_2 %>%
  pivot_longer(
    cols = 2:6,
    names_to = "quarter",
    values_to = "n_commuter_permits") %>%
  filter(X != "Andere") %>%
  mutate(
    source_country = recode(
      .x = X,
      Deutschland = "Germany",
      Frankreich = "France",
      Oesterreich = "Austria",
      Italien = "Italy"))
commuters_3$quarter <- gsub(
  commuters_3$quarter, pattern = "X", replacement = "")
commuters_final <- commuters_3 %>%
  separate(
    col = quarter, 
    into = c("year", "quarter"), 
    sep = "Q")

arrivals <- merge(
  x = tourists_final %>% select(-c(source_loc)),
  y = commuters_final %>% select(-c(X)),
  all = T)

arrivals$n_tourist_arrivals[is.na(arrivals$n_tourist_arrivals)] <- 0
arrivals$n_commuter_permits[is.na(arrivals$n_commuter_permits)] <- 0

arrivals <- arrivals %>%
  mutate(n_arrivals = n_tourist_arrivals + n_commuter_permits)

# Write out results
write.table(
  x = arrivals,
  file = paste(OUTDIR_DATA, "travel_per_country_month.txt", sep = "/"),
  row.names = F, col.names = T, quote = F, sep = "\t")

# Plot arrivals over time for sanity-check & supplemental
arrivals_2 <- arrivals %>%
  pivot_longer(
    cols = c("n_tourist_arrivals", "n_commuter_permits"),
    values_to = "n_persons", 
    names_to = "travel_type",
    names_prefix = "n_") %>%
  mutate(
    source_location = case_when(
      source_country %in% c("Italy", "Germany", "France", "Austria", "United States", "China", "United Kingdom") ~ source_country,
      T ~ "other"))

arrivals_2$source_location <- factor(
  x = arrivals_2$source_location,
  levels = c("France", "Germany", "Italy", "Austria", "other", "United States", "China", "United Kingdom"))

arrivals_2 <- arrivals_2 %>% 
  mutate(
    caveat = ifelse(
      test = (year = 2020 & month %in% (7 + 1):todays_month),
      yes = "data copied from last recorded month",
      no = "data from FSO"))

p <- ggplot(
  data = arrivals_2,
  aes(x = paste(year, month, sep = "-"), y = n_persons, alpha = caveat)) +
  scale_alpha_manual(values = c(0.3, 1)) + 
  geom_col(aes(fill = source_location)) + 
  facet_grid(. ~ travel_type) + 
  theme_bw() + 
  # scale_x_date(date_breaks = "1 month", date_labels = "%b") + 
  labs(x = element_blank(), y = "No. individuals") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

png(
  file = paste(OUTDIR, "travel_per_country_month.png", sep = "/"), 
  width = 6.5, height = 4, units = "in", res = 300)
show(p)
dev.off()

     