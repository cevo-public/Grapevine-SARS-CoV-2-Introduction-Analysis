# Create date to date_trun week lookup table for Tim.
# Will allow us to have the same Re breakpoints as sampling breakpoints.

source("database/R/utility.R")
source("utility_functions.R")

outdir <- "/Users/nadeaus/Repos/cov-swiss-phylogenetics/results_all/jan-dec_-01_max_sampling_1_travel_1_sim_context-sf_111_travel-wt/output/transmission_chain_alignments"

# Connect to local version of database
db_connection <- open_database_connection("local")

# Generate date data
min_date <- as.Date("2020-01-01")
max_date <- as.Date("2020-12-31")
dates <- data.frame(date = seq.Date(from = min_date, to = max_date, by = "day"))

# Create temporary table for date mapping
tbl_name <- "date_to_week"
DBI::dbCreateTable(
  conn = db_connection,
  name = tbl_name,
  fields = dates
)

DBI::dbWriteTable(
  conn = db_connection,
  name = tbl_name, 
  value = dates,
  overwrite = T
)

DBI::dbReadTable(
  conn = db_connection,
  name = tbl_name
)

# Generate mapping
dates_with_week <- dplyr::tbl(db_connection, tbl_name) %>%
  mutate(week = date_trunc('week', date)) %>%
  collect() %>%
  mutate(week_orig = week,
         week = format(week_orig, "%Y-%m-%d"))

# Write out mapping
write.csv(
  x = dates_with_week %>% select(date, week),
  file = paste(outdir, "date_to_week.csv", sep = "/"),
  row.names = F
)

# Remove table
DBI::dbRemoveTable(
  conn = db_connection,
  name = tbl_name
)
